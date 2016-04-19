#include "simple_remeshing.h"
#include "WeavingObject.h"

#include "Fenghui_Zhang_core/logging.h"
#include "Fenghui_Zhang_core//object_store.h"
#include "Fenghui_Zhang_core/vertex.h"
#include "Fenghui_Zhang_core/edge.h"
#include "Fenghui_Zhang_core/face.h"

#include "cyPoint.h"
#include <cmath>
#include <iostream>



// The following vector functions are more efficient than vector class, in terms
// of time and space. We aren't sure if we really need a math vector class.
// TODO: move these to a separate file so they can be shared.
void VectorAddition(const float* u, const float* v, float target[]) {
  for (int i = 0; i < 3; ++i) {
    target[i] = u[i] + v[i];
  }
}

bool VectorDivide(const float* u, float v, float p[]) {
  if (std::abs(v) < 0.000000001) return false;
  for (int i = 0; i < 3; ++i) {
    p[i] = u[i] / v;
  }
  return true;
}

// Triangulate a face by adding a vertex at the center of the face and connecting 
// the new vertex to all the vertices on the face.
void TriangulateFace(WeavingObject* obj, Face *f) {
  LOGLINE("Spliting face " << f->GetID());
  float *p = new float[3];
  p[0] = p[1] = p[2] = 0.0;
  Vertex* v =  obj->AddVertex();
  std::vector<Vertex*> face_vertices = f->GetVertices();
  std::vector<Edge*> face_edges = f->GetEdges();
  Edge* first_edge = NULL;
  for (int jj = 0; jj < face_vertices.size(); ++jj) {
    // p = p + obj->GetVertexCoordinates(face_vertices[jj]);
    VectorAddition(p, obj->GetVertexCoordinates(face_vertices[jj]), p);
    // We need edges from the face.
    Edge* new_edge = obj->AddEdge(
        face_vertices[jj],
        v,
        face_edges[jj == 0 ? face_edges.size() - 1 : jj - 1],
        first_edge);
    if (first_edge == NULL) first_edge = new_edge;
  }
  // Something is wrong if there is no vertex in the face.
  // We assume the face is at least a two-gon.
  if (face_vertices.size() < 2) {
    delete[] p;
    return;
  }
  // p = p / face_vertices.size();
  VectorDivide(p, face_vertices.size(), p);
  // VectorDivide(p, 0.8, p);

  // obj keeps p from now on.
  obj->SetVertexCoordinates(v, p);
}

// Triangulate all faces.
// TODO: write a test for this function.
void SimpleRemeshing(WeavingObject* obj) {
  // The caller has to make sure the faces stored in obj is up to date.
  std::set<Edge*> old_edges = obj->GetEdges();
  std::set<Face*> faces = obj->GetFaces();
  for (std::set<Face*>::iterator it = faces.begin(); it != faces.end(); ++it) {
    TriangulateFace(obj, *it);
  }
  /*
  // Remove old edges.
  for (std::set<Edge*>::iterator it =  old_edges.begin();
       it != old_edges.end(); ++it) {
    obj->RemoveEdge(*it);
  }
  */

  obj->ReComputeFaces();
  obj->ClearFaceNormals();
}

void DooSabin(WeavingObject* obj) {
  static float ratio = 0.5;
  float *p = new float[3];
  float *q;
  std::set<Face*> faces = obj->GetFaces();
  std::set<Edge*> old_edges = obj->GetEdges();
  std::set<Vertex*> old_vertices = obj->GetVertices();
  for (std::set<Face*>::iterator it = faces.begin(); it != faces.end(); ++it) {
    p[0] = p[1] = p[2] = 0.0;
    Face* f = *it;
    std::vector<Vertex*> face_vertices = f->GetVertices();
    std::vector<Vertex*> new_vertices;
    std::vector<Edge*> new_edges;
    // Get the centroid of the face.
    for (int j = 0; j < face_vertices.size(); ++j) {
      Vertex* u = face_vertices[j];
      VectorAddition(p, obj->GetVertexCoordinates(u), p);
      // Add a vertex to the obj.
      Vertex* v =  obj->AddVertex();
      new_vertices.push_back(v);
      obj->SetVertexCoordinates(v, p);
      LOGLINE("Pair: " << u->GetID() << "," << v->GetID());
    }
    VectorDivide(p, face_vertices.size(), p);
    // p is the centroid  of the face now.
    // LOGLINE("p:" << p[0] << "," << p[1] << "," << p[2]);
    int cnt = new_vertices.size();
    for (int j = 0; j < cnt; ++j) {
      Vertex* u = new_vertices[j];
      Vertex* v = face_vertices[j];
      // Create a cycle for the new vertices.
      Vertex* w = j < cnt - 1 ? new_vertices[j + 1] : new_vertices[0];
      new_edges.push_back(
          obj->AddEdge(
              u,
              w,
              j == 0 ? NULL : new_edges[j - 1],
              j < cnt - 1 ? NULL : new_edges[0]));
      // LOGLINE("AddEdge: " << u->GetID() << " - " << w->GetID());
      // q = -p;
      q = new float[3];
      VectorDivide(p, -1, q);
      // q = u + q;
      VectorAddition(q, obj->GetVertexCoordinates(v), q);
      // q = q / ratio;
      VectorDivide(q, 1 / ratio, q);
      // q = p + q;
      VectorAddition(p, q, q);
      // Now obj owns q.
      obj->SetVertexCoordinates(u, q);
    }
    
    // Connect new cycles withe paired faces (cycles). This will give us the
    // rotation of the newly added vertices around each old vertex.
    std::vector<Edge*> face_edges = f->GetEdges();
    for (int j = 0; j < cnt; ++j) {
      // Add edge between u and v.
      Vertex* u = new_vertices[j];
      Vertex* v = face_vertices[j];
      obj->AddEdge(u, v, new_edges[j],
          j == 0 ? face_edges[face_edges.size() - 1] : face_edges[j - 1]);
    }
  }

  /**
   * Each cycle has a "spin" when created, all we need to do is to make sure the
   * new edges conform to that spin.
   * We need to know who are paired with an old vertex, and their order wsp to
   * the rotation.
   */

  // Remove old edges.
  for (std::set<Edge*>::iterator it =  old_edges.begin();
       it != old_edges.end(); ++it) {
    obj->RemoveEdge(*it);
  }

  // Connect newly added cycles. 
  for (std::set<Vertex*>::iterator it =  old_vertices.begin();
       it != old_vertices.end(); ++it) {
    Vertex* u = *it;
    LOGLINE(u->GetID() << ":" << u->PrintRotation());
    std::vector<Edge*> rotation = u->GetRotation();
    for (int i = 0; i < rotation.size(); ++i) {
      Edge* ev = rotation[i];
      Vertex* v = ev->GetOtherEnd(u);
      // Edge* ew = rotation[i == rotation.size() - 1 ? 0 : i + 1];
      Edge* ew = u->GetPreviousEdgeInRotation(ev);
      Vertex* w = ew->GetOtherEnd(u);
      obj->AddEdge(v, w, ev, w->GetPreviousEdgeInRotation(ew));
      // obj->AddEdge(v, w, v->GetPreviousEdgeInRotation(ev), ew);
    }
  }

  // Remove old vertices
  for (std::set<Vertex*>::iterator it =  old_vertices.begin();
       it != old_vertices.end(); ++it) {
    obj->RemoveVertex(*it, true);
  }

  /*
  // Print vertex locations
  std::set<Vertex*> vertices = obj->GetVertices();
  for (std::set<Vertex*>::iterator it = vertices.begin();
       it != vertices.end(); ++it) {
    LOGLINE((*it)->GetID() << "["
              << obj->GetVertexCoordinates(*it)[0] << ","
              << obj->GetVertexCoordinates(*it)[1] << ","
              << obj->GetVertexCoordinates(*it)[2] << "]"
              );
  }
  // */
  
  delete[] p;

  obj->ReComputeFaces();
  obj->ClearFaceNormals();
}


/**
 * Apply Catmull_Clark subdivision on obj.
 */
void Catmull_Clark(WeavingObject* obj) {
	LinearSub(obj);
	Average(obj);
}

std::pair<Vertex*, Vertex*> make_orderless_pair(Vertex* v0, Vertex* v1) {
	if ( v0->GetID() < v1->GetID() ) {
		return std::make_pair(v0, v1);
	} 
	else {
		return std::make_pair(v1, v0);
	}
}

Edge* AddNewEdge(std::map<std::pair<Vertex*, Vertex*>, Edge*> &Verts_Edge,
				 Vertex* V0, Vertex* V1, WeavingObject* obj, int twist) {
	Edge* newEdge;
	std::map<std::pair<Vertex*, Vertex*>, Edge*>::iterator ei= Verts_Edge.find(make_orderless_pair(V0, V1));
	if ( Verts_Edge.end() == ei ) {		// edge not created yet
		// create new edge_edge
		newEdge = ObjectStore::GetInstance()->CreateEdge(V0, V1);
		obj->InsertEdges(newEdge);		// just insert to the edge set
		Verts_Edge[make_orderless_pair(V0, V1)] = newEdge;		// insert to harsh
		// set edge twist number
		obj->edgePro[newEdge].twist_ = twist;
	}
	else {
		newEdge = (*ei).second;
	}

	return newEdge;
}

void LinearSub(WeavingObject* obj) {
	Vertex* startVert;
	Vertex* endVert;
	std::map<std::pair<Vertex*, Vertex*>, Edge*> verts_edge;
	std::map<Edge*, EdgeInfo*> old_EdgeInfo;		// all old edges
	EdgeInfo* pre_EdgeInfo = NULL;
	EdgeInfo* this_EdgeInfo = NULL;

	// 1. add edge middle point   ---   old vertices remain
	std::set<Edge*> oldEdges = obj->GetEdges();
	Edge* thisEdge;
	Vertex* midVert;
	cyPoint3f startPos;
	cyPoint3f endPos;
	for (std::set<Edge*>::iterator ei = oldEdges.begin();
		ei != oldEdges.end(); ei++) {
			thisEdge = *ei;
			midVert = obj->AddVertex();
			startVert = thisEdge->GetStart();
			endVert = thisEdge->GetEnd();
			startPos.Set(obj->GetVertexCoordinates(startVert));
			endPos.Set(obj->GetVertexCoordinates(endVert));
			startPos = (startPos + endPos) / 2.0;
			float* pos = new float[3];
			startPos.GetValue(pos);
			obj->SetVertexCoordinates(midVert, pos);
			// insert to hash
			verts_edge[make_orderless_pair(startVert, endVert)] = thisEdge;
			this_EdgeInfo = new EdgeInfo(startVert, endVert, midVert, NULL, NULL, NULL, NULL, NULL, NULL, obj->edgePro[thisEdge].twist_);
			old_EdgeInfo[thisEdge] = this_EdgeInfo;

			// delete old edge
			ObjectStore::GetInstance()->DeleteEdge(thisEdge);
	}

	// 2. delete old edges + edge properties
	obj->ClearEdges();
	obj->edgePro.clear();
	// delete vertex rotation
	std::set<Vertex*> oldVerts = obj->GetVertices();
	for ( std::set<Vertex*>::iterator vi = oldVerts.begin();
		vi != oldVerts.end(); vi++ ) {
			Vertex* v = *vi;
			v->ClearRotations();
	}

	// 3. add new edges: edge_edge / face_edge
	std::set<Face*> faces = obj->GetFaces();
	Face* thisFace;
	Vertex* thisVert;
	Vertex* preVert;
	Vertex* nextVert;
	Vertex* Vc;		// face center
	float* centerPos;
	Edge* pre_OldEdge;
	Edge* this_OldEdge;
	Edge* newEdge;
	Edge* firstEdge;
	Edge* preEdge;
	int newEdge_twist;
	for (std::set<Face*>::iterator fi = faces.begin();
		fi != faces.end(); fi++) {
			thisFace = *fi;
			std::vector<Vertex*> faceVerts = thisFace->GetVertices();
			// add face center
			Vc = obj->AddVertex();
			centerPos = new float[3];
			(obj->GetFaceCenter(thisFace, false)).GetValue(centerPos);
			obj->SetVertexCoordinates(Vc, centerPos);
			// 1st vertex in this face
			preVert = faceVerts[faceVerts.size() - 1];
			thisVert = faceVerts[0];
			nextVert = faceVerts[1];
			// query old edge mid vertex + old edge twist number
			pre_OldEdge = verts_edge[make_orderless_pair(preVert, thisVert)];
			pre_EdgeInfo = old_EdgeInfo[pre_OldEdge];
			this_OldEdge = verts_edge[make_orderless_pair(thisVert, nextVert)];
			this_EdgeInfo = old_EdgeInfo[this_OldEdge];
			for (int i = 0; i < faceVerts.size(); i++) {
				// each face incident to a vertex is a quad
				// 1st edge
				// 1/2 of pre_OldEdge
				if (0 == pre_EdgeInfo->twist) {		// bug: not -1 twist case
					newEdge_twist = 1;
				} 
				else {
					newEdge_twist = pre_EdgeInfo->twist;
				}
				newEdge = AddNewEdge(verts_edge, pre_EdgeInfo->mid, thisVert, obj, newEdge_twist);
				// add orientation info for yarn start edge --- set twice
				if (obj->yarnStartEdges_set.find(pre_OldEdge) != obj->yarnStartEdges_set.end()) {
					if (thisVert->GetID() == pre_EdgeInfo->s->GetID()) {
						pre_EdgeInfo->e1 = newEdge;
						pre_EdgeInfo->right = Vc;
					} 
					else if (thisVert->GetID() == pre_EdgeInfo->t->GetID()) {
						pre_EdgeInfo->e2 = newEdge;
						pre_EdgeInfo->left = Vc;
					} 
					else {
						std::cout << "error: edge ends are not right" << std::endl;
					}
				}
				preEdge = newEdge;
				firstEdge = newEdge;		// remember 1st edge in the face incident to ith vertex
				// 2nd edge
				// 1/2 of this_OldEdge
				if (0 == this_EdgeInfo->twist) {
					newEdge_twist = 1;
				} 
				else {
					newEdge_twist = this_EdgeInfo->twist;
				}
				newEdge = AddNewEdge(verts_edge, thisVert, this_EdgeInfo->mid, obj, newEdge_twist);
				thisVert->InsertEdgeInRotation_load(newEdge, preEdge);		// add corner to rotation system
				preEdge = newEdge;
				// add orientation info for yarn start edge --- set twice
				if (obj->yarnStartEdges_set.find(this_OldEdge) != obj->yarnStartEdges_set.end()) {
					if (thisVert->GetID() == this_EdgeInfo->s->GetID()) {
						this_EdgeInfo->e1 = newEdge;
						this_EdgeInfo->left = Vc;
					} 
					else if (thisVert->GetID() == this_EdgeInfo->t->GetID()) {
						this_EdgeInfo->e2 = newEdge;
						this_EdgeInfo->right= Vc;
					} 
					else {
						std::cout << "error: edge ends are not right" << std::endl;
					}
				}
				// 3rd edge
				// 3/4 of this_OldEdge
				if (0 == this_EdgeInfo->twist) {
					newEdge_twist = 2;
				} 
				else {
					newEdge_twist = 1;
				}
				newEdge = AddNewEdge(verts_edge, this_EdgeInfo->mid, Vc, obj, newEdge_twist);
				this_EdgeInfo->mid->InsertEdgeInRotation_load(newEdge, preEdge);		// add corner to rotation system
				preEdge = newEdge;
				// add orientation info for yarn start edge --- set twice
				if (obj->yarnStartEdges_set.find(this_OldEdge) != obj->yarnStartEdges_set.end()) {
					if (thisVert->GetID() == this_EdgeInfo->s->GetID()) {
						this_EdgeInfo->e3 = newEdge;
					} 
					else if (thisVert->GetID() == this_EdgeInfo->t->GetID()) {
						this_EdgeInfo->e4 = newEdge;
					} 
					else {
						std::cout << "error: edge ends are not right" << std::endl;
					}
				}
				// 4th edge
				// 3/4 of pre_OldEdge
				if (0 == pre_EdgeInfo->twist) {
					newEdge_twist = 2;
				} 
				else {
					newEdge_twist = 1;
				}
				newEdge = AddNewEdge(verts_edge, Vc, pre_EdgeInfo->mid, obj, newEdge_twist);
				Vc->InsertEdgeInRotation_load(newEdge, preEdge);		// add corner to rotation system
				// add orientation info for yarn start edge --- set twice
				if (obj->yarnStartEdges_set.find(pre_OldEdge) != obj->yarnStartEdges_set.end()) {
					if (thisVert->GetID() == pre_EdgeInfo->s->GetID()) {
						pre_EdgeInfo->e4 = newEdge;
					} 
					else if (thisVert->GetID() == pre_EdgeInfo->t->GetID()) {
						pre_EdgeInfo->e3 = newEdge;
					} 
					else {
						std::cout << "error: edge ends are not right" << std::endl;
					}
				}
				// add corner of 1st and 4th edges
				pre_EdgeInfo->mid->InsertEdgeInRotation_load(firstEdge, newEdge);
				// update
				preVert = thisVert;
				thisVert = nextVert;
				nextVert = faceVerts[(i + 2) % faceVerts.size()];		// i has not incremented yet
				pre_OldEdge = this_OldEdge;
				pre_EdgeInfo = this_EdgeInfo;
				this_OldEdge = verts_edge[make_orderless_pair(thisVert, nextVert)];
				this_EdgeInfo = old_EdgeInfo[this_OldEdge];
			}	
	}

	// produce new yarn start corners
	std::vector<TraceCorner>* newCorners = new std::vector<TraceCorner>;
	obj->yarnStartEdges_set.clear();
	TraceCorner corner;
	for (std::vector<TraceCorner>::iterator ci = obj->yarnStartConters->begin();
		ci != obj->yarnStartConters->end(); ci++) {
			this_EdgeInfo = old_EdgeInfo[ci->e];
			if (ci->s->GetID() == this_EdgeInfo->s->GetID()) {		// same orientation
				if (0 != this_EdgeInfo->twist) {		// +1 / 2 twist
					// 1
					corner.e = this_EdgeInfo->e1;
					corner.s = this_EdgeInfo->s;
					corner.traceType = ci->traceType;
					newCorners->push_back(corner);
					obj->yarnStartEdges_set.insert(corner.e);
					// 2
					corner.e = this_EdgeInfo->e2;
					corner.s = this_EdgeInfo->mid;
					corner.traceType = ci->traceType;
					newCorners->push_back(corner);
					obj->yarnStartEdges_set.insert(corner.e);
				} 
				else {		// 0 twist
					if (1 == ci->traceType) {
						// 3
						corner.e = this_EdgeInfo->e3;
						corner.s = this_EdgeInfo->left;
						corner.traceType = -1;
						newCorners->push_back(corner);
						obj->yarnStartEdges_set.insert(corner.e);
						// 4
						corner.e = this_EdgeInfo->e4;
						corner.s = this_EdgeInfo->mid;
						corner.traceType = -1;
						newCorners->push_back(corner);
						obj->yarnStartEdges_set.insert(corner.e);
					} 
					else {
						// 4
						corner.e = this_EdgeInfo->e4;
						corner.s = this_EdgeInfo->right;
						corner.traceType = 1;
						newCorners->push_back(corner);
						obj->yarnStartEdges_set.insert(corner.e);
						// 3
						corner.e = this_EdgeInfo->e3;
						corner.s = this_EdgeInfo->mid;
						corner.traceType = 1;
						newCorners->push_back(corner);
						obj->yarnStartEdges_set.insert(corner.e);
					}
				}
			} 
			else if (ci->s->GetID() == this_EdgeInfo->t->GetID()) {	// ------ flip orientation
				if (0 != this_EdgeInfo->twist) {		// +1 / 2 twist
					// 2
					corner.e = this_EdgeInfo->e2;
					corner.s = this_EdgeInfo->t;
					corner.traceType = ci->traceType;
					newCorners->push_back(corner);
					obj->yarnStartEdges_set.insert(corner.e);
					// 1
					corner.e = this_EdgeInfo->e1;
					corner.s = this_EdgeInfo->mid;
					corner.traceType = ci->traceType;
					newCorners->push_back(corner);
					obj->yarnStartEdges_set.insert(corner.e);
				} 
				else {		// 0 twist
					if (1 == ci->traceType) {
						// 4
						corner.e = this_EdgeInfo->e4;
						corner.s = this_EdgeInfo->right;
						corner.traceType = -1;
						newCorners->push_back(corner);
						obj->yarnStartEdges_set.insert(corner.e);
						// 3
						corner.e = this_EdgeInfo->e3;
						corner.s = this_EdgeInfo->mid;
						corner.traceType = -1;
						newCorners->push_back(corner);
						obj->yarnStartEdges_set.insert(corner.e);
					} 
					else {
						// 3
						corner.e = this_EdgeInfo->e3;
						corner.s = this_EdgeInfo->left;
						corner.traceType = 1;
						newCorners->push_back(corner);
						obj->yarnStartEdges_set.insert(corner.e);
						// 4
						corner.e = this_EdgeInfo->e4;
						corner.s = this_EdgeInfo->mid;
						corner.traceType = 1;
						newCorners->push_back(corner);
						obj->yarnStartEdges_set.insert(corner.e);
					}
				}
			} 
			else {
				std::cout << "edge ends are wrong" << std::endl;
			}
	}
	std::cout << "num old yarns " << obj->yarnStartConters->size() << std::endl;
	obj->yarnStartConters->clear();
	delete obj->yarnStartConters;
	obj->yarnStartConters = newCorners;
	
	obj->ReComputeFaces();
	obj->ClearFaceNormals();
	obj->subdivisionLevel += 1;
	std::cout << "num new yarns " << obj->yarnStartConters->size() << std::endl;

	// delete old edge info
	for (std::map<Edge*, EdgeInfo*>::iterator mi = old_EdgeInfo.begin();
		mi != old_EdgeInfo.end(); mi++) {
			delete (*mi).second;
	}
	old_EdgeInfo.clear();
}

void Average(WeavingObject* obj) {
	std::set<Vertex*> vertices = obj->GetVertices();
	std::map<Vertex*, float> face_cout;
	std::map<Vertex*, cyPoint3f> vertPos;
	cyPoint3f faceCenter;
	Vertex* v = NULL;
	for ( std::set<Vertex*>::iterator vi = vertices.begin();
		vi != vertices.end(); vi++ ) {
		v = *vi;
		face_cout[v] = 0.0;
		vertPos[v] = faceCenter;	// zero vector
	}
	
	Face* f = NULL;
	std::set<Face*> faces = obj->GetFaces();
	std::vector<Vertex*> face_verts;
	for ( std::set<Face*>::iterator fi = faces.begin();
		fi != faces.end(); fi++) {
			f = *fi;
			faceCenter.Set(obj->GetFaceCenter(f, false));
			face_verts = f->GetVertices();
			for ( int i = 0; i < face_verts.size(); i++ ) {
				v = face_verts[i];
				vertPos[v] += faceCenter;
				face_cout[v] += 1.0;
			}
	}

	for ( std::set<Vertex*>::iterator vi = vertices.begin();
		vi != vertices.end(); vi++ ) {
			v = *vi;
			vertPos[v] /= face_cout[v];
			obj->ResetVertexCoordinates(v, vertPos[v].x, vertPos[v].y, vertPos[v].z);
	}

	obj->ClearFaceNormals();
}