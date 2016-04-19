#include "WeavingObject.h"

#include "Fenghui_Zhang_core/topology_object.h"
#include "Fenghui_Zhang_core/face.h"
#include "Fenghui_Zhang_core/edge.h"
#include "Fenghui_Zhang_core/vertex.h"
#include "Fenghui_Zhang_core/logging.h"
#include "Fenghui_Zhang_core/object_store.h"

#include "cyPoint.h"
#include "cyMatrix3.h"
#include "ColorSpace.h"
#include "Output.h"

#include <fl/gl.h>
#include "fl/glu.h"

#include <cstdlib>
#include <cstdio>
#include <ctime>
#include <cmath>
#include <string.h>

#include <iostream>
#include <map>
#include <vector>
#include <utility>
#include <algorithm>



// Multiplication - cross product
void CrossProduct(const float* u, const float* v, float p[]) {
	p[0] = u[1]*v[2] - u[2]*v[1];
	p[1] = u[2]*v[0] - u[0]*v[2];
	p[2] = u[0]*v[1] - u[1]*v[0];
}

// interpolate a0 and a1 to get b0 and b1 with constant c
void BiLinearInterpolate(const float c, const cyPoint3f a0, const cyPoint3f a1, cyPoint3f& b0, 
						 cyPoint3f& b1) {
	b0 = (1 - c) * a0 + c * a1;
	b1 = c * a0 + (1 - c) * a1;
}

WeavingObject::WeavingObject(void)
{
	isFaceDirty_ = true;
	weight4Center_ = 1.0;
	renderType_ = BEZIER_RIBBON;
	numOpenYarns_ = 0;
	numClosedYarns_ = 0;
	width_ = 0.5;
	curvature_ = 0.5;
	numSamples_ = 12;
	colorWheel_ = NULL;
	numColors_ = 0;
	displaceFactor_ = 0.05;
	numFacesObjFile_ = 0;
	subdivisionLevel = 0;
	yarnStartConters = NULL;
	isFlaten = false;
	tempWeavingCycle = NULL;
	tempNumCrossings = 0;
	weavingPieces = NULL;
	numWeavingPiece = 0;
	boundingBox_size = 5.0;
	holeRadius = 0.05;
	paperWidth = 16.5;
	paperHeight = 22;
}

WeavingObject::~WeavingObject(void)
{
	// release memory
	for ( int i = 0; i < 12; i++ ) {
		delete [] colorWheel_[i];
	}
	
	for ( std::map<Vertex*, VertexProperty>::iterator it = vertPro_.begin();
		it != vertPro_.end(); ++it ) {
			delete [] (*it).second.coords_;
	}

	for ( std::map<Face*, FaceProperty>::iterator it = facePro_.begin();
		it != facePro_.end(); ++it ) {
			delete [] (*it).second.normal_;
	}
}

bool WeavingObject::Clear() {
	// clear property maps
	for ( std::map<Vertex*, VertexProperty>::iterator it = vertPro_.begin();
		it != vertPro_.end(); ++it ) {
			delete [] (*it).second.coords_;
	}
	vertPro_.clear();
	edgePro.clear();
	for ( std::map<Face*, FaceProperty>::iterator it = facePro_.begin();
		it != facePro_.end(); ++it ) {
			delete [] (*it).second.normal_;
	}
	facePro_.clear();

	// clear vertices, edges, and faces
	for (std::set<Vertex*>::iterator it = vertices_.begin();
		it != vertices_.end(); ++it) {
			ObjectStore::GetInstance()->DeleteVertex(*it);
	}
	vertices_.clear();
	//vertPro_.clear();
	for (std::set<Edge*>::iterator it = edges_.begin();
		it != edges_.end(); ++it) {
			ObjectStore::GetInstance()->DeleteEdge(*it);
	}
	edges_.clear();
	//edgePro.clear();
	for (std::set<Face*>::iterator it = faces_.begin();
		it != faces_.end(); ++it) {
			ObjectStore::GetInstance()->DeleteFace(*it);
	}
	faces_.clear();
	//facePro_.clear();
	face_map_.clear();

	ObjectStore::GetInstance()->Compress();

	isFaceDirty_ = false;
	numOpenYarns_ = 0;
	numClosedYarns_ = 0;
	numFacesObjFile_ = 0;

	subdivisionLevel = 0;
	if (NULL != yarnStartConters) {
		delete yarnStartConters;
		yarnStartConters = NULL;
	}
	yarnStartEdges_set.clear();

	// children of 1st twist-1 edges
	childrenEdges.clear();

	// clear color wheel
	if (NULL != colorWheel_) 	{
		for (int i = 0; i <= numColors_; i++) {
			delete [] colorWheel_[i];
		}
		delete [] colorWheel_;
		colorWheel_ = NULL;
		numColors_ = 0;
	}

	return true;
}

void WeavingObject::init() {
	// initialize random seed
	srand ( time(NULL) );
}

// compute N colors to use
bool WeavingObject::ComputeColorWheel(int N, float minH, float maxH) {
	// clean
	if (NULL != colorWheel_) 	{
		for (int i = 0; i <= numColors_; i++) {
			delete [] colorWheel_[i];
		}
		delete [] colorWheel_;
		colorWheel_ = NULL;
	}
	
	// allocate memory
	numColors_ = N;
	colorWheel_ = new unsigned char*[numColors_+1];
	for ( int i = 0; i <= numColors_; i++ ) {
		colorWheel_[i] = new unsigned char[3];
	}

	float incrementH = (maxH - minH) / float(numColors_ - 1);
	float H = minH;		// S = 1		V = 1
	float R, G, B;
	for (int i = 0; i <= numColors_; i++) {
		if (HSVtoRGB(H, 1, 1, R, G, B)) {
			colorWheel_[i][0] = unsigned char(R);
			colorWheel_[i][1] = unsigned char(G);
			colorWheel_[i][2] = unsigned char(B);
			H += incrementH;
		} 
		else 	{
			return false;
		}
	}
	return true;
// 	complementary colors
// 		// red
// 		colorWheel_[0][0] = 255;
// 		colorWheel_[0][1] = 0;
// 		colorWheel_[0][2] = 0;
// 		// cyan
// 		colorWheel_[1][0] = 0;
// 		colorWheel_[1][1] = 255;
// 		colorWheel_[1][2] = 255;
// 		// orange
// 		colorWheel_[2][0] = 255;
// 		colorWheel_[2][1] = 125;
// 		colorWheel_[2][2] = 0;
// 		// ocean
// 		colorWheel_[3][0] = 0;
// 		colorWheel_[3][1] = 125;
// 		colorWheel_[3][2] = 255;
// 		// yellow
// 		colorWheel_[4][0] = 255;
// 		colorWheel_[4][1] = 255;
// 		colorWheel_[4][2] = 0;
// 		// blue
// 		colorWheel_[5][0] = 0;
// 		colorWheel_[5][1] = 0;
// 		colorWheel_[5][2] = 255;
// 		// spring green
// 		colorWheel_[6][0] = 125;
// 		colorWheel_[6][1] = 255;
// 		colorWheel_[6][2] = 0;
// 		// violet
// 		colorWheel_[7][0] = 125;
// 		colorWheel_[7][1] = 0;
// 		colorWheel_[7][2] = 255;
// 		// green
// 		colorWheel_[8][0] = 0;
// 		colorWheel_[8][1] = 255;
// 		colorWheel_[8][2] = 0;
// 		// magenta
// 		colorWheel_[9][0] = 255;
// 		colorWheel_[9][1] = 0;
// 		colorWheel_[9][2] = 255;
// 		// turquoise
// 		colorWheel_[10][0] = 0;
// 		colorWheel_[10][1] = 255;
// 		colorWheel_[10][2] = 125;
// 		// raspberry
// 		colorWheel_[11][0] = 255;
// 		colorWheel_[11][1] = 0;
// 		colorWheel_[11][2] = 125;
}

std::set<Face*> WeavingObject::GetFaces() {
	return faces_;
}

// Traces all the faces and add them to faces_.
bool WeavingObject::ReComputeFaces() {
	for (std::set<Face*>::iterator it = faces_.begin();
		it != faces_.end(); ++it) {
			ObjectStore::GetInstance()->DeleteFace(*it);
	}
	std::cout << "num old faces: " << faces_.size() << std::endl;
	faces_.clear();
	face_map_.clear();

	std::set<Edge*> visited_edges;
	std::set<Edge*> edges = GetEdges();
// 	for (std::set<Edge*>::iterator it = edges_.begin();
// 		it != edges_.end(); ++it) {
// 			edges.insert((*it));
// 	}

	while(!edges.empty()) {
		Edge* e = *(edges.begin());
		Vertex* u = e->GetStart();
		Vertex* v = e->GetEnd();
		// Make sure u->v is the half edge we haven't visited yet.
		if (face_map_.find(u) != face_map_.end() &&
			face_map_[u].find(v) != face_map_[u].end()) {
				std::swap(u, v);
		}
		// open boundary
		if ( NULL == v->GetNextEdgeInRotation(e) ) {
			if (visited_edges.find(e) != visited_edges.end()) {
				edges.erase(e);
			} else {
				visited_edges.insert(e);
			}
			// The half edge u->v belongs to NULL
			face_map_[u][v] = NULL;
			continue;	// process next edge
		}
		// Trace a face and collect vertices along the way.
		Face* new_face = ObjectStore::GetInstance()->CreateFace();
		std::vector<Vertex*> face_vertices;
		std::vector<Edge*> face_edges;
		Vertex* s = u;
		Vertex* t = v;
		Edge* e_current = e;
		do {
			if (visited_edges.find(e_current) != visited_edges.end()) {
				edges.erase(e_current);
			} else {
				visited_edges.insert(e_current);
			}
			// In a face, the ith vertex is one end of the ith edge
			face_vertices.push_back(s);
			face_edges.push_back(e_current);
			// The half edge s->t belongs to new_face now.
			face_map_[s][t] = new_face;
			// Next edge in rotation after e_current.
			Edge* e_next = t->GetNextEdgeInRotation(e_current);
			// Next vertex after v.
			Vertex* w = e_next->GetOtherEnd(t);
			e_current = e_next;
			s = t;
			t = w;
		} while(s != u || t != v || e_current != e);
		new_face->SetVertices(face_vertices);
		new_face->SetEdges(face_edges);
		faces_.insert(new_face);
	}

	if ( faces_.size() == numFacesObjFile_ ) {
		std::cout << "num new faces: " << faces_.size() << std::endl;
		return true;
	} 
	else {
		//std::cout << "error***********************" << std::endl;
		//std::cout << "number of faces in .obj file: " << numFacesObjFile_ << std::endl;
		std::cout << "number of new faces: " << faces_.size() << std::endl;
		return false;
	}
}

void WeavingObject::SetFaceDirtyFlag() {
	isFaceDirty_ = true;
}

void WeavingObject::ClearFaceNormals() {
	facePro_.clear();
}

void WeavingObject::SetCenterWeight(float w) {
	weight4Center_ = w;
}

void WeavingObject::SetRenderType(RenderType rt)
{
	renderType_ = rt;
}

void WeavingObject::SetWidth(const float w) {
	width_ = w;
}

void WeavingObject::SetCurvature(const float c) {
	curvature_ = c;
}

void WeavingObject::SetDiplaceFactor(const float d) {
	displaceFactor_ = d;
}

void WeavingObject::SetNumSamples(const int s) {
	numSamples_ = s;
}

float* WeavingObject::GetVertexCoordinates(Vertex* v){
	if ( (NULL == v) || vertPro_.find(v) == vertPro_.end() ) {
		return NULL;
	}
	return vertPro_[v].coords_;
}

// Capture the ownership of the coords pointer.
// TODO: release the pointer when it's not needed any more, basically when
//   purge the vertex.
bool WeavingObject::SetVertexCoordinates(Vertex* v, float* coords) {
	if (v == NULL) return false;
	vertPro_[v].coords_ = coords;
	LOGLINE("set position:" << v->GetID()
		<< "[" << coords[0] << "," << coords[1] << "," << coords[2] << "]");
	return true;
}

// Capture the ownership of the coords pointer.
// TODO: release the pointer when it's not needed any more, basically when
//   purge the vertex.
bool WeavingObject::ResetVertexCoordinates(Vertex* v, float x, float y, float z) {
	if (v == NULL) return false;
	vertPro_[v].coords_[0] = x;
	vertPro_[v].coords_[1] = y;
	vertPro_[v].coords_[2] = z;
	LOGLINE("reset position:" << v->GetID()
		<< "[" << coords[0] << "," << coords[1] << "," << coords[2] << "]");
	return true;
}

/**
* Possible problems with normals:
*
*   (1) Face not planar.
*      Need normal at each corner
*   (2) Corner is winged (0 degree or 180)
*      Need to find the nearest nonWinged corner
*   (3) Face in a straight line.
*      Identify it and ignore it.
*
* It's easier just adding all the normals of the corners up, in that case we
* don't have to worry about if it is planar or not.
*
*/
float* WeavingObject::GetFaceNormal(Face* f) {
	if (facePro_.find(f) == facePro_.end()) {
		// Compute normal for the face.
		float U[3] = {0.0, 0.0, 0.0};
		float V[3] = {0.0, 0.0, 0.0};
		float* u_1 = GetVertexCoordinates(f->GetVertices()[0]);
		float* u_2 = GetVertexCoordinates(f->GetVertices()[1]);
		float* u_3 = GetVertexCoordinates(f->GetVertices()[2]);
		U[0] = u_2[0] - u_1[0];
		U[1] = u_2[1] - u_1[1];
		U[2] = u_2[2] - u_1[2];

		V[0] = u_3[0] - u_2[0];
		V[1] = u_3[1] - u_2[1];
		V[2] = u_3[2] - u_2[2];

		float* normal = new float[3];
		CrossProduct(U, V, normal);

		float length = sqrt(normal[0]*normal[0]+normal[1]*normal[1]+normal[2]*normal[2]);
		normal[0] /= length;
		normal[1] /= length;
		normal[2] /= length;
		facePro_[f].normal_ = normal;
	}
	return facePro_[f].normal_;
}

cyPoint3f WeavingObject::GetFaceCenter(Face* f, bool isNullEdgeWeighted) {
	cyPoint3f center;
	float sumWeight = 0;	// sum of weights
	Edge* e;
	Vertex* v;
	std::map<Vertex*, float> vertWeight;
	std::vector<Edge*> face_edges = f->GetEdges();
	if ( isNullEdgeWeighted ) {
		for ( int i = 0; i < face_edges.size(); i++ ) {
			e = face_edges[i];
			if ( 0 == edgePro[e].twist_ ) {
				vertWeight[e->GetStart()] = weight4Center_;
				vertWeight[e->GetEnd()] = weight4Center_;
			}
		}
	}

	std::vector<Vertex*> face_vertices = f->GetVertices();
	for (int i = 0; i < face_vertices.size(); i++) {
		v = face_vertices[i];
		if ( vertWeight.find(v) != vertWeight.end() ) {
			center += vertWeight[v] * cyPoint3f( GetVertexCoordinates(v) );
			sumWeight += vertWeight[v];
		} 
		else {
			center += cyPoint3f( GetVertexCoordinates(v) );
			sumWeight += 1;
		}
	}
	if ( !face_vertices.empty() ) {
		center /= sumWeight;
	}

	return center;
}

bool WeavingObject::ComputeControlPoints(Edge* e, Vertex* s, 
										 const int traceType, 
										 std::vector<cyPoint3f>& Q, 
										 cyPoint3f& edgeNormal)
{
	Vertex* t = e->GetOtherEnd(s);
	int displaceType = traceType * edgePro[e].twist_;
	if ( 2 == edgePro[e].twist_ ) {
		displaceType = 0;
	}
	cyPoint3f cf0, cf1;		// face centers
	cyPoint3f sCoords(GetVertexCoordinates(s));
	cyPoint3f tCoords(GetVertexCoordinates(t));
	float edgeLength = cyPoint3f(tCoords - sCoords).Length();
	cyPoint3f mCoords = (sCoords + tCoords) / 2.0;
	edgeNormal.Zero();

	// get edge normal and face centers
	Face* f = NULL;
	bool isLeftNULL = true;
	bool isRightNULL = true;
	f = face_map_[s][t];
	if ( NULL != f ) {
		edgeNormal += cyPoint3f(GetFaceNormal(f));
		cf0.Set(GetFaceCenter(f, true));
		isLeftNULL = false;
	}
	f = face_map_[t][s];
	if ( NULL != f ) {
		edgeNormal += cyPoint3f(GetFaceNormal(f));
		cf1.Set(GetFaceCenter(f, true));
		isRightNULL = false;
	}
	if ( edgeNormal.Length() > 0) {
		edgeNormal.Normalize();
		// boundary edge: mirror reflect edge region
		if ( isLeftNULL ) {
			cf0 = 2 * mCoords - cf1;
		}
		if ( isRightNULL ) {
			cf1 = 2 * mCoords - cf0;
		}
	}
	else {
		return false;	// edge normal is zero vector
	}

	// displace to get the projection plane
	mCoords += float(displaceType) * edgeNormal * edgeLength * displaceFactor_;

	// project edge area to plane
	cyPoint3f q;
	q = sCoords - mCoords;
	sCoords -= q.Dot(edgeNormal) * edgeNormal;
	q = tCoords - mCoords;
	tCoords -= q.Dot(edgeNormal) * edgeNormal;
	q = cf0 - mCoords;
	cf0 -= q.Dot(edgeNormal) * edgeNormal;
	q = cf1 - mCoords;
	cf1 -= q.Dot(edgeNormal) * edgeNormal;

	// linear interpolate to get control points
	cyPoint3f a0, a1, b0, b1;
	if ( 0 == edgePro[e].twist_ ) {	// untwist edge
		if ( 1 == traceType ) {
			BiLinearInterpolate((1.0-width_)/2, cf0, sCoords, a0, a1);
			BiLinearInterpolate((1.0-width_)/2, cf0, tCoords, b0, b1);
			BiLinearInterpolate((1.0-curvature_)/2, a0, b0, Q[0], Q[4]);
			BiLinearInterpolate((1.0-curvature_)/2, a1, b1, Q[1], Q[5]);
		} 
		else {		// -1 == traceType
			BiLinearInterpolate((1.0-width_)/2, sCoords, cf1, a0, a1);
			BiLinearInterpolate((1.0-width_)/2, tCoords, cf1, b0, b1);
			BiLinearInterpolate((1.0-curvature_)/2, a0, b0, Q[0], Q[4]);
			BiLinearInterpolate((1.0-curvature_)/2, a1, b1, Q[1], Q[5]);
		}
	} 
	else if ( 2 == edgePro[e].twist_ ) {	// NULL edge
		if ( 1 == traceType ) {
			BiLinearInterpolate((1.0-width_)/2, cf0, sCoords, a0, a1);
			BiLinearInterpolate((1.0-width_)/2, cf1, sCoords, b0, b1);
			BiLinearInterpolate((1.0-curvature_)/2, a0, b0, Q[0], Q[4]);
			BiLinearInterpolate((1.0-curvature_)/2, a1, b1, Q[1], Q[5]);
		} 
		else {		// -1 == traceType
			BiLinearInterpolate((1.0-width_)/2, sCoords, cf1, a0, a1);
			BiLinearInterpolate((1.0-width_)/2, sCoords, cf0, b0, b1);
			BiLinearInterpolate((1.0-curvature_)/2, a0, b0, Q[0], Q[4]);
			BiLinearInterpolate((1.0-curvature_)/2, a1, b1, Q[1], Q[5]);
		}
	} 
	else {		// twist number: +1 / -1
		if ( 1 == traceType ) {
			BiLinearInterpolate((1.0-width_)/2, cf0, sCoords, a0, a1);
			BiLinearInterpolate((1.0-width_)/2, tCoords, cf1, b0, b1);
			BiLinearInterpolate((1.0-curvature_)/2, a0, b0, Q[0], Q[4]);
			BiLinearInterpolate((1.0-curvature_)/2, a1, b1, Q[1], Q[5]);
		} 
		else {		// -1 == traceType
			BiLinearInterpolate((1.0-width_)/2, sCoords, cf1, a0, a1);
			BiLinearInterpolate((1.0-width_)/2, cf0, tCoords, b0, b1);
			BiLinearInterpolate((1.0-curvature_)/2, a0, b0, Q[0], Q[4]);
			BiLinearInterpolate((1.0-curvature_)/2, a1, b1, Q[1], Q[5]);
		}
	}
	Q[2] = (Q[0] + Q[4]) / 2.0;
	Q[3] = (Q[1] + Q[5]) / 2.0;

	return true;
}

bool WeavingObject::LoadObjFile(const char *objFileName) {
	// clear object
	if ( !(vertices_.empty()) ) {
		Clear();
	}
	// save file name
	strcpy(fileName_, objFileName);

	// ------------------   load .obj file   ------------------
	FILE* inFile = fopen (objFileName, "r");	// open file
	if ( !inFile ) {
		std::cout << "Can not open file: " << fileName_ << std::endl;
	}
	else {
		std::cout << "model file: " << fileName_ << std::endl;
	}
	// 1st pass - read all the vertices
	char buffer[512] = {0};
	float* coord;
	Vertex* newVertex;
	std::map<int, Vertex*> id_vertex;	// vert ID starts from 1
	while ( fscanf( inFile, "%s", buffer) != EOF ) {
		if ( 'v' == buffer[0] ) {
			switch ( buffer[1] ) {
				case '\0': // v
					coord = new float[3];
					if ( fscanf(inFile, "%f %f %f", &(coord[0]), &(coord[1]), &(coord[2])) != 3 ) {
						delete [] coord;
						std::cout << "vertex " << vertices_.size() << " wrong" << std::endl;
						return false;
					}
					else {
						newVertex = AddVertex();
						SetVertexCoordinates(newVertex, coord);
						id_vertex[vertices_.size()] = newVertex;
						break;
					}

				case 'n':
				case 't':
					break;
			}
		}
	}
	std::cout << "number of vertices: " << vertices_.size() << std::endl;
	// 2nd pass - add edges
	numFacesObjFile_ = 0;
	int tempID, startID, endID;
	Edge* preEdge;
	Edge* newEdge;
	std::pair<int, int> edgeIDpair;	// start <= end
	std::pair<int, int> firstPair;	// 1st edge pair
	std::map<std::pair<int, int>, Edge*> vertex_id_edge;	// vert ID starts from 1
	std::map<std::pair<int, int>, Edge*>::iterator it_edge; 
	rewind(inFile);		// rewind to the beginning of the file
	// process each face
	while ( fscanf(inFile, "%s", buffer) != EOF) {
		if ( 'f' == buffer[0] && '\0' == buffer[1] ) {
			std::vector<int> fVertices;
			while ( fscanf(inFile, "%s", buffer) != EOF ) {
				if ( 1 == sscanf(buffer, "%d", &tempID) ) {
					fVertices.push_back(tempID);
				} 
				else {		// finish one face
					numFacesObjFile_++;
					if ( fVertices.size() < 2 ) {
						std::cout << "face only has " << fVertices.size() << " vertices" << std::endl;
						return false;
					}
					// 1st edge
					startID = fVertices[0];
					endID = fVertices[1];
					if ( startID > endID ) {
						std::swap(startID, endID);
					}
					edgeIDpair.first = startID;
					edgeIDpair.second = endID;
					firstPair = edgeIDpair;
					it_edge = vertex_id_edge.find(edgeIDpair);
					if ( vertex_id_edge.end() == it_edge ) {	// edge not created
						newEdge = ObjectStore::GetInstance()->CreateEdge(id_vertex[startID], id_vertex[endID]);
						edges_.insert(newEdge);
						vertex_id_edge[edgeIDpair] = newEdge;
					}
					else {
						newEdge = it_edge->second;
					}
					preEdge = newEdge;
					for ( int i = 1; i < fVertices.size(); i++ ) {
						startID = fVertices[i];
						endID = fVertices[(i+1 == fVertices.size() ? 0 : i+1)];
						if ( startID > endID ) {
							std::swap(startID, endID);
						}
						edgeIDpair.first = startID;
						edgeIDpair.second = endID;
						it_edge = vertex_id_edge.find(edgeIDpair);
						if ( vertex_id_edge.end() == it_edge) {		// edge not created
							newEdge = ObjectStore::GetInstance()->CreateEdge(id_vertex[startID], id_vertex[endID]);
							edges_.insert(newEdge);
							vertex_id_edge[edgeIDpair] = newEdge;
						}	
						else {
							newEdge = it_edge->second;
						}
						// add to rotation
						id_vertex[fVertices[i]]->InsertEdgeInRotation_load(newEdge, preEdge);
						preEdge = newEdge;
					}
					// last & 1st edges in rotation
					newEdge = vertex_id_edge[firstPair];
					id_vertex[fVertices[0]]->InsertEdgeInRotation_load(newEdge, preEdge);

					fVertices.clear();	// start a new face
					if ( 'f' != buffer[0] ) {
						break;
					}
				} 
			}
		}
	}
	fclose(inFile);		// close file
	std::cout << "number of edges: " << edges_.size() << std::endl;
	// compute faces
	std::cout << "number of faces in .obj file: " << numFacesObjFile_ << std::endl;
	ReComputeFaces();
	isFaceDirty_ = false;
	std::cout << "number of faces after face tracing: " << faces_.size() << std::endl;

	// ------------------   load .twist file   ------------------
	char* pt = strchr(fileName_, '.');
	strcpy(pt, ".twist");
	int twistNumber;
	tempID = 2 * edges_.size();		// edge number
	inFile  = fopen(fileName_, "r");
	if ( !inFile ) {
		std::cout << "Can not open file: " << fileName_ << std::endl;
		RandomTwist();
	}
	else {
		std::cout << "twist file: " << fileName_ << std::endl;
		//
		// edges are saved twice
		// vertex ID starts from 1
		//
		while (!feof(inFile)) {
			if ( fscanf(inFile, "%d %d %d", &startID, &endID, &twistNumber) != 3 ) {
				std::cout << "edge " << startID << " " << endID << " reading is wrong" << std::endl;
				//return false;
			}
			else {
				if ( startID > endID ) {
					std::swap(startID, endID);
				}
				edgeIDpair.first = startID;
				edgeIDpair.second = endID;
				it_edge = vertex_id_edge.find(edgeIDpair);
				if ( vertex_id_edge.end() == it_edge ) {	// edge not found
					std::cout << "edge " << startID << " " << endID << " is not found" << std::endl;
					return false;
				}
				else {
					edgePro[it_edge->second].twist_ = twistNumber;
				}
			}
		}
		fclose(inFile);		// close file
	}

	// 1st 1-twist edge
	std::set<Edge*> edges = GetEdges();
	for ( std::set<Edge*>::iterator it = edges.begin();
		it != edges.end() ; ++it ) {
			Edge* e = *it;
			if (1 == edgePro[e].twist_) {
				childrenEdges.insert(e);
				break;
			}
	}
	//std::cout << "1st twist-1 edge " << childrenEdges_.size() << std::endl;

	return true;	// loading succeeds
}

void WeavingObject::RandomTwist() {
	std::cout << "edge twist numbers are randomly generated" << std::endl;
	int i = 0;
	for ( std::set<Edge*>::iterator ei = edges_.begin();
		ei != edges_.end(); ei++) {
			i++;
//			if ( 10 == i || 9 == i || 11 == i || 12 == i) {
// 			if ( 2 == i || 5 == i) {
// 				edgePro_[*ei].twist_ = 0;
// 			} 
// 			else {
// 				edgePro_[*ei].twist_ = 1;//rand() % 3;		// 0 / 1 / 2
// 			}
// 			i = rand() % 4;		// 0 / 1 / 2
// 			if (i > 0) {
// 				edgePro[*ei].twist_ = 1;
// 			} 
// 			else {
// 				edgePro[*ei].twist_ = 0;
// 			}
			edgePro[*ei].twist_ = 1;//rand() % 3;
	}
}

bool WeavingObject::SaveFiles() {
	// ------------------   write .obj file   ------------------
	char* pt = strchr(fileName_, '.');
	strcpy(pt, "_out.obj");	// out file name: original_out.obj
	FILE * outFile = fopen(fileName_, "wb");
	if ( NULL == outFile ) {
		printf("can not open the outFile %s\n", fileName_);
		return false;
	}
	// vertex list
	std::set<Vertex*> vertices_out = GetVertices();
	int vertID = 1;
	std::map<Vertex*, int> vertex_out_ID;	// vert ID starts from 1
	for (std::set<Vertex*>::iterator it = vertices_out.begin();
		it != vertices_out.end(); ++it) {
			Vertex* v = *it;
			float* coords = GetVertexCoordinates(v);
			fprintf(outFile, "v %f %f %f\n", coords[0], coords[1], coords[2]);
			vertex_out_ID[v] = vertID++;
	}
	fprintf(outFile, "# vertices %d\n", vertices_out.size());
	// face list
	std::set<Face*> faces_out = GetFaces();
	for ( std::set<Face*>::iterator it = faces_out.begin();
		it != faces_out.end(); ++it)
	{
		fprintf(outFile, "f ");
		Face* f = *it;
		std::vector<Vertex*> face_vertices = f->GetVertices();
		for (int j = 0; j < face_vertices.size(); ++j) {
			Vertex* u = face_vertices[j];
			fprintf(outFile, "%d ", vertex_out_ID[u]);
		}
		fprintf(outFile, "\n");
	}
	fprintf(outFile, "# faces %d\n", faces_out.size());
	fclose(outFile);
	std::cout << "save the mesh in " << fileName_ << std::endl;
	
	// ------------   write .twist file   ------------
	// edges are saved twice
	//  vertex ID starts from 1
	pt = strchr(fileName_, '.');
	strcpy(pt, ".twist");	// out file name: original_out.twist
	outFile = fopen(fileName_, "wb");
	if ( NULL == outFile ) {
		printf("can not open the outfile %s\n", fileName_);
		return false;
	}
	// edge list
	std::set<Edge*> edges_out = GetEdges();
	int startVertID, endVert1ID;
	for ( std::set<Edge*>::iterator it = edges_out.begin();
		it != edges_out.end() ; ++it ) {
		Edge* e = *it;
		startVertID = vertex_out_ID[e->GetStart()];
		endVert1ID = vertex_out_ID[e->GetEnd()];
		fprintf(outFile, "%d %d %d\n", startVertID, endVert1ID, edgePro[e]);
		fprintf(outFile, "%d %d %d\n", endVert1ID, startVertID, edgePro[e]);
	}
	std::cout << "save edge twist info in " << fileName_ << std::endl;
	fclose(outFile);

	// ------------   write .start file   ------------
	// children 1st 1-twist edge
	//  vertex ID starts from 1
	pt = strchr(fileName_, '.');
	strcpy(pt, ".start");	// out file name: original_out.twist
	outFile = fopen(fileName_, "wb");
	if ( NULL == outFile ) {
		printf("can not open the outfile %s\n", fileName_);
		return false;
	}
	// edge list
	//std::cout << "1st twist-1 edge " << childrenEdges_.size() << std::endl;
	fprintf(outFile, "%d\n", childrenEdges.size());
	for ( std::set<Edge*>::iterator it = childrenEdges.begin();
		it != childrenEdges.end() ; ++it ) {
			Edge* e = *it;
			startVertID = vertex_out_ID[e->GetStart()];
			endVert1ID = vertex_out_ID[e->GetEnd()];
			fprintf(outFile, "%d %d\n", startVertID, endVert1ID);
	}
	std::cout << "save starting info in " << fileName_ << std::endl;
	fclose(outFile);

	return true;
}

// save weaving geometry to .obj file
bool WeavingObject::OutputObj() {
	char* pt = strchr(fileName_, '.');
	if (CYLINDER == renderType_) 	{
		strcpy(pt, "__cylinder.obj");	// out file name: original_cylinder.obj
	} 
	else if (BEZIER_RIBBON == renderType_) {
		strcpy(pt, "__ribbon.obj");	// out file name: original_ribbon.obj
	} 
	else 	{
		// to do
	}
	std::cout << "weaving geometry is saved to " << fileName_ << std::endl;
	FILE * outFile = fopen(fileName_, "wb");
	if ( NULL == outFile ) {
		printf("can not open the outFile %s\n", fileName_);
		return false;
	}
	// vertex list
	long numQuads = 0;
	Weave(true, numQuads, outFile);
	fprintf(outFile, "# vertex\n");
	// face list
	std::cout << "num of quads : " << numQuads << std::endl;
	long temp = 0;
	for (long i = 0; i < numQuads; i++)
	{
		temp = i * 4;
		fprintf(outFile, "f %d %d %d %d\n", temp+1, temp+2, temp+3, temp+4);
	}
	fprintf(outFile, "# obj file writer by Qing Xing\n");
	fclose(outFile);

	return true;
}

void WeavingObject::DrawBaseMesh() {
	// draw vertices
// 	glDisable(GL_LIGHTING);
// 	glDisable(GL_TEXTURE_2D);
// 	glColor4f(0.0, 0.0, 1.0, 1.0);	// opaque blue
// 	std::set<Vertex*> vertices = GetVertices();
// 	glBegin(GL_POINTS);
// 	for ( std::set<Vertex*>::iterator vi = vertices.begin();
// 		vi != vertices.end(); vi++ ) {
// 		glVertex3fv( GetVertexCoordinates(*vi) );
// 	}
// 	glEnd();

	// draw edges
	glDisable(GL_LIGHTING);
	glDisable(GL_TEXTURE_2D);
	std::set<Edge*> edges = GetEdges();
	glBegin(GL_LINES);
	for ( std::set<Edge*>::iterator ei = edges.begin();
		ei != edges.end(); ei++ ) {
			Edge* e = *ei;
			switch ( edgePro[e].twist_ ) {
				case 1:
					glLineWidth(3.0);
					glColor4f(1.0, 0.0, 0.0, 1.0);	// opaque red
					break;
				case -1:
					glColor4f(0.0, 0.0, 1.0, 1.0);	// opaque blue
					break;;
				case 0:
					glLineWidth(3.0);
					glColor4f(0.0, 0.0, 0.0, 1.0);	// opaque black
					break;;
				case 2:
					glLineWidth(3.0);
					glColor4f(0.0, 1.0, 0.0, 1.0);	// opaque green
					break;
			}
			// yarn starting edges
// 			if (yarnStartEdges_set.find(e) != yarnStartEdges_set.end()) {
// 				glLineWidth(3.0);
// 				glColor4f(0.0, 0.0, 1.0, 1.0);	// opaque blue
// 			}
			Vertex* vert = e->GetStart();
			glVertex3fv( GetVertexCoordinates(vert) );
			vert = e->GetEnd();
			glVertex3fv( GetVertexCoordinates(vert) );
	}
	glEnd();

	// draw faces
// 	if ( isFaceDirty_ ) {
// 		ReComputeFaces();
// 		isFaceDirty_ = false;
// 	}
// 	glEnable(GL_LIGHTING);
// 	glDisable(GL_TEXTURE_2D);
// 	glColor4f(1.0, 0.5, 0.2, 0.3);	// transparent gray
// 	std::set<Face*> faces = GetFaces();
// 	for(std::set<Face*>::iterator fi = faces.begin();
// 		fi != faces.end(); 
// 		++fi) {
// 			glBegin(GL_POLYGON);
// 			Face* f = *fi;
// 			glNormal3fv( GetFaceNormal(f) );
// 			std::vector<Vertex*> vertices = f->GetVertices();
// 			for( std::vector<Vertex*>::iterator vi= vertices.begin();
// 				vi != vertices.end(); vi++) {
// 					glVertex3fv( GetVertexCoordinates(*vi) );
// 			}
// 			glEnd();
// 	}

	// draw weaving pieces
	if ( isFlaten ) {
		glDisable(GL_LIGHTING);
		glDisable(GL_TEXTURE_2D);
		glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
		for (int i = 0; i < numWeavingPiece ; i++)
		{
			// paper strip
			glColor4f(0.0, 1.0, 0.0, 1.0);	// opaque green
			glBegin(GL_QUAD_STRIP);
			for (int b = 0; b < 33; b++)
			{
				glVertex3fv( &( weavingPieces[i].Apoints[b].x ) );
				glVertex3fv( &( weavingPieces[i].Bpoints[b].x ) );
			}
			glEnd();

			// points
			glColor4f(1.0, 1.0, 1.0, 1.0);	// white
			cyPoint3f tempP;
			glBegin(GL_POINTS);
				// left hole - 4
				tempP = (weavingPieces[i].Apoints[4] + weavingPieces[i].Bpoints[4]) / 2.0;
				glVertex3fv( &(tempP.x) );
				// e1 - 8
				tempP = (weavingPieces[i].Apoints[8] + weavingPieces[i].Bpoints[8]) / 2.0;
				glVertex3fv( &(tempP.x) );
				// e2 - 12
				tempP = (weavingPieces[i].Apoints[12] + weavingPieces[i].Bpoints[12]) / 2.0;
				glVertex3fv( &(tempP.x) );
				// middle hole - 16
				tempP = (weavingPieces[i].Apoints[16] + weavingPieces[i].Bpoints[16]) / 2.0;
				glVertex3fv( &(tempP.x) );
				// e4 - 20
				tempP = (weavingPieces[i].Apoints[20] + weavingPieces[i].Bpoints[20]) / 2.0;
				glVertex3fv( &(tempP.x) );
				// e3 - 24
				tempP = (weavingPieces[i].Apoints[24] + weavingPieces[i].Bpoints[24]) / 2.0;
				glVertex3fv( &(tempP.x) );
				// right hole - 28
				tempP = (weavingPieces[i].Apoints[28] + weavingPieces[i].Bpoints[28]) / 2.0;
				glVertex3fv( &(tempP.x) );
			glEnd();
			// edge ID
		}
		glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
	}

	// xyz frame
	glDisable(GL_LIGHTING);
	glDisable(GL_TEXTURE_2D);
	glBegin(GL_LINES);
	// x - red
	glColor4f(1.0, 0.0, 0.0, 1.0);
	glVertex3i(0, 0, 0);
	glVertex3i(10, 0, 0);
	// y - green
	glColor4f(0.0, 1.0, 0.0, 1.0);
	glVertex3i(0, 0, 0);
	glVertex3i(0, 10, 0);
	// z - blue
	glColor4f(0.0, 0.0, 1.0, 1.0);
	glVertex3i(0, 0, 0);
	glVertex3i(0, 0, 10);
	glEnd();
}

void WeavingObject::DrawBezierPiece(const cyPoint3f& P0, const cyPoint3f& P1, 
									const cyPoint3f& P2, const cyPoint3f& P3,
									const cyPoint3f& PNormal, 
									const cyPoint3f& Q0, const cyPoint3f& Q1, 
									const cyPoint3f& Q2, const cyPoint3f& Q3, 
									const cyPoint3f& QNormal,
									bool isOutput, long& numQuads, FILE* file)
{
	//  A0---A1---A2---A3
	//	|    |    |    | ---> weaving yarn
	//	B0---B1---B2---B3
	std::vector<cyPoint3f> A, AN, B, BN;
	A.reserve(4);
	AN.reserve(4);
	B.reserve(4);
	BN.reserve(4);
	A.push_back(P0);	// A0
	A.push_back(P2);	// A1
	A.push_back(Q0);	// A2
	A.push_back(Q2);	// A3
	B.push_back(P1);	// B0
	B.push_back(P3);	// B1
	B.push_back(Q1);	// B2
	B.push_back(Q3);	// B3
	AN.push_back(PNormal);	// AN0
	// AN1 : area-weighted normal
	cyPoint3f temp0 = A[0] - A[1];	// P02
	cyPoint3f temp1 = B[1] - A[1];	// P32
	cyPoint3f tempN = temp0.Cross(temp1);
	temp0 = A[2] - A[1];	// Q0P2
	tempN += temp1.Cross(temp0);
	AN.push_back( tempN.GetNormalized() );
	// AN2
	temp0 = A[1] - A[2];	// P2Q0
	temp1 = B[2] - A[2];	// Q10
	tempN = temp0.Cross(temp1);
	temp0 = A[3] - A[2];	// Q20
	tempN += temp1.Cross(temp0);
	AN.push_back( tempN.GetNormalized() );
	AN.push_back(QNormal);	// AN3
	BN.push_back(PNormal);	// BN0
	// BN1
	temp0 = B[0] - B[1];	// P13
	temp1 = A[1] - B[1];	// P23
	tempN = temp1.Cross(temp0);
	temp0 = B[2] - B[1];	// Q1P3
	tempN += temp0.Cross(temp1);
	BN.push_back( tempN.GetNormalized() );
	// BN2
	temp0 = B[1] - B[2];	// P3Q1
	temp1 = A[2] - B[2];	// Q01
	tempN = temp1.Cross(temp0);
	temp0 = B[3] - B[2];	// Q31
	tempN += temp0.Cross(temp1);
	BN.push_back( tempN.GetNormalized() );
	BN.push_back(QNormal);	// BN3

	float t, it, tstep = 1.0 / float(numSamples_);
	float t0, t1, t2, t3;
	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
	glBegin(GL_QUADS);
	for ( int i = 0; i < numSamples_; i++ ) {
		t = float(i) * tstep;
		it = 1.0 - t;
		t0 = it * it * it;
		t1= 3.0 * t * it *it;
		t2 = 3.0* t * t * it;
		t3 = t * t * t;
		// vertex Ai
		temp0 = t0 * A[0] + t1 * A[1] + t2 * A[2] + t3 * A[3];
		tempN = t0 * AN[0] + t1 * AN[1] + t2 * AN[2] + t3 * AN[3];
		glNormal3fv( &(tempN.x) );
		glTexCoord2f(t, 0.98);
		glVertex3fv( &(temp0.x) );
		if (isOutput) {
			fprintf(file, "v %f %f %f\n", temp0.x, temp0.y, temp0.z);
		}
		if ( isFlaten ) {
			tempWeavingCycle[tempNumCrossings].Apoints[i] = temp0;
		}
		// vertex Bi
		temp0 = t0 * B[0] + t1 * B[1] + t2 * B[2] + t3 * B[3];
		tempN = t0 * BN[0] + t1 * BN[1] + t2 * BN[2] + t3 * BN[3];
		glNormal3fv( &(tempN.x) );
		glTexCoord2f(t, 0.02);
		glVertex3fv( &(temp0.x) );
		if (isOutput) {
			fprintf(file, "v %f %f %f\n", temp0.x, temp0.y, temp0.z);
		}
		if ( isFlaten ) {
			tempWeavingCycle[tempNumCrossings].Bpoints[i] = temp0;
		}
		t += tstep;
		it = 1.0 - t;
		t0 = it * it * it;
		t1= 3.0 * t * it *it;
		t2 = 3.0* t * t * it;
		t3 = t * t * t;
		// vertex Bi+1
		temp0 = t0 * B[0] + t1 * B[1] + t2 * B[2] + t3 * B[3];
		tempN = t0 * BN[0] + t1 * BN[1] + t2 * BN[2] + t3 * BN[3];
		glNormal3fv( &(tempN.x) );
		glTexCoord2f(t, 0.02);
		glVertex3fv( &(temp0.x) );
		if (isOutput) {
			fprintf(file, "v %f %f %f\n", temp0.x, temp0.y, temp0.z);
		}
		// vertex Ai+1
		temp0 = t0 * A[0] + t1 * A[1] + t2 * A[2] + t3 * A[3];
		tempN = t0 * AN[0] + t1 * AN[1] + t2 * AN[2] + t3 * AN[3];
		glNormal3fv( &(tempN.x) );
		glTexCoord2f(t, 0.98);
		glVertex3fv( &(temp0.x) );
		if (isOutput) {
			fprintf(file, "v %f %f %f\n", temp0.x, temp0.y, temp0.z);
			numQuads++;
		}
	}
	glEnd();
}

void WeavingObject::DrawControlPolygon(const cyPoint3f& P0, const cyPoint3f& P1, 
									   const cyPoint3f& P2, const cyPoint3f& P3,
									   const cyPoint3f& PNormal, 
									   const cyPoint3f& Q0, const cyPoint3f& Q1, 
									   const cyPoint3f& Q2, const cyPoint3f& Q3, 
									   const cyPoint3f& QNormal,
									   bool isOutput, long& numQuads, FILE* file)
{
	//  A0---A1---A2---A3
	//	|    |    |    | ---> weaving yarn
	//	B0---B1---B2---B3
	std::vector<cyPoint3f> A, AN, B, BN;
	A.reserve(4);
	AN.reserve(4);
	B.reserve(4);
	BN.reserve(4);
	A.push_back(P0);	// A0
	A.push_back(P2);	// A1
	A.push_back(Q0);	// A2
	A.push_back(Q2);	// A3
	B.push_back(P1);	// B0
	B.push_back(P3);	// B1
	B.push_back(Q1);	// B2
	B.push_back(Q3);	// B3
	AN.push_back(PNormal);	// AN0
	// AN1 : area-weighted normal
	cyPoint3f temp0 = A[0] - A[1];	// P02
	cyPoint3f temp1 = B[1] - A[1];	// P32
	cyPoint3f tempN = temp0.Cross(temp1);
	temp0 = A[2] - A[1];	// Q0P2
	tempN += temp1.Cross(temp0);
	AN.push_back( tempN.GetNormalized() );
	// AN2
	temp0 = A[1] - A[2];	// P2Q0
	temp1 = B[2] - A[2];	// Q10
	tempN = temp0.Cross(temp1);
	temp0 = A[3] - A[2];	// Q20
	tempN += temp1.Cross(temp0);
	AN.push_back( tempN.GetNormalized() );
	AN.push_back(QNormal);	// AN3
	BN.push_back(PNormal);	// BN0
	// BN1
	temp0 = B[0] - B[1];	// P13
	temp1 = A[1] - B[1];	// P23
	tempN = temp1.Cross(temp0);
	temp0 = B[2] - B[1];	// Q1P3
	tempN += temp0.Cross(temp1);
	BN.push_back( tempN.GetNormalized() );
	// BN2
	temp0 = B[1] - B[2];	// P3Q1
	temp1 = A[2] - B[2];	// Q01
	tempN = temp1.Cross(temp0);
	temp0 = B[3] - B[2];	// Q31
	tempN += temp0.Cross(temp1);
	BN.push_back( tempN.GetNormalized() );
	BN.push_back(QNormal);	// BN3

	glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
	glTexCoord2f(0.5, 0.5);		// texture is the same
	glBegin(GL_QUAD_STRIP);
	for ( int i = 0; i < 4; i++ ) {
		glNormal3fv( &(AN[i].x) );
		glVertex3fv( &(A[i].x) );

		glNormal3fv( &(BN[i].x) );
		glVertex3fv( &(B[i].x) );
	}
	glEnd();
}

void WeavingObject::DrawCylinder(const cyPoint3f& P0, const cyPoint3f& P1, 
									   const cyPoint3f& P2, const cyPoint3f& P3,
									   const cyPoint3f& PNormal, 
									   const cyPoint3f& Q0, const cyPoint3f& Q1, 
									   const cyPoint3f& Q2, const cyPoint3f& Q3, 
									   const cyPoint3f& QNormal,
									   bool isOutput, long& numQuads, FILE* file)
{
	//	A---B---C---D
	//	|     |      |     | ---> weaving yarn
	//	A---B---C---D
	cyPoint3f A[8], AN[8], B[8], BN[8], C[8], CN[8], D[8], DN[8];
	cyPoint3f center = 0.5 * (P0 + P1);
	A[0] = P0;
	AN[0] = P0 - center;
	float radius = AN[0].Length();
	AN[0].Normalize();
	A[4] = P1;
	AN[4] = -1 * AN[0];
	A[2] = center + radius * PNormal;
	AN[2] = PNormal;
	A[6] = center - radius * PNormal;
	AN[6] = -1 * PNormal;
	AN[1] = AN[0] + AN[2];		AN[1].Normalize();
	A[1] = center + radius * AN[1];
	AN[5] = -1 * AN[1];
	A[5] = center + radius * AN[5];
	AN[3] = AN[2] + AN[4];		AN[3].Normalize();
	A[3] = center + radius * AN[3];
	AN[7] = -1 * AN[3];
	A[7] = center + radius * AN[7];

	center = 0.5 * (P2 + P3);
	B[0] = P2;
	BN[0] = P2 - center;
	radius = BN[0].Length();
	BN[0].Normalize();
	B[4] = P3;
	BN[4] = -1 * BN[0];
	B[2] = center + radius * PNormal;
	BN[2] = PNormal;
	B[6] = center - radius * PNormal;
	BN[6] = -1 * PNormal;
	BN[1] = BN[0] + BN[2];		BN[1].Normalize();
	B[1] = center + radius * BN[1];
	BN[5] = -1 * BN[1];
	B[5] = center + radius * BN[5];
	BN[3] = BN[2] + BN[4];		BN[3].Normalize();
	B[3] = center + radius * BN[3];
	BN[7] = -1 * BN[3];
	B[7] = center + radius * BN[7];

	center = 0.5 * (Q0 +Q1);
	C[0] =Q0;
	CN[0] = Q0 - center;
	radius = CN[0].Length();
	CN[0].Normalize();
	C[4] = Q1;
	CN[4] = -1 * CN[0];
	C[2] = center + radius * QNormal;
	CN[2] = QNormal;
	C[6] = center - radius * QNormal;
	CN[6] = -1 * QNormal;
	CN[1] = CN[0] + CN[2];		CN[1].Normalize();
	C[1] = center + radius * CN[1];
	CN[5] = -1 * CN[1];
	C[5] = center + radius * CN[5];
	CN[3] = CN[2] + CN[4];		CN[3].Normalize();
	C[3] = center + radius * CN[3];
	CN[7] = -1 * CN[3];
	C[7] = center + radius * CN[7];

	center = 0.5 * (Q2 + Q3);
	D[0] = Q2;
	DN[0] = Q2 - center;
	radius = DN[0].Length();
	DN[0].Normalize();
	D[4] = Q3;
	DN[4] = -1 * DN[0];
	D[2] = center + radius * QNormal;
	DN[2] = QNormal;
	D[6] = center - radius * QNormal;
	DN[6] = -1 * QNormal;
	DN[1] = DN[0] + DN[2];		DN[1].Normalize();
	D[1] = center + radius * DN[1];
	DN[5] = -1 * DN[1];
	D[5] = center + radius * DN[5];
	DN[3] = DN[2] + DN[4];		DN[3].Normalize();
	D[3] = center + radius * DN[3];
	DN[7] = -1 * DN[3];
	D[7] = center + radius * DN[7];

	//glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
	glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
	glBegin(GL_QUADS);
	for ( int i = 7; i >= 1; i-- ) {
		// A---B
		glNormal3fv( &(AN[i].x) );
		glVertex3fv( &(A[i].x) );
		if (isOutput) {
			fprintf(file, "v %f %f %f\n", A[i].x, A[i].y, A[i].z);
		}
		glNormal3fv( &(BN[i].x) );
		glVertex3fv( &(B[i].x) );
		if (isOutput) {
			fprintf(file, "v %f %f %f\n", B[i].x, B[i].y, B[i].z);
		}
		glNormal3fv( &(BN[i-1].x) );
		glVertex3fv( &(B[i-1].x) );
		if (isOutput) {
			fprintf(file, "v %f %f %f\n", B[i-1].x, B[i-1].y, B[i-1].z);
		}
		glNormal3fv( &(AN[i-1].x) );
		glVertex3fv( &(A[i-1].x) );
		if (isOutput) {
			fprintf(file, "v %f %f %f\n", A[i-1].x, A[i-1].y, A[i-1].z);
			numQuads++;
		}
		// B---C
		glNormal3fv( &(BN[i].x) );
		glVertex3fv( &(B[i].x) );
		if (isOutput) {
			fprintf(file, "v %f %f %f\n", B[i].x, B[i].y, B[i].z);
		}
		glNormal3fv( &(CN[i].x) );
		glVertex3fv( &(C[i].x) );
		if (isOutput) {
			fprintf(file, "v %f %f %f\n", C[i].x, C[i].y, C[i].z);
		}
		glNormal3fv( &(CN[i-1].x) );
		glVertex3fv( &(C[i-1].x) );
		if (isOutput) {
			fprintf(file, "v %f %f %f\n", C[i-1].x, C[i-1].y, C[i-1].z);
		}
		glNormal3fv( &(BN[i-1].x) );
		glVertex3fv( &(B[i-1].x) );
		if (isOutput) {
			fprintf(file, "v %f %f %f\n", B[i-1].x, B[i-1].y, B[i-1].z);
			numQuads++;
		}
		// C---D
		glNormal3fv( &(CN[i].x) );
		glVertex3fv( &(C[i].x) );
		if (isOutput) {
			fprintf(file, "v %f %f %f\n", C[i].x, C[i].y, C[i].z);
		}
		glNormal3fv( &(DN[i].x) );
		glVertex3fv( &(D[i].x) );
		if (isOutput) {
			fprintf(file, "v %f %f %f\n", D[i].x, D[i].y, D[i].z);
		}
		glNormal3fv( &(DN[i-1].x) );
		glVertex3fv( &(D[i-1].x) );
		if (isOutput) {
			fprintf(file, "v %f %f %f\n", D[i-1].x, D[i-1].y, D[i-1].z);
		}
		glNormal3fv( &(CN[i-1].x) );
		glVertex3fv( &(C[i-1].x) );
		if (isOutput) {
			fprintf(file, "v %f %f %f\n", C[i-1].x, C[i-1].y, C[i-1].z);
			numQuads++;
		}
	}
	glNormal3fv( &(AN[0].x) );
	glVertex3fv( &(A[0].x) );
	if (isOutput) {
		fprintf(file, "v %f %f %f\n", A[0].x, A[0].y, A[0].z);
	}
	glNormal3fv( &(BN[0].x) );
	glVertex3fv( &(B[0].x) );
	if (isOutput) {
		fprintf(file, "v %f %f %f\n", B[0].x, B[0].y, B[0].z);
	}
	glNormal3fv( &(BN[7].x) );
	glVertex3fv( &(B[7].x) );
	if (isOutput) {
		fprintf(file, "v %f %f %f\n", B[7].x, B[7].y, B[7].z);
	}
	glNormal3fv( &(AN[7].x) );
	glVertex3fv( &(A[7].x) );
	if (isOutput) {
		fprintf(file, "v %f %f %f\n", A[7].x, A[7].y, A[7].z);
		numQuads++;
	}

	glNormal3fv( &(BN[0].x) );
	glVertex3fv( &(B[0].x) );
	if (isOutput) {
		fprintf(file, "v %f %f %f\n", B[0].x, B[0].y, B[0].z);
	}
	glNormal3fv( &(CN[0].x) );
	glVertex3fv( &(C[0].x) );
	if (isOutput) {
		fprintf(file, "v %f %f %f\n", C[0].x, C[0].y, C[0].z);
	}
	glNormal3fv( &(CN[7].x) );
	glVertex3fv( &(C[7].x) );
	if (isOutput) {
		fprintf(file, "v %f %f %f\n", C[7].x, C[7].y, C[7].z);
	}
	glNormal3fv( &(BN[7].x) );
	glVertex3fv( &(B[7].x) );
	if (isOutput) {
		fprintf(file, "v %f %f %f\n", B[7].x, B[7].y, B[7].z);
		numQuads++;
	}

	glNormal3fv( &(CN[0].x) );
	glVertex3fv( &(C[0].x) );
	if (isOutput) {
		fprintf(file, "v %f %f %f\n", C[0].x, C[0].y, C[0].z);
	}
	glNormal3fv( &(DN[0].x) );
	glVertex3fv( &(D[0].x) );
	if (isOutput) {
		fprintf(file, "v %f %f %f\n", D[0].x, D[0].y, D[0].z);
	}
	glNormal3fv( &(DN[7].x) );
	glVertex3fv( &(D[7].x) );
	if (isOutput) {
		fprintf(file, "v %f %f %f\n", D[7].x, D[7].y, D[7].z);
	}
	glNormal3fv( &(CN[7].x) );
	glVertex3fv( &(C[7].x) );
	if (isOutput) {
		fprintf(file, "v %f %f %f\n", C[7].x, C[7].y, C[7].z);
		numQuads++;
	}
	glEnd();
}

void WeavingObject::DrawB_SplinePiece(const cyPoint3f& P0, const cyPoint3f& P1, 
									   const cyPoint3f& P2, const cyPoint3f& P3,
									   const cyPoint3f& PNormal, 
									   const cyPoint3f& Q0, const cyPoint3f& Q1, 
									   const cyPoint3f& Q2, const cyPoint3f& Q3, 
									   const cyPoint3f& QNormal,
									   bool isOutput, long& numQuads, FILE* file)
{
										   //  A0---A1---A2---A3
										   //	|    |    |    | ---> weaving yarn
										   //	B0---B1---B2---B3
										   std::vector<cyPoint3f> A, AN, B, BN;
										   A.reserve(4);
										   AN.reserve(4);
										   B.reserve(4);
										   BN.reserve(4);
										   A.push_back(P0);	// A0
										   A.push_back(P2);	// A1
										   A.push_back(Q0);	// A2
										   A.push_back(Q2);	// A3
										   B.push_back(P1);	// B0
										   B.push_back(P3);	// B1
										   B.push_back(Q1);	// B2
										   B.push_back(Q3);	// B3
										   AN.push_back(PNormal);	// AN0
										   // AN1 : area-weighted normal
										   cyPoint3f temp0 = A[0] - A[1];	// P02
										   cyPoint3f temp1 = B[1] - A[1];	// P32
										   cyPoint3f tempN = temp0.Cross(temp1);
										   temp0 = A[2] - A[1];	// Q0P2
										   tempN += temp1.Cross(temp0);
										   AN.push_back( tempN.GetNormalized() );
										   // AN2
										   temp0 = A[1] - A[2];	// P2Q0
										   temp1 = B[2] - A[2];	// Q10
										   tempN = temp0.Cross(temp1);
										   temp0 = A[3] - A[2];	// Q20
										   tempN += temp1.Cross(temp0);
										   AN.push_back( tempN.GetNormalized() );
										   AN.push_back(QNormal);	// AN3
										   BN.push_back(PNormal);	// BN0
										   // BN1
										   temp0 = B[0] - B[1];	// P13
										   temp1 = A[1] - B[1];	// P23
										   tempN = temp1.Cross(temp0);
										   temp0 = B[2] - B[1];	// Q1P3
										   tempN += temp0.Cross(temp1);
										   BN.push_back( tempN.GetNormalized() );
										   // BN2
										   temp0 = B[1] - B[2];	// P3Q1
										   temp1 = A[2] - B[2];	// Q01
										   tempN = temp1.Cross(temp0);
										   temp0 = B[3] - B[2];	// Q31
										   tempN += temp0.Cross(temp1);
										   BN.push_back( tempN.GetNormalized() );
										   BN.push_back(QNormal);	// BN3

										   glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
										   glBegin(GL_QUAD_STRIP);
										   for ( int i = 0; i < 4; i++ ) {
											   glNormal3fv( &(AN[i].x) );
											   glVertex3fv( &(A[i].x) );

											   glNormal3fv( &(BN[i].x) );
											   glVertex3fv( &(B[i].x) );
										   }
										   glEnd();
}

// trace from edge e, starting vertex s, tracType
// return: true - open yarn
//		   false - closed yarn
// displaceType * edge_normal = displacement
bool WeavingObject::YarnTrace(Edge* e0, Vertex* s0, const int traceType0, 
							  std::map< Edge*, std::pair<Vertex*, int> >& visited_edges,
							  std::set<Edge*>& edges,
							  bool isOutput, long& numQuads, FILE* file) 
{
	Edge* e = e0;
	Vertex* s = s0;		// standing vertex
	Vertex* t = e0->GetOtherEnd(s0);
	int traceType = traceType0;		// trace type: flip between -1 / +1
	bool isOneHalfDone = false;
	cyPoint3f edgeNormal;
	std::vector<cyPoint3f> P;	P.reserve(5);		// Bezier piece control points
	std::vector<cyPoint3f> Q;	Q.reserve(6);
	// generate vector elements
	for ( int i = 0; i < 6; i++ ) {
		Q.push_back(edgeNormal);
	}
	// for flatten
	tempNumCrossings = 0;

	do {
		// ------ finish process this edge ------
		if ( visited_edges.find(e) == visited_edges.end() ) {
			// 1st time to visit this edge
			visited_edges[e] = std::make_pair(s, traceType);
		} 
		else {
			// 2nd time to visit this edge
			edges.erase(e);
			visited_edges.erase(e);
		}

		// ------ projection method ------
		ComputeControlPoints(e, s, traceType, Q, edgeNormal);
		if ( !P.empty() ) {
			// for flatten
			if ( isFlaten ) {
				tempWeavingCycle[tempNumCrossings].edgeID = e->GetID();
				tempWeavingCycle[tempNumCrossings].displaceType = traceType * edgePro[e].twist_;
				tempWeavingCycle[tempNumCrossings].edgeNormal = edgeNormal;
			}
			switch ( renderType_ ) {
				case BEZIER_RIBBON:
					DrawBezierPiece(P[0], P[1], P[2], P[3], P[4], Q[0], Q[1], Q[2], Q[3], edgeNormal, isOutput, numQuads, file);
					break;
				case CONTROL_POLYGON:
					DrawControlPolygon(P[0], P[1], P[2], P[3], P[4], Q[0], Q[1], Q[2], Q[3], edgeNormal, isOutput, numQuads, file);
					break;
				case CYLINDER:
					DrawCylinder(P[0], P[1], P[2], P[3], P[4], Q[0], Q[1], Q[2], Q[3], edgeNormal, isOutput, numQuads, file);
					break;
				case B_SPLINE_RIBBON:
					DrawB_SplinePiece(P[0], P[1], P[2], P[3], P[4], Q[0], Q[1], Q[2], Q[3], edgeNormal, isOutput, numQuads, file);
					break;
				default:
					std::cout << "render type is wrong" << std::endl;
			}
			// for flatten
			if ( isFlaten ) {
				tempNumCrossings++;
			}
		}
		else {
			// generate vector elements
			for ( int i = 0; i < 5; i++ ) {
				P.push_back(edgeNormal);
			}
		}
		// P[5] for next Bezier piece
		P[0] = Q[2];
		P[1] = Q[3];
		P[2] = Q[4];
		P[3] = Q[5];
		P[4] = edgeNormal;
		// ------ proceed to next edge ------
		if ( 1 == edgePro[e].twist_  || -1 == edgePro[e].twist_ ) {
			// 0_twist: same
			// 2_twist: same
			// +1_twist / -1_twist: flip sign
			traceType *= -1;
		}
		if ( edgePro[e].twist_ != 2 ) {
			s = t;			// switch standing vertex
		}			
		// NULL edge : skip until find the x_next edge
		e = s->Get_X_NextEdgeInRotation(e, traceType);
		
		if ( NULL == e ) {
			// open yarn end
			if ( isOneHalfDone ) {
				return true;	// open yarn : e0 is in the middle
			} 
			else {
				isOneHalfDone = true;
				// trace the other half open yarn from e0
				s = s0;
				e = e0;
				traceType =  -1 * traceType0;
				e = s->Get_X_NextEdgeInRotation(e, traceType);

				if ( NULL == e ) {
					return true;	// open yarn : e0 is the open end
				} 
				else {
					ComputeControlPoints(e0, s0, traceType0, Q, edgeNormal);
					// P[5] for next Bezier piece
					P[0] = Q[3];	// opposite yarn direction
					P[1] = Q[2];
					P[2] = Q[1];
					P[3] = Q[0];
					P[4] = edgeNormal;
					// ------ proceed to next edge ------
					t = e->GetOtherEnd(s);
				}					
			}
		} 
		else {
			t = e->GetOtherEnd(s);
		}
	} while( e != e0 || s != s0 || traceType != traceType0 );

	// closed yarn
	// draw the Bezier ribbon between last edge and e0
	ComputeControlPoints(e0, s0, traceType0, Q, edgeNormal);
	// for flatten
	if ( isFlaten ) {
		tempWeavingCycle[tempNumCrossings].edgeID = e0->GetID();
		tempWeavingCycle[tempNumCrossings].displaceType = traceType0 * edgePro[e0].twist_;
		tempWeavingCycle[tempNumCrossings].edgeNormal = edgeNormal;
	}
	switch ( renderType_ ) {
		case BEZIER_RIBBON:
			DrawBezierPiece(P[0], P[1], P[2], P[3], P[4], Q[0], Q[1], Q[2], Q[3], edgeNormal, isOutput, numQuads, file);
			break;
		case CONTROL_POLYGON:
			DrawControlPolygon(P[0], P[1], P[2], P[3], P[4], Q[0], Q[1], Q[2], Q[3], edgeNormal, isOutput, numQuads, file);
			break;
		case CYLINDER:
			DrawCylinder(P[0], P[1], P[2], P[3], P[4], Q[0], Q[1], Q[2], Q[3], edgeNormal, isOutput, numQuads, file);
			break;
		case B_SPLINE_RIBBON:
			DrawB_SplinePiece(P[0], P[1], P[2], P[3], P[4], Q[0], Q[1], Q[2], Q[3], edgeNormal, isOutput, numQuads, file);
			break;
		default:
			std::cout << "render type is wrong" << std::endl;
	}
	// for flatten
	if ( isFlaten ) {
		tempNumCrossings++;
	}

	return false;
}

bool WeavingObject::Weave(bool isOutput, long& numQuads, FILE* file) {
	// < edge, <start vertex, trace type> >
	std::map< Edge*, std::pair<Vertex*, int> > visited_edges;		// half-visited
	std::set<Edge*> edges;
 	for (std::set<Edge*>::iterator it = edges_.begin();
		it != edges_.end(); ++it) {
			edges.insert((*it));
	}
	bool isOpenYarn;
	int colorID;
	Edge* e;
	TraceCorner corner;

	if (NULL == yarnStartConters) {		// ------ 1st time weave
		if ( !edges.empty() ) {		// weaving object has non-NULL edges
			yarnStartConters = new std::vector<TraceCorner>;
			numClosedYarns_ = 0;
			numOpenYarns_ = 0;
			while ( !edges.empty() ) {
				// 			if (NULL != colorWheel_) 	{
				// 				colorID = (numOpenYarns_ + numClosedYarns_ + numColors_ / 2) % numColors_;
				// 				glColor3ubv(colorWheel_[colorID]);
				// 			}
				e = *(edges.begin());
				if ( visited_edges.find(e) == visited_edges.end() ) {
					// 1st time to visit the edge
					corner.e = e;
					corner.s = e->GetStart();
					corner.traceType = 1;
					isOpenYarn = YarnTrace(corner.e, corner.s, corner.traceType, visited_edges, edges, isOutput, numQuads, file);
				} 
				else {
					// one corner has been processed : switch standing vertex
					corner.e = e;
					corner.s = e->GetOtherEnd(visited_edges[e].first);		// switch standing vertex
					corner.traceType = visited_edges[e].second;
					if((1 == edgePro[e].twist_) ||  (-1 == edgePro[e].twist_)) {
						corner.traceType *= -1;
					}
					isOpenYarn = YarnTrace(corner.e, corner.s, corner.traceType, visited_edges, edges, isOutput, numQuads, file);
				}
				if ( isOpenYarn )
					++numOpenYarns_;
				else
					++numClosedYarns_;

				// add the yarn starting edge to the map
				yarnStartConters->push_back(corner);
				yarnStartEdges_set.insert(corner.e);
			}

			if ( !visited_edges.empty() || !edges.empty()) {
				std::cout << "error: edge untwisted left" << std::endl;
				return false;
			}

			std::cout << "number of open yarns: " << numOpenYarns_ << std::endl;
			std::cout << "number of closed yarns: " << numClosedYarns_ << std::endl;

			// compute colors
			if (numOpenYarns_ + numClosedYarns_ > 0 && NULL == colorWheel_) {
				ComputeColorWheel( numOpenYarns_ + numClosedYarns_, 10, 300);
				std::cout << "---------   recompute colors   ---------" << std::endl;
			}
		}
	} 
	else {		// ------ trace yarns in the preset order
		/*
		std::map<> use binary tree structure, so the order of elements is different from the order that elements are inserted
		*/
		int yarnID = 0;
		int groupID = 0;
		unsigned int r, g, b;
		numClosedYarns_ = 0;
		numOpenYarns_ = 0;
		float temp = pow(2.0, subdivisionLevel);		// number of yarn in one subdivision level
		for (std::vector<TraceCorner>::iterator ci = yarnStartConters->begin();
			ci != yarnStartConters->end(); ci++) 
		{
					colorID = (yarnID / (int)temp) % numColors_;
					groupID = yarnID % (int)temp;
					glColor3ubv(colorWheel_[colorID]);
// 					r = unsigned int(colorWheel_[colorID][0]) + groupID * 30;
// 					if (r > 255) {
// 						r = 255;
// 					}
// 					g = unsigned int(colorWheel_[colorID][1]) + groupID * 30;
// 					if (g > 255) {
// 						g = 255;
// 					}
// 					b = unsigned int(colorWheel_[colorID][2]) + groupID * 30;
// 					if (b > 255) {
// 						b = 255;
// 					}
// 					glColor3ub(unsigned char(r), unsigned char(g), unsigned char(b));
					isOpenYarn = YarnTrace((*ci).e, (*ci).s, (*ci).traceType, visited_edges, edges, isOutput, numQuads, file);
					++yarnID;
					if ( isOpenYarn )
						++numOpenYarns_;
					else
						++numClosedYarns_;
					// for flatten
					if ( isFlaten ) {
						int fisrtDownCrossing = 0;
						std::pair<int, int> tempPair;
						for (int i = 0; i < tempNumCrossings ; i++)
						{
							if ( -1 == tempWeavingCycle[i].displaceType ) {
								fisrtDownCrossing = i;
								break;
							}
						}
						for (int j = 0; j < tempNumCrossings ; j += 2)
						{
							int index = (fisrtDownCrossing + j) % tempNumCrossings;
							// e1
							weavingPieces[numWeavingPiece].e1 = tempWeavingCycle[index].edgeID;
							for (int k = 0; k < 4 ; k++)
							{
								weavingPieces[numWeavingPiece].Apoints[k] = tempWeavingCycle[index].Apoints[8 + k];
								weavingPieces[numWeavingPiece].Bpoints[k] = tempWeavingCycle[index].Bpoints[8 + k];
							}
							// e2
							index = (index + 1) % tempNumCrossings;
							tempPair.first = tempWeavingCycle[index].edgeID;;
							weavingPieces[numWeavingPiece].e2 = tempPair.first;
							for (int k = 0; k < numSamples_ ; k++)
							{
								weavingPieces[numWeavingPiece].Apoints[4 + k] = tempWeavingCycle[index].Apoints[k];
								weavingPieces[numWeavingPiece].Bpoints[4 + k] = tempWeavingCycle[index].Bpoints[k];
							}
							weavingPieces[numWeavingPiece].upEdgeNormal = tempWeavingCycle[index].edgeNormal;
							// e3
							index = (index + 1) % tempNumCrossings;
							tempPair.second = tempWeavingCycle[index].edgeID;
							weavingPieces[numWeavingPiece].e3 = tempPair.second;
							for (int k = 0; k < numSamples_ ; k++)
							{
								weavingPieces[numWeavingPiece].Apoints[16 + k] = tempWeavingCycle[index].Apoints[k];
								weavingPieces[numWeavingPiece].Bpoints[16 + k] = tempWeavingCycle[index].Bpoints[k];
							}
							weavingPieces[numWeavingPiece].e4 = edgeCorner[tempPair];
							// after e3
							index = (index + 1) % tempNumCrossings;
							// 32 segments; 33 points
							for (int k = 0; k < 5 ; k++)
							{
								weavingPieces[numWeavingPiece].Apoints[28 + k] = tempWeavingCycle[index].Apoints[k];
								weavingPieces[numWeavingPiece].Bpoints[28 + k] = tempWeavingCycle[index].Bpoints[k];
							}

							numWeavingPiece++;
						}
					}
		}
		std::cout << "number of open yarns: " << numOpenYarns_ << std::endl;
		std::cout << "number of closed yarns: " << numClosedYarns_ << std::endl;
		std::cout << "number of yarns: " << yarnID << std::endl;

		if ( !visited_edges.empty() || !edges.empty()) {
			std::cout << "error: edge untwisted left" << std::endl;
			return false;
		}
	}

	return true;
}

bool WeavingObject::Flatten()
{
	// construct mapping from (e2, e3) -> e4
	Face* tempFace = NULL;
	std::set<Face*> faces = GetFaces();
	int e2, e3, e4;
	for(std::set<Face*>::iterator fi = faces.begin();
		fi != faces.end(); ++fi)
	{
			tempFace = *fi;
			std::vector<Edge*> faceEdges = tempFace->GetEdges();
			int numFaceEdges = faceEdges.size();
			for (int i = 0; i < numFaceEdges ; i++)
			{
				e2 = faceEdges[i]->GetID();
				e3 = faceEdges[(i+1)%numFaceEdges]->GetID();
				e4 = faceEdges[(i+2)%numFaceEdges]->GetID();
				edgeCorner[std::make_pair(e3, e2)] = e4;	// face order vs crossing order
			}
	}

	// save the weaving piece geometry
	isFlaten = true;
	tempNumCrossings = 2 * GetEdges().size();
	tempWeavingCycle = new EdgeBeizerPiece[tempNumCrossings];
	weavingPieces = new WeavingPiece[tempNumCrossings];
	long temp = 0;
	FILE* tempFile = NULL;
	Weave(false, temp, tempFile);
	delete [] tempWeavingCycle;
	edgeCorner.clear();

	// flatten in space
	cyPoint3f flatApoints[33];
	cyPoint3f flatBpoints[33];
	cyPoint3f X, Y, Z, P;
	float a, b;
	for (int i = 0; i < numWeavingPiece ; i++)
	{
		Z = weavingPieces[i].upEdgeNormal;
		Z.Normalize();
		flatApoints[16] = weavingPieces[i].Apoints[16];
		flatBpoints[16] = weavingPieces[i].Bpoints[16];
		// 16 ~ 33
		for (int k = 16; k <= 31 ; k++)
		{
			// Apoints
			Y = weavingPieces[i].Apoints[k] - weavingPieces[i].Bpoints[k];
			Y.Normalize();
			a = Y.Dot(weavingPieces[i].Apoints[k + 1] - weavingPieces[i].Apoints[k]);
			P = weavingPieces[i].Apoints[k] + a * Y;
			b = (weavingPieces[i].Apoints[k + 1] - P).Length();
			Y = flatApoints[k] - flatBpoints[k];
			Y.Normalize();
			X = Y.Cross(Z);
			X.Normalize();
			flatApoints[k + 1] = flatApoints[k] + a * Y + b * X;
			// Bpoints
			Y = weavingPieces[i].Apoints[k + 1] - weavingPieces[i].Bpoints[k];
			Y.Normalize();
			a = Y.Dot(weavingPieces[i].Bpoints[k + 1] - weavingPieces[i].Bpoints[k]);
			P = weavingPieces[i].Bpoints[k] + a * Y;
			b = (weavingPieces[i].Bpoints[k + 1] - P).Length();
			Y = flatApoints[k + 1] - flatBpoints[k];
			Y.Normalize();
			X = Y.Cross(Z);
			X.Normalize();
			flatBpoints[k + 1] = flatBpoints[k] + a * Y + b * X;
		}

		// 16 ~ 0
		for (int k = 16; k > 0; k--)
		{
			// Apoints
			Y = weavingPieces[i].Apoints[k] - weavingPieces[i].Bpoints[k];
			Y.Normalize();
			a = Y.Dot(weavingPieces[i].Apoints[k - 1] - weavingPieces[i].Apoints[k]);
			P = weavingPieces[i].Apoints[k] + a * Y;
			b = (weavingPieces[i].Apoints[k - 1] - P).Length();
			Y = flatApoints[k] - flatBpoints[k];
			Y.Normalize();
			X = Y.Cross(Z);
			X.Normalize();
			flatApoints[k - 1] = flatApoints[k] + a * Y - b * X;
			// Bpoints
			Y = weavingPieces[i].Apoints[k - 1] - weavingPieces[i].Bpoints[k];
			Y.Normalize();
			a = Y.Dot(weavingPieces[i].Bpoints[k - 1] - weavingPieces[i].Bpoints[k]);
			P = weavingPieces[i].Bpoints[k] + a * Y;
			b = (weavingPieces[i].Bpoints[k - 1] - P).Length();
			Y = flatApoints[k - 1] - flatBpoints[k];
			Y.Normalize();
			X = Y.Cross(Z);
			X.Normalize();
			flatBpoints[k - 1] = flatBpoints[k] + a * Y - b * X;
		}

		// update
		for (int k = 0; k < 33 ; k++)
		{
			weavingPieces[i].Apoints[k] = flatApoints[k];
			weavingPieces[i].Bpoints[k] = flatBpoints[k];
		}
	}

	return true;
}


// place weaving piece in x-y plane
void WeavingObject::Placing2D()
{
	cyPoint3f C;	// center of weaving piece
	cyPoint3f X, Y, Z;		// unit vectors
	cyMatrix3f rotationMatrix;

	for (int i = 0; i < numWeavingPiece ; i++)
	{
		C = (weavingPieces[i].Apoints[16] + weavingPieces[i].Bpoints[16]) / 2.0;
		// construction rotation matrix
		Z = weavingPieces[i].upEdgeNormal;
		Y = weavingPieces[i].Apoints[16] - weavingPieces[i].Bpoints[16];
		Y.Normalize();
		X = Y.Cross(Z);
		X.Normalize();
		rotationMatrix.Set(X, Y, Z);
		rotationMatrix.Transpose();

		for (int k = 0; k < 33 ; k++)
		{
			// step 1: translate center of face to origin
			weavingPieces[i].Apoints[k] -= C;
			weavingPieces[i].Bpoints[k] -= C;
			// step 2: rotation 
			weavingPieces[i].Apoints[k] = rotationMatrix * weavingPieces[i].Apoints[k] ;
			weavingPieces[i].Bpoints[k] = rotationMatrix * weavingPieces[i].Bpoints[k];
		}
	}
}


//	find bounding box for each paper strip
//  get the max size
float WeavingObject::GetMaxBoundingBoxSize(WeavingPiece &wp)
{
	float xMin = wp.Apoints[0].x;
	float xMax = xMin;
	float yMax = wp.Apoints[0].y;
	float yMin = yMax;
	float xSize = 0, ySize = 0;
	cyPoint3f tempP;
	for (int i = 0; i < 33 ; i++)
	{
		tempP = wp.Apoints[i];
		if ( xMin > tempP.x) {
			xMin = tempP.x;
		}
		if ( xMax < tempP.x) {
			xMax = tempP.x;
		}
		if ( yMin > tempP.y) {
			yMin = tempP.y;
		}
		if ( yMax < tempP.y) {
			yMax = tempP.y;
		}
		tempP = wp.Bpoints[i];
		if ( xMin > tempP.x) {
			xMin = tempP.x;
		}
		if ( xMax < tempP.x) {
			xMax = tempP.x;
		}
		if ( yMin > tempP.y) {
			yMin = tempP.y;
		}
		if ( yMax < tempP.y) {
			yMax = tempP.y;
		}
	}

	xSize = xMax - xMin;
	ySize = yMax - yMin;
	// update
	if ( xSize > ySize ) {
		return xSize;
	} 
	else {
		return ySize;
	}
}


// bounding box is boundingBox_X*boundingBox_X
void WeavingObject::ScalePaperStrips()
{
	// get the maxSize for all paper strips
	float maxSize = 0;
	float temp = 0;
	for (int i = 0; i < numWeavingPiece ; i++)
	{
		temp = GetMaxBoundingBoxSize(weavingPieces[i]);
		if ( maxSize < temp ) {
			maxSize = temp;
		}
	}
	maxSize = boundingBox_size / maxSize;	// uniform scale factor	

	// scale all paper strips so that the max bounding box size is boundingBox_X*boundingBox_X
	for (int i = 0; i < numWeavingPiece ; i++)
	{
		for (int k = 0; k < 33 ; k++)
		{
			weavingPieces[i].Apoints[k] *= maxSize;
			weavingPieces[i].Bpoints[k] *= maxSize;
		}
	}
}


// save paper strip to .eps file
bool WeavingObject::OutputVertexPaperStrips(){
	ScalePaperStrips();

	Output out;
	out.width = paperWidth;
	out.height = paperHeight;
	out.scale = boundingBox_size;
	out.holeRadius = holeRadius;
	bool result = out.OutputPS(fileName_, weavingPieces, numWeavingPiece);

	return result;

}