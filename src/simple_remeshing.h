#ifndef SIMPLE_REMESHING__
#define SIMPLE_REMESHING__

class Vertex;
class Edge;
class Face;
class WeavingObject;
/**
 * Put a vertex at the center of each face and add an edge between this new
 * vertex and each vertex of the original face.
 */
void TriangulateFace(WeavingObject* objt, Face *f);

/**
 * Apply SplitFace to all faces.
 */
void SimpleRemeshing(WeavingObject* obj);

/**
 * Apply DooSabin subdivision on obj.
 */
void DooSabin(WeavingObject* obj);

/**
 * Apply Catmull_Clark subdivision on obj.
 */
void Catmull_Clark(WeavingObject* obj);
void LinearSub(WeavingObject* obj);
void Average(WeavingObject* obj);

struct EdgeInfo 
{
	Vertex* s;
	Vertex* t;
	Vertex* mid;
	Vertex* left;
	Vertex* right;
	Edge* e1;
	Edge* e2;
	Edge* e3;
	Edge* e4;
	int twist;		// edge twist number

	EdgeInfo(Vertex* start, Vertex* end, Vertex* midPoint, Vertex* l, Vertex* r, Edge* c1, Edge* c2, Edge* c3, Edge* c4, int edgeTwist) {
		s = start;
		t = end;
		mid = midPoint;
		left = l;
		right = r;
		e1 = c1;
		e2 = c2;
		e3 = c3;
		e4 = c4;
		twist = edgeTwist;
	}
};


#endif // SIMPLE_REMESHING__
