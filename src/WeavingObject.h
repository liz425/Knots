/********************************************************************
	created:	2009/09/30
	created:	30:9:2009   18:01
	filename: 	d:\Research\Knots_v2.0_Fenghui_core\WeavingObject.h
	file path:	d:\Research\Knots_v2.0_Fenghui_core
	file base:	WeavingObject
	file ext:	h
	author:		Qing Xing
	
	purpose:	
*********************************************************************/


#ifndef WEAVING_OBJECT__
#define WEAVING_OBJECT__

#include "Fenghui_Zhang_core/topology_object.h"

#include "cyPoint.h"

#include <map>

class Vertex;
class Edge;
class Face;

// TODO: when the vertex/face is being removed, we have to delete memory of pointers.
struct VertexProperty {
	float* coords_;
};

struct EdgeProperty {
	int twist_;		// twist number
						// NULL edge : 2
};

struct FaceProperty {
	float* normal_;
};

struct TraceCorner 
{
	Edge* e;
	Vertex* s;		// standing vertex
	int traceType;
};

struct EdgeBeizerPiece 
{
	int edgeID;
	int displaceType;
	cyPoint3f edgeNormal;
	cyPoint3f Apoints[12];
	cyPoint3f Bpoints[12];
};

struct WeavingPiece 
{
	int e1, e2, e3, e4;
	cyPoint3f upEdgeNormal;
	cyPoint3f Apoints[33];	// 32 segments; 33 points
	cyPoint3f Bpoints[33];
};

enum RenderType {BEZIER_RIBBON, CONTROL_POLYGON, CYLINDER, B_SPLINE_RIBBON};

class WeavingObject : public TopologyObject {
public:
	WeavingObject(void);
	~WeavingObject(void);

	bool Clear();		// clear all the vertices, edges, faces
	void init();
	bool ComputeColorWheel(int N, float minH, float maxH);		// compute N colors to use
	std::set<Face*> GetFaces();
	bool ReComputeFaces();
	void SetFaceDirtyFlag();
	void ClearFaceNormals();
	void SetCenterWeight(float w);
	void SetRenderType(RenderType rt);
	void SetWidth(const float w);
	void SetCurvature(const float c);
	void SetDiplaceFactor(const float d);
	void SetNumSamples(const int s);
	float* GetVertexCoordinates(Vertex* v);
	bool SetVertexCoordinates(Vertex* v, float* coords);
	bool ResetVertexCoordinates(Vertex* v, float x, float y, float z);
	float* GetFaceNormal(Face* f);
	cyPoint3f GetFaceCenter(Face* f, bool isNullEdgeWeighted);
	bool ComputeControlPoints(Edge* e, Vertex* s, const int traceType,
		std::vector<cyPoint3f>& Q, cyPoint3f& edgeNormal);		// project edge area
	bool LoadObjFile(const char* objFileName);
	void RandomTwist();		// randomly generate twist number for edge
	bool SaveFiles();			// save object to .obj file and .twist file
	bool OutputObj();			// save weaving geometry to .obj file
	void DrawBaseMesh();	// wire frame
	void DrawBezierPiece(const cyPoint3f& P0, const cyPoint3f& P1, const cyPoint3f& P2, const cyPoint3f& P3,
		const cyPoint3f& PNormal, const cyPoint3f& Q0, const cyPoint3f& Q1, const cyPoint3f& Q2,
		const cyPoint3f& Q3, const cyPoint3f& QNormal, bool isOutput, long& numQuads, FILE* file);
	void DrawControlPolygon(const cyPoint3f& P0, const cyPoint3f& P1, const cyPoint3f& P2, const cyPoint3f& P3,
		const cyPoint3f& PNormal, const cyPoint3f& Q0, const cyPoint3f& Q1, const cyPoint3f& Q2,
		const cyPoint3f& Q3, const cyPoint3f& QNormal, bool isOutput, long& numQuads, FILE* file);
	void DrawCylinder(const cyPoint3f& P0, const cyPoint3f& P1, const cyPoint3f& P2, const cyPoint3f& P3,
		const cyPoint3f& PNormal, const cyPoint3f& Q0, const cyPoint3f& Q1, const cyPoint3f& Q2,
		const cyPoint3f& Q3, const cyPoint3f& QNormal, bool isOutput, long& numQuads, FILE* file);
	void DrawB_SplinePiece(const cyPoint3f& P0, const cyPoint3f& P1, const cyPoint3f& P2, const cyPoint3f& P3,
		const cyPoint3f& PNormal, const cyPoint3f& Q0, const cyPoint3f& Q1, const cyPoint3f& Q2,
		const cyPoint3f& Q3, const cyPoint3f& QNormal, bool isOutput, long& numQuads, FILE* file);
	bool YarnTrace(Edge* e0, Vertex* s0, const int traceType0,
					std::map< Edge*, std::pair<Vertex*, int> >& visited_edges,
					std::set<Edge*>& edges, bool isOutput, long& numQuads, FILE* file);
	bool Weave(bool isOutput, long& numQuads, FILE* file);
	inline void ClearEdges() { edges_.clear(); }
	inline void InsertEdges(Edge* e) { edges_.insert(e); }
	bool Flatten();
	void Placing2D();	// place weaving piece in x-y plane
	float GetMaxBoundingBoxSize(WeavingPiece &wp);
	void ScalePaperStrips();	// bounding box is 1*1
	bool OutputVertexPaperStrips();			// save paper strip to .esp file

	std::map<Edge*, EdgeProperty> edgePro;
	std::set<Edge*> childrenEdges;	// children 1st 1-twist edge
	int subdivisionLevel;
	std::vector<TraceCorner>* yarnStartConters;	// pointer to array of yarn trace starting corners
	std::set<Edge*> yarnStartEdges_set;		// find the yarn starting edge quickly
	bool isFlaten;
	EdgeBeizerPiece* tempWeavingCycle;
	long tempNumCrossings;
	WeavingPiece* weavingPieces;	// array of weaving pieces
	long numWeavingPiece;
	std::map<std::pair<int, int>, int> edgeCorner;
	float boundingBox_size;
	double holeRadius;
	double paperWidth;
	double paperHeight;

private: 
	std::set<Face*> faces_;
	std::map<Vertex*, VertexProperty> vertPro_;
	//std::map<Edge*, EdgeProperty> edgePro_;
	std::map<Face*, FaceProperty> facePro_;
	// vertex_0 -> vertex_1 -> face_ID.
	std::map< Vertex*, std::map<Vertex*, Face*> > face_map_;
	bool isFaceDirty_;		// need to recompute face
	float weight4Center_;	// vertex weight to compute face center
	RenderType renderType_;
	int numOpenYarns_, numClosedYarns_;
	float width_, curvature_;		// Bezier piece control polygon
	float displaceFactor_;
	int numSamples_;		// number of samples in a Bezier piece
	unsigned char** colorWheel_;	// complementary colors / colors in a rang on the hue wheel
	int numColors_;
	int numFacesObjFile_;	// number of faces in the .obj file
	char fileName_[128];
};

#endif // WEAVING_OBJECT__
