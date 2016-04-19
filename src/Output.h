#include "cyPoint.h"
#include "cstdio"
#include "Fenghui_Zhang_core/topology_object.h"
#include "Fenghui_Zhang_core/face.h"
#include "Fenghui_Zhang_core/edge.h"
#include "Fenghui_Zhang_core/vertex.h"
#include "Fenghui_Zhang_core/logging.h"
#include "Fenghui_Zhang_core/object_store.h"

#include "cyMatrix3.h"
#include "cyPoint.h"
#include "ColorSpace.h"
#include "WeavingObject.h"

#include <cstdlib>
#include <cstdio>

#include <ctime>
#include <cmath>
#include <string>
#include <cstring>


#include <iostream>
#include <map>
#include <vector>
#include <utility>
#include <algorithm>


class Output {

public:
	Output();

	char* syCreateNewPage(char* fileName, int pagenum);
	int syLinetoPostScript(float sx, float sy, float tx, float ty);
	int syCircletoPostScript(float cx, float cy, float radius);
	int syTexttoPostScript(int num, float cx, float cy);
	int syTexttoPostScript_corner(int num, float cx, float cy);
	bool syOutputOneWeavingPiece(WeavingPiece &wp);
	void sySetOrigin(WeavingPiece &wp, cyPoint3f origin);
	cyPoint3f syGetMinimumPoint(WeavingPiece &wp);
	void syArrageStrips(WeavingPiece &wp, int num);
	bool OutputPS(char* fileName, WeavingPiece* weavingPieces, long numWPs);
	int syDrawBoundray();

	double width;
	double height;
	double scale;
	double holeRadius;

	WeavingPiece tempWP;
	int numofRows;
	int numofColum;

	FILE * outFile;
};
