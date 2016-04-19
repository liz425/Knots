#include "WeavingObject.h"

#include "Fenghui_Zhang_core/topology_object.h"
#include "Fenghui_Zhang_core/face.h"
#include "Fenghui_Zhang_core/edge.h"
#include "Fenghui_Zhang_core/vertex.h"
#include "Fenghui_Zhang_core/logging.h"
#include "Fenghui_Zhang_core/object_store.h"

#include "cyMatrix3.h"
#include "cyPoint.h"
#include "ColorSpace.h"
#include "Output.h"

#include <fl/gl.h>
#include <fl/glu.h>

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


Output::Output()
{
	width = 16.5;
	height = 22;
	scale = 5;
	holeRadius = 0.05;
}

char* Output::syCreateNewPage(char* fileName, int pagenum){
	std::string stringname;

	char* num;
	char mynum[4]; /* should be long enough */

	sprintf(mynum, "%d", pagenum);
	std::string path = fileName;
 	size_t found = path.find_last_of("/\\"); 
// 	path.erase(0, found);   // erase the path

	found = path.find_last_of(".");
	path.erase(found, path.size());   // erase the initial file type
//	stringname = "./results";
	stringname += path;
	stringname += mynum;
	stringname += ".eps";

	char* charname = new char[128];
	strcpy (charname, stringname.c_str());
	outFile = fopen(charname, "wb");    

	//outFile = fopen(stringname.c_str(), "wb");    

	//    return fileName;
	return charname;
}

int Output::syDrawBoundray(){

	float zero = 0;
	fprintf(outFile, "0 1 0 setrgbcolor\n");
	fprintf(outFile, "%f in %f in moveto\n", zero, zero);
	fprintf(outFile, "%f in %f in lineto\n", zero, zero);
	fprintf(outFile, "%f in %f in lineto\n", height, zero);
	fprintf(outFile, "stroke\n\n");

	fprintf(outFile, "0 1 0 setrgbcolor\n");
	fprintf(outFile, "%f in %f in moveto\n", height, zero);
	fprintf(outFile, "%f in %f in lineto\n", height, zero);
	fprintf(outFile, "%f in %f in lineto\n", height, width);
	fprintf(outFile, "stroke\n\n");

	fprintf(outFile, "0 1 0 setrgbcolor\n");
	fprintf(outFile, "%f in %f in moveto\n", height, width);
	fprintf(outFile, "%f in %f in lineto\n", height, width);
	fprintf(outFile, "%f in %f in lineto\n", zero, width);
	fprintf(outFile, "stroke\n\n");

	fprintf(outFile, "0 1 0 setrgbcolor\n");
	fprintf(outFile, "%f in %f in moveto\n", zero, zero);
	fprintf(outFile, "%f in %f in lineto\n", zero, zero);
	fprintf(outFile, "%f in %f in lineto\n", zero, width);
	fprintf(outFile, "stroke\n\n");

	return 0;
}


int Output::syLinetoPostScript(float sx, float sy, float tx, float ty){
	fprintf(outFile, "0 0 1 setrgbcolor\n");
	fprintf(outFile, "%f in %f in moveto\n", sx, sy);
	fprintf(outFile, "%f in %f in lineto\n", sx, sy);
	fprintf(outFile, "%f in %f in lineto\n", tx, ty);
	fprintf(outFile, "stroke\n\n");

	return 0;
}

int Output::syCircletoPostScript(float cx, float cy, float radius){
	fprintf(outFile, "0 0 1 setrgbcolor\n");
	fprintf(outFile, "%f in %f in moveto\n", cx, cy);
	fprintf(outFile, "%f in %f in %f in 0 360 arc closepath\n", cx, cy, radius);
	fprintf(outFile, "stroke\n\n");

	return 0;
}

int Output::syTexttoPostScript(int num, float cx, float cy){
	fprintf(outFile, "1 1 0 setrgbcolor\n");	// yellow
	fprintf(outFile, "%f in %f in moveto\n", cx, cy);
	fprintf(outFile, "(%d) show\n", num);
	fprintf(outFile, "stroke\n\n");

	return 0;
}

int Output::syTexttoPostScript_corner(int num, float cx, float cy){
	fprintf(outFile, "1 1 0 setrgbcolor\n");	// yellow
	fprintf(outFile, "%f in %f in moveto\n", cx, cy);
	fprintf(outFile, "(___%d) show\n", num);
	fprintf(outFile, "stroke\n\n");

	return 0;
}


bool Output::syOutputOneWeavingPiece(WeavingPiece &wp)
{
	// Ai---Ai+1	Bi---Bi+1
	for (int i = 0; i <= 31 ; i++)
	{
		syLinetoPostScript(wp.Apoints[i].x, wp.Apoints[i].y, wp.Apoints[i + 1].x, wp.Apoints[i + 1].y);
		syLinetoPostScript(wp.Bpoints[i].x, wp.Bpoints[i].y, wp.Bpoints[i + 1].x, wp.Bpoints[i + 1].y);
	}
	// ends: A0---B0	A32---B32
	syLinetoPostScript(wp.Apoints[0].x, wp.Apoints[0].y, wp.Bpoints[0].x, wp.Bpoints[0].y);
	syLinetoPostScript(wp.Apoints[32].x, wp.Apoints[32].y, wp.Bpoints[32].x, wp.Bpoints[32].y);

//	// test: a square
// 	syLinetoPostScript(0, 0, 10, 0 );
// 	syLinetoPostScript(10, 0, 10, 10 );
// 	syLinetoPostScript(10, 10, 0, 10 );
// 	syLinetoPostScript(0, 10, 0, 0 );


	// holes
	// left hole - 4
	cyPoint3f tempP = (wp.Apoints[4] + wp.Bpoints[4]) / 2.0;
	syCircletoPostScript(tempP.x, tempP.y, holeRadius);
	// middle hole - 16
	tempP = (wp.Apoints[16] + wp.Bpoints[16]) / 2.0;
	syCircletoPostScript(tempP.x, tempP.y, holeRadius);
	// right hole - 28
	tempP = (wp.Apoints[28] + wp.Bpoints[28]) / 2.0;
	syCircletoPostScript(tempP.x, tempP.y, holeRadius);

	// edge ID
	// e1 - 8
	tempP = (wp.Apoints[8] + wp.Bpoints[8]) / 2.0;
	syTexttoPostScript(wp.e1, tempP.x, tempP.y);
	// e2 - 12
	tempP = (wp.Apoints[12] + wp.Bpoints[12]) / 2.0;
	syTexttoPostScript(wp.e2, tempP.x, tempP.y);
	// e3 - 24
	tempP = (wp.Apoints[24] + wp.Bpoints[24]) / 2.0;
	syTexttoPostScript(wp.e3, tempP.x, tempP.y);
	// e4 - 20
	tempP = (wp.Apoints[20] + wp.Bpoints[20]) / 2.0;
	syTexttoPostScript_corner(wp.e4, tempP.x, tempP.y);

	return true;
}


void Output::sySetOrigin(WeavingPiece &wp, cyPoint3f origin)
{
// 	origin.x -= 0.5;
// 	origin.y -= 0.5;
	for (int i = 0; i < 33; i++){
		wp.Apoints[i].x -= origin.x;
		wp.Apoints[i].y -= origin.y;
		wp.Bpoints[i].x -= origin.x;
		wp.Bpoints[i].y -= origin.y;
	}
}

cyPoint3f Output::syGetMinimumPoint(WeavingPiece &wp)
{
	float minx = 10000.0;
	float miny = 10000.0;

	for (int i = 0; i < 33 ; i++)
	{
		if ( wp.Apoints[i].x < minx) {
			minx = wp.Apoints[i].x;
		}
		if ( wp.Apoints[i].y < miny) {
			miny = wp.Apoints[i].y;
		}

		if ( wp.Bpoints[i].x < minx) {
			minx = wp.Bpoints[i].x;
		}
		if ( wp.Bpoints[i].y < miny) {
			miny = wp.Bpoints[i].y;
		}
	}

	cyPoint3f minPoint(minx, miny, 0);
	return minPoint;
}


void Output::syArrageStrips(WeavingPiece &wp, int num)
{
	int heightNum = numofRows;
	int widthNum = numofColum;

	int x_shift_num = num % widthNum;
	double x_shift = x_shift_num * scale;

	int y_shift_num = num / widthNum;
	double y_shift = y_shift_num * scale;

	cyPoint3f origin(0, 0, 0);
	origin.x = origin.x - x_shift;
	origin.y = origin.y - y_shift;

	sySetOrigin(wp, origin);
}

bool Output::OutputPS(char* fileName, WeavingPiece* weavingPieces, long numWPs)
{
	numofRows = width / scale;
	numofColum = height / scale;
	int totalStripsinaPage = numofRows * numofColum;
	int pagecount = 1;

	for (int i = 0; i < numWPs ; i++)
	{
		if(i % totalStripsinaPage == 0){  // Open a new page
			char* outputName = syCreateNewPage(fileName, pagecount);
			if ( NULL == outFile ) {
				printf("can not open the outFile %s\n", outputName);
				return false;
			}
			printf("paper strip information is saved to %s\n", outputName);
			//	    std:: string outputName_str = std::string(outputName);
			//	    std::cout << "paper strip information is saved to " << outputName_str << std::endl;
			fprintf(outFile, "%%!PS-Adobe-3.0 EPSF-3.0 \n%%%%\BoundingBox: 0 0 4000 4000\n%%%%EndComments \n\n");
			fprintf(outFile, "/in {72 mul} def\n\n");
			fprintf(outFile, "/Times-Roman findfont \n 8 scalefont\n setfont\n newpath\n\n");

			syDrawBoundray();
			pagecount++;
		}

		tempWP = weavingPieces[i];
		cyPoint3f origin = syGetMinimumPoint(tempWP);	// Find the minimum x and y 
		sySetOrigin(tempWP, origin);  // Make the minimum x and y as origin point for this strip

		syArrageStrips(tempWP, i % totalStripsinaPage);  // Put the strip at a suitable place in the page

		syOutputOneWeavingPiece(tempWP); // Output the positions of this strip into file

		if(i % totalStripsinaPage == (totalStripsinaPage-1)  // Close this pace, if the page is full
			|| i == (numWPs-1))
		{
			fclose(outFile);
		}
	}

	return true;
}
