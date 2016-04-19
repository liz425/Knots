#include "Gl_Window.h"

#include "glErrorUtil.h"
#include "loadImage.h"
#include "Text_File.h"
#include "simple_remeshing.h"

#include <fl/gl.h>
#include <fl/glu.h>

#include <cstdio>

#include <iostream>
#include <map>
#include <vector>

bool Gl_Window::init() 
{
	/* ----- Checking for OpenGL 2 --------------- */
	const GLubyte* string;
	string = glGetString(GL_VENDOR);
	printf("Vendor: %s\n", string);
	string = glGetString(GL_RENDERER);
	printf("Renderer: %s\n", string);
	string = glGetString(GL_VERSION);
	printf("OpenGL Version: %s\n", string);
	string = glGetString(GL_SHADING_LANGUAGE_VERSION);
	printf("GLSL Version: %s\n", string);
	// GLSL initialization
	glewInit();		// An OpenGL context must already exist
	if (glewIsSupported("GL_VERSION_2_0")) {
		printf("Ready for OpenGL 2.0\n");
		//setShaders();	// GLSL set up
	}
	else {
		printf("OpenGL 2.0 not supported\n");
		exit(1);
	}

	GLfloat values[2];
	glGetFloatv(GL_POINT_SIZE_GRANULARITY, values);
	printf("GL_ALIASED_POINT_SIZE_GRANULARITY value is %3.1f\n", values[0]);
	glGetFloatv(GL_ALIASED_POINT_SIZE_RANGE, values);
	printf("GL_ALIASED_POINT_SIZE_RANGE values are %3.1f %3.1f\n", values[0], values[1]);
	glGetFloatv(GL_LINE_WIDTH_GRANULARITY, values);
	printf("GL_ALIASED_LINE_WIDTH_GRANULARITY value is %3.1f\n", values[0]);
	glGetFloatv(GL_ALIASED_LINE_WIDTH_RANGE, values);
	printf("GL_ALIASED_LINE_WIDTH_RANGE values are %3.1f %3.1f\n", values[0], values[1]);

	// OpenGL state machine settings
	glClearColor(1.0, 1.0, 1.0, 1.0);
	glClearDepth(1.0);
	glEnable(GL_DEPTH_TEST);
	// no culling for transparent faces
	//glEnable(GL_CULL_FACE);	// cull back faces
	glDisable(GL_CULL_FACE);
	glCullFace(GL_BACK);
	glShadeModel(GL_SMOOTH);	//Enable smooth shading
	glColorMaterial(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE);
	glEnable(GL_COLOR_MATERIAL);	// color material
	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	glHint(GL_LINE_SMOOTH_HINT, GL_DONT_CARE);
	glEnable(GL_LINE_SMOOTH);	// aliased lines
	glLineWidth(3.0);
	glEnable(GL_POINT_SMOOTH);	// round points
	glPointSize(9.0);

	// lighting setup
	glLightModelf(GL_LIGHT_MODEL_LOCAL_VIEWER, 1.0);
	glLightModelf(GL_LIGHT_MODEL_TWO_SIDE, 0.0);
	//glEnable(GL_NORMALIZE);
	//global ambient light
	GLfloat ambientColor[] = {0.2f, 0.2f, 0.2f, 0.0f};
	glLightModelfv(GL_LIGHT_MODEL_AMBIENT, ambientColor);
	// Lighting parameters:
	GLfloat light_pos[] = {2.0, 4.0, 2.0, 1.0};
	GLfloat light_Ka[]  = {0.2, 0.2, 0.2, 1.0};
	GLfloat light_Kd[]  = {0.6, 0.6, 0.6, 1.0};
	GLfloat light_Ks[]  = {1.0, 1.0, 1.0, 1.0};
	glLightfv(GL_LIGHT0, GL_POSITION, light_pos);
	glLightfv(GL_LIGHT0, GL_AMBIENT, light_Ka);
	glLightfv(GL_LIGHT0, GL_DIFFUSE, light_Kd);
	glLightfv(GL_LIGHT0, GL_SPECULAR, light_Ks);
	glEnable(GL_LIGHT0);
	glEnable(GL_LIGHTING);

	// Material parameters:
	//GLfloat material_Ka[] = {0.6, 0.0, 0.0, 1.0};
	//GLfloat material_Kd[] = {0.6, 0.4, 0.4, 1.0};
	GLfloat material_Ke[] = {0.0, 0.0, 0.0, 1.0};
	GLfloat material_Ks[] = {0.3, 0.3, 0.3, 1.0};
	GLfloat material_Se	  = 20.0;
	//glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, material_Ka);
	//glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, material_Kd);
	glMaterialfv(GL_FRONT_AND_BACK, GL_EMISSION, material_Ke);
	glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, material_Ks);
	glMaterialf (GL_FRONT_AND_BACK, GL_SHININESS, material_Se);

	// default textures
	setTexture("D:/Images/lines/lines0a.jpg", 0);	// environment texture - 0
// 	setTexture("D:/Images/env/white.png", 0);	// environment texture - 0
// 	setTexture("D:/Images/jpg/lines.bmp", 1);	// object texture - 1

	// generate display list after the OpenGL context is created

	// BUG FIXED by Z.L.: according to https://www.opengl.org/discussion_boards/showthread.php/169001-glCallList-not-working 
	// In order to avoid list definition before OpenGL is fully initialized, 
	// we should move glGenLists(1) into the draw function, right before glNewList(), 
	//displayList_Weaving = glGenLists(1);
	//displayList_BaseMesh = glGenLists(1);

	mesh->init();

	return true;
}

void Gl_Window::setShaders()
{
	char* vs = NULL;
	char* fs = NULL;

	vs = textFileRead("Texture_Phong.vert");
	fs = textFileRead("Texture_Phong.frag");

	const char* pv = vs;
	const char* pf = fs;

	vertShader_ = glCreateShader(GL_VERTEX_SHADER);
	fragShader_ = glCreateShader(GL_FRAGMENT_SHADER);

	glShaderSource(vertShader_, 1, &pv, NULL);
	glShaderSource(fragShader_, 1, &pf, NULL);

	glCompileShader(vertShader_);
	glCompileShader(fragShader_);

	// check errors
	printShaderInfoLog(vertShader_);
	printShaderInfoLog(fragShader_);

	shaderProgram_Texture_Phong_ = glCreateProgram();

	glAttachShader(shaderProgram_Texture_Phong_, vertShader_);
	glAttachShader(shaderProgram_Texture_Phong_, fragShader_);

	glLinkProgram(shaderProgram_Texture_Phong_);
	// As a result of a successful link operation, all active user-defined 
	// uniform variables belonging to program will be initialized to 0, 
	// and each of the program object's active uniform variables will be assigned 
	// a location that can be queried by calling glGetUniformLocation. 
	// Also, any active user-defined attribute variables that have not been bound 
	// to a generic vertex attribute index will be bound to one at this time.

	// check errors
	printProgramInfoLog(shaderProgram_Texture_Phong_);

	// free memory - shaders text files are large
	free(vs);	
	free(fs);
	glDeleteShader(vertShader_);
	glDeleteShader(fragShader_);

	// bind the shader program
 	currentShaderProgram_ = shaderProgram_Texture_Phong_;
	currentShaderProgram_ = 0;
 	glUseProgram(currentShaderProgram_);

	loc_ = glGetUniformLocation(shaderProgram_Texture_Phong_, "isCylinder");	// must be called after glLinkProgram()
	// You must bind the shader before calling glUniform1i(). Binding is done with glUseProgram().
	glUniform1i(loc_, GL_FALSE);		// default render type : Bezier ribbon
}

// load and set up texture
void Gl_Window::setTexture(char *filename, GLuint texID)
{
	if (filename != NULL)
	{
		// if texture was already created, destroy it
		// do not create again, reuse the texture - the initial texture size must
		// be large enough to contain all later textures
		if (glIsTexture(textureID_[texID]))
			glDeleteTextures(1, &(textureID_[texID]));

		// set up texture
		glGenTextures(1, &(textureID_[texID]));
		glActiveTexture(GL_TEXTURE0 + texID);
		glBindTexture(GL_TEXTURE_2D, textureID_[texID]);

		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
		glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MIN_FILTER,GL_NEAREST);	
		glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MAG_FILTER,GL_NEAREST); // Disable Filtering

		// hook texture unit to shader uniform variables
// 		if (0 == texID)
// 		{
// 			loc_ = glGetUniformLocation(currentShaderProgram_, "envMap");
// 			glUniform1i(loc_, texID);	// assign to texture unit_0
// 		} 
// 		else if (1 == texID)
// 		{
// 			loc_ = glGetUniformLocation(currentShaderProgram_, "baseMap");
// 			glUniform1i(loc_, texID);	// assign to texture unit_1
// 		}

		// send texture date to GPU
		unsigned char* data;
		int iWidth, iHeight, iDepth;
		loadImage(filename, data, iWidth, iHeight, iDepth);
		if (iDepth == 3)
			glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB8, iWidth, iHeight, 0, GL_RGB, GL_UNSIGNED_BYTE, data);
			//gluBuild2DMipmaps(GL_TEXTURE_2D, GL_RGB8, iWidth, iHeight, GL_RGB, GL_UNSIGNED_BYTE, data);
		else if (iDepth == 4)
			glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA8, iWidth, iHeight, 0, GL_RGBA, GL_UNSIGNED_BYTE, data);
			//gluBuild2DMipmaps(GL_TEXTURE_2D, GL_RGBA8, iWidth, iHeight, GL_RGBA, GL_UNSIGNED_BYTE, data);
		else
			std::cout << "image file not right" << std::endl;

		delete [] data;		// texture has been set to GPU

		glBindTexture(GL_TEXTURE_2D, textureID_[texID]);	// activate it
	} 
	else
		std::cout << "no texture file name given" << std::endl;
}

void Gl_Window::draw() {
	if (!valid()) {
		valid(1);
		init();
	}
	if(h() == 0)	h(1);	// Prevent a divide by zero, when window is too short

	glViewport(0, 0, w(), h());
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	
	// draw back ground
// 	glMatrixMode(GL_PROJECTION);
// 	glLoadIdentity();
// 	gluOrtho2D(0, 1, 0, 1);
// 
// 	glMatrixMode(GL_MODELVIEW);
// 	glLoadIdentity();
// 	gluLookAt(0.0, 0.0, 0.0, 0.0, 0.0, -1.0, 0.0, 1.0, 0.0);

	//glUseProgram(0);	// default OpenGL
// 	glActiveTexture(GL_TEXTURE0);
// 	glBindTexture(GL_TEXTURE_2D, textureID_[0]);	// environment map texture
// 	glEnable(GL_TEXTURE_2D);
// 	glDisable(GL_LIGHTING);
// 
// 	glColor4f(1.0, 1.0, 1.0, 1.0);	// blend with texture color
// 	glBegin(GL_QUADS);
// 	{
// 		glTexCoord2f(0.15, 0.85); glVertex2i(0, 0);
// 		glTexCoord2f(0.85, 0.85); glVertex2i(1, 0);
// 		glTexCoord2f(0.85, 0.15); glVertex2i(1, 1);
// 		glTexCoord2f(0.15, 0.15); glVertex2i(0, 1);
// 	}
// 	glEnd();

	// reset viewport, viewing, project matrix
	camera_->PerspectiveDisplay(w(), h());

	// draw on top of the back ground
	//glClear(GL_DEPTH_BUFFER_BIT);

	glPushMatrix();
	glRotatef(rotateY, 0.0, 1.0, 0.0);
	// draw base mesh
	if ( isBaseMesh ) {
		//glUseProgram(0);	// default OpenGL
		if ( isNewBaseMesh ) {
			isNewBaseMesh = false;
			if ( glIsList(displayList_BaseMesh) ) {
				glDeleteLists(displayList_BaseMesh, 1);
			}
			displayList_BaseMesh = glGenLists(1);
			glNewList(displayList_BaseMesh, GL_COMPILE_AND_EXECUTE);
			mesh->DrawBaseMesh();
			glEndList();
		} 
		else {
			glCallList(displayList_BaseMesh);
		}
	}

		// draw weaving yarns
		//glUseProgram(currentShaderProgram_);
		// additional clipping planes
	glEnable(GL_LIGHTING);
	glEnable(GL_TEXTURE_2D);
	
	if ( isNewWeaving ) {
		isNewWeaving = false;	// flip flag
		if ( glIsList(displayList_Weaving) ) {
			glDeleteLists(displayList_Weaving, 1);
		}
		displayList_Weaving = glGenLists(1);
		glNewList(displayList_Weaving, GL_COMPILE_AND_EXECUTE);
		long temp = 0;
		FILE* tempFile = NULL;
		mesh->Weave(false, temp, tempFile);
		glEndList();
	} 
	else {
		glCallList(displayList_Weaving);
	}

	glPopMatrix();

	// draw titles
// 	gl_color(FL_RED);
// 	drawText_screen("+1 twist edge", 10, 980);
// 	gl_color(FL_BLACK);
// 	drawText_screen("0 twist edge", 10, 960);
// 	gl_color(FL_GREEN);
// 	drawText_screen("NULL edge", 10, 940);
}


// x and y are screen coordinates
inline void Gl_Window::drawText_screen(const char* str, const int x, const int y)
{
	// text
	gl_font(FL_HELVETICA_BOLD, 16);

	glPushAttrib(GL_ENABLE_BIT); // save OpenGL states
	glDisable(GL_LIGHTING);
	glDisable(GL_TEXTURE_2D);
	glDisable(GL_DEPTH_TEST); // disable depth buffer while drawing text,so text draws /over/ object.

	// save 3D and setup 2D projection
	glMatrixMode(GL_PROJECTION);
	glPushMatrix(); // save 3D projection
	glLoadIdentity();
	gluOrtho2D(0, w(), 0, h());

	// save 3D and setup 2D model-view transformation
	glMatrixMode(GL_MODELVIEW);
	glPushMatrix();
	glLoadIdentity();

	glRasterPos2i(x, y);
	gl_draw(str, strlen(str));

	// restore 3D state
	glPopAttrib(); // pop the attribute stack

	glPopMatrix(); // pop GL_MODELVIEW

	glMatrixMode(GL_PROJECTION);
	glPopMatrix(); // pop GL_PROJECTION
}

// handle events
int Gl_Window::handle(int event)
{
	if (event == FL_SHOW) {
		show();
	}
	else if (event == FL_SHORTCUT)
	{
		switch (Fl::event_key()) {
		case 'r':
			camera_->Reset();
			rotateY = 0;
			redraw();
			return 1;

		case  's':		// SimpleRemeshing
			SimpleRemeshing(mesh);
			std::cout << "apply SimpleRemeshing" << std::endl;
			mesh->RandomTwist();
			isNewBaseMesh = true;
			isNewWeaving = true;
			redraw();
			return 1;

		case  'd':		// DooSabin
			DooSabin(mesh);
			std::cout << "apply DooSabin subdivision" << std::endl;
			mesh->RandomTwist();
			isNewBaseMesh = true;
			isNewWeaving = true;
			redraw();
			return 1;

		case  'c':		// Catmull-Clark
			Catmull_Clark(mesh);
			std::cout << "apply Catmull-Clark subdivision" << std::endl;
			//std::cout << "1st twist-1 edge " << mesh->childrenEdges_.size() << std::endl;
			//mesh->RandomTwist();
			isNewBaseMesh = true;
			isNewWeaving = true;
			redraw();
			return 1;

		case  'l':		// Catmull-Clark: linear subdivision
			LinearSub(mesh);
			std::cout << "apply Catmull-Clark linear subdivision" << std::endl;
			//std::cout << "1st twist-1 edge " << mesh->childrenEdges_.size() << std::endl;
			//mesh->RandomTwist();
			isNewBaseMesh = true;
			isNewWeaving = true;
			redraw();
			return 1;

		case  'a':		// Catmull-Clark: averaging
			Average(mesh);
			std::cout << "apply Catmull-Clark averaging" << std::endl;
			isNewBaseMesh = true;
			isNewWeaving = true;
			redraw();
			return 1;

		case 'o':	// output weaving geometry
			mesh->OutputObj();
			return 1;

		case 'w':	// save .obj and .twist files
			mesh->SaveFiles();
			return 1;

		case 'y':
			glEnable(GL_CULL_FACE);
			redraw();
			return 1;

		case 'u':
			glDisable(GL_CULL_FACE);
			redraw();
			return 1;

		default:
			return Fl_Widget::handle(event);	
		}
	} 
	else{	// mouse & others
		if (camera_->HandleMouseEvent(event)) {
			//std::cout << "HandleMouseEvent executed!!!" << std::endl;
			redraw();
			return 1;
		}
		else {
			//std::cout << "FL_Widget::handle executed!!!" << std::endl;
			return Fl_Widget::handle(event);
		}
	}	
}
