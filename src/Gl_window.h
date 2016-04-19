#ifndef Gl_Window_h
#define Gl_Window_h

#include "WeavingObject.h"
#include "Camera.h"

#include "GLSL.h"

#include <fl/Fl.H>
#include <fl/Fl_Gl_Window.H>

#define TIMESTEP 0.01	// for animation

class Gl_Window : public Fl_Gl_Window {
public:
	Gl_Window(int X,int Y,int W,int H,const char*L=0) : Fl_Gl_Window(X,Y,W,H,L) {
		mesh= new WeavingObject;
		isNewWeaving = true;
		isNewBaseMesh = true;
		isBaseMesh = true;
		rotateY = 0;

		camera_ = new Camera(cyPoint3f(0,0,4), cyPoint3f(0,0,0), cyPoint3f(0,1,0), 
							0.1, 100, 45);
		vertShader_ = 0;
		fragShader_ = 0;
		shaderProgram_Texture_Phong_ = 0;
		currentShaderProgram_ = 0;
		loc_ = 0;
		textureID_[0] = 0;	// 0 is reserved by OpenGL texturing
		textureID_[1] = 0;
		displayList_Weaving = 0;	// must be generated before use
		displayList_BaseMesh = 0;	// 0 is for no list or error
	}
	virtual ~Gl_Window(void) {}

	bool init();	// initialization
	void setShaders();
	void setTexture(char* filename, GLuint texID);	// load and set up texture
	void draw();	// draw the widget
	void drawText_screen(const char*, const GLint x, const GLint y); // x and y are screen coordinates
	int handle(int event);		// handle events

	WeavingObject* mesh;
	bool isNewWeaving;		// flag : new display list for the weaving yarns
	bool isNewBaseMesh;		// flag : new display list for the base mesh
	bool isBaseMesh;		// flag : draw base mesh
	int rotateY;			// rotation angle around y axis in degree

private:
	Camera* camera_;
	GLuint vertShader_;
	GLuint fragShader_;
	GLuint shaderProgram_Texture_Phong_;
	GLint currentShaderProgram_;
	GLint loc_;		// location in shader
	GLuint textureID_[2];	// 0: environment map
							// 1: object texture
	GLuint displayList_Weaving;
	GLuint displayList_BaseMesh;
};


#endif