#ifndef Camera_H
#define Camera_H

#include "cyPoint.h"

#include "fl/Fl.H"
#include "fl/Fl_Gl_Window.H"

#define INACTIVE 0
#define TRANSLATE 1
#define ROTATE 2
#define ZOOM 3

class Camera {
private:
	// all variables starting with 'Default' hold the initial camera values
	// these values are used if the camera is reset to its initial state
	cyPoint3f DefaultPos;
	cyPoint3f DefaultAim;
	cyPoint3f DefaultUp;

	float DefaultAzim;
	float DefaultElev;

	float CurrentAzim;
	float CurrentElev;

	float NearPlane;
	float FarPlane;
	float Fov;
	float defaultFov;

	void Initialize();

public:
	cyPoint3f Pos;
	cyPoint3f Aim; 
	cyPoint3f Up;

	// constructors

	// default constructor
	Camera();

	// constructor setting up camera orientation
	// P is position in 3D, A is the aim coordinate in 3D, U is the up vector
	Camera(cyPoint3f& P, cyPoint3f& A, cyPoint3f& U);

	// constructor setting up camera orientation and view volume
	// P is position in 3D, A is aim coordinate in 3D, U is up vector
	// Near is near clipping plane, Far is far clipping plane, 
	// ViewAngle is field of view angle in degrees
	Camera(cyPoint3f& P, cyPoint3f& A, cyPoint3f& U, 
		float Near, float Far, float ViewAngle);

	// sets the clipping planes of the view volume
	void SetClippingPlanes(float Near, float Far);

	// sets the FOV, ViewAngle should be in degrees
	void SetFOV(float ViewAngle);	

	// set routines for Pos, Aim, and Up vector
	void SetPos(cyPoint3f& P);
	void SetAim(cyPoint3f& A);
	void SetUp(cyPoint3f& U);

	// reset the camera to its initial position
	void Reset();

	// focus camera to some input aim position
	void SetCenterOfFocus(cyPoint3f& NewAim);

	// function to use the camera as the OpenGL camera
	// W and H are the width and height of the window
	void PerspectiveDisplay(int W, int H);

	// function that handles mouse events
	int HandleMouseEvent(int event);

	const Camera& operator=(const Camera& cam);
};

#endif



