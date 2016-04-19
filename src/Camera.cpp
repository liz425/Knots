#include "Camera.h"

#include "cyMatrix4.h"

#include <fl/gl.h>
#include <fl/glu.h>

#include <cstdio>
#include <cmath>

#include <iostream>

/* Angle Conversions & Constants */
#define PI 3.1415926535897
#define RAD2DEG (180.0/PI)
#define DEG2RAD (PI/180.0)
#define DegToRad(x) ((x)*DEG2RAD)
#define RadToDeg(x) ((x)*RAD2DEG)

float DeltaAzim;
float DeltaElev;
float LocalDeltaAzim;
float LocalDeltaElev;

int MouseStartX;
int MouseStartY;
int MousePrevX;
int MousePrevY;

float epsilon = 0.0001;

GLdouble MvMatrix[16];
GLdouble ProjMatrix[16];
GLint ViewPort[4];

int CameraMode = INACTIVE;

cyPoint3f PrevMousePos;

void RotateX(cyPoint3f& v, float degree){
	float c = cos( DegToRad(degree) );
	float s = sin( DegToRad(degree) );
	float v1 = v.y * c - v.z * s;
	float v2 = v.y * s + v.z * c;
	v.Set(v.x, v1, v2);
}

void RotateY(cyPoint3f& v, float degree){
	float c = cos( DegToRad(degree) );
	float s = sin( DegToRad(degree) );
	float v0 = v.x * c + v.z * s;
	float v2 = -v.x * s + v.z * c;
	v.Set(v0, v.y, v2);
}

/* 
* ArbitraryRotate() - rotate around an arbitrary coordinate system specified by
*                     U, V, & W
*/
void ArbitraryRotate(cyPoint3f& U, cyPoint3f& V, cyPoint3f& W, 
					 float degreeX, float degreeY, 
					 cyPoint3f& point, cyPoint3f& aim) {
	float cx = cos( DegToRad(degreeX) );
	float sx = sin( DegToRad(degreeX) );
	float cy = cos( DegToRad(degreeY) );
	float sy = sin( DegToRad(degreeY) ); 

	cyMatrix4f trans(1, 0, 0, -aim.x,
					0, 1, 0, -aim.y,
					0, 0, 1, -aim.z,
					0, 0, 0, 1);
	trans.SetTranspose();

	cyMatrix4f mat(U.x, U.y, U.z, 0,
		V.x, V.y, V.z, 0,
		W.x, W.y, W.z, 0,
		0, 0, 0, 1);
	mat.SetTranspose();

	cyMatrix4f rot;
	cyPoint4f pos(point.x, point.y, point.z, 1);

	pos = trans * pos;

	pos = mat * pos;

	rot.Set(1,   0,  0, 0,
		0,  cx, sx, 0,
		0, -sx, cx, 0,
		0,   0,  0, 1);
	rot.SetTranspose();

	pos = rot * pos;

	rot.Set( cy, 0, sy, 0,
		0, 1,  0, 0,
		-sy, 0, cy, 0,
		0, 0,  0, 1);
	rot.SetTranspose();

	pos = rot * pos;

	pos = mat.GetInverse() * pos;

	pos = trans.GetInverse()  *pos;

	point.Set(pos.x, pos.y, pos.z);
}



/* constructors */

// default constructor... sets position to 0, 0, 5, aimed at the origin
// with the up vector set to the y axis
Camera::Camera() {
	Pos.Set(0, 0, 10);
	Aim.Set(0, 0, 0);
	Up.Set(0, 1, 0);

	// set default view volume
	NearPlane = 0.1;
	FarPlane = 1000.0;
	Fov = 60.0;
	defaultFov = 60.0;

	Initialize();
}

/* 
* constructor to set a camera to a desired orientation
* P is position in 3D
* A is the aim coordinate
* U is the up vector 
*/
Camera::Camera(cyPoint3f& P, cyPoint3f& A, cyPoint3f& U) {
	cyPoint3f dir = (A - P).GetNormalized();
	cyPoint3f up = U.GetNormalized();
	float dot = dir.Dot(up);

	// up vector and aim vector aren't perpendicular
	if (abs(dot) > epsilon) {
		fprintf (stderr, "Improper camera orientation. Can't create camera!\n");
		this->~Camera();
		return;
	} else if (dir.z != 0.0) {
		up.z = -(up.x*dir.x + up.y*dir.y)/dir.z;
	}

	Pos = P;
	Aim = A;
	Up = up;

	// set default view volume
	NearPlane = 0.1;
	FarPlane = 1000.0;
	Fov = 60.0;
	defaultFov = 60.0;

	Initialize();
}

/*
* Constructor setting up all values
*/
Camera::Camera(cyPoint3f& P, cyPoint3f& A, cyPoint3f& U,
			   float Near, float Far, float ViewAngle) {
	// orientation of camera isn't perpendicular
	if (((A-P).GetNormalized()).Dot(U.GetNormalized()) > epsilon) {
		fprintf (stderr, "Improper camera orientation. Can't create camera!\n");
		this->~Camera();
		return;
	}

	Pos = P;
	Aim = A;
	Up = U.GetNormalized();

	NearPlane = Near;
	FarPlane = Far;
	Fov = ViewAngle;
	defaultFov = ViewAngle;

	Initialize();
}

// Initialize routine setting up defaults
void Camera::Initialize() {
	cyPoint3f tmp, tmp1, tmp2;
	cyPoint3f axisOrigin, updatePos;
	float dist;

	DefaultPos = Pos;
	DefaultAim = Aim;
	DefaultUp = Up;

	// find the angle around the x axis 
	updatePos = Pos - Aim;
	axisOrigin.Set(updatePos.x, 0, 0);
	dist = (axisOrigin-updatePos).Length();
	tmp1.Set(updatePos.x, 0, dist);

	tmp = updatePos.GetNormalized();
	tmp1 = tmp1.GetNormalized();

	CurrentElev = RadToDeg(acos(tmp.Dot(tmp1)));

	// find the angle around the y axis
	axisOrigin.Set(0, updatePos.y, 0);
	dist = (axisOrigin-updatePos).Length();

	tmp2.Set(0, updatePos.y, dist);
	tmp2.Normalize();

	CurrentAzim = 360.0 - RadToDeg(acos(tmp2.Dot(tmp)));

	DefaultElev = CurrentElev;
	DefaultAzim = CurrentAzim;
}

// set functions for the Pos, Aim, and Up vectors....
// be careful with these because if you set either of them which causes the
// orientation of the camera to be un-orthogonal, then you'll see problems...

// just remember that (Aim - Pos).normalize() % Up == 0, or you'll see problems
void Camera::SetPos(cyPoint3f& P) {
	Pos = P;
}

void Camera::SetAim(cyPoint3f& A) {
	Aim = A;
}

void Camera::SetUp(cyPoint3f& U) {
	Up = U;
}


/*
* sets the near and far clipping planes for the camera view
*/
void Camera::SetClippingPlanes(float Near, float Far) {
	NearPlane = Near;
	FarPlane = Far;
}

/*
* sets the field of view of the camera, ViewAngle is in degrees
*/
void Camera::SetFOV(float ViewAngle) {
	Fov = ViewAngle;
}

/*
* resets the camera to its original orientation
*/
void Camera::Reset() {
	Pos = DefaultPos;
	Aim = DefaultAim;
	Up = DefaultUp;
	Fov = defaultFov;

	CurrentElev = DefaultElev;
	CurrentAzim = DefaultAzim;
}

/*
* sets the camera's aim to be the given vector v
*/
void Camera::SetCenterOfFocus(cyPoint3f& NewAim) {
	cyPoint3f& dif = NewAim - Aim;

	Aim = NewAim;
	Pos = Pos + dif;
}

/*
* draws an OpenGL window with the camera orientation
* W and H are the width and height of the window respectively
*/
void Camera::PerspectiveDisplay(int W, int H) {

	// set up the projection matrix
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluPerspective(Fov, float(W)/float(H), NearPlane, FarPlane);

	// set up the modelview matrix
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();

	// The matrix generated by gluLookAt postmultiplies the current matrix
	gluLookAt(Pos.x, Pos.y, Pos.z, 
		Aim.x, Aim.y, Aim.z,
		Up.x, Up.y, Up.z);
	
	// reset the view post to the whole window
	glViewport(0, 0, W, H);
}

/*
* mouse event handler function... should be called in the 
* mouse event handler function of your own code
*/
int Camera::HandleMouseEvent(int event) {
	GLdouble realy, wx, wy, wz, tempX, tempY, tempZ;
	int mouse_dx, mouse_dy, d;
	float z;
	cyPoint3f MousePos, dir;
	cyPoint3f WindowX, WindowY, WindowZ;

	switch(event) {
		case FL_PUSH:
			// set the new mouse state
			MouseStartX = MousePrevX = Fl::event_x();
			MouseStartY = MousePrevY = Fl::event_y();

			// alt key and mouse button have been pressed, camera will move
			switch (Fl::event_button()) {
				case FL_LEFT_MOUSE:	// rotating camera
					CameraMode = ROTATE;
					break;

				case FL_MIDDLE_MOUSE:	// translating camera:
					CameraMode = TRANSLATE;
					// get the modelview and projection matrices for projection
					// of the mouse's cursor into 3D
					glGetIntegerv(GL_VIEWPORT, ViewPort);
					glGetDoublev(GL_MODELVIEW_MATRIX, MvMatrix);
					glGetDoublev(GL_PROJECTION_MATRIX, ProjMatrix);

					// viewport[3] is height of window in pixels
					realy = ViewPort[3] - Fl::event_y() - 1;

					// project the aim of the camera into window coordinates
					// only concerned about getting the depth (wz) from here
					gluProject(Aim.x, Aim.y, Aim.z, 
						MvMatrix, ProjMatrix, ViewPort,
						&wx, &wy, &wz);

					// from the depth found from the previous call, project the
					// mouse coordinates into 3D coordinates
					gluUnProject(GLdouble(Fl::event_x()), GLdouble(realy), wz,
						MvMatrix, ProjMatrix, ViewPort, 
						&tempX, &tempY, &tempZ);
					PrevMousePos.x = tempX;
					PrevMousePos.y = tempY;
					PrevMousePos.z = tempZ;
					break;

				case FL_RIGHT_MOUSE:	// zooming camera:
					CameraMode = ZOOM;
					break;
			}			  
			return 1;

		case FL_DRAG:
			if (CameraMode != INACTIVE) {
				// find the greatest change in mouse position 
				mouse_dx = Fl::event_x() - MousePrevX;
				mouse_dy = Fl::event_y() - MousePrevY;

				if (abs(mouse_dx) > abs(mouse_dy)) 
					d = mouse_dx;
				else
					d = mouse_dy;

				switch (CameraMode) {
					case ZOOM:	// camera is zooming in
						z = float(d) / 100.0;
						dir = Aim - Pos;
						if (dir.Length() < 0.1 && z > 0) {
							// move the aim position too when you get in really close
							z *= 10.0;
							Aim = Aim + z*dir;
						}
						// update the new position
						Pos = Pos + z*dir;
						break;
					
					case ROTATE:	// camera is rotating
						// get rate of change in screen coordinates from when the 
						// mouse was first pressed
						DeltaAzim = float(Fl::event_x() - MouseStartX) / 5.0;
						DeltaElev = float(Fl::event_y() - MouseStartY) / 5.0;

						// get rate of change in screen coordinate from prev mouse pos
						LocalDeltaAzim = float(mouse_dx) / 5.0;
						LocalDeltaElev = float(mouse_dy) / 5.0;

						// rotate the window coordinate system by the rate of change
						// from the onset of the mouse event

						// got this small section of code from Dr. House
						WindowX.Set(1, 0, 0);
						WindowY.Set(0, 1, 0);

						RotateX(WindowX, CurrentElev+DeltaElev);
						RotateY(WindowX, CurrentAzim+DeltaAzim);
						WindowX.z = -WindowX.z;

						RotateX(WindowY, CurrentElev+DeltaElev);
						RotateY(WindowY, CurrentAzim+DeltaAzim);
						WindowY.z = -WindowY.z;

						WindowZ = (WindowX.Cross(WindowY)).GetNormalized();

						ArbitraryRotate(WindowX, WindowY, WindowZ, 
							LocalDeltaElev, 0, Pos, Aim);

						ArbitraryRotate(cyPoint3f(1, 0, 0), cyPoint3f(0, 1, 0), 
							cyPoint3f(0, 0, 1), 0, -LocalDeltaAzim, Pos, Aim);

						Up = WindowY.GetNormalized();
						break;
						
					case TRANSLATE:	// camera is translating
						realy = ViewPort[3] - Fl::event_y() - 1;

						gluProject(Aim.x, Aim.y, Aim.z, 
							MvMatrix, ProjMatrix, ViewPort, 
							&wx, &wy, &wz);

						gluUnProject(GLdouble(Fl::event_x()), GLdouble(realy), wz,
							MvMatrix, ProjMatrix, ViewPort, 
							&tempX, &tempY, &tempZ);
						MousePos.x = tempX;
						MousePos.y = tempY;
						MousePos.z = tempZ;

						// move both the camera position and its aim coordinate
						dir = MousePos - PrevMousePos;
						Pos = Pos - dir;
						Aim = Aim - dir;

						PrevMousePos = MousePos;
						break;
				}
				MousePrevX = Fl::event_x();
				MousePrevY = Fl::event_y();
			}
			return 1;

		case FL_RELEASE:
			if (CameraMode != INACTIVE) {
				// update the elevation and roll of the camera
				CurrentElev += DeltaElev;
				CurrentAzim += DeltaAzim;

				//printf("%f %f\n", CurrentElev, CurrentAzim);

				// reset the change in elevation and roll of the camera
				DeltaElev = DeltaAzim = 0.0;

				CameraMode = INACTIVE;
			}
			return 1;

		case FL_MOUSEWHEEL:
			SetFOV(Fov+Fl::event_dy());	// vertical wheel
			return 1;

		default:
			return 0;
			//return Fl_Widget::handle(event);
	}
}

// equals operator
const Camera& Camera::operator=(const Camera& Cam) {
	Aim = Cam.Aim;
	Pos = Cam.Pos;
	Up = Cam.Up;

	return *this;
}
