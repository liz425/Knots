#ifndef GLSL_h
#define GLSL_h

#include "gl/glew.h"
#define printOpenGLError() printOglError(__FILE__, __LINE__)

int printOglError(char *file, int line);

// vertex / fragment shader
void printShaderInfoLog(GLuint obj);
void printProgramInfoLog(GLuint obj);

#endif