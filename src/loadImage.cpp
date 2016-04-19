#include <cstdlib>
#include <cstdio>
#include <cstring>	// mmecpy()
#include <fl/Fl_Shared_Image.H>
#include "loadImage.h"

using namespace std;

int loadImage ( const char *filename, unsigned char * &image, int &imageWidth, int &imageHeight , int &imageDepth)
{
	fl_register_images();
	Fl_Shared_Image *img = Fl_Shared_Image::get(filename);
	if ( img->d() == 0 ) {
		perror(filename);
		exit(1);
	}
	printf("image file: %s\n", filename);

	imageWidth = img->w();
	imageHeight = img->h();
	imageDepth = img->d();
	size_t imageSize = imageWidth * imageHeight * imageDepth;	// in byte
	printf("image depth %d\n", imageDepth);
	image = new unsigned char[imageSize];

	if (img->count() == 1)	// bitmap and color images
		memcpy(image, img->data()[0], imageSize);
	else
		printf("Not supported: count=%d\n", img->count());

	img->release();
	return 1;
}
