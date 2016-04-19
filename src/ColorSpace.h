#ifndef COLOR_SPACE__
#define COLOR_SPACE__


// hue: 0~360		saturation: 0~1		value: 0~1
// red, green, blue: 0~255
bool HSVtoRGB(float H, float S, float V, float& R, float& G, float& B);
bool RGBtoHSV(float R, float G, float B, float& H, float& S, float& V);

#endif // COLOR_SPACE__
