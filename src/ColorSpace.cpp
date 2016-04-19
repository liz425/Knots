#include "ColorSpace.h"

#include <cstdio>
#include <cmath>



/*
HSV to RGB
Get some easy work done:
If Value V = 0, then we are done, color is black set R,G,B to 0
If Saturation S = 0, no color is dominant, set to some gray color
Hue is valued from 0 to 360, we chunk the space into 60 degree increments. At each 60 degrees we use a slightly different formula.
In general we will assign and set R,G and B exclusively as:
We set the most dominant color:
If H is 300 -> 60 , set R = V
If H is 60 -> 180, set G = V
If H is 180 -> 300, set B = V
The least dominant color is set as: pv = Value * ( 1 - Saturation )
The last remaining color is set as either:
qv = Value * ( 1 - Saturation * (Hue/60) - floor(Hue/60))
tv = Value * ( 1 - Saturation * ( 1 - ((Hue/60) - floor(Hue/60))))
Clean up, here we allow for i to be -1 or 6 just in case we have a very small floating point error otherwise, we have an undefined input.
Normalize R,G and B from 0 to 255.
Note: S and V are normalized between 0 and 1 in this example. H is its normal 0 to 360.
*/
bool HSVtoRGB(float H, float S, float V, float& R, float& G, float& B)
{
	if( V == 0 ) {
		R = 0;	
		G = 0;	
		B = 0;
	}
	else if( S == 0 ) {
		R = V;
		G = V;
		B = V;
	}
	else {
		const double hf = H / 60.0;
		const int       i  = (int) floor( hf );
		const double f  = hf - i;
		const double pv  = V * ( 1 - S );
		const double qv  = V * ( 1 - S * f );
		const double tv  = V * ( 1 - S * ( 1 - f ) );
		switch( i ) {
			// Red is the dominant color
			case 0:
				R = V;
				G = tv;                                                      
				B = pv;                                                      
				break;
			case 5:                                                         
				R = V;                                                        
				G = pv;                                                       
				B = qv;                                                      
				break;  

			// Green is the dominant color
			case 1:                                                        
				R = qv;                                                      
				G = V;                                                       
				B = pv;                                                      
				break;                                                       
			case 2:                                                        
				R = pv;                                                      
				G = V;                                                       
				B = tv;                                                      
				break;                                                       

			// Blue is the dominant color
			case 3:                                                        
				R = pv;                                                       
				G = qv;                                                      
				B = V;                                                       
				break;                                                       
			case 4:                                                        																
				R = tv;                                                      
				G = pv;                                                      
				B = V;                                                      
				break;                                                                                                    

			// Just in case we overshoot on our math by a little, we put these here. Since its a switch it won't slow us down at all to put these here.
			case 6:                                                       
				R = V;                                                       
				G = tv;                                                     
				B = pv;                                                      
				break;                                                       
			case -1:                                                       
				R = V;                                                        
				G = pv;                                                      
				B = qv;                                                       
				break;                                                       

			// The color is not defined, we should throw an error.
			default:                                                        
				printf("i Value error in Pixel conversion, Value is %d\n", i);   
				return false;                                                      
		}                                                               
	}                                                                  
	R *= 255.0F;                                                      
	G *= 255.0F;                                                        
	B *= 255.0F;
	return true;
}



bool RGBtoHSV(float R, float G, float B, float& H, float& S, float& V)
{
	return true;
}