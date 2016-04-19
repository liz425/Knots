uniform sampler2D envMap;
uniform sampler2D baseMap;
uniform bool isCylinder;		// flag

varying vec3 N;
varying vec3 v;   
varying vec3 L;
varying vec4 diffuse, ambient, specular;
varying vec2 TexcoordEnv, Texcoord;

void main (void)  
{   
   vec3 NN = normalize(N); // normalize the interpolated normal vector
   vec3 LL = normalize(L); // normalize the interpolated light direction vector
   vec3 E = normalize(-v); // we are in Eye Coordinates, so EyePos is (0,0,0) 
   vec3 R = normalize(-reflect(LL,NN));    

   vec4  envColor = texture2D( envMap, TexcoordEnv );
   vec4  baseColor = texture2D( baseMap, Texcoord );
   
   if(isCylinder) {
		baseColor = envColor;
   }
   else {
		baseColor = envColor * baseColor;
   }
   
   // calculate Diffuse Term:  
   vec4 Idiff = diffuse * max((dot(NN,LL)+1.0)/2.0, 0.0);    
   Idiff = clamp(Idiff, 0.0, 1.0); 
   
   // calculate Specular Term:
   vec4 Ispec = specular * pow(max(dot(R,E),0.0), 0.3*gl_FrontMaterial.shininess);
	
   // write Total Color:  
   gl_FragColor = ambient *baseColor + Idiff * baseColor + Ispec;  
   //gl_FragColor = vec4(baseColor.rgb, 1.0); 
}
