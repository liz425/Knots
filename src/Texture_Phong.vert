varying vec3 N;		// vertx normal in eye space
varying vec3 v;		// vertex position in eye space
varying vec3 L;		// light direction in eye space
varying vec4 diffuse, ambient, specular;
varying vec2 TexcoordEnv, Texcoord;

void main(void)
{          
   N = normalize(gl_NormalMatrix * gl_Normal);
   v = vec3(gl_ModelViewMatrix * gl_Vertex);
   L = normalize(gl_LightSource[0].position.xyz - v);
   
   diffuse	= gl_FrontMaterial.diffuse * gl_LightSource[0].diffuse;
   ambient	= gl_FrontMaterial.ambient * gl_LightSource[0].ambient;
   specular	= gl_FrontMaterial.specular * gl_LightSource[0].specular;
   		
   TexcoordEnv.x = (0.9 * dot(N, vec3(1.0, 0.0, 0.0)) + 1.0) / 2.0;
   TexcoordEnv.y = (0.9 * dot(N, vec3(0.0, -1.0, 0.0)) + 1.0) / 2.0;
	   
   Texcoord = gl_MultiTexCoord0.xy;
   
   gl_Position = ftransform();
}