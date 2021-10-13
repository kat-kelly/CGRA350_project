#version 330 core

// uniform data
uniform mat4 uProjectionMatrix;
uniform mat4 uModelViewMatrix;
uniform vec3 uColor;

// viewspace data (this must match the output of the fragment shader)
//in VertexData {
//	vec3 position;
//	vec3 normal;
//	vec2 textureCoord;
//} f_in;
in GeometryData {
	vec3 position;
	vec3 normal;
	vec2 textureCoord;
} g_in;

// framebuffer output
out vec4 fb_color;

void main() {
	// calculate lighting (hack)
	vec3 eye = normalize(-g_in.position);
	float light = abs(dot(normalize(g_in.normal), eye));
	vec3 color = mix(uColor / 4, uColor, light);

	// output to the frambuffer
	fb_color = vec4(color, 1);
}