#version 410 core

//in VertexData {
//	vec3 position;
//	vec3 normal;
//	vec2 textureCoord;
//} f_in[gl_MaxPatchVertices];

patch in float gl_TessLevelOuter[4];
patch in float gl_TessLevelInner[2];

layout (quads) in;

void main(void)
{
    // interpolate bottom
    vec4 p1 = mix(gl_in[0].gl_Position, gl_in[1].gl_Position, gl_TessCoord.x);
    // interpolate top
    vec4 p2 = mix(gl_in[2].gl_Position, gl_in[3].gl_Position, gl_TessCoord.x);

    gl_Position = mix(p1, p2, gl_TessCoord.y);
}