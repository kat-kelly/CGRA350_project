#version 410 core
#extension GL_ARB_tessellation_shader : enable
#extension GL_EXT_shader_io_blocks : enable

uniform vec4 cameraPos;

// output level of detail
patch out float gl_TessLevelOuter[4];
patch out float gl_TessLevelInner[2];

void main() {
//distance(gl_in[1].position
//clamp(distance(aPosition, cameraPos)/10, 1, 50);
	// outer tesselation level
	gl_TessLevelOuter[0] = clamp(distance(gl_in[0].gl_Position, cameraPos)/10, 1, 50);
	gl_TessLevelOuter[1] = clamp(distance(gl_in[1].gl_Position, cameraPos)/10, 1, 50);
	gl_TessLevelOuter[2] = clamp(distance(gl_in[2].gl_Position, cameraPos)/10, 1, 50);
	gl_TessLevelOuter[3] = clamp(distance(gl_in[3].gl_Position, cameraPos)/10, 1, 50);

	// inner tesselation level
	gl_TessLevelInner[0] = 0.5 * (gl_TessLevelOuter[0] + gl_TessLevelOuter[3]);
	gl_TessLevelInner[1] = 0.5 * (gl_TessLevelOuter[2] + gl_TessLevelOuter[1]);
}