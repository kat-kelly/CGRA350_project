#version 410 core
#extension GL_ARB_geometry_shader4 : enable

layout(lines_adjacency) in;
layout(triangle_strip, max_vertices = 50) out;

void main() {
	// lower triangle
	gl_Position = gl_in[0].gl_Position;
	EmitVertex();

	gl_Position = gl_in[1].gl_Position;
	EmitVertex();

	gl_Position = gl_in[2].gl_Position;
	EmitVertex();

	EndPrimitive();

	// upper triangle
	gl_Position = gl_in[0].gl_Position;
	EmitVertex();

	gl_Position = gl_in[2].gl_Position;
	EmitVertex();

	gl_Position = gl_in[3].gl_Position;
	EmitVertex();

	EndPrimitive();
}