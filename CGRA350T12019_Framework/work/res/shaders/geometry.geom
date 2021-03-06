#version 410 core
#extension GL_ARB_geometry_shader4 : enable

layout(lines_adjacency) in;
layout(triangle_strip, max_vertices = 50) out;

in TesselationEData {
	vec3 position;
	vec3 normal;
	vec2 textureCoord;
} g_in[];

out GeometryData {
	vec3 position;
	vec3 normal;
	vec2 textureCoord;
} g_out;

void main() {

	// lower triangle
	g_out.position = g_in[0].position;
	g_out.normal = g_in[0].normal;
	g_out.textureCoord = g_in[0].textureCoord;
	gl_Position = gl_in[0].gl_Position;
	EmitVertex();

	g_out.position = g_in[1].position;
	g_out.normal = g_in[1].normal;
	g_out.textureCoord = g_in[1].textureCoord;
	gl_Position = gl_in[1].gl_Position;
	EmitVertex();

	g_out.position = g_in[2].position;
	g_out.normal = g_in[2].normal;
	g_out.textureCoord = g_in[2].textureCoord;
	gl_Position = gl_in[2].gl_Position;
	EmitVertex();

	EndPrimitive();

	// upper triangle
	g_out.position = g_in[0].position;
	g_out.normal = g_in[0].normal;
	g_out.textureCoord = g_in[0].textureCoord;
	gl_Position = gl_in[0].gl_Position;
	EmitVertex();

	g_out.position = g_in[2].position;
	g_out.normal = g_in[2].normal;
	g_out.textureCoord = g_in[2].textureCoord;
	gl_Position = gl_in[2].gl_Position;
	EmitVertex();

	g_out.position = g_in[3].position;
	g_out.normal = g_in[3].normal;
	g_out.textureCoord = g_in[3].textureCoord;
	gl_Position = gl_in[3].gl_Position;
	EmitVertex();

	EndPrimitive();
}