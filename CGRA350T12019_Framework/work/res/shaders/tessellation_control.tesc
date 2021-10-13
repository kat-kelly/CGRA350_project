#version 410 core
#extension GL_ARB_tessellation_shader : enable

// uniform data
uniform vec3 cameraPos;

layout (vertices = 4) out;

in VertexData {
	vec3 position;
	vec3 normal;
	vec2 textureCoord;
} tc_in[];

out TesselationCData {
	vec3 position;
	vec3 normal;
	vec2 textureCoord;
} tc_out[];

void main()
{
    tc_out[gl_InvocationID].position = tc_in[gl_InvocationID].position;
    tc_out[gl_InvocationID].normal = tc_in[gl_InvocationID].normal;
    tc_out[gl_InvocationID].textureCoord = tc_in[gl_InvocationID].textureCoord;
    if(gl_InvocationID == 0)
    {
	// distance between point on grass and camera
        float dist = distance(gl_in[0].gl_Position.xyz, gl_in[1].gl_Position.xyz);
        float camDist = distance(gl_in[0].gl_Position.xyz, cameraPos);

        gl_TessLevelOuter[0] = 1.0;
        gl_TessLevelOuter[1] = max(dist/camDist, 1.0); // will increase further from camera
    }
    
    gl_out[gl_InvocationID].gl_Position  = gl_in[gl_InvocationID].gl_Position;
}