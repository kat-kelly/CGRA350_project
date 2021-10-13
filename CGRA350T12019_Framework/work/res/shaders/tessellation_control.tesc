#version 430 core
#extension GL_ARB_tessellation_shader : enable
//#extension GL_NV_gpu_program5 : enable

//uniform vec3 cameraPos;
// uniform data
//uniform mat4 uProjectionMatrix;
//uniform mat4 uModelViewMatrix;
//uniform vec3 uColor;

layout (vertices = 4) out;

void main()
{
    if(gl_InvocationID == 0)
    {
        float dist = length(gl_in[0].gl_Position.xyz - gl_in[1].gl_Position.xyz);

        gl_TessLevelOuter[0] = 1.0;
        gl_TessLevelOuter[1] = dist;
    }

    gl_out[gl_InvocationID].gl_Position = gl_in[gl_InvocationID].gl_Position;
}