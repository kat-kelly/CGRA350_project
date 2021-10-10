#pragma once
// glm
#include <glm/glm.hpp>

// project
#include "opengl.hpp"
using namespace glm;

// Wind_field data which can be used by other parts of the prject
class wind_field {
private:
	int N;
	int iter = 10;

	// timesetp viscocity diffusion
	float dt;
	float diffusion;
	float viscosity;

	// density previouse velocities 
	float* s;
	float* density;
	float* Vx0;
	float* Vy0;
	float* Vz0;

	// Method for stepping though our simulation and updating the wind_field
	void set_Boundaries(int b, float* x);
	void lin_solve(int b, float* x, float* x0, float a, float c);

	void diffuse(int b, float* x, float* x0, float diff);
	void advect(int b, float* d, float* d0, float* velocX, float* velocY, float* velocZ);
	void project(float* velocX, float* velocY, float* velocZ, float* p, float* div);

public:
	void addVelocity(vec3 index, vec3 amount);
	void addDen(vec3 index, float amount);

	void step();

	// current velocitys
	float* Vx;
	float* Vy;
	float* Vz;

	wind_field(int size, float diffusion, float vicosity, float dt);
};

class wind_model {
private:
	// display
	GLuint shader = 0;
	vec3 color = vec3(0.6, 0.6, 0.6);
	float S = 1;

	// fields parameter
	int N = 14;
	float visc = 0.1;
	float diff = 0;
	float d_time = 0.05;

	//sim speed
	double old_t = glfwGetTime();
	double frameCap = 0.01;

	bool D = true;
	float addon = 0;

	// stepping though
	void simulate();

	// for display visuals
	void draw(const mat4& view, const mat4& proj);

public:
	// Fluid wind field which has our vectors.
	wind_field* w_field = new wind_field(N, visc, diff, d_time);

	// 0 = Off
	// 1 = On
	// 2 = visualize
	int display = 0;

	float w_angle = 0;
	float w_strength = 1;
	float pulse = 0;

	// Setup
	wind_model();

	//Step method used for simulate and updating the wind_field vectors with options for displaying wind visualization.
	void run(const mat4& view, const mat4& proj);
};