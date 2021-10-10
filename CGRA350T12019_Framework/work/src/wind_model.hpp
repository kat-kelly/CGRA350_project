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
	int iter = 1;

	float dt;
	float diffusion;
	float viscosity;

	float* s;
	float* density;
	float* Vx0;
	float* Vy0;
	float* Vz0;

	void set_bnd(int b, float* x);
	void lin_solve(int b, float* x, float* x0, float a, float c);

	void diffuse(int b, float* x, float* x0, float diff);
	void advect(int b, float* d, float* d0, float* velocX, float* velocY, float* velocZ);
	void project(float* velocX, float* velocY, float* velocZ, float* p, float* div);

	//// Method for stepping though our simulation and updating the wind_field
	//void step();

	//// Functions used for calculating the simplified Navier-Stokes equation used for fluid simulation
	//void addForces(int x, int y, int z, vec3 amount);

	//void diffuse(int b, vec3* x, vec3* x0, float diff);
	//void diffuse(int b, float* x, float* x0, float diff);

	//void lin_solve(int b, vec3* x, vec3* x0, float a, float c);
	//void lin_solve(int b, float* x, float* x0, float a, float c);

	//void project(vec3* veloc, vec3* pdiv);
	//void advect(int b, float* d, float* d0, vec3* veloc);

	//void set_bnd(int b, vec3* x);
	//void set_bnd(int b, float* x);

public:
	void addVelocity(vec3 index, vec3 amount);
	void addDen(vec3 index, float amount);

	void step();

	float* Vx;
	float* Vy;
	float* Vz;

	wind_field(int size, float diffusion, float vicosity, float dt);
};

class wind_model {
private:
	// display
	bool w_show = true;
	GLuint shader = 0;
	vec3 color = vec3(0, 1, 0);
	float S = 1;

	//sim speed
	double old_t = glfwGetTime();
	double frameCap = 0.01;

	// fields parameter
	int N = 10;
	float visc = 0;
	float diff = 0.00000000001;
	float d_time = 0.1;

	// external force
	int aa = 1;

public:
	// Fluid wind field which has our vectors.
	wind_field* w_field = new wind_field(N, visc, diff, d_time);

	// Setup
	wind_model();

	//Step method used for simulate and updating the wind_field vectors with options for displaying wind visualization.
	void draw(const mat4& view, const mat4& proj);
};