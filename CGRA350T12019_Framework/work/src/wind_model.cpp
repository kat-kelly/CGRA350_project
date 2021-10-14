#define GLM_ENABLE_EXPERIMENTAL
#define INDEX(x, y, z) ((x) + (y) * N + (z) * N * N)

// C++
#include <iostream>
#include <math.h>


// glm
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>

//experimental
#include <glm/gtx/rotate_vector.hpp>
#include <glm/gtx/euler_angles.hpp>

// project
#include "cgra/cgra_geometry.hpp"
#include "wind_model.hpp"
#include <cgra/cgra_shader.hpp>

using namespace std;
using namespace glm;
using namespace cgra;



wind_model::wind_model()
{
	shader_builder sb;
	sb.set_shader(GL_VERTEX_SHADER, CGRA_SRCDIR + string("//res//shaders//color_vert.glsl"));
	sb.set_shader(GL_FRAGMENT_SHADER, CGRA_SRCDIR + string("//res//shaders//color_frag_wind.glsl"));
	shader = sb.build();


}
void wind_model::run(const mat4& view, const mat4& proj) {
	if (display < 1) return;
	simulate();

	if (display < 2) return;
	draw(view, proj);
}

void wind_model::simulate() {

	double t = glfwGetTime();
	if (t - old_t > frameCap) {
		for (int i = 3; i < N - 3; ++i) {
			vec3 index1 = vec3(2, 3, i);
			mat4 r = rotate(mat4(1), radians(w_angle), vec3(0, 1, 0));
			vec3 vel = r * vec4(w_strength, 0, 0, 0);

			if (pulse > 0) {
				float cV = w_strength + addon;
				vel = r * vec4(cV, 0, 0, 0);
				if (D && w_strength + addon >= w_strength + 1) D = false;
				if (!D && w_strength + addon <= -0) D = true;
				addon = (D) ? addon + pulse : addon - pulse;
			}

			w_field->addVelocity(index1, vel);
			w_field->addDensity(index1, length(vel) / 20);
		}

		w_field->step();
		old_t = t;
	}

}


void wind_model::draw(const mat4& view, const mat4& proj)
{
	glUseProgram(shader);
	glUniformMatrix4fv(glGetUniformLocation(shader, "uProjectionMatrix"), 1, false, value_ptr(proj));
	vec3 vec; // place holders
	mat4 mv;

	for (int k = 0; k < N; ++k) {
		for (int j = 0; j < N; ++j) {
			for (int i = 0; i < N; ++i) {
				float x = w_field->Vx[INDEX(i, j, k)];
				float y = w_field->Vy[INDEX(i, j, k)];
				float z = w_field->Vz[INDEX(i, j, k)];
				float d = w_field->density[INDEX(i, j, k)];
				if (d > 0) {
					//w_field->addDensity(vec3(i, j, k), -0.00001);
				}

				vec = vec3(x, y, z);// setup
				mv = view;

				vec3 norm = normalize(vec);
				float len = length(vec);

				mv = scale(mv, vec3(S));// scale everything

				mv = translate(mv, vec3(i + 0.5 - N / 2, j + 0.5, k + 0.5 - N / 2));	 // translate
				mv = mv * orientation(norm, vec3(0, 0, 1));		// rotation
				mv = scale(mv, vec3(0.01, 0.01, len));				// scale


				glUniform3fv(glGetUniformLocation(shader, "uColor"), 1, value_ptr(vec3(1 - len, len, 0)));
				glUniformMatrix4fv(glGetUniformLocation(shader, "uModelViewMatrix"), 1, false, value_ptr(mv));
				drawCone();

				if (balls) {
					mv = view;
					mv = scale(mv, vec3(S));// scale everything
					mv = translate(mv, vec3(i + 0.5 - N / 2, j + 0.5, k + 0.5 - N / 2));	 // translate
					mv = scale(mv, vec3(d, d, d));// scale
					glUniform3fv(glGetUniformLocation(shader, "uColor"), 1, value_ptr(vec3(1-d*10, 1-d*10, 1-d*10)));
					glUniformMatrix4fv(glGetUniformLocation(shader, "uModelViewMatrix"), 1, false, value_ptr(mv));
					drawSphere();
				}
				else {
					w_field->density[INDEX(i, j, k)] = 0;
				}
			}
		}

	}


}



wind_field::wind_field(int N, float diffusion, float viscosity, float dt)
{
	this->N = N;
	this->dt = dt;
	this->diffusion = diffusion;
	this->viscosity = viscosity;

	this->s = new float[N * N * N];
	this->density = new float[N * N * N];

	this->Vx = new float[N * N * N];
	this->Vy = new float[N * N * N];
	this->Vz = new float[N * N * N];

	this->Vx0 = new float[N * N * N];
	this->Vy0 = new float[N * N * N];
	this->Vz0 = new float[N * N * N];

	for (int i = 0; i < N * N * N; ++i) {
		s[i] = 0;
		density[i] = 0;

		Vx[i] = 0;
		Vy[i] = 0;
		Vz[i] = 0;

		Vx0[i] = 0;
		Vy0[i] = 0;
		Vz0[i] = 0;
	}


}

void wind_field::step()
{
	diffuse(1, Vx0, Vx, viscosity);
	diffuse(2, Vy0, Vy, viscosity);
	diffuse(3, Vz0, Vz, viscosity);

	project(Vx0, Vy0, Vz0, Vx, Vy);

	advect(1, Vx, Vx0, Vx0, Vy0, Vz0);
	advect(2, Vy, Vy0, Vx0, Vy0, Vz0);
	advect(3, Vz, Vz0, Vx0, Vy0, Vz0);

	project(Vx, Vy, Vz, Vx0, Vy0);

	diffuse(0, s, density, diffusion);
	advect(0, density, s, Vx, Vy, Vz);
}

void wind_field::addVelocity(vec3 index, vec3 amount)
{
	int i = INDEX(index.x, index.y, index.z);

	Vx[i] += amount.x;
	Vy[i] += amount.y;
	Vz[i] += amount.z;
}
void wind_field::addDensity(vec3 index, float amount)
{
	int i = INDEX(index.x, index.y, index.z);
	density[i] += amount;
}

void wind_field::diffuse(int b, float* x, float* x0, float diff)
{
	float a = dt * diff * (N - 2) * (N - 2);
	lin_solve(b, x, x0, a, 1 + 6 * a);
}

void wind_field::advect(int b, float* d, float* d0, float* velocX, float* velocY, float* velocZ)
{
	float i0, i1, j0, j1, k0, k1;

	float dtx = dt * (N - 2);
	float dty = dt * (N - 2);
	float dtz = dt * (N - 2);

	float s0, s1, t0, t1, u0, u1;
	float tmp1, tmp2, tmp3, x, y, z;

	float Nfloat = N;
	float ifloat, jfloat, kfloat;
	int i, j, k;

	for (k = 1, kfloat = 1; k < N - 1; ++k, kfloat++) {
		for (j = 1, jfloat = 1; j < N - 1; ++j, jfloat++) {
			for (i = 1, ifloat = 1; i < N - 1; ++i, ifloat++) {
				tmp1 = dtx * velocX[INDEX(i, j, k)];
				tmp2 = dty * velocY[INDEX(i, j, k)];
				tmp3 = dtz * velocZ[INDEX(i, j, k)];
				x = ifloat - tmp1;
				y = jfloat - tmp2;
				z = kfloat - tmp3;

				if (x < 0.5f) x = 0.5f;
				if (x > Nfloat + 0.5f) x = Nfloat + 0.5f;
				i0 = floorf(x);
				i1 = i0 + 1.0f;
				if (y < 0.5f) y = 0.5f;
				if (y > Nfloat + 0.5f) y = Nfloat + 0.5f;
				j0 = floorf(y);
				j1 = j0 + 1.0f;
				if (z < 0.5f) z = 0.5f;
				if (z > Nfloat + 0.5f) z = Nfloat + 0.5f;
				k0 = floorf(z);
				k1 = k0 + 1.0f;

				s1 = x - i0;
				s0 = 1.0f - s1;
				t1 = y - j0;
				t0 = 1.0f - t1;
				u1 = z - k0;
				u0 = 1.0f - u1;

				int i0i = i0;
				int i1i = i1;
				int j0i = j0;
				int j1i = j1;
				int k0i = k0;
				int k1i = k1;

				d[INDEX(i, j, k)] =

					s0 * (t0 * (u0 * d0[INDEX(i0i, j0i, k0i)]
						+ u1 * d0[INDEX(i0i, j0i, k1i)])
						+ (t1 * (u0 * d0[INDEX(i0i, j1i, k0i)]
							+ u1 * d0[INDEX(i0i, j1i, k1i)])))
					+ s1 * (t0 * (u0 * d0[INDEX(i1i, j0i, k0i)]
						+ u1 * d0[INDEX(i1i, j0i, k1i)])
						+ (t1 * (u0 * d0[INDEX(i1i, j1i, k0i)]
							+ u1 * d0[INDEX(i1i, j1i, k1i)])));
			}
		}
	}
	set_Boundaries(b, d);
}

void wind_field::project(float* velocX, float* velocY, float* velocZ, float* p, float* div)
{
	for (int k = 1; k < N - 1; ++k) {
		for (int j = 1; j < N - 1; ++j) {
			for (int i = 1; i < N - 1; ++i) {
				div[INDEX(i, j, k)] = -0.5f * (
					velocX[INDEX(i + 1, j, k)]
					- velocX[INDEX(i - 1, j, k)]
					+ velocY[INDEX(i, j + 1, k)]
					- velocY[INDEX(i, j - 1, k)]
					+ velocZ[INDEX(i, j, k + 1)]
					- velocZ[INDEX(i, j, k - 1)]
					) / N;
				p[INDEX(i, j, k)] = 0;
			}
		}
	}
	set_Boundaries(0, div);
	set_Boundaries(0, p);
	lin_solve(0, p, div, 1, 6);

	for (int k = 1; k < N - 1; ++k) {
		for (int j = 1; j < N - 1; ++j) {
			for (int i = 1; i < N - 1; ++i) {
				velocX[INDEX(i, j, k)] -= 0.5f * (p[INDEX(i + 1, j, k)]
					- p[INDEX(i - 1, j, k)]) * N;
				velocY[INDEX(i, j, k)] -= 0.5f * (p[INDEX(i, j + 1, k)]
					- p[INDEX(i, j - 1, k)]) * N;
				velocZ[INDEX(i, j, k)] -= 0.5f * (p[INDEX(i, j, k + 1)]
					- p[INDEX(i, j, k - 1)]) * N;
			}
		}
	}
	set_Boundaries(1, velocX);
	set_Boundaries(2, velocY);
	set_Boundaries(3, velocZ);
}

void wind_field::lin_solve(int b, float* x, float* x0, float a, float c)
{
	float cc = 1.0 / c;
	for (int k = 0; k < iter; ++k) {
		for (int m = 1; m < N - 1; ++m) {
			for (int j = 1; j < N - 1; ++j) {
				for (int i = 1; i < N - 1; ++i) {
					x[INDEX(i, j, m)] =
						(x0[INDEX(i, j, m)]
							+ a * (x[INDEX(i + 1, j, m)]
								+ x[INDEX(i - 1, j, m)]
								+ x[INDEX(i, j + 1, m)]
								+ x[INDEX(i, j - 1, m)]
								+ x[INDEX(i, j, m + 1)]
								+ x[INDEX(i, j, m - 1)]
								)) * cc;
				}
			}
		}
		set_Boundaries(b, x);
	}
}
void wind_field::set_Boundaries(int b, float* x)
{
	for (int j = 1; j < N - 1; ++j) {
		for (int i = 1; i < N - 1; ++i) {
			x[INDEX(i, j, 0)] = b == 3 ? -x[INDEX(i, j, 1)] : x[INDEX(i, j, 1)];
			x[INDEX(i, j, N - 1)] = b == 3 ? -x[INDEX(i, j, N - 2)] : x[INDEX(i, j, N - 2)];
		}
	}
	for (int k = 1; k < N - 1; ++k) {
		for (int i = 1; i < N - 1; ++i) {
			x[INDEX(i, 0, k)] = b == 2 ? -x[INDEX(i, 1, k)] : x[INDEX(i, 1, k)];
			x[INDEX(i, N - 1, k)] = b == 2 ? -x[INDEX(i, N - 2, k)] : x[INDEX(i, N - 2, k)];
		}
	}
	for (int k = 1; k < N - 1; ++k) {
		for (int j = 1; j < N - 1; ++j) {
			x[INDEX(0, j, k)] = b == 1 ? -x[INDEX(1, j, k)] : x[INDEX(1, j, k)];
			x[INDEX(N - 1, j, k)] = b == 1 ? -x[INDEX(N - 2, j, k)] : x[INDEX(N - 2, j, k)];
		}
	}

	x[INDEX(0, 0, 0)] = 0.33f * (x[INDEX(1, 0, 0)]
		+ x[INDEX(0, 1, 0)]
		+ x[INDEX(0, 0, 1)]);
	x[INDEX(0, N - 1, 0)] = 0.33f * (x[INDEX(1, N - 1, 0)]
		+ x[INDEX(0, N - 2, 0)]
		+ x[INDEX(0, N - 1, 1)]);
	x[INDEX(0, 0, N - 1)] = 0.33f * (x[INDEX(1, 0, N - 1)]
		+ x[INDEX(0, 1, N - 1)]
		+ x[INDEX(0, 0, N)]);
	x[INDEX(0, N - 1, N - 1)] = 0.33f * (x[INDEX(1, N - 1, N - 1)]
		+ x[INDEX(0, N - 2, N - 1)]
		+ x[INDEX(0, N - 1, N - 2)]);
	x[INDEX(N - 1, 0, 0)] = 0.33f * (x[INDEX(N - 2, 0, 0)]
		+ x[INDEX(N - 1, 1, 0)]
		+ x[INDEX(N - 1, 0, 1)]);
	x[INDEX(N - 1, N - 1, 0)] = 0.33f * (x[INDEX(N - 2, N - 1, 0)]
		+ x[INDEX(N - 1, N - 2, 0)]
		+ x[INDEX(N - 1, N - 1, 1)]);
	x[INDEX(N - 1, 0, N - 1)] = 0.33f * (x[INDEX(N - 2, 0, N - 1)]
		+ x[INDEX(N - 1, 1, N - 1)]
		+ x[INDEX(N - 1, 0, N - 2)]);
	x[INDEX(N - 1, N - 1, N - 1)] = 0.33f * (x[INDEX(N - 2, N - 1, N - 1)]
		+ x[INDEX(N - 1, N - 2, N - 1)]
		+ x[INDEX(N - 1, N - 1, N - 2)]);
}

//void wind_field::step() {
//	diffuse(1, V0, V, viscosity);
//
//	project(V0, V);
//
//	advect(1, Vx, Vx0, Vx0, Vy0, Vz0);
//	advect(2, Vy, Vy0, Vx0, Vy0, Vz0);
//	advect(3, Vz, Vz0, Vx0, Vy0, Vz0);
//
//	project(V,V0);
//
//	diffuse(0, s, density, diffusion);
//	advect(0, density, s, Vx, Vy, Vz);
//}
//
//void wind_field::addForces(int x, int y, int z, vec3 amount)
//{
//	int N = N;
//	V[IX(x, y, z)] += amount;
//}
//
//wind_field::wind_field(int N, int diffusion, int viscosity, float dt)
//{
//
//	this->N = N;
//	this->dt = dt;
//	this->diffusion = diffusion;
//	this->viscosity = viscosity;
//
//	this->s = new float[N * N * N];
//	this->density = new float[N * N * N];
//
//	this->V = new vec3[N * N * N];
//	this->V0 = new vec3[N * N * N];
//}
//
//void wind_field::diffuse(int b, vec3* x, vec3* x0, float diff)
//{
//	float a = dt * diff * (N - 2) * (N - 2);
//	lin_solve(b, x, x0, a, 1 + 6 * a);
//}
//
//void wind_field::diffuse(int b, float* x, float* x0, float diff)
//{
//	float a = dt * diff * (N - 2) * (N - 2);
//	lin_solve(b, x, x0, a, 1 + 6 * a);
//}
//
//void wind_field::lin_solve(int b, vec3* x, vec3* x0, float a, float c) {
//	float cRecip = 1.0 / c;
//	for (int k = 0; k < iter; ++k) {
//		for (int m = 1; m < N - 1; ++m) {
//			for (int j = 1; j < N - 1; ++j) {
//				for (int i = 1; i < N - 1; ++i) {
//
//					x[IX(i, j, m)].x =
//						(x0[IX(i, j, m)].x
//							+ a * (x[IX(i + 1, j, m)].x
//								+ x[IX(i - 1, j, m)].x
//								+ x[IX(i, j + 1, m)].x
//								+ x[IX(i, j - 1, m)].x
//								+ x[IX(i, j, m + 1)].x
//								+ x[IX(i, j, m - 1)].x
//								)) * cRecip;
//
//					x[IX(i, j, m)].y =
//						(x0[IX(i, j, m)].y
//							+ a * (x[IX(i + 1, j, m)].y
//								+ x[IX(i - 1, j, m)].y
//								+ x[IX(i, j + 1, m)].y
//								+ x[IX(i, j - 1, m)].y
//								+ x[IX(i, j, m + 1)].y
//								+ x[IX(i, j, m - 1)].y
//								)) * cRecip;
//
//					x[IX(i, j, m)].z =
//						(x0[IX(i, j, m)].z
//							+ a * (x[IX(i + 1, j, m)].z
//								+ x[IX(i - 1, j, m)].z
//								+ x[IX(i, j + 1, m)].z
//								+ x[IX(i, j - 1, m)].z
//								+ x[IX(i, j, m + 1)].z
//								+ x[IX(i, j, m - 1)].z
//								)) * cRecip;
//				}
//			}
//		}
//		//set_bnd(b, x);
//	}
//}
//
//void wind_field::lin_solve(int b, float* x, float* x0, float a, float c) {
//	float cRecip = 1.0 / c;
//	for (int k = 0; k < iter; ++k) {
//		for (int m = 1; m < N - 1; ++m) {
//			for (int j = 1; j < N - 1; ++j) {
//				for (int i = 1; i < N - 1; ++i) {
//					x[IX(i, j, m)] =
//						(x0[IX(i, j, m)]
//							+ a * (x[IX(i + 1, j, m)]
//								+ x[IX(i - 1, j, m)]
//								+ x[IX(i, j + 1, m)]
//								+ x[IX(i, j - 1, m)]
//								+ x[IX(i, j, m + 1)]
//								+ x[IX(i, j, m - 1)]
//								)) * cRecip;
//				}
//			}
//		}
//		set_bnd(b, x);
//	}
//}
//void wind_field::project(vec3* veloc, vec3* pdiv) {
//	for (int k = 1; k < N - 1; ++k) {
//		for (int j = 1; j < N - 1; ++j) {
//			for (int i = 1; i < N - 1; ++i) {
//				pdiv[IX(i, j, k)].y = -0.5f * (
//					veloc[IX(i + 1, j, k)].x
//					- veloc[IX(i - 1, j, k)].x
//					+ veloc[IX(i, j + 1, k)].y
//					- veloc[IX(i, j - 1, k)].y
//					+ veloc[IX(i, j, k + 1)].z
//					- veloc[IX(i, j, k - 1)].z
//					) / N;
//				pdiv[IX(i, j, k)].x = 0;
//			}
//		}
//	}
//	//set_bnd(0, div);
//	//set_bnd(0, p);
//	lin_solve(0, p, div, 1, 6);
//
//	for (int k = 1; k < N - 1; ++k) {
//		for (int j = 1; j < N - 1; ++j) {
//			for (int i = 1; i < N - 1; ++i) {
//				veloc[IX(i, j, k)].x -= 0.5f * (pdiv[IX(i + 1, j, k)].x
//					- pdiv[IX(i - 1, j, k)].x) * N;
//				veloc[IX(i, j, k)].y -= 0.5f * (pdiv[IX(i, j + 1, k)].x
//					- pdiv[IX(i, j - 1, k)].x) * N;
//				veloc[IX(i, j, k)].z -= 0.5f * (pdiv[IX(i, j, k + 1)].x
//					- pdiv[IX(i, j, k - 1)].x) * N;
//			}
//		}
//	}
//	//set_bnd(1, velocX);
//	//set_bnd(2, velocY);
//	//set_bnd(3, velocZ);
//}
//
//void wind_field::advect(int b, float* d, float* d0, vec3* veloc)
//{
//	vec3 ijk0, ijk1;
//
//	float ddt = dt * (N - 2);
//
//	vec3 stu0, stu1;
//
//	vec3 tmp, dim;
//
//	float Nfloat = N;
//	float ifloat, jfloat, kfloat;
//	int i, j, k;
//
//	for (k = 1, kfloat = 1; k < N - 1; ++k, kfloat++) {
//		for (j = 1, jfloat = 1; j < N - 1; ++j, jfloat++) {
//			for (i = 1, ifloat = 1; i < N - 1; ++i, ifloat++) {
//				tmp = ddt * veloc[IX(i, j, k)];
//				dim = vec3(ifloat - tmp.x, jfloat - tmp.y, kfloat - tmp.z);
//
//				if (dim.x < 0.5f)dim.x = 0.5f;
//				if (dim.x > Nfloat + 0.5f) dim.x = Nfloat + 0.5f;
//
//				if (dim.y < 0.5f) dim.y = 0.5f;
//				if (dim.y > Nfloat + 0.5f) dim.y = Nfloat + 0.5f;
//
//				if (dim.z < 0.5f) dim.z = 0.5f;
//				if (dim.z > Nfloat + 0.5f) dim.z = Nfloat + 0.5f;
//
//				ijk0 = vec3(floorf(dim.x), floorf(dim.y), floorf(dim.z));
//				ijk1 = vec3(ijk0.x + 1.0f, ijk0.y + 1.0f, ijk0.z + 1.0f);
//
//				stu1 = vec3(dim.x - ijk0.x, dim.y - ijk0.y, dim.z - ijk0.z);
//				stu0 = vec3(1.0f - stu1.x, 1.0f - stu1.y, 1.0f - stu1.z);
//
//				int i0i = ijk0.x;
//				int i1i = ijk1.x;
//				int j0i = ijk0.y;
//				int j1i = ijk1.y;
//				int k0i = ijk0.z;
//				int k1i = ijk1.z;
//
//				d[IX(i, j, k)] =
//
//					stu0.x * (stu0.y * (stu0.z * d0[IX(i0i, j0i, k0i)]
//										+ stu1.z * d0[IX(i0i, j0i, k1i)])
//						+ (stu1.y * (stu0.z * d0[IX(i0i, j1i, k0i)]
//										+ stu1.z * d0[IX(i0i, j1i, k1i)])))
//					+ stu1.x * (stu0.y * (stu0.z * d0[IX(i1i, j0i, k0i)]
//										+ stu1.z * d0[IX(i1i, j0i, k1i)])
//						+ (stu1.y * (stu0.z * d0[IX(i1i, j1i, k0i)]
//										+ stu1.z * d0[IX(i1i, j1i, k1i)])));
//			}
//		}
//	}
//	//set_bnd(b, d);
//}
//
//void wind_field::set_bnd(int b, vec3* x)
//{
//	for (int j = 1; j < N - 1; ++j) {
//		for (int i = 1; i < N - 1; ++i) {
//			x[IX(i, j, 0)] = b == 3 ? -x[IX(i, j, 1)] : x[IX(i, j, 1)];
//			x[IX(i, j, N - 1)] = b == 3 ? -x[IX(i, j, N - 2)] : x[IX(i, j, N - 2)];
//		}
//	}
//	for (int k = 1; k < N - 1; ++k) {
//		for (int i = 1; i < N - 1; ++i) {
//			x[IX(i, 0, k)] = b == 2 ? -x[IX(i, 1, k)] : x[IX(i, 1, k)];
//			x[IX(i, N - 1, k)] = b == 2 ? -x[IX(i, N - 2, k)] : x[IX(i, N - 2, k)];
//		}
//	}
//	for (int k = 1; k < N - 1; ++k) {
//		for (int j = 1; j < N - 1; ++j) {
//			x[IX(0, j, k)] = b == 1 ? -x[IX(1, j, k)] : x[IX(1, j, k)];
//			x[IX(N - 1, j, k)] = b == 1 ? -x[IX(N - 2, j, k)] : x[IX(N - 2, j, k)];
//		}
//	}
//
//	x[IX(0, 0, 0)] = 0.33f * (x[IX(1, 0, 0)]
//		+ x[IX(0, 1, 0)]
//		+ x[IX(0, 0, 1)]);
//	x[IX(0, N - 1, 0)] = 0.33f * (x[IX(1, N - 1, 0)]
//		+ x[IX(0, N - 2, 0)]
//		+ x[IX(0, N - 1, 1)]);
//	x[IX(0, 0, N - 1)] = 0.33f * (x[IX(1, 0, N - 1)]
//		+ x[IX(0, 1, N - 1)]
//		+ x[IX(0, 0, N)]);
//	x[IX(0, N - 1, N - 1)] = 0.33f * (x[IX(1, N - 1, N - 1)]
//		+ x[IX(0, N - 2, N - 1)]
//		+ x[IX(0, N - 1, N - 2)]);
//	x[IX(N - 1, 0, 0)] = 0.33f * (x[IX(N - 2, 0, 0)]
//		+ x[IX(N - 1, 1, 0)]
//		+ x[IX(N - 1, 0, 1)]);
//	x[IX(N - 1, N - 1, 0)] = 0.33f * (x[IX(N - 2, N - 1, 0)]
//		+ x[IX(N - 1, N - 2, 0)]
//		+ x[IX(N - 1, N - 1, 1)]);
//	x[IX(N - 1, 0, N - 1)] = 0.33f * (x[IX(N - 2, 0, N - 1)]
//		+ x[IX(N - 1, 1, N - 1)]
//		+ x[IX(N - 1, 0, N - 2)]);
//	x[IX(N - 1, N - 1, N - 1)] = 0.33f * (x[IX(N - 2, N - 1, N - 1)]
//		+ x[IX(N - 1, N - 2, N - 1)]
//		+ x[IX(N - 1, N - 1, N - 2)]);
//}
//
//void wind_field::set_bnd(int b, float* x)
//{
//	for (int j = 1; j < N - 1; ++j) {
//		for (int i = 1; i < N - 1; ++i) {
//			x[IX(i, j, 0)] = b == 3 ? -x[IX(i, j, 1)] : x[IX(i, j, 1)];
//			x[IX(i, j, N - 1)] = b == 3 ? -x[IX(i, j, N - 2)] : x[IX(i, j, N - 2)];
//		}
//	}
//	for (int k = 1; k < N - 1; ++k) {
//		for (int i = 1; i < N - 1; ++i) {
//			x[IX(i, 0, k)] = b == 2 ? -x[IX(i, 1, k)] : x[IX(i, 1, k)];
//			x[IX(i, N - 1, k)] = b == 2 ? -x[IX(i, N - 2, k)] : x[IX(i, N - 2, k)];
//		}
//	}
//	for (int k = 1; k < N - 1; ++k) {
//		for (int j = 1; j < N - 1; ++j) {
//			x[IX(0, j, k)] = b == 1 ? -x[IX(1, j, k)] : x[IX(1, j, k)];
//			x[IX(N - 1, j, k)] = b == 1 ? -x[IX(N - 2, j, k)] : x[IX(N - 2, j, k)];
//		}
//	}
//
//	x[IX(0, 0, 0)] = 0.33f * (x[IX(1, 0, 0)]
//		+ x[IX(0, 1, 0)]
//		+ x[IX(0, 0, 1)]);
//	x[IX(0, N - 1, 0)] = 0.33f * (x[IX(1, N - 1, 0)]
//		+ x[IX(0, N - 2, 0)]
//		+ x[IX(0, N - 1, 1)]);
//	x[IX(0, 0, N - 1)] = 0.33f * (x[IX(1, 0, N - 1)]
//		+ x[IX(0, 1, N - 1)]
//		+ x[IX(0, 0, N)]);
//	x[IX(0, N - 1, N - 1)] = 0.33f * (x[IX(1, N - 1, N - 1)]
//		+ x[IX(0, N - 2, N - 1)]
//		+ x[IX(0, N - 1, N - 2)]);
//	x[IX(N - 1, 0, 0)] = 0.33f * (x[IX(N - 2, 0, 0)]
//		+ x[IX(N - 1, 1, 0)]
//		+ x[IX(N - 1, 0, 1)]);
//	x[IX(N - 1, N - 1, 0)] = 0.33f * (x[IX(N - 2, N - 1, 0)]
//		+ x[IX(N - 1, N - 2, 0)]
//		+ x[IX(N - 1, N - 1, 1)]);
//	x[IX(N - 1, 0, N - 1)] = 0.33f * (x[IX(N - 2, 0, N - 1)]
//		+ x[IX(N - 1, 1, N - 1)]
//		+ x[IX(N - 1, 0, N - 2)]);
//	x[IX(N - 1, N - 1, N - 1)] = 0.33f * (x[IX(N - 2, N - 1, N - 1)]
//		+ x[IX(N - 1, N - 2, N - 1)]
//		+ x[IX(N - 1, N - 1, N - 2)]);
//}

	// ___________________________________		ITEMS __________________________________________________//
	//vec = vec3(N, N, N);// setup
	//mv = view;

	//mv = translate(mv, vec3(-N / 2, 0, -N / 2));	 // translate
	//mv = mv * orientation(normalize(vec), vec3(0, 0, 1));		// rotation
	//mv = scale(mv, vec3(0.1, 0.1, length(vec)));				// scale

	//glUniformMatrix4fv(glGetUniformLocation(shader, "uModelViewMatrix"), 1, false, value_ptr(mv));
	//drawCylinder();

	// ___________________________________		ITEMS idk __________________________________________________//
	//float x = w_field->Vy[IX(0, 0, 0)];
	//cout << x;


	//_________________________		Wind  testing      _________________________________________________//
	//vec = vec3(3, 5, 5); // setup
	//mv = view;

	//mv = mv * orientation(normalize(vec), vec3(0, 0, 1));	// direction
	//mv = scale(mv, vec3(1, 1, length(vec)));				 // strength

	//glUniformMatrix4fv(glGetUniformLocation(shader, "uModelViewMatrix"), 1, false, value_ptr(mv));
	//drawCone();