// ERROR: UNUSED

#include "grass_model.hpp"
#include "grass_model.hpp"

// std
#include <iostream>
#include <string>
#include <chrono>

// project
#include "opengl.hpp"
#include "grass_model.hpp"
#include "cgra_mesh.hpp"

using namespace std;
using namespace cgra;
using namespace glm;

grass_model::grass_model(glm::vec3 p0, glm::vec3 p1, glm::vec3 p2, glm::vec3 p3)
{
}

/*grass_model::grass_model()
{
}*/

vec3* grass_model::getControlPoints() {
	// return array containing the four control points
	return { control_point0, control_point1, control_point2, control_point3 };
}
void grass_model::setControlPoints(vec3 p0, vec3 p1, vec3 p2, vec3 p3) {
	control_point0 = p0;
	control_point1 = p1;
	control_point2 = p2;
	control_point3 = p3;
}

gl_mesh grass_model::getSpline() {
	return spline;
}
gl_mesh grass_model::getCurve() {
	return curve;
}

void grass_model::setSpline() {
	// set spline
	mesh_builder mb = mesh_builder();
	vec3 cp[] = getControlPoints();
	for (int i = 0; i < 4; i++) {
		mesh_vertex mv;
		mv.pos = cp[i];
		mb.push_vertex(mv);
		mb.push_index(i);
	}
	spline = mb.build();
}

void grass_model::setCurve() {

}

// drawing and rendering
void grass_model::drawPolyline(const glm::mat4& view, const glm::mat4 proj) {
	//spline.color = vec3(1, 0, 0);
	//spline.mesh.mode = GL_LINE_STRIP;

	/*mat4 modelview = view;

	glUseProgram(shader); // load shader and variables
	glUniformMatrix4fv(glGetUniformLocation(shader, "uProjectionMatrix"), 1, false, value_ptr(proj));
	glUniformMatrix4fv(glGetUniformLocation(shader, "uModelViewMatrix"), 1, false, value_ptr(modelview));
	glUniform3fv(glGetUniformLocation(shader, "uColor"), 1, value_ptr(color));*/

	// TODO: draw spline
	spline.draw();
}

void grass_model::drawCurve(const glm::mat4& view, const glm::mat4 proj) {
	// TODO: draw control points
	// TODO: draw straight lines between control points (in red?)
}

glm::vec3 grass_model::interpolate(float t) {
	// 0 <= t <= 1
	return pow((1.0f - t), 3.0f) * control_point0
		+ 3.0f * t * pow((1.0f - t), 2.0f) * control_point1
		+ 3.0f * pow(t, 2.0f) * (1.0f - t) * control_point2
		+ pow(t, 3.0f) * control_point3; // FIXME: cant use * for vec3
}