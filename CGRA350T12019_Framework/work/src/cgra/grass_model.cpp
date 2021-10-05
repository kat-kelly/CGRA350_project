
// std
#include <iostream>
#include <string>
#include <chrono>

// project
#include <opengl.hpp>
#include "grass_model.hpp"
#include "cgra_mesh.hpp"

using namespace std;
using namespace cgra;
using namespace glm;


grass_model::grass_model(glm::vec3 p0, glm::vec3 p1, glm::vec3 p2, glm::vec3 p3) {

}

vec3* grass_model::getControlPoints() {
	// return array containing the four control points
	return { control_point0, control_point1, control_point2, control_point3 };
}
void grass_model::setControlPoints(glm::vec3 p0, glm::vec3 p1, glm::vec3 p2, glm::vec3 p3) {
	control_point0 = p0;
	control_point1 = p1;
	control_point2 = p2;
	control_point3 = p3;
}

// drawing and rendering
void grass_model::drawPolyline() {
	// set spline
	mesh_builder mb = mesh_builder();
	vec3[] cp = getControlPoints();
	for (int i = 0; i < 4; i++) {
		mesh_vertex mv;
		mv.pos = cp[i];
		mb.push_vertex(mv);
		mb.push_index(i);
	}
	spline.mesh = mb.build();
	spline.color = vec3(1, 0, 0);
	spline.mesh.mode = GL_LINE_STRIP;

	// TODO: draw spline
	spline.draw();
}

void grass_model::drawCurve() {
	// TODO: draw control points
	// TODO: draw straight lines between control points (in red?)
}

glm::vec3 grass_model::interpolate(float t) {
	// 0 <= t <= 1
	return pow((1 - t), 3) * control_point0
		+ 3 * t * pow((1 - t), 2) * control_point1
		+ 3 * pow(t, 2) * (1 - t) * control_point2
		+ pow(t, 3) * control_point3
}