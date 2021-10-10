#pragma once

// ERROR: UNUSED

// std
#include <iostream>
#include <string>
#include <chrono>
#include <vector>

// glm
#include <glm/glm.hpp>
//#include <glm/gtc/type_ptr.hpp>

// project
//#include "opengl.hpp"
//#include "cgra/cgra_mesh.hpp"

/*struct basic_model {
	GLuint shader = 0;
	cgra::gl_mesh mesh;
	glm::vec3 color{ 0.7 };
	glm::mat4 modelTransform{ 1.0 };
	GLuint texture;

	void draw(const glm::mat4& view, const glm::mat4 proj);
};*/

// TODO: make grass a namespace?

/*struct grass_blade {
	std::vector<glm::vec3> control_points;
	cgra::gl_mesh curve;
	cgra::gl_mesh spline;
};*/

class grass_model {
private:
	// four control points
	glm::vec3 control_point0 = glm::vec3(0, 0, 0);
	glm::vec3 control_point1 = glm::vec3(1, 1, 0);
	glm::vec3 control_point2 = glm::vec3(0, 2, 0);
	glm::vec3 control_point3 = glm::vec3(0, 3, 0);

	// meshes
	cgra::gl_mesh curve; // FIXME: missing specifier
	cgra::gl_mesh spline; // FIXME: missing specifier
public:
	grass_model(glm::vec3 p0, glm::vec3 p1, glm::vec3 p2, glm::vec3 p3); // constructor
	grass_model() = default;

	// getters and setters
	glm::vec3 * getControlPoints();
	void setControlPoints(glm::vec3 p0, glm::vec3 p1, glm::vec3 p2, glm::vec3 p3);
	cgra::gl_mesh getSpline();
	cgra::gl_mesh getCurve();
	void setSpline();
	void setCurve();

	// drawing and rendering
	void drawPolyline(const glm::mat4& view, const glm::mat4 proj);
	void drawCurve(const glm::mat4& view, const glm::mat4 proj);
	glm::vec3 interpolate(float t);
};