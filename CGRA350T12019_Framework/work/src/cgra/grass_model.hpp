#pragma once

// glm
#include <glm/glm.hpp>
#include <glm/gtc/type_ptr.hpp>

// project
#include <opengl.hpp>
#include "application.hpp"

class grass_model {
private:
	// four control points
	glm::vec3 control_point0;
	glm::vec3 control_point1;
	glm::vec3 control_point2;
	glm::vec3 control_point3;

	// meshes
	basic_model curve; // FIXME: missing specifier
	basic_model spline; // FIXME: missing specifier
public:
	grass_model(glm::vec3 p0, glm::vec3 p1, glm::vec3 p2, glm::vec3 p3); // constructor

	// getters and setters
	glm::vec3 * getControlPoints();
	void setControlPoints(glm::vec3 p0, glm::vec3 p1, glm::vec3 p2, glm::vec3 p3);

	// drawing and rendering
	void drawPolyline();
	void drawCurve();
	glm::vec3 interpolate(float t);
};