
#pragma once

// glm
#include <glm/glm.hpp>
#include <glm/gtc/type_ptr.hpp>

// project
#include "opengl.hpp"
#include "cgra/cgra_mesh.hpp"
#include "skeleton_model.hpp"
#include "wind_model.hpp"
//#include "cgra/grass_model.hpp"


// Basic model that holds the shader, mesh and transform for drawing.
// Can be copied and modified for adding in extra information for drawing
// including textures for texture mapping etc.
struct basic_model {
	GLuint shader = 0;
	cgra::gl_mesh mesh;
	glm::vec3 color{0.7};
	glm::mat4 modelTransform{1.0};
	GLuint texture;

	void draw(const glm::mat4 &view, const glm::mat4 proj);
};

struct grass_model {
	GLuint shader = 0;
	cgra::gl_mesh spline_mesh;
	cgra::gl_mesh curve_mesh;
	glm::vec3 color{ 0.7 };
	glm::mat4 modelTransform{ 1.0 };
	glm::vec3 controlPts[4]{ glm::vec3(0, 0, 0), glm::vec3(1, 1, 0), glm::vec3(1, 2, 0), glm::vec3(0, 3, 0) }; // FIXME: make changable

	void drawSpline(const glm::mat4& view, const glm::mat4 proj);
	void drawCurve(const glm::mat4& view, const glm::mat4 proj);
	void setControlPts(glm::vec3 cp[4]);
	void setMeshes(GLuint shader);
	glm::vec3 interpolateBezier(float t);
};

// Main application class
//
class Application {
private:
	// window
	glm::vec2 m_windowsize = glm::vec2(10,10);
	GLFWwindow *m_window;

	// oribital camera
	float m_pitch = .86;
	float m_yaw = -.86;
	float m_distance = 20;

	// last input
	bool m_leftMouseDown = false;
	glm::vec2 m_mousePosition = glm::vec2(0,0);

	// drawing flags
	bool m_show_axis = true;
	bool m_show_grid = true;
	bool m_showWireframe = false;

	// geometry
	basic_model m_model;
	wind_model w_model;
	//basic_model m_model;
	grass_model grass;// = grass_model(glm::vec3(0, 0, 0), glm::vec3(1, 1, 0), glm::vec3(0, 2, 0), glm::vec3(0, 3, 0));

public:
	// setup
	Application(GLFWwindow *);

	// disable copy constructors (for safety)
	Application(const Application&) = delete;
	Application& operator=(const Application&) = delete;

	// rendering callbacks (every frame)
	void render();
	void renderGUI();

	// input callbacks
	void cursorPosCallback(double xpos, double ypos);
	void mouseButtonCallback(int button, int action, int mods);
	void scrollCallback(double xoffset, double yoffset);
	void keyCallback(int key, int scancode, int action, int mods);
	void charCallback(unsigned int c);

	// convert grass model into basic model
	basic_model meshToModel(cgra::gl_mesh mesh);
	basic_model toRenderMesh(grass_model blade);
};