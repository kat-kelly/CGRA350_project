
// std
#include <iostream>
#include <string>
#include <chrono>

// glm
#include <glm/gtc/constants.hpp>
#include <glm/gtc/matrix_transform.hpp>

// project
#include "application.hpp"
#include "cgra/cgra_geometry.hpp"
#include "cgra/cgra_gui.hpp"
#include "cgra/cgra_image.hpp"
#include "cgra/cgra_shader.hpp"
#include "cgra/cgra_wavefront.hpp"
//#include "cgra/grass_model.hpp"


using namespace std;
using namespace cgra;
using namespace glm;



void basic_model::draw(const glm::mat4 &view, const glm::mat4 proj) { // UNUSED
	mat4 modelview = view * modelTransform;

	glUseProgram(shader); // load shader and variables
	glUniformMatrix4fv(glGetUniformLocation(shader, "uProjectionMatrix"), 1, false, value_ptr(proj));
	glUniformMatrix4fv(glGetUniformLocation(shader, "uModelViewMatrix"), 1, false, value_ptr(modelview));
	glUniform3fv(glGetUniformLocation(shader, "uColor"), 1, value_ptr(color));
	mesh.draw(); // draw
}

// grass model methods
void grass_model::drawSpline(const glm::mat4& view, const glm::mat4 proj) {
	mat4 modelview = view * modelTransform;

	glUseProgram(shader); // load shader and variables
	glUniformMatrix4fv(glGetUniformLocation(shader, "uProjectionMatrix"), 1, false, value_ptr(proj));
	glUniformMatrix4fv(glGetUniformLocation(shader, "uModelViewMatrix"), 1, false, value_ptr(modelview));
	glUniform3fv(glGetUniformLocation(shader, "uColor"), 1, value_ptr(vec3(1,0,0))); // color is red

	spline_mesh.draw(); // draw
}

void grass_model::drawCurve(const glm::mat4& view, const glm::mat4 proj) {
	mat4 modelview = view * modelTransform;

	glUseProgram(shader); // load shader and variables
	glUniformMatrix4fv(glGetUniformLocation(shader, "uProjectionMatrix"), 1, false, value_ptr(proj));
	glUniformMatrix4fv(glGetUniformLocation(shader, "uModelViewMatrix"), 1, false, value_ptr(modelview));
	glUniform3fv(glGetUniformLocation(shader, "uColor"), 1, value_ptr(vec3(0,1,0))); // color is green

	curve_mesh.draw(); // draw
}

void grass_model::setMeshes(GLuint shad) {
	// spline mesh
	mesh_builder mb_spline = mesh_builder();
	for (int i = 0; i < 4; i++) {
		mesh_vertex mv;
		mv.pos = controlPts[i];
		mv.norm = vec3(1, 0, 0);
		mb_spline.push_vertex(mv);
		mb_spline.push_index(i);
	}
	mb_spline.mode = GL_LINE_STRIP;
	spline_mesh = mb_spline.build();

	// curve mesh
	mesh_builder mb_curve = mesh_builder();
	int count = 0;
	for (float i = 0; i <= 1; i = i + 0.2) {
		mesh_vertex mv;
		mv.pos = interpolateBezier(i);
		mv.norm = vec3(1, 0, 0);
		mb_curve.push_vertex(mv);
		mb_curve.push_index(count);
		count++;
	}
	mb_curve.mode = GL_LINE_STRIP;
	curve_mesh = mb_curve.build();
	shader = shad;
}

void grass_model::setControlPts(vec3 cp[4]) {
	*controlPts = *cp;
}

vec3 grass_model::interpolateBezier(float t) {
	return pow((1.0f - t), 3.0f) * controlPts[0]
		+ 3.0f * t * pow((1.0f - t), 2.0f) * controlPts[1]
		+ 3.0f * pow(t, 2.0f) * (1.0f - t) * controlPts[2]
		+ pow(t, 3.0f) * controlPts[3];
}

Application::Application(GLFWwindow* window) : m_window(window) {

	shader_builder sb;
	sb.set_shader(GL_VERTEX_SHADER, CGRA_SRCDIR + std::string("//res//shaders//color_vert.glsl"));
	sb.set_shader(GL_FRAGMENT_SHADER, CGRA_SRCDIR + std::string("//res//shaders//color_frag.glsl"));
	GLuint shader = sb.build();

	// set up grass blade
	vec3 arr[4] = { vec3(0, 0, 0), vec3(1, 1, 0), vec3(1, 2, 0), vec3(0, 3, 0) }; // new control points
	grass.setControlPts(arr);
	grass.setMeshes(shader);
}


void Application::render() {

	// retrieve the window hieght
	int width, height;
	glfwGetFramebufferSize(m_window, &width, &height);

	m_windowsize = vec2(width, height); // update window size
	glViewport(0, 0, width, height); // set the viewport to draw to the entire window

	// clear the back-buffer
	glClearColor(0.3f, 0.3f, 0.4f, 1.0f);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	// enable flags for normal/forward rendering
	glEnable(GL_DEPTH_TEST);
	glDepthFunc(GL_LESS);

	// projection matrix
	mat4 proj = perspective(1.f, float(width) / height, 0.1f, 1000.f);

	// view matrix
	mat4 view = translate(mat4(1), vec3(0, 0, -m_distance))
		* rotate(mat4(1), m_pitch, vec3(1, 0, 0))
		* rotate(mat4(1), m_yaw, vec3(0, 1, 0));


	// helpful draw options
	if (m_show_grid) drawGrid(view, proj);
	if (m_show_axis) drawAxis(view, proj);
	glPolygonMode(GL_FRONT_AND_BACK, (m_showWireframe) ? GL_LINE : GL_FILL);
  
  // wind
  w_model.run(view, proj);

	// draw grass blade
	grass.drawSpline(view, proj);
	grass.drawCurve(view, proj);
}

basic_model Application::meshToModel(gl_mesh mesh) { // UNUSED
	basic_model model;
	model.mesh = mesh;
	model.color = vec3(1, 0, 0); // TODO: make changable
	return model;
}

// added methods
basic_model Application::toRenderMesh(grass_model grass_blade) { // UNUSED
	// TODO: determine level of detail for this grass blade
	// TODO: construct mesh from line segments
	basic_model grass_bm;
	//grass.mesh = grass_blade.
	//grass.color = vec3(0, 1, 0);
	return grass_bm;
}


void Application::renderGUI() {

	// setup window
	ImGui::SetNextWindowPos(ImVec2(5, 5), ImGuiSetCond_Once);
	ImGui::SetNextWindowSize(ImVec2(300, 200), ImGuiSetCond_Once);
	ImGui::Begin("Options", 0);

	// display current camera parameters
	ImGui::Text("Application %.3f ms/frame (%.1f FPS)", 1000.0f / ImGui::GetIO().Framerate, ImGui::GetIO().Framerate);
	ImGui::SliderFloat("Pitch", &m_pitch, -pi<float>() / 2, pi<float>() / 2, "%.2f");
	ImGui::SliderFloat("Yaw", &m_yaw, -pi<float>(), pi<float>(), "%.2f");
	ImGui::SliderFloat("Distance", &m_distance, 0, 100, "%.2f", 2.0f);

	// helpful drawing options
	ImGui::Checkbox("Show axis", &m_show_axis);
	ImGui::SameLine();
	ImGui::Checkbox("Show grid", &m_show_grid);
	ImGui::Checkbox("Wireframe", &m_showWireframe);
	ImGui::SameLine();
	if (ImGui::Button("Screenshot")) rgba_image::screenshot(true);

	ImGui::Separator();
	ImGui::Text("Wind settings");
	if (ImGui::RadioButton("Panic button", w_model.display == -1)) {
		w_model.display = -1;
		w_model = wind_model();
	}
	if (ImGui::RadioButton("Off", w_model.display == 0)) w_model.display = 0;	ImGui::SameLine();
	if (ImGui::RadioButton("On", w_model.display == 1)) w_model.display = 1;	ImGui::SameLine();
	if (ImGui::RadioButton("Visualize", w_model.display == 2)) w_model.display = 2;
	if (w_model.display > 0) {
		ImGui::SliderFloat("Wind Strength", &w_model.w_strength, 0, 5, "%.1f");
		ImGui::SliderFloat("Wind Yaw", &w_model.w_angle, -20, 20, "%.1f");
		ImGui::SliderFloat("Wind Pulse", &w_model.pulse, 0.0, 0.05, "%.2f");
	}
	ImGui::Separator();


	// finish creating window
	ImGui::End();
}


void Application::cursorPosCallback(double xpos, double ypos) {
	if (m_leftMouseDown) {
		vec2 whsize = m_windowsize / 2.0f;

		// clamp the pitch to [-pi/2, pi/2]
		m_pitch += float(acos(glm::clamp((m_mousePosition.y - whsize.y) / whsize.y, -1.0f, 1.0f))
			- acos(glm::clamp((float(ypos) - whsize.y) / whsize.y, -1.0f, 1.0f)));
		m_pitch = float(glm::clamp(m_pitch, -pi<float>() / 2, pi<float>() / 2));

		// wrap the yaw to [-pi, pi]
		m_yaw += float(acos(glm::clamp((m_mousePosition.x - whsize.x) / whsize.x, -1.0f, 1.0f))
			- acos(glm::clamp((float(xpos) - whsize.x) / whsize.x, -1.0f, 1.0f)));
		if (m_yaw > pi<float>()) m_yaw -= float(2 * pi<float>());
		else if (m_yaw < -pi<float>()) m_yaw += float(2 * pi<float>());
	}

	// updated mouse position
	m_mousePosition = vec2(xpos, ypos);
}


void Application::mouseButtonCallback(int button, int action, int mods) {
	(void)mods; // currently un-used

	// capture is left-mouse down
	if (button == GLFW_MOUSE_BUTTON_LEFT)
		m_leftMouseDown = (action == GLFW_PRESS); // only other option is GLFW_RELEASE
}


void Application::scrollCallback(double xoffset, double yoffset) {
	(void)xoffset; // currently un-used
	m_distance *= pow(1.1f, -yoffset);
}


void Application::keyCallback(int key, int scancode, int action, int mods) {
	(void)key, (void)scancode, (void)action, (void)mods; // currently un-used
}


void Application::charCallback(unsigned int c) {
	(void)c; // currently un-used
}
