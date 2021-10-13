
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

// grass model methods by Katrina
void grass_model::drawSpline(const glm::mat4& view, const glm::mat4 proj) {
	mat4 modelview = view * modelTransform;

	glUseProgram(shader); // load shader and variables
	glUniformMatrix4fv(glGetUniformLocation(shader, "uProjectionMatrix"), 1, false, value_ptr(proj));
	glUniformMatrix4fv(glGetUniformLocation(shader, "uModelViewMatrix"), 1, false, value_ptr(modelview));
	glUniform3fv(glGetUniformLocation(shader, "uColor"), 1, value_ptr(vec3(1,0,0))); // color is red

	spline_mesh.draw(); // draw
}

void grass_model::drawCurve(const glm::mat4& view, const glm::mat4 proj) { // TODO: get LOD from distance
	//GLuint sh;

	//set_shader(GL_TESS_CONTROL_SHADER, CGRA_SRCDIR + std::string("//res//shaders//tessellation_control.glsl"));

	mat4 modelview = view * modelTransform;

	glUseProgram(shader); // load shader and variables
	glUniformMatrix4fv(glGetUniformLocation(shader, "uProjectionMatrix"), 1, false, value_ptr(proj));
	glUniformMatrix4fv(glGetUniformLocation(shader, "uModelViewMatrix"), 1, false, value_ptr(modelview));
	glUniform3fv(glGetUniformLocation(shader, "uColor"), 1, value_ptr(vec3(0,1,0))); // color is green

	curve_mesh.draw(); // draw
}

void grass_model::drawBlade(const glm::mat4& view, const glm::mat4 proj, glm::vec3 camera) {
	mat4 modelview = view * modelTransform;

	glUseProgram(shader); // load shader and variables
	glUniformMatrix4fv(glGetUniformLocation(shader, "uProjectionMatrix"), 1, false, value_ptr(proj));
	glUniformMatrix4fv(glGetUniformLocation(shader, "uModelViewMatrix"), 1, false, value_ptr(modelview));
	glUniform3fv(glGetUniformLocation(shader, "uColor"), 1, value_ptr(vec3(0, 1, 0))); // color is green
	// shaders added by Katrina
	glUniform3fv(glGetUniformLocation(shader, "cameraPos"), 1, value_ptr(camera)); // tessellation control

	spline_mesh.draw(); // error
}

void grass_model::setMeshes(GLuint shad) {
	// spline mesh
	mesh_builder mb_spline = mesh_builder();
	vec3 parentTrans = vec3(0); // heirarchical
	for (int i = 0; i < 4; i++) {
		mesh_vertex mv;
		mv.pos = parentTrans + controlPts[i];
		mv.norm = vec3(1, 0, 0);
		mb_spline.push_vertex(mv);
		mb_spline.push_index(i);
		parentTrans = parentTrans + controlPts[i];
	}
	mb_spline.mode = GL_LINE_STRIP;
	spline_mesh = mb_spline.build();

	// curve mesh
	mesh_builder mb_curve = mesh_builder();
	int count = 0;
	for (float i = 0; i <= 1; i = i + (1/lod)) {
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
	// heirarchical points
	vec3 p0 = controlPts[0];
	vec3 p1 = p0 + controlPts[1];
	vec3 p2 = p1 + controlPts[2];
	vec3 p3 = p2 + controlPts[3];

	// get point along curve
	return pow((1.0f - t), 3.0f) * p0
		+ 3.0f * t * pow((1.0f - t), 2.0f) * p1
		+ 3.0f * pow(t, 2.0f) * (1.0f - t) * p2
		+ pow(t, 3.0f) * p3;
}

Application::Application(GLFWwindow* window) : m_window(window) {
	setGrass();
}

void Application::setGrass()
{
	// set shaders
	shader_builder sb;
	sb.set_shader(GL_VERTEX_SHADER, CGRA_SRCDIR + std::string("//res//shaders//color_vert.glsl"));
	if (show_rendered_grass) {
		// added shaders by Katrina
		sb.set_shader(GL_TESS_CONTROL_SHADER, CGRA_SRCDIR + std::string("//res//shaders//tessellation_control.tesc"));
		sb.set_shader(GL_TESS_EVALUATION_SHADER, CGRA_SRCDIR + std::string("//res//shaders//tessellation_eval.tese"));
		sb.set_shader(GL_GEOMETRY_SHADER, CGRA_SRCDIR + std::string("//res//shaders//geometry.geom"));
	}
	sb.set_shader(GL_FRAGMENT_SHADER, CGRA_SRCDIR + std::string("//res//shaders//color_frag.glsl"));
	GLuint shader = sb.build();

	// clear existing patch
	grass_patch.clear();

	// set up grass blade/s
	vec3 arr[4] = { vec3(0, 0, 0), vec3(0, 1, 0), vec3(0, 1, 0), vec3(0, 1, 0) }; // grass blade in center
	grass_model centerBlade;
	centerBlade.setControlPts(arr);
	centerBlade.setMeshes(shader);
	grass_patch.push_back(centerBlade);

	// set up the rest of the patch
	for (int i = 1; i < blade_count; i++) {
		// create random location for blade
		vec3 cp[4] = { vec3(rand() % 20 - 10, 0, rand() % 20 - 10), vec3(0, 1, 0), vec3(0, 1, 0), vec3(0, 1, 0) };
		grass_model grass;
		grass.setControlPts(cp);
		grass.setMeshes(shader);
		grass_patch.push_back(grass);
	}

	m_model.mesh = load_wavefront_data(CGRA_SRCDIR + std::string("/res//assets//teapot.obj")).build();
	m_model.shader = shader;
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

	for (int i = 0; i < grass_patch.size(); i++) {
		grass_model grass = grass_patch.at(i);
		if (!show_rendered_grass) {
			// draw grass blade as bezier and spline
			grass.drawSpline(view, proj);
			grass.drawCurve(view, proj);
		}
		else {
			// draw grass blade rendered
			vec3 camera_position = vec3((m_distance * cos(m_pitch) * cos(m_yaw)), (m_distance * cos(m_pitch) * sin(m_yaw)), (m_distance * sin(m_pitch)));
			grass.drawBlade(view, proj, camera_position);
			//m_model.draw(view, proj);
		}
	}
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

	// add or decrease number of grass blades
	if (ImGui::SliderInt("Grass blades", &blade_count, 1, 100)) {
		// update grass patch vector
		setGrass();
	}

	// show grass or lines
	if (ImGui::Checkbox("Grass render", &show_rendered_grass)) {
		// update grass patch vector
		setGrass();
	}

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
