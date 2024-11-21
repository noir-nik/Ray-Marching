#include <vector>
#include <chrono>
#include <omp.h>
#include <glad/gl.h>

#include "fragment.h"
#include "settings.h"
#include "window.h"
#include "save_bmp.h"

#define RENDER_DIR "output/"

struct ShaderReader
{
	ShaderReader(const std::string &filePath)
	{
		std::ifstream fileStream(filePath, std::ios::in | std::ios::binary);
		if (fileStream.is_open())
		{
			fileStream.seekg(0, std::ios::end);
			size_t size = fileStream.tellg();
			fileStream.seekg(0, std::ios::beg);

			source.resize(size + 1);
			fileStream.read(source.data(), size);
			source[size] = '\0';

			data = source.data();
			this->size = size;

			fileStream.close();
		}
		else
		{
			std::cerr << "Unable to open file: " << filePath << std::endl;
			exit(1);
		}
	}

	std::vector<GLchar> source;
	const GLchar *data;
	GLint size;
};

// Vertex shader
const char *vertexShaderSource = R"(
    #version 330
    layout (location = 0) in vec2 aPos;
    void main() {
        gl_Position = vec4(aPos.x, aPos.y, 0.0, 1.0);
    }
)";

std::vector<uint32_t> framebuffer;

void clear_framebuffer(uint32_t color) {
  for (size_t i = 0; i < WINDOW_WIDTH * WINDOW_HEIGHT; i++) {
    framebuffer[i] = color;
  }
}

void draw_pixel(uint16_t x, uint16_t y, glm::vec4 color_vec) {
    uint32_t color;
    // ARGB8888
    color = glm::clamp(int(color_vec.b * 255.0f), 0, 255);
    color += glm::clamp(int(color_vec.g * 255.0f), 0, 255) << 8;
    color += glm::clamp(int(color_vec.r * 255.0f), 0, 255) << 16;
    color += glm::clamp(int(color_vec.a * 255.0f), 0, 255) << 24;
    framebuffer[(WINDOW_WIDTH * y) + x] = color;
}

int main(int argc, char *argv[])
{
	int success;
	char infoLog[512];
	std::chrono::high_resolution_clock::time_point start, start_copy, end;
	std::chrono::nanoseconds elapsed;

	framebuffer.resize(WINDOW_WIDTH* WINDOW_HEIGHT);

	// Clear CPU framebuffer
	clear_framebuffer(0x00000000);

#ifdef RENDER_SINGLE_THREAD
	// Single thread render
	printf("CPU render with 1 thread starts...\n");
	start = std::chrono::high_resolution_clock::now();
	for (size_t i = 0; i < WINDOW_HEIGHT; i++) {
		for (size_t j = 0; j < WINDOW_WIDTH; j++){
			glm::vec4 color_vec = fragment({j, i}, {WINDOW_WIDTH, WINDOW_HEIGHT});
			draw_pixel(j, i, color_vec);
		}
	}
	end = std::chrono::high_resolution_clock::now();
	elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start);
	printf("CPU render - 1 thread: %.3lf ms\n", (double)elapsed.count() / 1000000.0);
#endif

	// Multi thread render
	printf("CPU render with %d threads starts...\n", omp_get_max_threads());
	start = std::chrono::high_resolution_clock::now();

#pragma omp parallel for
	for (size_t i = 0; i < WINDOW_HEIGHT; i++) {
		for (size_t j = 0; j < WINDOW_WIDTH; j++) {
			glm::vec4 color_vec = fragment({j, i}, {WINDOW_WIDTH, WINDOW_HEIGHT});
			draw_pixel(j, i, color_vec);
		}
	}

	end = std::chrono::high_resolution_clock::now();
	elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start);
	printf("CPU render - %d threads: %.3lf ms\n", omp_get_max_threads(), (double)elapsed.count() / 1000000.0);
	generateBitmapImage(framebuffer.data(), WINDOW_HEIGHT, WINDOW_WIDTH, RENDER_DIR"out_cpu.bmp", true);
	printf("CPU render result saved to %s\n", RENDER_DIR"out_cpu.bmp");

	// GPU render
	if (!create_window()) {
		std::cerr << "Failed to create window" << std::endl;
		return -1;
	}

	// Vertex shader
	GLuint vertexShader = glCreateShader(GL_VERTEX_SHADER);
	glShaderSource(vertexShader, 1, &vertexShaderSource, nullptr);
	glCompileShader(vertexShader);

	glGetShaderiv(vertexShader, GL_COMPILE_STATUS, &success);
	if (!success)
	{
		glGetShaderInfoLog(vertexShader, 512, NULL, infoLog);
		std::cout << "ERROR::SHADER::VERTEX::COMPILATION_FAILED\n"
				  << infoLog << std::endl;
		return -1;
	}

	// Fragment shader
	GLuint fragmentShader = glCreateShader(GL_FRAGMENT_SHADER);
	ShaderReader sr(FRAGMENT_SOURCE);
	glShaderSource(fragmentShader, 1, &sr.data, &sr.size);
	glCompileShader(fragmentShader);
	glGetShaderiv(fragmentShader, GL_COMPILE_STATUS, &success);
	if (!success)
	{
		glGetShaderInfoLog(fragmentShader, 512, NULL, infoLog);
		std::cout << "ERROR::SHADER::FRAGMENT::COMPILATION_FAILED\n"
				  << infoLog << std::endl;
		return -1;
	}

	// Shader program
	GLuint shaderProgram = glCreateProgram();

	glAttachShader(shaderProgram, vertexShader);

	glAttachShader(shaderProgram, fragmentShader);

	glLinkProgram(shaderProgram);

	glDetachShader(shaderProgram, vertexShader);
	glDetachShader(shaderProgram, fragmentShader);
	glDeleteShader(vertexShader);
	glDeleteShader(fragmentShader);
	glGetProgramiv(shaderProgram, GL_LINK_STATUS, &success);
	if (!success)
	{
		glGetProgramInfoLog(shaderProgram, 512, NULL, infoLog);
		std::cout << "ERROR::SHADER::PROGRAM::COMPILATION_FAILED\n"
				  << infoLog << std::endl;
		return -1;
	}

	printf("GPU render starts...\n");
	glFinish();
	start_copy = std::chrono::high_resolution_clock::now();

	// Uniforms
	GLint resolutionLocation = glGetUniformLocation(shaderProgram, "iResolution");
	GLint timeLocation = glGetUniformLocation(shaderProgram, "iTime");
	GLint maxStepsLocation = glGetUniformLocation(shaderProgram, "uMaxSteps");
	GLint maxDistLocation = glGetUniformLocation(shaderProgram, "uMaxDist");
	GLint epsLocation = glGetUniformLocation(shaderProgram, "uEps");
	GLint AALocation = glGetUniformLocation(shaderProgram, "uAa");
	GLint reflectLocation = glGetUniformLocation(shaderProgram, "uReflect");

	// Camera
	GLint posLocation = glGetUniformLocation(shaderProgram, "uCamPos");
	GLint lookLocation = glGetUniformLocation(shaderProgram, "uLookAt");
	GLint FOVLocation = glGetUniformLocation(shaderProgram, "uFov");

	// Vertex data for two triangles forming the screen
	float vertices[] = {
		-1.0f, 1.0f,  // top left
		-1.0f, -1.0f, // bottom left
		1.0f, -1.0f,  // bottom right
		1.0f, 1.0f	  // top right
	};

	unsigned int indices[] = {
		0, 1, 2, // first triangle
		0, 2, 3	 // second triangle
	};

	// VAO, VBO, EBO
	GLuint VAO, VBO, EBO;
	glGenVertexArrays(1, &VAO);
	glGenBuffers(1, &VBO);
	glGenBuffers(1, &EBO);

	glBindVertexArray(VAO);

	glBindBuffer(GL_ARRAY_BUFFER, VBO);
	glBufferData(GL_ARRAY_BUFFER, sizeof(vertices), vertices, GL_STATIC_DRAW);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, EBO);
	glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(indices), indices, GL_STATIC_DRAW);

	glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, 2 * sizeof(float), (void *)0);
	glEnableVertexAttribArray(0);

	glBindVertexArray(0);

	// Copy uniforms
	glUseProgram(shaderProgram);
	glBindVertexArray(VAO);

	glUniform2f(resolutionLocation, WINDOW_WIDTH, WINDOW_HEIGHT);
	glUniform1i(AALocation, AA);
	glUniform1i(maxStepsLocation, MAX_STEPS);
	glUniform1f(maxDistLocation, MAX_DIST);
	glUniform1f(epsLocation, EPS);
	glUniform1i(reflectLocation, REFLECT);

	glUniform3f(posLocation, POS[0], POS[1], POS[2]);
	glUniform3f(lookLocation, LOOKAT[0], LOOKAT[1], LOOKAT[2]);
	glUniform1f(FOVLocation, FOV);

	glFinish();
	start = std::chrono::high_resolution_clock::now();

	// Render
	glDrawElements(GL_TRIANGLES, 6, GL_UNSIGNED_INT, 0);

	glFinish();
	end = std::chrono::high_resolution_clock::now();

	auto elapsed_render = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start);
	auto elapsed_copy = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start_copy);
	auto copy_time = std::chrono::duration_cast<std::chrono::nanoseconds>(start - start_copy);
	printf("GPU Render time: %.3lf ms\n", (double)elapsed_render.count() / 1000000.0);
	printf("Copy + GPU render time: %.3lf ms\n", (double)elapsed_copy.count() / 1000000.0);
	printf("Copy time: %.3lf ms\n", (double)copy_time.count() / 1000000.0);

	// Read pixels
	glReadBuffer(GL_BACK);
	glPixelStorei(GL_PACK_ALIGNMENT, 1);
	glReadPixels(0, 0, WINDOW_WIDTH, WINDOW_HEIGHT, GL_BGRA, GL_UNSIGNED_BYTE, framebuffer.data());

	// Save image
	generateBitmapImage(framebuffer.data(), WINDOW_HEIGHT, WINDOW_WIDTH, RENDER_DIR"out_gpu.bmp", false);
	printf("GPU render result saved to: %s\n", RENDER_DIR"out_gpu.bmp");

	// Cleanup
	glDeleteVertexArrays(1, &VAO);
	glDeleteBuffers(1, &VBO);
	glDeleteBuffers(1, &EBO);
	glDeleteProgram(shaderProgram);

	destroy_window();
	return 0;
}
