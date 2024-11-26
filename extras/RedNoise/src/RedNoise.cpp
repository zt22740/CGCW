#include <CanvasTriangle.h>
#include <DrawingWindow.h>
#include <Utils.h>
#include <fstream>
#include <vector>
#include <iostream>
#include <glm/glm.hpp>
#include <CanvasPoint.h>
#include <Colour.h>
#include <TextureMap.h>
#include <ModelTriangle.h>
#include <Utils.h>
#include <unordered_map>
#include <sstream>
#include <Interpolate.h>
#include <Draw.h>
#include <LoadModel.h>
#include <RayTriangleIntersection.h>
#include <map>
#include <Calculations.h>
#include <Rasterised.h>

#define WIDTH 320
#define HEIGHT 240

// void initializeDepthBuffer(std::vector<std::vector<float>>& depthBuffer) {
//     for (auto& row : depthBuffer) {
//         std::fill(row.begin(), row.end(), std::numeric_limits<float>::lowest());
//     }
// }

CanvasTriangle generateRandomTriangle(int width, int height) {
    // Generate three random points within the given width and height
    CanvasPoint v0(rand() % width, rand() % height);
    CanvasPoint v1(rand() % width, rand() % height);
    CanvasPoint v2(rand() % width, rand() % height);

    // Return a CanvasTriangle constructed from these random points
    return CanvasTriangle(v0, v1, v2);
}

void handleEvent(SDL_Event event, DrawingWindow &window, glm::vec3 &cameraPosition, glm::mat3 &cameraOrientation) {
    float pitch = 0.0f;  // Rotation around X-axis
    float yaw = 0.0f;    // Rotation around Y-axis

	if (event.type == SDL_KEYDOWN) {
		if (event.key.keysym.sym == SDLK_LEFT) {
            cameraPosition.x += 0.1f;
            }
		else if (event.key.keysym.sym == SDLK_RIGHT){
            cameraPosition.x -= 0.1f;
        }
		else if (event.key.keysym.sym == SDLK_UP){
            cameraPosition.y -= 0.1f;
        }
		else if (event.key.keysym.sym == SDLK_DOWN){
            cameraPosition.y += 0.1f;
        }
        else if (event.key.keysym.sym == SDLK_f){
            cameraPosition.z -= 0.1f;
        }
        else if (event.key.keysym.sym == SDLK_b){
            cameraPosition.z += 0.1f;
        }
        else if (event.key.keysym.sym == SDLK_j){
            yaw += 1.0f;
        }
        else if (event.key.keysym.sym == SDLK_l){
            yaw -= 1.0f;
        }
        else if (event.key.keysym.sym == SDLK_i){
            pitch += 1.0f;
        }
        else if (event.key.keysym.sym == SDLK_k){
            pitch -= 1.0f;
        }
		else if (event.key.keysym.sym == SDLK_u) {
            // Generate random vertices for the triangle
            CanvasTriangle triangle = generateRandomTriangle(window.width, window.height);

            // Generate random RGB color
            Colour color(rand() % 256, rand() % 256, rand() % 256);

            // Call the triangle drawing function
            Draw::drawTriangle(triangle, color, window);
        }
		// else if (event.key.keysym.sym == SDLK_f){
		// 	CanvasTriangle triangle = generateRandomTriangle(window.width, window.height);
		// 	Colour color(rand() % 256, rand() % 256, rand() % 256);
		// 	//drawFilledTriangle(triangle, color, window);
		// 	drawTriangle(triangle, Colour(255,255,255), window);
		// }
        Rasterised::applyCameraRotation(cameraPosition, cameraOrientation, pitch, yaw);
	} else if (event.type == SDL_MOUSEBUTTONDOWN) {
		window.savePPM("output.ppm");
		window.saveBMP("output.bmp");
	}
}


int main(int argc, char *argv[]) {
	DrawingWindow window = DrawingWindow(WIDTH, HEIGHT, false);
	SDL_Event event;
	std::vector<float> result;
	result = Interpolate::interpolateSingleFloats(2.2, 8.5, 7);
	for(size_t i=0; i<result.size(); i++) std::cout << result[i] << " ";
	std::cout << std::endl;
	Colour colour = Colour{255, 255, 255};
	// CanvasPoint topLeft = CanvasPoint(0, 0);
	// CanvasPoint topRight = CanvasPoint(WIDTH, 0);
	// CanvasPoint centre = CanvasPoint(WIDTH/2, HEIGHT/2);
    TextureMap textureMap = TextureMap("texture.ppm");
    // CanvasPoint v0(160, 10); 
    // CanvasPoint v1(300, 230);
    // CanvasPoint v2(10, 150);

    // // Set the corresponding texture points
    // v0.texturePoint = TexturePoint(195.0f / textureMap.width, 5.0f / textureMap.height);
    // v1.texturePoint = TexturePoint(395.0f / textureMap.width, 380.0f / textureMap.height);
    // v2.texturePoint = TexturePoint(65.0f / textureMap.width, 330.0f / textureMap.height);
    // CanvasTriangle triangle(v0, v1, v2);
    //Draw the textured triangle
    // Draw::fillTexturedTriangle(triangle, textureMap, window);
    std::string objFilePath = "models/cornell-box.obj";
    std::string mtlFilePath = "models/cornell-box.mtl";
    std::string TobjFilePath = "models/textured-cornell-box.obj";
    std::string TmtlFilePath = "models/textured-cornell-box.mtl";
    float scaleFactor = 0.35;
    float imageScaleFactor = 160.0f;

    std::vector<ModelTriangle> model = LoadModel::loadOBJ(TobjFilePath, TmtlFilePath, scaleFactor);

    Calculations::computeTriangleNormals(model);
    std::vector<std::vector<glm::vec3>> vertexNormals;
    Calculations::computeVertexNormals(model, vertexNormals);
    glm::vec3 cameraPosition(0.0f, 0.0f, 4.0f); // Camera is positioned at (0, 0, 4)
    glm::mat3 cameraOrientation = glm::mat3(
    glm::vec3(1, 0, 0),  // right (X-axis)
    glm::vec3(0, 1, 0),  // up (Y-axis)
    glm::vec3(0, 0, 1)   // forward (Z-axis)
);
    float focalLength = 2.0f; // Adjust based on your requirements
    float orbitAngle = 0.0f;
    std::vector<std::vector<float>> depthBuffer(HEIGHT, std::vector<float>(WIDTH, 0.0f));
    // std::map<int, float> reflectiveTriangles = {
    //     {0, 1.0f}, // Fully reflective triangle
    //     {1, 0.5f}, // Partially reflective triangle
    //     // Add more mappings as needed
    // };
	while (true) {
		// We MUST poll for events - otherwise the window will freeze !
		if (window.pollForInputEvents(event)) handleEvent(event, window, cameraPosition, cameraOrientation);
		// Draw::draw(window);
		// Draw::drawLine(topLeft, centre, colour, window);
		// Draw::drawLine(topRight, centre, colour, window);
		// Draw::drawLine(CanvasPoint(WIDTH/2, 0), CanvasPoint(WIDTH/2, HEIGHT), colour, window);
		// Draw::drawLine(CanvasPoint(WIDTH/3, HEIGHT/2), CanvasPoint((WIDTH/3)*2, HEIGHT/2), colour, window);
        Rasterised::renderPointcloud(window, model, cameraPosition, focalLength, imageScaleFactor, cameraOrientation);
        // Rasterised::renderWireframe(window, model, cameraPosition, focalLength, imageScaleFactor, depthBuffer, cameraOrientation);
        // Rasterised::drawOrbit(window, model, cameraPosition, cameraOrientation, orbitAngle, focalLength, imageScaleFactor, depthBuffer);
        // Draw::drawRayTracedScene(cameraPosition, focalLength, imageScaleFactor, model, window, textureMap);
        // Draw::drawRayTracedSceneGouraud(cameraPosition, focalLength, imageScaleFactor, model, vertexNormals, window, textureMap);

		// Need to render the frame at the end, or nothing actually gets shown on the screen !
		window.renderFrame();
	}
}