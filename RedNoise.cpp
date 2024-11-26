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
//#include <Render.h>

#define WIDTH 320
#define HEIGHT 240

// void initializeDepthBuffer(std::vector<std::vector<float>>& depthBuffer) {
//     for (auto& row : depthBuffer) {
//         std::fill(row.begin(), row.end(), std::numeric_limits<float>::lowest());
//     }
// }

void initializeDepthBuffer(std::vector<std::vector<float>>& depthBuffer) {
    for (auto& row : depthBuffer) {
        std::fill(row.begin(), row.end(), 0.0f); // Start with 0, signifying no object at any pixel
    }
}

CanvasTriangle generateRandomTriangle(int width, int height) {
    // Generate three random points within the given width and height
    CanvasPoint v0(rand() % width, rand() % height);
    CanvasPoint v1(rand() % width, rand() % height);
    CanvasPoint v2(rand() % width, rand() % height);

    // Return a CanvasTriangle constructed from these random points
    return CanvasTriangle(v0, v1, v2);
}

// void drawFilledTriangle(CanvasTriangle triangle, Colour colour, DrawingWindow &window, std::vector<std::vector<float>>& depthBuffer) {
//     CanvasTriangle sorted = sortVertices(triangle);
    
//     int yStart = (int)sorted[0].y;
//     int yEnd = (int)sorted[2].y;
//     int yMid = (int)sorted[1].y;

//     if (yStart < 0) yStart = 0;
//     if (yEnd >= HEIGHT) yEnd = HEIGHT - 1;

//     for (int y = yStart; y <= yEnd; y++) {
//         float t1, t2;
//         int x1, x2;
//         float z1, z2;

//         if (y < sorted[1].y) {
//             float yDiff1 = sorted[1].y - sorted[0].y;
//             float yDiff2 = sorted[2].y - sorted[0].y;
//             if (yDiff1 == 0){
//                 t1 = 0;
//             }
//             else{
//                 t1 = (y - sorted[0].y)/yDiff1;
//             }
//             if (yDiff2 == 0){
//                 t2 = 0;
//             }
//             else{
//                 t2 = (y - sorted[0].y)/yDiff2;
//             }

//             x1 = (int)(sorted[0].x + t1 * (sorted[1].x - sorted[0].x));
//             x2 = (int)(sorted[0].x + t2 * (sorted[2].x - sorted[0].x));

//             z1 = 1.0f / (1.0f / sorted[0].depth + t1 * (1.0f / sorted[1].depth - 1.0f / sorted[0].depth));
//             z2 = 1.0f / (1.0f / sorted[0].depth + t2 * (1.0f / sorted[2].depth - 1.0f / sorted[0].depth));
//         } else {
//             float yDiff1 = sorted[2].y - sorted[1].y;
//             float yDiff2 = sorted[2].y - sorted[0].y;

//             if (yDiff1 == 0){
//                 t1 = 0;
//             }
//             else{
//                 t1 = (y - sorted[1].y)/yDiff1;
//             }
//             if (yDiff2 == 0){
//                 t2 = 0;
//             }
//             else{
//                 t2 = (y - sorted[0].y)/yDiff2;
//             }

//             x1 = (int)(sorted[1].x + t1 * (sorted[2].x - sorted[1].x));
//             x2 = (int)(sorted[0].x + t2 * (sorted[2].x - sorted[0].x));

//             z1 = 1.0f / (1.0f / sorted[1].depth + t1 * (1.0f / sorted[2].depth - 1.0f / sorted[1].depth));
//             z2 = 1.0f / (1.0f / sorted[0].depth + t2 * (1.0f / sorted[2].depth - 1.0f / sorted[0].depth));
//         }

//         if (x1 > x2) {
//             std::swap(x1, x2);
//             std::swap(z1, z2);
//         }

//         int xStart = x1;
//         int xEnd = x2;

//         if (xStart < 0) xStart = 0;
//         if (xEnd >= WIDTH) xEnd = WIDTH - 1;

//         for (int x = xStart; x <= xEnd; x++) {
//             float t;
//             float xDiff = x2 - x1;
//             if (xDiff == 0){
//                 t = 0;
//             }
//             else{
//                 t = (x - x1) / xDiff;
//             }
//             //float t = (x - x1) / (xDiff == 0 ? 1 : xDiff);
//             float z = 1.0f / (1.0f / z1 + t * (1.0f / z2 - 1.0f / z1));

//             if (z > depthBuffer[y][x]) {
//                 depthBuffer[y][x] = z;
//                 uint32_t fillcolour = (255 << 24) + (uint32_t(colour.red) << 16) + (uint32_t(colour.green) << 8) + uint32_t(colour.blue);
//                 window.setPixelColour(x, y, fillcolour);
//             }
//             //std::cout << "y: " << y << ", x1: " << x1 << ", x2: " << x2 << ", z1: " << z1 << ", z2: " << z2 << std::endl;

//         }
//     }
// }

void applyCameraRotation(glm::vec3 &cameraPosition, glm::mat3 &cameraOrientation, float pitch, float yaw) {
    // Rotation matrix for pitch (X-axis rotation)
    glm::mat3 rotationMatrixX = glm::mat3(
        glm::vec3(1, 0, 0),
        glm::vec3(0, cos(glm::radians(pitch)), -sin(glm::radians(pitch))),
        glm::vec3(0, sin(glm::radians(pitch)), cos(glm::radians(pitch)))
    );
    
    // Rotation matrix for yaw (Y-axis rotation)
    glm::mat3 rotationMatrixY = glm::mat3(
        glm::vec3(cos(glm::radians(yaw)), 0, sin(glm::radians(yaw))),
        glm::vec3(0, 1, 0),
        glm::vec3(-sin(glm::radians(yaw)), 0, cos(glm::radians(yaw)))
    );
    
    // Apply the rotations to the camera orientation
    cameraOrientation = rotationMatrixY * cameraOrientation;
    cameraOrientation = rotationMatrixX * cameraOrientation;

    // Apply the rotations to the camera position
    cameraPosition = rotationMatrixY * cameraPosition;
    cameraPosition = rotationMatrixX * cameraPosition;
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
        applyCameraRotation(cameraPosition, cameraOrientation, pitch, yaw);
	} else if (event.type == SDL_MOUSEBUTTONDOWN) {
		window.savePPM("output.ppm");
		window.saveBMP("output.bmp");
	}
}

CanvasPoint projectVertexOntoCanvasPoint(const glm::vec3& cameraPosition, float focalLength, const glm::vec3& vertexPosition, const glm::mat3& cameraOrientation) {
    // Transform to camera coordinates
    glm::vec3 cameraCoords = (vertexPosition - cameraPosition) * cameraOrientation;

    // Check if the vertex is behind the camera
    if (cameraCoords.z > 0) {
        return CanvasPoint(std::numeric_limits<float>::quiet_NaN(), std::numeric_limits<float>::quiet_NaN());
    }

    // Project the point
    float u = -(focalLength * cameraCoords.x) / cameraCoords.z + (WIDTH / 2.0f);
    float v = (focalLength * cameraCoords.y) / cameraCoords.z + (HEIGHT / 2.0f);
    //std::cout << "vertexX" << u << "y" << v << "\n";
    // Create and return the CanvasPoint
    return CanvasPoint(u, v, -1.0f/cameraCoords.z);
}

// Update the renderPointcloud function
void renderPointcloud(DrawingWindow &window, const std::vector<ModelTriangle>& model, const glm::vec3& cameraPosition, float focalLength, float scaleFactor, const glm::mat3& cameraOrientation) {
    int pointsDrawn = 0;
    int totalPoints = 0;

    for (const auto& triangle : model) {
        for (int i = 0; i < 3; i++) {
            totalPoints++;
            CanvasPoint projectedPoint = projectVertexOntoCanvasPoint(cameraPosition, focalLength, triangle.vertices[i], cameraOrientation);

            if (projectedPoint.x != -1 && projectedPoint.y != -1) {
                // Apply scaling factor
                projectedPoint.x = (projectedPoint.x - window.width / 2) * scaleFactor + window.width / 2;
                projectedPoint.y = (projectedPoint.y - window.height / 2) * scaleFactor + window.height / 2;
                
                // Check if the point is within the window bounds
                if (projectedPoint.x >= 0 && projectedPoint.x < window.width && 
                    projectedPoint.y >= 0 && projectedPoint.y < window.height) {
                    Draw::drawPoint(window, projectedPoint, Colour(255, 255, 255)); // White color for all points
                    pointsDrawn++;
                }
            }
        }
    }
}

void renderWireframe(DrawingWindow &window, const std::vector<ModelTriangle>& model,
                     const glm::vec3& cameraPosition, float focalLength, float scaleFactor,
                     std::vector<std::vector<float>>& depthBuffer, const glm::mat3& cameraOrientation) {
    window.clearPixels();
    initializeDepthBuffer(depthBuffer);
    
    for (const auto& triangle : model) {
        CanvasTriangle canvasTriangle;
        for (int i = 0; i < 3; i++) {
            CanvasPoint projectedPoint = projectVertexOntoCanvasPoint(cameraPosition, focalLength, triangle.vertices[i], cameraOrientation);
            
            // Apply scaling factor
            projectedPoint.x = (projectedPoint.x - window.width / 2) * scaleFactor + window.width / 2;
            projectedPoint.y = (projectedPoint.y - window.height / 2) * scaleFactor + window.height / 2;
            
            // Calculate depth (inverse of z)
           // projectedPoint.depth = 1.0f / glm::length(cameraPosition - triangle.vertices[i]);
            
            canvasTriangle.vertices[i] = projectedPoint;
        }
            //drawTriangle(canvasTriangle, triangle.colour, window);
            Draw::drawFilledTriangle(canvasTriangle, triangle.colour, window, depthBuffer);
    }
}

void drawOrbit(DrawingWindow &window, std::vector<ModelTriangle> &model, glm::vec3 &cameraPosition, glm::mat3 &cameraOrientation, float &orbitAngle, float focalLength, float imageScaleFactor, std::vector<std::vector<float>> &depthBuffer) {
    window.clearPixels();

    // Camera orbit logic (updates camera position and orientation)
    orbitAngle = 0.1f;  // Adjust the orbit speed by changing this value

    // Apply the rotation to camera orientation (only yaw needed for orbit)
    float pitch = 0.0f;  // No pitch change
    applyCameraRotation(cameraPosition, cameraOrientation, pitch, orbitAngle);
    renderWireframe(window, model, cameraPosition, focalLength, imageScaleFactor, depthBuffer, cameraOrientation);
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

    Draw::computeTriangleNormals(model);
    std::vector<std::vector<glm::vec3>> vertexNormals;
    Draw::computeVertexNormals(model, vertexNormals);
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
        // renderPointcloud(window, model, cameraPosition, focalLength, imageScaleFactor, cameraOrientation);
        // renderWireframe(window, model, cameraPosition, focalLength, imageScaleFactor, depthBuffer, cameraOrientation);
        //drawOrbit(window, model, cameraPosition, cameraOrientation, orbitAngle, focalLength, imageScaleFactor, depthBuffer);
        Draw::drawRayTracedScene(cameraPosition, focalLength, imageScaleFactor, model, window, textureMap);
        // Draw::drawRayTracedSceneGouraud(cameraPosition, focalLength, imageScaleFactor, model, vertexNormals, window, textureMap);

		// Need to render the frame at the end, or nothing actually gets shown on the screen !
		window.renderFrame();
	}
}
