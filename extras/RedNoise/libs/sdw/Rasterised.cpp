#include "Rasterised.h"

void initializeDepthBuffer(std::vector<std::vector<float>>& depthBuffer) {
    for (auto& row : depthBuffer) {
        std::fill(row.begin(), row.end(), 0.0f); // Start with 0, signifying no object at any pixel
    }
}

void Rasterised::applyCameraRotation(glm::vec3 &cameraPosition, glm::mat3 &cameraOrientation, float pitch, float yaw) {
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
CanvasPoint Rasterised::projectVertexOntoCanvasPoint(const glm::vec3& cameraPosition, float focalLength, const glm::vec3& vertexPosition, const glm::mat3& cameraOrientation) {
    // Transform to camera coordinates
    glm::vec3 cameraCoords = (vertexPosition - cameraPosition) * cameraOrientation;

    // Check if the vertex is behind the camera
    if (cameraCoords.z > 0) {
        return CanvasPoint(std::numeric_limits<float>::quiet_NaN(), std::numeric_limits<float>::quiet_NaN());
    }

    // Project the point
    float u = -(focalLength * cameraCoords.x) / cameraCoords.z + (WIDTH / 2.0f);
    float v = (focalLength * cameraCoords.y) / cameraCoords.z + (HEIGHT / 2.0f);
    // Create and return the CanvasPoint
    return CanvasPoint(u, v, -1.0f/cameraCoords.z);
}

// Update the renderPointcloud function
void Rasterised::renderPointcloud(DrawingWindow &window, const std::vector<ModelTriangle>& model, const glm::vec3& cameraPosition, float focalLength, float scaleFactor, const glm::mat3& cameraOrientation) {
    int pointsDrawn = 0;
    int totalPoints = 0;

    for (const auto& triangle : model) {
        for (int i = 0; i < 3; i++) {
            totalPoints++;
            CanvasPoint projectedPoint = Rasterised::projectVertexOntoCanvasPoint(cameraPosition, focalLength, triangle.vertices[i], cameraOrientation);

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

void Rasterised::renderWireframe(DrawingWindow &window, const std::vector<ModelTriangle>& model,
                     const glm::vec3& cameraPosition, float focalLength, float scaleFactor,
                     std::vector<std::vector<float>>& depthBuffer, const glm::mat3& cameraOrientation) {
    window.clearPixels();
    initializeDepthBuffer(depthBuffer);
    
    for (const auto& triangle : model) {
        CanvasTriangle canvasTriangle;
        for (int i = 0; i < 3; i++) {
            CanvasPoint projectedPoint = Rasterised::projectVertexOntoCanvasPoint(cameraPosition, focalLength, triangle.vertices[i], cameraOrientation);
            
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

void Rasterised::drawOrbit(DrawingWindow &window, std::vector<ModelTriangle> &model, glm::vec3 &cameraPosition, glm::mat3 &cameraOrientation, float &orbitAngle, float focalLength, float imageScaleFactor, std::vector<std::vector<float>> &depthBuffer) {
    window.clearPixels();

    // Camera orbit logic (updates camera position and orientation)
    orbitAngle = 0.1f;  // Adjust the orbit speed by changing this value

    // Apply the rotation to camera orientation (only yaw needed for orbit)
    float pitch = 0.0f;  // No pitch change
    applyCameraRotation(cameraPosition, cameraOrientation, pitch, orbitAngle);
    Rasterised::renderWireframe(window, model, cameraPosition, focalLength, imageScaleFactor, depthBuffer, cameraOrientation);
}