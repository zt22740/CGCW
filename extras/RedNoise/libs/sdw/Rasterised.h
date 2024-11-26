#ifndef RASTERIZED_H
#define RASTERIZED_H

#include <vector>
#include <glm/glm.hpp>
#include "DrawingWindow.h"
#include "ModelTriangle.h"
#include "CanvasPoint.h"
#include "CanvasTriangle.h"
#include "Colour.h"
#include "TextureMap.h"
#include "Draw.h"

#define WIDTH 320
#define HEIGHT 240
class Rasterised {
public:
    static CanvasPoint projectVertexOntoCanvasPoint(const glm::vec3& cameraPosition, float focalLength, const glm::vec3& vertexPosition, const glm::mat3& cameraOrientation);

    static void renderPointcloud(DrawingWindow &window, const std::vector<ModelTriangle>& model,const glm::vec3& cameraPosition, float focalLength, float scaleFactor, const glm::mat3& cameraOrientation);

    static void renderWireframe(DrawingWindow &window, const std::vector<ModelTriangle>& model,const glm::vec3& cameraPosition, float focalLength, float scaleFactor, std::vector<std::vector<float>>& depthBuffer, const glm::mat3& cameraOrientation);

    static void drawOrbit(DrawingWindow &window, std::vector<ModelTriangle> &model, glm::vec3 &cameraPosition, glm::mat3 &cameraOrientation, float &orbitAngle, float focalLength, float imageScaleFactor, std::vector<std::vector<float>> &depthBuffer);

    static void applyCameraRotation(glm::vec3 &cameraPosition, glm::mat3 &cameraOrientation, float pitch, float yaw);
};

#endif // RASTERIZED_H
