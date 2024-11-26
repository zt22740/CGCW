#ifndef DRAW_H
#define DRAW_H

#include <vector>
#include <cstdint>
#include "DrawingWindow.h" // Include your drawing window class
#include "CanvasTriangle.h" // Include the CanvasTriangle definition
#include "CanvasPoint.h"    // Include the CanvasPoint definition
#include "TextureMap.h"     // Include your TextureMap definition
#include "Colour.h"
#include "Interpolate.h"
#include "RayTriangleIntersection.h"
#include <map>
#include "Calculations.h"

#define WIDTH 320
#define HEIGHT 240
class Draw {
public:
    static void drawLine(CanvasPoint from, CanvasPoint to, Colour colour, DrawingWindow &window);
    static void drawTriangle(const CanvasTriangle triangle, const Colour colour, DrawingWindow &window);
    static void drawFilledTriangle(CanvasTriangle triangle, Colour colour, DrawingWindow &window, std::vector<std::vector<float>>& depthBuffer);
    static void fillTexturedTriangle(CanvasTriangle triangle, TextureMap &textureMap, DrawingWindow &window);
    static void draw(DrawingWindow &window);
    static void drawPoint(DrawingWindow &window, const CanvasPoint &point, Colour colour);
    static void drawRayTracedScene(const glm::vec3& cameraPosition, float focalLength, float scaleFactor,  std::vector<ModelTriangle>& triangles, DrawingWindow &window, TextureMap& texture);
    static void computeTriangleNormals(std::vector<ModelTriangle>& triangles);
    static void computeVertexNormals(std::vector<ModelTriangle>& triangles, std::vector<std::vector<glm::vec3>>& vertexNormals);
    static void drawRayTracedSceneGouraud(const glm::vec3& cameraPosition, float focalLength, float scaleFactor, std::vector<ModelTriangle>& triangles, std::vector<std::vector<glm::vec3>>& vertexNormals, DrawingWindow &window, TextureMap &textureMap);
    static glm::vec3 traceRay(const glm::vec3& origin, const glm::vec3& direction, const std::vector<ModelTriangle>& triangles, const glm::vec3& lightPosition, int depth, int maxDepth, const std::map<int, float>& reflectiveTriangles);
};

#endif // DRAW_H
