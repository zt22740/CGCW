#ifndef RAYTRACED_H
#define RAYTRACED_H

#include <vector>
#include <glm/glm.hpp>
#include "ModelTriangle.h"
#include "DrawingWindow.h"
#include "TextureMap.h"
#include "RayTriangleIntersection.h"
#include "Calculations.h"
#include <bits/algorithmfwd.h>
#include "Draw.h"
#include <map>
#include <algorithm>

#define WIDTH 320
#define HEIGHT 240
class Raytraced {
public:
static void drawRayTracedScene(const glm::vec3& cameraPosition, float focalLength, float scaleFactor, std::vector<ModelTriangle>& triangles, DrawingWindow &window, TextureMap& textureMap, int maxReflectionDepth);
static void drawRayTracedSceneGouraud(const glm::vec3& cameraPosition, float focalLength, float scaleFactor, std::vector<ModelTriangle>& triangles, std::vector<std::vector<glm::vec3>>& vertexNormals, DrawingWindow &window, TextureMap &textureMap);
};

#endif // RAYTRACED_H