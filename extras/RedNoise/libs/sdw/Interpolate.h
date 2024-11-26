#ifndef INTERPOLATE_H
#define INTERPOLATE_H

#include <vector>
#include <glm/glm.hpp> // For glm::vec3
#include "TexturePoint.h" // Include your TexturePoint definition
#include "CanvasPoint.h" // Include your CanvasPoint definition
class Interpolate {
public:
// Function to interpolate single float values
static std::vector<float> interpolateSingleFloats(float from, float to, int numberOfValues);

// Function to interpolate three-element vectors (glm::vec3)
static std::vector<glm::vec3> interpolateThreeElementValues(glm::vec3 from, glm::vec3 to, int numberOfValues);

// Function to interpolate TexturePoints
static std::vector<TexturePoint> interpolateTexturePoints(TexturePoint from, TexturePoint to, float numberOfValues);

// Function to interpolate CanvasPoints
static std::vector<CanvasPoint> interpolateCanvasPoints(CanvasPoint from, CanvasPoint to, float numberOfValues);
};

#endif // INTERPOLATE_H
