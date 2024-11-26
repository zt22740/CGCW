#ifndef CALCULATIONS_H
#define CALCULATIONS_H

#include <vector>
#include <map>
#include <glm/glm.hpp>
#include "ModelTriangle.h" // Include your custom ModelTriangle definition here

// Class to handle drawing-related calculations
class Calculations {
public:
    // Compute triangle normals for a list of triangles
    static void computeTriangleNormals(std::vector<ModelTriangle>& triangles);

    // Compute vertex normals for a list of triangles
    static void computeVertexNormals(std::vector<ModelTriangle>& triangles, std::vector<std::vector<glm::vec3>>& vertexNormals);
    static glm::vec3 calculateBarycentricCoords(const glm::vec3 &p, const glm::vec3 &a, const glm::vec3 &b, const glm::vec3 &c);
};

#endif // CALCULATIONS_H
