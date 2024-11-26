#include "Calculations.h"
void Calculations::computeTriangleNormals(std::vector<ModelTriangle>& triangles) {
    for (ModelTriangle& triangle : triangles) {
        glm::vec3 e0 = triangle.vertices[1] - triangle.vertices[0];
        glm::vec3 e1 = triangle.vertices[2] - triangle.vertices[0];
        triangle.normal = glm::normalize(glm::cross(e0, e1));
    }
}
struct Vec3Comparator {
    bool operator()(const glm::vec3& a, const glm::vec3& b) const {
        if (a.x < b.x) return true;
        if (a.x > b.x) return false;
        if (a.y < b.y) return true;
        if (a.y > b.y) return false;
        return a.z < b.z;
    }
};
void Calculations::computeVertexNormals(std::vector<ModelTriangle>& triangles, std::vector<std::vector<glm::vec3>>& vertexNormals) {
    // Clear and resize the vertexNormals vector
    vertexNormals.clear();
    vertexNormals.resize(triangles.size());

    // First, map to store vertex positions to their connected triangle normals
    std::map<glm::vec3, std::vector<glm::vec3>, Vec3Comparator> vertexToNormalsMap;

    // Compute face normals and group them by vertex
    for (size_t i = 0; i < triangles.size(); ++i) {
        ModelTriangle& triangle = triangles[i];

        // Compute the face normal
        glm::vec3 e0 = triangle.vertices[1] - triangle.vertices[0];
        glm::vec3 e1 = triangle.vertices[2] - triangle.vertices[0];
        glm::vec3 faceNormal = glm::normalize(glm::cross(e0, e1));

        // Add the face normal to each vertex's list of normals
        for (int j = 0; j < 3; ++j) {
            vertexToNormalsMap[triangle.vertices[j]].push_back(faceNormal);
        }
    }

    // Compute vertex normals by averaging the connected face normals
    vertexNormals.resize(triangles.size());
    for (size_t i = 0; i < triangles.size(); ++i) {
        ModelTriangle& triangle = triangles[i];
        vertexNormals[i].clear();

        for (int j = 0; j < 3; ++j) {
            // Get all normals connected to this vertex
            std::vector<glm::vec3>& connectedNormals = vertexToNormalsMap[triangle.vertices[j]];
            
            // Average the normals
            glm::vec3 avgNormal(0.0f);
            for (const auto& normal : connectedNormals) {
                avgNormal += normal;
            }
            avgNormal = glm::normalize(avgNormal / static_cast<float>(connectedNormals.size()));
            
            vertexNormals[i].push_back(avgNormal);
        }
    }
}

glm::vec3 Calculations::calculateBarycentricCoords(const glm::vec3& p, const glm::vec3& a, const glm::vec3& b, const glm::vec3& c) {
    glm::vec3 v0 = b - a;
    glm::vec3 v1 = c - a;
    glm::vec3 v2 = p - a;

    float d00 = glm::dot(v0, v0);
    float d01 = glm::dot(v0, v1);
    float d11 = glm::dot(v1, v1);
    float d20 = glm::dot(v2, v0);
    float d21 = glm::dot(v2, v1);

    float denom = d00 * d11 - d01 * d01;

    // Check if denom is very close to zero
    if (std::abs(denom) < 1e-6f) {
        return glm::vec3(-1.0f, -1.0f, -1.0f); // return invalid coordinates if triangle is degenerate
    }

    float v = (d11 * d20 - d01 * d21) / denom;
    float w = (d00 * d21 - d01 * d20) / denom;
    float u = 1.0f - v - w;

    return glm::vec3(u, v, w);
}
