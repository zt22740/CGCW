#include "Draw.h"
#include <unordered_map>
#include <map>
#include <algorithm>
#include <iostream> // Include this if not already included
#include <bits/algorithmfwd.h>
#include <algorithm>

CanvasTriangle sortVertices(const CanvasTriangle& triangle) {
    CanvasTriangle sorted = triangle;
    if (sorted[0].y > sorted[1].y) std::swap(sorted[0], sorted[1]);
    if (sorted[1].y > sorted[2].y) std::swap(sorted[1], sorted[2]);
    if (sorted[0].y > sorted[1].y) std::swap(sorted[0], sorted[1]);
    return sorted;
}


// Draw a line
void Draw::drawLine(CanvasPoint from, CanvasPoint to, Colour colour, DrawingWindow &window) {
    float xDiff = to.x - from.x;
    float yDiff = to.y - from.y;
    float steps = fmax(abs(xDiff), abs(yDiff));
    float xStepSize = xDiff / steps;
    float yStepSize = yDiff / steps;
    for (float i = 0.0; i < steps; i++) {
        float x = from.x + (xStepSize * i);
        float y = from.y + (yStepSize * i);
        uint32_t ncolour = (255 << 24) + (uint32_t(colour.red) << 16) + (uint32_t(colour.green) << 8) + uint32_t(colour.blue);
        window.setPixelColour(x, y, ncolour);
    }
}

// Draw a triangle outline
void Draw::drawTriangle(const CanvasTriangle triangle, Colour colour, DrawingWindow &window) {
    drawLine(triangle[0], triangle[1], colour, window);
    drawLine(triangle[0], triangle[2], colour, window);
    drawLine(triangle[1], triangle[2], colour, window);
}
// Fill a top triangle
void fill_top_triangle(CanvasPoint v0, CanvasPoint v1, CanvasPoint v2, Colour colour, DrawingWindow &window, std::vector<std::vector<float>>& depthBuffer) {
    float stepX1 = (v1.x - v0.x) / (v1.y - v0.y);
    float stepX2 = (v2.x - v0.x) / (v2.y - v0.y);

    for (float y = v0.y; y <= v1.y; y++) {
        float x1 = v0.x + (y - v0.y) * stepX1;
        float x2 = v0.x + (y - v0.y) * stepX2;

        float z1 = (v0.depth + (y - v0.y) * ((v1.depth - v0.depth) / (v1.y - v0.y)));
        float z2 = (v0.depth + (y - v0.y) * ((v2.depth - v0.depth) / (v2.y - v0.y)));
        if (y < 0 || y >= HEIGHT) continue;
        if (x1 > x2) {
            std::swap(x1, x2);
            std::swap(z1, z2);
        }
        for (float x = std::max(0.0f, std::min(x1, x2)); x <= std::min(float(WIDTH - 1), std::max(x1, x2)); x++) {
            if (x < 0 || x > WIDTH) continue;
            // Interpolate Z across the horizontal line
            float t = (x - x1) / (x2 - x1);
            float zInterpolated = (z1 + t * (z2 - z1));
            
            // Check and set pixel if depth test passes
            if (zInterpolated > depthBuffer[y][x]) {
                depthBuffer[y][x] = zInterpolated;
                window.setPixelColour(x, y, (255 << 24) + (uint32_t(colour.red) << 16) + (uint32_t(colour.green) << 8) + uint32_t(colour.blue));
            }
        }
    }
}

// Fill a lower triangle
void fill_lower_triangle(CanvasPoint v0, CanvasPoint v1, CanvasPoint v2, Colour colour, DrawingWindow &window, std::vector<std::vector<float>>& depthBuffer) {
    float stepX1 = (v2.x - v0.x) / (v2.y - v0.y);
    float stepX2 = (v2.x - v1.x) / (v2.y - v1.y);

    for (float y = v0.y; y <= v2.y; y++) {
        float x1 = v0.x + (y - v0.y) * stepX1;
        float x2 = v1.x + (y - v1.y) * stepX2;

        float z1 = (v0.depth + (y - v0.y) * ((v2.depth - v0.depth) / (v2.y - v0.y)));
        float z2 = (v1.depth + (y - v1.y) * ((v2.depth - v1.depth) / (v2.y - v1.y)));
        if (y < 0 || y >= HEIGHT) continue;
        if (x1 > x2) {
            std::swap(x1, x2);
            std::swap(z1, z2);
        }
        for (float x = std::max(0.0f, std::min(x1, x2)); x <= std::min(float(WIDTH - 1), std::max(x1, x2)); x++) {
            if (x < 0 || x >= WIDTH) continue;
            // Interpolate Z across the horizontal line
            float t = (x - x1) / (x2 - x1);
            float zInterpolated = (z1 + t * (z2 - z1));

            // Check and set pixel if depth test passes
            if (zInterpolated > depthBuffer[y][x]) {
                depthBuffer[y][x] = zInterpolated;
                window.setPixelColour(x, y, (255 << 24) + (uint32_t(colour.red) << 16) + (uint32_t(colour.green) << 8) + uint32_t(colour.blue));
            }
        }
    }
}
// Fill a triangle
void Draw::drawFilledTriangle(CanvasTriangle triangle, Colour colour, DrawingWindow &window, std::vector<std::vector<float>>& depthBuffer) {
    // Sort the vertices by y-coordinate (ascending)
    CanvasTriangle sorted_triangle = sortVertices(triangle);
    CanvasPoint v0 = sorted_triangle.v0();
    CanvasPoint v1 = sorted_triangle.v1();
    CanvasPoint v2 = sorted_triangle.v2();

    // Case 1: Top-flat triangle
    if (v1.y == v0.y) {
        fill_lower_triangle(v0, v1, v2, colour, window, depthBuffer);
    }
    // Case 2: Bottom-flat triangle
    else if (v2.y == v1.y) {
        fill_top_triangle(v0, v1, v2, colour, window, depthBuffer);
    }
    // Case 3: General triangle, need to split it
    else {
        // Calculate the proportion of the height difference for interpolation
        float t = (v1.y - v0.y) / (v2.y - v0.y);

        // Interpolate the x and depth coordinates to find the split point (v3)
        float v3_x = v0.x + t * (v2.x - v0.x);
        // float v3_depth = 1.0f / (1.0f / v0.depth + t * (1.0f / v2.depth - 1.0f / v0.depth));
        float v3_depth = v0.depth + t * (v2.depth - v0.depth);

        // v3 will have the same y-coordinate as v1
        CanvasPoint v3 = {v3_x, v1.y, v3_depth};

        // Now fill the two smaller triangles formed by the split
        fill_top_triangle(v0, v1, v3, colour, window, depthBuffer);
        fill_lower_triangle(v1, v3, v2, colour, window, depthBuffer);
    }
}
uint32_t getTextureColour(TexturePoint texturePoint, TextureMap &textureMap) {
    int texX = texturePoint.x * (textureMap.width - 1);
    int texY = texturePoint.y * (textureMap.height - 1);

    int index = texY * textureMap.width + texX;

    return textureMap.pixels[index];
}

void drawTexturedLine(CanvasPoint from, CanvasPoint to, TexturePoint fromTP, TexturePoint toTP, TextureMap &textureMap, DrawingWindow &window) {
    float xDiff = to.x -from.x;
    float yDiff = to.y -from.y;

    float numberOfSteps = fmax(abs(xDiff), abs(yDiff));

    std::vector<CanvasPoint> canvasPoints = Interpolate::interpolateCanvasPoints(from, to, numberOfSteps);

    std::vector<TexturePoint> texturePoints = Interpolate::interpolateTexturePoints(fromTP, toTP, numberOfSteps);

    for (float i = 0.0; i < canvasPoints.size(); i++) {
        CanvasPoint canvasPoint = canvasPoints[i];
        TexturePoint texturePoint = texturePoints[i];

        uint32_t colour = getTextureColour(texturePoint, textureMap);

        window.setPixelColour(canvasPoint.x, canvasPoint.y, colour);
    }
}

// void Draw::fillTexturedTriangle(CanvasTriangle triangle, TextureMap &textureMap, DrawingWindow &window) {
//     CanvasTriangle sorted = sortVertices(triangle);
    
//     for (int y = sorted[0].y; y <= sorted[2].y; y++) {
//         // Calculate x1, x2 for the current scanline
//         float t1, t2;
//         int x1, x2;
//         TexturePoint tp1, tp2;

//         if (y < sorted[1].y) {
//             // Interpolate along edges (v0 -> v1) and (v0 -> v2) for the upper part
//             t1 = (y - sorted[0].y) / float(sorted[1].y - sorted[0].y);
//             t2 = (y - sorted[0].y) / float(sorted[2].y - sorted[0].y);

//             x1 = sorted[0].x + t1 * (sorted[1].x - sorted[0].x);
//             x2 = sorted[0].x + t2 * (sorted[2].x - sorted[0].x);

//             tp1.x = sorted[0].texturePoint.x + t1 * (sorted[1].texturePoint.x - sorted[0].texturePoint.x);
//             tp1.y = sorted[0].texturePoint.y + t1 * (sorted[1].texturePoint.y - sorted[0].texturePoint.y);
//             tp2.x = sorted[0].texturePoint.x + t2 * (sorted[2].texturePoint.x - sorted[0].texturePoint.x);
//             tp2.y = sorted[0].texturePoint.y + t2 * (sorted[2].texturePoint.y - sorted[0].texturePoint.y);
//         } else {
//             // Interpolate along edges (v1 -> v2) and (v0 -> v2) for the lower part
//             t1 = (y - sorted[1].y) / float(sorted[2].y - sorted[1].y);
//             t2 = (y - sorted[0].y) / float(sorted[2].y - sorted[0].y);

//             x1 = sorted[1].x + t1 * (sorted[2].x - sorted[1].x);
//             x2 = sorted[0].x + t2 * (sorted[2].x - sorted[0].x);

//             tp1.x = sorted[1].texturePoint.x + t1 * (sorted[2].texturePoint.x - sorted[1].texturePoint.x);
//             tp1.y = sorted[1].texturePoint.y + t1 * (sorted[2].texturePoint.y - sorted[1].texturePoint.y);
//             tp2.x = sorted[0].texturePoint.x + t2 * (sorted[2].texturePoint.x - sorted[0].texturePoint.x);
//             tp2.y = sorted[0].texturePoint.y + t2 * (sorted[2].texturePoint.y - sorted[0].texturePoint.y);
//         }

//         // Ensure x1 is always less than x2 for left-to-right scanlines
//         if (x1 > x2) {
//             std::swap(x1, x2);
//             std::swap(tp1, tp2);
//         }

//         // Draw textured line for this scanline
//         CanvasPoint from(x1, y);
//         CanvasPoint to(x2, y);
//         drawTexturedLine(from, to, tp1, tp2, textureMap, window);
//     }
// }
uint32_t getTextureColourPerspective(float alpha, float beta, float gamma, const CanvasTriangle& triangle, TextureMap& textureMap) {
    // Get the depth values at each vertex
    float w1 = 1.0f / std::max(0.001f, triangle[0].depth);
float w2 = 1.0f / std::max(0.001f, triangle[1].depth);
float w3 = 1.0f / std::max(0.001f, triangle[2].depth);

float w = 1.0f / (alpha * w1 + beta * w2 + gamma * w3);

float u = w * (alpha * w1 * triangle[0].texturePoint.x +
               beta * w2 * triangle[1].texturePoint.x +
               gamma * w3 * triangle[2].texturePoint.x);

float v = w * (alpha * w1 * triangle[0].texturePoint.y +
               beta * w2 * triangle[1].texturePoint.y +
               gamma * w3 * triangle[2].texturePoint.y);

int texX = u * (textureMap.width - 1);
int texY = v * (textureMap.height - 1);

// Manual clamping
if (texX < 0) texX = 0;
if (texX >= textureMap.width) texX = textureMap.width - 1;

if (texY < 0) texY = 0;
if (texY >= textureMap.height) texY = textureMap.height - 1;

return textureMap.pixels[texY * textureMap.width + texX];


    
    return textureMap.pixels[texY * textureMap.width + texX];
}

void Draw::fillTexturedTriangle(CanvasTriangle triangle, TextureMap& textureMap, DrawingWindow& window) {
    // Find bounding box of the triangle
    int minX = std::floor(triangle[0].x);
    int maxX = std::ceil(triangle[0].x);
    int minY = std::floor(triangle[0].y);
    int maxY = std::ceil(triangle[0].y);
    
    // Update min/max with remaining vertices
    for(int i = 1; i < 3; i++) {
        minX = std::floor(std::min(minX, (int)triangle[i].x));
        maxX = std::ceil(std::max(maxX, (int)triangle[i].x));
        minY = std::floor(std::min(minY, (int)triangle[i].y));
        maxY = std::ceil(std::max(maxY, (int)triangle[i].y));
    }
    
    // Clip against screen boundaries
    minX = std::max(minX, 0);
    maxX = std::min(maxX, static_cast<int>(WIDTH - 1));
    minY = std::max(minY, 0);
    maxY = std::min(maxY, static_cast<int>(HEIGHT - 1));

    // Precompute some values for barycentric coordinates
    float x1 = triangle[0].x, y1 = triangle[0].y;
    float x2 = triangle[1].x, y2 = triangle[1].y;
    float x3 = triangle[2].x, y3 = triangle[2].y;
    float denominator = (y2 - y3) * (x1 - x3) + (x3 - x2) * (y1 - y3);
    
    // Iterate over each pixel in the bounding box
    for (int y = minY; y <= maxY; y++) {
        for (int x = minX; x <= maxX; x++) {
            // Calculate barycentric coordinates
            float alpha = ((y2 - y3) * (x - x3) + (x3 - x2) * (y - y3)) / denominator;
            float beta = ((y3 - y1) * (x - x3) + (x1 - x3) * (y - y3)) / denominator;
            float gamma = 1.0f - alpha - beta;
            
            // Check if point is inside triangle
            if (alpha >= 0 && beta >= 0 && gamma >= 0) {
                uint32_t colour = getTextureColourPerspective(alpha, beta, gamma, triangle, textureMap);
                window.setPixelColour(x, y, colour);
            }
        }
    }
}
uint32_t packColor(glm::vec3 color) {
    uint32_t packedColor = (255 << 24) + (int(color.x) << 16) + (int(color.y) << 8) + int(color.z);
    return packedColor;
}
void Draw::draw(DrawingWindow &window) {
	window.clearPixels();
	glm::vec3 topLeft(255, 0, 0);        // red 
	glm::vec3 topRight(0, 0, 255);       // blue 
	glm::vec3 bottomRight(0, 255, 0);    // green 
	glm::vec3 bottomLeft(255, 255, 0);   // yellow
	for (size_t y = 0; y < window.height; y++) {
		for (size_t x = 0; x < window.width; x++) {
			glm::vec3 leftColour = Interpolate::interpolateThreeElementValues(topLeft, bottomLeft, window.height)[y];
			glm::vec3 rightColour = Interpolate::interpolateThreeElementValues(topRight, bottomRight, window.height)[y];

			std::vector<glm::vec3> rowColours = Interpolate::interpolateThreeElementValues(leftColour, rightColour, window.width);
			uint32_t packedColor = packColor(rowColours[x]);
            window.setPixelColour(x, y, packedColor);
		}
	}
}

void Draw::drawPoint(DrawingWindow &window, const CanvasPoint &point, Colour colour) {
    uint32_t c = (255 << 24) + (int(colour.red) << 16) + (int(colour.green) << 8) + int(colour.blue);
    window.setPixelColour(point.x, point.y, c);
}

// Function signature for getClosestIntersection
RayTriangleIntersection getClosestIntersection(const glm::vec3& cameraPosition, const glm::vec3& rayDirection, const std::vector<ModelTriangle>& triangles) {
    size_t triangleIndex = 0;
    float closestT = std::numeric_limits<float>::max(); // Initialize with max value
    RayTriangleIntersection closestIntersection;
    bool intersectionFound = false;

    // Iterate over each triangle in the scene
    for (const ModelTriangle& triangle : triangles) {
        // Step 1: Calculate edges and the vector from the camera to the triangle vertex
        glm::vec3 e0 = triangle.vertices[1] - triangle.vertices[0];
        glm::vec3 e1 = triangle.vertices[2] - triangle.vertices[0];
        glm::vec3 SPVector = cameraPosition - triangle.vertices[0];

        // Step 2: Set up the determinant matrix
        glm::mat3 DEMatrix(-rayDirection, e0, e1);
        glm::vec3 possibleSolution = glm::inverse(DEMatrix) * SPVector;

        // Extract t, u, v from the solution
        float t = possibleSolution.x;
        float u = possibleSolution.y;
        float v = possibleSolution.z;

        // Step 3: Check if this is a valid intersection
        if ((u >= 0.0f && u <= 1.0f) &&
            (v >= 0.0f && v <= 1.0f) &&
            (u + v <= 1.0f) &&
            t > 0.0f) {
            // Step 4: Update the closest intersection if this one is closer
            if (t < closestT) {
                closestT = t;
                intersectionFound = true;

                // Step 5: Calculate the intersection point
                glm::vec3 intersectionPoint = cameraPosition + t * rayDirection;

                // Set intersection details
                closestIntersection = RayTriangleIntersection(
                    intersectionPoint, 
                    t, 
                    triangle, 
                    triangleIndex
                );
            }
        }
        triangleIndex++;
    }

    // Return the closest intersection if found; otherwise return an empty or default intersection
    if (!intersectionFound) {
        return RayTriangleIntersection();
    } 
    return closestIntersection;
}

void Draw::computeTriangleNormals(std::vector<ModelTriangle>& triangles) {
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
void Draw::computeVertexNormals(std::vector<ModelTriangle>& triangles, std::vector<std::vector<glm::vec3>>& vertexNormals) {
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

glm::vec3 calculateBarycentricCoords(const glm::vec3& p, const glm::vec3& a, const glm::vec3& b, const glm::vec3& c) {
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

void Draw::drawRayTracedScene(const glm::vec3& cameraPosition, float focalLength, float scaleFactor,
                              std::vector<ModelTriangle>& triangles, DrawingWindow &window, TextureMap& textureMap) {
    // Light and rendering parameters
    const glm::vec3 lightPosition(0.0, 0.7, 0.0);
    const float shadowBias = 0.0001f;
    const float ambientStrength = 0.2f;
    const float specularStrength = 0.5f;
    const float specularExponent = 256.0f;

    for (int y = 0; y < HEIGHT; y++) {
        for (int x = 0; x < WIDTH; x++) {
            glm::vec3 rayDirection = glm::normalize(glm::vec3(
                ((x - WIDTH / 2.0f) * (1.0f / focalLength)) / scaleFactor,
                -((y - HEIGHT / 2.0f) * (1.0f / focalLength)) / scaleFactor,
                -1.0f
            ));

            RayTriangleIntersection intersection = getClosestIntersection(cameraPosition, rayDirection, triangles);

            if (intersection.distanceFromCamera > 0) {
                // Calculate lighting
                glm::vec3 lightDir = glm::normalize(lightPosition - intersection.intersectionPoint);
                glm::vec3 normal = intersection.intersectedTriangle.normal;
                float diffuse = std::max(glm::dot(normal, lightDir), 0.0f);

                // Specular reflection
                glm::vec3 viewDir = glm::normalize(cameraPosition - intersection.intersectionPoint);
                glm::vec3 reflectDir = glm::reflect(-lightDir, normal);
                float specular = std::pow(std::max(glm::dot(viewDir, reflectDir), 0.0f), specularExponent);

                // Shadow calculation
                glm::vec3 shadowRayOrigin = intersection.intersectionPoint + normal * shadowBias;
                RayTriangleIntersection shadowIntersection = getClosestIntersection(shadowRayOrigin, lightDir, triangles);
                float shadowFactor = 1.0f;
                if (shadowIntersection.distanceFromCamera > 0 &&
                    shadowIntersection.distanceFromCamera < glm::length(lightPosition - intersection.intersectionPoint)) {
                    shadowFactor = 0.5f;
                }

                // Combine lighting components
                float lightIntensity = ambientStrength + shadowFactor * (diffuse + specularStrength * specular);

                // Clamp intensity
                lightIntensity = glm::clamp(lightIntensity, 0.0f, 1.0f);

                // Compute barycentric coordinates
                glm::vec3 A = intersection.intersectedTriangle.vertices[0];
                glm::vec3 B = intersection.intersectedTriangle.vertices[1];
                glm::vec3 C = intersection.intersectedTriangle.vertices[2];
                glm::vec3 barycentricCoords = calculateBarycentricCoords(intersection.intersectionPoint, A, B, C);

                // Construct CanvasTriangle for texture mapping
                

                bool hasValidTexture = !intersection.intersectedTriangle.texturePoints.empty() && textureMap.width > 0 && textureMap.height > 0;
                bool hasNoTexture = std::all_of(intersection.intersectedTriangle.texturePoints.begin(), intersection.intersectedTriangle.texturePoints.end(),
                                                [](const TexturePoint& tp) { return tp.x == 0.0f && tp.y == 0.0f; });
                std::cout << "texture points: ";
                for (const auto& tp : intersection.intersectedTriangle.texturePoints) {
                    std::cout << "(" << tp.x << ", " << tp.y << ") ";
                }
                std::cout << std::endl;
                uint32_t pixelColor;

                if (hasValidTexture) {
                    // Interpolate texture coordinates
                    TexturePoint texPoint;
                    texPoint.x = barycentricCoords.x * intersection.intersectedTriangle.texturePoints[0].x + 
                                 barycentricCoords.y * intersection.intersectedTriangle.texturePoints[1].x + 
                                 barycentricCoords.z * intersection.intersectedTriangle.texturePoints[2].x;
                    texPoint.y = barycentricCoords.x * intersection.intersectedTriangle.texturePoints[0].y + 
                                 barycentricCoords.y * intersection.intersectedTriangle.texturePoints[1].y + 
                                 barycentricCoords.z * intersection.intersectedTriangle.texturePoints[2].y;

                    // Clamp texture coordinates
                    texPoint.x = std::max(0.0f, std::min(texPoint.x, static_cast<float>(textureMap.width - 1)));
                    texPoint.y = std::max(0.0f, std::min(texPoint.y, static_cast<float>(textureMap.height - 1)));

                    // Get texture color
                    uint32_t texColor = getTextureColour(texPoint, textureMap);

                    float texR = ((texColor >> 16) & 0xFF) / 255.0f;
                    float texG = ((texColor >> 8) & 0xFF) / 255.0f;
                    float texB = (texColor & 0xFF) / 255.0f;

                    // Apply texture color and lighting
                    pixelColor = (255 << 24) | 
                                 (int(std::min(texR * lightIntensity, 1.0f) * 255) << 16) | 
                                 (int(std::min(texG * lightIntensity, 1.0f) * 255) << 8) | 
                                 int(std::min(texB * lightIntensity, 1.0f) * 255);
                } else if (hasNoTexture) {
                    std::cout << "No texture found for triangle " << std::endl;
                    const Colour& triangleColor = intersection.intersectedTriangle.colour;
                    pixelColor = (255 << 24) |
                                 ((int(std::min(triangleColor.red * lightIntensity, 255.0f))) << 16) |
                                 ((int(std::min(triangleColor.green * lightIntensity, 255.0f))) << 8) |
                                 (int(std::min(triangleColor.blue * lightIntensity, 255.0f)));
                } else {
                    // If no texture and no color, use a default color (e.g., black)
                    pixelColor = 0; // Or use some other fallback
                }

                // Set the pixel color on the window
                window.setPixelColour(x, y, pixelColor);
            } else {
                // Background color (e.g., black if no intersection)
                window.setPixelColour(x, y, 0);
            }
        }
    }
}





// void Draw::drawRayTracedScene(const glm::vec3& cameraPosition, float focalLength, float scaleFactor, 
//                              std::vector<ModelTriangle>& triangles, DrawingWindow &window) {
//     // Light position from the model
//     const glm::vec3 lightPosition(0.0, 0.7, 0.0);
//     const float shadowBias = 0.0001f;
//     const float shadowSoftness = 0.4f;
//     // Enhanced ambient lighting parameters
//     const float ambientStrength = 0.2f;     // Base ambient light level
//     const float ambientOcclusion = 0.3f;    // Additional ambient in corners/crevices
//     const float minLightThreshold = 0.15f;  // Minimum light level for any surface
    
//     // Specular parameters
//     const float specularStrength = 0.5f;
//     const float specularExponent = 256.0f;
    
//     for (int y = 0; y < HEIGHT; y++) {
//         for (int x = 0; x < WIDTH; x++) {
//             glm::vec3 rayDirection = glm::normalize(glm::vec3(
//                 ((x - WIDTH / 2.0f) * (1.0f / focalLength)) / scaleFactor,
//                 -((y - HEIGHT / 2.0f) * (1.0f / focalLength)) / scaleFactor,
//                 -1.0f
//             ));

//             RayTriangleIntersection intersection = getClosestIntersection(cameraPosition, rayDirection, triangles);

//             if (intersection.distanceFromCamera > 0) {
//                 // Calculate lighting vectors
//                 glm::vec3 lightDir = lightPosition - intersection.intersectionPoint;
//                 float distanceToLight = glm::length(lightDir);
//                 lightDir = glm::normalize(lightDir);
                

//                 glm::vec3& normal = intersection.intersectedTriangle.normal;

//                 // Calculate diffuse lighting
//                 float diffuse = std::max(glm::dot(normal, lightDir), 0.0f);
                
//                 // Enhanced ambient lighting calculation
//                 // float ambientFactor = ambientStrength;
                
//                 // // Add extra ambient light in corners and crevices
//                 // // This simulates light bouncing in confined spaces
//                 // float occlusionFactor = 1.0f - std::abs(glm::dot(normal, glm::vec3(0, 1, 0)));
//                 // ambientFactor += ambientOcclusion * occlusionFactor;

//                 // Enhanced ambient lighting calculation
//                 float ambientFactor = ambientStrength; // Base ambient light

//                 // Optional: Refine occlusionFactor based on geometry or light direction
//                 float occlusionFactor = 1.0f - glm::max(glm::dot(normal, glm::vec3(0, 1, 0)), 0.0f); // Simulated occlusion
//                 ambientFactor *= (1.0f - ambientOcclusion * occlusionFactor);
                
//                 // Calculate proximity lighting (inverse square law)
//                 float proximityFactor = 1.0f / (4.0f * M_PI * distanceToLight * distanceToLight);
//                 proximityFactor = std::min(proximityFactor, 1.0f);
                
//                 // Specular calculation
//                 glm::vec3 viewDir = glm::normalize(cameraPosition - intersection.intersectionPoint);
//                 glm::vec3 reflectionDir = glm::normalize(lightDir - 2.0f * glm::dot(normal, lightDir) * normal);
//                 float specular = std::pow(std::max(glm::dot(viewDir, reflectionDir), 0.0f), specularExponent);
                
//                 // Shadow calculation
//                 glm::vec3 shadowRayOrigin = intersection.intersectionPoint + normal * shadowBias;
//                 RayTriangleIntersection shadowIntersection = getClosestIntersection(shadowRayOrigin, lightDir, triangles);

//                 float shadowFactor = 1.0f;
//                 if (shadowIntersection.distanceFromCamera > 0 && 
//                     shadowIntersection.distanceFromCamera < distanceToLight) {
//                     shadowFactor = shadowSoftness;
//                 }

//                 // Combine all lighting components
//                 float lightIntensity = (proximityFactor * shadowFactor + diffuse * shadowFactor) * (1.0f - ambientFactor) + ambientFactor + specularStrength *specular * shadowFactor;
                
                
//                 // Apply minimum light threshold to simulate indirect bounces
//                 lightIntensity = std::max(lightIntensity, minLightThreshold);
                
//                 // Ensure final intensity doesn't exceed 1.0
//                 lightIntensity = std::min(lightIntensity, 1.0f);

//                 // Apply color with all lighting effects
//                 const Colour& triangleColor = intersection.intersectedTriangle.colour;
//                 uint32_t pixelColor = (255 << 24) | 
//                     ((int(triangleColor.red * lightIntensity) & 0xFF) << 16) | 
//                     ((int(triangleColor.green * lightIntensity) & 0xFF) << 8) | 
//                     (int(triangleColor.blue * lightIntensity) & 0xFF);
                
//                 window.setPixelColour(x, y, pixelColor);
//             } else {
//                 window.setPixelColour(x, y, 0);
//             }
//         }
//     }
// }

// void Draw::drawRayTracedScene(const glm::vec3& cameraPosition, float focalLength, float scaleFactor,
//                               std::vector<ModelTriangle>& triangles, DrawingWindow &window) {
//     // Reflectivity map: Associate reflectivity with specific triangles
//     std::map<int, float> reflectiveTriangles = {
//         {0, 1.0f}, // Triangle 0 is fully reflective
//         {1, 0.5f}, // Triangle 1 is partially reflective
//         // Add more entries as needed
//     };

//     const glm::vec3 lightPosition(0.0, 0.7, 0.0);
//     const int maxReflectionDepth = 3; // Maximum recursion depth for reflections

//     for (int y = 0; y < HEIGHT; y++) {
//         for (int x = 0; x < WIDTH; x++) {
//             // Calculate primary ray direction
//             glm::vec3 rayDirection = glm::normalize(glm::vec3(
//                 ((x - WIDTH / 2.0f) * (1.0f / focalLength)) / scaleFactor,
//                 -((y - HEIGHT / 2.0f) * (1.0f / focalLength)) / scaleFactor,
//                 -1.0f
//             ));

//             // Trace the ray and get the final color
//             glm::vec3 pixelColor = traceRay(cameraPosition, rayDirection, triangles, lightPosition, 0, maxReflectionDepth, reflectiveTriangles);

//             // Convert to 8-bit color and set the pixel
//             uint32_t finalColor = (255 << 24) | 
//                 ((int(pixelColor.r * 255.0f) & 0xFF) << 16) |
//                 ((int(pixelColor.g * 255.0f) & 0xFF) << 8) |
//                 (int(pixelColor.b * 255.0f) & 0xFF);

//             window.setPixelColour(x, y, finalColor);
//         }
//     }
// }


// glm::vec3 Draw::traceRay(const glm::vec3& origin, const glm::vec3& direction, 
//                          const std::vector<ModelTriangle>& triangles, const glm::vec3& lightPosition, 
//                          int depth, int maxDepth, const std::map<int, float>& reflectiveTriangles) {
//     if (depth > maxDepth) return glm::vec3(0.0f); // Base case for recursion

//     RayTriangleIntersection intersection = getClosestIntersection(origin, direction, triangles);
//     if (intersection.distanceFromCamera < 0) return glm::vec3(0.0f); // No intersection

//     glm::vec3 normal = intersection.intersectedTriangle.normal;
//     glm::vec3 intersectionPoint = intersection.intersectionPoint;
//     const Colour& triangleColor = intersection.intersectedTriangle.colour;

//     // Lighting calculations (as before)
//     glm::vec3 lightDir = glm::normalize(lightPosition - intersectionPoint);
//     float diffuse = std::max(glm::dot(normal, lightDir), 0.0f);

//     glm::vec3 viewDir = glm::normalize(origin - intersectionPoint);
//     glm::vec3 reflectionDir = glm::normalize(lightDir - 2.0f * glm::dot(normal, lightDir) * normal);
//     float specular = std::pow(std::max(glm::dot(viewDir, reflectionDir), 0.0f), 256.0f);

//     float ambientFactor = 0.2f;
//     float occlusionFactor = 1.0f - glm::max(glm::dot(normal, glm::vec3(0, 1, 0)), 0.0f);
//     ambientFactor *= (1.0f - 0.3f * occlusionFactor);
//     float distanceToLight = glm::length(lightPosition - intersectionPoint);
//     float proximityFactor = 1.0f / (4.0f * M_PI * distanceToLight * distanceToLight);
//     proximityFactor = std::min(proximityFactor, 1.0f);

//     float lightIntensity = (proximityFactor + diffuse) * (1.0f - ambientFactor) + ambientFactor + 0.5f * specular;

//     glm::vec3 surfaceColor = glm::vec3(
//         triangleColor.red * lightIntensity,
//         triangleColor.green * lightIntensity,
//         triangleColor.blue * lightIntensity
//     ) / 255.0f;

//     // Reflection handling
//     int triangleIndex = intersection.triangleIndex; // Assuming you have an index in the intersection
//     if (reflectiveTriangles.count(triangleIndex) > 0 && depth < maxDepth) {
//         float reflectionIntensity = reflectiveTriangles.at(triangleIndex); // Get reflection value
//         glm::vec3 reflectionDir = glm::normalize(direction - 2.0f * glm::dot(direction, normal) * normal);
//         glm::vec3 reflectionColor = Draw::traceRay(intersectionPoint + normal * 0.0001f, reflectionDir, triangles, lightPosition, depth + 1, maxDepth, reflectiveTriangles);
//         surfaceColor = glm::mix(surfaceColor, reflectionColor, reflectionIntensity);
//     }

//     return glm::clamp(surfaceColor, 0.0f, 1.0f);
// }

void Draw::drawRayTracedSceneGouraud(const glm::vec3& cameraPosition, float focalLength, float scaleFactor, 
                                     std::vector<ModelTriangle>& triangles, std::vector<std::vector<glm::vec3>>& vertexNormals, 
                                     DrawingWindow &window, TextureMap &textureMap) {
    // Keep existing lighting constants
    const glm::vec3 lightPosition(0.0, 0.7, 0.6);
    const float shadowBias = 0.0001f;
    const float shadowSoftness = 0.6f;
    const float ambientStrength = 0.1f;
    const float ambientOcclusion = 0.2f;
    const float minLightThreshold = 0.05f;
    const float specularStrength = 0.5f;
    const float specularExponent = 256.0f;

    for (int y = 0; y < HEIGHT; y++) {
        for (int x = 0; x < WIDTH; x++) {
            glm::vec3 rayDirection = glm::normalize(glm::vec3(
                ((x - WIDTH / 2.0f) * (1.0f / focalLength)) / scaleFactor,
                -((y - HEIGHT / 2.0f) * (1.0f / focalLength)) / scaleFactor,
                -1.0f
            ));

            RayTriangleIntersection intersection = getClosestIntersection(cameraPosition, rayDirection, triangles);

            if (intersection.distanceFromCamera > 0) {
                const ModelTriangle& triangle = intersection.intersectedTriangle;
                const glm::vec3& p = intersection.intersectionPoint;
                const glm::vec3& a = triangle.vertices[0];
                const glm::vec3& b = triangle.vertices[1];
                const glm::vec3& c = triangle.vertices[2];

                // Calculate barycentric coordinates using the new function
                glm::vec3 barycentricCoords = calculateBarycentricCoords(p, a, b, c);

                if (barycentricCoords.x < 0.0f) {
                    continue;  // Skip degenerate triangle
                }

                float u = barycentricCoords.x;
                float v = barycentricCoords.y;
                float w = barycentricCoords.z;

                // Calculate vertex intensities (same as before)
                float vertexIntensities[3];
                for(int i = 0; i < 3; i++) {
                    glm::vec3 normal = vertexNormals[intersection.triangleIndex][i];
                    glm::vec3 vertexPos = triangle.vertices[i];
                    glm::vec3 lightDir = lightPosition - vertexPos;
                    float distanceToLight = glm::length(lightDir);
                    lightDir = glm::normalize(lightDir);

                    float diffuse = 0.0f;
                    float ambientFactor = ambientStrength;
                    float occlusionFactor = 1.0f - std::abs(glm::dot(normal, glm::vec3(0, 1, 0)));
                    ambientFactor += ambientOcclusion * occlusionFactor;
                    float proximityFactor = 1.0f / (4.0f * M_PI * distanceToLight * distanceToLight);
                    proximityFactor = std::min(proximityFactor, 1.0f);
                    float specular = 0.0f;

                    diffuse = std::max(glm::dot(normal, lightDir), 0.0f);

                    if (diffuse > 0.0f) {
                        glm::vec3 viewDir = glm::normalize(cameraPosition - vertexPos);
                        glm::vec3 reflectionDir = glm::normalize(lightDir - 2.0f * glm::dot(normal, lightDir) * normal);
                        specular = std::pow(std::max(glm::dot(viewDir, reflectionDir), 0.0f), specularExponent);
                    }

                    glm::vec3 shadowRayOrigin = vertexPos + normal * shadowBias;
                    RayTriangleIntersection shadowIntersection = getClosestIntersection(shadowRayOrigin, lightDir, triangles);

                    float shadowFactor = 1.0f;
                    if (shadowIntersection.distanceFromCamera > 0 && 
                        shadowIntersection.distanceFromCamera < distanceToLight) {
                        shadowFactor = shadowSoftness;
                    }

                    float diffuseIntensity = (proximityFactor * shadowFactor + diffuse * shadowFactor) * (1.0f - ambientFactor) + ambientFactor;
                    float vertexIntensity = diffuseIntensity + specularStrength * specular * shadowFactor;
                    vertexIntensity = std::max(vertexIntensity, minLightThreshold);
                    vertexIntensity = std::min(vertexIntensity, 1.0f);
                    
                    vertexIntensities[i] = vertexIntensity;
                }

                float interpolatedIntensity = u * vertexIntensities[0] + 
                                            v * vertexIntensities[1] + 
                                            w * vertexIntensities[2];

                // Check for valid texture map
                bool hasValidTexture = !triangle.texturePoints.empty() && textureMap.width > 0 && textureMap.height > 0;
                bool hasNoTexture = std::any_of(triangle.texturePoints.begin(), triangle.texturePoints.end(),
                                 [](const TexturePoint& tp) { return tp.x == 0.0f && tp.y == 0.0f; });

                uint32_t pixelColor;
                if (hasValidTexture) {
                    // Interpolate texture coordinates
                    TexturePoint texPoint;
                    texPoint.x = u * triangle.texturePoints[0].x + 
                                v * triangle.texturePoints[1].x + 
                                w * triangle.texturePoints[2].x;
                    texPoint.y = u * triangle.texturePoints[0].y + 
                                v * triangle.texturePoints[1].y + 
                                w * triangle.texturePoints[2].y;

                    // Clamp texture coordinates
                    texPoint.x = std::max(0.0f, std::min(texPoint.x, static_cast<float>(textureMap.width - 1)));
                    texPoint.y = std::max(0.0f, std::min(texPoint.y, static_cast<float>(textureMap.height - 1)));

                    uint32_t texColor;
                    texColor = getTextureColour(texPoint, textureMap);

                    if (hasValidTexture) {
                        float texR = ((texColor >> 16) & 0xFF) / 255.0f;
                        float texG = ((texColor >> 8) & 0xFF) / 255.0f;
                        float texB = (texColor & 0xFF) / 255.0f;
                        pixelColor = (255 << 24) | 
                                    (int(std::min(texR * interpolatedIntensity, 1.0f) * 255) << 16) | 
                                    (int(std::min(texG * interpolatedIntensity, 1.0f) * 255) << 8) | 
                                    int(std::min(texB * interpolatedIntensity, 1.0f) * 255);
                    }
                }
                // std::cout << "before !!!" << std::endl;
                if (hasNoTexture) {
                    const Colour& triangleColor = triangle.colour;
                    pixelColor = (255 << 24) |
                                 ((int(std::min(triangleColor.red * interpolatedIntensity, 255.0f))) << 16) |
                                 ((int(std::min(triangleColor.green * interpolatedIntensity, 255.0f))) << 8) |
                                 (int(std::min(triangleColor.blue * interpolatedIntensity, 255.0f)));
                }
                window.setPixelColour(x, y, pixelColor);
            }
        }
    }
}



// wihtout texture mapping
// void Draw::drawRayTracedSceneGouraud(const glm::vec3& cameraPosition, float focalLength, float scaleFactor, 
//                                      std::vector<ModelTriangle>& triangles, std::vector<std::vector<glm::vec3>>& vertexNormals, 
//                                      DrawingWindow &window) {
//     // Light position from the model
//     const glm::vec3 lightPosition(0.0, 0.7, 0.6);
//     const float shadowBias = 0.0001f;
//     const float shadowSoftness = 0.6f;
    
//     // Enhanced ambient lighting parameters
//     const float ambientStrength = 0.1f;
//     const float ambientOcclusion = 0.2f;
//     const float minLightThreshold = 0.05f;
    
//     // Specular parameters
//     const float specularStrength = 0.5f;
//     const float specularExponent = 256.0f;
    
//     for (int y = 0; y < HEIGHT; y++) {
//         for (int x = 0; x < WIDTH; x++) {
//             glm::vec3 rayDirection = glm::normalize(glm::vec3(
//                 ((x - WIDTH / 2.0f) * (1.0f / focalLength)) / scaleFactor,
//                 -((y - HEIGHT / 2.0f) * (1.0f / focalLength)) / scaleFactor,
//                 -1.0f
//             ));

//             RayTriangleIntersection intersection = getClosestIntersection(cameraPosition, rayDirection, triangles);

//             if (intersection.distanceFromCamera > 0) {
//                 // Calculate barycentric coordinates for the intersection point
//                 const ModelTriangle& triangle = intersection.intersectedTriangle;
//                 const glm::vec3& p = intersection.intersectionPoint;
//                 const glm::vec3& a = triangle.vertices[0];
//                 const glm::vec3& b = triangle.vertices[1];
//                 const glm::vec3& c = triangle.vertices[2];

//                 // Calculate vectors
//                 glm::vec3 v0 = b - a;
//                 glm::vec3 v1 = c - a;
//                 glm::vec3 v2 = p - a;

//                 // Calculate dot products
//                 float d00 = glm::dot(v0, v0);
//                 float d01 = glm::dot(v0, v1);
//                 float d11 = glm::dot(v1, v1);
//                 float d20 = glm::dot(v2, v0);
//                 float d21 = glm::dot(v2, v1);

//                 // Calculate barycentric coordinates
//                 float denom = d00 * d11 - d01 * d01;
//                 float v = (d11 * d20 - d01 * d21) / denom;
//                 float w = (d00 * d21 - d01 * d20) / denom;
//                 float u = 1.0f - v - w;

//                 // Calculate lighting using vertex normals
//                 float vertexIntensities[3];
//                 for(int i = 0; i < 3; i++) {
//                     // Get the vertex normal from the vertexNormals array
//                     glm::vec3 normal = vertexNormals[intersection.triangleIndex][i];
                    
//                     // Lighting calculations remain the same as before
//                     glm::vec3 vertexPos = triangle.vertices[i];
//                     glm::vec3 lightDir = lightPosition - vertexPos;
//                     float distanceToLight = glm::length(lightDir);
//                     lightDir = glm::normalize(lightDir);

//                     float diffuse = 0.0f;
//                     float ambientFactor = ambientStrength;
//                     float occlusionFactor = 1.0f - std::abs(glm::dot(normal, glm::vec3(0, 1, 0)));
//                     ambientFactor += ambientOcclusion * occlusionFactor;
//                     float proximityFactor = 1.0f / (4.0f * M_PI * distanceToLight * distanceToLight);
//                     proximityFactor = std::min(proximityFactor, 1.0f);
//                     float specular = 0.0f;

//                     diffuse = std::max(glm::dot(normal, lightDir), 0.0f);

//                     if (diffuse > 0.0f) {
//                         glm::vec3 viewDir = glm::normalize(cameraPosition - vertexPos);
//                         glm::vec3 reflectionDir = glm::normalize(lightDir - 2.0f * glm::dot(normal, lightDir) * normal);
//                         specular = std::pow(std::max(glm::dot(viewDir, reflectionDir), 0.0f), specularExponent);
//                     }
                    
//                     glm::vec3 shadowRayOrigin = vertexPos + normal * shadowBias;
//                     RayTriangleIntersection shadowIntersection = getClosestIntersection(shadowRayOrigin, lightDir, triangles);

//                     float shadowFactor = 1.0f;
//                     if (shadowIntersection.distanceFromCamera > 0 && 
//                         shadowIntersection.distanceFromCamera < distanceToLight) {
//                         shadowFactor = shadowSoftness;
//                     }

//                     float diffuseIntensity = (proximityFactor * shadowFactor + diffuse * shadowFactor) * (1.0f - ambientFactor) + ambientFactor;
//                     float vertexIntensity = diffuseIntensity + specularStrength * specular * shadowFactor;
//                     vertexIntensity = std::max(vertexIntensity, minLightThreshold);
//                     vertexIntensity = std::min(vertexIntensity, 1.0f);
                    
//                     vertexIntensities[i] = vertexIntensity;
//                 }

//                 // Interpolate the vertex intensities using calculated barycentric coordinates
//                 float interpolatedIntensity = u * vertexIntensities[0] + 
//                                             v * vertexIntensities[1] + 
//                                             w * vertexIntensities[2];

//                 // Apply color with interpolated lighting
//                 const Colour& triangleColor = intersection.intersectedTriangle.colour;
//                 uint32_t pixelColor = (255 << 24) | 
//                     ((int(triangleColor.red * interpolatedIntensity) & 0xFF) << 16) | 
//                     ((int(triangleColor.green * interpolatedIntensity) & 0xFF) << 8) | 
//                     (int(triangleColor.blue * interpolatedIntensity) & 0xFF);
                
//                 window.setPixelColour(x, y, pixelColor);
//             } else {
//                 window.setPixelColour(x, y, 0);
//             }
//         }
//     }
// }