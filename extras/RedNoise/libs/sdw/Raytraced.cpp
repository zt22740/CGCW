#include "Raytraced.h"

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

// Calculate diffuse lighting
float calculateDiffuse(const glm::vec3& normal, const glm::vec3& lightDir) {
    return std::max(glm::dot(normal, lightDir), 0.0f);
}

// Calculate ambient lighting with enhanced occlusion effect
float calculateAmbient(const glm::vec3& normal, float ambientStrength, float ambientOcclusion) {
    float ambientFactor = ambientStrength; // Base ambient light
    float occlusionFactor = 1.0f - glm::max(glm::dot(normal, glm::vec3(0, 1, 0)), 0.0f); // Simulated occlusion
    ambientFactor *= (1.0f - ambientOcclusion * occlusionFactor);
    return ambientFactor;
}

// Calculate proximity lighting based on the inverse square law
float calculateProximityLighting(float distanceToLight) {
    float proximityFactor = 1.0f / (4.0f * M_PI * distanceToLight * distanceToLight);
    return std::min(proximityFactor, 1.0f);
}

// Calculate specular lighting based on view direction and reflection direction
float calculateSpecular(const glm::vec3& viewDir, const glm::vec3& lightDir, const glm::vec3& normal, float specularExponent) {
    glm::vec3 reflectionDir = glm::normalize(lightDir - 2.0f * glm::dot(normal, lightDir) * normal);
    return std::pow(std::max(glm::dot(viewDir, reflectionDir), 0.0f), specularExponent);
}

// Calculate shadow factor based on intersection with a shadow ray
float calculateShadow(const glm::vec3& shadowRayOrigin, const glm::vec3& lightDir, const std::vector<ModelTriangle>& triangles, float shadowSoftness, float distanceToLight) {
    RayTriangleIntersection shadowIntersection = getClosestIntersection(shadowRayOrigin, lightDir, triangles);
    if (shadowIntersection.distanceFromCamera > 0 && shadowIntersection.distanceFromCamera < distanceToLight) {
        return shadowSoftness;
    }
    return 1.0f; // No shadow
}

void Raytraced::drawRayTracedScene(const glm::vec3& cameraPosition, float focalLength, float scaleFactor,
                              std::vector<ModelTriangle>& triangles, DrawingWindow &window, TextureMap& textureMap) {
    // Light and light position from the model
    const glm::vec3 lightPosition(0.0, 0.7, 0.0);
    const float shadowBias = 0.0001f;
    const float shadowSoftness = 0.4f;
    const float ambientStrength = 0.2f;     // Base ambient light level
    const float ambientOcclusion = 0.3f;    // Additional ambient in corners/crevices
    const float minLightThreshold = 0.15f;   // Minimum light intensity
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
                // Lighting vectors
                glm::vec3 lightDir = lightPosition - intersection.intersectionPoint;
                float distanceToLight = glm::length(lightDir);
                lightDir = glm::normalize(lightDir);
                glm::vec3& normal = intersection.intersectedTriangle.normal;

                // Calculate lighting components
                float diffuse = calculateDiffuse(normal, lightDir);
                float ambientFactor = calculateAmbient(normal, ambientStrength, ambientOcclusion);
                float proximityFactor = calculateProximityLighting(distanceToLight);
                glm::vec3 viewDir = glm::normalize(cameraPosition - intersection.intersectionPoint);
                float specular = calculateSpecular(viewDir, lightDir, normal, specularExponent);

                // Shadow calculation
                glm::vec3 shadowRayOrigin = intersection.intersectionPoint + normal * shadowBias;
                float shadowFactor = calculateShadow(shadowRayOrigin, lightDir, triangles, shadowSoftness, distanceToLight);

                // Combine all lighting components
                float lightIntensity = (proximityFactor * shadowFactor + diffuse * shadowFactor) * (1.0f - ambientFactor) + ambientFactor + specularStrength * specular * shadowFactor;

                // Apply minimum light threshold to simulate indirect bounces
                lightIntensity = std::max(lightIntensity, minLightThreshold);
                
                // Ensure final intensity doesn't exceed 1.0
                lightIntensity = std::min(lightIntensity, 1.0f);

                // Compute barycentric coordinates
                glm::vec3 A = intersection.intersectedTriangle.vertices[0];
                glm::vec3 B = intersection.intersectedTriangle.vertices[1];
                glm::vec3 C = intersection.intersectedTriangle.vertices[2];
                glm::vec3 barycentricCoords = Calculations::calculateBarycentricCoords(intersection.intersectionPoint, A, B, C);

                bool hasValidTexture = !intersection.intersectedTriangle.texturePoints.empty() && textureMap.width > 0 && textureMap.height > 0;
                bool hasNoTexture = std::any_of(intersection.intersectedTriangle.texturePoints.begin(), intersection.intersectedTriangle.texturePoints.end(),
                                                [](const TexturePoint& tp) { return tp.x == 0.0f && tp.y == 0.0f; });

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
                    uint32_t texColor = Draw::getTextureColour(texPoint, textureMap);

                    float texR = ((texColor >> 16) & 0xFF) / 255.0f;
                    float texG = ((texColor >> 8) & 0xFF) / 255.0f;
                    float texB = (texColor & 0xFF) / 255.0f;

                    // Apply texture color and lighting
                    pixelColor = (255 << 24) | 
                                 (int(std::min(texR * lightIntensity, 1.0f) * 255) << 16) | 
                                 (int(std::min(texG * lightIntensity, 1.0f) * 255) << 8) | 
                                 int(std::min(texB * lightIntensity, 1.0f) * 255);
                } 
                if (hasNoTexture) {
                    const Colour& triangleColor = intersection.intersectedTriangle.colour;
                    pixelColor = (255 << 24) |
                                 ((int(std::min(triangleColor.red * lightIntensity, 255.0f))) << 16) |
                                 ((int(std::min(triangleColor.green * lightIntensity, 255.0f))) << 8) |
                                 (int(std::min(triangleColor.blue * lightIntensity, 255.0f)));
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


void Raytraced::drawRayTracedSceneGouraud(const glm::vec3& cameraPosition, float focalLength, float scaleFactor, 
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
                glm::vec3 barycentricCoords = Calculations::calculateBarycentricCoords(p, a, b, c);

                if (barycentricCoords.x < 0.0f) {
                    continue;  // Skip degenerate triangle
                }

                float u = barycentricCoords.x;
                float v = barycentricCoords.y;
                float w = barycentricCoords.z;

                // Calculate vertex intensities (same as before)
                float vertexIntensities[3];
                for (int i = 0; i < 3; i++) {
                    glm::vec3 normal = vertexNormals[intersection.triangleIndex][i];
                    glm::vec3 vertexPos = triangle.vertices[i];
                    glm::vec3 lightDir = lightPosition - vertexPos;
                    float distanceToLight = glm::length(lightDir);
                    lightDir = glm::normalize(lightDir);

                    // Calculate lighting components using the new functions
                    float diffuse = calculateDiffuse(normal, lightDir);
                    float ambientFactor = calculateAmbient(normal, ambientStrength, ambientOcclusion);
                    float proximityFactor = calculateProximityLighting(distanceToLight);

                    glm::vec3 viewDir = glm::normalize(cameraPosition - vertexPos);
                    float specular = calculateSpecular(viewDir, lightDir, normal, specularExponent);

                    glm::vec3 shadowRayOrigin = vertexPos + normal * shadowBias;
                    float shadowFactor = calculateShadow(shadowRayOrigin, lightDir, triangles, shadowSoftness, distanceToLight);

                    // Combine the lighting components
                    float diffuseIntensity = (proximityFactor * shadowFactor + diffuse * shadowFactor) * (1.0f - ambientFactor) + ambientFactor;
                    float vertexIntensity = diffuseIntensity + specularStrength * specular * shadowFactor;
                    vertexIntensity = std::max(vertexIntensity, minLightThreshold);
                    vertexIntensity = std::min(vertexIntensity, 1.0f);
                    
                    vertexIntensities[i] = vertexIntensity;
                }

                // Interpolate vertex intensities using barycentric coordinates
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
                    texColor = Draw::getTextureColour(texPoint, textureMap);

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
