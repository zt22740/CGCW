#include "Draw.h"
#include <unordered_map>
#include <map>
#include <algorithm>
#include <iostream>
#include <bits/algorithmfwd.h>

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
uint32_t Draw::getTextureColour(TexturePoint texturePoint, TextureMap &textureMap) {
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

        uint32_t colour = Draw::getTextureColour(texturePoint, textureMap);

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