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
#include "TexturePoint.h"

#define WIDTH 320
#define HEIGHT 240
class Draw {
public:
    static uint32_t getTextureColour(TexturePoint texturePoint, TextureMap &textureMap);
    static void drawLine(CanvasPoint from, CanvasPoint to, Colour colour, DrawingWindow &window);
    static void drawTriangle(const CanvasTriangle triangle, const Colour colour, DrawingWindow &window);
    static void drawFilledTriangle(CanvasTriangle triangle, Colour colour, DrawingWindow &window, std::vector<std::vector<float>>& depthBuffer);
    static void fillTexturedTriangle(CanvasTriangle triangle, TextureMap &textureMap, DrawingWindow &window);
    static void draw(DrawingWindow &window);
    static void drawPoint(DrawingWindow &window, const CanvasPoint &point, Colour colour);
};

#endif // DRAW_H
