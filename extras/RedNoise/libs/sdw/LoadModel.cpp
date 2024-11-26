#include "LoadModel.h"
#include <fstream>
#include <sstream>
#include <iostream>
#include "Colour.h"
#include "ModelTriangle.h"
#include "Utils.h"
#include "Draw.h"

std::unordered_map<std::string, Colour> LoadModel::loadMaterials(
    const std::string& filename, 
    std::unordered_map<std::string, TextureMap*>& textureMaps) 
{
    std::unordered_map<std::string, Colour> colourPalette;
    std::ifstream file(filename);
    std::string line;

    if (!file.is_open()) {
        std::cerr << "Could not open materials file: " << filename << std::endl;
        return colourPalette;
    }

    std::string currentMaterial;

    while (std::getline(file, line)) {
        std::istringstream iss(line);
        std::string command;
        iss >> command;

        if (command == "newmtl") {  // Define a new material
            iss >> currentMaterial;
        } else if (command == "Kd") {  // Diffuse color
            float r, g, b;
            if (iss >> r >> g >> b) {
                // Convert float (0-1) to int (0-255)
                int red = static_cast<int>(r * 255);
                int green = static_cast<int>(g * 255);
                int blue = static_cast<int>(b * 255);
                colourPalette[currentMaterial] = Colour(currentMaterial, red, green, blue);
            }
        } else if (command == "map_Kd") {  // Texture map
            std::string textureFile;
            if (iss >> textureFile) {
                textureMaps[currentMaterial] = new TextureMap(textureFile);  // Load the texture
            }
        }
    }

    file.close();
    return colourPalette;
}


std::vector<ModelTriangle> LoadModel::loadOBJ(const std::string &filePath, const std::string &mtlPath, float scale) {
    std::vector<glm::vec3> vertices;
    std::vector<TexturePoint> textureCoords;
    std::vector<ModelTriangle> triangles;
    std::unordered_map<std::string, TextureMap*> textureMaps;
    std::unordered_map<std::string, Colour> colourPalette = loadMaterials(mtlPath, textureMaps);

    // std::unordered_map<std::string, Colour> colourPalette = loadMaterials(mtlPath);

    std::ifstream objFile(filePath);
    if (!objFile.is_open()) {
        std::cerr << "Error: Could not open file " << filePath << std::endl;
        return triangles;
    }

    Colour defaultColour = Colour("Default", 255, 255, 255);
    std::string currentMaterial;
    std::string line;

    while (getline(objFile, line)) {
        if (line.empty()) continue;

        std::vector<std::string> tokens = split(line, ' ');

        if (tokens[0] == "v") {  // Vertex
            float x = std::stof(tokens[1]) * scale;
            float y = std::stof(tokens[2]) * scale;
            float z = std::stof(tokens[3]) * scale;
            vertices.push_back(glm::vec3(x, y, z));
        } else if (tokens[0] == "vt") {  // Texture coordinate
            float u = std::stof(tokens[1]);
            float v = std::stof(tokens[2]);
            textureCoords.emplace_back(u, v);
        } else if (tokens[0] == "usemtl") {
            currentMaterial = tokens[1];
        } else if (tokens[0] == "f") {  // Face
            std::vector<int> vertexIndices;
            std::vector<int> textureIndices;

            for (size_t i = 1; i < tokens.size(); i++) {
                std::vector<std::string> indices = split(tokens[i], '/');
                vertexIndices.push_back(std::stoi(indices[0]) - 1);  // Vertex index

                // Optional texture index
                if (indices.size() > 1 && !indices[1].empty()) {
                    textureIndices.push_back(std::stoi(indices[1]) - 1);
                }
            }

            if (vertexIndices.size() >= 3) {
                // Create triangle
                ModelTriangle triangle(
                    vertices[vertexIndices[0]],
                    vertices[vertexIndices[1]],
                    vertices[vertexIndices[2]],
                    colourPalette.count(currentMaterial) ? colourPalette[currentMaterial] : defaultColour
                );

                // Assign texture coordinates if available
                if (textureIndices.size() >= 3) {
                    triangle.texturePoints[0] = textureCoords[textureIndices[0]];
                    triangle.texturePoints[1] = textureCoords[textureIndices[1]];
                    triangle.texturePoints[2] = textureCoords[textureIndices[2]];
                }

                triangles.push_back(triangle);
            } else {
                std::cerr << "Warning: Face with less than 3 vertices encountered in " << filePath << std::endl;
            }
        }
    }

    objFile.close();
    return triangles;
}