#ifndef LOADMODEL_H
#define LOADMODEL_H

#include <unordered_map>
#include <string>
#include <vector>
#include <glm/vec3.hpp>
#include "Colour.h"
#include "ModelTriangle.h"
#include "TextureMap.h"

class LoadModel {
public:
    static std::unordered_map<std::string, Colour> loadMaterials(const std::string& filename, std::unordered_map<std::string, TextureMap*>& textureMaps);
    static std::vector<ModelTriangle> loadOBJ(const std::string& filePath, const std::string& mtlPath, float scale);
};

#endif // LOADMODEL_H
