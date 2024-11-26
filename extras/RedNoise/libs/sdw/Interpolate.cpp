#include <vector>
std::vector<float> interpolateSingleFloats(float from, float to, int numberOfValues) {
    // Create a vector to store the result
    std::vector<float> result;
    
    // Calculate the step size
    float step = (to - from) / (numberOfValues - 1);
    
    // Fill the vector with evenly spaced values
    for (int i = 0; i < numberOfValues; i++) {
        result.push_back(from + i * step);
    }
    
    return result;
}

std::vector<glm::vec3> interpolateThreeElementValues(glm::vec3 from, glm::vec3 to, int numberOfValues){
	std::vector<glm::vec3> result;

	float step_X = (to.x - from.x) / (numberOfValues - 1);
	float step_Y = (to.y - from.y) / (numberOfValues - 1);
	float step_Z = (to.z - from.z) / (numberOfValues - 1);

	for (int i = 0; i < numberOfValues; i++) {
		float interpolated_X = from.x + i * step_X;
		float interpolated_Y = from.y + i * step_Y;
		float interpolated_Z = from.z + i * step_Z;

		glm::vec3 interpolatedVec(interpolated_X, interpolated_Y, interpolated_Z);
        result.push_back(interpolatedVec);
    }
    
    return result;
}
std::vector<TexturePoint> interpolateTexturePoints(TexturePoint from, TexturePoint to, float numberOfValues) {
    std::vector<TexturePoint> result;

    float stepX = (to.x - from.x) / (numberOfValues - 1);
    float stepY = (to.y - from.y) / (numberOfValues - 1);

    for (int i = 0; i < numberOfValues; ++i) {
        TexturePoint texturePoint;
        texturePoint.x = from.x + stepX * i;
        texturePoint.y = from.y + stepY * i;

        result.push_back(texturePoint);
    }
    return result;
}
std::vector<CanvasPoint> interpolateCanvasPoints(CanvasPoint from, CanvasPoint to, float numberOfValues) {
    std::vector<CanvasPoint> result;

    float stepX = (to.x - from.x) / (numberOfValues - 1);
    float stepY = (to.y - from.y) / (numberOfValues - 1);

    for (int i = 0; i < numberOfValues; ++i) {
        CanvasPoint canvasPoints;
        canvasPoints.x = from.x + stepX * i;
        canvasPoints.y = from.y + stepY * i;


        result.push_back(canvasPoints);
    }
    return result;
}