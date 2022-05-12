/**
 * @file skinning.hpp
 * @brief tree skinning
 * @date 2022-04-24
 * 
 */
#ifndef PTREE_SKINNING_H
#define PTREE_SKINNING_H

#include <procedural_tree.hpp>

namespace ptree {

class Colorizer {
public:
    virtual ~Colorizer() = default;
    virtual glm::vec3 color_depth(int d) = 0;
};

class DefaultColorizer : public Colorizer {
public:
    DefaultColorizer(glm::vec3 c0, glm::vec3 c1, int max_depth) : c0{c0}, c1{c1}, max_depth{max_depth} {}

    glm::vec3 color_depth(int d) override {
        float m = (float)d / max_depth;
        return (1 - m) * c0 + m * c1;
    }

    int max_depth = 0;
    glm::vec3 c0, c1;
};

// skin using global orientation vectors and ignore joint transforms
void Skin_GO(int faces, const Skeleton& skeleton, std::vector<Vertex>& vertices, std::vector<uint32_t>& indices, bool hard_normals, float thickness, Colorizer* colorizer);

}

#endif // PTREE_SKINNING_H