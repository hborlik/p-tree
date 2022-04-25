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

// skin using global orientation vectors and ignore joint transforms
void Skin_GO(int faces, const Skeleton& skeleton, std::vector<Vertex>& vertices, std::vector<uint32_t>& indices);

}

#endif // PTREE_SKINNING_H