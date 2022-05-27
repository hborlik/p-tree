#include <skinning.hpp>

namespace ptree {

/**
 * @brief returns a vertor that is non parallel to input vector and lies in the y axis plane
 * 
 * @param v 
 * @return glm::vec3 
 */
glm::vec3 nonParallelVectorInHoriz(const glm::vec3& v) {
    if (glm::dot(v, {1, 0, 0}) > 0.999)
        return {0, 0, 1};
    return {1, 0, 0};
}

/**
 * @brief harden mesh normals, ccw front facing
 * 
 * @param vertices 
 * @param indices 
 */
void harden_normals(std::vector<Vertex>& vertices, std::vector<uint32_t>& indices) {
    for (int i = 0; i < indices.size(); i+=3) {
        Vertex &A = vertices[indices[i+0]];
        Vertex &B = vertices[indices[i+1]];
        Vertex &C = vertices[indices[i+2]];
        glm::vec3 X = C.pos - B.pos;
        glm::vec3 Y = A.pos - B.pos;

        glm::vec3 N = glm::normalize(glm::cross(X, Y));
        A.normal = N;
        B.normal = N;
        C.normal = N;
    }
}


void Skin_GO(int faces, const Skeleton& skeleton, std::vector<Vertex>& vertices, std::vector<uint32_t>& indices, bool hard_normals, float thickness, Colorizer* colorizer)
{
    vertices.clear();
    indices.clear();

    const std::vector<Joint> &skeleton_joints = skeleton.joints;
    const std::vector<uint32_t> &skeleton_indices = skeleton.indices;
    
    vertices.reserve(skeleton_indices.size() * faces * 2 * (hard_normals ? 2 : 1));
    indices.reserve(skeleton_indices.size() * 6);

    std::vector<glm::vec3> v_pos;
    glm::vec3 I = {0.5, 0, 0};
    for(int i = 0; i < faces; i++) {
        float rad = (float)M_PI * 2.0f / faces * i; // radians to rotate
        v_pos.push_back(glm::rotate(glm::mat4{1.0f}, rad, glm::vec3{0, 0, 1}) * glm::vec4(I, .0f));
    }

    const uint32_t NVerts = v_pos.size();
    for (int i = 0; i < skeleton_indices.size(); i += 2) {

        const uint32_t  vi_a = skeleton_indices[i],
                        vi_b = skeleton_indices[i + 1];
        const Joint& joint_a = skeleton_joints[vi_a];
        const Joint& joint_b = skeleton_joints[vi_b];
        const glm::vec3 sk_va_pos = glm::vec3(joint_a.position);
        const glm::vec3 sk_vb_pos = glm::vec3(joint_b.position);

        const int VertFirstInd = vertices.size(); // initial vertex buffer size
        const int IndexFirstInd = indices.size(); // initial index buffer

        const float a_scale = joint_a.width_scale * thickness;
        const float b_scale = joint_b.width_scale * thickness;

        const glm::vec3 bone_tangent_a = glm::normalize(joint_a.tangent);
        const glm::vec3 bone_tangent_b = glm::normalize(joint_b.tangent);

        const glm::vec3 horiz_a = nonParallelVectorInHoriz(bone_tangent_a);
        const glm::vec3 horiz_b = nonParallelVectorInHoriz(bone_tangent_b);

        // normal vector cross forward and x axis
        const glm::vec3 a_normal = glm::normalize(glm::cross(glm::vec3{bone_tangent_a}, horiz_a));
        glm::vec3 b_normal = glm::normalize(glm::cross(glm::vec3{bone_tangent_b}, horiz_b));

        const glm::vec3 a_binormal = glm::normalize(glm::cross(glm::vec3{bone_tangent_a}, a_normal));
        glm::vec3 b_binormal = glm::normalize(glm::cross(glm::vec3{bone_tangent_b}, b_normal));

        if (glm::dot(a_binormal, b_binormal) < 0) {
            b_normal = -b_normal;
            b_binormal = -b_binormal;
        }

        const glm::mat4 transform_a = glm::mat4{
            glm::vec4{a_normal, 0.0f},
            glm::vec4{a_binormal, 0.0f}, 
            glm::vec4{glm::normalize(bone_tangent_a), 0.0f}, 
            glm::vec4{sk_va_pos, 1}
        } * glm::scale(glm::mat4{1.0f}, {a_scale, a_scale, 1.0f});
        const glm::mat4 transform_b = glm::mat4{
            glm::vec4{b_normal, 0.0f},
            glm::vec4{b_binormal, 0.0f}, 
            glm::vec4{glm::normalize(bone_tangent_b), 0.0f},
            glm::vec4{sk_vb_pos, 1}
        } * glm::scale(glm::mat4{1.0f}, {b_scale, b_scale, 1.0f});
        
        const glm::mat3 normal_transform_a = glm::transpose(glm::inverse(glm::mat3(transform_a)));
        const glm::mat3 normal_transform_b = glm::transpose(glm::inverse(glm::mat3(transform_b)));

        Vertex nv{};
        nv.color = glm::vec4(1.0, 1.0, 0.5, 1.0);

        // push two circles of vertices to vertex buffer
        for (int v = 0; v < NVerts; ++v) {
            // vertex set a
            nv.pos = transform_a * glm::vec4(v_pos[v], 1.0f);
            nv.normal = glm::normalize(normal_transform_a * v_pos[v]);
            nv.color = colorizer ? colorizer->color_depth(joint_a.depth) : glm::vec3{};
            vertices.push_back(nv);
            if (hard_normals) {
                vertices.push_back(nv);
            }
            // vertex set b
            nv.pos = transform_b * glm::vec4(v_pos[v], 1.0f);
            nv.normal = glm::normalize(normal_transform_b * v_pos[v]);
            nv.color = colorizer ? colorizer->color_depth(joint_b.depth) : glm::vec3{};
            vertices.push_back(nv);
            if (hard_normals) {
                vertices.push_back(nv);
            }
        }

        std::size_t n_pushed_verts = vertices.size() - VertFirstInd;

        // stitch vertices into triangles
        for (int f = 0; f < NVerts; f += 1) {
            // each quad has two triangles
            const int inds = 2 * f;
            if (hard_normals) {
                indices.push_back(VertFirstInd + 2 * inds);
                indices.push_back(VertFirstInd + (2 * (inds+2) + 1) % (n_pushed_verts));
                indices.push_back(VertFirstInd + 2 * (inds+1));

                indices.push_back(VertFirstInd + (2 * (inds+2) + 1) % (n_pushed_verts));
                indices.push_back(VertFirstInd + (2 * (inds+3) + 1) % (n_pushed_verts));
                indices.push_back(VertFirstInd + 2 * (inds+1));
            } else {
                indices.push_back(VertFirstInd + inds);
                indices.push_back(VertFirstInd + (inds+2) % (n_pushed_verts));
                indices.push_back(VertFirstInd + inds+1);

                indices.push_back(VertFirstInd + (inds+2) % (n_pushed_verts));
                indices.push_back(VertFirstInd + (inds+3) % (n_pushed_verts));
                indices.push_back(VertFirstInd + inds+1);
            }
        }
    }

    if (hard_normals)
        harden_normals(vertices, indices);
}

} // namespace ptree