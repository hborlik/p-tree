#include <skinning.hpp>

namespace ptree {

void Skin(int faces, const Skeleton& skeleton, std::vector<Vertex>& vertices, std::vector<uint32_t>& indices)
{
    const std::vector<Joint> &skeleton_joints = skeleton.joints;
    const std::vector<uint32_t> &skeleton_indices = skeleton.indices;

    std::vector<glm::vec3> v_pos;
    glm::vec3 I = {0.5, 0, 0};
    glm::mat3 rot = glm::rotate(glm::mat4{1.0f}, (float)M_PI * 2.0f / faces, glm::vec3{0, 0, 1});
    for(int i = 0; i < faces; i++) {
        v_pos.push_back(I);
        I = rot * I;
    }

    const uint32_t NVerts = v_pos.size();
    for (int i = 0; i < skeleton_indices.size(); i += 2) {

        const uint32_t  vi_a = skeleton_indices[i],
                        vi_b = skeleton_indices[i + 1];
        const Joint& joint_a = skeleton_joints[vi_a];
        const Joint& joint_b = skeleton_joints[vi_b];
        const glm::vec3 sk_va_pos = glm::vec3(joint_a.tr.position);
        const glm::vec3 sk_vb_pos = glm::vec3(joint_b.tr.position);
        const float bone_length = glm::length(sk_vb_pos - sk_va_pos);

        const int VertFirstInd = vertices.size(); // initial vertex buffer size

        const glm::mat4 bone_basis_a = joint_a.transform();
        const glm::mat4 bone_basis_b = joint_b.transform();

        const glm::mat4 transform_a = bone_basis_a;
        const glm::mat4 transform_b = bone_basis_b;
        
        const glm::mat3 normal_transform_a = glm::transpose(glm::inverse(glm::mat3(transform_a)));
        const glm::mat3 normal_transform_b = glm::transpose(glm::inverse(glm::mat3(transform_b)));

        Vertex nv{};
        nv.color = glm::vec4(1.0, 1.0, 0.5, 1.0);

        for (int v = 0; v < NVerts; ++v) {
            // vertex set a
            nv.pos = transform_a * glm::vec4(v_pos[v], 1.0f);
            nv.normal = glm::normalize(normal_transform_a * v_pos[v]);
            vertices.push_back(nv);
            // vertex set b
            nv.pos = transform_b * glm::vec4(v_pos[v], 1.0f);
            nv.normal = glm::normalize(normal_transform_b * v_pos[v]);
            vertices.push_back(nv);
        }

        for (int f = 0; f < NVerts; ++f) {
            const int inds = 2 * f;
            indices.push_back(VertFirstInd + inds);
            indices.push_back(VertFirstInd + (inds+2) % (NVerts * 2));
            indices.push_back(VertFirstInd + inds+1);

            indices.push_back(VertFirstInd + (inds+2) % (NVerts * 2));
            indices.push_back(VertFirstInd + (inds+3) % (NVerts * 2));
            indices.push_back(VertFirstInd + inds+1);
        }
    }
}

void Skin_GO(int faces, const Skeleton& skeleton, std::vector<Vertex>& vertices, std::vector<uint32_t>& indices)
{
    const std::vector<Joint> &skeleton_joints = skeleton.joints;
    const std::vector<uint32_t> &skeleton_indices = skeleton.indices;

    std::vector<glm::vec3> v_pos;
    glm::vec3 I = {0.5, 0, 0};
    glm::mat3 rot = glm::rotate(glm::mat4{1.0f}, (float)M_PI * 2.0f / faces, glm::vec3{0, 0, 1});
    for(int i = 0; i < faces; i++) {
        v_pos.push_back(I);
        I = rot * I;
    }

    const uint32_t NVerts = v_pos.size();
    for (int i = 0; i < skeleton_indices.size(); i += 2) {

        const uint32_t  vi_a = skeleton_indices[i],
                        vi_b = skeleton_indices[i + 1];
        const Joint& joint_a = skeleton_joints[vi_a];
        const Joint& joint_b = skeleton_joints[vi_b];
        const glm::vec3 sk_va_pos = glm::vec3(joint_a.tr.position);
        const glm::vec3 sk_vb_pos = glm::vec3(joint_b.tr.position);

        const int VertFirstInd = vertices.size(); // initial vertex buffer size

        const glm::mat4 bone_basis_a = joint_a.transform();
        const glm::mat4 bone_basis_b = joint_b.transform();

        const float a_scale = joint_a.width_scale;
        const float b_scale = joint_b.width_scale;

        // normal vector cross forward and x axis
        const glm::vec3 a_normal = glm::normalize(glm::cross(glm::vec3{bone_basis_a[2]}, {1, 0, 0}));
        const glm::vec3 b_normal = glm::normalize(glm::cross(glm::vec3{bone_basis_b[2]}, {1, 0, 0}));

        const glm::vec3 a_binormal = glm::normalize(glm::cross(glm::vec3{bone_basis_a[2]}, a_normal));
        const glm::vec3 b_binormal = glm::normalize(glm::cross(glm::vec3{bone_basis_b[2]}, b_normal));


        const glm::mat4 transform_a = glm::mat4{
            glm::vec4{a_normal, 0.0f},
            glm::vec4{a_binormal, 0.0f}, 
            glm::normalize(bone_basis_a[2]), 
            glm::vec4{sk_va_pos, 1}
        } * glm::scale(glm::mat4{1.0f}, {a_scale, a_scale, 1.0f});
        const glm::mat4 transform_b = glm::mat4{
            glm::vec4{b_normal, 0.0f},
            glm::vec4{b_binormal, 0.0f}, 
            glm::normalize(bone_basis_b[2]),
            glm::vec4{sk_vb_pos, 1}
        } * glm::scale(glm::mat4{1.0f}, {b_scale, b_scale, 1.0f});
        
        const glm::mat3 normal_transform_a = glm::transpose(glm::inverse(glm::mat3(transform_a)));
        const glm::mat3 normal_transform_b = glm::transpose(glm::inverse(glm::mat3(transform_b)));

        Vertex nv{};
        nv.color = glm::vec4(1.0, 1.0, 0.5, 1.0);

        for (int v = 0; v < NVerts; ++v) {
            // vertex set a
            nv.pos = transform_a * glm::vec4(v_pos[v], 1.0f);
            nv.normal = glm::normalize(normal_transform_a * v_pos[v]);
            vertices.push_back(nv);
            // vertex set b
            nv.pos = transform_b * glm::vec4(v_pos[v], 1.0f);
            nv.normal = glm::normalize(normal_transform_b * v_pos[v]);
            vertices.push_back(nv);
        }

        for (int f = 0; f < NVerts; ++f) {
            const int inds = 2 * f;
            indices.push_back(VertFirstInd + inds);
            indices.push_back(VertFirstInd + (inds+2) % (NVerts * 2));
            indices.push_back(VertFirstInd + inds+1);

            indices.push_back(VertFirstInd + (inds+2) % (NVerts * 2));
            indices.push_back(VertFirstInd + (inds+3) % (NVerts * 2));
            indices.push_back(VertFirstInd + inds+1);
        }
    }
}

} // namespace ptree