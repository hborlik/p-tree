#include <procedural_tree.hpp>
#include <spline.hpp>

#include <math.h>

namespace ptree {

constexpr float degToRad(const float deg) noexcept {
    return deg / 180.f * M_PI;
};

constexpr uint32_t S_A      = ptree::TurtleCommands::S_A;
constexpr uint32_t S_B      = S_A + 1;
constexpr uint32_t S_C      = S_B + 1;

namespace sympodial {

/*
 * Productions based on page 50 of Algorithmic Botany Book
 * 
 */

template<typename T>
struct P_1 : public Production<T> {

    /**
     * @brief Construct a new p 1 object
     * 
     * @param a matching character
     */
    P_1() : Production<T>{1.0f, S_A} {}

    bool matches(const SymbolN<T>& sym) const override {
        return sym.center()->RepSym == this->A;
    }

    SymbolString<T> translate(const SymbolN<T>& sym) const override {
        SymbolString<T> ret{};
        if (matches(sym)) {
            // F(1)[+A][-A]
            ret.push_back({TurtleCommands::S_forward, 1});
            ret.push_back({TurtleCommands::S_push});
            ret.push_back({TurtleCommands::S_yaw, M_PI_2 * 0.9});
            ret.push_back({S_A});
            ret.push_back({TurtleCommands::S_pop});
            ret.push_back({TurtleCommands::S_push});
            ret.push_back({TurtleCommands::S_yaw, -M_PI_2 * 0.9});
            ret.push_back({S_A});
            ret.push_back({TurtleCommands::S_pop});
        }
        return ret;
    }
};

template<typename T>
struct P_2 : public Production<T> {
    const float _R;

    explicit P_2(float R) : Production<T>{1.0f, TurtleCommands::S_forward}, _R{R} {}

    bool matches(const SymbolN<T>& sym) const override {
        return sym.center()->RepSym;
    } 

    SymbolString<T> translate(const SymbolN<T>& sym) const override {
        SymbolString<T> ret{};
        if (matches(sym)) {
            // F(s * R)
            ret.push_back({TurtleCommands::S_forward, sym.center()->value * _R});
        }
        return ret;
    }
};

}

/**
 * @brief Monopodial tree-like structures of Honda
 *  based on pg. 56 of Algorithmic Botany
 * 
 */
namespace monopodial {


struct Pair {
    float l = 0,w = 0;

    operator float() const noexcept {return l;}
};

// example a
// constexpr const float R_1   = 0.9f;             /* contraction ratio for the trunk */
// constexpr const float R_2   = 0.6f;             /* contraction ratio for the branches */
// constexpr const float a_0   = degToRad(45);     /* branching angle from the trunk */
// constexpr const float a_2   = degToRad(45);     /* branching angle for the lateral axes */
// constexpr const float d     = degToRad(137.5f); /* divergence angle */
// constexpr const float w_r   = 0.707f;           /* width decrease rate */

// // example b
// constexpr const float R_1   = 0.9f;             /* contraction ratio for the trunk */
// constexpr const float R_2   = 0.9f;             /* contraction ratio for the branches */
// constexpr const float a_0   = degToRad(45);     /* branching angle from the trunk */
// constexpr const float a_2   = degToRad(45);     /* branching angle for the lateral axes */
// constexpr const float d     = degToRad(137.5f); /* divergence angle */
// constexpr const float w_r   = 0.707f;           /* width decrease rate */

// // example c
// constexpr const float R_1   = 0.9f;             /* contraction ratio for the trunk */
// constexpr const float R_2   = 0.8f;             /* contraction ratio for the branches */
// constexpr const float a_0   = degToRad(45);     /* branching angle from the trunk */
// constexpr const float a_2   = degToRad(45);     /* branching angle for the lateral axes */
// constexpr const float d     = degToRad(137.5f); /* divergence angle */
// constexpr const float w_r   = 0.707f;           /* width decrease rate */

// example c
// constexpr const float R_1   = 0.9f;             /* contraction ratio for the trunk */
// constexpr const float R_2   = 0.7f;             /* contraction ratio for the branches */
// constexpr const float a_0   = degToRad(30);     /* branching angle from the trunk */
// constexpr const float a_2   = degToRad(-30);     /* branching angle for the lateral axes */
// constexpr const float d     = degToRad(137.5f); /* divergence angle */
// constexpr const float w_r   = 0.707f;           /* width decrease rate */

constexpr Symbol<Pair> Axiom = {S_A, {2, 0.5}};

struct MonopodialProduction : public Production<Pair> {

    const float R_1   = 0.9f;             /* contraction ratio for the trunk */
    const float R_2   = 0.7f;             /* contraction ratio for the branches */
    const float a_0   = degToRad(30);     /* branching angle from the trunk */
    const float a_2   = degToRad(-30);     /* branching angle for the lateral axes */
    const float d     = degToRad(137.5f); /* divergence angle */
    const float w_r   = 0.707f;           /* width decrease rate */

    MonopodialProduction() = default;

    MonopodialProduction(float p, uint32_t sym) : Production<Pair>{p, sym} {}

    MonopodialProduction(const std::map<std::string, float> &param, float p, uint32_t sym) : Production<Pair>{p, sym},
        R_1{param.at("R_1")},
        R_2{param.at("R_2")},
        a_0{param.at("a_0")}, 
        a_2{param.at("a_2")},
        d{param.at("d")},
        w_r{param.at("d_r")}
    {
    }

    bool matches(const SymbolN<Pair>& sym) const override {
        return sym.center()->RepSym == this->A;
    }
};


struct P_1 : public MonopodialProduction {

    P_1() : MonopodialProduction{1.0f, S_A} {}

    SymbolString<Pair> translate(const SymbolN<Pair>& sym) const override {
        const Pair& value = sym.center()->value;
        float L = value.l;
        float W = value.w;
        SymbolString<Pair> ret{};
        if (matches(sym)) {
            // !(w) F(l) [ &(a0) B(l * R_2, w * w_r) ] /(d) A(l * R_1, w * w_r)
            ret.push_back({TurtleCommands::SetWidth, {W}});
            ret.push_back({TurtleCommands::S_forward, {L}});
            ret.push_back({TurtleCommands::S_push});
            ret.push_back({TurtleCommands::S_pitch, {a_0}});
            ret.push_back({S_B, {L * R_2, W * w_r}});
            ret.push_back({TurtleCommands::S_pop});
            ret.push_back({TurtleCommands::S_roll, {d}});
            ret.push_back({S_A, {L * R_1, W * w_r}});
        }
        return ret;
    }
};

struct P_2 : public MonopodialProduction {
    P_2() : MonopodialProduction{1.0f, S_B} {}

    SymbolString<Pair> translate(const SymbolN<Pair>& sym) const override {
        const Pair& value = sym.center()->value;
        float L = value.l;
        float W = value.w;
        SymbolString<Pair> ret{};
        if (matches(sym)) {
            // !(w) F(L) [ -(a_2) $ C(l * R_2, w * w_r) ] C(l * R_1, w * w_r)
            ret.push_back({TurtleCommands::SetWidth, {W}});
            ret.push_back({TurtleCommands::S_forward, {L}});
            ret.push_back({TurtleCommands::S_push});
            ret.push_back({TurtleCommands::S_yaw, {-a_2}});
            ret.push_back({TurtleCommands::S_Dollar});
            ret.push_back({S_C, {L * R_2, W * w_r}});
            ret.push_back({TurtleCommands::S_pop});
            ret.push_back({S_C, {L * R_1, W * w_r}});
        }
        return ret;
    }
};

struct P_3 : public MonopodialProduction {
    P_3() : MonopodialProduction{1.0f, S_C} {}

    SymbolString<Pair> translate(const SymbolN<Pair>& sym) const override {
        const Pair& value = sym.center()->value;
        float L = value.l;
        float W = value.w;
        SymbolString<Pair> ret{};
        if (matches(sym)) {
            // !(w) F(L) [ +(a_2) $ B(l * R_2, w * w_r) ] B(l * R_1, w * w_r)
            ret.push_back({TurtleCommands::SetWidth, {W}});
            ret.push_back({TurtleCommands::S_forward, {L}});
            ret.push_back({TurtleCommands::S_push});
            ret.push_back({TurtleCommands::S_yaw, {a_2}});
            ret.push_back({TurtleCommands::S_Dollar});
            ret.push_back({S_B, {L * R_2, W * w_r}});
            ret.push_back({TurtleCommands::S_pop});
            ret.push_back({S_B, {L * R_1, W * w_r}});
        }
        return ret;
    }
};

}

Skeleton CreateSkeleton(int iterations) {
    using namespace monopodial;
    LSystemTr<Pair> lsys{};

    P_1 p1{};
    P_2 p2{};
    P_3 p3{};

    lsys.add_rule(&p1);
    lsys.add_rule(&p2);
    lsys.add_rule(&p3);

    SymbolString<Pair> str{Axiom};

    // PrintSymbolString<float>(Tree::Library, str);

    for (int i = 0; i < iterations; ++i) {
        str = lsys.evaluate(str);
        // PrintSymbolString<float>(Tree::Library, str);
    }

    Tree tree{};

    auto sk = tree.str_to_skeleton(str);
    if (sk) {
        // tree.simple_skeleton(10);
        return *sk;
    }
    return {};
}

void Skin(int faces, Skeleton& skeleton, std::vector<Vertex>& vertices, std::vector<uint32_t>& indices)
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
        const glm::vec3 sk_va_pos = glm::vec3(joint_a.position);
        const glm::vec3 sk_vb_pos = glm::vec3(joint_b.position);
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

void Skin_GO(int faces, Skeleton& skeleton, std::vector<Vertex>& vertices, std::vector<uint32_t>& indices)
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
        const glm::vec3 sk_va_pos = glm::vec3(joint_a.position);
        const glm::vec3 sk_vb_pos = glm::vec3(joint_b.position);

        const int VertFirstInd = vertices.size(); // initial vertex buffer size

        const glm::mat4 bone_basis_a = joint_a.transform();
        const glm::mat4 bone_basis_b = joint_b.transform();

        const float a_scale = joint_a.scale;
        const float b_scale = joint_b.scale;

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

Turtle::Turtle() : rotation{glm::quatLookAt(GravityDir, {1, 0, 0})}, position{} {

}

void Turtle::yaw(float rad) {
    rotation = glm::rotate(rotation, rad, {0, 1, 0});
}

void Turtle::roll(float rad) {
    rotation = glm::rotate(rotation, rad, {0, 0, 1});
}

void Turtle::pitch(float rad) {
    rotation = glm::rotate(rotation, rad, {1, 0, 0});
}

void Turtle::forward(float distance, std::vector<Joint>& joints, std::vector<uint32_t>& indices) {
    push_edge(joints, indices);
    position += heading() * distance;
}

void Turtle::skip(float distance, std::vector<Joint>& joints, std::vector<uint32_t>& indices) {
    push_edge(joints, indices);
    reset_line();
    position += heading() * distance;
}

void Turtle::level() {
    // orients the turtle level to the ground
    glm::vec3 h = heading();
    glm::vec3 l = left();
    glm::vec3 u = up();

    l = glm::cross(-GravityDir, h);
    l = l * (1.f / glm::length(l));

    u = glm::cross(h, l);

    rotation = glm::mat3(l, u, h);
}

void Turtle::push_edge(std::vector<Joint>& joints, std::vector<uint32_t>& indices) {
    Joint tr = joint_transform();
    if (joint_index != -1) {
        // add previous joint
        indices.push_back(joint_index);
        indices.push_back(joints.size());

        // transform rotation is blend between last and current headings
        tr.rotation = glm::slerp(joints[joint_index].rotation, tr.rotation, 0.5f);
    }
    // add current position
    joint_index = joints.size();
    joints.push_back(tr);
}

void Tree::apply_tropism(Turtle& turtle, const glm::vec3& T, float F, float b_l) {
    const glm::vec3 hxt = glm::cross(turtle.heading(), T);
    const float alpha = F * hxt.length() * b_l * b_l / (2.0f * turtle.width);
    turtle.rotation = glm::rotate(glm::quat{glm::mat4{1.0f}}, alpha, hxt) * turtle.rotation;
}

void Tree::eval_turtle_step(uint32_t sym, float value, uint32_t depth, Turtle& turtle, std::stack<Turtle>& turtle_stack, Skeleton& sk) {
    switch(sym) {
        case TurtleCommands::S_forward:
        {
            if (depth > 0) {
                apply_tropism(turtle, GravityDir, 1.0f / (depth + 1.0f) * 0.02f, value);
            }
            turtle.forward(value, sk.joints, sk.indices);
            break;
        }
        case TurtleCommands::S_skip:
            turtle.skip(value, sk.joints, sk.indices);
            break;
        case TurtleCommands::S_yaw:
            turtle.yaw(value);
            break;
        case TurtleCommands::S_pitch:
            turtle.pitch(value);
            break;
        case TurtleCommands::S_roll:
            turtle.roll(value);
            break;
        case TurtleCommands::S_push:
            turtle_stack.push(turtle);
            turtle.reset_line();
            break;
        case TurtleCommands::S_pop:
            turtle.push_edge(sk.joints, sk.indices);
            turtle = turtle_stack.top();
            turtle_stack.pop();
            break;
        case TurtleCommands::S_Dollar:
            turtle.level();
            break;
        case TurtleCommands::SetWidth:
            turtle.width = value;
            break;
        default:
            break; // no op for undefined symbols
    }
}

} // namespace ptree