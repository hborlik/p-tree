#include <procedural_tree.hpp>

#include <math.h>

namespace ptree {

constexpr float degToRad(const float deg) noexcept {
    return deg / 180.f * M_PI;
};

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
            ret.push_back({S_forward, 1});
            ret.push_back({S_push});
            ret.push_back({S_yaw, M_PI_2 * 0.9});
            ret.push_back({S_A});
            ret.push_back({S_pop});
            ret.push_back({S_push});
            ret.push_back({S_yaw, -M_PI_2 * 0.9});
            ret.push_back({S_A});
            ret.push_back({S_pop});
        }
        return ret;
    }
};

template<typename T>
struct P_2 : public Production<T> {
    const float _R;

    explicit P_2(float R) : Production<T>{1.0f, S_forward}, _R{R} {}

    bool matches(const SymbolN<T>& sym) const override {
        return sym.center()->RepSym;
    } 

    SymbolString<T> translate(const SymbolN<T>& sym) const override {
        SymbolString<T> ret{};
        if (matches(sym)) {
            // F(s * R)
            ret.push_back({S_forward, sym.center()->value * _R});
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

constexpr uint32_t S_A      = ptree::S_A;
constexpr uint32_t S_B      = S_A + 1;
constexpr uint32_t S_C      = S_B + 1;

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
constexpr const float R_1   = 0.9f;             /* contraction ratio for the trunk */
constexpr const float R_2   = 0.7f;             /* contraction ratio for the branches */
constexpr const float a_0   = degToRad(30);     /* branching angle from the trunk */
constexpr const float a_2   = degToRad(-30);     /* branching angle for the lateral axes */
constexpr const float d     = degToRad(137.5f); /* divergence angle */
constexpr const float w_r   = 0.707f;           /* width decrease rate */

constexpr Symbol<Pair> Axiom = {S_A, {1, 10}};


struct P_1 : public Production<Pair> {

    P_1() : Production<Pair>{1.0f, S_A} {}

    bool matches(const SymbolN<Pair>& sym) const override {
        return sym.center()->RepSym == this->A;
    }

    SymbolString<Pair> translate(const SymbolN<Pair>& sym) const override {
        const Pair& value = sym.center()->value;
        float L = value.l;
        float W = value.w;
        SymbolString<Pair> ret{};
        if (matches(sym)) {
            // !(w) F(l) [ &(a0) B(l * R_2, w * w_r) ] /(d) A(l * R_1, w * w_r)
            ret.push_back({S_Bang, {W}});
            ret.push_back({S_forward, {L}});
            ret.push_back({S_push});
            ret.push_back({S_pitch, {a_0}});
            ret.push_back({S_B, {L * R_2, W * w_r}});
            ret.push_back({S_pop});
            ret.push_back({S_roll, {d}});
            ret.push_back({S_A, {L * R_1, W * w_r}});
        }
        return ret;
    }
};

struct P_2 : public Production<Pair> {
    P_2() : Production<Pair>{1.0f, S_B} {}

    bool matches(const SymbolN<Pair>& sym) const override {
        return sym.center()->RepSym == this->A;
    } 

    SymbolString<Pair> translate(const SymbolN<Pair>& sym) const override {
        const Pair& value = sym.center()->value;
        float L = value.l;
        float W = value.w;
        SymbolString<Pair> ret{};
        if (matches(sym)) {
            // !(w) F(L) [ -(a_2) $ C(l * R_2, w * w_r) ] C(l * R_1, w * w_r)
            ret.push_back({S_Bang, {W}});
            ret.push_back({S_forward, {L}});
            ret.push_back({S_push});
            ret.push_back({S_yaw, {-a_2}});
            ret.push_back({S_Dollar});
            ret.push_back({S_C, {L * R_2, W * w_r}});
            ret.push_back({S_pop});
            ret.push_back({S_C, {L * R_1, W * w_r}});
        }
        return ret;
    }
};

struct P_3 : public Production<Pair> {
    P_3() : Production<Pair>{1.0f, S_C} {}

    bool matches(const SymbolN<Pair>& sym) const override {
        return sym.center()->RepSym == this->A;
    } 

    SymbolString<Pair> translate(const SymbolN<Pair>& sym) const override {
        const Pair& value = sym.center()->value;
        float L = value.l;
        float W = value.w;
        SymbolString<Pair> ret{};
        if (matches(sym)) {
            // !(w) F(L) [ +(a_2) $ B(l * R_2, w * w_r) ] B(l * R_1, w * w_r)
            ret.push_back({S_Bang, {W}});
            ret.push_back({S_forward, {L}});
            ret.push_back({S_push});
            ret.push_back({S_yaw, {a_2}});
            ret.push_back({S_Dollar});
            ret.push_back({S_B, {L * R_2, W * w_r}});
            ret.push_back({S_pop});
            ret.push_back({S_B, {L * R_1, W * w_r}});
        }
        return ret;
    }
};

}

void CreateSkeleton(int iterations, std::vector<glm::mat4>& joints, std::vector<uint32_t>& indices) {
    using namespace monopodial;
    LSystemTr<Pair> lsys{};

    P_1 p1{};
    P_2 p2{};
    P_3 p3{};

    lsys.add_rule(&p1);
    lsys.add_rule(&p2);
    lsys.add_rule(&p3);

    SymbolString<Pair> str{Axiom};

    // PrintSymbolString<float>(ProceduralTree::Library, str);

    for (int i = 0; i < iterations; ++i) {
        str = lsys.evaluate(str);
        // PrintSymbolString<float>(ProceduralTree::Library, str);
    }

    ProceduralTree eval{};

    eval.to_skeleton(str, joints, indices);
}

void Skin(const std::vector<glm::mat4> &skeleton_joints, const std::vector<uint32_t> &skeleton_indices,
          std::vector<Vertex>& vertices, std::vector<uint32_t>& indices)
{
    for (int i = 0; i < skeleton_indices.size(); i += 2) {

        const uint32_t  vi_a = skeleton_indices[i], 
                        vi_b = skeleton_indices[i + 1];
        const glm::vec3 sk_va_pos = glm::vec3(skeleton_joints[vi_a][3]);
        const glm::vec3 sk_vb_pos = glm::vec3(skeleton_joints[vi_b][3]);
        const glm::mat4& bone_basis = skeleton_joints[vi_a];
        const float bone_length = glm::length(sk_vb_pos - sk_va_pos);

        const int l_v_pos = vertices.size(); // initial vertex buffer size

        glm::mat4 transform = glm::scale(glm::translate(bone_basis, {0, 0, bone_length / 2}), {0.01, 0.01, bone_length});

        Vertex nv{};
        nv.color = glm::vec4(1.0, 1.0, 0.5, 1.0);

        glm::vec3 v_pos[] = {
            {0.5, 0.5, 0.5},
            {0.5, -0.5, -0.5},
            {0.5, -0.5, 0.5},
            {0.5, 0.5, -0.5},

            {-0.5, 0.5, 0.5},
            {-0.5, -0.5, -0.5},
            {-0.5, -0.5, 0.5},
            {-0.5, 0.5, -0.5}
        };

        for (int v = 0; v < sizeof(v_pos) / sizeof(glm::vec3); ++v) {
            // nv.pos = sk_v_pos + v_pos[v];
            nv.pos = transform * glm::vec4(v_pos[v], 1.0f);
            nv.normal = glm::normalize(skeleton_joints[i] * glm::vec4(v_pos[v], 0.f));
            vertices.push_back(nv);
        }


        indices.push_back(0 + l_v_pos);
        indices.push_back(2 + l_v_pos);
        indices.push_back(1 + l_v_pos);
        indices.push_back(0 + l_v_pos);
        indices.push_back(1 + l_v_pos);
        indices.push_back(3 + l_v_pos);

        indices.push_back(4 + l_v_pos);
        indices.push_back(5 + l_v_pos);
        indices.push_back(6 + l_v_pos);
        indices.push_back(4 + l_v_pos);
        indices.push_back(7 + l_v_pos);
        indices.push_back(5 + l_v_pos);

        indices.push_back(4 + l_v_pos);
        indices.push_back(0 + l_v_pos);
        indices.push_back(3 + l_v_pos);
        indices.push_back(4 + l_v_pos);
        indices.push_back(3 + l_v_pos);
        indices.push_back(7 + l_v_pos);

        indices.push_back(6 + l_v_pos);
        indices.push_back(2 + l_v_pos);
        indices.push_back(0 + l_v_pos);
        indices.push_back(6 + l_v_pos);
        indices.push_back(0 + l_v_pos);
        indices.push_back(4 + l_v_pos);

        indices.push_back(3 + l_v_pos);
        indices.push_back(1 + l_v_pos);
        indices.push_back(5 + l_v_pos);
        indices.push_back(3 + l_v_pos);
        indices.push_back(5 + l_v_pos);
        indices.push_back(7 + l_v_pos);

        indices.push_back(5 + l_v_pos);
        indices.push_back(2 + l_v_pos);
        indices.push_back(6 + l_v_pos);
        indices.push_back(5 + l_v_pos);
        indices.push_back(1 + l_v_pos);
        indices.push_back(2 + l_v_pos);
    }
}

const char* ProceduralTree::Library[] = {
    "Nil",
    "F",
    "f", 
    "+", 
    "&", 
    "\\", 
    "[", 
    "]", 
    "A" 
};

Turtle::Turtle() : rotation{glm::identity<glm::quat>()}, position{} {

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

void Turtle::forward(float distance) {
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

void ProceduralTree::apply_tropism(Turtle& turtle, const glm::vec3& T, float F, float b_l) {
    const glm::vec3 hxt = glm::cross(turtle.heading(), T);
    const float alpha = F * hxt.length() * b_l * b_l / (2.0f * turtle.width);
    turtle.rotation = glm::rotate(glm::quat{glm::mat4{1.0f}}, alpha, hxt) * turtle.rotation;
}

} // namespace ptree