#include <procedural_tree.hpp>
#include <spline.hpp>

#include <math.h>

namespace ptree {

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

constexpr Symbol<TreeSymbol> Axiom = {S_A, {2, 0.5}};

struct MonopodialProduction : public Production<TreeSymbol> {

    const float R_1   = 0.9f;             /* contraction ratio for the trunk */
    const float R_2   = 0.7f;             /* contraction ratio for the branches */
    const float a_0   = degToRad(30);     /* branching angle from the trunk */
    const float a_2   = degToRad(-30);     /* branching angle for the lateral axes */
    const float d     = degToRad(137.5f); /* divergence angle */
    const float w_r   = 0.707f;           /* width decrease rate */

    MonopodialProduction() = default;

    MonopodialProduction(float p, uint32_t sym) : Production<TreeSymbol>{p, sym} {}

    MonopodialProduction(const std::map<std::string, float> &param, float p, uint32_t sym) : Production<TreeSymbol>{p, sym},
        R_1{param.at("R_1")},
        R_2{param.at("R_2")},
        a_0{param.at("a_0")}, 
        a_2{param.at("a_2")},
        d{param.at("d")},
        w_r{param.at("d_r")}
    {
    }

    bool matches(const SymbolN<TreeSymbol>& sym) const override {
        return sym.center()->RepSym == this->A;
    }
};


struct P_1 : public MonopodialProduction {

    P_1() : MonopodialProduction{1.0f, S_A} {}

    SymbolString<TreeSymbol> translate(const SymbolN<TreeSymbol>& sym) const override {
        const TreeSymbol& value = sym.center()->value;
        float L = value.l;
        float W = value.w;
        SymbolString<TreeSymbol> ret{};
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

    SymbolString<TreeSymbol> translate(const SymbolN<TreeSymbol>& sym) const override {
        const TreeSymbol& value = sym.center()->value;
        float L = value.l;
        float W = value.w;
        SymbolString<TreeSymbol> ret{};
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

    SymbolString<TreeSymbol> translate(const SymbolN<TreeSymbol>& sym) const override {
        const TreeSymbol& value = sym.center()->value;
        float L = value.l;
        float W = value.w;
        SymbolString<TreeSymbol> ret{};
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
    LSystemTr<TreeSymbol> mlsys{};

    P_1 p1{};
    P_2 p2{};
    P_3 p3{};

    mlsys.add_rule(&p1);
    mlsys.add_rule(&p2);
    mlsys.add_rule(&p3);

    SymbolString<TreeSymbol> str{Axiom};

    for (int i = 0; i < iterations; ++i) {
        str = mlsys.evaluate(str);
    }

    SplineSkeleton ssk;
    HermiteSpline<Transform> spline{};
    spline.points.push_back(Transform{glm::quatLookAt(-glm::normalize(glm::vec3{1, 1, 1}), {0, 1, 0}), {}});
    spline.points.push_back(Transform{glm::quatLookAt(-glm::normalize(glm::vec3{1, 1, 1}), {0, 1, 0}), {1, 1, 1}});
    spline.points.push_back(Transform{glm::quatLookAt(-glm::normalize(glm::vec3{1, 1, 1}), {0, 1, 0}), {1.2, 1.2, 1.2}});
    spline.points.push_back(Transform{glm::quatLookAt(-glm::normalize(glm::vec3{1, 1, 1}), {0, 1, 0}), {2, 2, 2}});
    spline.points.push_back(Transform{glm::quatLookAt(-glm::normalize(glm::vec3{0, 0.8, 0.05}), {0, 1, 0}), {3, 5, 3}});
    ssk.sections.push_back(spline);

    // return ssk.toSkeleton(50);

    Tree tree{};
    if (tree.from_symbol_string(str)) {
        // tree.simple_skeleton(5);
        return tree.to_skeleton();
    }
    return {};
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

void Turtle::forward(float distance) {
    push_vertex();
    position += heading() * distance;
}

void Turtle::skip(float distance) {
    push_vertex();
    branch();
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

void Turtle::branch() {
    Branch b{};
    b.parent = current_branch;

    const int ind = current_branch->children.size();
    current_branch->children.push_back(b);
    current_branch = &(current_branch->children[ind]);
}

void Turtle::push_vertex() {
    Joint joint = joint_transform();
    if (current_branch != nullptr) {        
        current_branch->joints.push_back(joint);
    }
    depth++;
}

void Tree::apply_tropism(Turtle& turtle, const glm::vec3& T, float F, float b_l) {
    const glm::vec3 hxt = glm::cross(turtle.heading(), T);
    const float alpha = F * hxt.length() * b_l * b_l / (2.0f * turtle.width);
    turtle.rotation = glm::rotate(glm::quat{glm::mat4{1.0f}}, alpha, hxt) * turtle.rotation;
}

void Tree::eval_turtle_step(uint32_t sym, float value, uint32_t depth, Turtle& turtle, std::stack<Turtle>& turtle_stack) {
    switch(sym) {
        case TurtleCommands::S_forward:
        {
            if (depth > 0) {
                apply_tropism(turtle, GravityDir, 1.0f / (depth + 1.0f) * 0.02f, value);
            }
            turtle.forward(value);
            break;
        }
        case TurtleCommands::S_skip:
            turtle.skip(value);
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
            turtle.branch();
            break;
        case TurtleCommands::S_pop:
            turtle.push_vertex();
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