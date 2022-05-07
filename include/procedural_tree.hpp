#ifndef OL_EVAL_H
#define OL_EVAL_H

#include <stack>
#include <map>
#include <optional>

#define GLM_FORCE_RADIANS
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>

#include <ol_system.hpp>
#include <spline.hpp>

namespace ptree {

struct TreeSymbol {
    float l = 0,w = 0;

    operator float() const noexcept {return l;}
};

struct Vertex {
    glm::vec3 pos;
    glm::vec3 normal;
    glm::vec2 uv;
    glm::vec4 color;
};

struct Transform {
    glm::quat rotation;
    glm::vec3 position;

    Transform(const glm::quat& rotation, const glm::vec3& position) : rotation{rotation}, position{position} {}

    glm::mat4 transform() const noexcept {
        glm::mat4 t = glm::mat4{rotation};
        t[3] = glm::vec4{position, 1.0f};
        return t;
    }

private:
    friend Transform operator*(float scalar, const Transform& t) noexcept {
        return t * scalar;
    }

    friend Transform operator*(const Transform& t, float scalar) noexcept {
        return Transform{t.rotation * scalar, t.position * scalar};
    }

    friend Transform operator-(const Transform& t) noexcept {
        return t * -1.0f;
    }
    friend Transform operator+(const Transform& ta, const Transform& tb) noexcept {
        return Transform{ta.rotation + tb.rotation, ta.position + tb.position};
    }
    friend Transform operator-(const Transform& ta, const Transform& tb) noexcept {
        return ta + -tb;
    }
};

struct Joint {
    glm::vec3 tangent;
    glm::vec3 position;
    float width_scale;

    Joint(const glm::vec3& tangent, const glm::vec3& position, float width_scale) : tangent{tangent}, position{position}, width_scale{width_scale} {}
};

struct Branch {
    std::vector<Joint> joints;
    std::vector<Branch> children;
    Branch *parent = nullptr;

    /**
     * @brief The last joint added to vertices
     * 
     * @return const Joint& 
     */
    const Joint& current_joint() const {
        return joints[joints.size() - 1];
    }

    void add_joint(Joint joint) {
        joints.emplace_back(std::move(joint));
    }
};

struct Skeleton {
    std::vector<Joint> joints;
    std::vector<uint32_t> indices;
};

/**
 * @brief skeleton made of spline sections
 * 
 */
struct SplineSkeleton {
    std::vector<HermiteSpline<Transform>> sections;

    Skeleton toSkeleton(int curve_samples) const {
        Skeleton interm_skel;
        for (const auto& section : sections) {
            int last_joint_ind = -1;
            // sample spline on section and add joints for it
            for (int i = 0; i < curve_samples; ++i) {
                const float x = float(i) / (curve_samples - 1);
                if (last_joint_ind != -1) {
                    interm_skel.indices.push_back(last_joint_ind);
                    interm_skel.indices.push_back(interm_skel.joints.size());
                }
                last_joint_ind = interm_skel.joints.size();
                interm_skel.joints.emplace_back<Joint>({section.eval_dt(x).position, section.eval(x).position, 1.0f});
            }
        }
        return interm_skel;
    }
};

constexpr glm::vec3 GravityDir  = {0, -1, 0};
constexpr glm::vec3 UpDir       = -GravityDir;

using FSymbol = Symbol<float>;

Skeleton CreateSkeleton(int iterations);

struct Turtle {
    glm::quat rotation;
    glm::vec3 position;
    float width = 0.1f;

    Branch *current_branch = nullptr; // current branch, used for forward commands

    Turtle();

    glm::vec3 heading() const noexcept  {return glm::mat3(rotation)[2];}
    glm::vec3 up() const noexcept       {return glm::mat3(rotation)[1];}
    glm::vec3 left() const noexcept     {return glm::mat3(rotation)[0];}

    void yaw(float rad);
    void roll(float rad);
    void pitch(float rad);

    void forward(float distance);
    void skip(float distance);

    void level();

    Joint joint_transform() const noexcept {
        return {
            heading(),
            position,
            width
        };
    }

    /**
     * @brief split this turtle from parent branch by creating a new branch and adding it to previous as a child
     * 
     */
    void branch() {
        Branch b{};
        b.parent = current_branch;

        const int ind = current_branch->children.size();
        current_branch->children.push_back(b);
        current_branch = &(current_branch->children[ind]);
    }

    // push edge by adding position to the branch
    void push_vertex();
};

static const char* Library[] = {
    "Nil",
    "F",
    "f", 
    "+", 
    "&", 
    "\\",
    "[", 
    "]", 
    "A",
    "B",
    "C"
};

struct TurtleCommands {
    enum CMD {
        S_forward = 1,
        S_skip    ,     
        S_yaw     ,     // +-
        S_pitch   ,     // &
        S_roll    ,     // /
        S_Dollar  ,     // special turtle command
        S_push    ,     
        S_pop     ,     
        SetWidth    ,     
        S_A             // other intermediate symbols should be above this value
    };
};


class Tree {
public:
    Branch root;

    void set_parameter(const std::string& name, float value) {parameters[name] = value;}
    float get_parameter(const std::string& name) const {return parameters.at(name);}

    /**
     * @brief creates tree branch structure from symbol string
     * 
     * @param ss symbol string created from an l-system evaluation
     * @param joints 
     * @param indices 
     */
    template<typename T>
    bool from_symbol_string(const SymbolString<T>& ss) {
        std::stack<Turtle> turtle_stack;
        Turtle turtle{};
        root = Branch{}; // reset tree
        turtle.current_branch = &root;

        for (const auto& symbol : ss) {
            const uint32_t depth = turtle_stack.size();
            // operator determination, symbol to turtle command
            eval_turtle_step(symbol.RepSym, symbol.value, depth, turtle, turtle_stack);
        }
        return turtle_stack.size() == 0;
    }

    Skeleton to_skeleton() const {
        Skeleton sk;
        std::stack<const Branch*> branch_stack;
        branch_stack.push(&root);

        while(!branch_stack.empty()) {
            const Branch *cur = branch_stack.top();
            branch_stack.pop();
            for (auto &c : cur->children) {
                branch_stack.push(&c);
            }
            int last_joint = -1;
            for (auto &j : cur->joints) {
                if (last_joint != -1) {
                    // start adding edges once there are at least two new vertices
                    sk.indices.push_back(last_joint);
                    sk.indices.push_back(sk.joints.size());
                }
                // index of this new joint
                last_joint = sk.joints.size();
                sk.joints.push_back(j);
            }
        }
        return sk;
    }

    void simple_skeleton(int len) {
        root = Branch{}; // reset tree
        Turtle turtle{};
        turtle.current_branch = &root;
        turtle.width = 2.0f;
        for (int i = 0; i < len; i++) {
            // apply_tropism(turtle, GravityDir, 4.22f, 1.0f);
            turtle.forward(0.5f);
            // turtle.yaw(0.18f);
        }
        Turtle turtle_a = turtle;
        turtle.branch();
        turtle.yaw(0.7f);
        for (int i = 0; i < len; i++) {
            // apply_tropism(turtle, GravityDir, 4.22f, 1.0f);
            turtle.forward(0.5f);
            // turtle.yaw(-0.09f);
        }
        turtle_a.forward(1.5f);
    }

private:
    std::map<std::string, float> parameters;

    void apply_tropism(Turtle& turtle, const glm::vec3& T, float F, float b_l);
    void eval_turtle_step(uint32_t sym, float value, uint32_t depth, Turtle& turtle, std::stack<Turtle>& turtle_stack);
};


} // namespace ptree

#endif // OL_EVAL_H