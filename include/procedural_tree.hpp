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
    Transform tr;
    float width_scale;

    Joint(const Transform& tr, float width_scale) : tr{tr}, width_scale{width_scale} {}

    glm::mat4 transform() const noexcept {
        auto t = tr.transform() * glm::scale(glm::mat4{1.0f}, {width_scale, width_scale, 1.0f});
        return t;
    }
};

struct Skeleton {
    std::vector<ptree::Joint> joints;
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
                interm_skel.joints.emplace_back<Joint>({section.eval(x), 1.0f});
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

    int joint_index = -1; // joint associated with current turtle position, used for forward ops

    Turtle();

    glm::vec3 heading() const noexcept  {return glm::mat3(rotation)[2];}
    glm::vec3 up() const noexcept       {return glm::mat3(rotation)[1];}
    glm::vec3 left() const noexcept     {return glm::mat3(rotation)[0];}

    void yaw(float rad);
    void roll(float rad);
    void pitch(float rad);

    void forward(float distance, std::vector<Joint>& joints, std::vector<uint32_t>& indices);
    void skip(float distance, std::vector<Joint>& joints, std::vector<uint32_t>& indices);

    void level();

    Joint joint_transform() const noexcept {
        return {
            Transform{
                rotation,
                position},
            width
        };
    }

    void reset_line() {
        joint_index = -1;
    }

    // push edge by adding the last and current position to the skeleton
    void push_edge(std::vector<Joint>& joints, std::vector<uint32_t>& indices);
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


struct Branch {
    std::vector<glm::vec3> vertices;
    std::vector<Branch> children;
    Branch *parent = nullptr;

};

class Tree {
public:
    Branch root;

    void set_parameter(const std::string& name, float value) {parameters[name] = value;}
    float get_parameter(const std::string& name) const {return parameters.at(name);}

    /**
     * @brief converts a given symbol string to a 3d representation. Line primitives are places in index array
     * 
     * @param ss 
     * @param joints 
     * @param indices 
     */
    template<typename T>
    std::optional<Skeleton> from_symbol_string(const SymbolString<T>& ss) {
        Skeleton sk;
        std::stack<Turtle> turtle_stack;
        Turtle turtle{};

        for (const auto& symbol : ss) {
            const uint32_t depth = turtle_stack.size();
            // operator determination, symbol to turtle command
            eval_turtle_step(symbol.RepSym, symbol.value, depth, turtle, turtle_stack, sk);
        }
        if (turtle_stack.size() == 0)
            return {sk};
        return {};
    }

    Skeleton simple_skeleton(int len) {
        Skeleton skeleton;
        Turtle turtle{};
        turtle.width = 2.0f;
        for (int i = 0; i < len; i++) {
            // apply_tropism(turtle, GravityDir, 4.22f, 1.0f);
            turtle.forward(0.5f, skeleton.joints, skeleton.indices);
            turtle.yaw(0.2f);
        }
        return skeleton;
    }

private:
    std::map<std::string, float> parameters;

    void apply_tropism(Turtle& turtle, const glm::vec3& T, float F, float b_l);
    void eval_turtle_step(uint32_t sym, float value, uint32_t depth, Turtle& turtle, std::stack<Turtle>& turtle_stack, Skeleton& sk);
};


} // namespace ptree

#endif // OL_EVAL_H