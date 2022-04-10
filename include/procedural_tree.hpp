#ifndef OL_EVAL_H
#define OL_EVAL_H

#include <stack>
#include <map>

#define GLM_FORCE_RADIANS
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>

#include <ol_system.hpp>

namespace ptree {

struct Vertex {
    glm::vec3 pos;
    glm::vec3 normal;
    glm::vec2 uv;
    glm::vec4 color;
};

constexpr glm::vec3 GravityDir  = {0, -1, 0};
constexpr glm::vec3 UpDir       = -GravityDir;

using FSymbol = Symbol<float>;

void CreateSkeleton(int iterations, std::vector<glm::mat4>& joints, std::vector<uint32_t>& indices);

void Skin(const std::vector<glm::mat4>& skeleton_joints, const std::vector<uint32_t>& skeleton_indices, 
        std::vector<Vertex>& vertices, std::vector<uint32_t>& indices);

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

    void forward(float distance, std::vector<glm::mat4>& joints, std::vector<uint32_t>& indices);
    void skip(float distance);

    void level();

    glm::mat4 transform() const {
        glm::mat4 tr = glm::mat4{rotation} * scaling();
        tr[3] = glm::vec4(position, 1.0f);
        return tr;
    }

    glm::mat4 linear_transform() const {
        glm::mat4 tr = glm::mat4{rotation};
        tr[3] = glm::vec4(position, 1.0f);
        return tr;
    }

    glm::mat4 scaling() const noexcept {
        return glm::scale(glm::mat4{1.0f}, {width, width, 1.0f});
    }

    void reset_line() {
        joint_index = -1;
    }

    // push edge by adding the last and current position to the skeleton
    void push_edge(std::vector<glm::mat4>& joints, std::vector<uint32_t>& indices);
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

class ProceduralTree {
public:
    std::vector<glm::mat4> joints;
    std::vector<uint32_t> indices;

    static const char* Library[];

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
    bool str_to_skeleton(const SymbolString<T>& ss) {
        std::stack<Turtle> turtle_stack;
        Turtle turtle{};

        for (const auto& symbol : ss) {
            const uint32_t depth = turtle_stack.size();
            // operator determination, symbol to turtle command
            eval_turtle_step(symbol.RepSym, symbol.value, depth, turtle, turtle_stack);
        }
        return turtle_stack.size() == 0;
    }

    void simple_skeleton(int len) {
        Turtle turtle{};
        turtle.width = 2.0f;
        for (int i = 0; i < len; i++) {
            // apply_tropism(turtle, GravityDir, 4.22f, 1.0f);
            turtle.forward(0.5f, joints, indices);
            turtle.yaw(0.2f);
        }
    }

private:
    std::map<std::string, float> parameters;

    void apply_tropism(Turtle& turtle, const glm::vec3& T, float F, float b_l);
    void eval_turtle_step(uint32_t sym, float value, uint32_t depth, Turtle& turtle, std::stack<Turtle>& turtle_stack);
};


} // namespace ptree

#endif // OL_EVAL_H