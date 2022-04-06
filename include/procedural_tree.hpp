#ifndef OL_EVAL_H
#define OL_EVAL_H

#include <stack>

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

constexpr glm::vec3 GravityDir = {0, 0, -1};

using FSymbol = Symbol<float>;

void CreateSkeleton(int iterations, std::vector<glm::mat4>& joints, std::vector<uint32_t>& indices);

void Skin(const std::vector<glm::mat4>& skeleton_joints, const std::vector<uint32_t>& skeleton_indices, 
        std::vector<Vertex>& vertices, std::vector<uint32_t>& indices);

struct Turtle {
    glm::quat rotation;
    glm::vec3 position;
    float width = 0.1f;

    Turtle();

    glm::vec3 heading() const noexcept  {return glm::mat3(rotation)[2];}
    glm::vec3 up() const noexcept       {return glm::mat3(rotation)[1];}
    glm::vec3 left() const noexcept     {return glm::mat3(rotation)[0];}

    void yaw(float rad);
    void roll(float rad);
    void pitch(float rad);

    void forward(float distance);

    void level();
};

constexpr uint32_t S_forward = 1;
constexpr uint32_t S_skip    = 2;
constexpr uint32_t S_yaw     = 3; // +-
constexpr uint32_t S_pitch   = 4; // &
constexpr uint32_t S_roll    = 5; // /
constexpr uint32_t S_Dollar  = 6; // special turtle command
constexpr uint32_t S_push    = 7;
constexpr uint32_t S_pop     = 8;
constexpr uint32_t S_Bang    = 9;
constexpr uint32_t S_A       = 10; // other intermediate symbols should be above this value

class ProceduralTree {
public:

    static const char* Library[];

    void apply_tropism(Turtle& turtle, const glm::vec3& T, float F, float b_l);

    /**
     * @brief converts a given symbol string to a 3d representation. Line primitives are places in index array
     * 
     * @param ss 
     * @param joints 
     * @param indices 
     */
    template<typename T>
    bool to_skeleton(const SymbolString<T>& ss, std::vector<glm::mat4>& joints, std::vector<uint32_t>& indices) {
        std::stack<Turtle> turtle_stack;
        Turtle turtle;
        for (const auto& symbol : ss) {
            int depth = turtle_stack.size();
            // operator determination, symbol to turtle command
            switch(symbol.RepSym) {
                case S_forward:
                {
                    glm::mat4 scaling = glm::scale(glm::mat4{1.0f}, {turtle.width, turtle.width, 1.0f});

                    if (depth > 0)
                        apply_tropism(turtle, GravityDir, 1.0f / (depth + 1.0f) * 4.22f, symbol.value);

                    indices.push_back(joints.size());
                    glm::mat4 tr = glm::mat4{turtle.rotation} * scaling;
                    tr[3] = glm::vec4(turtle.position, 1.0f);
                    joints.push_back(tr);

                    turtle.forward(symbol.value);

                    indices.push_back(joints.size());
                    tr = glm::mat4{turtle.rotation} * scaling;
                    tr[3] = glm::vec4(turtle.position, 1.0f);
                    joints.push_back(tr);
                    break;
                }
                case S_skip:
                {
                    turtle.forward(symbol.value);
                    break;
                }
                case S_yaw:
                    turtle.yaw(symbol.value);
                    break;
                case S_pitch:
                    turtle.pitch(symbol.value);
                    break;
                case S_roll:
                    turtle.roll(symbol.value);
                    break;
                case S_push:
                    turtle_stack.push(turtle);
                    break;
                case S_pop:
                    turtle = turtle_stack.top();
                    turtle_stack.pop();
                    break;
                case S_Dollar:
                    turtle.level();
                    break;
                case S_Bang:
                    turtle.width = symbol.value;
                    break;
                default:
                    break; // no op for undefined symbols
            }
        }
        return turtle_stack.size() == 0;
    }
};

} // namespace ptree

#endif