/**
 * @file spline.h
 * @brief 
 * @date 2022-04-22
 * 
 * from https://github.com/recp/cglm/blob/8f09cc858377b3692f3931dfd1d3fa3719551112/include/cglm/bezier.h
 */
#ifndef PTREE_SPLINE_H
#define PTREE_SPLINE_H

#include <glm/glm.hpp>

namespace ptree {

constexpr glm::mat4 GLM_BEZIER_MAT = {{-1.0f,  3.0f, -3.0f,  1.0f},
                                    { 3.0f, -6.0f,  3.0f,  0.0f},
                                    {-3.0f,  3.0f,  0.0f,  0.0f},
                                    { 1.0f,  0.0f,  0.0f,  0.0f}};


constexpr glm::mat4 GLM_HERMITE_MAT = {{ 2.0f, -3.0f,  0.0f,  1.0f},
                                    {-2.0f,  3.0f,  0.0f,  0.0f},
                                    { 1.0f, -2.0f,  1.0f,  0.0f},
                                    { 1.0f, -1.0f,  0.0f,  0.0f}};

/**
 * @brief cubic bezier interpolation function
 * 
 * @tparam T 
 * @param s interpolation value [0, 1]
 * @param p0 point 0
 * @param c0 control point 0
 * @param c1 control point 1
 * @param p1 point 1
 * @return T 
 */
template<typename T>
T bezier(float s, const T& p0, const T& c0, const T& c1, const T& p1) noexcept {
    const glm::vec4 wv = glm::vec4{s*s*s, s*s, s, 1} * GLM_BEZIER_MAT; // {(1-s)^3, s(1-s)^2, s^2(1-s), s^3}
    return p0 * wv[0] + c0 * wv[1] + c1 * wv[2] + p1 * wv[3];
}

template<typename T>
T hermite(float s, const T& p0, const T& t0, const T& t1, const T& p1) noexcept {
    float ss, d, a, b, c, e, f;

    ss = s  * s;
    a  = ss + ss;
    c  = a  + ss;
    b  = a  * s;
    d  = s  * ss;
    f  = d  - ss;
    e  = b  - c;

    return p0 * (e + 1.0f) + t0 * (f - ss + s) + t1 * f - p1 * e;
}

/**
 * @brief interpolate data [n-1, n, n+1, n+2] from B to C
 * 
 * @tparam T 
 * @param A n-1
 * @param B n
 * @param C n+1
 * @param D n+2
 * @param t [0, 1]
 * @return T 
 */
template<typename T>
T cubic_hermite (const T& A, const T& B, const T& C, const T& D, float t)
{
    T a = 0.5f * (-A + 3.0f*B - 3.0f*C + D);
    T b = A - 0.5f * (5.0f*B) + 2.0f*C - 0.5f * D;
    T c = 0.5f * (-A + C);
    T d = B;
 
    return a*t*t*t + b*t*t + c*t + d;
}

template<typename T>
T cubic_hermite_dt (const T& A, const T& B, const T& C, const T& D, float t)
{
    T a = 0.5f * (-A + 3.0f*B - 3.0f*C + D);
    T b = A - (5.0f*B)/2.0f + 2.0f*C - D / 2.0f;
    T c = 0.5f * (-A + C);
 
    return 3.0f * a*t*t + 2.0f * b*t + c;
}

template<typename T>
struct HermiteSpline {
    std::vector<T> points;

    T eval(float t) const noexcept {
        const float x = (points.size() - 1) * t;
        const int index = x;
        
        return cubic_hermite<T>(
            points[glm::clamp<int>(index - 1, 0, points.size() - 1)],
            points[glm::clamp<int>(index    , 0, points.size() - 1)],
            points[glm::clamp<int>(index + 1, 0, points.size() - 1)],
            points[glm::clamp<int>(index + 2, 0, points.size() - 1)],
            x - floor(x)
        );
    }

    T eval_dt(float t) const noexcept {
        const float x = (points.size() - 1) * t;
        const int index = x;
        
        return cubic_hermite_dt<T>(
            points[glm::clamp(index - 1, 0, points.size() - 1)],
            points[glm::clamp(index    , 0, points.size() - 1)],
            points[glm::clamp(index + 1, 0, points.size() - 1)],
            points[glm::clamp(index + 2, 0, points.size() - 1)],
            x - floor(x)
        );
    }
};

}

#endif // PTREE_SPLINE_H