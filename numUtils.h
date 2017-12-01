//
// Created by lejonmcgowan on 9/24/17.
//

#ifndef JET_NUMUTILS_H
#define JET_NUMUTILS_H
#include <eigen3/Eigen/Dense>

typedef Eigen::Array<long,3,1> Size3;

typedef Eigen::Vector3d Vec3;
typedef Eigen::Vector2d Vec2;

typedef Eigen::Vector4d Vec4;

typedef Eigen::Matrix4d Mat4;

enum Medium
{
    FLUID = 0,
    AIR = 1,
    BOUNDARY = 2
};

class Ray
{
public:
    Ray(Vec3 origin, Vec3 direction)
        : o(origin), d(direction.normalized())
    {}
    inline Vec3 origin()
    { return o; }
    inline Vec3 direction()
    { return d; }
    inline Vec3 timeAlongRay(double t)
    { return o + t * d; }
private:

    const Vec3 o, d;

};

#endif //JET_NUMUTILS_H
