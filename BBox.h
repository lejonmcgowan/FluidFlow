//
// Created by lejonmcgowan on 9/24/17.
//

#ifndef JET_BBOX3D_H
#define JET_BBOX3D_H

#include <eigen3/Eigen/Dense>
#include "numUtils.h"

class BBox
{
    Eigen::Vector3d min, max;

public:
    BBox(const Eigen::Vector3d min, const Eigen::Vector3d max)
        : min(min), max(max)
    {

    }

    BBox()
        : min(Vec3(0, 0, 0)), max(Vec3(0, 0, 0))
    {

    }

    const Eigen::Vector3d &getMin() const
    {
        return min;
    }

    const Eigen::Vector3d &getMax() const
    {
        return max;
    }
};

#endif //JET_BBOX3D_H
