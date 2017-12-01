//
// Created by lejonmcgowan on 9/26/17.
//

#ifndef JET_MYVECTORBUILDER_H
#define JET_MYVECTORBUILDER_H

#include <memory>
#include "grids/VectorGrid.h"

class MyVectorBuilder
{
public:
    MyVectorBuilder();

    virtual ~MyVectorBuilder();

    virtual std::shared_ptr<MyVectorBuilder> build(const Size3 resolution,const Eigen::Vector3d spacing, Eigen::Vector3d origin, Eigen::Vector3d initialValue) = 0;
};


#endif //JET_MYVECTORBUILDER_H
