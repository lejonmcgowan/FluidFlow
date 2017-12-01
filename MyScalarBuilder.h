//
// Created by lejonmcgowan on 9/26/17.
//

#ifndef JET_MYSCALARBUILDER_H
#define JET_MYSCALARBUILDER_H

#include <memory>
#include "grids/ScalarGrid.h"
class MyScalarBuilder
{
public:
    MyScalarBuilder();

    virtual ~MyScalarBuilder();

    virtual std::shared_ptr<ScalarGrid> build(const Size3 resolution,const Eigen::Vector3d spacing, Eigen::Vector3d origin, double initialValue) = 0;
};
#endif //JET_MYSCALARBUILDER_H
