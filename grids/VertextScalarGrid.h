//
// Created by lejonmcgowan on 9/24/17.
//

#ifndef JET_VERTEXTCELLGRID_H
#define JET_VERTEXTCELLGRID_H
#include "ScalarGrid.h"
class VertexScalarGrid: public ScalarGrid
{
public:
    VertexScalarGrid()
    {}

    VertexScalarGrid(const Size3 resolution,
                     const Eigen::Vector3d origin = Vec3(0, 0, 0),
                     const Eigen::Vector3d gridSpacing = Vec3(1, 1, 1),
                     double initialValue = 0.0)
    {
        setSize(resolution,origin,gridSpacing);
        this->initialValue = initialValue;
        const Size3 dimSize = dataSize();

        data.reserve(static_cast<unsigned long>(3 * dimSize[0] * dimSize[1] * dimSize[2]));
        std::fill(data.begin(),data.end(),initialValue);
    }

    VertexScalarGrid(const VertexScalarGrid &other)
        : ScalarGrid(other)
    {

    }

    Size3 dataSize() const override
    {
        return res() + Size3(1,1,1);
    }

    Vec3 dataOrigin() const override
    {
        return origin();
    }
};
#endif //JET_VERTEXTCELLGRID_H
