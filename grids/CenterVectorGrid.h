//
// Created by lejonmcgowan on 9/24/17.
//

#ifndef JET_CENTERVECTORGRID_H
#define JET_CENTERVECTORGRID_H

#include "grids/CollocatedGrid.h"
class CenterVectorGrid: public CollocatedGrid
{
public:
    Size3 dataSize() const override
    {
        return res();
    }
    Eigen::Vector3d dataOrigin() const override
    {
        return origin() + 0.5 * spacing();
    }
};

#endif //JET_CENTERVECTORGRID_H
