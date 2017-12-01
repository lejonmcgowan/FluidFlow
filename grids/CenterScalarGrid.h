//
// Created by lejonmcgowan on 9/24/17.
//

#ifndef JET_CENTERCELLGRID_H
#define JET_CENTERCELLGRID_H

#include "grids/ScalarGrid.h"
class CenterScalarGrid: public ScalarGrid
{
public:
    CenterScalarGrid()
    {

    }

    CenterScalarGrid(const Size3 resolution,
                     const Eigen::Vector3d origin,
                     const Eigen::Vector3d gridSpacing,
                     double initialValue)
    {
        setSize(resolution, origin, gridSpacing);
        this->initialValue = initialValue;
        const Size3 dimSize = dataSize();

        unsigned long size = static_cast<unsigned long>(dimSize[0] * dimSize[1] * dimSize[2]);
        data.clear();


        for (unsigned int i = 0; i < size; i++)
            data.push_back(initialValue);

    }

    CenterScalarGrid(const CenterScalarGrid &other)
        : ScalarGrid(other)
    {

    }

    virtual Size3 dataSize() const override
    {
        return res();
    }

    virtual Eigen::Vector3d dataOrigin() const override
    {
        return origin() + 0.5 * spacing();
    }
};


#endif //JET_CENTERCELLGRID_H
