//
// Created by lejonmcgowan on 9/24/17.
//

#ifndef JET_MYGRID_H
#define JET_MYGRID_H

#include <eigen3/Eigen/Dense>
#include "BBox.h"
#include "numUtils.h"
class MyGrid
{
public:

    MyGrid(){};
    virtual ~MyGrid(){};

    inline const Size3 res() const{ return gridRes; };
    inline const Vec3 origin() const{ return gridOrigin;};
    inline const Vec3 spacing() const{ return gridSpacing;};
    inline const BBox getBounds() const
    { return bounds; };

    Vec3 cellCenterPosition(Size3 indices) const
    {
        Vec3 centerIndices(indices[0] + 0.5,indices[1] + 0.5,indices[2] + 0.5);
        Vec3 offsets = centerIndices.cwiseProduct(gridSpacing);
        return gridOrigin + offsets;
    };

    virtual void setSize(const Size3 res, const Vec3 origin, const Vec3 spacing)
    {
        this->gridRes = res;
        this->gridOrigin = origin;
        this->gridSpacing = spacing;

        Vec3 dGridRes(gridRes[0],gridRes[1],gridRes[2]);
        dGridRes = dGridRes.cwiseProduct(gridSpacing);

        bounds = BBox(origin, origin + dGridRes);
    }

    void setSizeDefault(Size3 res)
    {
        setSize(res,Vec3(0,0,0),Vec3(1,1,1));
    }

private:
    Size3 gridRes;
    Vec3 gridSpacing = Vec3(1, 1, 1);
    //NOT the dataOrigin, which is where the point in a grid. This is the lower bounds of the bounds
    Vec3 gridOrigin = Vec3(0,0,0);
    BBox bounds = BBox(Vec3(1, 1, 1), Vec3(-1, -1, -1));

};

#endif //JET_MYGRID_H
