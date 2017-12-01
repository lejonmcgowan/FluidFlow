//
// Created by lejonmcgowan on 9/24/17.
//

#ifndef JET_MYVECTORFIELD_H
#define JET_MYVECTORFIELD_H

#include "MyGrid.h"
#include "numUtils.h"

class VectorGrid: public MyGrid
{
public:
    VectorGrid()
    {};
    virtual ~VectorGrid()
    {};


    virtual void setSize(const Size3 res, const Vec3 origin, const Vec3 spacing) override
    {
        MyGrid::setSize(res, origin, spacing);
        onResize(res, spacing, origin, Vec3(0,0,0));
    }

    void setSize(const Size3 res, const Vec3 origin, const Vec3 spacing, Vec3 initValue)
    {
        MyGrid::setSize(res, origin, spacing);
        onResize(res, spacing, origin, initValue);
    }

    virtual Vec3 sample(const Vec3 point) = 0;

protected:
    virtual void onResize(const Size3& resolution, Eigen::Vector3d gridSpacing,
                          Eigen::Vector3d origin, Eigen::Vector3d inittialValue) = 0;


};
#endif //JET_MYVECTORFIELD_H
