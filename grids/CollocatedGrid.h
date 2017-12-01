//
// Created by lejonmcgowan on 9/24/17.
//

#ifndef JET_COLLOCATEDGRID_H
#define JET_COLLOCATEDGRID_H
#include "VectorGrid.h"
class CollocatedGrid: public VectorGrid
{
public:
    CollocatedGrid(){};
    virtual ~CollocatedGrid(){};

    CollocatedGrid(const CollocatedGrid &other)
    {
        setSize(other.res(), other.origin(), other.spacing());
        for (unsigned int i = 0; i < other.data.size(); i++)
        {
            data.push_back(other.data[i]);
        }
    }

    virtual Size3 dataSize() const = 0;
    virtual Vec3 dataOrigin() const = 0;

    Vec3 at(int i, int j, int k);
    const Vec3 at(int i, int j, int k) const;
    const Vec3 at(Size3 indices) const;
    Vec3 at(Size3 indices);

    double divergenceAt(long i, long j, long k) const;
    double divergenceAt(Size3 indices) const;
    double divergence(const Vec3 x) const;

    Vec3 curlAt(Size3 indices) const;
    Vec3 curlAt(long i, long j, long k) const;

    void setData(std::function<Vec3(double,double,double)> predicate);
    void setData(std::function<Vec3(Vec3)> predicate);
    void setData(Vec3 value);

    void getData(std::vector<Vec3>& result);

    Vec3 dataCenterPosition(Size3 indices) const
    {
        Vec3 indexPositions(indices[0], indices[1], indices[2]);
        Vec3 offsets = indexPositions.cwiseProduct(spacing());
        return dataOrigin() + offsets;
    }

    Vec3 dataCenterPosition(long i, long j, long k) const
    {
        Vec3 indexPositions(i,j,k);
        Vec3 offsets = indexPositions.cwiseProduct(spacing());
        return dataOrigin() + offsets;
    }

    Vec3 sample(const Vec3 point) override;

private:
    void onResize(const Size3 &resolution, Vec3 gridSpacing, Vec3 origin, Vec3 initVal) override;

private:

    std::vector<Vec3> data;
};
#endif //JET_COLLOCATEDGRID_H
