//
// Created by lejonmcgowan on 9/24/17.
//

#ifndef JET_MYSCALRFIELD_H
#define JET_MYSCALRFIELD_H

#include "MyGrid.h"
#include "numUtils.h"

class ScalarGrid: public MyGrid
{
public:
    ScalarGrid()
    {};

    ScalarGrid(const ScalarGrid &other)
    {
        setSize(other.res(), other.origin(), other.spacing());
        this->initialValue = other.initialValue;

        for (unsigned int i = 0; i < other.data.size(); i++)
            data.push_back(other.data[i]);
    }

    virtual ~ScalarGrid()
    {};


    virtual Size3 dataSize() const = 0;

    virtual Eigen::Vector3d dataOrigin() const = 0;

    double cubicSample(const Vec3 point) const;
    double sample(const Vec3 point) const;

    inline long getFlatIndex(Size3 indices) const {
        const Size3 dimSize = dataSize();
        return indices[0] * dimSize[1] * dimSize[2] +
                indices[1] * dimSize[2] + indices[2];}

    double at(long i, long j, long k) const;

    double at(Size3 indices) const;


    Vec3 gradient(const Eigen::Vector3d x) const;

    Vec3 gradient(const double i, const double j, const double k) const;


    double laplacian(Vec3 indices) const;

    double laplacian(int i, int j, int k) const;

    Vec3 gradientAt(long i, long j, long k) const;

    Vec3 gradientAt(const Size3 indices) const;

    double laplacianAt(Size3 indices) const;

    double laplacianAt(int i, int j, int k) const;

    //just to ensure, the entire vector parameter passed is will be OVERWRITTEN with the data from the field
    void getData(std::vector<double> *data) const;

    void fillData(const std::function<double(Vec3)> &predicate);
    void fillData(const std::function<double(double, double, double)> &predicate);
    void fillData(double value);
    void setDataAt(Size3 indices, double value);

    void foreach(std::function<void(double)> predicate) const;
    void foreachIndex(std::function<void(long, long, long)> predicate) const;

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

    void setInitValue(double value)
    { initialValue = value; }
protected:
    double initialValue;
    std::vector<double> data;
};

#endif //JET_MYSCALRFIELD_H
