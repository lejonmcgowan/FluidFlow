//
// Created by lejonmcgowan on 9/28/17.
//

#include "CollocatedGrid.h"
#include "mathUtils.h"
#include <algorithm>

double CollocatedGrid::divergenceAt(Size3 indices) const
{
    return divergenceAt(indices[0],indices[1],indices[2]);
}

double CollocatedGrid::divergenceAt(long i, long j, long k) const
{
    Size3 left, right, front, back, up, down;
    const Size3 gridRes = dataSize();

    left  << std::max(i - 1,(long)0),j,k;
    right << std::min(i + 1,gridRes[0]),j,k;

    down << i, std::max(j - 1, (long) 0), k;
    up << i, std::min(j + 1, gridRes[1]), k;

    back  << i,j,std::max(k - 1,(long)0);
    front << i,j,std::min(k + 1,gridRes[2]);

    double x, y, z;

    x = 0.5 * (at(right)[0] - at(left)[0]) / spacing()[0];
    y = 0.5 * (at(up)[1] - at(down)[1]) / spacing()[1];
    z = 0.5 * (at(front)[2] - at(back)[2]) / spacing()[2];


    return x + y + z;
}

double CollocatedGrid::divergence(const Vec3 x) const
{
    double result = 0;
    std::array<Size3,8> indices;
    std::array<double,8> weights;

    Utils::getCoordinatesAndWeights(Size3(), spacing(), Vec3(), x, indices, weights, nullptr);

    for(int i = 0; i < 8; i++)
    {
        result += weights[i] * divergenceAt(indices[i]);
    }

    return result;

}

Vec3 CollocatedGrid::at(int i, int j, int k)
{
    long idx = k + j * dataSize()[2] + i * dataSize()[1] * dataSize()[2];
    return data[idx];
}

const Vec3 CollocatedGrid::at(int i, int j, int k) const
{
    long idx = k + j * dataSize()[2] + i * dataSize()[1] * dataSize()[2];
    return data[idx];
}

const Vec3 CollocatedGrid::at(Size3 indices) const
{
    return at(indices[0],indices[1], indices[2]);
}

Vec3 CollocatedGrid::at(Size3 indices)
{
    return at(indices[0],indices[1], indices[2]);
}
Vec3 CollocatedGrid::curlAt(Size3 indices) const
{
    return curlAt(indices[0],indices[1],indices[2]);
}

Vec3 CollocatedGrid::curlAt(long i, long j, long k) const
{
    Size3 left, right, front, back, up, down;
    const Size3 gridRes = res();
    const Vec3 gridSpacing = spacing();

    left  << std::max(i - 1,(long)0),j,k;
    right << std::min(i + 1,gridRes[0]),j,k;

    up   << i,std::max(j - 1,(long)0),k;
    down << i,std::min(j + 1,gridRes[1]),k;

    back  << i,j,std::max(k - 1,(long)0);
    front << i,j,std::min(k + 1,gridRes[2]);

    double x,y,z;

    x = 0.5 * (up[2] - down[2]) / gridSpacing[1] - 0.5 * (front[1] - back[1]) / gridSpacing[2];
    y = 0.5 * (front[0] - back[0]) / gridSpacing[2] - 0.5 * (left[2] - right[2]) / gridSpacing[0];
    z = 0.5 * (left[1] - right[1]) / gridSpacing[0] - 0.5 * (up[0] - down[0]) / gridSpacing[1];

    return Vec3(x,y,z);
}

void CollocatedGrid::onResize(const Size3 &resolution, const Eigen::Vector3d gridSpacing, const Eigen::Vector3d origin,
                              Eigen::Vector3d initVal)
{
    this->data.resize(dataSize()[0] * dataSize()[1] * dataSize()[2],initVal);
}

void CollocatedGrid::setData(std::function<Vec3(double, double, double)> predicate)
{
    Size3 size = dataSize();

    data.clear();
    for(int i = 0; i < size[0]; i++)
    {
        for(int j = 0; j < size[1]; j++)
        {
            for(int k = 0; k < size[2]; k++)
            {
                Vec3 position = dataCenterPosition(Size3(i,j,k));

                data.push_back(predicate(position[0],position[1],position[2]));
            }
        }
    }
}

void CollocatedGrid::setData(std::function<Vec3(Vec3)> predicate)
{
    Size3 size = dataSize();

    data.clear();
    for(int i = 0; i < size[0]; i++)
    {
        for(int j = 0; j < size[1]; j++)
        {
            for(int k = 0; k < size[2]; k++)
            {
                Vec3 position = dataCenterPosition(Size3(i,j,k));

                data.push_back(predicate(position));
            }
        }
    }
}

void CollocatedGrid::setData(Vec3 value)
{
    Size3 size = dataSize();

    data.clear();
    for(int i = 0; i < size[0]; i++)
    {
        for(int j = 0; j < size[1]; j++)
        {
            for(int k = 0; k < size[2]; k++)
            {
                data.push_back(value);
            }
        }
    }
}

void CollocatedGrid::getData(std::vector<Vec3> &result)
{
    result.clear();
    for (unsigned int i = 0; i < data.size(); i++)
    {
        result.push_back(data[i]);
    }
}

Vec3 CollocatedGrid::sample(const Vec3 point)
{
    Vec3 result(0, 0, 0);
    std::array<Size3, 8> indices;
    std::array<double, 8> weights;

    Utils::getCoordinatesAndWeights(dataSize(), spacing(), dataOrigin(), point, indices, weights, nullptr);

    for (int i = 0; i < 8; i++)
    {
        result += weights[i] * at(indices[i]);
    }

    return result;
}
