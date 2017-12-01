//
// Created by lejonmcgowan on 9/26/17.
//

#include "ScalarGrid.h"
#include "mathUtils.h"
#include <algorithm>
#include <iostream>

double ScalarGrid::cubicSample(const Vec3 point) const
{
    return Utils::cubicInterpolateFromGrid(*this,point);
}

double ScalarGrid::sample(const Vec3 point) const
{
    double result = 0.0;
    std::array<Size3,8> indices;
    std::array<double,8> weights;

    Utils::getCoordinatesAndWeights(dataSize(),spacing(),dataOrigin(),point,indices,weights);

    for(int i = 0; i < 8; i++)
    {
        result += weights[i] * at(indices[i]);
    }

    return result;
}

Vec3 ScalarGrid::gradientAt(long i, long j, long k) const
{
    const Size3 &dimSize = dataSize();
    assert(i >= 0 && i < dimSize[0]);
    assert(j >= 0 && j < dimSize[1]);
    assert(k >= 0 && k < dimSize[2]);

    Eigen::Vector3i is,js,ks;

    //extrapolate values outside grid to the index itself
    is << std::max(i-1,(long)0),i,std::min(i+1,dimSize[0] - 1);
    js << std::max(j-1,(long)0),j,std::min(j+1,dimSize[1] - 1);
    ks << std::max(k-1,(long)0),k,std::min(k+1,dimSize[2] - 1);

    double left   = at(is[0],j,k);
    double right  = at(is[2],j,k);

    double bottom = at(i,js[0],k);
    double top    = at(i,js[2],k);

    double back   = at(i,j,ks[0]);
    double front  = at(i,j,ks[2]);

    //using central difference
    double x = (right - left) / 2 / spacing()[0];
    double y = (top - bottom) / 2 / spacing()[1];
    double z = (front - back) / 2 / spacing()[2];

    return Vec3(x,y,z);
}
Vec3 ScalarGrid::gradientAt(Size3 indices) const
{
    return gradientAt(indices[0],indices[1],indices[2]);
}
Vec3 ScalarGrid::gradient(const Vec3 x) const
{
    Vec3 result(0,0,0);
    std::array<Size3,8> indices;
    std::array<double,8> weights;

    Utils::getCoordinatesAndWeights(Size3(), spacing(), Vec3(), x, indices, weights, nullptr);

    for(int i = 0; i < 8; i++)
    {
        result += weights[i] * gradientAt(indices[i]);
    }

    return result;

}
Vec3 ScalarGrid::gradient(double i, double j, double k) const
{
    return gradient(Vec3(i,j,k));
}
double ScalarGrid::laplacianAt(Size3 indices) const
{
    return laplacianAt(indices[0],indices[1],indices[2]);
}
double ScalarGrid::laplacianAt(int i, int j, int k) const
{
    const double center = at(i,j,k);
    double dLeft,dRight,dUp,dDown,dFront,dBack;
    //any samples outside the domain (samples taken at dge) becomes 0
    dLeft = dRight = dUp = dDown = dFront = dBack = 0;

    if(i > 0)
    {
        dLeft = at(i-1,j,k);
    }
    if(i + 1 < dataSize()[0])
    {
        dRight = at(i+1,j,k);
    }

    if(j > 0)
    {
        dDown = at(i,j-1,k);
    }
    if(j + 1 < dataSize()[1])
    {
        dUp = at(i,j+1,k);
    }

    if(k > 0)
    {
        dBack = at(i,j,k-1);
    }
    if(k + 1 < dataSize()[2])
    {
        dFront = at(i,j,k+1);
    }
    double ppx = (dRight + dLeft - 2 * center);
    double ppy = (dDown + dUp    - 2 * center);
    double ppz = (dFront + dBack - 2 * center);
    return ppx / spacing()[0] / spacing()[0] +
           ppy / spacing()[1] / spacing()[1] +
           ppz / spacing()[2] / spacing()[2];
}
double ScalarGrid::laplacian(Vec3 x) const
{
    double result = 0.0;
    std::array<Size3,8> indices;
    std::array<double,8> weights;

    Utils::getCoordinatesAndWeights(Size3(), spacing(), Vec3(), x, indices, weights, nullptr);

    for(int i = 0; i < 8; i++)
    {
        result += weights[i] * laplacianAt(indices[i]);
    }

    return result;

}
double ScalarGrid::laplacian(int i, int j, int k) const
{
    return laplacian(Vec3(i,j,k));
}
double ScalarGrid::at(long i, long j, long k) const
{
    return data[getFlatIndex(Size3(i,j,k))];
}
double ScalarGrid::at(Size3 indices) const
{
    return at(indices[0],indices[1],indices[2]);
}


void ScalarGrid::getData(std::vector<double> *data) const
{
    Size3 dSize = dataSize();
    unsigned long numElements = dSize[0] * dSize[1] * dSize[2];

    data->clear();
    data->resize(numElements);

    std::copy(this->data.begin(),this->data.end(),data->begin());

}

void ScalarGrid::fillData(const std::function<double(Vec3)> &predicate)
{
    Size3 dSize = dataSize();
    long numElements = dSize[0] * dSize[1] * dSize[2];

    this->data.clear();
    this->data.reserve(numElements);
    for(int i = 0; i < dSize[0]; i++)
    {
        for(int j = 0; j < dSize[1]; j++)
        {
            for(int k = 0; k < dSize[2]; k++)
            {
                data.push_back(predicate(dataCenterPosition(Size3(i,j,k))));
            }
        }
    }
}

void ScalarGrid::fillData(const std::function<double(double, double, double)> &predicate)
{
    Size3 dSize = dataSize();
    long numElements = dSize[0] * dSize[1] * dSize[2];

    this->data.clear();
    this->data.reserve(numElements);
    for(int i = 0; i < dSize[0]; i++)
    {
        for(int j = 0; j < dSize[1]; j++)
        {
            for(int k = 0; k < dSize[2]; k++)
            {
                const Vec3 &dataPosition = dataCenterPosition(Size3(i, j, k));
                data.push_back(predicate(dataPosition[0],dataPosition[1],dataPosition[2]));
            }
        }
    }
}

void ScalarGrid::fillData(double value)
{
    Size3 dSize = dataSize();
    long numElements = dSize[0] * dSize[1] * dSize[2];

    data.resize(numElements);
    std::fill(data.begin(),data.end(),value);
}

void ScalarGrid::setDataAt(Size3 indices, double value)
{
    assert(indices[0] >= 0 && indices[0] < dataSize()[0]);
    assert(indices[1] >= 0 && indices[1] < dataSize()[1]);
    assert(indices[2] >= 0 && indices[2] < dataSize()[2]);

    data[getFlatIndex(indices)] = value;
}
void ScalarGrid::foreach(std::function<void(double)> predicate) const
{
    for(unsigned int i = 0; i < data.size();i++)
        predicate(data[i]);
}
void ScalarGrid::foreachIndex(std::function<void(long,long,long)> predicate) const
{
    Size3 size = dataSize();

    for(int i = 0; i < size[0]; i++)
    {
        for(int j = 0; j < size[1]; j++)
        {
            for(int k = 0; k < size[2]; k++)
            {
                predicate(i,j,k);
            }
        }
    }
}
