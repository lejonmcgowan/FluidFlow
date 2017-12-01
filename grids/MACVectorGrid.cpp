//
// Created by lejonmcgowan on 9/30/17.
//


#include "MACVectorGrid.h"
#include "mathUtils.h"

Vec3 MACVectorGrid::sample(const Vec3 point)
{
   return Utils::trilinearInterpolateFromGrid(*this,point);
}

Vec3 MACVectorGrid::cubicSample(const Vec3 point)
{
    return Utils::cubicInterpolateFromGrid(*this,point);
}

double MACVectorGrid::divergenceAt(long i, long j, long k) const
{
    double left, right, front, back, up, down;

    left = u(i, j, k);
    right = u(i + 1, j, k);

    down = v(i, j, k);
    up = v(i, j + 1, k);

    back = w(i, j, k);
    front = w(i, j, k + 1);

    double x, y, z;

    x = (right - left) / spacing()[0];
    y = (up - down) / spacing()[1];
    z = (front - back) / spacing()[2];


    return x + y + z;
}

double MACVectorGrid::divergenceAt(Size3 indices) const
{
    return divergenceAt(indices[0], indices[1], indices[2]);
}

double MACVectorGrid::divergence(const Vec3 x) const
{
    double result = 0.0;
    std::array<Size3, 8> indices;
    std::array<double, 8> weights;

    Utils::getCoordinatesAndWeights(res() - Size3(1, 1, 1),
                                    spacing(),
                                    origin() + 0.5 * spacing(),
                                    x,
                                    indices,
                                    weights,
                                    nullptr);

    for (int i = 0; i < 8; i++)
    {
        result += weights[i] * divergenceAt(indices[i]);
    }

    return result;

}

Vec3 MACVectorGrid::curlAt(Size3 indices) const
{
    return curlAt(indices[0], indices[1], indices[2]);
}

Vec3 MACVectorGrid::curlAt(long i, long j, long k) const
{
    Vec3 left, right, front, back, up, down;
    const Size3 gridRes = res();
    const Vec3 gridSpacing = spacing();

    left = atCenter(std::max(i - 1, (long) 0), j, k);
    right = atCenter(std::min(i + 1, gridRes[0] - 1), j, k);

    up = atCenter(i, std::max(j - 1, (long) 0), k);
    down = atCenter(i, std::min(j + 1, gridRes[1] - 1), k);

    back = atCenter(i, j, std::max(k - 1, (long) 0));
    front = atCenter(i, j, std::min(k + 1, gridRes[2] - 1));

    double x, y, z;

    x = 0.5 * (up[2] - down[2]) / gridSpacing[1] - 0.5 * (front[1] - back[1]) / gridSpacing[2];
    y = 0.5 * (front[0] - back[0]) / gridSpacing[2] - 0.5 * (left[2] - right[2]) / gridSpacing[0];
    z = 0.5 * (left[1] - right[1]) / gridSpacing[0] - 0.5 * (up[0] - down[0]) / gridSpacing[1];

    return Vec3(x, y, z);
}

Vec3 MACVectorGrid::curl(const Vec3 point) const
{
    Vec3 result(0, 0, 0);
    std::array<Size3, 8> indices;
    std::array<double, 8> weights;

    Utils::getCoordinatesAndWeights(res() - Size3(1, 1, 1),
                                    spacing(),
                                    origin() + 0.5 * spacing(),
                                    point,
                                    indices,
                                    weights,
                                    nullptr);

    for (int i = 0; i < 8; i++)
    {
        result += weights[i] * curlAt(indices[i]);
    }

    return result;
}

Vec3 MACVectorGrid::atCenter(long i, long j, long k) const
{
    const Vec3 sample = Vec3(u(i, j, k) + u(i + 1, j, k), v(i, j, k) + v(i, j + 1, k),
                             w(i, j, k) + w(i, j, k + 1));
    return 0.5 * sample;
}

void MACVectorGrid::onResize(const Size3 &resolution, const Eigen::Vector3d gridSpacing, const Eigen::Vector3d origin,
                             Eigen::Vector3d initialValue)
{
    dataU.clear();
    dataV.clear();
    dataW.clear();

    unsigned long sizeX = static_cast<unsigned long>(getUDims()[0] * getUDims()[1] * getUDims()[2]);
    unsigned long sizeY = static_cast<unsigned long>(getVDims()[0] * getVDims()[1] * getVDims()[2]);
    unsigned long sizeZ = static_cast<unsigned long>(getWDims()[0] * getWDims()[1] * getWDims()[2]);

    dataU.resize(sizeX);
    dataV.resize(sizeY);
    dataW.resize(sizeZ);

    std::fill(dataU.begin(),dataU.end(),initialValue[0]);
    std::fill(dataV.begin(),dataV.end(),initialValue[1]);
    std::fill(dataW.begin(),dataW.end(),initialValue[2]);

    double x = gridSpacing[0];
    double y = gridSpacing[1];
    double z = gridSpacing[2];

    originU = origin + 0.5 * Vec3(0, y, z);
    originV = origin + 0.5 * Vec3(x, 0, z);
    originW = origin + 0.5 * Vec3(x, y, 0);
}

double MACVectorGrid::divergence(long i, long j, long k) const
{
    return divergence(Vec3(i, j, k));
}

void MACVectorGrid::fillData(std::function<Vec3(Vec3)> predicate)
{
    dataU.clear();
    dataV.clear();
    dataW.clear();

    const Size3 r = res() + Size3(1, 1, 1);
    for (int i = 0; i < r[0]; i++)
    {
        for (int j = 0; j < r[1]; j++)
        {
            for (int k = 0; k < r[2]; k++)
            {
                Vec3 positionU = dataCenterPosU(i, j, k);
                if(j < r[1] - 1 && k < r[2] - 1)
                    dataU.push_back(predicate(positionU)[0]);

                Vec3 positionV = dataCenterPosV(i, j, k);
                if(i < r[0] - 1 && k < r[2] - 1)
                    dataV.push_back(predicate(positionV)[1]);

                Vec3 positionW = dataCenterPosW(i, j, k);
                if(i < r[0] - 1 && j < r[1] - 1)
                    dataW.push_back(predicate(positionW)[2]);

            }
        }
    }
}

void MACVectorGrid::fillData(std::function<Vec3(double, double, double)> predicate)
{
    dataU.clear();
    dataV.clear();
    dataW.clear();

    const Size3 r = res() + Size3(1, 1, 1);
    for (int i = 0; i < r[0]; i++)
    {
        for (int j = 0; j < r[1]; j++)
        {
            for (int k = 0; k < r[2]; k++)
            {
                Vec3 positionU = dataCenterPosU(i, j, k);
                if(j < r[1] - 1 && k < r[2] - 1)
                    dataU.push_back(predicate(positionU[0], positionU[1], positionU[2])[0]);

                Vec3 positionV = dataCenterPosV(i, j, k);
                if(i < r[0] - 1 && k < r[2] - 1)
                    dataV.push_back(predicate(positionV[0], positionV[1], positionV[2])[1]);

                Vec3 positionW = dataCenterPosW(i, j, k);
                if(i < r[0] - 1 && j < r[1] - 1)
                    dataW.push_back(predicate(positionW[0], positionW[1], positionW[2])[2]);

            }
        }
    }
}

void MACVectorGrid::fillData(double value)
{
    dataU.clear();
    dataV.clear();
    dataW.clear();

    const Size3 r = res() + Size3(1, 1, 1);
    for (int i = 0; i < r[0]; i++)
    {
        for (int j = 0; j < r[1]; j++)
        {
            for (int k = 0; k < r[2]; k++)
            {
                if(j < r[1] - 1 && k < r[2] - 1)
                    dataU.push_back(value);

                if(i < r[0] - 1 && k < r[2] - 1)
                    dataV.push_back(value);

                if(i < r[0] - 1 && j < r[1] - 1)
                    dataW.push_back(value);

            }
        }
    }
}

void MACVectorGrid::getDataU(std::vector<double> &result) const
{
    result.clear();
    for (double value: dataU)
        result.push_back(value);
}

void MACVectorGrid::getDataV(std::vector<double> &result) const
{
    result.clear();
    for (double value: dataV)
        result.push_back(value);
}

void MACVectorGrid::getDataW(std::vector<double> &result) const
{
    result.clear();
    for (double value: dataW)
        result.push_back(value);
}

Vec3 MACVectorGrid::dataCenterPosW(const Vec3 indices) const
{
    return originW + Vec3(spacing().array() * indices.array());
}

Vec3 MACVectorGrid::dataCenterPosW(long i, long j, long k) const
{
    return dataCenterPosW(Vec3(i, j, k));
}

Vec3 MACVectorGrid::dataCenterPosV(const Vec3 indices) const
{
    return originV + Vec3(spacing().array() * indices.array());
}

Vec3 MACVectorGrid::dataCenterPosV(long i, long j, long k) const
{
    return dataCenterPosV(Vec3(i, j, k));
}

Vec3 MACVectorGrid::dataCenterPosU(const Vec3 indices) const
{
    return originU + Vec3(spacing().array() * indices.array());
}

Vec3 MACVectorGrid::dataCenterPosU(long i, long j, long k) const
{
    return dataCenterPosU(Vec3(i, j, k));
}
double MACVectorGrid::sizeU()
{
    return dataU.size();
}

double MACVectorGrid::sizeV()
{
    return dataV.size();
}

double MACVectorGrid::sizeW()
{
    return dataW.size();
}

const double MACVectorGrid::u(long i, long j, long k) const
{
    return dataU[flatIndexU(i,j,k)];
}

const double MACVectorGrid::v(long i, long j, long k) const
{
    return dataV[flatIndexV(i,j,k)];
}

const double MACVectorGrid::w(long i, long j, long k) const
{
    return dataW[flatIndexW(i,j,k)];
}

double MACVectorGrid::u(long i, long j, long k)
{
    return dataU[flatIndexU(i,j,k)];
}

double MACVectorGrid::v(long i, long j, long k)
{
    return dataV[flatIndexV(i,j,k)];
}

double MACVectorGrid::w(long i, long j, long k)
{
    return dataW[flatIndexW(i,j,k)];
}
void MACVectorGrid::setUDataAt(long i, long j, long k, double value)
{
    long flatIdx = flatIndexU(i, j, k);
    dataU[flatIdx] = value;
}

void MACVectorGrid::setVDataAt(long i, long j, long k, double value)
{
    long flatIdx = flatIndexV(i, j, k);
    dataV[flatIdx] = value;
}

void MACVectorGrid::setWDataAt(long i, long j, long k, double value)
{
    long flatIdx = flatIndexW(i, j, k);;
    dataW[flatIdx] = value;
}

void MACVectorGrid::setUDataAt(Size3 indices, double value)
{
    setUDataAt(indices[0], indices[1], indices[2], value);
}

void MACVectorGrid::setVDataAt(Size3 indices, double value)
{
    setVDataAt(indices[0], indices[1], indices[2], value);
}

void MACVectorGrid::setWDataAt(Size3 indices, double value)
{
    setWDataAt(indices[0], indices[1], indices[2], value);
}

void MACVectorGrid::forEachU(std::function<void(long, long, long)> operation)
{
    for (int i = 0; i < getUDims()[0]; i++)
    {
        for (int j = 0; j < getUDims()[1]; j++)
        {
            for (int k = 0; k < getUDims()[2]; k++)
            {
                operation(i, j, k);
            }
        }
    }
}

void MACVectorGrid::forEachV(std::function<void(long, long, long)> operation)
{
    for (int i = 0; i < getVDims()[0]; i++)
    {
        for (int j = 0; j < getVDims()[1]; j++)
        {
            for (int k = 0; k < getVDims()[2]; k++)
            {
                operation(i, j, k);
            }
        }
    }
}

void MACVectorGrid::forEachW(std::function<void(long, long, long)> operation)
{
    for (int i = 0; i < getWDims()[0]; i++)
    {
        for (int j = 0; j < getWDims()[1]; j++)
        {
            for (int k = 0; k < getWDims()[2]; k++)
            {
                operation(i, j, k);
            }
        }
    }
}



