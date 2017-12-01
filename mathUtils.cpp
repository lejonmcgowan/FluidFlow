//
// Created by lejonmcgowan on 10/1/17.
//

#include "mathUtils.h"
#include <grids/ScalarGrid.h>
#include <grids/MACVectorGrid.h>

void Utils::getCoordinatesAndWeights(const Size3 res,
                                     const Vec3 spacing,
                                     const Vec3 origin,
                                     const Eigen::Vector3d x,
                                     std::array<Size3, 8> &indices,
                                     std::array<double, 8> &weights,
                                     Vec3 *triweights)
{
    long i, j, k;
    double fx, fy, fz;

    assert(spacing[0] > 0.0 && spacing[1] > 0.0 && spacing[2] > 0.0);

    Eigen::Vector3d normalizedX = (x - origin);

    normalizedX[0] /= spacing[0];
    normalizedX[1] /= spacing[1];
    normalizedX[2] /= spacing[2];

    const ssize_t iSize = static_cast<ssize_t>(res[0]);
    const ssize_t jSize = static_cast<ssize_t>(res[1]);
    const ssize_t kSize = static_cast<ssize_t>(res[2]);

    getBarycentric(normalizedX[0], 0, iSize - 1, i, fx);
    getBarycentric(normalizedX[1], 0, jSize - 1, j, fy);
    getBarycentric(normalizedX[2], 0, kSize - 1, k, fz);

    if (triweights != nullptr)
    {
        *triweights << fx, fy, fz;
    }

    const ssize_t ip1 = std::min(i + 1, iSize - 1);
    const ssize_t jp1 = std::min(j + 1, jSize - 1);
    const ssize_t kp1 = std::min(k + 1, kSize - 1);

    indices[0] = Size3(i, j, k);
    indices[1] = Size3(i, j, kp1);
    indices[2] = Size3(i, jp1, k);
    indices[3] = Size3(i, jp1, kp1);
    indices[4] = Size3(ip1, j, k);
    indices[5] = Size3(ip1, j, kp1);
    indices[6] = Size3(ip1, jp1, k);
    indices[7] = Size3(ip1, jp1, kp1);

    weights[0] = (1 - fx) * (1 - fy) * (1 - fz);
    weights[1] = (1 - fx) * (1 - fy) * fz;
    weights[2] = (1 - fx) * fy * (1 - fz);
    weights[3] = (1 - fx) * fy * fz;
    weights[4] = fx * (1 - fy) * (1 - fz);
    weights[5] = fx * (1 - fy) * fz;
    weights[6] = fx * fy * (1 - fz);
    weights[7] = fx * fy * fz;
}

double Utils::trilinearInterpolateFromGrid(const ScalarGrid &grid, const Vec3 point)
{
    std::array<Size3, 8> indices{};
    std::array<double, 8> weights{};
    std::array<double, 8> values{};
    Vec3 triweights;
    getCoordinatesAndWeights(grid.dataSize(), grid.spacing(), grid.dataOrigin(), point, indices, weights, &triweights);

    for (int i = 0; i < 8; i++)
    {
        values[i] = grid.at(indices[i]);
    }

    return trilinearInterpolate(values, triweights);
}

Vec3 Utils::trilinearInterpolateFromGrid(const MACVectorGrid &grid, const Vec3 point)
{
    std::array<double, 8> values{};
    Vec3 uTriWeights,vTriWeights,wTriWeights;

    std::array<Size3, 8> indicesU, indicesV, indicesW;
    std::array<double, 8> weightsU, weightsV, weightsW;

    Utils::getCoordinatesAndWeights(grid.getUDims(), grid.spacing(), grid.getUOrigin(), point, indicesU, weightsU, &uTriWeights);
    Utils::getCoordinatesAndWeights(grid.getVDims(), grid.spacing(), grid.getVOrigin(), point, indicesV, weightsV, &vTriWeights);
    Utils::getCoordinatesAndWeights(grid.getWDims(), grid.spacing(), grid.getWOrigin(), point, indicesW, weightsW, &wTriWeights);

    for (int i = 0; i < 8; i++)
    {
        values[i] = grid.u(indicesU[i][0],indicesU[i][1],indicesU[i][2]);
    }
    double x = trilinearInterpolate(values, uTriWeights);

    for (int i = 0; i < 8; i++)
    {
        values[i] = grid.v(indicesV[i][0],indicesV[i][1],indicesV[i][2]);
    }
    double y = trilinearInterpolate(values, vTriWeights);

    for (int i = 0; i < 8; i++)
    {
        values[i] = grid.w(indicesW[i][0],indicesW[i][1],indicesW[i][2]);
    }

    double z = trilinearInterpolate(values, wTriWeights);

    return Vec3(x,y,z);
}

double Utils::trilinearInterpolate(std::array<double, 8> values, Vec3 weights)
{
    std::array<double, 4> values1, values2;
    Vec2 biweights(weights[2], weights[1]);
    for (int i = 0; i < 4; i++)
    {
        values1[i] = values[i];
    }
    for (int i = 0; i < 4; i++)
    {
        values2[i] = values[i + 4];
    }

    return lerp(bilinearInterpolate(values1, biweights), bilinearInterpolate(values2, biweights), weights[0]);
}
double Utils::bilinearInterpolate(std::array<double, 4> values, Vec2 weights)
{
    return lerp(lerp(values[0], values[1], weights[0]),
                lerp(values[2], values[3], weights[0]),
                weights[1]);
}
double Utils::lerp(double x, double y, double weight)
{
    return (1 - weight) * x + weight * y;
}
void Utils::extrapolateToRegion(const std::vector<double> &input,
                                const Size3 dims,
                                const std::vector<bool> &validMatrix,
                                unsigned int numIters,
                                std::vector<double> &output)
{
    assert(static_cast<unsigned int>(dims[0] * dims[1] * dims[2]) == input.size());

    std::vector<bool> valid0(validMatrix.size());
    std::vector<bool> valid1(validMatrix.size());

    for (unsigned int i = 0; i < validMatrix.size(); i++)
    {
        valid0[i] = validMatrix[i];
        output[i] = input[i];
    }

    const long iIter = dims[1] * dims[2];
    const long jIter = dims[2];
    const long kIter = 1;

    for (unsigned int iter = 0; iter < numIters; ++iter)
    {
        for (int i = 0; i < dims[0]; i++)
        {
            for (int j = 0; j < dims[1]; j++)
            {
                for (int k = 0; k < dims[2]; k++)
                {
                    long flatIndex = i * dims[1] * dims[2] + j * dims[2] + k;

                    double sum = 0;

                    unsigned int count = 0;

                    if (!valid0[flatIndex])
                    {
                        if (i + 1 < dims[0] && valid0[flatIndex + iIter])
                        {
                            sum += output[flatIndex + iIter];
                            ++count;
                        }

                        if (i > 0 && valid0[flatIndex - iIter])
                        {
                            sum += output[flatIndex - iIter];
                            ++count;
                        }

                        if (j + 1 < dims[1] && valid0[flatIndex + jIter])
                        {
                            sum += output[flatIndex + jIter];
                            ++count;
                        }

                        if (j > 0 && valid0[flatIndex - jIter])
                        {
                            sum += output[flatIndex - jIter];
                            ++count;
                        }

                        if (k + 1 < dims[2] && valid0[flatIndex + kIter])
                        {
                            sum += output[flatIndex + kIter];
                            ++count;
                        }

                        if (k > 0 && valid0[flatIndex - kIter])
                        {
                            sum += output[flatIndex - kIter];
                            ++count;
                        }

                        if (count > 0)
                        {
                            output[flatIndex] = sum / (double) count;
                            valid1[flatIndex] = true;
                        }
                    }
                    else
                    {
                        valid1[flatIndex] = true;
                    }
                }

            }
        }

        std::swap(valid0, valid1);
    }
}
double Utils::pseudoTimeStep(const ScalarGrid &sdf, const Vec3 gridSpacing, double maxCFL)
{
    //const Size3 size = sdf.dataSize();

    const double hMax = std::max(gridSpacing[0], std::max(gridSpacing[1], gridSpacing[2]));
    const double hMin = std::min(gridSpacing[0], std::min(gridSpacing[1], gridSpacing[2]));

    double initS = sdf.at(0, 0, 0);
    double maxS = initS / std::sqrt(initS * initS + hMin * hMin);
    double dtau = maxCFL * hMax;
    sdf.foreach([&](double value)
                {
                    maxS = std::max(maxS, value / std::sqrt(value * value + hMin * hMin));
                });


    while (dtau * maxS / hMax > maxCFL)
    {
        dtau *= 0.5;
    }

    return dtau;
}
double Utils::cubicInterpolateFromGrid(const ScalarGrid &grid, const Vec3 point)
{

    std::array<Size3, 8> indices{};
    std::array<double, 8> weights{};
    std::array<double, 8> values{};
    Vec3 triweights;
    getCoordinatesAndWeights(grid.dataSize(), grid.spacing(), grid.dataOrigin(), point, indices, weights, &triweights);


    long i,j,k;

    //indices to start the sample from are at the lowest corner of bounding cube
    i = indices[0][0];
    j = indices[0][1];
    k = indices[0][2];
    Size3 size = grid.dataSize();
    long sizeI = size[0];
    long sizeJ = size[1];
    long sizeK = size[2];
    //setup the indices needed for catmull-rom
    std::array<long,4> iPoints = {
        std::max(i - 1, (long)0),
        i,
        std::min(i + 1, sizeI - 1),
        std::min(i + 2, sizeI - 1)
    };
    std::array<long,4> jPoints = {
        std::max(j - 1, (long)0),
        j,
        std::min(j + 1, sizeJ - 1),
        std::min(j + 2, sizeJ - 1)
    };
    std::array<long,4> kPoints = {
        std::max(k - 1, (long)0),
        k,
        std::min(k + 1, sizeK - 1),
        std::min(k + 2, sizeK - 1)
    };

    std::array<double, 4> kValues{};

    //perform 3-dimensional catmull-rom interpolation
    for (int catK = 0; catK < 4; catK++)
    {
        std::array<double, 4> jValues{};


        for (int catJ = 0; catJ < 4; catJ++)
        {
            std::array<double,4> iValues = {
                grid.at(iPoints[0], jPoints[catJ], kPoints[catK]),
                grid.at(iPoints[1], jPoints[catJ], kPoints[catK]),
                grid.at(iPoints[2], jPoints[catJ], kPoints[catK]),
                grid.at(iPoints[3], jPoints[catJ], kPoints[catK]),
            };
            jValues[catJ] = catmullRom(iValues,triweights[0]);
        }

        kValues[catK] = catmullRom(jValues,triweights[1]);
    }

    return catmullRom(kValues, triweights[2]);
}

double helperMacCubicInterpolate(Size3 res, Vec3 spacing, Vec3 origin, const std::vector<double>& data, Vec3 point)
{
    std::array<Size3, 8> indices{};
    std::array<double, 8> weights{};
    std::array<double, 8> values{};
    Vec3 triweights;
    Utils::getCoordinatesAndWeights(res, spacing, origin, point, indices, weights, &triweights);


    long i,j,k;

    //indices to start the sample from are at the lowest corner of bounding cube
    i = indices[0][0];
    j = indices[0][1];
    k = indices[0][2];
    long sizeI = res[0];
    long sizeJ = res[1];
    long sizeK = res[2];
    //setup the indices needed for catmull-rom
    std::array<long,4> iPoints = {
        std::max(i - 1, (long)0),
        i,
        std::min(i + 1, sizeI - 1),
        std::min(i + 2, sizeI - 1)
    };
    std::array<long,4> jPoints = {
        std::max(j - 1, (long)0),
        j,
        std::min(j + 1, sizeJ - 1),
        std::min(j + 2, sizeJ - 1)
    };
    std::array<long,4> kPoints = {
        std::max(k - 1, (long)0),
        k,
        std::min(k + 1, sizeK - 1),
        std::min(k + 2, sizeK - 1)
    };

    std::array<double, 4> kValues{};

    const int kIter = 1;
    const int jIter = res[2];
    const int iIter = res[2] * res[1];
    //perform 3-dimensional catmull-rom interpolation
    for (int catK = 0; catK < 4; catK++)
    {
        std::array<double, 4> jValues{};


        for (int catJ = 0; catJ < 4; catJ++)
        {
            double flatIndexm1 = iPoints[0] * iIter + jPoints[catJ] * jIter + kPoints[catK] * kIter;
            double flatIndex0  = iPoints[1] * iIter + jPoints[catJ] * jIter + kPoints[catK] * kIter;
            double flatIndex1  = iPoints[2] * iIter + jPoints[catJ] * jIter + kPoints[catK] * kIter;
            double flatIndex2  = iPoints[3] * iIter + jPoints[catJ] * jIter + kPoints[catK] * kIter;

            std::array<double,4> iValues = {
                data[flatIndexm1],
                data[flatIndex0],
                data[flatIndex1],
                data[flatIndex2],
            };
            jValues[catJ] = Utils::catmullRom(iValues,triweights[0]);
        }

        kValues[catK] = Utils::catmullRom(jValues,triweights[1]);
    }

    return Utils::catmullRom(kValues, triweights[2]);
}

Vec3 Utils::cubicInterpolateFromGrid(const MACVectorGrid& grid, const Vec3 point)
{
    Vec3 spacing = grid.spacing();

    std::vector<double> uData, vData, wData;
    grid.getDataU(uData);
    grid.getDataV(vData);
    grid.getDataW(wData);

    double x = helperMacCubicInterpolate(grid.getUDims(),spacing,grid.getUOrigin(),uData,point);
    double y = helperMacCubicInterpolate(grid.getVDims(),spacing,grid.getVOrigin(),vData,point);
    double z = helperMacCubicInterpolate(grid.getWDims(),spacing,grid.getWOrigin(),wData,point);
    return Vec3(x,y,z);
}

double Utils::catmullRom(std::array<double, 4> values, double x)
{
    double d0 = (values[2] - values[0]) / 2.0;
    double d1 = (values[3] - values[1]) / 2.0;
    double midDiff = values[2] - values[1];

    //take into account monotonicity by settting derivatives at extremes (when the slope changes) to 0

    if(std::fabs(midDiff) < std::numeric_limits<double>::epsilon())
    {
        d0 = d1 = 0;
    }
    if(midDiff * d1 < 0.0)
    {
        d1 = 0;
    }
    if(midDiff * d0 < 0.0)
    {
        d0 = 0;
    }

    double a,b,c,d;

    a = d1 + d0 - 2 * midDiff;
    b = 3 * midDiff - 2 * d0 - d1;
    c = d0;
    d = values[1];

    return a * x * x * x + b * x * x + c * x + d;
}
