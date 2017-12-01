#include <gtest/gtest.h>
#include <numUtils.h>

#include <grids/ScalarGrid.h>
#include <grids/VertextScalarGrid.h>
#include <grids/CenterScalarGrid.h>
#include <grids/CenterScalarGrid.h>
#include <grids/MACVectorGrid.h>
#include <grids/VertexVectorGrid.h>

#include <mathUtils.h>
#include <grids/CenterVectorGrid.h>

void COMPARE_VEC_INDICES(Size3 vec1, Size3 vec2, std::string name = "", double delta = 1e-10)
{
    EXPECT_TRUE(vec1.isApprox(vec2, delta)) << name << ": check between [" << vec1.transpose() << "] and ["
                                            << vec2.transpose() << "] is not approximate " << std::endl;
}

void COMPARE_VEC(Vec3 vec1, Vec3 vec2, std::string name = "", double delta = 1e-10)
{
    EXPECT_TRUE(vec1.isApprox(vec2, delta)) << name << ": check between [" << vec1.transpose() << "] and ["
                                            << vec2.transpose() << "] is not approximate " << std::endl;
}

// Test Math utils
TEST(MathTest, Barycentric)
{
    long index;
    double weight;

    //between ranges
    double value = 6.5;
    long ilow = 1, ihigh = 9;
    Utils::getBarycentric(value, ilow, ihigh, index, weight);

    EXPECT_DOUBLE_EQ(7, index);
    EXPECT_DOUBLE_EQ(0.5, weight);

    value = 4.3;
    ilow = 0, ihigh = 5;
    Utils::getBarycentric(value, ilow, ihigh, index, weight);

    EXPECT_DOUBLE_EQ(4, index);
    EXPECT_DOUBLE_EQ(0.3, weight);

    //higher
    value = 3.1;
    ilow = 3, ihigh = 5;
    Utils::getBarycentric(value, ilow, ihigh, index, weight);

    EXPECT_DOUBLE_EQ(4, index);
    EXPECT_DOUBLE_EQ(1, weight);

    //lower
    value = -0.9;
    ilow = 2, ihigh = 8;
    Utils::getBarycentric(value, ilow, ihigh, index, weight);

    EXPECT_DOUBLE_EQ(2, index);
    EXPECT_DOUBLE_EQ(0, weight);

    //clamped range
    value = 0.5;
    ilow = 2, ihigh = 2;
    Utils::getBarycentric(value, ilow, ihigh, index, weight);

    EXPECT_DOUBLE_EQ(2, index);
    EXPECT_DOUBLE_EQ(0, weight);
}

TEST(MathTest, CoordsAndWeightsSimple)
{
    Size3 res(3, 3, 3);
    Vec3 origin(0, 0, 0);
    Vec3 spacing(1, 1, 1);
    MyGrid *grid = new VertexScalarGrid(res, origin, spacing, 0);

    std::array<double, 8> weights;
    std::array<Size3, 8> indices;

    //do a smaple with equal weights
    Vec3 sample(1.5, 1.5, 1.5);
    Utils::getCoordinatesAndWeights(Size3(), *grid, Vec3(), sample, indices, weights, nullptr);

    std::array<Size3,8> answerIndices = {
            Size3(1,1,1),
            Size3(1,1,2),
            Size3(1,2,1),
            Size3(1,2,2),
            Size3(2,1,1),
            Size3(2,1,2),
            Size3(2,2,1),
            Size3(2,2,2)
    };

    for(int i = 0; i < 8; i++)
    {
        EXPECT_DOUBLE_EQ(0.125,weights[i]);
        COMPARE_VEC_INDICES(answerIndices[i],indices[i]);
    }



}

TEST(MathTest, CoordsAndWeightsOffset)
{
    Size3 res(3, 3, 3);
    //offset origin, this will change the sample result
    Vec3 origin(1, 1, 1);
    Vec3 spacing(1, 1, 1);
    MyGrid *grid = new VertexScalarGrid(res, origin, spacing, 0);

    std::array<double, 8> weights;
    std::array<Size3, 8> indices;

    //do a sample with equal weights
    Vec3 sample(1.5, 1.5, 1.5);
    Utils::getCoordinatesAndWeights(Size3(), *grid, Vec3(), sample, indices, weights, nullptr);

    std::array<Size3,8> answerIndices = {
            Size3(0,0,0),
            Size3(0,0,1),
            Size3(0,1,0),
            Size3(0,1,1),
            Size3(1,0,0),
            Size3(1,0,1),
            Size3(1,1,0),
            Size3(1,1,1)
    };

    for(int i = 0; i < 8; i++)
    {
        EXPECT_DOUBLE_EQ(0.125,weights[i]);
        COMPARE_VEC_INDICES(answerIndices[i],indices[i]);
    }

}

TEST(MathTest, CoordsAndWeightsEdge)
{
    Size3 res(3, 3, 3);
    Vec3 origin(0, 0, 0);
    Vec3 spacing(1, 1, 1);
    MyGrid *grid = new VertexScalarGrid(res, origin, spacing, 0);

    std::array<double, 8> weights;
    std::array<Size3, 8> indices;

    std::array<Size3,8> answerIndices = {
            Size3(0,0,0),
            Size3(0,0,1),
            Size3(0,1,0),
            Size3(0,1,1),
            Size3(1,0,0),
            Size3(1,0,1),
            Size3(1,1,0),
            Size3(1,1,1)
    };


    for (int i = 0; i < 2; i++)
    {
        for(int j = 0; j < 2; j++)
        {
            for(int k = 0; k < 2; k++)
            {
                //do a smaple on exact coordinates
                Vec3 sample(i,j,k);
                int flatIndex = i * 4 + j * 2 + k;
                Utils::getCoordinatesAndWeights(Size3(), *grid, Vec3(), sample, indices, weights, nullptr);

                EXPECT_DOUBLE_EQ(1,weights[0]) << "incorrect weight on flatindex " << flatIndex;
                COMPARE_VEC_INDICES(answerIndices[flatIndex],indices[0]);
            }
        }
    }
}



//Test inintialization
TEST(GridTest, InitVertexCell)
{
    std::vector<double> data;
    Size3 size(8, 8, 8);

    VertexScalarGrid field;
    field.setSizeDefault(size);

    COMPARE_VEC_INDICES(field.res(), size, "resolution init");
    COMPARE_VEC_INDICES(field.dataSize(), size + Size3(1, 1, 1), "Data Size init");
    COMPARE_VEC(field.origin(), Vec3(0, 0, 0), "Origin init");
    COMPARE_VEC(field.spacing(), Vec3(1, 1, 1), "Spacing init");

    field.setSize(size, Vec3(-1, -1, -1), Vec3(0.5, 0.5, 0.5));

    for (int i = 0; i < size[0] + 1; i++)
    {
        for (int j = 0; j < size[1] + 1; j++)
        {
            for (int k = 0; k < size[2] + 1; k++)
            {
                double di, dj, dk;

                di = field.spacing()[0] * i + field.origin()[0];
                dj = field.spacing()[1] * j + field.origin()[1];
                dk = field.spacing()[2] * k + field.origin()[2];

                Vec3 comTest = field.dataCenterPosition(Size3(i, j, k));
                comTest = comTest;

                data.push_back(std::sin(di) * std::sin(dj) * std::sin(dk));
            }
        }
    }

    COMPARE_VEC_INDICES(field.res(), size, "resolution set");
    COMPARE_VEC_INDICES(field.dataSize(), size + Size3(1, 1, 1), "Data Size set");
    COMPARE_VEC(field.origin(), Vec3(-1, -1, -1), "Origin set");
    COMPARE_VEC(field.spacing(), Vec3(0.5, 0.5, 0.5), "Spacing set");

    std::vector<double> copyData;
    auto func = [](Vec3 indices)
    { return std::sin(indices[0]) * std::sin(indices[1]) * std::sin(indices[2]); };
    field.fillData(func);
    field.getData(&copyData);

    EXPECT_EQ(data.size(), copyData.size()) << "retrieved data in Vertex Scalar grid is not the same size";

    for (int i = 0; i < size[0] + 1; i++)
    {
        for (int j = 0; j < size[1] + 1; j++)
        {
            for (int k = 0; k < size[2] + 1; k++)
            {
                Vec3 index(i, j, k);
                long flatIndex = k + (size[2] + 1) * ((size[1] + 1) * i + j);
                EXPECT_EQ(flatIndex, field.getFlatIndex(Size3(i, j, k)))
                                    << "flat intdex not consistent with grid implementation";
                EXPECT_EQ(data[flatIndex], copyData[flatIndex]) << "retreved data in Vertex Scalar grid at index " <<
                                                                index.transpose() << " is not the same";
                EXPECT_EQ(data[flatIndex], field.at(i, j, k))
                                    << "directly sampled data in Vertex Scalar grid at index " <<
                                    index.transpose() << " is not the same";
            }
        }
    }
}

TEST(GridTest, InitCenterCell)
{
    std::vector<double> data;
    Size3 size(8, 8, 8);

    CenterScalarGrid field;
    field.setSizeDefault(size);

    COMPARE_VEC_INDICES(field.res(), size, "resolution init");
    COMPARE_VEC_INDICES(field.dataSize(), size, "Data Size init");
    COMPARE_VEC(field.origin(), Vec3(0, 0, 0), "Origin init");
    COMPARE_VEC(field.spacing(), Vec3(1, 1, 1), "Spacing init");

    field.setSize(size, Vec3(-1, -1, -1), Vec3(0.5, 0.5, 0.5));

    COMPARE_VEC_INDICES(field.res(), size, "resolution set");
    COMPARE_VEC_INDICES(field.dataSize(), size, "Data Size set");
    COMPARE_VEC(field.origin(), Vec3(-1, -1, -1), "Origin set");
    COMPARE_VEC(field.spacing(), Vec3(0.5, 0.5, 0.5), "Spacing set");

    for (int i = 0; i < size[0]; i++)
    {
        for (int j = 0; j < size[1]; j++)
        {
            for (int k = 0; k < size[2]; k++)
            {
                double di, dj, dk;

                //add 0.5 to offset for cell center
                di = field.spacing()[0] * i + field.origin()[0] + 0.5 * field.spacing()[0];
                dj = field.spacing()[1] * j + field.origin()[1] + 0.5 * field.spacing()[1];
                dk = field.spacing()[2] * k + field.origin()[2] + 0.5 * field.spacing()[2];

                Vec3 comTest = field.dataCenterPosition(Size3(i, j, k));
                comTest = comTest;

                data.push_back(std::sin(di) * std::sin(dj) * std::sin(dk));
            }
        }
    }

    std::vector<double> copyData;
    field.fillData([](Vec3 indices)
                   { return std::sin(indices[0]) * std::sin(indices[1]) * std::sin(indices[2]); });
    field.getData(&copyData);

    EXPECT_EQ(data.size(), copyData.size()) << "retrieved data in centered Scalar grid is not the same size";

    for (int i = 0; i < size[0]; i++)
    {
        for (int j = 0; j < size[1]; j++)
        {
            for (int k = 0; k < size[2]; k++)
            {
                Vec3 index(i, j, k);
                unsigned long flatIndex = k + size[2] * (size[1] * i + j);
                EXPECT_EQ(data[flatIndex], copyData[flatIndex]) << "retreved data in centered Scalar grid at index " <<
                                                                index.transpose() << " is not the same";
                EXPECT_EQ(data[flatIndex], field.at(i, j, k))
                                    << "directly sampled data in centered Scalar grid at index " <<
                                    index.transpose() << " is not the same";
            }
        }
    }
}

TEST(GridTest, InitVertexVector)
{
    Size3 size(8, 8, 8);
    Size3 actualSize = size + Size3(1,1,1);
    std::vector<Vec3> data, copyData;

    VertexVectorGrid field;
    field.setSizeDefault(size);

    COMPARE_VEC_INDICES(field.res(), size, "resolution init");
    COMPARE_VEC_INDICES(field.dataSize(), actualSize, "Data Size init");
    COMPARE_VEC(field.origin(), Vec3(0, 0, 0), "Origin init");
    COMPARE_VEC(field.spacing(), Vec3(1, 1, 1), "Spacing init");

    field.setSize(size, Vec3(-1, -1, -1), Vec3(0.5, 0.5, 0.5));

    COMPARE_VEC_INDICES(field.res(), size, "resolution set");
    COMPARE_VEC_INDICES(field.dataSize(), actualSize, "Data Size set");
    COMPARE_VEC(field.origin(), Vec3(-1, -1, -1), "Origin set");
    COMPARE_VEC(field.spacing(), Vec3(0.5, 0.5, 0.5), "Spacing set");

    field.setData([](double i, double j, double k) -> Vec3
                  {
                     return Vec3(i,j,k);
                  });

    for(int i = 0; i < actualSize[0]; i++)
    {
        for(int j = 0; j < actualSize[1]; j++)
        {
            for(int k = 0; k < actualSize[2]; k++)
            {
                Vec3 result = field.origin() + Vec3(i,j,k).cwiseProduct(field.spacing());
                data.push_back(result);
            }
        }
    }

    field.getData(copyData);

    for(int i = 0; i < actualSize[0]; i++)
    {
        for(int j = 0; j < actualSize[1]; j++)
        {
            for(int k = 0; k < actualSize[2]; k++)
            {
                int flatIndex = k + j * actualSize[2] + i * actualSize[2] * actualSize[1];
                COMPARE_VEC(data[flatIndex],copyData[flatIndex],"VERTEX VECTOR compare retrieved data and exptected result");
                COMPARE_VEC(copyData[flatIndex], field.at(i,j,k),"VERTEX VECTOR compare retrieved data and queried result");
            }
        }
    }

}

TEST(GridTest, InitCenterVector)
{
    Size3 size(8, 8, 8);
    Size3 actualSize = size;
    std::vector<Vec3> data, copyData;

    CenterVectorGrid field;
    field.setSizeDefault(size);

    COMPARE_VEC_INDICES(field.res(), size, "resolution init");
    COMPARE_VEC_INDICES(field.dataSize(), actualSize, "Data Size init");
    COMPARE_VEC(field.origin(), Vec3(0, 0, 0), "Origin init");
    COMPARE_VEC(field.spacing(), Vec3(1, 1, 1), "Spacing init");

    field.setSize(size, Vec3(-1, -1, -1), Vec3(0.5, 0.5, 0.5));

    COMPARE_VEC_INDICES(field.res(), size, "resolution set");
    COMPARE_VEC_INDICES(field.dataSize(), actualSize, "Data Size set");
    COMPARE_VEC(field.origin(), Vec3(-1, -1, -1), "Origin set");
    COMPARE_VEC(field.spacing(), Vec3(0.5, 0.5, 0.5), "Spacing set");

    field.setData([](Vec3 indices) -> Vec3
                  {
                      return indices;
                  });

    const Vec3 spacing = field.spacing();
    for(int i = 0; i < actualSize[0]; i++)
    {
        for(int j = 0; j < actualSize[1]; j++)
        {
            for(int k = 0; k < actualSize[2]; k++)
            {
                Vec3 result = field.origin() + spacing/2 +  Vec3(i,j,k).cwiseProduct(spacing);
                data.push_back(result);
            }
        }
    }

    field.getData(copyData);

    for(int i = 0; i < actualSize[0]; i++)
    {
        for(int j = 0; j < actualSize[1]; j++)
        {
            for(int k = 0; k < actualSize[2]; k++)
            {
                int flatIndex = k + j * actualSize[2] + i * actualSize[2] * actualSize[1];
                COMPARE_VEC(data[flatIndex],copyData[flatIndex],"CENTER VECTOR compare retrieved data and exptected result");
                COMPARE_VEC(copyData[flatIndex], field.at(i,j,k),"CENTER VECTOR compare retrieved data and queried result");
            }
        }
    }

}

TEST(GridTest, InitMACVector)
{
    Size3 size(3, 4, 5);

    double actualSizeU = 4 * 4 * 5;
    double actualSizeV = 3 * 5 * 5;
    double actualSizeW = 3 * 4 * 6;

    std::vector<double> dataU,dataV,dataW, copyDataU,copyDataV,copyDataW;

    MACVectorGrid field;
    field.setSizeDefault(size);

    field.fillData([](Vec3 indices) -> Vec3
                   {
                       return indices;
                   });


    COMPARE_VEC_INDICES(field.res(), size, "resolution init");
    EXPECT_DOUBLE_EQ(field.sizeU(), actualSizeU) <<  "MAC GRID: Data Size U init did not match";
    EXPECT_DOUBLE_EQ(field.sizeV(), actualSizeV) <<  "MAC GRID: Data Size V init did not match";
    EXPECT_DOUBLE_EQ(field.sizeW(), actualSizeW) <<  "MAC GRID: Data Size W init did not match";
    COMPARE_VEC(field.origin(), Vec3(0, 0, 0), "Origin init");
    COMPARE_VEC(field.spacing(), Vec3(1, 1, 1), "Spacing init");

    const Vec3 spacing = field.spacing();
    Size3 dims = field.getUDims();

    field.forEachU([&](long i, long j , long k)
                   {
                      dataU.push_back(i);
                   });
    field.forEachV([&](long i, long j , long k)
                   {
                       dataV.push_back(j);
                   });
    field.forEachW([&](long i, long j , long k)
                   {
                       dataW.push_back(k);
                   });

    field.getDataU(copyDataU);
    field.getDataV(copyDataV);
    field.getDataW(copyDataW);

    for(int i = 0; i < size[0] + 1; i++)
    {
        for(int j = 0; j < size[1]; j++)
        {
            for(int k = 0; k < size[2]; k++)
            {
                long flatIndex = field.flatIndexU(i,j,k);
                EXPECT_DOUBLE_EQ(dataU[flatIndex],copyDataU[flatIndex]) << "MAC GRID: U data retrieved not equal at index " << Vec3(i,j,k).transpose();
                EXPECT_DOUBLE_EQ(copyDataU[flatIndex],field.u(i,j,k))   << "MAC GRID: U data queried not equal at index " << Vec3(i,j,k).transpose();
            }
        }
    }

    for(int i = 0; i < size[0]; i++)
    {
        for(int j = 0; j < size[1] + 1; j++)
        {
            for(int k = 0; k < size[2]; k++)
            {
                long flatIndex = field.flatIndexV(i,j,k);
                EXPECT_DOUBLE_EQ(dataV[flatIndex],copyDataV[flatIndex]) << "MAC GRID: V data retrieved not equal at index " << Vec3(i,j,k).transpose();
                EXPECT_DOUBLE_EQ(copyDataV[flatIndex],field.v(i,j,k))   << "MAC GRID: V data queried not equal at index " << Vec3(i,j,k).transpose();
            }
        }
    }

    for(int i = 0; i < size[0]; i++)
    {
        for(int j = 0; j < size[1]; j++)
        {
            for(int k = 0; k < size[2] + 1; k++)
            {
                long flatIndex = field.flatIndexW(i,j,k);
                EXPECT_DOUBLE_EQ(dataW[flatIndex],copyDataW[flatIndex]) << "MAC GRID: W data retrieved not equal at index " << Vec3(i,j,k).transpose();
                EXPECT_DOUBLE_EQ(copyDataW[flatIndex],field.w(i,j,k))   << "MAC GRID: W data queried not equal at index " << Vec3(i,j,k).transpose();
            }
        }
    }

}

TEST(GridTest, InitMACVector2)
{
    Size3 size(3, 4, 5);

    double actualSizeU = 4 * 4 * 5;
    double actualSizeV = 3 * 5 * 5;
    double actualSizeW = 3 * 4 * 6;

    std::vector<double> dataU,dataV,dataW, copyDataU,copyDataV,copyDataW;

    MACVectorGrid field;
    field.setSize(size,Vec3(-1,-1,-1),Vec3(2,3,0.5));

    field.fillData([](double i, double j, double k) -> Vec3
                   {
                       return Vec3(i, j, k);
                   });


    COMPARE_VEC_INDICES(field.res(), size, "resolution init");
    EXPECT_DOUBLE_EQ(field.sizeU(), actualSizeU) <<  "MAC GRID: Data Size U init did not match";
    EXPECT_DOUBLE_EQ(field.sizeV(), actualSizeV) <<  "MAC GRID: Data Size V init did not match";
    EXPECT_DOUBLE_EQ(field.sizeW(), actualSizeW) <<  "MAC GRID: Data Size W init did not match";
    COMPARE_VEC(field.origin(), Vec3(-1,-1,-1), "Origin init");
    COMPARE_VEC(field.spacing(), Vec3(2, 3, 0.5), "Spacing init");


    field.forEachU([&](long i, long j, long k)
                   {
                       dataU.push_back(i * 2 - 1);
                   });

    field.forEachV([&](long i, long j, long k)
                   {
                       dataV.push_back(j * 3 - 1);
                   });

    field.forEachW([&](long i, long j, long k)
                   {
                       dataW.push_back(k / 2.0 - 1);
                   });

    field.getDataU(copyDataU);
    field.getDataV(copyDataV);
    field.getDataW(copyDataW);


    Size3 dims = field.getUDims();
    for(int i = 0; i < dims[0]; i++)
    {
        for(int j = 0; j < dims[1]; j++)
        {
            for(int k = 0; k < dims[2]; k++)
            {
                long flatIndex = field.flatIndexU(i,j,k);
                EXPECT_DOUBLE_EQ(dataU[flatIndex],copyDataU[flatIndex]) << "MAC GRID: U data retrieved not equal at index " << Vec3(i,j,k).transpose();
                EXPECT_DOUBLE_EQ(copyDataU[flatIndex],field.u(i,j,k))   << "MAC GRID: U data queried not equal at index " << Vec3(i,j,k).transpose();
            }
        }
    }

    dims = field.getVDims();

    for(int i = 0; i < dims[0]; i++)
    {
        for(int j = 0; j < dims[1]; j++)
        {
            for(int k = 0; k < dims[2]; k++)
            {
                long flatIndex = field.flatIndexV(i,j,k);
                EXPECT_DOUBLE_EQ(dataV[flatIndex],copyDataV[flatIndex]) << "MAC GRID: V data retrieved not equal at index " << Vec3(i,j,k).transpose();
                EXPECT_DOUBLE_EQ(copyDataV[flatIndex],field.v(i,j,k))   << "MAC GRID: V data queried not equal at index " << Vec3(i,j,k).transpose();
            }
        }
    }

    dims = field.getWDims();
    for(int i = 0; i < dims[0]; i++)
    {
        for(int j = 0; j < dims[1]; j++)
        {
            for(int k = 0; k < dims[2]; k++)
            {
                long flatIndex = field.flatIndexW(i,j,k);
                EXPECT_DOUBLE_EQ(dataW[flatIndex],copyDataW[flatIndex]) << "MAC GRID: W data retrieved not equal at index " << Vec3(i,j,k).transpose();
                EXPECT_DOUBLE_EQ(copyDataW[flatIndex],field.w(i,j,k))   << "MAC GRID: W data queried not equal at index " << Vec3(i,j,k).transpose();
            }
        }
    }

}

TEST(GridTest, MACCenter)
{
    Size3 size(3, 4, 5);

    double actualSizeU = 3 * 5 * 6;
    double actualSizeV = 4 * 4 * 6;
    double actualSizeW = 4 * 5 * 5;

    std::vector<double> dataU, dataV, dataW, copyDataU, copyDataV, copyDataW;

    MACVectorGrid field;
    field.setSizeDefault(size);

    field.fillData(1);


    for (int i = 0; i < size[0] - 1; i++)
    {
        for (int j = 0; j < size[1] - 1; j++)
        {
            for (int k = 0; k < size[2] - 1; k++)
            {
                int flatIndex = k + j * (size[2] + 1) + i * (size[2] + 1) * (size[1] + 1);
                COMPARE_VEC(Vec3(1, 1, 1), field.atCenter(i, j, k), "MAC GRID: Center sample not correct");
            }
        }
    }

}


//Test Scalar calculus
TEST(GridTest, GradientCenterCell)
{
    Size3 size(5, 8, 6);
    Vec3 dx(2, 3, 1.5);

    CenterScalarGrid field;
    field.setSize(size, Vec3(0, 0, 0), dx);
    field.fillData([](Vec3 indices) -> double
                   { return indices[0] + 2 * indices[1] - 3 * indices[2]; });

    //homogeneous grid, all [non-extrapolated] data point should have the same gradient: <1,2,3>
    for (int i = 0; i < size[0]; i++)
    {
        for (int j = 0; j < size[1]; j++)
        {
            for (int k = 0; k < size[2]; k++)
            {
                //expected partial given normal circumstance
                Vec3 ans(1, 2, -3);
                bool edge = false;
                /**
                 * extrapolated values, since they handle border cases by using the respective end border,
                 * which is equal to the partial derivative of the respective component / 2
                */
                if (i == 0 || i == size[0] - 1)
                {
                    ans[0] /= 2;
                    edge = true;
                }
                if (j == 0 || j == size[1] - 1)
                {
                    ans[1] /= 2;
                    edge = true;
                }
                if (k == 0 || k == size[2] - 1)
                {
                    ans[2] /= 2;
                    edge = true;
                }

                COMPARE_VEC(field.gradientAt(i, j, k), ans, "CELL-CENTERED SCALAR: gradient comparison");
                //trivial laplacian test; should always be zero
                if (!edge)
                    EXPECT_EQ(field.laplacianAt(i, j, k), 0) << "VERTEX-CENTERED SCALAR: trivial laplacian comparison";
            }
        }
    }
}

TEST(GridTest, GradientVertextCell)
{
    Size3 size(5, 8, 6);
    Vec3 dx(2, 3, 1.5);

    VertexScalarGrid field;
    field.setSize(size, Vec3(0, 0, 0), dx);
    field.fillData([](Vec3 indices)
                   { return indices[0] + 2 * indices[1] - 3 * indices[2]; });

    //homogeneous grid, all [non-extrapolated] data point should have the same gradient: <1,2,3>
    for (int i = 0; i < size[0] + 1; i++)
    {
        for (int j = 0; j < size[1] + 1; j++)
        {
            for (int k = 0; k < size[2] + 1; k++)
            {
                //expected partial given normal circumstance
                Vec3 ans(1, 2, -3);
                bool edge = false;
                /**
                 * extrapolated values, since they handle border cases by using the respective end border,
                 * which is equal to the partial derivative of the respective component / 2
                */
                if (i == 0 || i == size[0])
                {
                    ans[0] /= 2;
                    edge = true;
                }
                if (j == 0 || j == size[1])
                {
                    ans[1] /= 2;
                    edge = true;
                }
                if (k == 0 || k == size[2])
                {
                    ans[2] /= 2;
                    edge = true;
                }

                COMPARE_VEC(field.gradientAt(i, j, k), ans, "VERTEX-CENTERED SCALAR: gradient comparison");
                //trivial laplacian test; should always be zero
                if (!edge)
                    EXPECT_EQ(field.laplacianAt(i, j, k), 0) << "VERTEX-CENTERED SCALAR: trivial laplacian comparison";
            }
        }
    }
}

TEST(GridTest, LaplacianVertextCell)
{
    Size3 size(7, 4, 5);

    VertexScalarGrid field1, field2, field3, field4;
    Vec3 dx(0.75, 1, 2);
    Vec3 orig(2, -3.4, 1);
    field1.setSize(size, orig, dx);
    field2.setSize(size, orig, dx);
    field3.setSize(size, orig, dx);
    field4.setSize(size, orig, dx);

    //test separate components
    //2x^3
    field1.fillData([](double i, double j, double k)
                    {
                        return 2 * (i * i * i);
                    });
    //3y^3
    field2.fillData([](double i, double j, double k)
                    {
                        return 3 * (j * j * j);
                    });
    //-2z^3
    field3.fillData([](double i, double j, double k)
                    {
                        return -2 * (k * k * k);
                    });

    //test that takes all components into account

    //x^3 + 4y^3 - 2z^3
    field4.fillData([](double i, double j, double k)
                    {
                        return i * i * i + 4 * j * j * j - 2 * k * k * k;
                    });

    //homogeneous grid, all [non-extrapolated] data point should have the same laplacian:
    //12x, the ccoordinates are offset by the center spacing
    //todo once again, no check for boundary points, that get ugly and I'm too lazy to calulate it.
    for (int i = 1; i < size[0] - 1; i++)
    {
        for (int j = 1; j < size[1] - 1; j++)
        {
            for (int k = 1; k < size[2] - 1; k++)
            {
                //same size and origin, so this should apply to all graphs
                Vec3 position = field1.dataCenterPosition(i, j, k);

                double x = orig[0] + i * dx[0];
                double y = orig[1] + j * dx[1];
                double z = orig[2] + k * dx[2];

                //expected partial given normal circumstance
                double ans1 = 12 * x;
                double ans2 = 18 * y;
                double ans3 = -12 * z;
                double ans4 = 6 * x + 24 * y - 12 * z;

                double laplacian1 = field1.laplacianAt(i, j, k);
                double laplacian2 = field2.laplacianAt(i, j, k);
                double laplacian3 = field3.laplacianAt(i, j, k);
                double laplacian4 = field4.laplacianAt(i, j, k);


                EXPECT_NEAR(laplacian1, ans1, 1e-6)
                                    << "VERTEX-CENTERED SCALAR: X   laplacian at " << Size3(i, j, k).transpose();
                EXPECT_NEAR(laplacian2, ans2, 1e-6)
                                    << "VERTEX-CENTERED SCALAR: Y   laplacian at " << Size3(i, j, k).transpose();
                EXPECT_NEAR(laplacian3, ans3, 1e-6)
                                    << "VERTEX-CENTERED SCALAR: Z   laplacian at " << Size3(i, j, k).transpose();
                EXPECT_NEAR(laplacian4, ans4, 1e-6)
                                    << "VERTEX-CENTERED SCALAR: ALL laplacian at " << Size3(i, j, k).transpose();
            }
        }
    }
}

TEST(GridTest, LaplacianCenterCell)
{
    Size3 size(6, 5, 7);

    CenterScalarGrid field1, field2, field3, field4;
    Vec3 dx(0.75, 1, 2);
    Vec3 orig(2, -3.4, 1);
    field1.setSize(size, orig, dx);
    field2.setSize(size, orig, dx);
    field3.setSize(size, orig, dx);
    field4.setSize(size, orig, dx);

    //test separate components
    //2x^3
    field1.fillData([](double i, double j, double k)
                    {
                        return 2 * (i * i * i);
                    });
    //3y^3
    field2.fillData([](double i, double j, double k)
                    {
                        return 3 * (j * j * j);
                    });
    //-2z^3
    field3.fillData([](double i, double j, double k)
                    {
                        return -2 * (k * k * k);
                    });

    //test that takes all components into account

    //x^3 + 4y^3 - 2z^3
    field4.fillData([](double i, double j, double k)
                    {
                        return i * i * i + 4 * j * j * j - 2 * k * k * k;
                    });

    //homogeneous grid, all [non-extrapolated] data point should have the same laplacian:
    //12x, the ccoordinates are offset by the center spacing
    //todo once again, no check for boundary points, that get ugly and I'm too lazy to calulate it.
    for (int i = 1; i < size[0] - 1; i++)
    {
        for (int j = 1; j < size[1] - 1; j++)
        {
            for (int k = 1; k < size[2] - 1; k++)
            {
                //same size and origin, so this should apply to all graphs
                Vec3 position = field1.dataCenterPosition(i, j, k);

                double x = orig[0] + i * dx[0] + 0.5 * dx[0];
                double y = orig[1] + j * dx[1] + 0.5 * dx[1];
                double z = orig[2] + k * dx[2] + 0.5 * dx[2];

                //expected partial given normal circumstance
                double ans1 = 12 * x;
                double ans2 = 18 * y;
                double ans3 = -12 * z;
                double ans4 = 6 * x + 24 * y - 12 * z;

                double laplacian1 = field1.laplacianAt(i, j, k);
                double laplacian2 = field2.laplacianAt(i, j, k);
                double laplacian3 = field3.laplacianAt(i, j, k);
                double laplacian4 = field4.laplacianAt(i, j, k);


                EXPECT_NEAR(laplacian1, ans1, 1e-6)
                                    << "CELL-CENTERED SCALAR: X   laplacian at " << Size3(i, j, k).transpose();
                EXPECT_NEAR(laplacian2, ans2, 1e-6)
                                    << "CELL-CENTERED SCALAR: Y   laplacian at " << Size3(i, j, k).transpose();
                EXPECT_NEAR(laplacian3, ans3, 1e-6)
                                    << "CELL-CENTERED SCALAR: Z   laplacian at " << Size3(i, j, k).transpose();
                EXPECT_NEAR(laplacian4, ans4, 1e-6)
                                    << "CELL-CENTERED SCALAR: ALL laplacian at " << Size3(i, j, k).transpose();
            }
        }
    }
}

TEST(GridTest, DivergenceVertex)
{
    Size3 size(4, 3, 5);

    VertexVectorGrid field;
    field.setSizeDefault(size);

    field.setData([](double i, double j, double k) -> Vec3
                  {
                      return Vec3(i, j, k);
                  }
    );


    for (int i = 1; i < size[0] - 1; i++)
    {
        for (int j = 1; j < size[1] - 1; j++)
        {
            for (int k = 1; k < size[2] - 1; k++)
            {
                EXPECT_DOUBLE_EQ(3, field.divergenceAt(i, j, k));
            }
        }
    }
}

TEST(GridTest, DivergenceCenter)
{
    Size3 size(4, 3, 5);

    CenterVectorGrid field;
    field.setSizeDefault(size);

    field.setData([](double i, double j, double k) -> Vec3
                  {
                      return Vec3(i, j, k);
                  }
    );


    for (int i = 1; i < size[0] - 1; i++)
    {
        for (int j = 1; j < size[1] - 1; j++)
        {
            for (int k = 1; k < size[2] - 1; k++)
            {
                EXPECT_DOUBLE_EQ(3, field.divergenceAt(i, j, k));
            }
        }
    }
}

TEST(GridTest, DivergenceMAC)
{
    Size3 size(4, 3, 5);

    MACVectorGrid field;
    field.setSizeDefault(size);

    field.fillData([](double i, double j, double k) -> Vec3
                   {
                       return Vec3(i, j, k);
                   }
    );


    for (int i = 1; i < size[0] - 1; i++)
    {
        for (int j = 1; j < size[1] - 1; j++)
        {
            for (int k = 1; k < size[2] - 1; k++)
            {
                EXPECT_DOUBLE_EQ(3, field.divergenceAt(i, j, k));
            }
        }
    }
}

TEST(GridTest, CurlCenter)
{
    Size3 size(4, 3, 5);

    CenterVectorGrid field;
    field.setSizeDefault(size);

    //trivial curl: curl of this shuold be zero
    field.setData([](double i, double j, double k) -> Vec3
                  {
                      return Vec3(i, 2 * j, 3 * k);
                  }
    );


    for (int i = 0; i < size[0]; i++)
    {
        for (int j = 0; j < size[1]; j++)
        {
            for (int k = 0; k < size[2]; k++)
            {
                COMPARE_VEC(Vec3(0, 0, 0), field.curlAt(i, j, k));
            }
        }
    }
}

TEST(GridTest, CurlMAC)
{
    Size3 size(4, 3, 5);

    MACVectorGrid field;
    field.setSizeDefault(size);

    //trivial test; curl of zero
    field.fillData([](double i, double j, double k) -> Vec3
                   {
                       return Vec3(i, j, k);
                   }
    );

    //not checking boundaries here
    for (int i = 1; i < size[0] - 1; i++)
    {
        for (int j = 1; j < size[1] - 1; j++)
        {
            for (int k = 1; k < size[2] - 1; k++)
            {
                COMPARE_VEC(Vec3(0, 0, 0), field.curlAt(i, j, k));
            }
        }
    }
}


