//
// Created by lejonmcgowan on 10/23/17.
//
#include "PressureSolver.h"

#include<eigen3/Eigen/IterativeLinearSolvers>

PressureSolver::PressureSolver()
{

}
PressureSolver::~PressureSolver()
{

}

void PressureSolver::solve(const MACVectorGrid &input,
                           const ScalarGrid &boundarySDF,
                           const ScalarGrid &fluidSDF,
                           MACVectorGrid &output)
{

    setupMarkers(input, boundarySDF, fluidSDF);
    setupSystem(input);

//    std::cout << "System b" << std::endl;
//    for (int k = 0; k < input.res()[2]; k++)
//    {
//        for (int j = 0; j < input.res()[1]; j++)
//        {
//            for (int i = 0; i < input.res()[0]; i++)
//            {
//                long flatIndex = i * input.res()[2] * input.res()[1] + j * input.res()[2] + k;
//                std::cout << "(" << i << "," << j << "," << k << "): " << x[flatIndex] << std::endl;
//            }
//        }
//    }


    //Symmetric matrix, so the solver only has to resolve the upper part (right, up, front). bottom half isn't even
    //filled out
    Eigen::ConjugateGradient<Eigen::SparseMatrix<double>,Eigen::Upper > solver;
    x = solver.compute(A).solve(b);



    applyPressureGradient(input, output);

}
void PressureSolver::setupMarkers(const MACVectorGrid &input,
                                  const ScalarGrid &boundary,
                                  const ScalarGrid &fluid)
{
    Size3 size = input.res();

    markers.clear();
    markers.resize(size[0] * size[1] * size[2]);
    for (int k = 0; k < size[2]; k++)
    {
        for (int j = 0; j < size[1]; j++)
        {
            for (int i = 0; i < size[0]; i++)
            {
                long flatIndex = k + j * size[2] + i * size[2] * size[1];

                Vec3 pos = input.cellCenterPosition(Size3(i, j, k));
                double fluidSample = fluid.sample(pos);
                if (boundary.sample(pos) < 0.0)
                    markers[flatIndex] = (BOUNDARY);
                else if (fluidSample < 0.0)
                {
                    markers[flatIndex] = (FLUID);
                }
                else
                    markers[flatIndex] = (AIR);
            }
        }
    }
//    std::cout << "Markers" << std::endl;
//    for (int k = 0; k < size[2]; k++)
//    {
//        for (int j = 0; j < size[1]; j++)
//        {
//            for (int i = 0; i < size[0]; i++)
//            {
//                int flatIndex = i * size[2] * size[1] + j * size[2] + k;
//                std::cout << "(" << i << "," << j << "," << k << "): " << markers[flatIndex] << std::endl;
//            }
//        }
//    }
//    std::cout << std::endl;
}
void PressureSolver::setupSystem(const MACVectorGrid &input)
{
    typedef Eigen::Triplet<double> DataPoint;

    Size3 res = input.res();
    long size = res[0] * res[1] * res[2];
    A = Eigen::SparseMatrix<double>(size, size);
    std::vector<DataPoint> Adata;
    b = Eigen::VectorXd::Zero(size);
    x.resize(size);

    Vec3 c(1 / (input.spacing()[0] * input.spacing()[0]),
           1 / (input.spacing()[1] * input.spacing()[1]),
           1 / (input.spacing()[2] * input.spacing()[2]));


    int count = 0;
    //buil matrix along diagonal
    for (int diagIndex = 0; diagIndex < A.rows(); diagIndex++)
    {
        int tempI = diagIndex;
        //unflatten indices
        long unflatI = tempI / (res[1] * res[2]);
        tempI -= unflatI * res[1] * res[2];
        long unflatJ = tempI / res[2];
        long unflatK = tempI % res[2];

        Size3 gridIndices(unflatI, unflatJ, unflatK);
        long kIter = 1, jIter = res[2], iIter = res[2] * res[1];



        if (markers[diagIndex] == FLUID)
        {
            count++;
            double value = input.divergenceAt(gridIndices);
            b[diagIndex] = value;

            //in x dimension
            if (unflatI + 1 < res[0] && markers[diagIndex + iIter] != BOUNDARY)
            {
                Adata.emplace_back(diagIndex, diagIndex, c[0]);
                if (markers[diagIndex + iIter] == FLUID)
                {
                    Adata.emplace_back(diagIndex, diagIndex + iIter, -c[0]);
                }
            }

            if (unflatI > 0 && markers[diagIndex - iIter] != BOUNDARY)
            {
                Adata.emplace_back(diagIndex, diagIndex, c[0]);
            }
            //in y dimension
            if (unflatJ + 1 < res[1] && markers[diagIndex + jIter] != BOUNDARY)
            {
                Adata.emplace_back(diagIndex, diagIndex, c[1]);
                if (markers[diagIndex + jIter] == FLUID)
                {
                    Adata.emplace_back(diagIndex, diagIndex + jIter, -c[1]);
                }
            }

            if (unflatJ > 0 && markers[diagIndex - jIter] != BOUNDARY)
            {
                Adata.emplace_back(diagIndex, diagIndex, c[1]);
            }
            //in z dimension
            if (unflatK + 1 < res[2] && markers[diagIndex + kIter] != BOUNDARY)
            {
                Adata.emplace_back(diagIndex, diagIndex, c[2]);
                if (markers[diagIndex + kIter] == FLUID)
                {
                    Adata.emplace_back(diagIndex, diagIndex + kIter, -c[2]);
                }
            }

            if (unflatK > 0 && markers[diagIndex - kIter] != BOUNDARY)
            {
                Adata.emplace_back(diagIndex, diagIndex, c[2]);

            }

        }
        else
        {
            Adata.emplace_back(diagIndex, diagIndex, 1);
        }
    }

    A.setFromTriplets(Adata.begin(),Adata.end());


}

void PressureSolver::applyPressureGradient(const MACVectorGrid &input, MACVectorGrid &output)
{
    Size3 res = input.res();
    Vec3 invSpacing(1 / input.spacing()[0], 1 / input.spacing()[1], 1 / input.spacing()[2]);

    for (int i = 0; i < res[0]; i++)
    {
        for (int j = 0; j < res[1]; j++)
        {
            for (int k = 0; k < res[2]; k++)
            {
                long flatIndex = k + j * res[2] + i * res[2] * res[1];
                double kIter = 1, jIter = res[2], iIter = res[2] * res[1];
                if (markers[flatIndex] == FLUID)
                {
                    if (i + 1 < res[0] && markers[flatIndex + iIter] != BOUNDARY)
                    {
                        double uData = input.u(i + 1, j, k);
                        double x1 = x[flatIndex + iIter];
                        double x0 = x[flatIndex];
                        output.setUDataAt(i + 1, j, k, uData + invSpacing[0] * (x1 - x0));
                    }

                    if (j + 1 < res[1] && markers[flatIndex + jIter] != BOUNDARY)
                    {
                        double vData = input.v(i, j + 1, k);
                        double value = vData + invSpacing[1] * (x[flatIndex + jIter] - x[flatIndex]);
                        output.setVDataAt(i, j + 1, k, value);
                    }

                    if (k + 1 < res[2] && markers[flatIndex + kIter] != BOUNDARY)
                    {
                        double wData = input.w(i, j, k + 1);
                        double additionalPressW = invSpacing[2] * (x[flatIndex + kIter] - x[flatIndex]);
                        output.setWDataAt(i, j, k + 1, wData + additionalPressW);
                    }
                }
            }
        }
    }
}

