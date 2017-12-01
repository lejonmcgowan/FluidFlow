//
// Created by lejonmcgowan on 10/15/17.
//

#include "DiffusionSolver.h"
void DiffusionSolver::solve(const std::shared_ptr<MACVectorGrid> &input,
                            double kDiffusion,
                            double dt,
                            const std::shared_ptr<ScalarGrid> &fluidSdf,
                            const std::shared_ptr<ScalarGrid> &boundarySDF,
                            std::shared_ptr<MACVectorGrid> &output)
{
    Vec3 h = input->spacing();
    Vec3 c = Vec3(1 / h[0] / h[0], 1 / h[1] / h[1], 1 / h[2] / h[2]);
    c *= dt * kDiffusion;
    Eigen::ConjugateGradient<Eigen::SparseMatrix<double>,Eigen::Upper > solver;

    //solve u
    setupMarkers(*input, *boundarySDF, *fluidSdf);
    setupSystem(input, c, input->getUDims(),0);


    x = solver.compute(A).solve(b);

    Size3 size = input->getUDims();
    for (int i = 0; i < size[0]; i++)
    {
        for (int j = 0; j < size[1]; j++)
        {
            for (int k = 0; k < size[2]; k++)
            {
                double value = input->u(i, j, k);
                output->setUDataAt(Size3(i, j, k), x[input->flatIndexU(i,j,k)]);
            }
        }
    }

    //solve v
    setupMarkers(*input, *boundarySDF, *fluidSdf);
    setupSystem(input, c, input->getVDims(),1);


    x = solver.compute(A).solve(b);

    size = input->getVDims();
    for (int i = 0; i < size[0]; i++)
    {
        for (int j = 0; j < size[1]; j++)
        {
            for (int k = 0; k < size[2]; k++)
            {
                double value = input->u(i, j, k);
                output->setUDataAt(Size3(i, j, k), x[input->flatIndexV(i,j,k)]);
            }
        }
    }

    //solve w
    setupMarkers(*input, *boundarySDF, *fluidSdf);
    setupSystem(input, c, input->getWDims(),2);


    x = solver.compute(A).solve(b);

    size = input->getWDims();
    for (int i = 0; i < size[0]; i++)
    {
        for (int j = 0; j < size[1]; j++)
        {
            for (int k = 0; k < size[2]; k++)
            {
                double value = input->u(i, j, k);
                output->setUDataAt(Size3(i, j, k), x[input->flatIndexW(i,j,k)]);
            }
        }
    }
}
void DiffusionSolver::setupSystem(const std::shared_ptr<MACVectorGrid> &input, Vec3 c, Size3 dims, int type)
{
    typedef Eigen::Triplet<double> DataPoint;

    long size = dims[0] * dims[1] * dims[2];
    A = Eigen::SparseMatrix<double>(size, size);
    std::vector<DataPoint> Adata;
    b = Eigen::VectorXd::Zero(size);
    x.resize(size);

    //build matrix along diagonal
    for (int diagIndex = 0; diagIndex < A.rows(); diagIndex++)
    {
        int tempI = diagIndex;
        //unflatten indices
        long unflatI = tempI / (dims[1] * dims[2]);
        tempI -= unflatI * dims[1] * dims[2];
        long unflatJ = tempI / dims[2];
        long unflatK = tempI % dims[2];

        Size3 gridIndices(unflatI, unflatJ, unflatK);
        long kIter = 1, jIter = dims[2], iIter = dims[2] * dims[1];


        if (markers[diagIndex] == FLUID)
        {

            b[diagIndex] = 1.0;

            //in x dimension
            if (unflatI + 1 < dims[0])
            {
                Adata.emplace_back(diagIndex, diagIndex, c[0]);
                Adata.emplace_back(diagIndex, diagIndex + iIter, -c[0]);
            }

            if (unflatI > 0)
            {
                Adata.emplace_back(diagIndex, diagIndex, c[0]);
            }
            //in y dimension
            if (unflatJ + 1 < dims[1])
            {
                Adata.emplace_back(diagIndex, diagIndex, c[1]);
                Adata.emplace_back(diagIndex, diagIndex + jIter, -c[1]);
            }

            if (unflatJ > 0)
            {
                Adata.emplace_back(diagIndex, diagIndex, c[1]);
            }
            //in z dimension
            if (unflatK + 1 < dims[2])
            {
                Adata.emplace_back(diagIndex, diagIndex, c[2]);
                Adata.emplace_back(diagIndex, diagIndex + kIter, -c[2]);
            }

            if (unflatK > 0)
            {
                Adata.emplace_back(diagIndex, diagIndex, c[2]);

            }

        }
        else
        {
            Adata.emplace_back(diagIndex, diagIndex, 1);
        }
    }

    A.setFromTriplets(Adata.begin(), Adata.end());

    //setup solution vector
    assert(type >= 0 && type <= 2);
    for (int i = 0; i < size; i++)
    {
        int tempI = i;
        //unflatten indices
        long unflatI = tempI / (dims[1] * dims[2]);
        tempI -= unflatI * dims[1] * dims[2];
        long unflatJ = tempI / dims[2];
        long unflatK = tempI % dims[2];

        if(type == 0)
            b[i] = x[i] = input->u(unflatI, unflatJ, unflatK);
        else if(type == 1)
            b[i] = x[i] = input->v(unflatI, unflatJ, unflatK);
        else
            b[i] = x[i] = input->w(unflatI, unflatJ, unflatK);
    }
}

void DiffusionSolver::setupMarkers(const MACVectorGrid &input, const ScalarGrid &boundary, const ScalarGrid &fluid)
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
}
void DiffusionSolver::applyPressureGradient(const MACVectorGrid &input, MACVectorGrid &output)
{

}
