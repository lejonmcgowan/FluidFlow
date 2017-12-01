//
// Created by lejonmcgowan on 10/10/17.
//

#include <mathUtils.h>
#include <iomanip>
#include "GridFluidSolver.h"

void GridFluidSolver::printSDFContent(std::string message)
{
    auto velocity = grids().getVelocityGrid();
    const std::shared_ptr<ScalarGrid> &sdf = grids().getScalarGrid("MESH");

    Size3 size = sdf->res();
    std::cout << message << std::endl;

    for (int k = 0; k < size[2]; k++)
    {
        for (int j = 0; j < size[1]; j++)
        {
            for (int i = 0; i < size[0]; i++)
            {

                Vec3 vec = velocity->cellCenterPosition(Size3(i, j, k));
                double sample = sdf->sample(vec);
                
                if(sample != 0.0)
                    std::cout << "(" << i << "," << j << "," << k << "): " << sample << std::endl;
            }
        }
    }

    std::cout << std::endl << std::endl;
}

void GridFluidSolver::printYVelocityContent(std::string message)
{
    auto velocityGrid = std::dynamic_pointer_cast<MACVectorGrid>(solverGrids.getVelocityGrid());
    auto sdf = std::dynamic_pointer_cast<CenterScalarGrid>(solverGrids.getScalarGrid("MESH"));

    Size3 size = velocityGrid->getVDims();
    std::cout << message << std::endl;

    for (int k = 0; k < size[2]; k++)
    {
        for (int j = 0; j < size[1]; j++)
        {
            for (int i = 0; i < size[0]; i++)
            {
                if(velocityGrid->v(i,j,k) != 0.0)
                    std::cout << "(" << i << "," << j << "," << k << "): " << velocityGrid->v(i, j, k) << std::endl;
            }
        }
    }
    std::cout << std::endl << std::endl;
}

void GridFluidSolver::printXVelocityContent(std::string message)
{
    auto velocityGrid = std::dynamic_pointer_cast<MACVectorGrid>(solverGrids.getVelocityGrid());
    auto sdf = std::dynamic_pointer_cast<CenterScalarGrid>(solverGrids.getScalarGrid("MESH"));

    Size3 size = velocityGrid->getUDims();
    std::cout << message << std::endl;

    for (int k = 0; k < size[2]; k++)
    {
        for (int j = 0; j < size[1]; j++)
        {
            for (int i = 0; i < size[0]; i++)
            {
                if(velocityGrid->u(i,j,k) != 0.0)
                    std::cout << "(" << i << "," << j << "," << k << "): " << velocityGrid->u(i, j, k) << std::endl;
            }
        }
    }
    std::cout << std::endl << std::endl;
}

void GridFluidSolver::printZVelocityContent(std::string message)
{
    auto velocityGrid = std::dynamic_pointer_cast<MACVectorGrid>(solverGrids.getVelocityGrid());
    auto sdf = std::dynamic_pointer_cast<CenterScalarGrid>(solverGrids.getScalarGrid("MESH"));

    Size3 size = velocityGrid->getWDims();
    std::cout << message << std::endl;

    for (int k = 0; k < size[2]; k++)
    {
        for (int j = 0; j < size[1]; j++)
        {
            for (int i = 0; i < size[0]; i++)
            {
                if(velocityGrid->w(i,j,k) != 0.0)
                    std::cout << "(" << i << "," << j << "," << k << "): " << velocityGrid->w(i, j, k) << std::endl;
            }
        }
    }
    std::cout << std::endl << std::endl;
}

GridFluidSolver::GridFluidSolver(Size3 res, Vec3 origin, Vec3 spacing)
    :
    solverRes(res), solverOrigin(origin), solverSpacing(spacing), solverGrids(res, origin, spacing)
{
    CenterScalarBuilder meshGridBuilder;
    solverGrids.addScalarGrid("MESH", meshGridBuilder);

    //by default, no collider, so everything is infinitely far away
    colliderSDF = std::make_shared<CenterScalarGrid>(res, origin, spacing, std::numeric_limits<double>().max());
}

void GridFluidSolver::onAct(double delta)
{
    preAct(delta);

    calcExternalForces(delta);
    calcViscocity(delta);
    calcPressure(delta);
    calcAdvection(delta);

    //printYVelocityContent("Velocity after advection");

    postAct(delta);
}

void GridFluidSolver::calcExternalForces(double delta)
{
    applyGravity(delta);
}

void GridFluidSolver::calcAdvection(double delta)
{
    auto velocityGrid = std::dynamic_pointer_cast<MACVectorGrid>(solverGrids.getVelocityGrid());
    auto sdf = std::dynamic_pointer_cast<CenterScalarGrid>(solverGrids.getScalarGrid("MESH"));
    auto sdfInput = std::make_shared<CenterScalarGrid>(*sdf);

    extrapolateVelocityToAir();
    //printYVelocityContent("Velocity after extrapolation to air");

    std::shared_ptr<MACVectorGrid> velocityInput = std::make_shared<MACVectorGrid>(*velocityGrid);

    advectionSolver.advect(*sdfInput, velocityInput, delta, *sdf, colliderSDF);
    advectionSolver.advect(*velocityInput, velocityInput, delta, *velocityGrid, colliderSDF);

    applyBoundaryCondition();

}

void GridFluidSolver::calcViscocity(double delta)
{
    if (calcDiffusion)
    {
        double waterVisc = 1e-4;
        auto velocity = std::dynamic_pointer_cast<MACVectorGrid>(solverGrids.getVelocityGrid());
        auto velInput = std::make_shared<MACVectorGrid>(*velocity);
        auto fluidSDF =  solverGrids.getScalarGrid("MESH");

        diffusionSolver.solve(velInput,waterVisc,delta,fluidSDF,colliderSDF,velocity);
    }
    applyBoundaryCondition();
}

void GridFluidSolver::calcPressure(double delta)
{
    auto velocity = std::dynamic_pointer_cast<MACVectorGrid>(solverGrids.getVelocityGrid());
    auto velInput = std::make_shared<MACVectorGrid>(*velocity);
    auto sdf = std::dynamic_pointer_cast<CenterScalarGrid>(solverGrids.getScalarGrid("MESH"));

    pressureSolver.solve(*velInput, *colliderSDF, *sdf, *velocity);

    applyBoundaryCondition();

}

void GridFluidSolver::preAct(double delta)
{
    assert(grids().spacing()[0] == grids().getVelocityGrid()->spacing()[0]);
}

void GridFluidSolver::postAct(double delta)
{
    reinitializeMesh();
}
void GridFluidSolver::applyGravity(double delta)
{
    auto velocityGrid = std::dynamic_pointer_cast<MACVectorGrid>(solverGrids.getVelocityGrid());
    assert(velocityGrid != nullptr);

    velocityGrid->forEachU([&](long i, long j, long k)
                           {
                               double uValue = velocityGrid->u(i, j, k);
                               velocityGrid->setUDataAt(i, j, k, uValue + delta * gravity[0]);
                           });

    velocityGrid->forEachV([&](long i, long j, long k)
                           {
                               double vValue = velocityGrid->v(i, j, k);
                               velocityGrid->setVDataAt(i, j, k, vValue + delta * gravity[1]);
                           });

    velocityGrid->forEachW([&](long i, long j, long k)
                           {
                               double wValue = velocityGrid->w(i, j, k);
                               velocityGrid->setWDataAt(i, j, k, wValue + delta * gravity[2]);
                           });


    applyBoundaryCondition();

}

void GridFluidSolver::applyBoundaryCondition()
{
    auto velocity = std::dynamic_pointer_cast<MACVectorGrid>(solverGrids.getVelocityGrid());
    assert(velocity != nullptr);

    //for now, I am just satisfying the no flux condition on the domain boundaries

    //A.K.A set velocity to zero on the edges of the grid

    Size3 uDims = velocity->getUDims();
    Size3 vDims = velocity->getVDims();
    Size3 wDims = velocity->getWDims();
    //right and left side
    for (int j = 0; j < uDims[1]; j++)
    {
        for (int k = 0; k < uDims[2]; k++)
        {
            velocity->setUDataAt(0, j, k, 0);
            velocity->setUDataAt(uDims[0] - 1, j, k, 0);
        }
    }
    //top and bottom side
    for (int i = 0; i < vDims[0]; i++)
    {
        for (int k = 0; k < vDims[2]; k++)
        {
            velocity->setVDataAt(i, 0, k, 0);
            velocity->setVDataAt(i, vDims[1] - 1, k, 0);
        }
    }
    //front and back side
    for (int i = 0; i < wDims[0]; i++)
    {
        for (int j = 0; j < wDims[1]; j++)
        {
            velocity->setWDataAt(i, j, 0, 0);
            velocity->setWDataAt(i, j, wDims[2] - 1, 0);
        }
    }
}

Size3 GridFluidSolver::res() const
{ return solverRes; }
Vec3 GridFluidSolver::origin() const
{ return solverOrigin; }
Vec3 GridFluidSolver::spacing() const
{ return solverSpacing; }

GridManager &GridFluidSolver::grids()
{
    return solverGrids;
}

void GridFluidSolver::onRender()
{

}

void GridFluidSolver::reinitializeMesh(const ScalarGrid &inputMesh, double maxDistance, ScalarGrid &outputMesh)
{

    CenterScalarGrid tempGrid(inputMesh.res(), inputMesh.origin(), inputMesh.spacing(), 0);

    inputMesh.foreachIndex([&](long i, long j, long k)
                           {
                               outputMesh.setDataAt(Size3(i, j, k), inputMesh.at(i, j, k));
                               tempGrid.setDataAt(Size3(i, j, k), inputMesh.at(i, j, k));
                           });

    double pseudoT = Utils::pseudoTimeStep(inputMesh, inputMesh.spacing(), 0.5);
    int numIters = static_cast<int>(std::ceil(maxDistance / pseudoT));
    Vec3 spacing = inputMesh.spacing();
    double minSpacing = std::min(spacing[0], std::min(spacing[1], spacing[2]));
    for (int iter = 0; iter < numIters; iter++)
    {
        inputMesh.foreachIndex([&](long i, long j, long k)
                               {
                                   double sdf = outputMesh.at(i, j, k);
                                   double sign = sdf / std::sqrt(sdf * sdf + minSpacing * minSpacing);

                                   std::array<double, 2> minDx, minDy, minDz;

                                   //get derivatives
                                   getDerivatives(outputMesh, Size3(i, j, k), minDx, minDy, minDz);
                                   std::array<double, 2> maxDx = minDx, maxDy = minDy, maxDz = minDz;

                                   minDx[0] = std::max(0.0, minDx[0]);
                                   minDx[1] = std::min(0.0, minDx[1]);
                                   minDy[0] = std::max(0.0, minDy[0]);
                                   minDy[1] = std::min(0.0, minDy[1]);
                                   minDz[0] = std::max(0.0, minDz[0]);
                                   minDz[1] = std::min(0.0, minDz[1]);

                                   minDx[0] *= minDx[0];
                                   minDx[1] *= minDx[1];
                                   minDy[0] *= minDy[0];
                                   minDy[1] *= minDy[1];
                                   minDz[0] *= minDz[0];
                                   minDz[1] *= minDz[1];

                                   maxDx[0] = std::min(0.0, maxDx[0]);
                                   maxDx[1] = std::max(0.0, maxDx[1]);
                                   maxDy[0] = std::min(0.0, maxDy[0]);
                                   maxDy[1] = std::max(0.0, maxDy[1]);
                                   maxDz[0] = std::min(0.0, maxDz[0]);
                                   maxDz[1] = std::max(0.0, maxDz[1]);

                                   maxDx[0] *= maxDx[0];
                                   maxDx[1] *= maxDx[1];
                                   maxDy[0] *= maxDy[0];
                                   maxDy[1] *= maxDy[1];
                                   maxDz[0] *= maxDz[0];
                                   maxDz[1] *= maxDz[1];

                                   double minDistance =
                                       std::sqrt(
                                           minDx[0] + minDx[1] + minDy[0] + minDy[1] + minDz[0] + minDz[1]) - 1;
                                   double maxDistance =
                                       std::sqrt(
                                           maxDx[0] + maxDx[1] + maxDy[0] + maxDy[1] + maxDz[0] + maxDz[1]) - 1;

                                   double value = outputMesh.at(i,j,k)
                                       -pseudoT * std::max(sign, 0.0) * minDistance
                                       - pseudoT * std::min(sign, 0.0) * maxDistance;

                                   tempGrid.setDataAt(Size3(i, j, k), value);

                               });

        //swap grid values
        outputMesh.foreachIndex([&](long i, long j, long k)
                                {
                                    double tempValue = tempGrid.at(i, j, k);
                                    double value = outputMesh.at(i, j, k);

                                    outputMesh.setDataAt(Size3(i, j, k), tempValue);
                                    tempGrid.setDataAt(Size3(i, j, k), value);
                                });
    }
    std::cout << std::endl;

}
void GridFluidSolver::getDerivatives(ScalarGrid &grid,
                                     const Size3 indices,
                                     std::array<double, 2> &dx,
                                     std::array<double, 2> &dy,
                                     std::array<double, 2> &dz)
{
    auto upwindMethod = [](Vec3 values, double d) -> std::array<double,2>
    {
        double invd = 1 / d;
        std::array<double, 2> df{};
        df[0] = invd * (values[1] - values[0]);
        df[1] = invd * (values[2] - values[1]);
        return df;
    };

    long i = indices[0];
    long j = indices[1];
    long k = indices[2];

    Vec3 D0;
    Size3 size = grid.dataSize();
    Vec3 spacing = grid.spacing();

    const long im1 = (i < 1) ? 0 : i - 1;
    const long ip1 = std::min(i + 1, size[0] - 1);
    const long jm1 = (j < 1) ? 0 : j - 1;
    const long jp1 = std::min(j + 1, size[1] - 1);
    const long km1 = (k < 1) ? 0 : k - 1;
    const long kp1 = std::min(k + 1, size[2] - 1);

    D0[0] = grid.at(im1, j, k);
    D0[1] = grid.at(i, j, k);
    D0[2] = grid.at(ip1, j, k);
    dx = upwindMethod(D0, spacing[0]);

    D0[0] = grid.at(i, jm1, k);
    D0[1] = grid.at(i, j, k);
    D0[2] = grid.at(i, jp1, k);
    dy = upwindMethod(D0, spacing[1]);

    D0[0] = grid.at(i, j, km1);
    D0[1] = grid.at(i, j, k);
    D0[2] = grid.at(i, j, kp1);
    dz = upwindMethod(D0, spacing[2]);
}
void GridFluidSolver::reinitializeMesh()
{
    std::shared_ptr<CenterScalarGrid>
        fluid = std::dynamic_pointer_cast<CenterScalarGrid>(grids().getScalarGrid("MESH"));
    assert(fluid != nullptr);

    CenterScalarGrid oldSdf(*fluid);

    double h = std::min(oldSdf.spacing()[0],std::min(oldSdf.spacing()[1],oldSdf.spacing()[2]));
    double minDistance  = 10.0;
    reinitializeMesh(oldSdf, minDistance * h, *fluid);
}
void GridFluidSolver::extrapolateIntoCollider(ScalarGrid &grid)
{
    Size3 dims = grid.dataSize();
    unsigned long size = dims[0] * dims[1] * dims[2];

    //use to quickly determine which parts of the sdf grid is inside, with "true" meaning that the point is inside the sdf (a.k.a distance < 0)
    std::vector<bool> outsideMatrix(size);


    grid.foreachIndex([&](size_t i, size_t j, size_t k)
                      {
                          long flatIndex = k + dims[2] * (j + i * dims[1]);
                          //for now, evething is outside the boundaries, since there are no boundaries to worry about. update here later when you need
                          //actual boundaries
                          outsideMatrix[flatIndex] = true;
                      });

    std::vector<double> flatData;
    grid.getData(&flatData);
    std::vector<double> output(flatData.size());

    unsigned int depth = 1; //w/ CFL of 0.5 ceiling-ed up
    Utils::extrapolateToRegion(flatData, dims, outsideMatrix, depth, output);

    grid.foreachIndex([&](long i, long j, long k)
                      {
                          long flatIndex = k + dims[2] * (j + i * dims[1]);
                          grid.setDataAt(Size3(i, j, k), output[flatIndex]);
                      });
}

void GridFluidSolver::extrapolateVelocityToAir()
{
    double cfl = 5.0;
    auto sdf = grids().getScalarGrid("MESH");
    auto vel = std::dynamic_pointer_cast<MACVectorGrid>(grids().getVelocityGrid());
    assert(vel != nullptr);

    std::vector<bool> insideMatrixU(static_cast<unsigned long>(vel->sizeU()));
    std::vector<bool> insideMatrixV(static_cast<unsigned long>(vel->sizeV()));
    std::vector<bool> insideMatrixW(static_cast<unsigned long>(vel->sizeW()));


    vel->forEachU([&](size_t i, size_t j, size_t k)
                  {
                      const Vec3 &position = vel->dataCenterPosU(Vec3(i, j, k));
                      bool isInsideFluid = sdf->sample(position) < 0;
                      insideMatrixU[vel->flatIndexU(i, j, k)] = isInsideFluid;
                      if (!isInsideFluid)
                          vel->setUDataAt(i, j, k, 0);
                  });

    vel->forEachV([&](size_t i, size_t j, size_t k)
                  {
                      const Vec3 &position = vel->dataCenterPosV(Vec3(i, j, k));
                      bool isInsideFluid = sdf->sample(position) < 0;
                      insideMatrixV[vel->flatIndexV(i, j, k)] = isInsideFluid;
                      if (!isInsideFluid)
                          vel->setVDataAt(i, j, k, 0);
                  });
    vel->forEachW([&](size_t i, size_t j, size_t k)
                  {
                      const Vec3 &position = vel->dataCenterPosW(Vec3(i, j, k));
                      bool isInsideFluid = sdf->sample(position) < 0;
                      insideMatrixW[vel->flatIndexW(i, j, k)] = isInsideFluid;
                      if (!isInsideFluid)
                          vel->setWDataAt(i, j, k, 0);
                  });


    std::vector<double> uData, vData, wData;
    vel->getDataU(uData);
    vel->getDataV(vData);
    vel->getDataW(wData);

    std::vector<double> uOutput(uData.size()), vOutput(vData.size()), wOutput(wData.size());

    Utils::extrapolateToRegion(uData, vel->getUDims(), insideMatrixU, cfl, uOutput);
    Utils::extrapolateToRegion(vData, vel->getVDims(), insideMatrixV, cfl, vOutput);
    Utils::extrapolateToRegion(wData, vel->getWDims(), insideMatrixW, cfl, wOutput);

    vel->forEachU([&](size_t i, size_t j, size_t k)
                  {
                      double value = uOutput[vel->flatIndexU(i, j, k)];
                      vel->setUDataAt(i, j, k, value);
                  });

    vel->forEachV([&](size_t i, size_t j, size_t k)
                  {
                      double value = vOutput[vel->flatIndexV(i, j, k)];
                      vel->setVDataAt(i, j, k, value);
                  });
    vel->forEachW([&](size_t i, size_t j, size_t k)
                  {
                      double value = wOutput[vel->flatIndexW(i, j, k)];
                      vel->setWDataAt(i, j, k, value);
                  });

    applyBoundaryCondition();
}
BBox GridFluidSolver::bounds()
{
    return grids().getVelocityGrid()->getBounds();
}
