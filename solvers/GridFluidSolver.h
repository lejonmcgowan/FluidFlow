//
// Created by lejonmcgowan on 10/10/17.
//

#ifndef JET_GRIDFLUIDSOLVER_H
#define JET_GRIDFLUIDSOLVER_H

#include <GridManager.h>
#include <grids/CenterScalarGrid.h>
#include "FluidSolver.h"
#include "SemiLagrangianAdvectSolver.h"
#include "DiffusionSolver.h"
#include "PressureSolver.h"

class GridFluidSolver: public FluidSolver
{
public:
    GridFluidSolver(Size3 res, Vec3 origin = Vec3(0, 0, 0), Vec3 spacing = Vec3(1, 1, 1));
    virtual ~GridFluidSolver()
    {}

    Size3 res() const;
    Vec3 spacing() const;
    Vec3 origin() const;
    BBox bounds();
protected:
    void onRender() override;
public:
    GridManager &grids();
protected:
    void onAct(double delta) override;

    virtual void calcExternalForces(double delta);
    virtual void calcViscocity(double delta);
    virtual void calcPressure(double delta);
    virtual void calcAdvection(double delta);
private:
    Size3 solverRes;
    Vec3 solverOrigin;
    Vec3 solverSpacing;
    Vec3 gravity = Vec3(0, -9.8, 0);
    GridManager solverGrids;
    std::shared_ptr<CenterScalarGrid> colliderSDF;
    int count = 0;

    bool calcDiffusion;
    SemiLagrangianAdvectSolver advectionSolver;
    DiffusionSolver diffusionSolver;
    PressureSolver pressureSolver;

    void applyGravity(double delta);
    void applyBoundaryCondition();
    void preAct(double delta);
    void postAct(double delta);

    void reinitializeMesh();
    void extrapolateIntoCollider(ScalarGrid &grid);
    void reinitializeMesh(const ScalarGrid &inputMesh, double maxDistance, ScalarGrid &outputMesh);
    void getDerivatives(ScalarGrid &grid,
                        const Size3 indices,
                        std::array<double, 2> &dx,
                        std::array<double, 2> &dy,
                        std::array<double, 2> &dz);

    void extrapolateVelocityToAir();



    void printSDFContent(std::string message);
    void printYVelocityContent(std::string message);
    void printZVelocityContent(std::string message);
    void printXVelocityContent(std::string message);
};

#endif //JET_GRIDFLUIDSOLVER_H
