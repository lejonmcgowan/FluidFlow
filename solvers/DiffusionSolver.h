//
// Created by lejonmcgowan on 10/15/17.
//

#ifndef JET_DIFFUSIONSOLVER_H
#define JET_DIFFUSIONSOLVER_H


#include <memory>
#include <grids/ScalarGrid.h>
#include <grids/MACVectorGrid.h>

#include <eigen3/Eigen/Sparse>
class DiffusionSolver
{
public:
    void solve(const std::shared_ptr<MACVectorGrid> &input, double kDiffusion, double dt, const std::shared_ptr<ScalarGrid> &fluidSdf, const std::shared_ptr<ScalarGrid> &boundarySDF,std::shared_ptr<MACVectorGrid> &output);
private:
    std::vector<Medium> markers;

    void setupMarkers(const MACVectorGrid &input, const ScalarGrid &boundary, const ScalarGrid &fluid);
    void setupSystem(const std::shared_ptr<MACVectorGrid> &input, Vec3 matrix, Size3 dims, int type);
    void applyPressureGradient(const MACVectorGrid &input, MACVectorGrid &output);
    Eigen::SparseMatrix<double> A;
    Eigen::VectorXd x, b;
};


#endif //JET_DIFFUSIONSOLVER_H
