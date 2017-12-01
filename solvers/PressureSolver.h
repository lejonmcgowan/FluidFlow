//
// Created by lejonmcgowan on 10/23/17.
//

#ifndef JET_PRESSURESOLVER_H
#define JET_PRESSURESOLVER_H
#include "../grids/MACVectorGrid.h"
#include "../grids/ScalarGrid.h"

#include <eigen3/Eigen/Sparse>
class PressureSolver
{
public:
    PressureSolver();
    virtual ~PressureSolver();
    void solve(const MACVectorGrid &input,
               const ScalarGrid &boundarySDF,
               const ScalarGrid &fluidSDF,
               MACVectorGrid &output);

private:
    std::vector<Medium> markers;

    void setupMarkers(const MACVectorGrid &input, const ScalarGrid &boundary, const ScalarGrid &fluid);
    void setupSystem(const MACVectorGrid &input);
    void applyPressureGradient(const MACVectorGrid &input, MACVectorGrid &output);
    Eigen::SparseMatrix<double> A;
    Eigen::VectorXd x, b;
};
#endif //JET_PRESSURESOLVER_H
