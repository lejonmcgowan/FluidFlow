//
// Created by lejonmcgowan on 10/14/17.
//

#ifndef JET_SEMIGRANGESOLVER_H
#define JET_SEMIGRANGESOLVER_H


#include <grids/MACVectorGrid.h>
#include <memory>
class SemiLagrangianAdvectSolver
{
public:
    SemiLagrangianAdvectSolver();
    virtual ~SemiLagrangianAdvectSolver();

    void advect(MACVectorGrid &input,
                std::shared_ptr<MACVectorGrid> flow,
                double dt,
                MACVectorGrid &output,
                std::shared_ptr<ScalarGrid> colliderSDF);

    void advect(ScalarGrid &input,
                std::shared_ptr<MACVectorGrid> flow,
                double dt,
                ScalarGrid &output,
                std::shared_ptr<ScalarGrid> colliderSDF);
private:
    Vec3 backtrace(std::shared_ptr<MACVectorGrid> flow,
                   double dt,
                   double h,
                   Vec3 start,
                   std::shared_ptr<ScalarGrid> boundarySDF);
};


#endif //JET_SEMIGRANGESOLVER_H
