//
// Created by lejonmcgowan on 10/14/17.
//

#include <grids/ScalarGrid.h>
#include "SemiLagrangianAdvectSolver.h"
#include <cmath>

SemiLagrangianAdvectSolver::SemiLagrangianAdvectSolver()
{}

SemiLagrangianAdvectSolver::~SemiLagrangianAdvectSolver()
{

}
Vec3 SemiLagrangianAdvectSolver::backtrace(std::shared_ptr<MACVectorGrid> flow,
                                     double dt,
                                     double h,
                                     Vec3 start,
                                     std::shared_ptr<ScalarGrid> boundarySDF)
{
    Vec3 pt0 = start;
    Vec3 pt1 = start;

    double remaining = dt;

    //because floating point comparison
    while (remaining > std::numeric_limits<double>::epsilon())
    {
        //get midpoint
        Vec3 vel0 = flow->sample(pt0);
        double numSubSteps = std::max(std::ceil(vel0.norm() * remaining / h), 1.0);
        dt = remaining / numSubSteps;

        Vec3 midPoint = pt0 - 0.5 * dt * vel0;
        Vec3 midVelocity = flow->sample(midPoint);
        pt1 = pt0 - dt * midVelocity;

        //boundary handle
        double phi0 = boundarySDF->sample(pt0);
        double phi1 = boundarySDF->sample(pt1);

        if (phi0 * phi1 < 0.0)
        {
            double w = std::fabs(phi1) / (std::fabs(phi0) + std::fabs(phi1));
            pt1 = w * pt0 + (1.0 - w) * pt1;
            break;
        }

        remaining -= dt;
        pt0 = pt1;
    }

    return pt1;
}

void SemiLagrangianAdvectSolver::advect(MACVectorGrid &input,
                                  std::shared_ptr<MACVectorGrid> flow,
                                  double dt,
                                  MACVectorGrid &output,
                                  std::shared_ptr<ScalarGrid> colliderSDF)
{
    double h = std::min(input.spacing()[0], std::min(input.spacing()[1], input.spacing()[2]));
    input.forEachU([&](long i, long j, long k)
                   {
                       Vec3 inPos = input.dataCenterPosU(i, j, k);
                       Vec3 outPos = output.dataCenterPosU(i, j, k);
                       if (colliderSDF->sample(inPos) > 0.0)
                       {
                           Vec3 outPoint = backtrace(flow, dt, h, outPos, colliderSDF);
                           Vec3 sample = input.cubicSample(outPoint);
                           output.setUDataAt(i, j, k, sample[0]);
                       }
                   });

    input.forEachV([&](long i, long j, long k)
                   {
                       Vec3 inPos = input.dataCenterPosV(i, j, k);
                       Vec3 outPos = output.dataCenterPosV(i, j, k);
                       if (colliderSDF->sample(inPos) > 0.0)
                       {
                           Vec3 outPoint = backtrace(flow, dt, h, outPos, colliderSDF);
                           Vec3 sampleOut = input.cubicSample(outPoint);
                           output.setVDataAt(i, j, k, sampleOut[1]);
                       }
                   });

    input.forEachW([&](long i, long j, long k)
                   {
                       Vec3 inPos = input.dataCenterPosW(i, j, k);
                       Vec3 outPos = output.dataCenterPosW(i, j, k);
                       if (colliderSDF->sample(inPos) > 0.0)
                       {
                           Vec3 outPoint = backtrace(flow, dt, h, outPos, colliderSDF);
                           output.setWDataAt(i, j, k, input.cubicSample(outPoint)[2]);
                       }
                   });

}
void SemiLagrangianAdvectSolver::advect(ScalarGrid &input,
                                  std::shared_ptr<MACVectorGrid> flow,
                                  double dt,
                                  ScalarGrid &output,
                                  std::shared_ptr<ScalarGrid> colliderSDF)
{
    double h = std::min(output.spacing()[0], std::min(output.spacing()[1], output.spacing()[2]));

    Size3 size = input.dataSize();
    for (int k = 0; k < size[2]; k++)
    {
        for (int j = 0; j < size[1]; j++)
        {
            for (int i = 0; i < size[0]; i++)
            {

                Vec3 inPos = input.dataCenterPosition(i, j, k);
                Vec3 outPos = output.dataCenterPosition(i, j, k);
                double collideSample = colliderSDF->sample(inPos);
                //take into account NaN, not in a boundary if it does this
                collideSample = std::isnan(collideSample) ? 1 : collideSample;

                if (collideSample > 0.0)
                {
                    Vec3 outPoint = backtrace(flow, dt, h, outPos, colliderSDF);
                    double outPointSample = input.cubicSample(outPoint);
                    output.setDataAt(Size3(i, j, k), outPointSample);
                }

            }
        }
    }

}
