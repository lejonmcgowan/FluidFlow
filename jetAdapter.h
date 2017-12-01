//
// Created by lejonmcgowan on 10/19/17.
//

#ifndef JET_JETADAPTER_H
#define JET_JETADAPTER_H


#include <jet/triangle_mesh3.h>
#include <jet/marching_cubes.h>
#include <grids/ScalarGrid.h>

jet::TriangleMesh3 makeMeshFromGrid(ScalarGrid &grid)
{
    jet::TriangleMesh3 result;
    jet::Size3 res(grid.dataSize()[0], grid.dataSize()[1], grid.dataSize()[2]);
    jet::Vector3D spacing(grid.spacing()[0], grid.spacing()[1], grid.spacing()[2]);
    jet::Vector3D origin(grid.origin()[0], grid.origin()[1], grid.origin()[2]);

    //make jet's accessor from my grid
    std::vector<double> data;
    grid.getData(&data);
    jet::ConstArrayAccessor3<double> jetArray(res, data.data());

    int flag = jet::kDirectionAll & ~jet::kDirectionDown;
    jet::marchingCubes(jetArray, spacing, origin, &result, 0.0, flag);

    return result;
}
#endif //JET_JETADAPTER_H
