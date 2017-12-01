//
// Created by lejonmcgowan on 10/10/17.
//

#include <memory>
#include <numUtils.h>
#include <grids/VertexVectorGrid.h>
#include <grids/CenterVectorGrid.h>
#include <grids/MACVectorGrid.h>
#include "VectorBuilders.h"

std::shared_ptr<VectorGrid> VertexVectorBuilder::build(Size3 resolution, Vec3 origin, Vec3 spacing, Vec3 initValue)
{
    std::shared_ptr<VertexVectorGrid> grid = std::make_shared<VertexVectorGrid>();
    grid->setSize(resolution, origin, spacing, initValue);
    return grid;
}

std::shared_ptr<VectorGrid> CenterVectorBuilder::build(Size3 resolution, Vec3 origin, Vec3 spacing, Vec3 initValue)
{
    std::shared_ptr<CenterVectorGrid> grid = std::make_shared<CenterVectorGrid>();
    grid->setSize(resolution, origin, spacing, initValue);
    return grid;
}

std::shared_ptr<VectorGrid> MACGridBuilder::build(Size3 resolution, Vec3 origin, Vec3 spacing, Vec3 initValue)
{
    std::shared_ptr<MACVectorGrid> grid = std::make_shared<MACVectorGrid>();
    grid->setSize(resolution, origin, spacing, initValue);
    return grid;
}
