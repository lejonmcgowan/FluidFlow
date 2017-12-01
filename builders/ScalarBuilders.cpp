//
// Created by lejonmcgowan on 10/10/17.
//

#include <memory>
#include <numUtils.h>
#include <grids/VertextScalarGrid.h>
#include <grids/CenterScalarGrid.h>
#include "ScalarBuilders.h"

std::shared_ptr<ScalarGrid> VertexScalarBuilder::build(Size3 resolution, Vec3 origin, Vec3 spacing, double initValue)
{
    std::shared_ptr<VertexScalarGrid> grid = std::make_shared<VertexScalarGrid>();
    grid->setSize(resolution, origin, spacing);
    grid->setInitValue(initValue);

    return grid;
}

std::shared_ptr<ScalarGrid> CenterScalarBuilder::build(Size3 resolution, Vec3 origin, Vec3 spacing, double initValue)
{
    std::shared_ptr<CenterScalarGrid> grid = std::make_shared<CenterScalarGrid>();
    grid->setSize(resolution, origin, spacing);
    grid->setInitValue(initValue);

    return grid;
}
