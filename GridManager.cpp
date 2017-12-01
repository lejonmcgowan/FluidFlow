//
// Created by lejonmcgowan on 10/10/17.
//

#include "GridManager.h"

GridManager::GridManager(const Size3 resolution, const Vec3 gridOrigin, const Vec3 gridSpacing)
    :
    resolution(resolution),
    gridSpacing(gridSpacing),
    gridOrigin(gridOrigin)
{
    vectorGrids["VELOCITY"] = std::make_shared<MACVectorGrid>();
    auto velocityGrid = std::dynamic_pointer_cast<MACVectorGrid>(vectorGrids["VELOCITY"]);

    velocityGrid->setSize(resolution, gridOrigin, gridSpacing);
    gridBounds = velocityGrid->getBounds();
}

void GridManager::addScalarGrid(std::string name, ScalarGridBuilder &builder, double initValue)
{
    std::transform(name.begin(), name.end(), name.begin(), ::toupper);
    scalarGrids[name] = builder.build(resolution, origin(), spacing(), initValue);
}

void GridManager::addVectorGrid(std::string name, VectorGridBuilder &builder, Vec3 initValue)
{
    std::transform(name.begin(), name.end(), name.begin(), ::toupper);
    //no, you no overwrite velocity grid
    if (name != "VELOCITY")
        vectorGrids[name] = builder.build(resolution, origin(), spacing(), initValue);
}
std::shared_ptr<VectorGrid> GridManager::getVelocityGrid()
{
    return vectorGrids["VELOCITY"];
}
std::shared_ptr<VectorGrid> GridManager::getVectorGrid(std::string name)
{
    return vectorGrids[name];
}
std::shared_ptr<ScalarGrid> GridManager::getScalarGrid(std::string name)
{
    return scalarGrids[name];
}
