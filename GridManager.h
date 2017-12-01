//
// Created by lejonmcgowan on 10/10/17.
//

#ifndef JET_GRIDMANAGER_H
#define JET_GRIDMANAGER_H


#include <memory>
#include <map>
#include <builders/ScalarBuilders.h>
#include <builders/VectorBuilders.h>
#include "grids/MACVectorGrid.h"
#include "grids/ScalarGrid.h"

class GridManager
{
public:
    GridManager(Size3 resolution, Vec3 gridSpacing, Vec3 gridOrigin);
    virtual ~GridManager()
    {
        std::cout << "deleting gridmanager" << std::endl;
    }

    Size3 res() const
    { return resolution; }
    Vec3 spacing() const
    { return gridSpacing; }
    Vec3 origin() const
    { return gridOrigin; }
    BBox bounds() const
    { return gridBounds; }

    void addScalarGrid(std::string name, ScalarGridBuilder &builder, double initValue = 0);
    void addVectorGrid(std::string name, VectorGridBuilder &builder, Vec3 initValue = Vec3(0, 0, 0));
    std::shared_ptr<VectorGrid> getVectorGrid(std::string name);
    std::shared_ptr<ScalarGrid> getScalarGrid(std::string name);
    std::shared_ptr<VectorGrid> getVelocityGrid();

private:
    std::map<std::string, std::shared_ptr<ScalarGrid> > scalarGrids;
    std::map<std::string, std::shared_ptr<VectorGrid> > vectorGrids;

    Size3 resolution;
    Vec3 gridSpacing;
    Vec3 gridOrigin;
    BBox gridBounds;
};


#endif //JET_GRIDMANAGER_H
