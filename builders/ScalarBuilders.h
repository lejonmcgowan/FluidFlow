//
// Created by lejonmcgowan on 10/10/17.
//

#ifndef JET_SCALARGRIDBUILDER_H
#define JET_SCALARGRIDBUILDER_H

#include <grids/ScalarGrid.h>
#include <memory>

class ScalarGridBuilder
{
public:
    ScalarGridBuilder()
    {}
    virtual ~ScalarGridBuilder()
    {}
    virtual std::shared_ptr<ScalarGrid> build(Size3 resolution, Vec3 origin, Vec3 spacing, double initValue) = 0;

};

class VertexScalarBuilder: public ScalarGridBuilder
{
public:
    std::shared_ptr<ScalarGrid> build(Size3 resolution, Vec3 origin, Vec3 spacing, double initValue) override;

};

class CenterScalarBuilder: public ScalarGridBuilder
{
public:

    std::shared_ptr<ScalarGrid> build(Size3 resolution, Vec3 origin, Vec3 spacing, double initValue) override;

};
#endif //JET_SCALARGRIDBUILDER_H
