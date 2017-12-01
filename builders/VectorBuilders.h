//
// Created by lejonmcgowan on 10/10/17.
//

#ifndef JET_VECTORGRIDBUILDER_H
#define JET_VECTORGRIDBUILDER_H

#include <grids/VectorGrid.h>
#include <memory>

class VectorGridBuilder
{
public:
    VectorGridBuilder()
    {}
    virtual ~VectorGridBuilder()
    {}
    virtual std::shared_ptr<VectorGrid> build(Size3 resolution, Vec3 origin, Vec3 spacing, Vec3 initValue) = 0;

};

class VertexVectorBuilder: public VectorGridBuilder
{
public:
    std::shared_ptr<VectorGrid> build(Size3 resolution, Vec3 origin, Vec3 spacing, Vec3 initValue) override;

};

class CenterVectorBuilder: public VectorGridBuilder
{
public:
    std::shared_ptr<VectorGrid> build(Size3 resolution, Vec3 origin, Vec3 spacing, Vec3 initValue) override;

};

class MACGridBuilder: public VectorGridBuilder
{
public:
    std::shared_ptr<VectorGrid> build(Size3 resolution, Vec3 origin, Vec3 spacing, Vec3 initValue) override;

};
#endif //JET_VECTORGRIDBUILDER_H
