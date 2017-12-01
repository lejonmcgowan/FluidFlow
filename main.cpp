// Copyright (c) 2017 Doyub Kim
//
// I am making my contributions/submissions to this project solely in my
// personal capacity and am not conveying any rights to any intellectual
// property of any third parties.

#include <solvers/GridFluidSolver.h>
#include "grids/CenterScalarGrid.h"
#include "jetAdapter.h"
#include "mathUtils.h"
#include <jet/sphere3.h>
#include <jet/jet.h>
#include <iomanip>
#include <set>

void compareData(CenterScalarGrid& grid, MACVectorGrid& velocity, jet::ImplicitSurface3& surface);
void printMeshInfo(jet::TriangleMesh3 &mesh, int frame);
void testMatrix();

int zmain()
{
    using namespace jet;
    Vector3D pt(0.27499973385562959, 0.1744999669591201, 0.024999999999948189);

    jet::FaceCenteredGrid3 grid(jet::Size3(15,15,15),Vector3D(0.05,0.05,0.05));
    auto bounds = grid.boundingBox();
    Vector3D maxBound = 0.5 * bounds.upperCorner;
    const double radius = 0.35 * grid.resolution().x * grid.gridSpacing().x;
    maxBound[1] = radius;
    jet::Vector3D start(maxBound[0], maxBound[1], maxBound[2]);

    grid.fill([&](Vector3D vec) -> Vector3D
              {
                  double i = vec.x;
                  double j = vec.y;
                  double k = vec.z;
                  jet::Vector3D pos(i, j, k);
                  return pos;
              });

    auto sampler = grid.sampler();

    Vector3D s = sampler(pt);
    Vec3 sample(s.x,s.y,s.z);

    auto uSourceSampler = CubicArraySampler3<double, double>(
        grid.uConstAccessor(),
        grid.gridSpacing(),
        grid.uOrigin());
    auto vSourceSampler = CubicArraySampler3<double, double>(
        grid.vConstAccessor(),
        grid.gridSpacing(),
        grid.vOrigin());
    auto wSourceSampler = CubicArraySampler3<double, double>(
        grid.wConstAccessor(),
        grid.gridSpacing(),
        grid.wOrigin());
    auto cubicSampler = [uSourceSampler, vSourceSampler, wSourceSampler](const Vector3D& x) -> Vector3D
    {
            return Vector3D(
                uSourceSampler(x), vSourceSampler(x), wSourceSampler(x));
    };

    Vector3D c = cubicSampler(pt);
    Vec3 cubicSample(c.x,c.y,c.z);


    MACVectorGrid myGrid;
    myGrid.setSize(::Size3(15,15,15),Vec3(0,0,0),Vec3(0.05,0.05,0.05),Vec3(0,0,0));
    myGrid.fillData([&](double i, double j, double k) -> Vec3
                       {
                           Vec3 pos(i, j, k);
                           return pos;
                       });



    Vec3 mySample = myGrid.sample(Vec3(pt.x,pt.y,pt.z));
    Vec3 myCubeSample = myGrid.cubicSample(Vec3(pt.x,pt.y,pt.z));

    std::cout << "sample    (linear): " << sample.transpose()             << std::endl;
    std::cout << "my result (linear): " << mySample.transpose()    << std::endl << std::endl;
    std::cout << "sample    (cubic): " << cubicSample.transpose()             << std::endl;
    std::cout << "my result (cubic): " << myCubeSample.transpose() << std::endl;

    return 0;
}

int main()
{
    //testMatrix();
    int numFrames = 300;
    Vec3 spacing(0.05, 0.05, 0.05);
    GridFluidSolver solver(Size3(30,35,30), Vec3(0, 0, 0), spacing);

    auto meshGrid = std::dynamic_pointer_cast<CenterScalarGrid>(solver.grids().getScalarGrid("MESH"));
    auto velocity = std::dynamic_pointer_cast<MACVectorGrid>(solver.grids().getVelocityGrid());
    assert(meshGrid != nullptr);

    BBox bounds = solver.bounds();
    Vec3 maxBound = 0.5 * bounds.getMax();
    const double radius = 0.35 * solver.res()[0] * spacing[0];
    jet::Vector3D start(maxBound[0], maxBound[1], maxBound[2]);
    //start[1] = radius;
    auto sphere = std::make_shared<jet::Sphere3>(jet::Sphere3::Builder().withRadius(radius).withCenter(start).build());
    auto implicitSphere = jet::SurfaceToImplicit3(sphere);
    meshGrid->fillData([&](double i, double j, double k) -> double
                       {
                           jet::Vector3D pos(i, j, k);
                           return implicitSphere.signedDistance(pos);
                       });

    //render first; frame 0 should be the iniial mesh
    jet::TriangleMesh3 finalMesh = makeMeshFromGrid(*meshGrid);

    std::stringstream filename;
    filename << "frame_" << std::setfill('0') << std::setw(6) << 0 << ".obj";
    std::cout << "rendering frame 0" << std::endl;
    finalMesh.writeObj(filename.str());
    std::vector<double> origData;
    std::vector<double> frameMeshData;
    meshGrid->getData(&origData);
    printMeshInfo(finalMesh, 0);
    compareData(*meshGrid,*velocity,implicitSphere);
    for (int i = 1; i <= numFrames; i++)
    {
        if(i < 0)
            std::cout.setstate(std::ios_base::failbit);

        std::cout << "rendering frame " << i << std::endl;
        solver.update(i - 1);
        std::cout.clear();

        meshGrid->getData(&frameMeshData);
        meshGrid->getData(&origData);

        finalMesh = makeMeshFromGrid(*meshGrid);
        printMeshInfo(finalMesh, i);

        std::stringstream newFileName;
        newFileName << "frame_" << std::setfill('0') << std::setw(6) << i << ".obj";
        finalMesh.writeObj(newFileName.str());
    }

    return 0;
}

void compareData(CenterScalarGrid& grid, MACVectorGrid& velocity, jet::ImplicitSurface3& surface)
{
    auto oracleBuilder = jet::CellCenteredScalarGrid3::Builder();

    auto oracleGrid = oracleBuilder.withResolution(jet::Size3(grid.res()[0], grid.res()[1], grid.res()[2]))
                                                .withOrigin(grid.origin()[0],grid.origin()[1],grid.origin()[2])
                                                .withGridSpacing(grid.spacing()[0],grid.spacing()[1],grid.spacing()[2])
                                                .build();

    oracleGrid.fill([&](jet::Vector3D pos)
    {
        return surface.signedDistance(pos);
    });

    oracleGrid.forEachDataPointIndex([&](size_t i, size_t j, size_t k)
                                {
                                    double diff = std::fabs(grid.at(i,j,k) - oracleGrid(i, j, k));
                                    if(diff > 0.01)
                                        std::cout << "difference of " << diff << "at index " << Size3(i,j,k).transpose() << std::endl;
                                });
}

void printMeshInfo(jet::TriangleMesh3 &mesh, int frame)
{
    std::cout << std::endl << "Mesh " << frame << std::endl;
    std::cout << "-----------------------" << std::endl;
    std::cout << "num points  : " << mesh.numberOfPoints() << std::endl;
    std::cout << "num tris    : " << mesh.numberOfTriangles() << std::endl;
    std::cout << "num normals : " << mesh.numberOfNormals() << std::endl;
    std::cout << "num UVs     : " << mesh.numberOfUvs() << std::endl << std::endl;
}

void testMatrix()
{
    const int DIM = 8;
    Eigen::MatrixXd mat1D = Eigen::MatrixXd::Zero(DIM, DIM);

    for (int i = 0; i < DIM; i++)
    {
        if (i > 0)
        {
            mat1D(i, i)++;
            mat1D(i - 1, i)--;
        }
        if (i + 1 < DIM)
        {
            mat1D(i, i)++;
            mat1D(i + 1, i)--;
        }

    }

    std::cout << mat1D << std::endl;
}
