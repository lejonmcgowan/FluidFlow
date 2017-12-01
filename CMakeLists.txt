#
# Copyright (c) 2017 Doyub Kim
#
# I am making my contributions/submissions to this project solely in my personal
# capacity and am not conveying any rights to any intellectual property of any
# third parties.
#

# Target name
set(target sandbox)

# Includes
include_directories(${CMAKE_CURRENT_SOURCE_DIR})

# Sources
file(GLOB my_sources
        grids/CenterScalarGrid.h grids/VertextScalarGrid.h grids/CollocatedGrid.h grids/CenterVectorGrid.h grids/VertexVectorGrid.h
        grids/MACVectorGrid.h MyScalarBuilder.h MyVectorBuilder.h grids/ScalarGrid.cpp mathUtils.h grids/CollatedGrid.cpp
        grids/MACVectorGrid.cpp mathUtils.cpp
        MyGrid.h BBox.h grids/ScalarGrid.h numUtils.h grids/VectorGrid.h
        GridManager.cpp GridManager.h builders/ScalarBuilders.h
        builders/ScalarBuilders.cpp builders/VectorBuilders.cpp
        builders/VectorBuilders.h solvers/GridFluidSolver.h
        solvers/GridFluidSolver.cpp solvers/FluidSolver.cpp solvers/FluidSolver.h Collider.h
        Surface.h ImplicitSurface.h
        solvers/PressureSolver.h solvers/PressureSolver.cpp
        solvers/SemiLagrangianAdvectSolver.cpp solvers/SemiLagrangianAdvectSolver.h
        solvers/DiffusionSolver.cpp solvers/DiffusionSolver.h)

# Build executable
add_executable(${target} ${my_sources} main.cpp jetAdapter.h )


add_executable(${target}_test ${my_sources} test/GridCalcTests.cpp test/main.cpp)

#change requirements for test to not care about unused variables. Generally need these in test files for logging purposes
set_target_properties(${target}_test PROPERTIES COMPILE_FLAGS "-Wno-unused-variable -Wno-unused-but-set-variable")
set_target_properties(${target} PROPERTIES COMPILE_FLAGS "-Wno-unused-variable -Wno-unused-but-set-variable")

# Project options
set_target_properties(${target}
        PROPERTIES
        ${DEFAULT_PROJECT_OPTIONS})

set_target_properties(${target}_test
        PROPERTIES
        ${DEFAULT_PROJECT_OPTIONS})

# mostly so I don't have to worry about uninitialized variables (use those for debugging). I'll re-enable it later
# Compile options
target_compile_options(${target}
        PRIVATE
        PUBLIC
        ${DEFAULT_COMPILE_OPTIONS}
        INTERFACE)

target_compile_options(${target}_test
        PRIVATE
        PUBLIC
        ${DEFAULT_COMPILE_OPTIONS}
        INTERFACE)

# Link libraries
target_link_libraries(${target}
        PRIVATE
        ${DEFAULT_LINKER_OPTIONS}
        jet)

target_link_libraries(${target}_test
        PRIVATE
        ${DEFAULT_LINKER_OPTIONS}
        gtest)
