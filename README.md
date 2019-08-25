# FluidFlow

![example fluid simulated from the code](https://github.com/lejonmcgowan/FluidFlow/blob/master/fluinAnimFast.gif?raw=true)

![example fluid in slow motion](https://github.com/lejonmcgowan/FluidFlow/blob/master/fluid_large_sim.gif?raw=true)

Senior Project Fluid solver. Eulerian based, 3D, C++. A simulator made based on the research from Doyub Kim's *fluid engine development*: https://fluidenginedevelopment.org/


Work on this involved the implementation of various grid data strutcutres to store 3d spatial data of fluids into, and then creating various solvers to figure out the velocity, advection, and pressure of the grid every frame. I then took my prepared data and bridged into Kim's implementation of the Marching Cubes algortihm in the author's repository: https://github.com/doyubkim/fluid-engine-dev

in oder to be able to render my results into a visual frame (I intend to create my own implementation of Marching Cubes in the future so the program can render on its own without any dependencies)


