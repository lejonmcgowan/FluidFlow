//
// Created by lejonmcgowan on 10/10/17.
//

#include <iostream>
#include "FluidSolver.h"

const double FluidSolver::DEFAULT_TIMESTEP = 1 / 60.0;

void FluidSolver::update(int frameNumber)
{
    int numFrames = 0;
    if (frameNumber > frameIdx)
    {
        numFrames = frameNumber - frameIdx;
        for (int i = frameIdx; i < frameIdx + numFrames; i++)
        {
            act(timestep);
            std::cout << "Frame " << (i + 1) << " integrated" << std::endl;
        }
        frameIdx = frameNumber;
    }

}

void FluidSolver::act(double delta)
{
    //for now, let's just do 1 act per frame. we can do more substeps and get more precise answers later
    const int numSubSteps = 1;
    double subDelta = delta / (double) numSubSteps;

    for (int i = 0; i < numSubSteps; i++)
        onAct(subDelta);

}

void FluidSolver::render()
{

}
