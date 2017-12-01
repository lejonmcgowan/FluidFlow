//
// Created by lejonmcgowan on 10/10/17.
//

#ifndef JET_FLUIDSOLVER_H
#define JET_FLUIDSOLVER_H


class FluidSolver
{
public:
    const static double DEFAULT_TIMESTEP;

    FluidSolver()
    {}
    FluidSolver(double timestep)
        : timestep(timestep)
    {}
    virtual ~FluidSolver()
    {}
    void update(int frameNumber);

protected:
    virtual void onAct(double delta) = 0;
    virtual void onRender() = 0;
private:
    int frameIdx = -1;

    void act(double delta);
    void render();

    const double timestep = DEFAULT_TIMESTEP;
};


#endif //JET_FLUIDSOLVER_H
