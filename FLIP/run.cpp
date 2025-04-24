#include "fluid.hpp"

// entry point for the simulation
int main()
{
    cv::viz::Viz3d window("3D Plot");
    Fluid fluid(500,500,10, window);
    fluid.run();

    return 0;
}