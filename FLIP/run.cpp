#include "fluid.hpp"

// entry point for the simulation
int main()
{
    Fluid fluid(500,500);
    fluid.run();

    return 0;
}