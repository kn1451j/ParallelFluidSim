#include "grid.hpp"

/*
Perform a weighted particle to grid update for the velocities
*/
void Grid::transfer_to_grid(std::vector<Particle>& particles)
{
    return;
}

/*
Perform a FLIP fluid update
*/
void Grid::transfer_from_grid(std::vector<Particle>& particles)
{
    // for(size_t idx = 0; idx<particles.size(); idx++)
    // {
    //     particles[idx] = Particle();
    // }
}

void Grid::solve_pressure()
{
    return;
}