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

grid_idx_t Grid::get_particle_idx(Particle p)
{
    return std::tuple<size_t, size_t>(p.position.x/this->cell_width, p.position.y/this->cell_height);
}