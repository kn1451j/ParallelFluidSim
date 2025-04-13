#include <string.h>
#include "grid.hpp"


/*
Perform a weighted particle to grid update for the velocities
*/
void Grid::transfer_to_grid(std::vector<Particle>& particles)
{
    // interpolate particles based on grid location -> weighted average of particles based on distance from vertex (normalized)

    this->horizontal_velocity.zero();
    // first we want to compute for each cell the count of particles in range for normalization
    

    // then we want to linearly interpolate velocity and apply a particle mass weighting

    // then we set horizontal and vertical velocities based off of this for each grid vertex
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

void Grid::print_grid()
{
    for(size_t row_idx = 0; row_idx<ROW_NUM; row_idx++)
    {
        for(size_t col_idx = 0; col_idx<COL_NUM; col_idx++)
        {
            char cell_str = this->cells[row_idx][col_idx].type ==FLUID ? 'F' : cells[row_idx][col_idx].type == GAS ? 'G' : 'S';
            printf("[%zu, %zu] (%c):\n", row_idx, col_idx, cell_str);
            printf("    p_ij: %f\n", pressure_grid[row_idx][col_idx].value);
            printf("    u_i0, j: %f; u_i1, j: %f\n", horizontal_velocity[row_idx][2*col_idx].value, horizontal_velocity[row_idx][2*col_idx+1].value);
            printf("    u_i0, j: %f; u_i1, j: %f\n", vertical_velocity[2*row_idx][col_idx].value, vertical_velocity[2*row_idx+1][col_idx].value);
        }
    }
}