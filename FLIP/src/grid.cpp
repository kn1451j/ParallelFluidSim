#include <string.h>
#include "grid.hpp"


Grid::Grid(double width, double height){
    this->width = width;
    this->height = height;
    this->cell_width = width/COL_NUM;
    this->cell_height = height/ROW_NUM;

    // initialize the "poses" for each vertex
    for(size_t row_idx = 0; row_idx<ROW_NUM; row_idx++)
    {
        for(size_t col_idx = 0; col_idx<COL_NUM; col_idx++)
        {
            Vertex& vi0 = this->vertical_velocity[row_idx][col_idx];
            vi0.position = Point(this->cell_width*(col_idx+0.5), this->cell_height*row_idx, 0.0);

            Vertex& ui0 = this->horizontal_velocity[row_idx][col_idx];
            ui0.position = Point(this->cell_width*col_idx, this->cell_height*(row_idx+0.5), 0.0);
        }

        Vertex& ui0 = this->horizontal_velocity[row_idx][COL_NUM];
        ui0.position =  Point(this->cell_width*COL_NUM, this->cell_height*(row_idx+0.5), 0.0);
    }

    for(size_t col_idx = 0; col_idx<COL_NUM; col_idx++)
    { 
        Vertex& ui0 = this->vertical_velocity[ROW_NUM][col_idx];
        ui0.position = Point(this->cell_width*(col_idx+0.5), this->cell_height*ROW_NUM, 0.0);
    }

    this->A = new SparseMatrix(ROW_NUM, COL_NUM);
};

grid_idx_t Grid::get_grid_idx(Particle p)
{
    return std::pair<size_t, size_t>(
        static_cast<size_t>(p.position.x/this->cell_width), 
        static_cast<size_t>(p.position.y/this->cell_height)
    );
}

// returns the 4 horizontal velocity vertex indices bounding the particle
Neighbors get_horizontal_neighbors(Particle p)
{
    // TODO
}

Neighbors get_vertical_neighbors(Particle p)
{
    // TODO
}

/*
Perform a weighted particle to grid update for the velocities
*/
void Grid::transfer_to_grid(std::vector<Particle>& particles)
{
    // reset cells
    for(size_t row_idx = 0; row_idx<ROW_NUM; row_idx++)
    {
        for(size_t col_idx = 0; col_idx<COL_NUM; col_idx++)
        {
            // check if bounds - if so, solid. o.w. gas
            // if(row_idx == ROW_NUM - 1 || col_idx == 0 || col_idx == COL_NUM - 1)
            //     this->cells[row_idx][col_idx].type = SOLID;
            // else
                
            this->cells[row_idx][col_idx].type = GAS;
            this->pressure_grid[row_idx][col_idx].reset();
            this->horizontal_velocity[row_idx][col_idx].reset();
            this->vertical_velocity[row_idx][col_idx].reset();
        }

        // update horizontal velocity for every column
        this->horizontal_velocity[row_idx][COL_NUM].reset();
    }

    for(size_t col_idx = 0; col_idx<COL_NUM; col_idx++)
    { 
        this->vertical_velocity[ROW_NUM][col_idx].reset();
    }


    // interpolate particles based on grid location -> weighted average of particles based on distance from vertex (normalized)
    for(Particle p : particles){
        // get the neighboring cells
        auto cell_idx = this->get_grid_idx(p);

        this->cells[cell_idx.first][cell_idx.second].type = FLUID;

        Neighbors nh = this->get_horizontal_neighbors(p);
        for(grid_idx_t cell_idx : nh.neighbors){
            Vertex& ui0 = this->horizontal_velocity[cell_idx.first][cell_idx.second];

            // then we want to linearly interpolate velocity and apply a particle mass weighting
            // then we set horizontal and vertical velocities based off of this for each grid vertex
            ui0.value += p.velocity.x * p.position.l1_distance(ui0.position); // TODO -> should this be bilinear? just 4 closest neighbors?
            ui0.normalization += p.position.l1_distance(ui0.position);
        }

        // do the same for vertical
        Neighbors nv = this->get_vertical_neighbors(p);
        for(grid_idx_t cell_idx : nv.neighbors){
            Vertex& vi0 = this->vertical_velocity[cell_idx.first][cell_idx.second];
            vi0.value += p.velocity.y * p.position.l1_distance(vi0.position);
            vi0.normalization += p.position.l1_distance(vi0.position);
        }
    }

    // for each non-zero value cell, we normalize
    for(size_t row_idx = 0; row_idx<ROW_NUM; row_idx++)
    {
        for(size_t col_idx = 0; col_idx<COL_NUM; col_idx++)
        {
            Vertex& ui0 = this->horizontal_velocity[row_idx][col_idx];
            Vertex& vi0 = this->vertical_velocity[row_idx][col_idx];

            if(ui0.normalization>EPS) ui0.value /= ui0.normalization;
            if(vi0.normalization>EPS) vi0.value /= vi0.normalization;
        }

        // update horizontal velocity for every column
        Vertex& ui0 = this->horizontal_velocity[row_idx][COL_NUM];
        if(ui0.normalization>EPS) ui0.value /= ui0.normalization;
    }

    for(size_t col_idx = 0; col_idx<COL_NUM; col_idx++)
    { 
        Vertex& vi0 = this->vertical_velocity[ROW_NUM][col_idx];
        if(vi0.normalization>EPS) vi0.value /= vi0.normalization;
    }
}

/*
Perform a weighted PIC-FLIP fluid update by bilinearly interpolating grid vertices onto particles
*/
void Grid::transfer_from_grid(std::vector<Particle>& particles)
{
    double cell_size = this->cell_height * this->cell_width; //used as the denominator in bilinear interpolation

    for(Particle& p : particles)
    {
        p.velocity = Point();

        // get the vertical and horizontal neighbors for each particle, bilinearly interpolate their velocity onto the particle
        Neighbors nh = this->get_horizontal_neighbors(p);
        for(grid_idx_t cell_idx : nh.neighbors){
            Vertex& u = this->horizontal_velocity[cell_idx.first][cell_idx.second];
            p.velocity.x += PIC_WEIGHT * p.position.bilinear(u.position) * u.value / cell_size; // PIC update
            p.velocity.x += (1 - PIC_WEIGHT) * p.position.bilinear(u.position) * (u.value - u.prev_value) / cell_size; // FLIP update
        }

        Neighbors nv = this->get_vertical_neighbors(p);
        for(grid_idx_t cell_idx : nv.neighbors){
            Vertex& v = this->vertical_velocity[cell_idx.first][cell_idx.second];
            p.velocity.y += PIC_WEIGHT * p.position.bilinear(v.position) * v.value / cell_size; // PIC update
            p.velocity.y += (1 - PIC_WEIGHT) * p.position.bilinear(v.position) * (v.value - v.prev_value) / cell_size; // FLIP update
        }
    }
}

void Grid::solve_pressure()
{
    return;
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
            printf("    u, j: %f; v, j: %f\n", horizontal_velocity[row_idx][col_idx].value, vertical_velocity[row_idx][col_idx].value);
        }
    }
}