#include <string.h>
#include <cassert>
#include "grid.hpp"


Grid::Grid(double width, double height){
    this->width = width;
    this->height = height;
    this->cell_width = width/COL_NUM;
    this->cell_height = height/ROW_NUM;
    this->cell_volume = this->cell_width * this->cell_height;

    // initialize the "poses" for each vertex
    for(size_t row_idx = 0; row_idx<ROW_NUM; row_idx++)
    {
        for(size_t col_idx = 0; col_idx<COL_NUM; col_idx++)
        {
            // set the center position of each cell
            this->cells[row_idx][col_idx].position = Point((col_idx+0.5)*this->cell_width, (row_idx+0.5)*this->cell_height, 0.0);
            // printf(this->cells[row_idx][col_idx].position.print().c_str());

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

// returns -> [row_idx, col_idx]
grid_idx_t Grid::get_grid_idx(Particle p)
{
    return std::pair<int, int>(
        static_cast<int>(p.position.y/this->cell_height),
        static_cast<int>(p.position.x/this->cell_width)
    );
}

// returns the 4 horizontal cell vertex indices bounding the particle
Neighbors Grid::get_cell_neighbors(Particle p)
{
    int row_idx = static_cast<int>(p.position.y/this->cell_height - 0.5);
    int col_idx = static_cast<int>(p.position.x/this->cell_width - 0.5);

    Neighbors n {};
    n.type = CELL;

    // we dont check the value for negativity here, since negative implies solid boundary -> info we might need
    n.neighbors[0] = std::pair<int, int>(row_idx, col_idx);
    n.neighbors[1] = std::pair<int, int>(row_idx+1, col_idx);
    n.neighbors[2] = std::pair<int, int>(row_idx, col_idx+1);
    n.neighbors[3] = std::pair<int, int>(row_idx+1, col_idx+1);

    return n;
}

// returns the 4 horizontal velocity vertex indices bounding the particle
Neighbors Grid::get_horizontal_neighbors(Particle p)
{
    int row_idx = static_cast<int>(p.position.y/this->cell_height - 0.5);
    int col_idx = static_cast<int>(p.position.x/this->cell_width);

    Neighbors n {};
    n.type = HORIZONTAL;

    // we dont check the value for negativity here, since negative implies solid boundary -> info we might need
    n.neighbors[0] = std::pair<int, int>(row_idx, col_idx);
    n.neighbors[1] = std::pair<int, int>(row_idx+1, col_idx);
    n.neighbors[2] = std::pair<int, int>(row_idx, col_idx+1);
    n.neighbors[3] = std::pair<int, int>(row_idx+1, col_idx+1);

    return n;
}

Neighbors Grid::get_vertical_neighbors(Particle p)
{
    int row_idx = static_cast<int>(p.position.y/this->cell_height);
    int col_idx = static_cast<int>(p.position.x/this->cell_width - 0.5);

    Neighbors n {};
    n.type = VERTICAL;

    // we dont check the value for negativity here, since negative implies solid boundary -> info we might need
    n.neighbors[0] = std::pair<int, int>(row_idx, col_idx);
    n.neighbors[1] = std::pair<int, int>(row_idx+1, col_idx);
    n.neighbors[2] = std::pair<int, int>(row_idx, col_idx+1);
    n.neighbors[3] = std::pair<int, int>(row_idx+1, col_idx+1);

    return n;
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
                
            this->cells[row_idx][col_idx].reset();
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
        // printf("particle pos: %s\n", p.print().c_str());

        auto cell_idx = this->get_grid_idx(p);
        Cell cell = this->cells[cell_idx.first][cell_idx.second];
        // printf("cell (%d, %d)\n", cell_idx.first, cell_idx.second);

        // at least one particle in cell -> the cell becomes fluid
        cell.type = FLUID;
        // DENSITY CALCULATED LATER

        // get the neighboring horizontal velocity cells
        Neighbors nh = this->get_horizontal_neighbors(p);
        for(grid_idx_t cell_idx : nh.neighbors){
            // printf("cell (%d, %d)\n", cell_idx.first, cell_idx.second);
            // check if neighbor is border. if so, do nothing -> this will be a solid cell
            if(cell_idx.first < 0) continue;

            #ifdef DEBUG
            assert(cell_idx.second >= 0);
            #endif

            Vertex& ui0 = this->horizontal_velocity[cell_idx.first][cell_idx.second];

            // then we set horizontal and vertical velocities based off of this for each grid vertex
            ui0.value += p.mass * p.velocity.x * p.position.bilinear(ui0.position, this->cell_width, this->cell_height); // TODO -> should this be bilinear? just 4 closest neighbors?
            // ui0.normalization += p.position.bilinear(ui0.position, this->cell_width, this->cell_height);
            ui0.normalization += p.mass;
        }

        // do the same for vertical
        Neighbors nv = this->get_vertical_neighbors(p);
        for(grid_idx_t cell_idx : nv.neighbors){
            // check if neighbor is border. if so, do nothing -> this will be a solid cell
            if(cell_idx.second < 0) continue;

            #ifdef DEBUG
            assert(cell_idx.first >= 0);
            #endif

            Vertex& vi0 = this->vertical_velocity[cell_idx.first][cell_idx.second];
            vi0.value += p.mass * p.velocity.y * p.position.bilinear(vi0.position, this->cell_width, this->cell_height);

            // printf("dist: %f", p.position.ngp_distance(vi0.position, this->cell_width, this->cell_height));

            // printf("vi.value: %f\n", vi0.value);
            // vi0.normalization += p.position.bilinear(vi0.position, this->cell_width, this->cell_height);
            vi0.normalization += p.mass;
        }

        // calculate density by interpolating
        Neighbors n = this->get_cell_neighbors(p);
        for(grid_idx_t cell_idx : n.neighbors){
            // check if neighbor is border. if so, do nothing -> this will be a solid cell (no density)
            if(cell_idx.second < 0 || cell_idx.first < 0) continue;

            Cell& v = this->cells[cell_idx.first][cell_idx.second];
            v.density += (p.mass/this->cell_volume) * p.position.bilinear(v.position, this->cell_width, this->cell_height);
            // v.normalization += p.position.bilinear(v.position, this->cell_width, this->cell_height);
        }
    }

    // for each non-zero value cell (at least one associated particle), we normalize the velocity weighted by distance and
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
    for(Particle& p : particles)
    {
        p.velocity = Point();

        // get the vertical and horizontal neighbors for each particle, bilinearly interpolate their velocity onto the particle
        Neighbors nh = this->get_horizontal_neighbors(p);
        for(grid_idx_t cell_idx : nh.neighbors){
            Vertex& u = this->horizontal_velocity[cell_idx.first][cell_idx.second];
            p.velocity.x += PIC_WEIGHT * p.position.bilinear(u.position, this->cell_width, this->cell_height) * u.value; // PIC update
            p.velocity.x += (1 - PIC_WEIGHT) * p.position.bilinear(u.position, this->cell_width, this->cell_height) * (u.value - u.prev_value); // FLIP update
        }

        Neighbors nv = this->get_vertical_neighbors(p);
        for(grid_idx_t cell_idx : nv.neighbors){
            Vertex& v = this->vertical_velocity[cell_idx.first][cell_idx.second];
            p.velocity.y += PIC_WEIGHT * p.position.bilinear(v.position, this->cell_width, this->cell_height) * v.value; // PIC update
            p.velocity.y += (1 - PIC_WEIGHT) * p.position.bilinear(v.position, this->cell_width, this->cell_height) * (v.value - v.prev_value); // FLIP update
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
            printf("    density: %f\n", cells[row_idx][col_idx].density);
        }
    }
}