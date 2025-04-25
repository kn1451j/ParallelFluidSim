#include <string.h>
#include <cassert>

#include "grid.hpp"
// #include "profiler.hpp"

Grid::Grid(double width, double height, double depth, double density, Profiler *p = nullptr) {
    // initialize profiler
    this->profiler = p;

    this->width = width;
    this->height = height;
    this->depth = depth;
    this->cell_width = width / COL_NUM;
    this->cell_height = height / ROW_NUM;
    this->cell_depth = depth / DEPTH_NUM;
    this->cell_volume = this->cell_width * this->cell_height * this->cell_depth;
    this->density = density;

    // initialize the "poses" for each vertex
    for(size_t depth_idx = 0; depth_idx<DEPTH_NUM; depth_idx++)
    {
        for(size_t row_idx = 0; row_idx<ROW_NUM; row_idx++)
        {
            for(size_t col_idx = 0; col_idx<COL_NUM; col_idx++)
            {
                // set the center position of each cell
                double x = (col_idx+0.5)*this->cell_width;
                double y = (row_idx+0.5)*this->cell_height;
                double z = (depth_idx+0.5)*this->cell_depth;
                this->cells[row_idx][col_idx][depth_idx].position = Point(x, y, z);
                // printf(this->cells[row_idx][col_idx].position.print().c_str());

                Vertex& vi0 = this->vertical_velocity[row_idx][col_idx][depth_idx];
                double vx = this->cell_width*(col_idx+0.5);
                double vy = this->cell_height*row_idx;
                double vz = this->cell_depth*(depth_idx+0.5);
                vi0.position = Point(vx, vy, vz);

                Vertex& ui0 = this->horizontal_velocity[row_idx][col_idx][depth_idx];
                double ux = this->cell_width*col_idx;
                double uy = this->cell_height*(row_idx+0.5);
                double uz = this->cell_depth*(depth_idx+0.5);
                ui0.position = Point(ux, uy, uz);

                Vertex& zi0 = this->depth_velocity[row_idx][col_idx][depth_idx];
                double zx = this->cell_width*(col_idx+0.5);
                double zy = this->cell_height*(row_idx+0.5);
                double zz = this->cell_depth*depth_idx;
                zi0.position = Point(zx, zy, zz);
            }

            Vertex& ui0 = this->horizontal_velocity[row_idx][COL_NUM][depth_idx];
            double ux = this->cell_width*COL_NUM;
            double uy = this->cell_height*(row_idx+0.5);
            double uz = this->cell_depth*(depth_idx+0.5);
            ui0.position =  Point(ux, uy, uz);
        }

        for(size_t col_idx = 0; col_idx<COL_NUM; col_idx++)
        {
            Vertex& vi0 = this->vertical_velocity[ROW_NUM][col_idx][depth_idx];
            double vx = this->cell_width*(col_idx+0.5);
            double vy = this->cell_height*ROW_NUM;
            double vz = this->cell_depth*(depth_idx+0.5);
            vi0.position = Point(vx, vy, vz);
        }
    }
    
    for(size_t row_idx = 0; row_idx<ROW_NUM; row_idx++)
    {
        for(size_t col_idx = 0; col_idx<COL_NUM; col_idx++)
        {
            Vertex& zi0 = this->depth_velocity[row_idx][col_idx][DEPTH_NUM];
            double ux = this->cell_width*(col_idx+0.5);
            double uy = this->cell_height*(row_idx+0.5);
            double uz = this->cell_depth*DEPTH_NUM;
            zi0.position = Point(ux, uy, uz);
        }
    }

    // initialize pressure solver
    this->dV.resize(ROW_NUM * COL_NUM * DEPTH_NUM);
    this->sparseA.resize(ROW_NUM * COL_NUM * DEPTH_NUM, std::vector<double>(SPARSE_WIDTH));
    this->pVec.resize(ROW_NUM * COL_NUM * DEPTH_NUM);
    this->diagE.resize(ROW_NUM * COL_NUM * DEPTH_NUM);
};

// returns -> [row_idx, col_idx]
grid_idx_t Grid::get_grid_idx(Particle p)
{
    return {
        static_cast<int>(p.position.y/this->cell_height),
        static_cast<int>(p.position.x/this->cell_width),
        static_cast<int>(p.position.z/this->cell_depth),
    };
}

// returns the four cell vertex indices bounding the particle
Neighbors Grid::get_cell_neighbors(Particle p)
{
    grid_idx_t cell_idx = this->get_grid_idx(p);

    Neighbors n {};
    n.type = CELL;

    // we dont check the value for negativity here, since negative implies solid boundary -> info we might need
    n.neighbors[0] = {cell_idx.row - 1, cell_idx.col, cell_idx.depth};
    n.neighbors[1] = {cell_idx.row, cell_idx.col - 1, cell_idx.depth};
    n.neighbors[2] = {cell_idx.row, cell_idx.col + 1, cell_idx.depth};
    n.neighbors[3] = {cell_idx.row + 1, cell_idx.col, cell_idx.depth};
    n.neighbors[4] = {cell_idx.row, cell_idx.col, cell_idx.depth - 1};
    n.neighbors[5] = {cell_idx.row, cell_idx.col, cell_idx.depth + 1};

    return n;
}

// returns the 4 horizontal velocity vertex indices bounding the particle
Neighbors Grid::get_horizontal_neighbors(Particle p)
{
    int row_idx = static_cast<int>(p.position.y/this->cell_height - 0.5);
    int col_idx = static_cast<int>(p.position.x/this->cell_width);
    int depth_idx = static_cast<int>(p.position.z/this->cell_depth - 0.5);

    Neighbors n {};
    n.type = HORIZONTAL;

    // we dont check the value for negativity here, since negative implies solid boundary -> info we might need
    n.neighbors[0] = {row_idx, col_idx, depth_idx};
    n.neighbors[1] = {row_idx + 1, col_idx, depth_idx};
    n.neighbors[2] = {row_idx, col_idx + 1, depth_idx};
    n.neighbors[3] = {row_idx + 1, col_idx + 1, depth_idx};
    n.neighbors[4] = {row_idx, col_idx, depth_idx + 1};
    n.neighbors[5] = {row_idx + 1, col_idx, depth_idx + 1};
    n.neighbors[6] = {row_idx, col_idx + 1, depth_idx + 1};
    n.neighbors[7] = {row_idx + 1, col_idx + 1, depth_idx + 1};

    return n;
}

Neighbors Grid::get_vertical_neighbors(Particle p)
{
    int row_idx = static_cast<int>(p.position.y/this->cell_height);
    int col_idx = static_cast<int>(p.position.x/this->cell_width - 0.5);
    int depth_idx = static_cast<int>(p.position.z/this->cell_depth - 0.5);

    Neighbors n {};
    n.type = VERTICAL;

    // we dont check the value for negativity here, since negative implies solid boundary -> info we might need
    n.neighbors[0] = {row_idx, col_idx, depth_idx};
    n.neighbors[1] = {row_idx + 1, col_idx, depth_idx};
    n.neighbors[2] = {row_idx, col_idx + 1, depth_idx};
    n.neighbors[3] = {row_idx + 1, col_idx + 1, depth_idx};
    n.neighbors[4] = {row_idx, col_idx, depth_idx + 1};
    n.neighbors[5] = {row_idx + 1, col_idx, depth_idx + 1};
    n.neighbors[6] = {row_idx, col_idx + 1, depth_idx + 1};
    n.neighbors[7] = {row_idx + 1, col_idx + 1, depth_idx + 1};

    return n;
}

Neighbors Grid::get_depth_neighbors(Particle p)
{
    int row_idx = static_cast<int>(p.position.y/this->cell_height - 0.5);
    int col_idx = static_cast<int>(p.position.x/this->cell_width - 0.5);
    int depth_idx = static_cast<int>(p.position.z/this->cell_depth);

    Neighbors n {};
    n.type = DEPTH;

    // we dont check the value for negativity here, since negative implies solid boundary -> info we might need
    n.neighbors[0] = {row_idx, col_idx, depth_idx};
    n.neighbors[1] = {row_idx + 1, col_idx, depth_idx};
    n.neighbors[2] = {row_idx, col_idx + 1, depth_idx};
    n.neighbors[3] = {row_idx + 1, col_idx + 1, depth_idx};
    n.neighbors[4] = {row_idx, col_idx, depth_idx + 1};
    n.neighbors[5] = {row_idx + 1, col_idx, depth_idx + 1};
    n.neighbors[6] = {row_idx, col_idx + 1, depth_idx + 1};
    n.neighbors[7] = {row_idx + 1, col_idx + 1, depth_idx + 1};

    return n;
}

/*
Perform a weighted particle to grid update for the velocities
*/
void Grid::transfer_to_grid(std::vector<Particle>& particles)
{
    // reset cells
    for (size_t depth_idx = 0; depth_idx<DEPTH_NUM; depth_idx++){
        for(size_t row_idx = 0; row_idx<ROW_NUM; row_idx++)
        {
            for(size_t col_idx = 0; col_idx<COL_NUM; col_idx++)
            {  
                this->cells[row_idx][col_idx][depth_idx].reset();
                this->pressure_grid[row_idx][col_idx][depth_idx].reset();
                this->horizontal_velocity[row_idx][col_idx][depth_idx].reset();
                this->vertical_velocity[row_idx][col_idx][depth_idx].reset();
                this->depth_velocity[row_idx][col_idx][depth_idx].reset();
            }

            // update horizontal velocity for every column
            this->horizontal_velocity[row_idx][COL_NUM][depth_idx].reset();
        }
        
        for(size_t col_idx = 0; col_idx<COL_NUM; col_idx++)
        { 
            this->vertical_velocity[ROW_NUM][col_idx][depth_idx].reset();
        }
    }
    
    for(size_t row_idx = 0; row_idx<ROW_NUM; row_idx++)
    {
        for(size_t col_idx = 0; col_idx<COL_NUM; col_idx++)
        {
            this->depth_velocity[row_idx][col_idx][DEPTH_NUM].reset();
        }
    }

    // interpolate particles based on grid location -> weighted average of particles based on distance from vertex (normalized)
    for(Particle p : particles){
        // at least one particle in cell -> the cell becomes fluid
        grid_idx_t cell_idx = this->get_grid_idx(p);
        Cell& cell = this->cells[cell_idx.row][cell_idx.col][cell_idx.depth];
        cell.type = FLUID;
        double dist = p.position.bilinear(cell.position, this->cell_width, this->cell_height, this->cell_depth);
        cell.density += p.mass * dist;

        // get the neighboring horizontal velocity cells
        Neighbors nh = this->get_horizontal_neighbors(p);
        for(grid_idx_t cell_idx : nh.neighbors){
            // check if neighbor is border. if so, do nothing -> this will be a solid cell
            if(!this->_valid_hcell(cell_idx)) continue;

            Vertex& ui0 = this->horizontal_velocity[cell_idx.row][cell_idx.col][cell_idx.depth];

            // then we set horizontal and vertical velocities based off of this for each grid vertex
            double dist = p.position.bilinear(ui0.position, this->cell_width, this->cell_height, this->cell_depth);
            ui0.value += p.mass * dist * p.velocity.x; // TODO -> should this be bilinear? just 4 closest neighbors?
            ui0.normalization += dist;
        } 

        // do the same for vertical
        Neighbors nv = this->get_vertical_neighbors(p);
        for(grid_idx_t cell_idx : nv.neighbors){
            // check if neighbor is border. if so, do nothing -> this will be a solid cell
            if(!this->_valid_vcell(cell_idx)) continue;

            Vertex& vi0 = this->vertical_velocity[cell_idx.row][cell_idx.col][cell_idx.depth];

            double dist =  p.position.bilinear(vi0.position, this->cell_width, this->cell_height, this->cell_depth);
            vi0.value += p.mass * dist * p.velocity.y;
            vi0.normalization += dist;
        }

        Neighbors nz = this->get_depth_neighbors(p);
        for(grid_idx_t cell_idx : nz.neighbors){
            // check if neighbor is border. if so, do nothing -> this will be a solid cell
            if(!this->_valid_dcell(cell_idx)) continue;

            Vertex& di0 = this->depth_velocity[cell_idx.row][cell_idx.col][cell_idx.depth];

            // then we set horizontal and vertical velocities based off of this for each grid vertex
            double dist = p.position.bilinear(di0.position, this->cell_width, this->cell_height, this->cell_depth);
            di0.value += p.mass * dist * p.velocity.z; // TODO -> should this be bilinear? just 4 closest neighbors?
            di0.normalization += dist;
        } 

        // calculate density by interpolating
        Neighbors n = this->get_cell_neighbors(p);
        for(int idx = 0; idx < NUM_NEIGHBORS; idx++){
            grid_idx_t cell_idx = n.neighbors[idx];
            // check if neighbor is border. if so, do nothing -> this will be a solid cell (no density)
            if(!this->_valid_cell(cell_idx)) continue;

            Cell& cell = this->cells[cell_idx.row][cell_idx.col][cell_idx.depth];
            cell.type = FLUID;

            double dist = p.position.bilinear(cell.position, this->cell_width, this->cell_height, this->cell_depth);
            cell.density += (p.mass / this->cell_volume) * dist;
        }
    }

    for (size_t depth_idx = 0; depth_idx<DEPTH_NUM; depth_idx++){
        for(size_t row_idx = 0; row_idx<ROW_NUM; row_idx++)
        {
            for(size_t col_idx = 0; col_idx<COL_NUM; col_idx++)
            {  
                Vertex& ui0 =this->horizontal_velocity[row_idx][col_idx][depth_idx];
                Vertex& vi0 =this->vertical_velocity[row_idx][col_idx][depth_idx];
                Vertex& zi0 =this->depth_velocity[row_idx][col_idx][depth_idx];

                if(ui0.normalization>EPS) ui0.value /= ui0.normalization;
                if(vi0.normalization>EPS) vi0.value /= vi0.normalization;
                if(zi0.normalization>EPS) zi0.value /= zi0.normalization;
            }

            // update horizontal velocity for every column
            Vertex& ui0 = this->horizontal_velocity[row_idx][COL_NUM][depth_idx];
            if(ui0.normalization>EPS) ui0.value /= ui0.normalization;
        }
        
        for(size_t col_idx = 0; col_idx<COL_NUM; col_idx++)
        { 
            Vertex& vi0 = this->vertical_velocity[ROW_NUM][col_idx][depth_idx];
            if(vi0.normalization>EPS) vi0.value /= vi0.normalization;
        }
    }
    
    for(size_t row_idx = 0; row_idx<ROW_NUM; row_idx++)
    {
        for(size_t col_idx = 0; col_idx<COL_NUM; col_idx++)
        {
            Vertex& zi0 = this->depth_velocity[row_idx][col_idx][DEPTH_NUM];
            if(zi0.normalization>EPS) zi0.value /= zi0.normalization;
        }
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
        double velocity_norm = 0.0;
        for(grid_idx_t cell_idx : nh.neighbors){
            if(!this->_valid_hcell(cell_idx)) continue;

            Vertex& u = this->horizontal_velocity[cell_idx.row][cell_idx.col][cell_idx.depth];
            double dist = p.position.bilinear(u.position, this->cell_width, this->cell_height, this->cell_depth);
            p.velocity.x += PIC_WEIGHT * dist * u.value; // PIC update
            p.velocity.x += (1 - PIC_WEIGHT) * dist * (u.value - u.prev_value); // FLIP update
            velocity_norm += dist;
        }

        // normalize
        if(velocity_norm>0)
            p.velocity.x /= velocity_norm;

        Neighbors nv = this->get_vertical_neighbors(p);
        velocity_norm = 0.0;
        for(grid_idx_t cell_idx : nv.neighbors){
            if(!this->_valid_vcell(cell_idx)) continue;

            Vertex& v = this->vertical_velocity[cell_idx.row][cell_idx.col][cell_idx.depth];
            double dist = p.position.bilinear(v.position, this->cell_width, this->cell_height, this->cell_depth);
            p.velocity.y += PIC_WEIGHT * dist * v.value; // PIC update
            p.velocity.y += (1 - PIC_WEIGHT) * dist * (v.value - v.prev_value); // FLIP update
            velocity_norm += dist;
        }

        if(velocity_norm>0)
            p.velocity.y /= velocity_norm;

        Neighbors nz = this->get_depth_neighbors(p);
        velocity_norm = 0.0;
        for(grid_idx_t cell_idx : nz.neighbors){
            if(!this->_valid_dcell(cell_idx)) continue;

            Vertex& v = this->depth_velocity[cell_idx.row][cell_idx.col][cell_idx.depth];
            double dist = p.position.bilinear(v.position, this->cell_width, this->cell_height, this->cell_depth);
            p.velocity.z += PIC_WEIGHT * dist * v.value; // PIC update
            p.velocity.z += (1 - PIC_WEIGHT) * dist * (v.value - v.prev_value); // FLIP update
            velocity_norm += dist;
        }

        if(velocity_norm>0)
            p.velocity.z /= velocity_norm;
    }
}

/*
    Solves d * V = dt (d * dP) / density
        for every pressure cell i,j
    Then updates velocities s.t. d * V = 0
    
    Uses conjugate gradient method to solve Ap = dV
        where A is a first-order linear approx of dt * d / density
*/
void Grid::solve_pressure(double dt)
{
    // reset the pressure solver each time
    this->reset();

    // populate right hand-side divergence of current velocity field
    double dx = this->cell_width;
    double dy = this->cell_height;
    double dz = this->cell_depth;
    for (int depth_idx = 0; depth_idx < DEPTH_NUM; depth_idx ++)
    {    
        for (int row_idx = 0; row_idx < ROW_NUM; row_idx ++)
        {
            for (int col_idx = 0; col_idx < COL_NUM; col_idx ++)
            {
                if (this->_full_fluid_cell({row_idx, col_idx, depth_idx})){
                    // Calculate divergence correctly
                    double dU = (this->horizontal_velocity[row_idx][col_idx+1][depth_idx].value
                            - this->horizontal_velocity[row_idx][col_idx][depth_idx].value) / dx;
                    double dV = (this->vertical_velocity[row_idx+1][col_idx][depth_idx].value
                            - this->vertical_velocity[row_idx][col_idx][depth_idx].value) / dy;
                    double dZ = (this->depth_velocity[row_idx][col_idx][depth_idx + 1].value
                            - this->depth_velocity[row_idx][col_idx][depth_idx].value) / dz;

                    this->set_dV_idx({row_idx, col_idx, depth_idx}, -(dU + dV + dZ));
                
                    // Handle solid boundaries
                    if (this->_solid_hborder_left({row_idx, col_idx, depth_idx})){
                        double hv = this->horizontal_velocity[row_idx][col_idx][depth_idx].value;
                        this->add_to_dV({row_idx, col_idx, depth_idx}, - hv / dx);
                    }
                    if (this->_solid_hborder_right({row_idx, col_idx, depth_idx})){
                        double hv = this->horizontal_velocity[row_idx][col_idx + 1][depth_idx].value;
                        this->add_to_dV({row_idx, col_idx, depth_idx}, hv / dx);
                    }
                    if (this->_solid_vborder_top({row_idx, col_idx, depth_idx})){
                        double hv = this->vertical_velocity[row_idx + 1][col_idx][depth_idx].value;
                        this->add_to_dV({row_idx, col_idx, depth_idx}, hv / dy);
                    }
                    if (this->_solid_vborder_bot({row_idx, col_idx, depth_idx})){
                        double hv = this->vertical_velocity[row_idx][col_idx][depth_idx].value;
                        this->add_to_dV({row_idx, col_idx, depth_idx}, - hv / dy);
                    } 
                    if (this->_solid_dborder_front({row_idx, col_idx, depth_idx})){
                        double hv = this->depth_velocity[row_idx][col_idx][depth_idx].value;
                        this->add_to_dV({row_idx, col_idx, depth_idx}, - hv / dz);
                    }
                    if (this->_solid_dborder_back({row_idx, col_idx, depth_idx})){
                        double hv = this->depth_velocity[row_idx][col_idx][depth_idx + 1].value;
                        this->add_to_dV({row_idx, col_idx, depth_idx}, hv / dz);
                    } 
                }
            }
        }
    }

    // consider solid boundaries 
    // TODO -> our solids are only constant velocity rn

    // set entries of A matrix for the pressure solver
    double dt_dx2 = dt / (this->cell_width * this->cell_width);
    double dt_dy2 = dt / (this->cell_height * this->cell_height);
    double dt_dz2 = dt / (this->cell_depth * this->cell_depth);
    for (int depth_idx = 0; depth_idx < DEPTH_NUM; depth_idx ++)
    {    
        for (int row_idx = 0; row_idx < ROW_NUM; row_idx ++)
        {
            for (int col_idx = 0; col_idx < COL_NUM; col_idx ++)
            {
                if (this->_full_fluid_cell({row_idx, col_idx, depth_idx})){
                grid_idx_t my_idx = {row_idx, col_idx, depth_idx};
                // horizontal neighbors
                grid_idx_t lneighbor = _lneighbor(my_idx);
                grid_idx_t rneighbor = _rneighbor(my_idx);

                if(_full_fluid_cell(lneighbor)){
                    // should be valid index since a fluid cell
                    // double ldensity = this->cells[lneighbor.row][lneighbor.col].density; 
                    this->set_A_idx(my_idx, LEFT, -dt_dx2 / density);
                    this->add_to_A(my_idx, CENTER, dt_dx2 / density);
                }
                // still add to pressure solve (just no pressure from this dir)
                // else if(_air_cell(lneighbor)){
                //     this->add_to_A(my_idx, CENTER, dt_dx2 / density);
                // }

                if(_full_fluid_cell(rneighbor)){
                    // double rdensity = this->cells[rneighbor.row][rneighbor.col].density; 
                    this->set_A_idx(my_idx, RIGHT, -dt_dx2 / density);
                    this->add_to_A(my_idx, CENTER, dt_dx2 / density);
                }
                // else if(_air_cell(rneighbor)){
                //     this->add_to_A(my_idx, CENTER, dt_dx2 / density);
                // }

                // vertical neighbors
                grid_idx_t tneighbor = _tneighbor(my_idx);
                grid_idx_t bneighbor = _bneighbor(my_idx);
                if(_full_fluid_cell(tneighbor)){
                    // double tdensity = this->cells[tneighbor.row][tneighbor.col].density; 
                    this->set_A_idx(my_idx, TOP, -dt_dy2 / density);
                    this->add_to_A(my_idx, CENTER, dt_dy2 / density);
                }
                // else if(_air_cell(tneighbor)){
                //     this->add_to_A(my_idx, CENTER, dt_dy2 / density);
                // }

                if(_full_fluid_cell(bneighbor)){
                    // double bdensity = this->cells[bneighbor.row][bneighbor.col].density; 
                    this->set_A_idx(my_idx, BOTTOM, -dt_dy2 / density);
                    this->add_to_A(my_idx, CENTER, dt_dy2 / density);
                }
                // else if(_air_cell(bneighbor)){
                //     this->add_to_A(my_idx, CENTER, dt_dy2 / density);
                // }

                // depth neighbors
                grid_idx_t baneighbor = _baneighbor(my_idx);
                grid_idx_t fneighbor = _fneighbor(my_idx);
                if(_full_fluid_cell(baneighbor)){
                    // double tdensity = this->cells[tneighbor.row][tneighbor.col].density; 
                    this->set_A_idx(my_idx, BACK, -dt_dz2 / density);
                    this->add_to_A(my_idx, CENTER, dt_dz2 / density);
                }
                // else if(_air_cell(baneighbor)){
                //     this->add_to_A(my_idx, CENTER, dt_dz2 / density);
                // }

                if(_full_fluid_cell(fneighbor)){
                    // double bdensity = this->cells[bneighbor.row][bneighbor.col].density; 
                    this->set_A_idx(my_idx, FRONT, -dt_dz2 / density);
                    this->add_to_A(my_idx, CENTER, dt_dz2 / density);
                }
                // else if(_air_cell(fneighbor)){
                //     this->add_to_A(my_idx, CENTER, dt_dz2 / density);
                // }
                }
            }
        }
    }

    // compute the final p matrix
    if(!this->solve_with_PCG())
    {
        printf("ERROR SOLVING PCG\n");
    }

    // distribute pVec to the pressure grid
    // TODO -> idt this is actually necessary but might be nice for viz ?
    for (int depth_idx = 0; depth_idx < DEPTH_NUM; depth_idx ++)
    {    
        for (int row_idx = 0; row_idx < ROW_NUM; row_idx ++)
        {
            for (int col_idx = 0; col_idx < COL_NUM; col_idx ++)
            {
                // printf("%f\n", this->pVec[get_flat_idx({row_idx, col_idx, depth_idx})]);
                if(this->_full_fluid_cell({row_idx, col_idx, depth_idx})) // had non-zero density so could solve for its value
                {
                    this->pressure_grid[row_idx][col_idx][depth_idx].value = this->pVec[get_flat_idx({row_idx, col_idx, depth_idx})];
                }
                else
                {
                    this->pressure_grid[row_idx][col_idx][depth_idx].value = 0;
                }
            }
        }
    }

    // iterate and interpolate "unknown" gas ghost pressures
    for (int depth_idx = 0; depth_idx < DEPTH_NUM; depth_idx ++)
    {    
    for (int row_idx = 0; row_idx < ROW_NUM; row_idx ++)
    {
        for (int col_idx = 0; col_idx < COL_NUM; col_idx ++)
        {
            if(this->_air_cell({row_idx, col_idx, depth_idx}))
            {
                grid_idx_t my_idx = {row_idx, col_idx, depth_idx};

                // horizontal neighbors
                grid_idx_t lneighbor = _lneighbor(my_idx);
                grid_idx_t rneighbor = _rneighbor(my_idx);
                // vertical neighbors
                grid_idx_t tneighbor = _tneighbor(my_idx);
                grid_idx_t bneighbor = _bneighbor(my_idx);
                // depth neighbors
                grid_idx_t baneighbor = _baneighbor(my_idx);
                grid_idx_t fneighbor = _fneighbor(my_idx);

                int denom = _full_fluid_cell(lneighbor) + _full_fluid_cell(rneighbor) 
                + _full_fluid_cell(tneighbor) + _full_fluid_cell(bneighbor)
                + _full_fluid_cell(fneighbor) + _full_fluid_cell(baneighbor);

                if(_full_fluid_cell(lneighbor))
                    this->pressure_grid[row_idx][col_idx][depth_idx].value -= this->pressure_grid[lneighbor.row][lneighbor.col][lneighbor.depth].value/denom;
                if(_full_fluid_cell(rneighbor))
                    this->pressure_grid[row_idx][col_idx][depth_idx].value -= this->pressure_grid[rneighbor.row][rneighbor.col][rneighbor.depth].value/denom;
                if(_full_fluid_cell(tneighbor))
                    this->pressure_grid[row_idx][col_idx][depth_idx].value -= this->pressure_grid[tneighbor.row][tneighbor.col][tneighbor.depth].value/denom;
                if(_full_fluid_cell(bneighbor))
                    this->pressure_grid[row_idx][col_idx][depth_idx].value -= this->pressure_grid[bneighbor.row][bneighbor.col][bneighbor.depth].value/denom;
                if(_full_fluid_cell(baneighbor))
                    this->pressure_grid[row_idx][col_idx][depth_idx].value -= this->pressure_grid[baneighbor.row][baneighbor.col][baneighbor.depth].value/denom;
                if(_full_fluid_cell(fneighbor))
                    this->pressure_grid[row_idx][col_idx][depth_idx].value -= this->pressure_grid[fneighbor.row][fneighbor.col][fneighbor.depth].value/denom;
            }
        }
    }
    }

    // update velocities for every non-air cell
    double dt_dx = dt / (this->cell_width);
    double dt_dy = dt / (this->cell_height);
    double dt_dz = dt / (this->cell_depth);

    for (int depth_idx = 0; depth_idx < DEPTH_NUM; depth_idx ++)
    {
        for(int row_idx = 0; row_idx<ROW_NUM; row_idx++)
        {
            for(int col_idx = 0; col_idx<COL_NUM; col_idx++)
            {
                grid_idx_t cell_idx = {row_idx, col_idx, depth_idx};
                grid_idx_t lneigh = this->_lneighbor(cell_idx);
                grid_idx_t bneigh = this->_bneighbor(cell_idx);
                grid_idx_t fneigh = this->_fneighbor(cell_idx);

                Vertex& ui0 = this->horizontal_velocity[row_idx][col_idx][depth_idx];
                Vertex& vi0 = this->vertical_velocity[row_idx][col_idx][depth_idx];
                Vertex& zi0 = this->depth_velocity[row_idx][col_idx][depth_idx];

                // Update horizontal velocity
                if(this->_solid_hcell(cell_idx)) {
                    ui0.value = 0.0;
                } else if(_full_fluid_cell(cell_idx) || _full_fluid_cell(lneigh)) {
                    double hvalue = dt_dx * (pressure_grid[lneigh.row][lneigh.col][lneigh.depth].value - 
                        pressure_grid[row_idx][col_idx][depth_idx].value) / density;
                    ui0.value += hvalue; // Add pressure gradient
                } else {
                    ui0.value = 0.0;
                }

                // Update vertical velocity
                if(this->_solid_vcell(cell_idx)) {
                    vi0.value = 0.0;
                } else if(_full_fluid_cell(cell_idx) || _full_fluid_cell(bneigh)) {
                    double vvalue = dt_dy * (pressure_grid[bneigh.row][bneigh.col][bneigh.depth].value - 
                        pressure_grid[row_idx][col_idx][depth_idx].value) / density;
                    vi0.value += vvalue; // Add pressure gradient
                } else {
                    vi0.value = 0.0;
                }

                // Update depth velocity
                if(this->_solid_dcell(cell_idx)) {
                    zi0.value = 0.0;
                } else if(_full_fluid_cell(cell_idx) || _full_fluid_cell(fneigh)) {
                    double dvalue = dt_dz * (pressure_grid[fneigh.row][fneigh.col][fneigh.depth].value - 
                        pressure_grid[row_idx][col_idx][depth_idx].value) / density;
                    zi0.value += dvalue; // Add pressure gradient
                } else {
                    zi0.value = 0.0;
                }
            }

            // Set boundary velocities to zero
            Vertex& ui0 = this->horizontal_velocity[row_idx][COL_NUM][depth_idx];
            ui0.value = 0;
        }

        for(size_t col_idx = 0; col_idx<COL_NUM; col_idx++)
        { 
            Vertex& vi0 = this->vertical_velocity[ROW_NUM][col_idx][depth_idx];
            vi0.value = 0;
        }
    }

    for(size_t row_idx = 0; row_idx<ROW_NUM; row_idx++)
    {
        for(size_t col_idx = 0; col_idx<COL_NUM; col_idx++)
        {
            Vertex& zi0 = this->depth_velocity[row_idx][col_idx][DEPTH_NUM];
            zi0.value = 0;
        }
    }

    // print the vector
    #if DEBUG
    for (int depth_idx = 0; depth_idx < DEPTH_NUM; depth_idx ++)
    {
    for (int row_idx = 0; row_idx < ROW_NUM; row_idx ++)
    {
        printf("row %d: \n", row_idx);
        for (int col_idx = 0; col_idx < COL_NUM; col_idx ++)
        {
            grid_idx_t my_idx = {row_idx, col_idx, depth_idx};
            printf("[%d]: %f ", col_idx, this->pVec[this->get_flat_idx(my_idx)]);
        }
    }
    }

    printf("\n");
    #endif
}

void Grid::print_grid()
{
    for (size_t depth_idx = 0; depth_idx < DEPTH_NUM; depth_idx ++)
    {
    for(size_t row_idx = 0; row_idx<ROW_NUM; row_idx++)
    {
        for(size_t col_idx = 0; col_idx<COL_NUM; col_idx++)
        {
            char cell_str = this->cells[row_idx][col_idx][depth_idx].type ==FLUID ? 'F' : cells[row_idx][col_idx][depth_idx].type == GAS ? 'G' : 'S';
            printf("[%zu, %zu, %zu] (%c):\n", row_idx, col_idx, depth_idx, cell_str);
            printf("    pos: %s\n", cells[row_idx][col_idx][depth_idx].position.print().c_str());
            printf("    p_ij: %f\n", pressure_grid[row_idx][col_idx][depth_idx].value);
            printf("    u, j: %f; v, j: %f; z: %f \n", horizontal_velocity[row_idx][col_idx][depth_idx].value, vertical_velocity[row_idx][col_idx][depth_idx].value,
                depth_velocity[row_idx][col_idx][depth_idx].value);
            printf("    density: %f\n", cells[row_idx][col_idx][depth_idx].density);
        }
    }
    }
}