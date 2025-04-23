#include <string.h>
#include <cassert>
#include "grid.hpp"


Grid::Grid(double width, double height, double density = 1.0) {
    this->width = width;
    this->height = height;
    this->cell_width = width/COL_NUM;
    this->cell_height = height/ROW_NUM;
    this->cell_volume = this->cell_width * this->cell_height;
    this->density = density;

    // initialize the "poses" for each vertex
    for(size_t row_idx = 0; row_idx<ROW_NUM; row_idx++)
    {
        for(size_t col_idx = 0; col_idx<COL_NUM; col_idx++)
        {
            // set the center position of each cell
            double x = (col_idx+0.5)*this->cell_width;
            double y = (row_idx+0.5)*this->cell_height;
            this->cells[row_idx][col_idx].position = Point(x, y, 0.0);
            // printf(this->cells[row_idx][col_idx].position.print().c_str());

            Vertex& vi0 = this->vertical_velocity[row_idx][col_idx];
            double vx = this->cell_width*(col_idx+0.5);
            double vy = this->cell_height*row_idx;
            vi0.position = Point(vx, vy, 0.0);

            Vertex& ui0 = this->horizontal_velocity[row_idx][col_idx];
            double ux = this->cell_width*col_idx;
            double uy = this->cell_height*(row_idx+0.5);
            ui0.position = Point(ux, uy, 0.0);
        }

        Vertex& ui0 = this->horizontal_velocity[row_idx][COL_NUM];
        double ux = this->cell_width*COL_NUM;
        double uy = this->cell_height*(row_idx+0.5);
        ui0.position =  Point(ux, uy, 0.0);
    }

    for(size_t col_idx = 0; col_idx<COL_NUM; col_idx++)
    { 
        Vertex& ui0 = this->vertical_velocity[ROW_NUM][col_idx];
        double ux = this->cell_width*(col_idx+0.5);
        double uy = this->cell_height*ROW_NUM;
        ui0.position = Point(ux, uy, 0.0);
    }

    // initialize pressure solver
    this->dV.resize(ROW_NUM * COL_NUM);
    this->sparseA.resize(ROW_NUM * COL_NUM, std::vector<double>(SPARSE_WIDTH));
    this->pVec.resize(ROW_NUM * COL_NUM);
    this->diagE.resize(ROW_NUM * COL_NUM);
};

// returns -> [row_idx, col_idx]
grid_idx_t Grid::get_grid_idx(Particle p)
{
    return std::pair<int, int>(
        static_cast<int>(p.position.y/this->cell_height),
        static_cast<int>(p.position.x/this->cell_width)
    );
}

// returns the four cell vertex indices bounding the particle
Neighbors Grid::get_cell_neighbors(Particle p)
{
    // int row_start_idx = static_cast<int>((this->height - p.position.y)/this->cell_height - 0.5);
    // int col_startidx = static_cast<int>(p.position.x/this->cell_width - 0.5);

    // Neighbors n {};
    // n.type = CELL;

    // // we dont check the value for negativity here, since negative implies solid boundary -> info we might need
    // n.neighbors[0] = std::pair<int, int>(row_idx, col_idx);
    // n.neighbors[1] = std::pair<int, int>(row_idx+1, col_idx);
    // n.neighbors[2] = std::pair<int, int>(row_idx, col_idx+1);
    // n.neighbors[3] = std::pair<int, int>(row_idx+1, col_idx+1);

    // return n;

    grid_idx_t cell_idx = this->get_grid_idx(p);

    Neighbors n {};
    n.type = CELL;

    // we dont check the value for negativity here, since negative implies solid boundary -> info we might need
    n.neighbors[0] = std::pair<int, int>(cell_idx.first - 1, cell_idx.second);
    n.neighbors[1] = std::pair<int, int>(cell_idx.first, cell_idx.second - 1);
    n.neighbors[2] = std::pair<int, int>(cell_idx.first, cell_idx.second + 1);
    n.neighbors[3] = std::pair<int, int>(cell_idx.first + 1, cell_idx.second);

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

            // printf("(%d, %d): %f\n ", row_idx, col_idx, this->cells[row_idx][col_idx].density);
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
        // at least one particle in cell -> the cell becomes fluid
        grid_idx_t cell_idx = this->get_grid_idx(p);
        Cell& cell = this->cells[cell_idx.first][cell_idx.second];
        cell.type = FLUID;
        double dist = p.position.bilinear(cell.position, this->cell_width, this->cell_height);
        cell.density += (p.mass) * dist;

        // get the neighboring horizontal velocity cells
        Neighbors nh = this->get_horizontal_neighbors(p);
        for(grid_idx_t cell_idx : nh.neighbors){
            // printf("cell (%d, %d)\n", cell_idx.first, cell_idx.second);
            // check if neighbor is border. if so, do nothing -> this will be a solid cell
            if(!this->_valid_hcell(cell_idx)) continue;

            #ifdef DEBUG
            assert(cell_idx.second >= 0 && cell_idx.second <= COL_NUM);
            #endif

            Vertex& ui0 = this->horizontal_velocity[cell_idx.first][cell_idx.second];

            // if this is a solid cell, we keep its velocity 0 (later)
            // if(this->_solid_hcell(cell_idx)){ ui0.value = 0; continue;}

            // then we set horizontal and vertical velocities based off of this for each grid vertex
            double dist = p.position.bilinear(ui0.position, this->cell_width, this->cell_height);
            ui0.value += p.mass * dist * p.velocity.x; // TODO -> should this be bilinear? just 4 closest neighbors?
            // ui0.normalization += p.mass; // TODO -> normalize per particles
            // ui0.normalization = NUM_PARTICLES / (this->cell_width * this->cell_height);
            ui0.normalization += dist;
        } 

        // do the same for vertical
        Neighbors nv = this->get_vertical_neighbors(p);
        for(grid_idx_t cell_idx : nv.neighbors){
            // check if neighbor is border. if so, do nothing -> this will be a solid cell
            if(!this->_valid_vcell(cell_idx)) continue;

            #ifdef DEBUG
            assert(cell_idx.first >= 0 && cell_idx.first <= ROW_NUM);
            #endif

            Vertex& vi0 = this->vertical_velocity[cell_idx.first][cell_idx.second];

            // if this is a solid cell, we keep its velocity 0 (later)
            // if(this->_solid_vcell(cell_idx)){ vi0.value = 0; continue;}

            double dist =  p.position.bilinear(vi0.position, this->cell_width, this->cell_height);
            vi0.value += p.mass * dist * p.velocity.y;
            vi0.normalization += dist;
            // vi0.normalization += p.mass;
            // vi0.normalization = NUM_PARTICLES / (this->cell_width * this->cell_height);
        }

        // calculate density by interpolating
        Neighbors n = this->get_cell_neighbors(p);
        for(grid_idx_t cell_idx : n.neighbors){
            // check if neighbor is border. if so, do nothing -> this will be a solid cell (no density)
            if(!this->_valid_cell(cell_idx)) continue;

            Cell& cell = this->cells[cell_idx.first][cell_idx.second];
            cell.type = FLUID;

            // printf("[%zu, %zu] (%f, %f):\n", cell_idx.first, cell_idx.second, v.position.x, v.position.y);
            // printf("             (%f, %f):\n",  p.position.x, p.position.y);
            // printf("%f\n", p.position.bilinear(v.position, this->cell_width, this->cell_height));
            // assert(p.position.bilinear(cell.position, this->cell_width, this->cell_height) >= 0 &&  p.position.bilinear(cell.position, this->cell_width, this->cell_height) <= 1);
            double dist = p.position.bilinear(cell.position, this->cell_width, this->cell_height);
            cell.density += (p.mass / this->cell_volume) * dist;
            // cell.normalization += dist;
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
        double velocity_norm = 0.0;
        for(grid_idx_t cell_idx : nh.neighbors){
            if(!this->_valid_hcell(cell_idx)) continue;

            Vertex& u = this->horizontal_velocity[cell_idx.first][cell_idx.second];
            double dist = p.position.bilinear(u.position, this->cell_width, this->cell_height);
            p.velocity.x += PIC_WEIGHT * dist * u.value; // PIC update
            p.velocity.x += (1 - PIC_WEIGHT) * dist * (u.value - u.prev_value); // FLIP update
            velocity_norm += dist;
        }

        // normalize
        p.velocity.x /= velocity_norm;

        Neighbors nv = this->get_vertical_neighbors(p);
        velocity_norm = 0.0;
        for(grid_idx_t cell_idx : nv.neighbors){
            if(!this->_valid_vcell(cell_idx)) continue;

            Vertex& v = this->vertical_velocity[cell_idx.first][cell_idx.second];
            double dist = p.position.bilinear(v.position, this->cell_width, this->cell_height);
            p.velocity.y += PIC_WEIGHT * dist * v.value; // PIC update
            p.velocity.y += (1 - PIC_WEIGHT) * dist * (v.value - v.prev_value); // FLIP update
            velocity_norm += dist;
        }

        p.velocity.y /= velocity_norm;
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

    // printf("reset solver\n");

    // populate right hand-side divergence of current velocity field
    double dx = this->cell_width;
    double dy = this->cell_height;
    for (int row_idx = 0; row_idx < ROW_NUM; row_idx ++)
    {
        for (int col_idx = 0; col_idx < COL_NUM; col_idx ++)
        {
            // todo -> chat is this right? fluid or air?
            if (this->_full_fluid_cell({row_idx, col_idx})){
                // printf("FOUND A FLUID CELL\n");
                // this should always be okay since borders contain leeway in indices
                double dU = (this->horizontal_velocity[row_idx][col_idx+1].value - this->horizontal_velocity[row_idx][col_idx].value) / dx;
                double dV = (this->vertical_velocity[row_idx+1][col_idx].value - this->vertical_velocity[row_idx][col_idx].value) / dy;

                // printf("%f\n", (dU + dV));
                this->set_dV_idx({row_idx, col_idx}, -(dU + dV));
            
                // modify solid bordering fluid cells
                // TODO -> this assumes solid always moves at velocity 0
                if (this->_solid_hborder_left({row_idx, col_idx})){
                    double hv = this->horizontal_velocity[row_idx][col_idx].value;
                    this->add_to_dV({row_idx, col_idx}, - hv / dx);
                }
                if (this->_solid_hborder_right({row_idx, col_idx})){
                    double hv = this->horizontal_velocity[row_idx][col_idx + 1].value;
                    this->add_to_dV({row_idx, col_idx}, - hv / dx);
                }
                if (this->_solid_vborder_top({row_idx, col_idx})){
                    double hv = this->vertical_velocity[row_idx + 1][col_idx].value;
                    this->add_to_dV({row_idx, col_idx}, hv / dy);
                }
                if (this->_solid_vborder_bot({row_idx, col_idx})){
                    double hv = this->vertical_velocity[row_idx][col_idx].value;
                    this->add_to_dV({row_idx, col_idx}, - hv / dy);
                } 
            }
        }
    }

    // printf("populated right\n");

    // consider solid boundaries 
    // TODO -> our solids are only constant velocity rn

    // set entries of A matrix for the pressure solver
    double dt_dx2 = dt / (this->cell_width * this->cell_width);
    double dt_dy2 = dt / (this->cell_height * this->cell_height);
    for (int row_idx = 0; row_idx < ROW_NUM; row_idx ++)
    {
        for (int col_idx = 0; col_idx < COL_NUM; col_idx ++)
        {
            if (this->_full_fluid_cell({row_idx, col_idx})){
            grid_idx_t my_idx = {row_idx, col_idx};

            // TODO -> should this be an interpolation
            // double density = this->cells[row_idx][col_idx].density; 
            // horizontal neighbors
            grid_idx_t lneighbor = _lneighbor(my_idx);
            grid_idx_t rneighbor = _rneighbor(my_idx);
            if(_full_fluid_cell(lneighbor)){
                // should be valid index since a fluid cell
                // double ldensity = this->cells[lneighbor.first][lneighbor.second].density; 
                this->set_A_idx(my_idx, LEFT, -dt_dx2 / density);
                this->add_to_A(my_idx, CENTER, dt_dx2 / density);
            }
            // still add to pressure solve (just no pressure from this dir)
            else if(_air_cell(lneighbor)){
                this->add_to_A(my_idx, CENTER, dt_dx2 / density);
            }
            if(_full_fluid_cell(rneighbor)){
                // double rdensity = this->cells[rneighbor.first][rneighbor.second].density; 
                this->set_A_idx(my_idx, RIGHT, -dt_dx2 / density);
                this->add_to_A(my_idx, CENTER, dt_dx2 / density);
            }
            else if(_air_cell(rneighbor)){
                this->add_to_A(my_idx, CENTER, dt_dx2 / density);
            }

            // vertical neighbors
            grid_idx_t tneighbor = _tneighbor(my_idx);
            grid_idx_t bneighbor = _bneighbor(my_idx);
            if(_full_fluid_cell(tneighbor)){
                // double tdensity = this->cells[tneighbor.first][tneighbor.second].density; 
                this->set_A_idx(my_idx, TOP, -dt_dy2 / density);
                this->add_to_A(my_idx, CENTER, dt_dy2 / density);
            }
            else if(_air_cell(tneighbor)){
                this->add_to_A(my_idx, CENTER, dt_dy2 / density);
            }

            if(_full_fluid_cell(bneighbor)){
                // double bdensity = this->cells[bneighbor.first][bneighbor.second].density; 
                this->set_A_idx(my_idx, BOTTOM, -dt_dy2 / density);
                this->add_to_A(my_idx, CENTER, dt_dy2 / density);
            }
            else if(_air_cell(bneighbor)){
                this->add_to_A(my_idx, CENTER, dt_dy2 / density);
            }
            }
            // the pressure is just copied from its fluid cell neighbors (if exist, o.w. 0)
            // else if(this->_air_cell({row_idx, col_idx}))
            // {
            //     grid_idx_t my_idx = {row_idx, col_idx};

            //     // horizontal neighbors
            //     grid_idx_t lneighbor = _lneighbor(my_idx);
            //     grid_idx_t rneighbor = _rneighbor(my_idx);
            //     // vertical neighbors
            //     grid_idx_t tneighbor = _tneighbor(my_idx);
            //     grid_idx_t bneighbor = _bneighbor(my_idx);

            //     if(_full_fluid_cell(lneighbor)){
            //         double ldensity = this->cells[lneighbor.first][lneighbor.second].density; 
            //         this->set_A_idx(my_idx, LEFT, -dt_dx2 / ldensity);
            //     }
            //     if(_full_fluid_cell(rneighbor)){
            //         double rdensity = this->cells[rneighbor.first][rneighbor.second].density; 
            //         this->set_A_idx(my_idx, RIGHT, -dt_dx2 / rdensity);
            //     }
            //     if(_full_fluid_cell(tneighbor)){
            //         double tdensity = this->cells[tneighbor.first][tneighbor.second].density; 
            //         this->set_A_idx(my_idx, TOP, -dt_dy2 / tdensity);
            //     }
            //     if(_full_fluid_cell(bneighbor)){
            //         double bdensity = this->cells[bneighbor.first][bneighbor.second].density; 
            //         this->set_A_idx(my_idx, BOTTOM, -dt_dy2 / bdensity);
            //     }
            // }
        }
    }

    // printf("built A\n");

    // compute the final p matrix
    if(!this->solve_with_PCG())
    {
        printf("ERROR SOLVING PCG\n");
    }

    // distribute pVec to the pressure grid
    // TODO -> idt this is actually necessary but might be nice for viz ?
    for (int row_idx = 0; row_idx < ROW_NUM; row_idx ++)
    {
        for (int col_idx = 0; col_idx < COL_NUM; col_idx ++)
        {
            // printf("%f\n", this->pVec[get_flat_idx({row_idx, col_idx})]);
            if(this->_full_fluid_cell({row_idx, col_idx})) // had non-zero density so could solve for its value
            {
                this->pressure_grid[row_idx][col_idx].value = this->pVec[get_flat_idx({row_idx, col_idx})];
            }
            else
            {
                this->pressure_grid[row_idx][col_idx].value = 0;
            }
        }
    }

    // iterate and interpolate "unknown" gas ghost pressures
    for (int row_idx = 0; row_idx < ROW_NUM; row_idx ++)
    {
        for (int col_idx = 0; col_idx < COL_NUM; col_idx ++)
        {
            if(this->_air_cell({row_idx, col_idx}))
            {
                grid_idx_t my_idx = {row_idx, col_idx};

                // horizontal neighbors
                grid_idx_t lneighbor = _lneighbor(my_idx);
                grid_idx_t rneighbor = _rneighbor(my_idx);
                // vertical neighbors
                grid_idx_t tneighbor = _tneighbor(my_idx);
                grid_idx_t bneighbor = _bneighbor(my_idx);

                int denom = _full_fluid_cell(lneighbor) + _full_fluid_cell(rneighbor) 
                + _full_fluid_cell(tneighbor) + _full_fluid_cell(bneighbor);

                if(_full_fluid_cell(lneighbor))
                    this->pressure_grid[row_idx][col_idx].value -= this->pressure_grid[lneighbor.first][lneighbor.second].value/denom;
                if(_full_fluid_cell(rneighbor))
                    this->pressure_grid[row_idx][col_idx].value -= this->pressure_grid[rneighbor.first][rneighbor.second].value/denom;
                if(_full_fluid_cell(tneighbor))
                    this->pressure_grid[row_idx][col_idx].value -= this->pressure_grid[tneighbor.first][tneighbor.second].value/denom;
                if(_full_fluid_cell(bneighbor))
                    this->pressure_grid[row_idx][col_idx].value -= this->pressure_grid[bneighbor.first][bneighbor.second].value/denom;
            }
        }
    }

    // update velocities for every non-air cell
    double dt_dx = dt / (this->cell_width);
    double dt_dy = dt / (this->cell_height);

    for(size_t row_idx = 0; row_idx<ROW_NUM; row_idx++)
    {
        for(size_t col_idx = 0; col_idx<COL_NUM; col_idx++)
        {
            // TODO -> should this be an interpolation
            // double density = this->cells[row_idx][col_idx].density; 
            // double density = 1;

            grid_idx_t cell_idx = {row_idx, col_idx};
            grid_idx_t lneigh = this->_lneighbor(cell_idx);
            grid_idx_t bneigh = this->_bneighbor(cell_idx);

            Vertex& ui0 = this->horizontal_velocity[row_idx][col_idx];
            Vertex& vi0 = this->vertical_velocity[row_idx][col_idx];

            double hvalue = 0.0;
            // if horizontally borders a solid cell, horizontal velocity is always zero
            if(this->_solid_hcell(cell_idx)) 
                ui0.value = hvalue;
            else 
            if(_full_fluid_cell(cell_idx) || _full_fluid_cell(lneigh)){
                // double avgd = (density + this->cells[lneigh.first][lneigh.second].density)/2;
                hvalue = dt_dx * (pressure_grid[row_idx][col_idx].value - 
                    pressure_grid[lneigh.first][lneigh.second].value) / density;

                // printf("hvalue: %f\n", hvalue);
                ui0.value -= hvalue;
            }
            else 
                ui0.value = hvalue;

            double vvalue = 0.0;
            // if vertically borders a solid cell, vertical velocity is always zero
            if(this->_solid_vcell(cell_idx)) 
                vi0.value = hvalue;
            else 
            if(_full_fluid_cell(cell_idx) || _full_fluid_cell(bneigh)){
                // double avgd = (density + this->cells[tneigh.first][tneigh.second].density)/2;
                vvalue = dt_dy * (pressure_grid[row_idx][col_idx].value - 
                    pressure_grid[bneigh.first][bneigh.second].value) / density;
                    // printf("vvalue: %f\n", vvalue);
                vi0.value -= vvalue;
            }
            else 
                vi0.value = vvalue;
        }

        // set horizontal velocities bordering solid cells
        Vertex& ui0 = this->horizontal_velocity[row_idx][COL_NUM];
        ui0.value = 0;
    }

    for(size_t col_idx = 0; col_idx<COL_NUM; col_idx++)
    { 
        // vertical velocities at the bottom
        Vertex& vi0 = this->vertical_velocity[ROW_NUM][col_idx];
        vi0.value = 0;
    }

    // print the vector
    #ifdef DEBUG
    for (int row_idx = 0; row_idx < ROW_NUM; row_idx ++)
    {
        printf("row %d: \n", row_idx);
        for (int col_idx = 0; col_idx < COL_NUM; col_idx ++)
        {
            grid_idx_t my_idx = {row_idx, col_idx};
            printf("[%d]: %f ", col_idx, this->pVec[this->get_flat_idx(my_idx)]);
        }
    }

    printf("\n");
    #endif
}

void Grid::print_grid()
{
    for(size_t row_idx = 0; row_idx<ROW_NUM; row_idx++)
    {
        for(size_t col_idx = 0; col_idx<COL_NUM; col_idx++)
        {
            char cell_str = this->cells[row_idx][col_idx].type ==FLUID ? 'F' : cells[row_idx][col_idx].type == GAS ? 'G' : 'S';
            printf("[%zu, %zu] (%c):\n", row_idx, col_idx, cell_str);
            printf("    pos: %s\n", cells[row_idx][col_idx].position.print().c_str());
            printf("    p_ij: %f\n", pressure_grid[row_idx][col_idx].value);
            printf("    u, j: %f; v, j: %f\n", horizontal_velocity[row_idx][col_idx].value, vertical_velocity[row_idx][col_idx].value);
            printf("    density: %f\n", cells[row_idx][col_idx].density);
        }
    }
}