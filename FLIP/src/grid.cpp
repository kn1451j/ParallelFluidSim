#include <string.h>
#include <cassert>
#include "grid.hpp"


Grid::Grid(double width, double height) {
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
            double x = (col_idx+0.5)*this->cell_width;
            double y = this->height - (row_idx+0.5)*this->cell_height;
            this->cells[row_idx][col_idx].position = Point(x, y, 0.0);
            // printf(this->cells[row_idx][col_idx].position.print().c_str());

            Vertex& vi0 = this->vertical_velocity[row_idx][col_idx];
            double vx = this->cell_width*(col_idx+0.5);
            double vy = this->height - this->cell_height*row_idx;
            vi0.position = Point(vx, vy, 0.0);

            Vertex& ui0 = this->horizontal_velocity[row_idx][col_idx];
            double ux = this->cell_width*col_idx;
            double uy = this->height - this->cell_height*(row_idx+0.5);
            ui0.position = Point(ux, uy, 0.0);
        }

        Vertex& ui0 = this->horizontal_velocity[row_idx][COL_NUM];
        double ux = this->cell_width*COL_NUM;
        double uy = this->height - this->cell_height*(row_idx+0.5);
        ui0.position =  Point(ux, uy, 0.0);
    }

    for(size_t col_idx = 0; col_idx<COL_NUM; col_idx++)
    { 
        Vertex& ui0 = this->vertical_velocity[ROW_NUM][col_idx];
        double ux = this->cell_width*(col_idx+0.5);
        double uy = this->height - this->cell_height*ROW_NUM;
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
        static_cast<int>((this->height - p.position.y)/this->cell_height),
        static_cast<int>(p.position.x/this->cell_width)
    );
}

// returns the 4 horizontal cell vertex indices bounding the particle
Neighbors Grid::get_cell_neighbors(Particle p)
{
    int row_idx = static_cast<int>((this->height - p.position.y)/this->cell_height - 0.5);
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
    int row_idx = static_cast<int>((this->height - p.position.y)/this->cell_height - 0.5);
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
    int row_idx = static_cast<int>((this->height - p.position.y)/this->cell_height);
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
        // DENSITY CALCULATED LATER

        // get the neighboring horizontal velocity cells
        Neighbors nh = this->get_horizontal_neighbors(p);
        for(grid_idx_t cell_idx : nh.neighbors){
            // printf("cell (%d, %d)\n", cell_idx.first, cell_idx.second);
            // check if neighbor is border. if so, do nothing -> this will be a solid cell
            if(!this->_valid_cell(cell_idx)) continue;

            #ifdef DEBUG
            assert(cell_idx.second >= 0 && cell_idx.second <= COL_NUM);
            #endif

            Vertex& ui0 = this->horizontal_velocity[cell_idx.first][cell_idx.second];

            // then we set horizontal and vertical velocities based off of this for each grid vertex
            double dist = p.position.bilinear(ui0.position, this->cell_width, this->cell_height);
            ui0.value += p.mass * p.velocity.x * dist; // TODO -> should this be bilinear? just 4 closest neighbors?
            ui0.normalization += dist; // TODO -> normalize per particles
        } 

        // do the same for vertical
        Neighbors nv = this->get_vertical_neighbors(p);
        for(grid_idx_t cell_idx : nv.neighbors){
            // check if neighbor is border. if so, do nothing -> this will be a solid cell
            if(!this->_valid_cell(cell_idx)) continue;

            #ifdef DEBUG
            assert(cell_idx.first >= 0 && cell_idx.first <= ROW_NUM);
            #endif

            Vertex& vi0 = this->vertical_velocity[cell_idx.first][cell_idx.second];
            double dist =  p.position.bilinear(vi0.position, this->cell_width, this->cell_height);
            vi0.value += p.mass * p.velocity.y * dist;
            vi0.normalization += dist;
            // vi0.normalization += p.mass;
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
            cell.density += (p.mass/this->cell_volume) * dist;
            cell.normalization += dist;
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
            if(!this->_fluid_cell(cell_idx)) continue;

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
            if(!this->_fluid_cell(cell_idx)) continue;

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
    for (int row_idx = 0; row_idx < ROW_NUM; row_idx ++)
    {
        for (int col_idx = 0; col_idx < COL_NUM; col_idx ++)
        {
            if (this->_fluid_cell({row_idx, col_idx})){
                // printf("FOUND A FLUID CELL\n");
                // this should always be okay since borders contain leeway
                double dU = (this->horizontal_velocity[row_idx][col_idx+1].value - this->horizontal_velocity[row_idx][col_idx].value)/this->cell_width;
                double dV = (this->vertical_velocity[row_idx+1][col_idx].value - this->vertical_velocity[row_idx][col_idx].value)/this->cell_height;

                // printf("%f\n", (dU + dV));
                this->set_dV_idx({row_idx, col_idx}, (dU + dV));
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
            grid_idx_t my_idx = {row_idx, col_idx};

            // TODO -> should this be an interpolation
            double density = this->cells[row_idx][col_idx].density; 
            // horizontal neighbors
            grid_idx_t lneighbor = _lneighbor(my_idx);
            grid_idx_t rneighbor = _rneighbor(my_idx);
            if(_fluid_cell(lneighbor)){
                // should be valid index since a fluid cell
                // double density = this->cells[]
                this->set_A_idx(my_idx, LEFT, -dt_dx2 / density);
                this->add_to_A(my_idx, CENTER, dt_dx2 / density);
            }
            // still add to pressure solve (just no pressure from this dir)
            else if(_air_cell(lneighbor)){
                this->add_to_A(my_idx, CENTER, dt_dx2 / density);
            }
            if(_fluid_cell(rneighbor)){
                this->set_A_idx(my_idx, RIGHT, -dt_dx2 / density);
                this->add_to_A(my_idx, CENTER, dt_dx2 / density);
            }
            else if(_air_cell(rneighbor)){
                this->add_to_A(my_idx, CENTER, dt_dx2 / density);
            }

            // vertical neighbors
            grid_idx_t tneighbor = _tneighbor(my_idx);
            grid_idx_t bneighbor = _bneighbor(my_idx);
            if(_fluid_cell(tneighbor)){
                this->set_A_idx(my_idx, TOP, -dt_dy2 / density);
                this->add_to_A(my_idx, CENTER, dt_dy2 / density);
            }
            else if(_air_cell(tneighbor)){
                this->add_to_A(my_idx, CENTER, dt_dy2 / density);
            }

            if(_fluid_cell(bneighbor)){
                this->set_A_idx(my_idx, BOTTOM, -dt_dy2 / density);
                this->add_to_A(my_idx, CENTER, dt_dy2 / density);
            }
            else if(_air_cell(bneighbor)){
                this->add_to_A(my_idx, CENTER, dt_dy2 / density);
            }
        }
    }

    // printf("built A\n");

    // compute the final p matrix
    if(!this->solve_with_PCG())
    {
        printf("ERROR SOLVING PCG\n");
    }

    // distribute pVec to the pressure grid
    // TODO -> this isn't actually necessary but might be nice for viz ?
    for (int row_idx = 0; row_idx < ROW_NUM; row_idx ++)
    {
        for (int col_idx = 0; col_idx < COL_NUM; col_idx ++)
        {
            // printf("%f\n", this->pVec[get_flat_idx({row_idx, col_idx})]);
            this->pressure_grid[row_idx][col_idx].value = this->pVec[get_flat_idx({row_idx, col_idx})];
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
            double density = this->cells[row_idx][col_idx].density; 

            grid_idx_t cell_idx = {row_idx, col_idx};
            grid_idx_t lneigh = this->_lneighbor(cell_idx);
            grid_idx_t tneigh = this->_tneighbor(cell_idx);

            Vertex& ui0 = this->horizontal_velocity[row_idx][col_idx];
            Vertex& vi0 = this->vertical_velocity[row_idx][col_idx];

            double hvalue = 0.0;
            if(_fluid_cell(cell_idx) && _fluid_cell(lneigh)){
                hvalue = dt_dx * (pressure_grid[row_idx][col_idx].value - 
                    pressure_grid[lneigh.first][lneigh.second].value) / density;

                // printf("hvalue: %f\n", hvalue);
                ui0.value -= hvalue;
            }
            else 
                ui0.value = hvalue;

            double vvalue = 0.0;
            if(_fluid_cell(cell_idx) && _fluid_cell(tneigh)){
                vvalue = dt_dy * (pressure_grid[row_idx][col_idx].value - 
                    pressure_grid[tneigh.first][tneigh.second].value) / density;
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