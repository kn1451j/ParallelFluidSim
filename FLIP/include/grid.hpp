#include <vector>
#include <utility>
#include <set>
#include <utility>
#include <vector>
#include <algorithm>
#include <numeric>
#include <cmath>

#include "particle.hpp"

#define DEBUG 
// #define REALTIME 1

// Grid coarsness
#define ROW_NUM 4
#define COL_NUM 4
#define DEPTH_NUM 1

#define CLAMP 0.25
#define EPS 0.0001
#define TOL 0.001
#define MAX_ITER 100

#define PIC_WEIGHT 1.0

#define SPARSE_WIDTH 7
#define NUM_NEIGHBORS 6

enum Direction {
    LEFT = 0, 
    RIGHT = 1, 
    CENTER = 2,
    TOP = 3, 
    BOTTOM = 4,
    FRONT = 5,
    BACK = 6
};

enum SolverType
{
    PCG
};

// -1 can indicate that we're in "solid" region
struct grid_idx_t
{
    int row;
    int col;
    int depth;
};

// typedef Eigen::Matrix< double, 4, 5 > Matrix45d;

enum CELL_TYPE {FLUID, SOLID, GAS};
enum VELOCITY {HORIZONTAL, VERTICAL, CELL, DEPTH};


// TODO -> fix mass weighting
struct Cell
{
    CELL_TYPE type = GAS;

    // for density convolution
    Point position;
    double density = 0.0;
    double normalization = 0.0;

    void reset () {
        this->type = GAS;
        this->density = 0.0;
        this->normalization = 0.0;
    };
};

struct Vertex
{
    double value = 0.0;
    void reset () {
        this->prev_value = this->value;
        this->value = 0; 
        this->normalization = 0;
    };

    void set(double value) {
        this->value = value;
    };

    // used for computing particle to grid interpolation
    double normalization = 0.0;
    
    // used for flip velocity change
    double prev_value = 0.0;

    Point position = Point();
};

struct Neighbors
{
    VELOCITY type;
    
    grid_idx_t neighbors[8];
};

// Staggered Grid class for maintaining grid physics system
class Grid
{
    /*
        Need a way to store 
        1. Pressures
        2. Horizontal (x velocities)
        3. Vertical (y velocities)
    */
    public:
        // makes a staggered grid with ROW_NUM rows and COL_NUM cols
        Grid(double width, double height, double depth, double density);

        // exclusive owners of the pressure solver
        ~Grid() {};

        void transfer_to_grid(std::vector<Particle>& particles);
        void transfer_from_grid(std::vector<Particle>& particles);
        void solve_pressure(double dt);

        grid_idx_t get_grid_idx(Particle p);

        Neighbors get_horizontal_neighbors(Particle p);
        Neighbors get_vertical_neighbors(Particle p);
        Neighbors get_depth_neighbors(Particle p);
        Neighbors get_cell_neighbors(Particle p);

        /*

            Position coordinate system:

            [0][0]
            ^ ------------------  [0][COL_NUM]
            |                   |
            |                   |
            |                   |
            .-------------------> [ROW_NUM][COL_NUM]
            [ROW_NUM][0]

            [row_idx][col_idx] -> (x, y) : 
                x = col_idx * cell_width;
                y = height - row_idx * cell_height

            (x, y) -> [row_idx][col_idx] :
                row_idx = x/cell_width
                col_idx = (height - y)/cell_height
            
            idx [0][0] -> position (0, height)
            idx [ROW_NUM][0] -> position (0, 0)
            idx [0][COL_NUM] -> position (width, height)
            idx [ROW_NUM][COL_NUM] -> position (width, 0)
        */

        // TODO change for this for cache optimiality
        Vertex pressure_grid[ROW_NUM][COL_NUM][DEPTH_NUM] {};
        Vertex horizontal_velocity[ROW_NUM][COL_NUM+1][DEPTH_NUM] {};
        Vertex vertical_velocity[ROW_NUM+1][COL_NUM][DEPTH_NUM] {};
        Vertex depth_velocity[ROW_NUM][COL_NUM][DEPTH_NUM+1] {};
        Cell cells[ROW_NUM][COL_NUM][DEPTH_NUM]{};

        double width;
        double height;
        double depth;

        double cell_width;
        double cell_height;
        double cell_depth;
        double cell_volume;
        double density;

        void print_grid();

    private:
        void print_cell(grid_idx_t cell_idx) { printf("(%d, %d, %d) \n", cell_idx.row, cell_idx.col, cell_idx.depth);}

        bool _solid_vcell(grid_idx_t cell_idx)
        {
            return cell_idx.row==ROW_NUM || cell_idx.row==0;
        }

        bool _solid_hcell(grid_idx_t cell_idx)
        {
            return cell_idx.col==0 || cell_idx.col==COL_NUM;
        }

        bool _solid_dcell(grid_idx_t cell_idx)
        {
            return cell_idx.depth==0 || cell_idx.depth==DEPTH_NUM;
        }

        bool _solid_hborder_left(grid_idx_t cell_idx) {return cell_idx.col==0;};
        bool _solid_hborder_right(grid_idx_t cell_idx) {return cell_idx.col==COL_NUM - 1;};
        bool _solid_vborder_top(grid_idx_t cell_idx) {return cell_idx.row==ROW_NUM - 1;};
        bool _solid_vborder_bot(grid_idx_t cell_idx) {return cell_idx.row==0;};
        bool _solid_dborder_front(grid_idx_t cell_idx) {return cell_idx.depth==DEPTH_NUM - 1;};
        bool _solid_dborder_back(grid_idx_t cell_idx) {return cell_idx.depth==0;};
        
        bool _valid_cell(grid_idx_t cell_idx) {
            bool liquid_coord = (cell_idx.col >= 0 && cell_idx.row >= 0 && cell_idx.col < COL_NUM && cell_idx.row < ROW_NUM);
            bool liquid = liquid_coord && cell_idx.depth >= 0 && cell_idx.depth < DEPTH_NUM; 
            return liquid;
        }

        // valid horizontal velocity cell
        bool _valid_hcell(grid_idx_t cell_idx) {
            bool liquid_coord = (cell_idx.col >= 0 && cell_idx.row >= 0 && cell_idx.col <= COL_NUM && cell_idx.row < ROW_NUM);
            bool liquid = liquid_coord && cell_idx.depth >= 0 && cell_idx.depth < DEPTH_NUM; 
            return liquid;
        }

        // valid vertical velocity cell
        bool _valid_vcell(grid_idx_t cell_idx) {
            bool liquid_coord = (cell_idx.col >= 0 && cell_idx.row >= 0 && cell_idx.col < COL_NUM && cell_idx.row <= ROW_NUM);
            bool liquid = liquid_coord && cell_idx.depth >= 0 && cell_idx.depth < DEPTH_NUM; 
            return liquid;
        }

        // valid depth velocity cell
        bool _valid_dcell(grid_idx_t cell_idx) {
            bool liquid_coord = (cell_idx.col >= 0 && cell_idx.row >= 0 && cell_idx.col < COL_NUM && cell_idx.row < ROW_NUM);
            bool liquid = liquid_coord && cell_idx.depth >= 0 && cell_idx.depth <= DEPTH_NUM; 
            return liquid;
        }

        bool _fluid_cell(grid_idx_t cell_idx){
            // fluid and non-zero density (to avoid nans)
            return _valid_cell(cell_idx) && this->cells[cell_idx.row][cell_idx.col][cell_idx.depth].type==FLUID;
        }

        bool _full_fluid_cell(grid_idx_t cell_idx){
            return _fluid_cell(cell_idx) && this->cells[cell_idx.row][cell_idx.col][cell_idx.depth].density>0;
        }

        bool _air_cell(grid_idx_t cell_idx){
            return _valid_cell(cell_idx) && (!this->_full_fluid_cell(cell_idx));
        }

        grid_idx_t _lneighbor(grid_idx_t cell_idx){return {cell_idx.row, cell_idx.col - 1, cell_idx.depth};}
        grid_idx_t _rneighbor(grid_idx_t cell_idx){return {cell_idx.row, cell_idx.col + 1, cell_idx.depth};}
        grid_idx_t _tneighbor(grid_idx_t cell_idx){return {cell_idx.row + 1, cell_idx.col, cell_idx.depth};}
        grid_idx_t _bneighbor(grid_idx_t cell_idx){return {cell_idx.row - 1, cell_idx.col, cell_idx.depth};}
        grid_idx_t _fneighbor(grid_idx_t cell_idx){return {cell_idx.row, cell_idx.col, cell_idx.depth - 1};}
        grid_idx_t _baneighbor(grid_idx_t cell_idx){return {cell_idx.row, cell_idx.col, cell_idx.depth + 1};}

        int get_flat_idx(grid_idx_t idx) {return idx.col + COL_NUM * idx.row + ROW_NUM * COL_NUM * idx.depth;}

        /*
        Helper functions for solving pressure equations
        */
        void set_dV_idx(grid_idx_t cell, double entry){
            this->dV[this->get_flat_idx(cell)] = entry;
        };

        void add_to_dV(grid_idx_t cell, double entry){
            this->dV[this->get_flat_idx(cell)] += entry;
        };

        void set_A_idx(grid_idx_t cell, Direction dir, double entry){
            // printf("setting %d to ", dir);
            // print_cell(cell);
            // printf("%d size %d\n", this->get_flat_idx(cell), this->sparseA[this->get_flat_idx(cell)].size());
            this->sparseA[this->get_flat_idx(cell)][dir] = entry;
        };

        void add_to_A(grid_idx_t cell, Direction dir, double entry){
            // print_cell(cell);
            this->sparseA[this->get_flat_idx(cell)][dir] += entry;
        };

        void reset(){
            // printf("resetting solver...\n");
            for(int idx = 0; idx < this->dV.size(); idx++)
                this->dV[idx] = 0.0;

            for(int idx = 0; idx < this->pVec.size(); idx++)
                this->pVec[idx] = 0.0;
        

            for(int idx = 0; idx < this->diagE.size(); idx++)
                this->diagE[idx] = 0.0;
            
            for(std::vector<double>& vec : this->sparseA){
                for(int idx = 0; idx < vec.size(); idx++)
                    vec[idx] = 0.0;
            }
        };

        // pressure components
        std::vector<double> dV;
        // WARNING -> right now, because we only store sparseA as a NUM_CELLS x 5 length matrix, its not spatially correct
        std::vector<std::vector<double>> sparseA;
        std::vector<double> pVec;
        std::vector<double> diagE;

        // methods for solving PCG
        std::vector<double> apply_preconditioner(std::vector<double> res);
        std::vector<double> apply_A(std::vector<double> search);
        void build_preconditioner();
        bool solve_with_PCG();
};