#include <vector>
#include <tuple>
#include <set>

#include "particle.hpp"
#include "sparse_matrix.hpp"

// Grid coarsness
#define ROW_NUM 4
#define COL_NUM 4

typedef std::tuple<size_t, size_t> grid_idx_t;

// typedef Eigen::Matrix< double, 4, 5 > Matrix45d;

enum CELL_TYPE {FLUID, SOLID, GAS};

struct Cell
{
    CELL_TYPE type = GAS;
};

struct Vertex
{
    double value;
    void zero () {this->value = 0;};
    void set(double value) {this->value = value;};

    // used for computing particle to grid interpolation
    double normalization;
    std::set<size_t> neighbor_particles {};
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
        Grid(double width, double height){
            this->width = width;
            this->height = height;
            this->cell_width = width/COL_NUM;
            this->cell_height = height/ROW_NUM;

            this->A = new SparseMatrix(ROW_NUM, COL_NUM);
        };

        // exclusive owners of A
        ~Grid() {delete this->A;};

        void transfer_to_grid(std::vector<Particle>& particles);
        void transfer_from_grid(std::vector<Particle>& particles);
        void solve_pressure();

        grid_idx_t get_particle_idx(Particle p);

        // change for this for cache optimiality
        Vertex pressure_grid[ROW_NUM][COL_NUM] {};
        Vertex horizontal_velocity[ROW_NUM][2*COL_NUM] {};
        Vertex vertical_velocity[2*ROW_NUM][COL_NUM] {};
        Cell cells[ROW_NUM][COL_NUM] {};

        double width;
        double height;

        double cell_width;
        double cell_height;

        void print_grid();

    private:
        SparseMatrix *A;
        
};