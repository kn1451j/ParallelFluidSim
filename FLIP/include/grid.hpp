#include <vector>
#include <utility>
#include <set>

#include "particle.hpp"
#include "sparse_matrix.hpp"

// Grid coarsness
#define ROW_NUM 4
#define COL_NUM 4

#define EPS 0.000001

#define PIC_WEIGHT 0.9

typedef std::pair<size_t, size_t> grid_idx_t;

// typedef Eigen::Matrix< double, 4, 5 > Matrix45d;

enum CELL_TYPE {FLUID, SOLID, GAS};
enum VELOCITY {HORIZONTAL, VERTICAL};

struct Cell
{
    CELL_TYPE type = GAS;
};

struct Vertex
{
    double value;
    void reset () {
        this->prev_value = this->value;
        this->value = 0; 
        this->normalization = 0;
    };

    void set(double value) {
        this->value = value;
    };

    // used for computing particle to grid interpolation
    double normalization;
    
    // used for flip velocity change
    double prev_value;

    Point position;
};

struct Neighbors
{
    VELOCITY type;
    grid_idx_t neighbors[4];
}

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
        Grid(double width, double height);

        // exclusive owners of A
        ~Grid() {delete this->A;};

        void transfer_to_grid(std::vector<Particle>& particles);
        void transfer_from_grid(std::vector<Particle>& particles);
        void solve_pressure();

        grid_idx_t get_grid_idx(Particle p);

        Neighbors get_horizontal_neighbors(Particle p);
        Neighbors get_vertical_neighbors(Particle p);

        // change for this for cache optimiality
        Vertex pressure_grid[ROW_NUM][COL_NUM] {};
        Vertex horizontal_velocity[ROW_NUM][COL_NUM+1] {};
        Vertex vertical_velocity[ROW_NUM+1][COL_NUM] {};
        Cell cells[ROW_NUM][COL_NUM] {};

        double width;
        double height;

        double cell_width;
        double cell_height;

        void print_grid();

    private:
        SparseMatrix *A;
        
};