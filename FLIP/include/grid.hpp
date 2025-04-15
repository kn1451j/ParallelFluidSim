#include <vector>
#include <utility>
#include <set>

#include "particle.hpp"
#include "sparse_matrix.hpp"

#define DEBUG 1
// #define REALTIME 1

// Grid coarsness
#define ROW_NUM 4
#define COL_NUM 4

#define EPS 0.000001

#define PIC_WEIGHT 0.9

// -1 can indicate that we're in "solid" region
typedef std::pair<int, int> grid_idx_t;

// typedef Eigen::Matrix< double, 4, 5 > Matrix45d;

enum CELL_TYPE {FLUID, SOLID, GAS};
enum VELOCITY {HORIZONTAL, VERTICAL, CELL};


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
    
    // order: top left, top right, bottom right, bottom left
    grid_idx_t neighbors[4];
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
        Grid(double width, double height);

        // exclusive owners of A
        ~Grid() {delete this->A;};

        void transfer_to_grid(std::vector<Particle>& particles);
        void transfer_from_grid(std::vector<Particle>& particles);
        void solve_pressure();

        grid_idx_t get_grid_idx(Particle p);

        Neighbors get_horizontal_neighbors(Particle p);
        Neighbors get_vertical_neighbors(Particle p);
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
        Vertex pressure_grid[ROW_NUM][COL_NUM] {};
        Vertex horizontal_velocity[ROW_NUM][COL_NUM+1] {};
        Vertex vertical_velocity[ROW_NUM+1][COL_NUM] {};
        Cell cells[ROW_NUM][COL_NUM] {};

        double width;
        double height;

        double cell_width;
        double cell_height;
        double cell_volume;

        void print_grid();

    private:
        SparseMatrix *A;
        
};