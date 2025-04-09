#include <vector>
#include <tuple>
#include "particle.hpp"

// Grid coarsness
#define ROW_NUM 500
#define COL_NUM 500

typedef std::tuple<size_t, size_t> grid_idx_t;

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
        };

        void transfer_to_grid(std::vector<Particle>& particles);
        void transfer_from_grid(std::vector<Particle>& particles);
        void solve_pressure();

        grid_idx_t get_particle_idx(Particle p);

        double pressure_grid[ROW_NUM][COL_NUM] {};
        double horizontal_velocity[ROW_NUM][COL_NUM] {};
        double vertical_velocity[ROW_NUM][COL_NUM] {};

        double width;
        double height;

        double cell_width;
        double cell_height;
};