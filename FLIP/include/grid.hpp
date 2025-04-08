#include <vector>
#include "particle.hpp"

// Grid coarsness
#define ROW_NUM 500
#define COL_NUM 500

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
        };

        void transfer_to_grid(std::vector<Particle>& particles);
        void transfer_from_grid(std::vector<Particle>& particles);
        void solve_pressure();

    private:
        double pressure_grid[ROW_NUM][COL_NUM] {};
        double horizontal_velocity[ROW_NUM][COL_NUM] {};
        double vertical_velocity[ROW_NUM][COL_NUM] {};

        double width;
        double height;
};