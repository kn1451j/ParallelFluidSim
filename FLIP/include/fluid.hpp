#include <math.h>
#include <chrono>
#include <opencv2/opencv.hpp>

#include "grid.hpp"

#define PARTICLE_MASS 1.0

static Point GRAVITY = Point(0.0f, -9.8f, 0.0f);

class Fluid
{
    public:
        Fluid(double width, double height)
        : display()
        {
            this->particles.resize(NUM_PARTICLES);
            this->width = width;
            this->height = height;

            // initialize a grid
            this->grid = new Grid(width, height);

            // initialize particles uniformly across grid
            for (int row_idx=0; row_idx<SQRT_NUM_PARTICLES; row_idx++)
            {
                for (int col_idx=0; col_idx<SQRT_NUM_PARTICLES; col_idx++)
                {
                    double x = width*(col_idx+0.5)/SQRT_NUM_PARTICLES;
                    double y = height*(row_idx+0.5)/SQRT_NUM_PARTICLES;

                    Point pose = Point(x, y, 0.0f);
                    Point vel = Point();

                    size_t particle_idx = row_idx * SQRT_NUM_PARTICLES + col_idx;
                    this->particles[particle_idx] = Particle(particle_idx, pose, vel, PARTICLE_MASS);
                }
            }
        };

        ~Fluid(){
            delete this->grid;
        };

        void run();
        void timestep();
        void display_particles();
        void print_particles();

    private:
        Grid *grid;
        
        double width;
        double height;

        // TODO make sure this is memory safe in terms of ownership
        std::vector<Particle> particles;
        std::chrono::time_point<std::chrono::steady_clock> prev_timestep;

        void advection(double timestep);
        void external_forces(double timestep);

        // size_t get_particle_idx(int row_idx, int col_idx){};

        // for displaying a 2D slice
        cv::Mat display;
};