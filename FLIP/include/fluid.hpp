#include <math.h>
#include <chrono>
#include <opencv2/opencv.hpp>
#include <opencv2/viz.hpp>

#include "grid.hpp"

#define PARTICLE_MASS 1.0

static Point GRAVITY = Point(0.0f, -9.8f, 0.0f);

class Fluid
{
    public:
        Fluid(double width, double height, double depth)
        : display()
        {
            this->particles.resize(NUM_PARTICLES);
            this->width = width;
            this->height = height;
            this->depth = depth;
            double density = PARTICLE_MASS * NUM_PARTICLES / (width * height * depth);
            // double density = 1;

            // initialize a grid
            this->grid = new Grid(width, height, depth, density);

            // initialize particles randomly across grid
            for (int idx=0; idx<NUM_PARTICLES; idx++)
            {
                    float r = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
                    double x = width*r;

                    r = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
                    double y = height*r;

                    r = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
                    double z = depth*r;

                    Point pose = Point(x, y, z);
                    Point vel = Point();

                    this->particles[idx] = Particle(idx, pose, vel, PARTICLE_MASS);
            }

            // this->display.setWindowSize(cv::Size(500, 500));
            // this->display.setWindowPosition(cv::Point(150, 150));
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
        double depth;

        // TODO make sure this is memory safe in terms of ownership
        std::vector<Particle> particles;
        std::chrono::time_point<std::chrono::steady_clock> prev_timestep;

        void advection(double timestep);
        void external_forces(double timestep);

        // for displaying a 3D image
        // cv::viz::Viz3d display;
        cv::Mat display;
};