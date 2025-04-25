#include <math.h>
#include <chrono>
#include <opencv2/opencv.hpp>
#include <opencv2/viz.hpp>

#include "grid.hpp"

#define PARTICLE_MASS 1.0

static Point GRAVITY = Point(0.0f, -9.8f, 0.0f);

// class Profiler;

class Fluid
{
    public:
        Fluid(double width, double height, double depth);
        ~Fluid();

        void run();
        void timestep();
        void display_particles();
        void print_particles();

    private:
        Grid *grid;
        Profiler *profiler;
        
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