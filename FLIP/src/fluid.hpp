#include <math.h>
#include "grid.hpp"

#define NUM_PARTICLES 1024
#define SQRT_NUM_PARTICLES 32

struct Point
{
    public:
        Point(double x, double y, double z)
        {
            this->x = x;
            this->y = y;
            this->z = z;
        };

        Point()
        {
            this->x = 0;
            this->y = 0;
            this->z = 0;
        };

        double x;
        double y;
        double z;
};

struct Particle
{
    public:
        Particle(double x, double y, double z)
        {
            // TODO make sure this is memory safe in terms of ownership
            this->velocity = Point();
            this->position = Point(x, y, z);
        };

        Particle()
        {
            // TODO make sure this is memory safe
            this->velocity = Point();
            this->position = Point();
        };

        Point velocity;
        Point position;
};

class Fluid
{
    public:
        Fluid(double width, double height)
        {
            this->grid = new Grid(width, height);

            // initialize particles uniformly across grid
            for (int row_idx=0; row_idx<SQRT_NUM_PARTICLES; row_idx++)
            {
                for (int col_idx=0; col_idx<SQRT_NUM_PARTICLES; col_idx++)
                {
                    double x = width*row_idx/SQRT_NUM_PARTICLES;
                    double y = width*col_idx/SQRT_NUM_PARTICLES;
                    particles[row_idx*SQRT_NUM_PARTICLES+col_idx] = Particle(x, y, 0.0f);
                }
            }
        };

        ~Fluid() 
        {
            delete this->grid;
        };

        void timestep();

    private:
        Grid *grid;

        // TODO make sure this is memory safe in terms of ownership
        Particle particles[NUM_PARTICLES];

        void advection();
        void external_forces();
};