#include "fluid.hpp"

void Fluid::run()
{
    //initializes the display
    this->display_particles();

    // start timestep timer
    this->prev_timestep = std::chrono::high_resolution_clock::now();

    while(1){
        this->timestep();

        printf("Updated %zu particles\n", this->particles.size());

        // display a frame (for now just in 2D)
        this->display_particles();

        if(cv::waitKey(1)>0) break;
    }
}

void Fluid::timestep()
{
    auto now = std::chrono::high_resolution_clock::now();
    double dt = std::chrono::duration<double>(now - this->prev_timestep).count();

    // update the time
    this->prev_timestep = now;

    // iterate through
    this->advection(dt);
    this->external_forces(dt);
    this->grid->transfer_to_grid(this->particles);
    // TODO -> add ghost pressures
    this->grid->solve_pressure();
    this->grid->transfer_from_grid(this->particles);
}

/*
    For now, each particle just moves in its velocity dir by timestep
    Simple first-order update of the particle's position
*/
void Fluid::advection(double timestep)
{
    for(Particle &p : this->particles)
    {
        p.advect(timestep);
    }
}

/*
    Simple first-order update of graviational force
*/
void Fluid::external_forces(double timestep)
{
    for(Particle &p : this->particles)
    {
        // for now we just have a gravitational force
        p.apply_force(GRAVITY, timestep);
    }
}

/*
    Just shows a 2D front-facing slice of particles
*/
void Fluid::display_particles()
{
    this->display = cv::Mat(this->width, this->height, CV_8UC3, cv::Scalar(0, 0, 0));

    for(Particle p : this->particles)
    {
        // ensure inverting the height
        cv::circle(this->display, cv::Point(p.position.x, this->height - p.position.y), 1, cv::Scalar(0, 0, 255), 1);
    }

    cv::imshow("Fluid Sim Slice", this->display);
}