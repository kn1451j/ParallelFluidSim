#include "fluid.hpp"

void Fluid::run()
{
    //initializes the display
    this->display_particles();

    // start timestep timer
    this->prev_timestep = std::chrono::high_resolution_clock::now();

    while(1){
        this->timestep();

        // printf("Updated %zu particles\n", this->particles.size());

        // display a frame (for now just in 2D)
        this->display_particles();

        #ifndef DEBUG
        if(cv::waitKey(1)>0) break;
        #endif
        #ifdef DEBUG
        cv::waitKey(0); // frame by frame for debugging
        #endif
    }
}

void Fluid::timestep()
{
    // set dt to just be a constant for debugging
    #ifndef REALTIME 
    double dt = 0.5;
    #endif
    #ifdef REALTIME
    auto now = std::chrono::high_resolution_clock::now();
    double dt = std::chrono::duration<double>(now - this->prev_timestep).count();
    // update the time
    this->prev_timestep = now;
    #endif

    // iterate through
    this->advection(dt);
    this->external_forces(dt);

    // this->print_particles();

    this->grid->transfer_to_grid(this->particles);

    // this->grid->print_grid();

    // TODO -> add ghost pressures
    this->grid->solve_pressure(dt);

    this->grid->transfer_from_grid(this->particles);

    #ifdef DEBUG
    // this->print_particles();
    // this->grid->print_grid();
    #endif

}

/*
    For now, each particle just moves in its velocity dir by timestep
    Simple first-order update of the particle's position
*/
void Fluid::advection(double timestep)
{
    for(Particle &p : this->particles)
    {
        // advect
        p.advect(timestep);

        // clamp particle position to screen (TODO magic nums)
        p.position.x = std::clamp(p.position.x, CLAMP, this->width-CLAMP);
        p.position.y = std::clamp(p.position.y, CLAMP, this->height-CLAMP);
        p.position.z = std::clamp(p.position.z, CLAMP, this->depth-CLAMP);
        // add z
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
    // std::vector<cv::Point3f> points;
    // // Fill the 'points' matrix with your 3D data
    // for(Particle p : this->particles)
    // {
    //     // ensure inverting the height
    //     points.push_back(cv::Point3f(p.position.x, this->height - p.position.y, p.position.z));
    //     cv::circle(this->display, cv::Point(p.position.x, this->height - p.position.y), 5, cv::Scalar(0, 0, 255), 1, cv::FILLED);
    // }

    // cv::viz::WCloud cloud_widget(points);
    // this->display.showWidget("point_cloud", cloud_widget);

    this->display = cv::Mat(this->width, this->height, CV_8UC4, cv::Scalar(0, 0, 0));
    // for(size_t depth_idx = 0; depth_idx<DEPTH_NUM; depth_idx++)
    // {
    //     for(size_t row_idx = 0; row_idx<ROW_NUM; row_idx++)
    //     {
    //         for(size_t col_idx = 0; col_idx<COL_NUM; col_idx++)
    //         {
    //             Cell cell = this->grid->cells[row_idx][col_idx][depth_idx];

    //             double color = (depth - cell.position.z) / depth;
    //             if(cell.density>0)
    //                 cv::rectangle(this->display, cv::Point(cell.position.x - this->grid->cell_width/2, this->height - (cell.position.y - this->grid->cell_height/2)), 
    //                 cv::Point(cell.position.x + this->grid->cell_width/2, this->height - (cell.position.y + this->grid->cell_height/2)), 
    //                 cv::Scalar(255, 255*color, 255*color), cv::FILLED);
    //         }
    //     }
    // }

    for(Particle p : this->particles)
    {
        // ensure inverting the height
        double color = (depth - p.position.z) / depth;
        cv::circle(this->display, cv::Point(p.position.x, this->height - p.position.y), 5, cv::Scalar(255*color, 0, 255*(1 - color), 1), 1, cv::FILLED);
    }

    cv::imshow("Fluid Sim Slice", this->display);
}

void Fluid::print_particles()
{
    for(Particle p : this->particles)
    {
        printf("Particle %d at (%f, %f, %f)\n", p.id, p.position.x, p.position.y, p.position.z);
        printf("uv: (%f, %f, %f)\n", p.velocity.x, p.velocity.y, p.velocity.z);
    }
}