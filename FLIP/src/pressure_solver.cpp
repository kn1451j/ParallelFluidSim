// #include "grid.hpp"
// #include <cmath>
// #include <cassert>

// #include "profiler.hpp"

void Grid::apply_A(std::vector<double> search, std::vector<double> &dest)
{
    // printf("A:\n");
    for (int depth_idx = 0; depth_idx < DEPTH_NUM; depth_idx ++)
    {   
    for(int row_idx = 0; row_idx < ROW_NUM; row_idx++){
        for(int col_idx = 0; col_idx < COL_NUM; col_idx++)
        {
            int cell_idx = this->get_flat_idx({row_idx, col_idx, depth_idx});

            // multiply the of A with the corresponding pij values
            if(this->_fluid_cell({row_idx, col_idx, depth_idx}))
            {
                // iterate through 5 values and apply A
                // #if DEBUG
                // printf("(%d, %d, %d): %f, %f, %f, %f, %f, %f, %f\n", row_idx, col_idx,
                //     depth_idx,
                //     this->sparseA[cell_idx][LEFT], 
                //     this->sparseA[cell_idx][RIGHT], this->sparseA[cell_idx][CENTER], this->sparseA[cell_idx][TOP],
                //     this->sparseA[cell_idx][BOTTOM], this->sparseA[cell_idx][BACK], this->sparseA[cell_idx][FRONT]
                // );

                // printf("%f\n",
                //     search[cell_idx]
                // );
                // #endif
                
                dest[cell_idx] = this->sparseA[cell_idx][CENTER] * search[cell_idx];

                grid_idx_t lneigh = this->_lneighbor({row_idx, col_idx, depth_idx});
                grid_idx_t tneigh = this->_tneighbor({row_idx, col_idx, depth_idx});

                if(this->_valid_cell(lneigh))
                {
                    int lidx = this->get_flat_idx(lneigh);
                    dest[cell_idx] += this->sparseA[cell_idx][LEFT] * search[lidx];
                }
    
                if(this->_valid_cell(tneigh))
                {
                    int tidx = this->get_flat_idx(tneigh);
                    dest[cell_idx] += this->sparseA[cell_idx][TOP] * search[tidx];
                }

                grid_idx_t rneigh = this->_rneighbor({row_idx, col_idx, depth_idx});
                grid_idx_t bneigh = this->_bneighbor({row_idx, col_idx, depth_idx});
    
                if(this->_valid_cell(rneigh))
                {
                    int ridx = this->get_flat_idx(rneigh);
                    dest[cell_idx] += this->sparseA[cell_idx][RIGHT] * search[ridx];
                }
    
                if(this->_valid_cell(bneigh))
                {
                    int bidx = this->get_flat_idx(bneigh);
                    dest[cell_idx] += this->sparseA[cell_idx][BOTTOM] * search[bidx];
                }

                grid_idx_t fneigh = this->_fneighbor({row_idx, col_idx, depth_idx});
                grid_idx_t baneigh = this->_baneighbor({row_idx, col_idx, depth_idx});
    
                if(this->_valid_cell(fneigh))
                {
                    int fidx = this->get_flat_idx(fneigh);
                    dest[cell_idx] += this->sparseA[cell_idx][FRONT] * search[fidx];
                }
    
                if(this->_valid_cell(baneigh))
                {
                    int baidx = this->get_flat_idx(baneigh);
                    dest[cell_idx] += this->sparseA[cell_idx][BACK] * search[baidx];
                }

                assert(!std::isnan(dest[cell_idx]));
                // printf("result: %f\n", result[cell_idx]);
            }
            else{
                dest[cell_idx] = 0;
            }
        }
    }
    }
}

/*
    Uses preconditioned conjugate gradient method to solve for the pressure vector p
    Returns true if success, false otherwise
*/
bool Grid::solve_with_PCG()
{
    // solve Ap = d * V with PCG
    std::vector<double> residual = this->dV;

    // check convergence
    bool converged = true;
    for (double el : residual) 
        if(std::abs(el) > TOL) {
            converged = false;
            break;
        }
    if(converged) return true;

    // initialize necessary vectors
    std::vector<double> auxilary_vec = residual;
    std::vector<double> search_vec = auxilary_vec;

    // learning rate
    double sigma = std::inner_product(residual.begin(), residual.end(), residual.begin(), 0.0);

    for(int iter = 0; iter < MAX_ITER; iter++)
    {
        this->apply_A(search_vec, auxilary_vec);

        // how different the residual and search vec are
        double s_dot_a = std::inner_product(search_vec.begin(), search_vec.end(), auxilary_vec.begin(), 0.0);

        // printf("s_dot_a: %f\n");

        // A is singular! everything is ok i think
        if(abs(s_dot_a)<EPS) return true;

        double alpha = sigma / s_dot_a;

        // update estimates
        for (int idx = 0; idx < this->pVec.size(); idx++)
        {
            // printf("%f ", auxilary_vec[idx]);
            this->pVec[idx] += alpha * search_vec[idx];
            residual[idx] -= alpha * auxilary_vec[idx];
        }

        // check convergence
        converged = true;
        double max = 0.0;
        for (double el : residual) 
            if(std::abs(el) > TOL) {
                converged = false;
                max = std::max(max,std::abs(el));
                break;
            }

        // printf("MAX VALUE: %f\n", max);
        if(converged){ 
            // printf("\n");
            return true;}

        double sigma_new = std::inner_product(residual.begin(), residual.end(), residual.begin(), 0.0);

        // update search vector
        for (int idx = 0; idx < this->pVec.size(); idx++)
            search_vec[idx] = residual[idx] + (sigma_new / sigma) * search_vec[idx];
        
        
        sigma = sigma_new;
    }

    return false;
}

