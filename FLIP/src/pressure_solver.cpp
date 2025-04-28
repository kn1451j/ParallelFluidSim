#include "grid.hpp"
#include <cmath>
#include <cassert>


bool Grid::solve_with_PCG()
{
    // solve Ap = d * V with PCG
    std::vector<double> residual = this->dV;
    // initialize necessary vectors
    std::vector<double> auxilary_vec = residual;
    std::vector<double> search_vec = auxilary_vec;

    // check convergence
    bool converged = true;
    // double max = 0.0;
    #pragma omp parallel for reduction(&& : converged) schedule(static)
    for (double el : residual) 
        if(std::abs(el) > TOL) {
            converged &= false;
        }

    // learning rate
    double sigma = 0;
    #pragma omp parallel for reduction(+ : sigma) schedule(static)
    for(int idx = 0; idx < this->pVec.size(); idx++)
    {
        sigma += residual[idx] * residual[idx];
    }

    for(int iter = 0; iter < MAX_ITER; iter++)
    {

    #pragma omp parallel for schedule(dynamic)
    for (int depth_idx = 0; depth_idx < DEPTH_NUM; depth_idx ++)
    {   
        for(int row_idx = 0; row_idx < ROW_NUM; row_idx++){
            for(int col_idx = 0; col_idx < COL_NUM; col_idx++)
            {
            int cell_idx = this->get_flat_idx({row_idx, col_idx, depth_idx});

            // multiply the of A with the corresponding pij values
            if(this->_fluid_cell({row_idx, col_idx, depth_idx}))
            {   
                auxilary_vec[cell_idx] = this->sparseA[cell_idx][CENTER] * search_vec[cell_idx];

                grid_idx_t lneigh = this->_lneighbor({row_idx, col_idx, depth_idx});
                grid_idx_t tneigh = this->_tneighbor({row_idx, col_idx, depth_idx});

                if(this->_valid_cell(lneigh))
                {
                    int lidx = this->get_flat_idx(lneigh);
                    auxilary_vec[cell_idx] += this->sparseA[cell_idx][LEFT] * search_vec[lidx];
                }
    
                if(this->_valid_cell(tneigh))
                {
                    int tidx = this->get_flat_idx(tneigh);
                    auxilary_vec[cell_idx] += this->sparseA[cell_idx][TOP] * search_vec[tidx];
                }

                grid_idx_t rneigh = this->_rneighbor({row_idx, col_idx, depth_idx});
                grid_idx_t bneigh = this->_bneighbor({row_idx, col_idx, depth_idx});
    
                if(this->_valid_cell(rneigh))
                {
                    int ridx = this->get_flat_idx(rneigh);
                    auxilary_vec[cell_idx] += this->sparseA[cell_idx][RIGHT] * search_vec[ridx];
                }
    
                if(this->_valid_cell(bneigh))
                {
                    int bidx = this->get_flat_idx(bneigh);
                    auxilary_vec[cell_idx] += this->sparseA[cell_idx][BOTTOM] * search_vec[bidx];
                }

                grid_idx_t fneigh = this->_fneighbor({row_idx, col_idx, depth_idx});
                grid_idx_t baneigh = this->_baneighbor({row_idx, col_idx, depth_idx});
    
                if(this->_valid_cell(fneigh))
                {
                    int fidx = this->get_flat_idx(fneigh);
                    auxilary_vec[cell_idx] += this->sparseA[cell_idx][FRONT] * search_vec[fidx];
                }
    
                if(this->_valid_cell(baneigh))
                {
                    int baidx = this->get_flat_idx(baneigh);
                    auxilary_vec[cell_idx] += this->sparseA[cell_idx][BACK] * search_vec[baidx];
                }
            }
            else{
                auxilary_vec[cell_idx] = 0;
            }
        }
        }
        }

        // how different the residual and search vec are
        double s_dot_a = 0;
        #pragma omp parallel for reduction(+ : s_dot_a) schedule(static)
        for(int idx = 0; idx < this->pVec.size(); idx++)
        {
            s_dot_a += search_vec[idx] * auxilary_vec[idx];
        }

        // printf("s_dot_a: %f\n");

        // A is singular! everything is ok i think
        if(abs(s_dot_a)<EPS) return true;

        double alpha = sigma / s_dot_a;

        // update estimates
        #pragma omp parallel for schedule(static, 8) // cache line 64 bytes, these are doubles, should be 4?
        for (int idx = 0; idx < this->pVec.size(); idx++)
        {
            // printf("%f ", auxilary_vec[idx]);
            this->pVec[idx] += alpha * search_vec[idx];
            residual[idx] -= alpha * auxilary_vec[idx];
        }

        // check convergence
        bool converged = true;
        // double max = 0.0;
        #pragma omp parallel for reduction(&& : converged) schedule(static)
        for (double el : residual) 
            if(std::abs(el) > TOL) {
                converged &= false;
            }

        // printf("MAX VALUE: %f\n", max);
        if(converged){ 
            // printf("\n");
            return true;}

        double sigma_new = 0;
        #pragma omp parallel for reduction(+ : sigma_new) schedule(static)
        for(int idx = 0; idx < this->pVec.size(); idx++)
        {
            sigma_new += residual[idx] * residual[idx];
        }

        // update search vector
        #pragma omp parallel for schedule(static, 8)
        for (int idx = 0; idx < this->pVec.size(); idx++)
            search_vec[idx] = residual[idx] + (sigma_new / sigma) * search_vec[idx];
        
        
        sigma = sigma_new;

        // ensure everyone has the same copy of the neighboring search vec (will because of implicit barrier)
    }

    return false;
}