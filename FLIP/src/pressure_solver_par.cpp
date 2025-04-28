#include "grid.hpp"
#include <cmath>
#include <cassert>
#include <mpi.h>

/*
Each grid block owns a block of a grid:
Grid block i out of nproc will own the indices (i * width -> width - 1, j * width -> width - 1, k * width -> width - 1)
for width = nproc / i and i,j,k \in nproc
*/
class GridBlock
{
    public:
        GridBlock();
        // pressure components
        std::vector<double> dV;
        std::vector<std::vector<double>> sparseA;
        std::vector<double> pVec;
        std::vector<double> diagE;

        void reset(){
            // printf("resetting solver...\n");
            for(int idx = 0; idx < this->dV.size(); idx++)
                this->dV[idx] = 0.0;

            for(int idx = 0; idx < this->pVec.size(); idx++)
                this->pVec[idx] = 0.0;
        

            for(int idx = 0; idx < this->diagE.size(); idx++)
                this->diagE[idx] = 0.0;
            
            for(std::vector<double>& vec : this->sparseA){
                for(int idx = 0; idx < vec.size(); idx++)
                    vec[idx] = 0.0;
            }
        };
}


bool GridBlock::solve_with_PCG()
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

    for (int depth_idx = 0; depth_idx < DEPTH_NUM; depth_idx ++)
    {   
        for(int row_idx = 0; row_idx < ROW_NUM; row_idx++){
            for(int col_idx = 0; col_idx < COL_NUM; col_idx++)
            {
            int cell_idx = this->get_flat_idx({row_idx, col_idx, depth_idx});

            // multiply the of A with the corresponding pij values
            if(this->_fluid_cell({row_idx, col_idx, depth_idx}))
            {   
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
            }
            else{
                dest[cell_idx] = 0;
            }
        }
        }
        }

        // how different the residual and search vec are
        double s_dot_a = 0;
        for(int idx = 0; idx < GRID_SIZE; idx++)
        {
            s_dot_a += search_vec[idx] * auxilary_vec.size();
        }
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

        // ensure everyone has the same copy of the neighboring search vec
    }

    return false;
}