#include "grid.hpp"

void Grid::build_preconditioner(){
    // TODO MAKE IT SO THAT AIR CELLS CAN EXIST
    for (int depth_idx = 0; depth_idx < DEPTH_NUM; depth_idx ++)
    {   
    for(int row_idx = 0; row_idx < ROW_NUM; row_idx++){
        for(int col_idx = 0; col_idx < COL_NUM; col_idx++)
        {
            if(_fluid_cell({row_idx, col_idx, depth_idx})){

            int idx = this->get_flat_idx({row_idx, col_idx, depth_idx});

            grid_idx_t lneigh = this->_lneighbor({row_idx, col_idx, depth_idx});
            grid_idx_t tneigh = this->_tneighbor({row_idx, col_idx, depth_idx});

            double root = sparseA[idx][CENTER];

            // TODO -> replace air and solid cells w 0
            double t1,t2,t3,t4 = 0.0;
            if(this->_fluid_cell(lneigh))
            {
                int lidx = this->get_flat_idx(lneigh);
                double t2 = (sparseA[lidx][RIGHT]/diagE[lidx]) * (sparseA[lidx][RIGHT]/diagE[lidx]);
                double t4 = (sparseA[lidx][RIGHT]/diagE[lidx]) * (sparseA[lidx][BOTTOM]/diagE[lidx]);
            }

            if(this->_fluid_cell(tneigh))
            {
                int tidx = this->get_flat_idx(tneigh);
                double t1 = (sparseA[tidx][BOTTOM]/diagE[tidx]) * (sparseA[tidx][BOTTOM]/diagE[tidx]);
                double t3 = (sparseA[tidx][BOTTOM]/diagE[tidx]) * (sparseA[tidx][RIGHT]/diagE[tidx]);
            }

            this->diagE[idx] = std::sqrt(root - t1 - t2 - t3 - t4);

            }
        }
    }
    }
}

std::vector<double> Grid::apply_preconditioner(std::vector<double> res)
{
    /*
    Calculate q for Lq = res
     L = FE^-1 + E 
     F is lower triangular of A, E is diagonal preconditioner
     then can perform Gaussian Elimination:

        for each cell_idx:
            r = (FE^-1 + E) @ q 
            -> q[idx] = r - (FE^-1 + E)[idx] @ q
                = r - F @ (E^-1 @ q) - (diagE[cell_idx] * q[cell_idx] (recursive), so divide)
                = ( r - A[cell_idx][LEFT] * (q[lneighbor] / diagE[lneighbor]) 
                    - A[cell_idx][TOP] * (q[tneighbor] / diagE[tneighbor]) 
                    ) / diagE[cell_idx]
    */ 
   std::vector<double> q_est(ROW_NUM*COL_NUM*DEPTH_NUM);
   for (int depth_idx = 0; depth_idx < DEPTH_NUM; depth_idx ++)
   {   
    for(int row_idx = 0; row_idx < ROW_NUM; row_idx++){
        for(int col_idx = 0; col_idx < COL_NUM; col_idx++)
        {
            // find corresponding indices
            int cell_idx = this->get_flat_idx({row_idx, col_idx, depth_idx});

            grid_idx_t lneigh = this->_lneighbor({row_idx, col_idx, depth_idx});
            grid_idx_t tneigh = this->_tneighbor({row_idx, col_idx, depth_idx});

            double left, top = 0.0;
            if(this->_valid_cell(lneigh))
            {
                int lidx = this->get_flat_idx(lneigh);
                left = sparseA[cell_idx][LEFT] * q_est[lidx] / diagE[lidx];
            }

            if(this->_valid_cell(tneigh))
            {
                int tidx = this->get_flat_idx(tneigh);
                top = sparseA[cell_idx][TOP] * q_est[tidx] / diagE[tidx];
            }
            
            
            q_est[cell_idx] = (res[cell_idx] - left - top) / diagE[cell_idx];
        }
    }
    }

    // for (double el : q_est) 
    //     printf("%f ", el);
    // printf("\n");

    /*
    Calculate z for L^Tz = q -> Az \approx L(L^T)z = res
    L^T = F^TE^-1 + E
        for each cell_idx:
        q = (F^TE^-1 + E) @ z
        -> q[idx] = r - (F^TE^-1 + E)[idx] @ z
                = r - F^T @ (E^-1 @ q) - (diagE[cell_idx] * z[cell_idx] (recursive), so divide)
                = (r - A[cell_idx][RIGHT] * (q[rneighbor] / diagE[rneighbor]) 
                    - A[cell_idx][BOTTOM] * (q[bneighbor] / diagE[bneighbor]) 
                    ) / diagE[cell_idx]
    */
    std::vector<double> z_est(ROW_NUM*COL_NUM*DEPTH_NUM);
    for (int depth_idx = 0; depth_idx < DEPTH_NUM; depth_idx ++)
    {   
    for(int row_idx = ROW_NUM - 1; row_idx >= 0; row_idx--){
        for(int col_idx = COL_NUM - 1; col_idx >= 0; col_idx--)
        {
            // find corresponding indices
            int cell_idx = this->get_flat_idx({row_idx, col_idx, depth_idx});

            grid_idx_t rneigh = this->_rneighbor({row_idx, col_idx, depth_idx});
            grid_idx_t bneigh = this->_bneighbor({row_idx, col_idx, depth_idx});

            double right, bottom = 0.0;
            if(this->_valid_cell(rneigh))
            {
                int ridx = this->get_flat_idx(rneigh);
                right = sparseA[cell_idx][RIGHT] * z_est[ridx] / diagE[ridx];
            }

            if(this->_valid_cell(bneigh))
            {
                int bidx = this->get_flat_idx(bneigh);
                bottom = sparseA[cell_idx][BOTTOM] * z_est[bidx] / diagE[bidx];
            }
            
            
            z_est[cell_idx] = (q_est[cell_idx] - bottom - right) / diagE[cell_idx];
        }
    }
    }

    return z_est;
}

// 
std::vector<double> Grid::apply_A(std::vector<double> search)
{
    // printf("A:\n");
    std::vector<double> result(this->sparseA.size());
    for (int depth_idx = 0; depth_idx < DEPTH_NUM; depth_idx ++)
    {   
    for(int row_idx = 0; row_idx < ROW_NUM; row_idx++){
        for(int col_idx = 0; col_idx < COL_NUM; col_idx++)
        {
            // multiply the of A with the corresponding pij values
            if(this->_fluid_cell({row_idx, col_idx, depth_idx}))
            {
                int cell_idx = this->get_flat_idx({row_idx, col_idx, depth_idx});

                // iterate through 5 values and apply A
                #ifdef DEBUG
                printf("(%d, %d): %f, %f, %f, %f, %f\n", row_idx, col_idx,
                    this->sparseA[cell_idx][LEFT], 
                    this->sparseA[cell_idx][RIGHT], this->sparseA[cell_idx][CENTER], this->sparseA[cell_idx][TOP],
                    this->sparseA[cell_idx][BOTTOM]
                );

                printf("%f\n",
                    search[cell_idx]
                );
                #endif
                
                result[cell_idx] += this->sparseA[cell_idx][CENTER] * search[cell_idx];

                grid_idx_t lneigh = this->_lneighbor({row_idx, col_idx, depth_idx});
                grid_idx_t tneigh = this->_tneighbor({row_idx, col_idx, depth_idx});

                if(this->_valid_cell(lneigh))
                {
                    int lidx = this->get_flat_idx(lneigh);
                    result[cell_idx] += this->sparseA[cell_idx][LEFT] * search[lidx];
                }
    
                if(this->_valid_cell(tneigh))
                {
                    int tidx = this->get_flat_idx(tneigh);
                    result[cell_idx] += this->sparseA[cell_idx][TOP] * search[tidx];
                }

                grid_idx_t rneigh = this->_rneighbor({row_idx, col_idx, depth_idx});
                grid_idx_t bneigh = this->_bneighbor({row_idx, col_idx, depth_idx});
    
                if(this->_valid_cell(rneigh))
                {
                    int ridx = this->get_flat_idx(rneigh);
                    result[cell_idx] += this->sparseA[cell_idx][RIGHT] * search[ridx];
                }
    
                if(this->_valid_cell(bneigh))
                {
                    int bidx = this->get_flat_idx(bneigh);
                    result[cell_idx] += this->sparseA[cell_idx][BOTTOM] * search[bidx];
                }

                grid_idx_t fneigh = this->_rneighbor({row_idx, col_idx, depth_idx});
                grid_idx_t baneigh = this->_bneighbor({row_idx, col_idx, depth_idx});
    
                if(this->_valid_cell(fneigh))
                {
                    int fidx = this->get_flat_idx(fneigh);
                    result[cell_idx] += this->sparseA[cell_idx][FRONT] * search[fidx];
                }
    
                if(this->_valid_cell(baneigh))
                {
                    int baidx = this->get_flat_idx(baneigh);
                    result[cell_idx] += this->sparseA[cell_idx][BACK] * search[baidx];
                }

                // printf("result: %f\n", result[cell_idx]);
            }
        }
    }
    }

    return result;
}

/*
    Uses preconditioned conjugate gradient method to solve for the pressure vector p
    Returns true if success, false otherwise
*/
bool Grid::solve_with_PCG()
{
    // construct preconditioner
    this->build_preconditioner();

    // printf("built preconditioner\n");

    // for (double el : this->diagE) 
    //     printf("%f ", el);
    // printf("\n");

    // solve Ap = d * V with PCG
    std::vector<double> residual = this->dV;

    // printf("residual\n");

    // for (double el : residual) 
    //     printf("%f ", el);
    // printf("\n");


    // check convergence
    bool converged = true;
    for (double el : residual) 
        if(std::abs(el) > TOL) {
            converged = false;
            break;
        }
    if(converged) return true;

    // initialize necessary vectors
    std::vector<double> auxilary_vec;
    std::vector<double> search_vec = residual;

    // printf("applied preconditioner\n");

    // for (double el : auxilary_vec) 
    //     printf("%f ", el);
    // printf("\n");

    // learning rate
    // double sigma = std::inner_product(residual.begin(), residual.end(), auxilary_vec.begin(), 0.0);
    double sigma = std::inner_product(residual.begin(), residual.end(), residual.begin(), 0.0);

    // idk what any of this math does
    for(int iter = 0; iter < MAX_ITER; iter++)
    {
        auxilary_vec = this->apply_A(search_vec);

        // printf("applied A\n");

        // printf("res\n");
        // for (double el : residual) 
        //     printf("%f ", el);
        // printf("\n");
    
        // printf("auxilary\n");
        // for (double el : auxilary_vec) 
        //     printf("%f ", el);
        // printf("\n");

        // printf("search\n");
        // for (double el : search_vec) 
        //     printf("%f ", el);
        // printf("\n");

        // printf("pressure\n");
        // for (double el : pVec) 
        //     printf("%f ", el);
        // printf("\n");

        // how different the residual and search vec are
        double s_dot_a = std::inner_product(search_vec.begin(), search_vec.end(), auxilary_vec.begin(), 0.0);

        // printf("s_dot_a: %f\n");

        // A is singular! everything is ok i think
        if(s_dot_a<0.0000001) return true;

        double alpha = sigma / s_dot_a;

        // update estimates
        for (int idx = 0; idx < this->pVec.size(); idx++)
        {
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
                // break;
            }

        // printf("MAX VALUE: %f\n", max);
        // if(converged) return true;

        // printf("res\n");
        // for (double el : residual) 
        //     printf("%f ", el);
        // printf("\n");

        // auxilary_vec = this->apply_preconditioner(residual);
        double sigma_new = std::inner_product(residual.begin(), residual.end(), residual.begin(), 0.0);

        // printf("applied preconditioner\n");

        // update search vector
        for (int idx = 0; idx < this->pVec.size(); idx++)
            search_vec[idx] = residual[idx] + (sigma_new / sigma) * search_vec[idx];
        
        
        sigma = sigma_new;
    }

    return false;
}