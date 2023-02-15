#include <fixed_point_constraints.h>
#include <algorithm>
#include <iostream>
void fixed_point_constraints(Eigen::SparseMatrixd &P, unsigned int q_size, const std::vector<unsigned int> indices) {

//Input:
//  q_size - total number of scalar generalized coordinaes (3 times number of vertices in the mesh)
//  indices - indices (row ids in V) for fixed vertices
//Output:
//  P - mxn sparse matrix which projects out fixed vertices

//TODO: P is the selection matrix for non-moving particles
// should be 3m * 3n but treat valid one as 3*3 identity matrix
// m identity matrix in total
    
    
    /*
    // for debug only
    std::cout<<"start to test: "<<std::endl;
    for (unsigned idx = 0; idx < indices.size(); idx ++)
    {
        std::cout << idx << "th entry: " <<indices[idx] <<std::endl;
    }
    */
    std::vector<Eigen::Triplet<double>> trip_list;
    double n_moving = q_size - indices.size() * 3, n_gap = 0, counts = 0, row = n_moving / 3, col = q_size / 3;
    for (unsigned row_idx = 0, col_idx = 0; row_idx < row && col_idx < col; row_idx++, col_idx++)
    {
	while (col_idx == indices[n_gap] && col_idx < col)
	{
	    n_gap++;
	    col_idx++;
	}

	trip_list.push_back(Eigen::Triplet<double>(row_idx*3    , col_idx*3    , 1.));
	trip_list.push_back(Eigen::Triplet<double>(row_idx*3 + 1, col_idx*3 + 1, 1.));
	trip_list.push_back(Eigen::Triplet<double>(row_idx*3 + 2, col_idx*3 + 2, 1.));
    }
    P.resize(n_moving, q_size);
    P.setFromTriplets(trip_list.begin(), trip_list.end());

}

