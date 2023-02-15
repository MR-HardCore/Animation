#include <mass_matrix_particles.h>
#include <iostream>

void mass_matrix_particles(Eigen::SparseMatrixd &M, Eigen::Ref<const Eigen::VectorXd> q, double mass) {

//Input:
//  q - generalized coordiantes for the mass-spring system
//  mass - the mass of each particle in the mass-spring system.
//Output:
//  M - sparse mass matrix for mass-spring system

//TODO: create a diagonal matrix: all values are m (double)
// size_q * size_q diagonal matrix with m

    double size_q = q.size();
    std::vector<Eigen::Triplet<double>> trip_list;
    for (unsigned idx = 0; idx < size_q; idx++){ trip_list.push_back(Eigen::Triplet<double>(idx, idx, mass)); }     
    M.resize(size_q, size_q);
    M.setFromTriplets(trip_list.begin(), trip_list.end());

}
