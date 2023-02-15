#include <mass_matrix_mesh.h>
#include <iostream>
//Assemble the full mass matrix for the entire tetrahedral mesh.

//Input: 
//  q - generalized coordinates for the FEM system
//  V - the nx3 matrix of undeformed vertex positions
//  F - the mx3 matrix of triangle-vertex indices
//  density - the density of the cloth material
//  areas - the mx1 vector of undeformed triangle areas
//Output:
//  M - sparse mass matrix for the entire mesh
void mass_matrix_mesh(Eigen::SparseMatrixd &M, Eigen::Ref<const Eigen::VectorXd> q, 
                         Eigen::Ref<const Eigen::MatrixXd> V, Eigen::Ref<const Eigen::MatrixXi> F,
                         double density, Eigen::Ref<const Eigen::VectorXd> areas) {
    
    double size_F = F.rows(), size_q = q.size();
    std::vector<Eigen::Triplet<double>> trip_list;
    for (unsigned idx = 0; idx < size_F; idx++)
    {
        int v1 = F(idx, 0), v2 = F(idx, 1), v3 = F(idx, 2);
        double mass_fac = density * areas(idx) / 12.0;
        trip_list.push_back(Eigen::Triplet<double>(v1 * 3 + 0, v1 * 3 + 0, mass_fac * 2));
        trip_list.push_back(Eigen::Triplet<double>(v1 * 3 + 1, v1 * 3 + 1, mass_fac * 2));
        trip_list.push_back(Eigen::Triplet<double>(v1 * 3 + 2, v1 * 3 + 2, mass_fac * 2));
        
        trip_list.push_back(Eigen::Triplet<double>(v1 * 3 + 0, v2 * 3 + 0, mass_fac));
        trip_list.push_back(Eigen::Triplet<double>(v1 * 3 + 1, v2 * 3 + 1, mass_fac));
        trip_list.push_back(Eigen::Triplet<double>(v1 * 3 + 2, v2 * 3 + 2, mass_fac));

        trip_list.push_back(Eigen::Triplet<double>(v1 * 3 + 0, v3 * 3 + 0, mass_fac));
        trip_list.push_back(Eigen::Triplet<double>(v1 * 3 + 1, v3 * 3 + 1, mass_fac));
        trip_list.push_back(Eigen::Triplet<double>(v1 * 3 + 2, v3 * 3 + 2, mass_fac));

        trip_list.push_back(Eigen::Triplet<double>(v2 * 3 + 0, v1 * 3 + 0, mass_fac));
        trip_list.push_back(Eigen::Triplet<double>(v2 * 3 + 1, v1 * 3 + 1, mass_fac));
        trip_list.push_back(Eigen::Triplet<double>(v2 * 3 + 2, v1 * 3 + 2, mass_fac));
        
        trip_list.push_back(Eigen::Triplet<double>(v2 * 3 + 0, v2 * 3 + 0, mass_fac * 2));
        trip_list.push_back(Eigen::Triplet<double>(v2 * 3 + 1, v2 * 3 + 1, mass_fac * 2));
        trip_list.push_back(Eigen::Triplet<double>(v2 * 3 + 2, v2 * 3 + 2, mass_fac * 2));
        
        trip_list.push_back(Eigen::Triplet<double>(v2 * 3 + 0, v3 * 3 + 0, mass_fac));
        trip_list.push_back(Eigen::Triplet<double>(v2 * 3 + 1, v3 * 3 + 1, mass_fac));
        trip_list.push_back(Eigen::Triplet<double>(v2 * 3 + 2, v3 * 3 + 2, mass_fac));

        trip_list.push_back(Eigen::Triplet<double>(v3 * 3 + 0, v1 * 3 + 0, mass_fac));
        trip_list.push_back(Eigen::Triplet<double>(v3 * 3 + 1, v1 * 3 + 1, mass_fac));
        trip_list.push_back(Eigen::Triplet<double>(v3 * 3 + 2, v1 * 3 + 2, mass_fac));

        trip_list.push_back(Eigen::Triplet<double>(v3 * 3 + 0, v2 * 3 + 0, mass_fac));
        trip_list.push_back(Eigen::Triplet<double>(v3 * 3 + 1, v2 * 3 + 1, mass_fac));
        trip_list.push_back(Eigen::Triplet<double>(v3 * 3 + 2, v2 * 3 + 2, mass_fac));
        
        trip_list.push_back(Eigen::Triplet<double>(v3 * 3 + 0, v3 * 3 + 0, mass_fac * 2));
        trip_list.push_back(Eigen::Triplet<double>(v3 * 3 + 1, v3 * 3 + 1, mass_fac * 2));
        trip_list.push_back(Eigen::Triplet<double>(v3 * 3 + 2, v3 * 3 + 2, mass_fac * 2));
    }    
    M.resize(size_q, size_q);
    M.setFromTriplets(trip_list.begin(), trip_list.end());
    
}
 
