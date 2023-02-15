#include <mass_matrix_mesh.h>
#include <mass_matrix_linear_tetrahedron.h>
#include <iostream>

//Assemble the full mass matrix for the entire tetrahedral mesh.

//Input:
//  qdot - generalized velocity for the FEM system
//  T - the mx4 vertex indices for tet mesh
//  density - density of material
//  v0 - the undeformed tetrahedra volumes
//Output:
//  M - Sparse mass matrix for the whole mesh.

void mass_matrix_mesh(Eigen::SparseMatrixd &M, Eigen::Ref<const Eigen::VectorXd> qdot, Eigen::Ref<const Eigen::MatrixXi> T, double density, Eigen::Ref<const Eigen::VectorXd> v0) {

    int n_tet = T.rows();
    std::vector<Eigen::Triplet<double>> trip_list;
    for (unsigned idx = 0; idx < n_tet; idx++)
    {
	Eigen::Matrix1212d temp_M;
	temp_M.setZero();
	int v1 = T(idx, 0), v2 = T(idx, 1), v3 = T(idx, 2), v4 = T(idx, 3);
	mass_matrix_linear_tetrahedron(temp_M, qdot, Eigen::RowVector4i(v1, v2, v3, v4), density, v0(idx));
	// temp_m has 48 elements for 16 grids as
	/*  v1v1 v1v2 v1v3 v1v4 
	 *  v2v1 v2v2 v2v3 v2v4
	 *  v3v1 v3v2 v3v3 v3v4
	 *  v4v1 v4v2 v4v3 v4v4
	 */
	//each small grid is 3 * 3 diagonal matrix

	//v1v1
	trip_list.push_back(Eigen::Triplet<double>(3 * v1 + 0 , 3 * v1 + 0, temp_M(0 , 0)));
	trip_list.push_back(Eigen::Triplet<double>(3 * v1 + 1 , 3 * v1 + 1, temp_M(1 , 1)));
	trip_list.push_back(Eigen::Triplet<double>(3 * v1 + 2 , 3 * v1 + 2, temp_M(2 , 2)));

	//v1v2
	trip_list.push_back(Eigen::Triplet<double>(3 * v1 + 0 , 3 * v2 + 0, temp_M(0 , 3)));
	trip_list.push_back(Eigen::Triplet<double>(3 * v1 + 1 , 3 * v2 + 1, temp_M(1 , 4)));
	trip_list.push_back(Eigen::Triplet<double>(3 * v1 + 2 , 3 * v2 + 2, temp_M(2 , 5)));

	//v1v3
	trip_list.push_back(Eigen::Triplet<double>(3 * v1 + 0 , 3 * v3 + 0, temp_M(0 , 6)));
	trip_list.push_back(Eigen::Triplet<double>(3 * v1 + 1 , 3 * v3 + 1, temp_M(1 , 7)));
	trip_list.push_back(Eigen::Triplet<double>(3 * v1 + 2 , 3 * v3 + 2, temp_M(2 , 8)));

	//v1v4
	trip_list.push_back(Eigen::Triplet<double>(3 * v1 + 0 , 3 * v4 + 0, temp_M(0 , 9)));
	trip_list.push_back(Eigen::Triplet<double>(3 * v1 + 1 , 3 * v4 + 1, temp_M(1 , 10)));
	trip_list.push_back(Eigen::Triplet<double>(3 * v1 + 2 , 3 * v4 + 2, temp_M(2 , 11)));
//////////////////
	//v2v1
	trip_list.push_back(Eigen::Triplet<double>(3 * v2 + 0 , 3 * v1 + 0, temp_M(3 , 0)));
	trip_list.push_back(Eigen::Triplet<double>(3 * v2 + 1 , 3 * v1 + 1, temp_M(4 , 1)));
	trip_list.push_back(Eigen::Triplet<double>(3 * v2 + 2 , 3 * v1 + 2, temp_M(5 , 2)));

	//v2v2
	trip_list.push_back(Eigen::Triplet<double>(3 * v2 + 0 , 3 * v2 + 0, temp_M(3 , 3)));
	trip_list.push_back(Eigen::Triplet<double>(3 * v2 + 1 , 3 * v2 + 1, temp_M(4 , 4)));
	trip_list.push_back(Eigen::Triplet<double>(3 * v2 + 2 , 3 * v2 + 2, temp_M(5 , 5)));

	//v2v3
	trip_list.push_back(Eigen::Triplet<double>(3 * v2 + 0 , 3 * v3 + 0, temp_M(3 , 6)));
	trip_list.push_back(Eigen::Triplet<double>(3 * v2 + 1 , 3 * v3 + 1, temp_M(4 , 7)));
	trip_list.push_back(Eigen::Triplet<double>(3 * v2 + 2 , 3 * v3 + 2, temp_M(5 , 8)));

	//v2v4
	trip_list.push_back(Eigen::Triplet<double>(3 * v2 + 0 , 3 * v4 + 0, temp_M(3 , 9)));
	trip_list.push_back(Eigen::Triplet<double>(3 * v2 + 1 , 3 * v4 + 1, temp_M(4 , 10)));
	trip_list.push_back(Eigen::Triplet<double>(3 * v2 + 2 , 3 * v4 + 2, temp_M(5 , 11)));
////////////////////
	//v3v1
	trip_list.push_back(Eigen::Triplet<double>(3 * v3 + 0 , 3 * v1 + 0, temp_M(6 , 0)));
	trip_list.push_back(Eigen::Triplet<double>(3 * v3 + 1 , 3 * v1 + 1, temp_M(7 , 1)));
	trip_list.push_back(Eigen::Triplet<double>(3 * v3 + 2 , 3 * v1 + 2, temp_M(8 , 2)));

	//v3v2
	trip_list.push_back(Eigen::Triplet<double>(3 * v3 + 0 , 3 * v2 + 0, temp_M(6 , 3)));
	trip_list.push_back(Eigen::Triplet<double>(3 * v3 + 1 , 3 * v2 + 1, temp_M(7 , 4)));
	trip_list.push_back(Eigen::Triplet<double>(3 * v3 + 2 , 3 * v2 + 2, temp_M(8 , 5)));

	//v3v3
	trip_list.push_back(Eigen::Triplet<double>(3 * v3 + 0 , 3 * v3 + 0, temp_M(6 , 6)));
	trip_list.push_back(Eigen::Triplet<double>(3 * v3 + 1 , 3 * v3 + 1, temp_M(7 , 7)));
	trip_list.push_back(Eigen::Triplet<double>(3 * v3 + 2 , 3 * v3 + 2, temp_M(8 , 8)));

	//v3v4
	trip_list.push_back(Eigen::Triplet<double>(3 * v3 + 0 , 3 * v4 + 0, temp_M(6 , 9)));
	trip_list.push_back(Eigen::Triplet<double>(3 * v3 + 1 , 3 * v4 + 1, temp_M(7 , 10)));
	trip_list.push_back(Eigen::Triplet<double>(3 * v3 + 2 , 3 * v4 + 2, temp_M(8 , 11)));
//////////////////////
	//v4v1
	trip_list.push_back(Eigen::Triplet<double>(3 * v4 + 0 , 3 * v1 + 0, temp_M(9 , 0)));
	trip_list.push_back(Eigen::Triplet<double>(3 * v4 + 1 , 3 * v1 + 1, temp_M(10 , 1)));
	trip_list.push_back(Eigen::Triplet<double>(3 * v4 + 2 , 3 * v1 + 2, temp_M(11 , 2)));

	//v4v2
	trip_list.push_back(Eigen::Triplet<double>(3 * v4 + 0 , 3 * v2 + 0, temp_M(9 , 3)));
	trip_list.push_back(Eigen::Triplet<double>(3 * v4 + 1 , 3 * v2 + 1, temp_M(10 , 4)));
	trip_list.push_back(Eigen::Triplet<double>(3 * v4 + 2 , 3 * v2 + 2, temp_M(11 , 5)));

	//v4v3
	trip_list.push_back(Eigen::Triplet<double>(3 * v4 + 0 , 3 * v3 + 0, temp_M(9 , 6)));
	trip_list.push_back(Eigen::Triplet<double>(3 * v4 + 1 , 3 * v3 + 1, temp_M(10 , 7)));
	trip_list.push_back(Eigen::Triplet<double>(3 * v4 + 2 , 3 * v3 + 2, temp_M(11 , 8)));

	//v4v4
	trip_list.push_back(Eigen::Triplet<double>(3 * v4 + 0 , 3 * v4 + 0, temp_M(9 , 9)));
	trip_list.push_back(Eigen::Triplet<double>(3 * v4 + 1 , 3 * v4 + 1, temp_M(10 , 10)));
	trip_list.push_back(Eigen::Triplet<double>(3 * v4 + 2 , 3 * v4 + 2, temp_M(11 , 11)));

    }
    M.resize(qdot.size(), qdot.size());
    M.setFromTriplets(trip_list.begin(), trip_list.end());
}
