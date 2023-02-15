#include <assemble_stiffness.h>

//Assemble the global stiffness matrix for the finite element mesh.

//Input:
//  q - generalized coordinates for the FEM system
//  qdot - generalized velocity for the FEM system
//  V - the nx3 matrix of undeformed vertex positions. Each row is a single undeformed vertex position.
//  T - the mx4 tetrahedron connectivity matrix. Each row contains to indices into V that indicate a spring between those vertices.
//  v0 - the mx1 vector of undeformed tetrahedron volumes
//  C,D - material parameters for the Neo-Hookean model
//Output:
//  K - the sparse, global stiffness matrix
void assemble_stiffness(Eigen::SparseMatrixd &K, Eigen::Ref<const Eigen::VectorXd> q, Eigen::Ref<const Eigen::VectorXd> qdot, 
                     Eigen::Ref<const Eigen::MatrixXd> V, Eigen::Ref<const Eigen::MatrixXi> T, Eigen::Ref<const Eigen::VectorXd> v0, 
                     double C, double D) { 
        //d2V_linear_tetrahedron_dq2

    std::vector<Eigen::Triplet<double>> trip_list;
    int n_tet = T.rows();
    for (unsigned idx = 0; idx < n_tet; idx++)
    {
 	Eigen::Matrix1212d temp_H; 
	temp_H.setZero();
	int v1 = T(idx, 0), v2 = T(idx, 1), v3 = T(idx, 2), v4 = T(idx, 3);
	//(Eigen::Matrix1212d &H, Eigen::Ref<const Eigen::VectorXd> q,  Eigen::Ref<const Eigen::MatrixXd> V, Eigen::Ref<const Eigen::RowVectorXi> element, double volume, double C, double D)
	d2V_linear_tetrahedron_dq2(temp_H, q, V, Eigen::RowVector4i(v1, v2, v3, v4), v0(idx), C, D); 

	//TODO:set H to corresponding K
	//v1v1	
	trip_list.push_back(Eigen::Triplet<double>(v1 * 3    , v1 * 3    , - temp_H(0 , 0)));
	trip_list.push_back(Eigen::Triplet<double>(v1 * 3    , v1 * 3 + 1, - temp_H(0 , 1)));
	trip_list.push_back(Eigen::Triplet<double>(v1 * 3    , v1 * 3 + 2, - temp_H(0 , 2)));
	trip_list.push_back(Eigen::Triplet<double>(v1 * 3 + 1, v1 * 3    , - temp_H(1 , 0)));
	trip_list.push_back(Eigen::Triplet<double>(v1 * 3 + 1, v1 * 3 + 1, - temp_H(1 , 1)));
	trip_list.push_back(Eigen::Triplet<double>(v1 * 3 + 1, v1 * 3 + 2, - temp_H(1 , 2)));
	trip_list.push_back(Eigen::Triplet<double>(v1 * 3 + 2, v1 * 3    , - temp_H(2 , 0)));
	trip_list.push_back(Eigen::Triplet<double>(v1 * 3 + 2, v1 * 3 + 1, - temp_H(2 , 1)));
	trip_list.push_back(Eigen::Triplet<double>(v1 * 3 + 2, v1 * 3 + 2, - temp_H(2 , 2)));

	//v1v2	
	trip_list.push_back(Eigen::Triplet<double>(v1 * 3    , v2 * 3    , - temp_H(0 , 3)));
	trip_list.push_back(Eigen::Triplet<double>(v1 * 3    , v2 * 3 + 1, - temp_H(0 , 4)));
	trip_list.push_back(Eigen::Triplet<double>(v1 * 3    , v2 * 3 + 2, - temp_H(0 , 5)));
	trip_list.push_back(Eigen::Triplet<double>(v1 * 3 + 1, v2 * 3    , - temp_H(1 , 3)));
	trip_list.push_back(Eigen::Triplet<double>(v1 * 3 + 1, v2 * 3 + 1, - temp_H(1 , 4)));
	trip_list.push_back(Eigen::Triplet<double>(v1 * 3 + 1, v2 * 3 + 2, - temp_H(1 , 5)));
	trip_list.push_back(Eigen::Triplet<double>(v1 * 3 + 2, v2 * 3    , - temp_H(2 , 3)));
	trip_list.push_back(Eigen::Triplet<double>(v1 * 3 + 2, v2 * 3 + 1, - temp_H(2 , 4)));
	trip_list.push_back(Eigen::Triplet<double>(v1 * 3 + 2, v2 * 3 + 2, - temp_H(2 , 5)));

	//v1v3	
	trip_list.push_back(Eigen::Triplet<double>(v1 * 3    , v3 * 3    , - temp_H(0 , 6)));
	trip_list.push_back(Eigen::Triplet<double>(v1 * 3    , v3 * 3 + 1, - temp_H(0 , 7)));
	trip_list.push_back(Eigen::Triplet<double>(v1 * 3    , v3 * 3 + 2, - temp_H(0 , 8)));
	trip_list.push_back(Eigen::Triplet<double>(v1 * 3 + 1, v3 * 3    , - temp_H(1 , 6)));
	trip_list.push_back(Eigen::Triplet<double>(v1 * 3 + 1, v3 * 3 + 1, - temp_H(1 , 7)));
	trip_list.push_back(Eigen::Triplet<double>(v1 * 3 + 1, v3 * 3 + 2, - temp_H(1 , 8)));
	trip_list.push_back(Eigen::Triplet<double>(v1 * 3 + 2, v3 * 3    , - temp_H(2 , 6)));
	trip_list.push_back(Eigen::Triplet<double>(v1 * 3 + 2, v3 * 3 + 1, - temp_H(2 , 7)));
	trip_list.push_back(Eigen::Triplet<double>(v1 * 3 + 2, v3 * 3 + 2, - temp_H(2 , 8)));

	//v1v4	
	trip_list.push_back(Eigen::Triplet<double>(v1 * 3    , v4 * 3    , - temp_H(0 , 9)));
	trip_list.push_back(Eigen::Triplet<double>(v1 * 3    , v4 * 3 + 1, - temp_H(0 , 10)));
	trip_list.push_back(Eigen::Triplet<double>(v1 * 3    , v4 * 3 + 2, - temp_H(0 , 11)));
	trip_list.push_back(Eigen::Triplet<double>(v1 * 3 + 1, v4 * 3    , - temp_H(1 , 9)));
	trip_list.push_back(Eigen::Triplet<double>(v1 * 3 + 1, v4 * 3 + 1, - temp_H(1 , 10)));
	trip_list.push_back(Eigen::Triplet<double>(v1 * 3 + 1, v4 * 3 + 2, - temp_H(1 , 11)));
	trip_list.push_back(Eigen::Triplet<double>(v1 * 3 + 2, v4 * 3    , - temp_H(2 , 9)));
	trip_list.push_back(Eigen::Triplet<double>(v1 * 3 + 2, v4 * 3 + 1, - temp_H(2 , 10)));
	trip_list.push_back(Eigen::Triplet<double>(v1 * 3 + 2, v4 * 3 + 2, - temp_H(2 , 11)));

	//v2v1	
	trip_list.push_back(Eigen::Triplet<double>(v2 * 3    , v1 * 3    , - temp_H(3 , 0)));
	trip_list.push_back(Eigen::Triplet<double>(v2 * 3    , v1 * 3 + 1, - temp_H(3 , 1)));
	trip_list.push_back(Eigen::Triplet<double>(v2 * 3    , v1 * 3 + 2, - temp_H(3 , 2)));
	trip_list.push_back(Eigen::Triplet<double>(v2 * 3 + 1, v1 * 3    , - temp_H(4 , 0)));
	trip_list.push_back(Eigen::Triplet<double>(v2 * 3 + 1, v1 * 3 + 1, - temp_H(4 , 1)));
	trip_list.push_back(Eigen::Triplet<double>(v2 * 3 + 1, v1 * 3 + 2, - temp_H(4 , 2)));
	trip_list.push_back(Eigen::Triplet<double>(v2 * 3 + 2, v1 * 3    , - temp_H(5 , 0)));
	trip_list.push_back(Eigen::Triplet<double>(v2 * 3 + 2, v1 * 3 + 1, - temp_H(5 , 1)));
	trip_list.push_back(Eigen::Triplet<double>(v2 * 3 + 2, v1 * 3 + 2, - temp_H(5 , 2)));

	//v2v2	
	trip_list.push_back(Eigen::Triplet<double>(v2 * 3    , v2 * 3    , - temp_H(3 , 3)));
	trip_list.push_back(Eigen::Triplet<double>(v2 * 3    , v2 * 3 + 1, - temp_H(3 , 4)));
	trip_list.push_back(Eigen::Triplet<double>(v2 * 3    , v2 * 3 + 2, - temp_H(3 , 5)));
	trip_list.push_back(Eigen::Triplet<double>(v2 * 3 + 1, v2 * 3    , - temp_H(4 , 3)));
	trip_list.push_back(Eigen::Triplet<double>(v2 * 3 + 1, v2 * 3 + 1, - temp_H(4 , 4)));
	trip_list.push_back(Eigen::Triplet<double>(v2 * 3 + 1, v2 * 3 + 2, - temp_H(4 , 5)));
	trip_list.push_back(Eigen::Triplet<double>(v2 * 3 + 2, v2 * 3    , - temp_H(5 , 3)));
	trip_list.push_back(Eigen::Triplet<double>(v2 * 3 + 2, v2 * 3 + 1, - temp_H(5 , 4)));
	trip_list.push_back(Eigen::Triplet<double>(v2 * 3 + 2, v2 * 3 + 2, - temp_H(5 , 5)));

	//v2v3	
	trip_list.push_back(Eigen::Triplet<double>(v2 * 3    , v3 * 3    , - temp_H(3 , 6)));
	trip_list.push_back(Eigen::Triplet<double>(v2 * 3    , v3 * 3 + 1, - temp_H(3 , 7)));
	trip_list.push_back(Eigen::Triplet<double>(v2 * 3    , v3 * 3 + 2, - temp_H(3 , 8)));
	trip_list.push_back(Eigen::Triplet<double>(v2 * 3 + 1, v3 * 3    , - temp_H(4 , 6)));
	trip_list.push_back(Eigen::Triplet<double>(v2 * 3 + 1, v3 * 3 + 1, - temp_H(4 , 7)));
	trip_list.push_back(Eigen::Triplet<double>(v2 * 3 + 1, v3 * 3 + 2, - temp_H(4 , 8)));
	trip_list.push_back(Eigen::Triplet<double>(v2 * 3 + 2, v3 * 3    , - temp_H(5 , 6)));
	trip_list.push_back(Eigen::Triplet<double>(v2 * 3 + 2, v3 * 3 + 1, - temp_H(5 , 7)));
	trip_list.push_back(Eigen::Triplet<double>(v2 * 3 + 2, v3 * 3 + 2, - temp_H(5 , 8)));

	//v2v4	
	trip_list.push_back(Eigen::Triplet<double>(v2 * 3    , v4 * 3    , - temp_H(3 , 9)));
	trip_list.push_back(Eigen::Triplet<double>(v2 * 3    , v4 * 3 + 1, - temp_H(3 , 10)));
	trip_list.push_back(Eigen::Triplet<double>(v2 * 3    , v4 * 3 + 2, - temp_H(3 , 11)));
	trip_list.push_back(Eigen::Triplet<double>(v2 * 3 + 1, v4 * 3    , - temp_H(4 , 9)));
	trip_list.push_back(Eigen::Triplet<double>(v2 * 3 + 1, v4 * 3 + 1, - temp_H(4 , 10)));
	trip_list.push_back(Eigen::Triplet<double>(v2 * 3 + 1, v4 * 3 + 2, - temp_H(4 , 11)));
	trip_list.push_back(Eigen::Triplet<double>(v2 * 3 + 2, v4 * 3    , - temp_H(5 , 9)));
	trip_list.push_back(Eigen::Triplet<double>(v2 * 3 + 2, v4 * 3 + 1, - temp_H(5 , 10)));
	trip_list.push_back(Eigen::Triplet<double>(v2 * 3 + 2, v4 * 3 + 2, - temp_H(5 , 11)));

	//v3v1	
	trip_list.push_back(Eigen::Triplet<double>(v3 * 3    , v1 * 3    , - temp_H(6 , 0)));
	trip_list.push_back(Eigen::Triplet<double>(v3 * 3    , v1 * 3 + 1, - temp_H(6 , 1)));
	trip_list.push_back(Eigen::Triplet<double>(v3 * 3    , v1 * 3 + 2, - temp_H(6 , 2)));
	trip_list.push_back(Eigen::Triplet<double>(v3 * 3 + 1, v1 * 3    , - temp_H(7 , 0)));
	trip_list.push_back(Eigen::Triplet<double>(v3 * 3 + 1, v1 * 3 + 1, - temp_H(7 , 1)));
	trip_list.push_back(Eigen::Triplet<double>(v3 * 3 + 1, v1 * 3 + 2, - temp_H(7 , 2)));
	trip_list.push_back(Eigen::Triplet<double>(v3 * 3 + 2, v1 * 3    , - temp_H(8 , 0)));
	trip_list.push_back(Eigen::Triplet<double>(v3 * 3 + 2, v1 * 3 + 1, - temp_H(8 , 1)));
	trip_list.push_back(Eigen::Triplet<double>(v3 * 3 + 2, v1 * 3 + 2, - temp_H(8 , 2)));

	//v3v2	
	trip_list.push_back(Eigen::Triplet<double>(v3 * 3    , v2 * 3    , - temp_H(6 , 3)));
	trip_list.push_back(Eigen::Triplet<double>(v3 * 3    , v2 * 3 + 1, - temp_H(6 , 4)));
	trip_list.push_back(Eigen::Triplet<double>(v3 * 3    , v2 * 3 + 2, - temp_H(6 , 5)));
	trip_list.push_back(Eigen::Triplet<double>(v3 * 3 + 1, v2 * 3    , - temp_H(7 , 3)));
	trip_list.push_back(Eigen::Triplet<double>(v3 * 3 + 1, v2 * 3 + 1, - temp_H(7 , 4)));
	trip_list.push_back(Eigen::Triplet<double>(v3 * 3 + 1, v2 * 3 + 2, - temp_H(7 , 5)));
	trip_list.push_back(Eigen::Triplet<double>(v3 * 3 + 2, v2 * 3    , - temp_H(8 , 3)));
	trip_list.push_back(Eigen::Triplet<double>(v3 * 3 + 2, v2 * 3 + 1, - temp_H(8 , 4)));
	trip_list.push_back(Eigen::Triplet<double>(v3 * 3 + 2, v2 * 3 + 2, - temp_H(8 , 5)));

	//v3v3	
	trip_list.push_back(Eigen::Triplet<double>(v3 * 3    , v3 * 3    , - temp_H(6 , 6)));
	trip_list.push_back(Eigen::Triplet<double>(v3 * 3    , v3 * 3 + 1, - temp_H(6 , 7)));
	trip_list.push_back(Eigen::Triplet<double>(v3 * 3    , v3 * 3 + 2, - temp_H(6 , 8)));
	trip_list.push_back(Eigen::Triplet<double>(v3 * 3 + 1, v3 * 3    , - temp_H(7 , 6)));
	trip_list.push_back(Eigen::Triplet<double>(v3 * 3 + 1, v3 * 3 + 1, - temp_H(7 , 7)));
	trip_list.push_back(Eigen::Triplet<double>(v3 * 3 + 1, v3 * 3 + 2, - temp_H(7 , 8)));
	trip_list.push_back(Eigen::Triplet<double>(v3 * 3 + 2, v3 * 3    , - temp_H(8 , 6)));
	trip_list.push_back(Eigen::Triplet<double>(v3 * 3 + 2, v3 * 3 + 1, - temp_H(8 , 7)));
	trip_list.push_back(Eigen::Triplet<double>(v3 * 3 + 2, v3 * 3 + 2, - temp_H(8 , 8)));

	//v3v4	
	trip_list.push_back(Eigen::Triplet<double>(v3 * 3    , v4 * 3    , - temp_H(6 , 9)));
	trip_list.push_back(Eigen::Triplet<double>(v3 * 3    , v4 * 3 + 1, - temp_H(6 , 10)));
	trip_list.push_back(Eigen::Triplet<double>(v3 * 3    , v4 * 3 + 2, - temp_H(6 , 11)));
	trip_list.push_back(Eigen::Triplet<double>(v3 * 3 + 1, v4 * 3    , - temp_H(7 , 9)));
	trip_list.push_back(Eigen::Triplet<double>(v3 * 3 + 1, v4 * 3 + 1, - temp_H(7 , 10)));
	trip_list.push_back(Eigen::Triplet<double>(v3 * 3 + 1, v4 * 3 + 2, - temp_H(7 , 11)));
	trip_list.push_back(Eigen::Triplet<double>(v3 * 3 + 2, v4 * 3    , - temp_H(8 , 9)));
	trip_list.push_back(Eigen::Triplet<double>(v3 * 3 + 2, v4 * 3 + 1, - temp_H(8 , 10)));
	trip_list.push_back(Eigen::Triplet<double>(v3 * 3 + 2, v4 * 3 + 2, - temp_H(8 , 11)));

	//v4v1	
	trip_list.push_back(Eigen::Triplet<double>(v4 * 3    , v1 * 3    , - temp_H(9 , 0)));
	trip_list.push_back(Eigen::Triplet<double>(v4 * 3    , v1 * 3 + 1, - temp_H(9 , 1)));
	trip_list.push_back(Eigen::Triplet<double>(v4 * 3    , v1 * 3 + 2, - temp_H(9 , 2)));
	trip_list.push_back(Eigen::Triplet<double>(v4 * 3 + 1, v1 * 3    , - temp_H(10 , 0)));
	trip_list.push_back(Eigen::Triplet<double>(v4 * 3 + 1, v1 * 3 + 1, - temp_H(10 , 1)));
	trip_list.push_back(Eigen::Triplet<double>(v4 * 3 + 1, v1 * 3 + 2, - temp_H(10 , 2)));
	trip_list.push_back(Eigen::Triplet<double>(v4 * 3 + 2, v1 * 3    , - temp_H(11 , 0)));
	trip_list.push_back(Eigen::Triplet<double>(v4 * 3 + 2, v1 * 3 + 1, - temp_H(11 , 1)));
	trip_list.push_back(Eigen::Triplet<double>(v4 * 3 + 2, v1 * 3 + 2, - temp_H(11 , 2)));

	//v4v2	
	trip_list.push_back(Eigen::Triplet<double>(v4 * 3    , v2 * 3    , - temp_H(9 , 3)));
	trip_list.push_back(Eigen::Triplet<double>(v4 * 3    , v2 * 3 + 1, - temp_H(9 , 4)));
	trip_list.push_back(Eigen::Triplet<double>(v4 * 3    , v2 * 3 + 2, - temp_H(9 , 5)));
	trip_list.push_back(Eigen::Triplet<double>(v4 * 3 + 1, v2 * 3    , - temp_H(10 , 3)));
	trip_list.push_back(Eigen::Triplet<double>(v4 * 3 + 1, v2 * 3 + 1, - temp_H(10 , 4)));
	trip_list.push_back(Eigen::Triplet<double>(v4 * 3 + 1, v2 * 3 + 2, - temp_H(10 , 5)));
	trip_list.push_back(Eigen::Triplet<double>(v4 * 3 + 2, v2 * 3    , - temp_H(11 , 3)));
	trip_list.push_back(Eigen::Triplet<double>(v4 * 3 + 2, v2 * 3 + 1, - temp_H(11 , 4)));
	trip_list.push_back(Eigen::Triplet<double>(v4 * 3 + 2, v2 * 3 + 2, - temp_H(11 , 5)));

	//v4v3	
	trip_list.push_back(Eigen::Triplet<double>(v4 * 3    , v3 * 3    , - temp_H(9 , 6)));
	trip_list.push_back(Eigen::Triplet<double>(v4 * 3    , v3 * 3 + 1, - temp_H(9 , 7)));
	trip_list.push_back(Eigen::Triplet<double>(v4 * 3    , v3 * 3 + 2, - temp_H(9 , 8)));
	trip_list.push_back(Eigen::Triplet<double>(v4 * 3 + 1, v3 * 3    , - temp_H(10 , 6)));
	trip_list.push_back(Eigen::Triplet<double>(v4 * 3 + 1, v3 * 3 + 1, - temp_H(10 , 7)));
	trip_list.push_back(Eigen::Triplet<double>(v4 * 3 + 1, v3 * 3 + 2, - temp_H(10 , 8)));
	trip_list.push_back(Eigen::Triplet<double>(v4 * 3 + 2, v3 * 3    , - temp_H(11 , 6)));
	trip_list.push_back(Eigen::Triplet<double>(v4 * 3 + 2, v3 * 3 + 1, - temp_H(11 , 7)));
	trip_list.push_back(Eigen::Triplet<double>(v4 * 3 + 2, v3 * 3 + 2, - temp_H(11 , 8)));

	//v4v4	
	trip_list.push_back(Eigen::Triplet<double>(v4 * 3    , v4 * 3    , - temp_H(9 , 9)));
	trip_list.push_back(Eigen::Triplet<double>(v4 * 3    , v4 * 3 + 1, - temp_H(9 , 10)));
	trip_list.push_back(Eigen::Triplet<double>(v4 * 3    , v4 * 3 + 2, - temp_H(9 , 11)));
	trip_list.push_back(Eigen::Triplet<double>(v4 * 3 + 1, v4 * 3    , - temp_H(10 , 9)));
	trip_list.push_back(Eigen::Triplet<double>(v4 * 3 + 1, v4 * 3 + 1, - temp_H(10 , 10)));
	trip_list.push_back(Eigen::Triplet<double>(v4 * 3 + 1, v4 * 3 + 2, - temp_H(10 , 11)));
	trip_list.push_back(Eigen::Triplet<double>(v4 * 3 + 2, v4 * 3    , - temp_H(11 , 9)));
	trip_list.push_back(Eigen::Triplet<double>(v4 * 3 + 2, v4 * 3 + 1, - temp_H(11 , 10)));
	trip_list.push_back(Eigen::Triplet<double>(v4 * 3 + 2, v4 * 3 + 2, - temp_H(11 , 11)));

    }
    K.resize(q.size(), q.size());
    K.setFromTriplets(trip_list.begin(), trip_list.end());
        
};
