#include <assemble_stiffness.h>
#include <iostream>

void assemble_stiffness(Eigen::SparseMatrixd &K, Eigen::Ref<const Eigen::VectorXd> q, Eigen::Ref<const Eigen::VectorXd> qdot, 
                     Eigen::Ref<const Eigen::MatrixXd> V, Eigen::Ref<const Eigen::MatrixXi> E, Eigen::Ref<const Eigen::VectorXd> l0, 
                     double k) { 
//Input:
//  q - generalized coordinates for the mass-spring system
//  qdot - generalized velocity for the mass spring system
//  V - the nx3 matrix of undeformed vertex positions. Each row is a single undeformed vertex position.
//  E - the mx2 spring connectivity matrix. Each row contains to indices into V that indicate a spring between those vertices.
//  l0 - the mx1 vector of undeformed spring lengths
//  k - the stiffness of each spring in the mass-spring system
//Output:
//  K - the 3nx3n sparse stiffness matrix which is the negative hessian of the potential energy function.       
    
    //should use d2V_spring_particle_particle_dq2.cpp to assemble hessian 
    /* Hessian:
     * A | B *
     * C | D */
    //connection: E (v1, v2)
    
    //Block Location:
    // A - (v1, v1)
    // B - (v1, v2)
    // C - (v2, v1)
    // D - (v2, v2)
    
    std::vector<Eigen::Triplet<double>> trip_list;
    for (unsigned idx = 0; idx < E.rows(); idx++)
    {
    	int v1 = E(idx, 0), v2 = E(idx, 1);
	//TODO:declare var
        Eigen::Vector3d q0(q(3*v1), q(3*v1+1), q(3*v1+2));
	Eigen::Vector3d q1(q(3*v2), q(3*v2+1), q(3*v2+2));
	Eigen::Matrix66d H;
	d2V_spring_particle_particle_dq2(H, q0, q1, l0(idx), k); 

	//TODO:set H to corresponding K
	//A	
	trip_list.push_back(Eigen::Triplet<double>(v1 * 3    , v1 * 3    , - H(0 , 0)));
	trip_list.push_back(Eigen::Triplet<double>(v1 * 3    , v1 * 3 + 1, - H(0 , 1)));
	trip_list.push_back(Eigen::Triplet<double>(v1 * 3    , v1 * 3 + 2, - H(0 , 2)));
	trip_list.push_back(Eigen::Triplet<double>(v1 * 3 + 1, v1 * 3    , - H(1 , 0)));
	trip_list.push_back(Eigen::Triplet<double>(v1 * 3 + 1, v1 * 3 + 1, - H(1 , 1)));
	trip_list.push_back(Eigen::Triplet<double>(v1 * 3 + 1, v1 * 3 + 2, - H(1 , 2)));
	trip_list.push_back(Eigen::Triplet<double>(v1 * 3 + 2, v1 * 3    , - H(2 , 0)));
	trip_list.push_back(Eigen::Triplet<double>(v1 * 3 + 2, v1 * 3 + 1, - H(2 , 1)));
	trip_list.push_back(Eigen::Triplet<double>(v1 * 3 + 2, v1 * 3 + 2, - H(2 , 2)));

	//B	
	trip_list.push_back(Eigen::Triplet<double>(v1 * 3    , v2 * 3    , - H(0 , 3)));
	trip_list.push_back(Eigen::Triplet<double>(v1 * 3    , v2 * 3 + 1, - H(0 , 4)));
	trip_list.push_back(Eigen::Triplet<double>(v1 * 3    , v2 * 3 + 2, - H(0 , 5)));
	trip_list.push_back(Eigen::Triplet<double>(v1 * 3 + 1, v2 * 3    , - H(1 , 3)));
	trip_list.push_back(Eigen::Triplet<double>(v1 * 3 + 1, v2 * 3 + 1, - H(1 , 4)));
	trip_list.push_back(Eigen::Triplet<double>(v1 * 3 + 1, v2 * 3 + 2, - H(1 , 5)));
	trip_list.push_back(Eigen::Triplet<double>(v1 * 3 + 2, v2 * 3    , - H(2 , 3)));
	trip_list.push_back(Eigen::Triplet<double>(v1 * 3 + 2, v2 * 3 + 1, - H(2 , 4)));
	trip_list.push_back(Eigen::Triplet<double>(v1 * 3 + 2, v2 * 3 + 2, - H(2 , 5)));

	//C	
	trip_list.push_back(Eigen::Triplet<double>(v2 * 3    , v1 * 3    , - H(3 , 0)));
	trip_list.push_back(Eigen::Triplet<double>(v2 * 3    , v1 * 3 + 1, - H(3 , 1)));
	trip_list.push_back(Eigen::Triplet<double>(v2 * 3    , v1 * 3 + 2, - H(3 , 2)));
	trip_list.push_back(Eigen::Triplet<double>(v2 * 3 + 1, v1 * 3    , - H(4 , 0)));
	trip_list.push_back(Eigen::Triplet<double>(v2 * 3 + 1, v1 * 3 + 1, - H(4 , 1)));
	trip_list.push_back(Eigen::Triplet<double>(v2 * 3 + 1, v1 * 3 + 2, - H(4 , 2)));
	trip_list.push_back(Eigen::Triplet<double>(v2 * 3 + 2, v1 * 3    , - H(5 , 0)));
	trip_list.push_back(Eigen::Triplet<double>(v2 * 3 + 2, v1 * 3 + 1, - H(5 , 1)));
	trip_list.push_back(Eigen::Triplet<double>(v2 * 3 + 2, v1 * 3 + 2, - H(5 , 2)));

	//D	
	trip_list.push_back(Eigen::Triplet<double>(v2 * 3    , v2 * 3    , - H(3 , 3)));
	trip_list.push_back(Eigen::Triplet<double>(v2 * 3    , v2 * 3 + 1, - H(3 , 4)));
	trip_list.push_back(Eigen::Triplet<double>(v2 * 3    , v2 * 3 + 2, - H(3 , 5)));
	trip_list.push_back(Eigen::Triplet<double>(v2 * 3 + 1, v2 * 3    , - H(4 , 3)));
	trip_list.push_back(Eigen::Triplet<double>(v2 * 3 + 1, v2 * 3 + 1, - H(4 , 4)));
	trip_list.push_back(Eigen::Triplet<double>(v2 * 3 + 1, v2 * 3 + 2, - H(4 , 5)));
	trip_list.push_back(Eigen::Triplet<double>(v2 * 3 + 2, v2 * 3    , - H(5 , 3)));
	trip_list.push_back(Eigen::Triplet<double>(v2 * 3 + 2, v2 * 3 + 1, - H(5 , 4)));
	trip_list.push_back(Eigen::Triplet<double>(v2 * 3 + 2, v2 * 3 + 2, - H(5 , 5)));
    }
    K.resize(q.size(), q.size());
    K.setFromTriplets(trip_list.begin(), trip_list.end());

};
