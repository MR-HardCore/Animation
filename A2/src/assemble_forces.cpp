#include <assemble_forces.h>
#include <iostream>


void assemble_forces(Eigen::VectorXd &f, Eigen::Ref<const Eigen::VectorXd> q, Eigen::Ref<const Eigen::VectorXd> qdot, 
                     Eigen::Ref<const Eigen::MatrixXd> V, Eigen::Ref<const Eigen::MatrixXi> E, Eigen::Ref<const Eigen::VectorXd> l0, 
                     double mass, double k) {
//Input:
//  q - generalized coordinates for the mass-spring system
//  qdot - generalized velocity for the mass spring system
//  V - the nx3 matrix of undeformed vertex positions. Each row is a single undeformed vertex position.
//  E - the mx2 spring connectivity matrix. Each row contains to indices into V that indicate a spring between those vertices.
//  l0 - the mx1 vector of undeformed spring lengths
//  m - the mass of each particle in the mass-spring system
//  k - the stiffness of each spring in the mass-spring system
//Output:
//  f - the vector 3xn vector of forces acting on each node of the mass-spring system 



    //should use dV_gravity_particle_dq and dV_spring_particle_particle_dq
    /* Jacobian:
     * A *
     * B */
    //connection: E (v1, v2)
    
    //Block Location:
    // A - v1
    // B - v2
    
    f.resize(q.size());
    f.setZero();
    //for spring force 
    for (unsigned idx = 0; idx < E.rows(); idx++)
    {
	int v1 = E(idx, 0), v2 = E(idx, 1);
	//TODO:make variable work
        Eigen::Vector3d q0(q(3*v1), q(3*v1+1), q(3*v1+2));
	Eigen::Vector3d q1(q(3*v2), q(3*v2+1), q(3*v2+2));
	Eigen::Vector6d f_sp;
	dV_spring_particle_particle_dq(f_sp, q0, q1, l0(idx), k);

        //TODO:set f to proper position
	f(v1 * 3    ) -= f_sp(0);
	f(v1 * 3 + 1) -= f_sp(1);
	f(v1 * 3 + 2) -= f_sp(2);
	f(v2 * 3    ) -= f_sp(3);
	f(v2 * 3 + 1) -= f_sp(4);
	f(v2 * 3 + 2) -= f_sp(5);
	
    }

};
