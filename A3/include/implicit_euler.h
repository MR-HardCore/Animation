#include <Eigen/Dense>
#include <EigenTypes.h>
#include <newtons_method.h>
#include <iostream>

//Input:
//  q - generalized coordinates for the FEM system
//  qdot - generalized velocity for the FEM system
//  dt - the time step in seconds
//  mass - the mass matrix
//  energy(q, qdot) -  a function that computes the energy of the FEM system. This takes q and qdot as parameters, returns the energy value.
//  force(f, q, qdot) - a function that computes the force acting on the FEM system. This takes q and qdot as parameters, returns the force in f.
//  stiffness(K, q, qdot) - a function that computes the stiffness (negative second derivative of the potential energy). This takes q and qdot as parameters, returns the stiffness matrix in K.  
//  tmp_qdot - scratch space for storing velocities 
//  tmp_force - scratch space to collect forces
//  tmp_stiffness - scratch space to collect stiffness matrix
//Output:
//  q - set q to the updated generalized coordinate using linearly implicit time integration
//  qdot - set qdot to the updated generalized velocity using linearly implicit time integration
template<typename ENERGY, typename FORCE, typename STIFFNESS> 
inline void implicit_euler(Eigen::VectorXd &q, Eigen::VectorXd &qdot, double dt, 
                            const Eigen::SparseMatrixd &mass,  ENERGY &energy, FORCE &force, STIFFNESS &stiffness, 
                            Eigen::VectorXd &tmp_qdot, Eigen::VectorXd &tmp_force, Eigen::SparseMatrixd &tmp_stiffness) {
    
    //5 iterations to be considered
    Eigen::VectorXd temp_q = q;
    Eigen::VectorXd temp_qdot = qdot;
    tmp_force.setZero();
    tmp_stiffness.setZero();
    //std::cout<<"Energy: "<<energy(q)<<std::endl;
    //double diff = energy(temp_q) - 4.35047e+09;
    //diff = diff * diff;
    for (unsigned idx = 0; idx < 5 ; idx++) //&& diff > 1e+4
    {
	force(tmp_force, temp_q, temp_qdot); //step 1
        stiffness(tmp_stiffness, temp_q, temp_qdot); //step 1
	//g = M * qdot(k - 1) - M* qdot + dt * tmp_force;  
	tmp_force = mass * temp_qdot - mass * qdot - dt * tmp_force;
	//H = M + dt*dt*tmp_stiffness;
	tmp_stiffness = mass - dt * dt * tmp_stiffness;

	double alpha = newtons_method(temp_qdot, energy, force, stiffness, 10, tmp_force, tmp_stiffness); //step 2 with BLS

        //std::cout<<"Alpha is: "<<alpha <<std::endl;
  
	temp_q = q + dt * temp_qdot; //step 3
	
	//diff = energy(temp_q) - 4.35047e+09;
	//diff = diff * diff;
    }

    q = temp_q;
    qdot = temp_qdot;
}
