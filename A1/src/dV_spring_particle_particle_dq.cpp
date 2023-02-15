#include <dV_spring_particle_particle_dq.h>
#include <iostream>


//Input:
//  q - the current generalized coordinates for the mass-spring system
//  stiffness - the stiffness parameter (spring constant) for the mass-spring system

//Output:
// dV - the gradient of the potential energy with respect to the generalised coordinates.

// dv - actually is the force exerted on the object
// q(0) stands for x (displacement)
// q(1) stands for v (velocity)
void dV_spring_particle_particle_dq(Eigen::VectorXd &dV, const Eigen::VectorXd &q, double stiffness) {
    dV.resize(1);
    //std::cout<<"Here is FORCE"<<std::endl;
    dV(0) = q(0) * stiffness;
    //std::cout<<dV(0)<<std::endl;
}
