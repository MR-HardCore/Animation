#include <V_spring_particle_particle.h>
#include <iostream>

//the potential energy of a spring with 3D end points q0 and qd and undeformed length l0
void V_spring_particle_particle(double &V, Eigen ::Ref<const Eigen::Vector3d> q0,  Eigen::Ref<const Eigen::Vector3d> q1, double l0, double stiffness) {
//Input:
//  q0 - the generalized coordinates of the first node of the spring
//  q1 - the generalized coordinates of the second node of the spring
//  l0 - the undeformed length of the spring
//  stiffness - the stiffness constant for this spring
//Output:
//  V - the potential energy of this spring

    double dq1 = q1(0) - q0(0), 
	   dq2 = q1(1) - q0(1), 
	   dq3 = q1(2) - q0(2), 
           l1 = sqrt(dq1 * dq1 + dq2 * dq2 + dq3 * dq3);

    V = 0;
    V = (l1 - l0) * (l1 - l0) * stiffness;
    //std::cout<<"In VPP: "<<V << " "<<l1<<" "<<l0 <<std::endl;
}
