#include <T_particle.h>
#include <iostream>

void T_particle(double &T, Eigen::Ref<const Eigen::VectorXd> qdot, double mass) {
//Input:
//  qdot - generalized velocity for the mass spring system
//  mass - the mass of a particle
//Output:
//  T - kinetic energy of the all particles in the mass-spring system
    //std::cout<<"The size of qdot is: "<<qdot.size()<<std::endl;
    //std::cout<<"The dimension of this p: "<<qdot(0) << " " << qdot(1) << " "<< qdot(2) <<std::endl;
    T = mass * (qdot(0) * qdot(0) + qdot(1) * qdot(1) + qdot(2) * qdot(2)) / 2.0;
    //std::cout <<"Kinetic for this particle is: "<< T <<std::endl;
}
