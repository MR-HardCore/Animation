#include <V_gravity_particle.h>
#include <iostream>
void V_gravity_particle(double &V, Eigen::Ref<const Eigen::Vector3d> q,  double mass, Eigen::Ref<const Eigen::Vector3d> g) {
//Input:
//  q - generalized coordinate of a particle
//  mass - the mass of particles in the mass spring system
//  g - the gravity acceleration vector
//Output:
//  V - the potential energy due to gravity acting on this particle
    
    V = - mass * ( g(0) * q(0) + g(1) * q(1) + g(2) * q(2) );
}
