#include <dV_gravity_particle_dq.h>

void dV_gravity_particle_dq(Eigen::Ref<Eigen::Vector3d> f,  double mass, Eigen::Ref<const Eigen::Vector3d> g) {
//Input:
//  mass - mass of this particle
//  g - 3x1 gravity acceleration vector
//Output:
//  f - the 3x1 force of gravity vector acting on the particle

    f(0) = - mass * g(0);
    f(1) = - mass * g(1);
    f(2) = - mass * g(2);
}
