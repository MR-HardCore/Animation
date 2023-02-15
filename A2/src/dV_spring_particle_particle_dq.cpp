#include <dV_spring_particle_particle_dq.h>

void dV_spring_particle_particle_dq(Eigen::Ref<Eigen::Vector6d> f, Eigen::Ref<const Eigen::Vector3d> q0,  Eigen::Ref<const Eigen::Vector3d>     q1, double l0, double stiffness) {
//Input:
//  q0 - the generalized coordinates of the first node of the spring
//  q1 - the generalized coordinates of the second node of the spring
//  l0 - the undeformed length of the spring
//  stiffness - the stiffness constant for this spring
//Output:
//  f - the 6x1 per spring generalized force vector

    // q0(0) = x1, q0(1) = x2, q0(2) = x3, q1(0) = x4, q1(1) = x5, q1(2) = x6
    double x1 = q0(0), x2 = q0(1), x3 = q0(2), x4 = q1(0), x5 = q1(1), x6 = q1(2), k = stiffness;
    //force is negative gradient
    f(0) = -k*(x1*2.0-x4*2.0)*(l0-sqrt(pow(x1-x4,2.0)+pow(x2-x5,2.0)+pow(x3-x6,2.0)))*1.0/sqrt(pow(x1-x4,2.0)+pow(x2-x5,2.0)+pow(x3-x6,2.0));
    f(1) = -k*(x2*2.0-x5*2.0)*(l0-sqrt(pow(x1-x4,2.0)+pow(x2-x5,2.0)+pow(x3-x6,2.0)))*1.0/sqrt(pow(x1-x4,2.0)+pow(x2-x5,2.0)+pow(x3-x6,2.0));
    f(2) = -k*(x3*2.0-x6*2.0)*(l0-sqrt(pow(x1-x4,2.0)+pow(x2-x5,2.0)+pow(x3-x6,2.0)))*1.0/sqrt(pow(x1-x4,2.0)+pow(x2-x5,2.0)+pow(x3-x6,2.0));
    f(3) = k*(x1*2.0-x4*2.0)*(l0-sqrt(pow(x1-x4,2.0)+pow(x2-x5,2.0)+pow(x3-x6,2.0)))*1.0/sqrt(pow(x1-x4,2.0)+pow(x2-x5,2.0)+pow(x3-x6,2.0));
    f(4) = k*(x2*2.0-x5*2.0)*(l0-sqrt(pow(x1-x4,2.0)+pow(x2-x5,2.0)+pow(x3-x6,2.0)))*1.0/sqrt(pow(x1-x4,2.0)+pow(x2-x5,2.0)+pow(x3-x6,2.0));
    f(5) = k*(x3*2.0-x6*2.0)*(l0-sqrt(pow(x1-x4,2.0)+pow(x2-x5,2.0)+pow(x3-x6,2.0)))*1.0/sqrt(pow(x1-x4,2.0)+pow(x2-x5,2.0)+pow(x3-x6,2.0));

}

