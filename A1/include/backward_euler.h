#include <Eigen/Dense>

//Input:
//  q - generalized coordiantes for the mass-spring system
//  qdot - generalized velocity for the mass spring system
//  dt - the time step in seconds
//  mass - the mass
//  force(q, qdot) - a function that computes the force acting on the mass as a function. This takes q and qdot as parameters.
//  stiffness(q, qdot) - a function that computes the stiffness (negative second derivative of the potential energy). This takes q and qdot as parameters.  
//Output:
//  q - set q to the updated generalized coordinate using Backward Euler time integration
//  qdot - set qdot to the updated generalized velocity using Backward Euler time integration

template<typename FORCE, typename STIFFNESS> 
inline void backward_euler(Eigen::VectorXd &q, Eigen::VectorXd &qdot, double dt, double mass,  FORCE &force, STIFFNESS &stiffness) {
    //def: q = {x} qdot = {v}
    Eigen::VectorXd f; // '-' sign added
    Eigen::MatrixXd k; // '-' sign added
    force (f, q, qdot);
    stiffness(k, q, qdot);
    double xt = q(0);
    double vt = qdot(0);
    qdot(0) = (vt + f(0) * dt / mass) / (1 - dt * dt * k(0, 0) / mass);
    q(0) = xt + dt * qdot(0);
    //std::cout<<"BE iterated done"<<std::endl;
}
