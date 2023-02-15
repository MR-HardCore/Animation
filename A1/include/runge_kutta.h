//Input:
//  q - generalized coordiantes for the mass-spring system
//  qdot - generalized velocity for the mass spring system
//  dt - the time step in seconds
//  mass - the mass
//  force(q, qdot) - a function that computes the force acting on the mass as a function. This takes q and qdot as parameters.
//Output:
//  q - set q to the updated generalized coordinate using Runge-Kutta time integration
//  qdot - set qdot to the updated generalized velocity using Runge-Kutta time integration

template<typename FORCE> 
inline void runge_kutta(Eigen::VectorXd &q, Eigen::VectorXd &qdot, double dt, double mass,  FORCE &force) {
    //def: q = {x} qdot = {v}
    Eigen::VectorXd f;
    force (f, q, qdot); // '-' sign added
    double cv = f(0) * dt / ( q(0) * mass); // '-' sign added
    double cx = dt;
    double xt = q(0);
    double vt = qdot(0);
    
    //4th rk 
    double kv1 = cv * xt;
    double kx1 = cx * vt;
    double kv2 = cv * ( xt + kx1 / 2 );
    double kx2 = cx * ( vt + kv1 / 2 );
    double kv3 = cv * ( xt + kx2 / 2 );
    double kx3 = cx * ( vt + kv2 / 2 );
    double kv4 = cv * ( xt + kx3 );
    double kx4 = cx * ( vt + kv3 );
    qdot(0) = vt + ( kv1 + 2 * kv2 + 2 * kv3 + kv4 ) / 6;
    q(0) = xt + ( kx1 +2 * kx2 + 2 * kx3 + kx4 ) / 6;

    //std::cout<<"RK iteration done"<<std::endl;
}
