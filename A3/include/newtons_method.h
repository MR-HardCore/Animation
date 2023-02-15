#include <Eigen/Dense>
#include <EigenTypes.h>
#include <iostream>
//Input:
//  x0 - initial point for newtons search
//  f(x) - function that evaluates and returns the cost function at x
//  g(dx, x) - function that evaluates and returns the gradient of the cost function in dx
//  H(dH, x) - function that evaluates and returns the Hessian in dH (as a sparse matrix).
//  max steps - the maximum newton iterations to take
//  tmp_g and tmp_H are scratch space to store gradients and hessians
//Output: 
//  x0 - update x0 to new value
template<typename Objective, typename Jacobian, typename Hessian>
double newtons_method(Eigen::VectorXd &x0, Objective &f, Jacobian &g, Hessian &H, unsigned int maxSteps, Eigen::VectorXd &tmp_g, Eigen::SparseMatrixd &tmp_H) {
   //Implement Newton’s method with backtracking line search. Use the following parameter values: alpha (initial step length) = 1, p (scaling factor) = 0.5, c (ensure sufficient decrease) = 1e-8.

    //ENERGY FUNCTION as Objective f(x0)
    //x0 is qdot0
    //temp_g: gradient at this time
    //temp_H: hessian at this time
    //dV at this time = - H^(-1) * g
    double alpha = 1.0, p = 0.5, c = 1e-8, prev_e = f(x0);
    Eigen::SimplicialLDLT<Eigen::SparseMatrixd> solverLDLT;
    solverLDLT.compute(tmp_H);
    Eigen::VectorXd dV = solverLDLT.solve(Eigen::VectorXd(-tmp_g));
    //Eigen::VectorXd dV = - tmp_H^(-1) * tmp_g;
    //intial updated new_x0 = x0 + alpha*dV;
    Eigen::VectorXd temp_x0 = x0 + dV;
    double offset = c * tmp_g.transpose() * dV + prev_e; 
    //factoring alpha: alpha * p while f(new_x0) > f(x0) + c * temp_g.transpose() * dV
    for (unsigned idx = 0; idx < maxSteps && f(temp_x0) >= offset; idx++)
    {
	alpha = alpha * p;
	temp_x0 = x0 + alpha * dV;
    }

    x0 = temp_x0;
    return alpha;
}

