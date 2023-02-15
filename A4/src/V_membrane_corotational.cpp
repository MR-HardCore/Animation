#include <V_membrane_corotational.h>
#include <iostream>

//Potential energy for the cloth stretching force

//  q - generalized coordinates for the FEM system
//  dX - the 3x3 matrix containing dphi/dX
//  V - the nx3 matrix of undeformed vertex positions. Each row is a single undeformed vertex position.
//  element - the vertex indices of this triangle
//  area - the area of this triangle
//  mu,lambda - material parameters for the cloth material model
//Output:
//  energy- the per-triangle potential energy (the linear model described in the README).

//Allowed to use libigl SVD or Eigen SVD for this part
void V_membrane_corotational(double &energy, Eigen::Ref<const Eigen::VectorXd> q, Eigen::Ref<const Eigen::Matrix3d> dX, 
                          Eigen::Ref<const Eigen::MatrixXd> V, Eigen::Ref<const Eigen::RowVectorXi> element, double area, 
                          double mu, double lambda) {

    int v1 = element(0), v2 = element(1), v3 = element(2);
    
    Eigen::Matrix3d F, U, W, dx; 
    Eigen::Vector3d S;
    
    dx(0, 0) = q(3 * v1 + 0);
    dx(1, 0) = q(3 * v1 + 1);
    dx(2, 0) = q(3 * v1 + 2);
    
    dx(0, 1) = q(3 * v2 + 0);
    dx(1, 1) = q(3 * v2 + 1);
    dx(2, 1) = q(3 * v2 + 2);

    dx(0, 2) = q(3 * v3 + 0);
    dx(1, 2) = q(3 * v3 + 1);
    dx(2, 2) = q(3 * v3 + 2);

    F = dx * dX;

    Eigen::JacobiSVD<Eigen::Matrix3d> svd(F, Eigen::ComputeFullU | Eigen::ComputeFullV);
    U = svd.matrixU();
    W = svd.matrixV();
    S = svd.singularValues();

    double psi = mu * ((S(0) - 1) * (S(0) - 1) + (S(1) - 1) * (S(1) - 1) ) + lambda * (S(0) + S(1) - 2 ) * (S(0) + S(1) - 2 ) * 0.5;

    energy = psi * area;
}
