#include <dV_membrane_corotational_dq.h>
#include <iostream>

//Gradient of the cloth stretching energy.

//  q - generalized coordinates for the FEM system
//  dX - the 3x3 matrix containing dphi/dX
//  V - the nx3 matrix of undeformed vertex positions. Each row is a single undeformed vertex position.
//  element - the vertex indices of this triangle
//  area - the area of this triangle
//  mu,lambda - material parameters for the cloth material model
//Output:
//  f - the per-triangle gradient of the membrane potential energy (the linear model described in the README).

void dV_membrane_corotational_dq(Eigen::Vector9d &dV, Eigen::Ref<const Eigen::VectorXd> q, Eigen::Ref<const Eigen::Matrix3d> dX, 
                          Eigen::Ref<const Eigen::MatrixXd> V, Eigen::Ref<const Eigen::RowVectorXi> element, double area, 
                          double mu, double lambda) {

    //Deformation Gradient
    Eigen::Matrix3d dx; //deformed tangent matrix 
    Eigen::Matrix3d U;
    Eigen::Vector3d S; 
    Eigen::Matrix3d W; 
    
    //TODO: SVD Here
    int v1 = element(0), v2 = element(1), v3 = element(2);
    Eigen::Matrix3d F;
    F.setZero();

    //def space vertex
    dx(0, 0) = q(3 * v1 + 0);
    dx(1, 0) = q(3 * v1 + 1);
    dx(2, 0) = q(3 * v1 + 2);
    
    dx(0, 1) = q(3 * v2 + 0);
    dx(1, 1) = q(3 * v2 + 1);
    dx(2, 1) = q(3 * v2 + 2);

    dx(0, 2) = q(3 * v3 + 0);
    dx(1, 2) = q(3 * v3 + 1);
    dx(2, 2) = q(3 * v3 + 2);

    //F
    F = dx * dX;

    //SVD
    Eigen::JacobiSVD<Eigen::Matrix3d> svd(F, Eigen::ComputeFullU | Eigen::ComputeFullV);
    U = svd.matrixU();
    W = svd.matrixV();
    S = svd.singularValues();
    
    
    //check for failure
    if(std::fabs(S[2]) > 1e-8) {
        std::cout<<"SVD Fail \n";
        exit(1);
    }

    //Fix for inverted elements (thanks to Danny Kaufman)
    double det = S[0]*S[1];
    
     if(det <= -1e-10)
    {
        if(S[0] < 0) S[0] *= -1;
        if(S[1] < 0) S[1] *= -1;
        if(S[2] < 0) S[2] *= -1;
    }
    
    if(U.determinant() <= 0)
    {
        U(0, 2) *= -1;
        U(1, 2) *= -1;
        U(2, 2) *= -1;
    }
    
    if(W.determinant() <= 0)
    {
        W(0, 2) *= -1;
        W(1, 2) *= -1;
        W(2, 2) *= -1;
    }

    //TODO: energy model gradient
    Eigen::Matrix99d B;
    B.setZero();
    
    B(0, 0) = dX(0, 0);
    B(0, 3) = dX(1, 0);
    B(0, 6) = dX(2, 0);
    B(1, 0) = dX(0, 1);
    B(1, 3) = dX(1, 1);
    B(1, 6) = dX(2, 1);
    B(2, 0) = dX(0, 2);
    B(2, 3) = dX(1, 2);
    B(2, 6) = dX(2, 2);
    B(3, 1) = dX(0, 0);
    B(3, 4) = dX(1, 0);
    B(3, 7) = dX(2, 0);
    B(4, 1) = dX(0, 1);
    B(4, 4) = dX(1, 1);
    B(4, 7) = dX(2, 1);
    B(5, 1) = dX(0, 2);
    B(5, 4) = dX(1, 2);
    B(5, 7) = dX(2, 2);
    B(6, 2) = dX(0, 0);
    B(6, 5) = dX(1, 0);
    B(6, 8) = dX(2, 0);
    B(7, 2) = dX(0, 1);
    B(7, 5) = dX(1, 1);
    B(7, 8) = dX(2, 1);
    B(8, 2) = dX(0, 2);
    B(8, 5) = dX(1, 2);
    B(8, 8) = dX(2, 2);


    Eigen::Matrix3d dpsi;
    dpsi.setZero();
     
    //dpsi(0, 0) is dpsids0 dpsi(1, 1) is dpsids1
    dpsi(0, 0) = (lambda*(S(0)*2.0+S(1)*2.0-4.0))/2.0+mu*(S(0)*2.0-2.0);
    dpsi(1, 1) = (lambda*(S(0)*2.0+S(1)*2.0-4.0))/2.0+mu*(S(1)*2.0-2.0);
    
    Eigen::Matrix3d dpsidF = U * dpsi * W.transpose();

    Eigen::Vector9d grad;
    grad.setZero();

    //flaten dpsidF in row
    grad(0) = dpsidF(0, 0);
    grad(1) = dpsidF(0, 1);
    grad(2) = dpsidF(0, 2);
    grad(3) = dpsidF(1, 0);
    grad(4) = dpsidF(1, 1);
    grad(5) = dpsidF(1, 2);
    grad(6) = dpsidF(2, 0);
    grad(7) = dpsidF(2, 1);
    grad(8) = dpsidF(2, 2);

    dV = area * B.transpose() * grad;
}
