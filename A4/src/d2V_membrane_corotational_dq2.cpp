#include <d2V_membrane_corotational_dq2.h>
#include <iostream>

//Hessian matrix of the cloth stretching energy

//  q - generalized coordinates for the FEM system
//  dX - the 3x3 matrix containing dphi/dX
//  V - the nx3 matrix of undeformed vertex positions. Each row is a single undeformed vertex position.
//  element - the vertex indices of this triangle
//  area - the area of this triangle
//  mu,lambda - material parameters for the cloth material model
//Output:
//  H - the per-triangle Hessian of the potential energy (the linear model described in the README).

void d2V_membrane_corotational_dq2(Eigen::Matrix99d &H, Eigen::Ref<const Eigen::VectorXd> q, Eigen::Ref<const Eigen::Matrix3d> dX, 
                          Eigen::Ref<const Eigen::MatrixXd> V, Eigen::Ref<const Eigen::RowVectorXi> element, double area, 
                          double mu, double lambda) {
    

    //SVD = USW^T
    Eigen::Matrix3d U;
    Eigen::Vector3d S; 
    Eigen::Matrix3d W; 
    Eigen::Matrix3d F; //deformation gradient
    Eigen::Matrix3d dx;
    
    double tol = 1e-5;
    
    //Compute SVD of F here
    int v1 = element(0), v2 = element(1), v3 = element(2);
    
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

    //deal with singularity in the svd gradient
    if(std::fabs(S[0] - S[1]) < tol || std::fabs(S[1] - S[2]) < tol || std::fabs(S[0] - S[2]) < tol) {
        F += Eigen::Matrix3d::Random()*tol;
        Eigen::JacobiSVD<Eigen::Matrix3d> svd2(F, Eigen::ComputeFullU | Eigen::ComputeFullV);
        U = svd2.matrixU();
        W = svd2.matrixV();
        S = svd2.singularValues();
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

    //TODO: compute H, the hessian of the corotational energy


    //DSVD
    Eigen::Tensor3333d dU, dV;
    Eigen::Tensor333d dS;
    dsvd(dU, dS, dV, F);

    Eigen::Matrix3d dpsi, d2psidS2;
    dpsi.setZero();
    d2psidS2.setZero();
    dpsi(0, 0) = (lambda*(S(0)*2.0+S(1)*2.0-4.0))/2.0+mu*(S(0)*2.0-2.0);
    dpsi(1, 1) = (lambda*(S(0)*2.0+S(1)*2.0-4.0))/2.0+mu*(S(1)*2.0-2.0);
    
    d2psidS2(0, 0) = lambda + mu * 2.0;
    d2psidS2(0, 1) = lambda;
    d2psidS2(1, 0) = lambda;
    d2psidS2(1, 1) = lambda + mu * 2.0;

    Eigen::Matrix99d d2psidF;
    d2psidF.setZero();
    
    //Eigen::Tensor3333d d2psidF = dU * dpsi * W.transpose() + U * (d2psidS2 * dS) * W.transpose() + U * dpsi * dV;
    //setup hessian by hand
    // dU: U [3x3] x dF [3x3]
    // dV: V [3x3] x dF [3x3] notice V.T ij means ji in V [3x3]
    // dS: S [3] x dF [3x3]
    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            Eigen::Matrix3d d_dsij, Hij;
            Eigen::Vector3d dSdF, dsij;

            d_dsij.setZero();
            Hij.setZero();
            dSdF.setZero();
            dsij.setZero();

            Eigen::Matrix3d dUdF = dU[i][j];
            Eigen::Matrix3d dWdF = dV[i][j];  

            dSdF = dS[i][j];
            dsij = d2psidS2 * dSdF;

            d_dsij(0, 0) = dsij(0);
            d_dsij(1, 1) = dsij(1);
            d_dsij(2, 2) = dsij(2);
                
            Hij = dUdF * dpsi * W.transpose() + U * d_dsij * W.transpose() + U * dpsi * dWdF.transpose();

            d2psidF(0, i*3 + j) = Hij(0, 0);
            d2psidF(1, i*3 + j) = Hij(0, 1);
            d2psidF(2, i*3 + j) = Hij(0, 2);
            d2psidF(3, i*3 + j) = Hij(1, 0);
            d2psidF(4, i*3 + j) = Hij(1, 1);
            d2psidF(5, i*3 + j) = Hij(1, 2);
            d2psidF(6, i*3 + j) = Hij(2, 0);
            d2psidF(7, i*3 + j) = Hij(2, 1);
            d2psidF(8, i*3 + j) = Hij(2, 2);
        }
    }

    //B is dFdq
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


    H = area * B.transpose() * d2psidF * B;

    //fix errant eigenvalues
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix99d> es(H);
    
    Eigen::MatrixXd DiagEval = es.eigenvalues().real().asDiagonal();
    Eigen::MatrixXd Evec = es.eigenvectors().real();
    
    for (int i = 0; i < 9; ++i) {
        if (es.eigenvalues()[i]<1e-6) {
            DiagEval(i,i) = 1e-3;
        }
    }
    
    H = Evec * DiagEval * Evec.transpose();
   
}
