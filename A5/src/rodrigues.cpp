#include <rodrigues.h>
#include <cmath>
#include <iostream>

//The rodrigues formula for computing the marix exponential of a , skew-symmetric matrix.

//  R - rotation matrix
//  omega - angular velocity vector

void rodrigues(Eigen::Matrix3d &R, Eigen::Ref<const Eigen::Vector3d> omega) {

    // the norm of the vector is the angle
    double angle = sqrt(omega(0) * omega(0) + omega(1) * omega(1) + omega(2) * omega(2));
    
    // the normalized vector is the direction
    Eigen::Vector3d n_omega = omega / angle;
    
    /*
     * K Matrix:
     *           |   0 | -kz |  ky |
     *           |  kz |   0 | -kx |
     *           | -ky |  kx |   0 |
     */
    Eigen::Matrix3d K;
    
    K.setZero();
    K(1, 0) =  n_omega(2);
    K(2, 0) = -n_omega(1);
    K(0, 1) = -n_omega(2);
    K(2, 1) =  n_omega(0);
    K(0, 2) =  n_omega(1);
    K(1, 2) = -n_omega(0);

    R.setZero();
    R(0, 0) = 1.0;
    R(1, 1) = 1.0;
    R(2, 2) = 1.0;
    
    //R = I + sin(theta) * K + (1 - cos(theta)) * K^2
    R = R + sin(angle) * K + (1 - cos(angle)) * K * K;

}
