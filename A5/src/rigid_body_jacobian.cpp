#include <rigid_body_jacobian.h>
#include <iostream>

//The Jacobian of the rigid-to-world transform.

//Input:
//  R - rotation matrix for rigid body
//  p - world space position of center-of-mass
//  X -  undeformed position at which to compute the Jacobian.
//Output:
//  J - the rigid body jacobian acting on the undeformed space point X.

void rigid_body_jacobian(Eigen::Matrix36d &J, 
                         Eigen::Ref<const Eigen::Matrix3d> R, Eigen::Ref<const Eigen::Vector3d> p, 
                         Eigen::Ref<const Eigen::Vector3d> x) {
    //N(X) is rigid_body_jacobian!!!
    // N(X) = [R * [X].T * R.T    I ]
    // x = RX + p
    // X = R^(-1) * (x - p)
    
    /*
    //inverse matrix of R: R_inv
    double R1_1 = R(0, 0), R1_2 = R(0, 1), R1_3 = R(0, 2),
           R2_1 = R(1, 0), R2_2 = R(1, 1), R2_3 = R(1, 2),
           R3_1 = R(2, 0), R3_2 = R(2, 1), R3_3 = R(2, 2);
    
    Eigen::Matrix3d inv;
    inv.setZero();
    
    inv(0, 0) = (R2_2*R3_3-R2_3*R3_2)/(R1_1*R2_2*R3_3-R1_1*R2_3*R3_2-R1_2*R2_1*R3_3+R1_2*R2_3*R3_1+R1_3*R2_1*R3_2-R1_3*R2_2*R3_1);
    inv(0, 1) = -(R1_2*R3_3-R1_3*R3_2)/(R1_1*R2_2*R3_3-R1_1*R2_3*R3_2-R1_2*R2_1*R3_3+R1_2*R2_3*R3_1+R1_3*R2_1*R3_2-R1_3*R2_2*R3_1);
    inv(0, 2)= (R1_2*R2_3-R1_3*R2_2)/(R1_1*R2_2*R3_3-R1_1*R2_3*R3_2-R1_2*R2_1*R3_3+R1_2*R2_3*R3_1+R1_3*R2_1*R3_2-R1_3*R2_2*R3_1);
    inv(1, 0) = -(R2_1*R3_3-R2_3*R3_1)/(R1_1*R2_2*R3_3-R1_1*R2_3*R3_2-R1_2*R2_1*R3_3+R1_2*R2_3*R3_1+R1_3*R2_1*R3_2-R1_3*R2_2*R3_1);
    inv(1, 1) = (R1_1*R3_3-R1_3*R3_1)/(R1_1*R2_2*R3_3-R1_1*R2_3*R3_2-R1_2*R2_1*R3_3+R1_2*R2_3*R3_1+R1_3*R2_1*R3_2-R1_3*R2_2*R3_1);
    inv(1, 2) = -(R1_1*R2_3-R1_3*R2_1)/(R1_1*R2_2*R3_3-R1_1*R2_3*R3_2-R1_2*R2_1*R3_3+R1_2*R2_3*R3_1+R1_3*R2_1*R3_2-R1_3*R2_2*R3_1);
    inv(2, 0) = (R2_1*R3_2-R2_2*R3_1)/(R1_1*R2_2*R3_3-R1_1*R2_3*R3_2-R1_2*R2_1*R3_3+R1_2*R2_3*R3_1+R1_3*R2_1*R3_2-R1_3*R2_2*R3_1);
    inv(2, 1) = -(R1_1*R3_2-R1_2*R3_1)/(R1_1*R2_2*R3_3-R1_1*R2_3*R3_2-R1_2*R2_1*R3_3+R1_2*R2_3*R3_1+R1_3*R2_1*R3_2-R1_3*R2_2*R3_1);
    inv(2, 2) = (R1_1*R2_2-R1_2*R2_1)/(R1_1*R2_2*R3_3-R1_1*R2_3*R3_2-R1_2*R2_1*R3_3+R1_2*R2_3*R3_1+R1_3*R2_1*R3_2-R1_3*R2_2*R3_1);
    
    Eigen::Vector3d X = inv * (x - p);
    
    Eigen::Matrix3d c_X;
    c_X.setZero();
    c_X(1, 0) =  X(2);
    c_X(2, 0) = -X(1);
    c_X(0, 1) = -X(2);
    c_X(2, 1) =  X(0);
    c_X(0, 2) =  X(1);
    c_X(1, 2) = -X(0);
    */
    
    Eigen::Matrix3d c_X;
    c_X.setZero();
    c_X(1, 0) =  x(2);
    c_X(2, 0) = -x(1);
    c_X(0, 1) = -x(2);
    c_X(2, 1) =  x(0);
    c_X(0, 2) =  x(1);
    c_X(1, 2) = -x(0);
    J.setZero();
    //N(X) = [R * [X].T * R.T    I ]
    //might be wrong with block
    J.block(0, 0, 3, 3) = R * c_X.transpose() * R.transpose();
    
    J(0, 3) = 1.0;
    J(1, 4) = 1.0;
    J(2, 5) = 1.0;
    //exit(0);
}

