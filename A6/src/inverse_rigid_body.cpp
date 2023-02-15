#include <inverse_rigid_body.h>

//A method to transform a point from world (deformed) space to body (undeformed) space.

//Input:
//  x - world space position
//  R - rotation from undeformed to world space
//  p - center-of-mass translation
//Output:
//  X - undeformed position of point x
//  qdot - updated generalized velocities

void inverse_rigid_body(Eigen::Vector3d &x, Eigen::Ref<const Eigen::Vector3d> x_world, 
                        Eigen::Ref<const Eigen::Matrix3d> R, Eigen::Ref<const Eigen::Vector3d> p) {
    
    //x = R*X + p
    //X = R^(-1) * (x - p)
    Eigen::Matrix3d inv;
    inv.setZero();
    double R1_1 = R(0, 0), R1_2 = R(0, 1), R1_3 = R(0, 2),
           R2_1 = R(1, 0), R2_2 = R(1, 1), R2_3 = R(1, 2),
           R3_1 = R(2, 0), R3_2 = R(2, 1), R3_3 = R(2, 2);
    
    inv(0, 0) = (R2_2*R3_3-R2_3*R3_2)/(R1_1*R2_2*R3_3-R1_1*R2_3*R3_2-R1_2*R2_1*R3_3+R1_2*R2_3*R3_1+R1_3*R2_1*R3_2-R1_3*R2_2*R3_1);
    inv(0, 1) = -(R1_2*R3_3-R1_3*R3_2)/(R1_1*R2_2*R3_3-R1_1*R2_3*R3_2-R1_2*R2_1*R3_3+R1_2*R2_3*R3_1+R1_3*R2_1*R3_2-R1_3*R2_2*R3_1);
    inv(0, 2) = (R1_2*R2_3-R1_3*R2_2)/(R1_1*R2_2*R3_3-R1_1*R2_3*R3_2-R1_2*R2_1*R3_3+R1_2*R2_3*R3_1+R1_3*R2_1*R3_2-R1_3*R2_2*R3_1);
    inv(1, 0) = -(R2_1*R3_3-R2_3*R3_1)/(R1_1*R2_2*R3_3-R1_1*R2_3*R3_2-R1_2*R2_1*R3_3+R1_2*R2_3*R3_1+R1_3*R2_1*R3_2-R1_3*R2_2*R3_1);
    inv(1, 1) = (R1_1*R3_3-R1_3*R3_1)/(R1_1*R2_2*R3_3-R1_1*R2_3*R3_2-R1_2*R2_1*R3_3+R1_2*R2_3*R3_1+R1_3*R2_1*R3_2-R1_3*R2_2*R3_1);
    inv(1, 2) = -(R1_1*R2_3-R1_3*R2_1)/(R1_1*R2_2*R3_3-R1_1*R2_3*R3_2-R1_2*R2_1*R3_3+R1_2*R2_3*R3_1+R1_3*R2_1*R3_2-R1_3*R2_2*R3_1);
    inv(2, 0) = (R2_1*R3_2-R2_2*R3_1)/(R1_1*R2_2*R3_3-R1_1*R2_3*R3_2-R1_2*R2_1*R3_3+R1_2*R2_3*R3_1+R1_3*R2_1*R3_2-R1_3*R2_2*R3_1);
    inv(2, 1) = -(R1_1*R3_2-R1_2*R3_1)/(R1_1*R2_2*R3_3-R1_1*R2_3*R3_2-R1_2*R2_1*R3_3+R1_2*R2_3*R3_1+R1_3*R2_1*R3_2-R1_3*R2_2*R3_1);
    inv(2, 2) = (R1_1*R2_2-R1_2*R2_1)/(R1_1*R2_2*R3_3-R1_1*R2_3*R3_2-R1_2*R2_1*R3_3+R1_2*R2_3*R3_1+R1_3*R2_1*R3_2-R1_3*R2_2*R3_1);

    x = inv * (x_world - p);
}
