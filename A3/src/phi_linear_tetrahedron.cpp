#include <phi_linear_tetrahedron.h>
#include <iostream>


//Evaluate the linear shape functions for a tetrahedron. This function returns a 4D vector which contains the values of the shape functions for each vertex at the world space point x (assumed to be inside the element).

//Input:
//  V - the nx3 matrix of undeformed vertex positions. Each row is a single undeformed vertex position.
//  element - the 1x4 vertex indices for this tetrahedron
//  X - the position in the underformed space at which to compute the gradient
//Output:
//  phi - the 4x1 values the basis functions
void phi_linear_tetrahedron(Eigen::Vector4d &phi, Eigen::Ref<const Eigen::MatrixXd> V, Eigen::Ref<const Eigen::RowVectorXi> element, Eigen::Ref<const Eigen::Vector3d> x) {

    //V is undef vertex matrix nx3
    unsigned i1 = element(0), i2 = element(1), i3 = element(2), i4 = element(3);
    double X_11 = V(i1, 0), X_12 = V(i1, 1), X_13 = V(i1, 2), // X1 - undef vertex1
           X_21 = V(i2, 0), X_22 = V(i2, 1), X_23 = V(i2, 2), // X2 - undef vertex2
           X_31 = V(i3, 0), X_32 = V(i3, 1), X_33 = V(i3, 2), // X3 - undef vertex3
           X_41 = V(i4, 0), X_42 = V(i4, 1), X_43 = V(i4, 2), // X4 - undef vertex4
           X1 = x(0), X2 = x(1), X3 = x(2); // X - undef vertex

    phi.setZero();

    phi(0) = (X3*(X_21*X_32 - X_22*X_31 - X_21*X_42 + X_22*X_41 + X_31*X_42 - X_32*X_41))/(X_11*X_22*X_33 - X_11*X_23*X_32 - X_12*X_21*X_33 + X_12*X_23*X_31 + X_13*X_21*X_32 - X_13*X_22*X_31 - X_11*X_22*X_43 + X_11*X_23*X_42 + X_12*X_21*X_43 - X_12*X_23*X_41 - X_13*X_21*X_42 + X_13*X_22*X_41 + X_11*X_32*X_43 - X_11*X_33*X_42 - X_12*X_31*X_43 + X_12*X_33*X_41 + X_13*X_31*X_42 - X_13*X_32*X_41 - X_21*X_32*X_43 + X_21*X_33*X_42 + X_22*X_31*X_43 - X_22*X_33*X_41 - X_23*X_31*X_42 + X_23*X_32*X_41) - (X_21*X_32*X_43 - X_21*X_33*X_42 - X_22*X_31*X_43 + X_22*X_33*X_41 + X_23*X_31*X_42 - X_23*X_32*X_41)/(X_11*X_22*X_33 - X_11*X_23*X_32 - X_12*X_21*X_33 + X_12*X_23*X_31 + X_13*X_21*X_32 - X_13*X_22*X_31 - X_11*X_22*X_43 + X_11*X_23*X_42 + X_12*X_21*X_43 - X_12*X_23*X_41 - X_13*X_21*X_42 + X_13*X_22*X_41 + X_11*X_32*X_43 - X_11*X_33*X_42 - X_12*X_31*X_43 + X_12*X_33*X_41 + X_13*X_31*X_42 - X_13*X_32*X_41 - X_21*X_32*X_43 + X_21*X_33*X_42 + X_22*X_31*X_43 - X_22*X_33*X_41 - X_23*X_31*X_42 + X_23*X_32*X_41) - (X2*(X_21*X_33 - X_23*X_31 - X_21*X_43 + X_23*X_41 + X_31*X_43 - X_33*X_41))/(X_11*X_22*X_33 - X_11*X_23*X_32 - X_12*X_21*X_33 + X_12*X_23*X_31 + X_13*X_21*X_32 - X_13*X_22*X_31 - X_11*X_22*X_43 + X_11*X_23*X_42 + X_12*X_21*X_43 - X_12*X_23*X_41 - X_13*X_21*X_42 + X_13*X_22*X_41 + X_11*X_32*X_43 - X_11*X_33*X_42 - X_12*X_31*X_43 + X_12*X_33*X_41 + X_13*X_31*X_42 - X_13*X_32*X_41 - X_21*X_32*X_43 + X_21*X_33*X_42 + X_22*X_31*X_43 - X_22*X_33*X_41 - X_23*X_31*X_42 + X_23*X_32*X_41) + (X1*(X_22*X_33 - X_23*X_32 - X_22*X_43 + X_23*X_42 + X_32*X_43 - X_33*X_42))/(X_11*X_22*X_33 - X_11*X_23*X_32 - X_12*X_21*X_33 + X_12*X_23*X_31 + X_13*X_21*X_32 - X_13*X_22*X_31 - X_11*X_22*X_43 + X_11*X_23*X_42 + X_12*X_21*X_43 - X_12*X_23*X_41 - X_13*X_21*X_42 + X_13*X_22*X_41 + X_11*X_32*X_43 - X_11*X_33*X_42 - X_12*X_31*X_43 + X_12*X_33*X_41 + X_13*X_31*X_42 - X_13*X_32*X_41 - X_21*X_32*X_43 + X_21*X_33*X_42 + X_22*X_31*X_43 - X_22*X_33*X_41 - X_23*X_31*X_42 + X_23*X_32*X_41);
    
    phi(1) = (X_11*X_32*X_43 - X_11*X_33*X_42 - X_12*X_31*X_43 + X_12*X_33*X_41 + X_13*X_31*X_42 - X_13*X_32*X_41)/(X_11*X_22*X_33 - X_11*X_23*X_32 - X_12*X_21*X_33 + X_12*X_23*X_31 + X_13*X_21*X_32 - X_13*X_22*X_31 - X_11*X_22*X_43 + X_11*X_23*X_42 + X_12*X_21*X_43 - X_12*X_23*X_41 - X_13*X_21*X_42 + X_13*X_22*X_41 + X_11*X_32*X_43 - X_11*X_33*X_42 - X_12*X_31*X_43 + X_12*X_33*X_41 + X_13*X_31*X_42 - X_13*X_32*X_41 - X_21*X_32*X_43 + X_21*X_33*X_42 + X_22*X_31*X_43 - X_22*X_33*X_41 - X_23*X_31*X_42 + X_23*X_32*X_41) - (X3*(X_11*X_32 - X_12*X_31 - X_11*X_42 + X_12*X_41 + X_31*X_42 - X_32*X_41))/(X_11*X_22*X_33 - X_11*X_23*X_32 - X_12*X_21*X_33 + X_12*X_23*X_31 + X_13*X_21*X_32 - X_13*X_22*X_31 - X_11*X_22*X_43 + X_11*X_23*X_42 + X_12*X_21*X_43 - X_12*X_23*X_41 - X_13*X_21*X_42 + X_13*X_22*X_41 + X_11*X_32*X_43 - X_11*X_33*X_42 - X_12*X_31*X_43 + X_12*X_33*X_41 + X_13*X_31*X_42 - X_13*X_32*X_41 - X_21*X_32*X_43 + X_21*X_33*X_42 + X_22*X_31*X_43 - X_22*X_33*X_41 - X_23*X_31*X_42 + X_23*X_32*X_41) + (X2*(X_11*X_33 - X_13*X_31 - X_11*X_43 + X_13*X_41 + X_31*X_43 - X_33*X_41))/(X_11*X_22*X_33 - X_11*X_23*X_32 - X_12*X_21*X_33 + X_12*X_23*X_31 + X_13*X_21*X_32 - X_13*X_22*X_31 - X_11*X_22*X_43 + X_11*X_23*X_42 + X_12*X_21*X_43 - X_12*X_23*X_41 - X_13*X_21*X_42 + X_13*X_22*X_41 + X_11*X_32*X_43 - X_11*X_33*X_42 - X_12*X_31*X_43 + X_12*X_33*X_41 + X_13*X_31*X_42 - X_13*X_32*X_41 - X_21*X_32*X_43 + X_21*X_33*X_42 + X_22*X_31*X_43 - X_22*X_33*X_41 - X_23*X_31*X_42 + X_23*X_32*X_41) - (X1*(X_12*X_33 - X_13*X_32 - X_12*X_43 + X_13*X_42 + X_32*X_43 - X_33*X_42))/(X_11*X_22*X_33 - X_11*X_23*X_32 - X_12*X_21*X_33 + X_12*X_23*X_31 + X_13*X_21*X_32 - X_13*X_22*X_31 - X_11*X_22*X_43 + X_11*X_23*X_42 + X_12*X_21*X_43 - X_12*X_23*X_41 - X_13*X_21*X_42 + X_13*X_22*X_41 + X_11*X_32*X_43 - X_11*X_33*X_42 - X_12*X_31*X_43 + X_12*X_33*X_41 + X_13*X_31*X_42 - X_13*X_32*X_41 - X_21*X_32*X_43 + X_21*X_33*X_42 + X_22*X_31*X_43 - X_22*X_33*X_41 - X_23*X_31*X_42 + X_23*X_32*X_41);
    
    phi(2) = (X3*(X_11*X_22 - X_12*X_21 - X_11*X_42 + X_12*X_41 + X_21*X_42 - X_22*X_41))/(X_11*X_22*X_33 - X_11*X_23*X_32 - X_12*X_21*X_33 + X_12*X_23*X_31 + X_13*X_21*X_32 - X_13*X_22*X_31 - X_11*X_22*X_43 + X_11*X_23*X_42 + X_12*X_21*X_43 - X_12*X_23*X_41 - X_13*X_21*X_42 + X_13*X_22*X_41 + X_11*X_32*X_43 - X_11*X_33*X_42 - X_12*X_31*X_43 + X_12*X_33*X_41 + X_13*X_31*X_42 - X_13*X_32*X_41 - X_21*X_32*X_43 + X_21*X_33*X_42 + X_22*X_31*X_43 - X_22*X_33*X_41 - X_23*X_31*X_42 + X_23*X_32*X_41) - (X_11*X_22*X_43 - X_11*X_23*X_42 - X_12*X_21*X_43 + X_12*X_23*X_41 + X_13*X_21*X_42 - X_13*X_22*X_41)/(X_11*X_22*X_33 - X_11*X_23*X_32 - X_12*X_21*X_33 + X_12*X_23*X_31 + X_13*X_21*X_32 - X_13*X_22*X_31 - X_11*X_22*X_43 + X_11*X_23*X_42 + X_12*X_21*X_43 - X_12*X_23*X_41 - X_13*X_21*X_42 + X_13*X_22*X_41 + X_11*X_32*X_43 - X_11*X_33*X_42 - X_12*X_31*X_43 + X_12*X_33*X_41 + X_13*X_31*X_42 - X_13*X_32*X_41 - X_21*X_32*X_43 + X_21*X_33*X_42 + X_22*X_31*X_43 - X_22*X_33*X_41 - X_23*X_31*X_42 + X_23*X_32*X_41) - (X2*(X_11*X_23 - X_13*X_21 - X_11*X_43 + X_13*X_41 + X_21*X_43 - X_23*X_41))/(X_11*X_22*X_33 - X_11*X_23*X_32 - X_12*X_21*X_33 + X_12*X_23*X_31 + X_13*X_21*X_32 - X_13*X_22*X_31 - X_11*X_22*X_43 + X_11*X_23*X_42 + X_12*X_21*X_43 - X_12*X_23*X_41 - X_13*X_21*X_42 + X_13*X_22*X_41 + X_11*X_32*X_43 - X_11*X_33*X_42 - X_12*X_31*X_43 + X_12*X_33*X_41 + X_13*X_31*X_42 - X_13*X_32*X_41 - X_21*X_32*X_43 + X_21*X_33*X_42 + X_22*X_31*X_43 - X_22*X_33*X_41 - X_23*X_31*X_42 + X_23*X_32*X_41) + (X1*(X_12*X_23 - X_13*X_22 - X_12*X_43 + X_13*X_42 + X_22*X_43 - X_23*X_42))/(X_11*X_22*X_33 - X_11*X_23*X_32 - X_12*X_21*X_33 + X_12*X_23*X_31 + X_13*X_21*X_32 - X_13*X_22*X_31 - X_11*X_22*X_43 + X_11*X_23*X_42 + X_12*X_21*X_43 - X_12*X_23*X_41 - X_13*X_21*X_42 + X_13*X_22*X_41 + X_11*X_32*X_43 - X_11*X_33*X_42 - X_12*X_31*X_43 + X_12*X_33*X_41 + X_13*X_31*X_42 - X_13*X_32*X_41 - X_21*X_32*X_43 + X_21*X_33*X_42 + X_22*X_31*X_43 - X_22*X_33*X_41 - X_23*X_31*X_42 + X_23*X_32*X_41);

    phi(3) = (X_11*X_22*X_33 - X_11*X_23*X_32 - X_12*X_21*X_33 + X_12*X_23*X_31 + X_13*X_21*X_32 - X_13*X_22*X_31)/(X_11*X_22*X_33 - X_11*X_23*X_32 - X_12*X_21*X_33 + X_12*X_23*X_31 + X_13*X_21*X_32 - X_13*X_22*X_31 - X_11*X_22*X_43 + X_11*X_23*X_42 + X_12*X_21*X_43 - X_12*X_23*X_41 - X_13*X_21*X_42 + X_13*X_22*X_41 + X_11*X_32*X_43 - X_11*X_33*X_42 - X_12*X_31*X_43 + X_12*X_33*X_41 + X_13*X_31*X_42 - X_13*X_32*X_41 - X_21*X_32*X_43 + X_21*X_33*X_42 + X_22*X_31*X_43 - X_22*X_33*X_41 - X_23*X_31*X_42 + X_23*X_32*X_41) - (X3*(X_11*X_22 - X_12*X_21 - X_11*X_32 + X_12*X_31 + X_21*X_32 - X_22*X_31))/(X_11*X_22*X_33 - X_11*X_23*X_32 - X_12*X_21*X_33 + X_12*X_23*X_31 + X_13*X_21*X_32 - X_13*X_22*X_31 - X_11*X_22*X_43 + X_11*X_23*X_42 + X_12*X_21*X_43 - X_12*X_23*X_41 - X_13*X_21*X_42 + X_13*X_22*X_41 + X_11*X_32*X_43 - X_11*X_33*X_42 - X_12*X_31*X_43 + X_12*X_33*X_41 + X_13*X_31*X_42 - X_13*X_32*X_41 - X_21*X_32*X_43 + X_21*X_33*X_42 + X_22*X_31*X_43 - X_22*X_33*X_41 - X_23*X_31*X_42 + X_23*X_32*X_41) + (X2*(X_11*X_23 - X_13*X_21 - X_11*X_33 + X_13*X_31 + X_21*X_33 - X_23*X_31))/(X_11*X_22*X_33 - X_11*X_23*X_32 - X_12*X_21*X_33 + X_12*X_23*X_31 + X_13*X_21*X_32 - X_13*X_22*X_31 - X_11*X_22*X_43 + X_11*X_23*X_42 + X_12*X_21*X_43 - X_12*X_23*X_41 - X_13*X_21*X_42 + X_13*X_22*X_41 + X_11*X_32*X_43 - X_11*X_33*X_42 - X_12*X_31*X_43 + X_12*X_33*X_41 + X_13*X_31*X_42 - X_13*X_32*X_41 - X_21*X_32*X_43 + X_21*X_33*X_42 + X_22*X_31*X_43 - X_22*X_33*X_41 - X_23*X_31*X_42 + X_23*X_32*X_41) - (X1*(X_12*X_23 - X_13*X_22 - X_12*X_33 + X_13*X_32 + X_22*X_33 - X_23*X_32))/(X_11*X_22*X_33 - X_11*X_23*X_32 - X_12*X_21*X_33 + X_12*X_23*X_31 + X_13*X_21*X_32 - X_13*X_22*X_31 - X_11*X_22*X_43 + X_11*X_23*X_42 + X_12*X_21*X_43 - X_12*X_23*X_41 - X_13*X_21*X_42 + X_13*X_22*X_41 + X_11*X_32*X_43 - X_11*X_33*X_42 - X_12*X_31*X_43 + X_12*X_33*X_41 + X_13*X_31*X_42 - X_13*X_32*X_41 - X_21*X_32*X_43 + X_21*X_33*X_42 + X_22*X_31*X_43 - X_22*X_33*X_41 - X_23*X_31*X_42 + X_23*X_32*X_41);
    
}