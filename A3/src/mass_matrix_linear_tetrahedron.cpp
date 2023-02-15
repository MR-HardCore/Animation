 
 #include <mass_matrix_linear_tetrahedron.h>

//Compute the dense mass matrix for a single tetrahedron.

//Input:
//  qdot - generalized velocity for the FEM system
//  element - the 1x4 vertex indices for this tetrahedron
//  density - density of material
//  volume - the undeformed tetrahedron volume
//Output:
//  M - dense 12x12 per-tetrahedron mass matrix

 void mass_matrix_linear_tetrahedron(Eigen::Matrix1212d &M, Eigen::Ref<const Eigen::VectorXd> qdot, Eigen::Ref<const Eigen::RowVectorXi> element, double density, double volume) {

    
    M.setZero();
    M(0, 0) = 2;
    M(0, 3) = 1;
    M(0, 6) = 1;
    M(0, 9) = 1;

    M(1, 1) = 2;
    M(1, 4) = 1;
    M(1, 7) = 1;
    M(1, 10) = 1;

    M(2, 2) = 2;
    M(2, 5) = 1;
    M(2, 8) = 1;
    M(2, 11) = 1;      
  
    M(3, 0) = 1;
    M(3, 3) = 2;
    M(3, 6) = 1;
    M(3, 9) = 1;

    M(4, 1) = 1;
    M(4, 4) = 2;
    M(4, 7) = 1;
    M(4, 10) = 1;

    M(5, 2) = 1;
    M(5, 5) = 2;
    M(5, 8) = 1;
    M(5, 11) = 1;

    M(6, 0) = 1;
    M(6, 3) = 1;
    M(6, 6) = 2;
    M(6, 9) = 1;

    M(7, 1) = 1;
    M(7, 4) = 1;
    M(7, 7) = 2;
    M(7, 10) = 1;

    M(8, 2) = 1;
    M(8, 5) = 1;
    M(8, 8) = 2;
    M(8, 11) = 1;

    M(9, 0) = 1;
    M(9, 3) = 1;
    M(9, 6) = 1;
    M(9, 9) = 2;

    M(10, 1) = 1;
    M(10, 4) = 1;
    M(10, 7) = 1;
    M(10, 10) = 2;

    M(11, 2) = 1;
    M(11, 5) = 1;
    M(11, 8) = 1;
    M(11, 11) = 2;

    double M_fac = density * volume / 20.0;
    M = M * M_fac;


 }
