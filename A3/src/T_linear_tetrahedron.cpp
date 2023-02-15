#include <T_linear_tetrahedron.h>
#include <mass_matrix_linear_tetrahedron.h>
#include <iostream>

//Compute the kinetic energy of a single tetrahedron.

//Input:
// qdot - generalied velocity of FEM system
// element - vertex indices of the element
// density - material density
// volume - volume of tetrahedron
//Output:
//  T - kinetic energy of tetrahedron

void T_linear_tetrahedron(double &T, Eigen::Ref<const Eigen::VectorXd> qdot, Eigen::Ref<const Eigen::RowVectorXi> element, double density, double volume) {

    double qdot11 = qdot(3 * element(0)), qdot12 = qdot(3 * element(0) + 1), qdot13 = qdot(3 * element(0) + 2),
           qdot21 = qdot(3 * element(1)), qdot22 = qdot(3 * element(1) + 1), qdot23 = qdot(3 * element(1) + 2),
           qdot31 = qdot(3 * element(2)), qdot32 = qdot(3 * element(2) + 1), qdot33 = qdot(3 * element(2) + 2),
           qdot41 = qdot(3 * element(3)), qdot42 = qdot(3 * element(3) + 1), qdot43 = qdot(3 * element(3) + 2);
    
    double qdot_square = qdot11 * qdot11 + qdot12 * qdot12 + qdot13 * qdot13 + qdot21 * qdot21 + qdot22 * qdot22 + qdot23 * qdot23 +
                         qdot31 * qdot31 + qdot32 * qdot32 + qdot33 * qdot33 + qdot41 * qdot41 + qdot42 * qdot42 + qdot43 * qdot43;
    

    T = 0;
    T = density * volume * qdot_square;
}
