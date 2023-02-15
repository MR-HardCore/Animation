#include <d2V_spring_particle_particle_dq2.h>
#include <iostream>

//Input:
//  q - the current generalized coordinates for the mass-spring system
//  stiffness - the stiffness parameter (spring constant) for the mass-spring system

//Output:
// H - the second derivtive of the potential energy

// H - is the stiffness of string?
// I believe this can be a single value as well?
void d2V_spring_particle_particle_dq2(Eigen::MatrixXd &H, const Eigen::VectorXd &q, double stiffness) {
    H.resize(1,1);
    //std::cout<<"Here is STIFF"<<std::endl;
    H(0, 0) = stiffness;
}
