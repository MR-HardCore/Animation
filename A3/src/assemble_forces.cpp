#include <assemble_forces.h>
#include <iostream>

//Assemble the global force vector for the finite element mesh.

//Input:
//  q - generalized coordinates for the FEM system
//  qdot - generalized velocity for the FEM system
//  V - the nx3 matrix of undeformed vertex positions. Each row is a single undeformed vertex position.
//  T - the mx4 tetrahedron connectivity matrix. Each row contains to indices into V that indicate a spring between those vertices.
//  v0 - the mx1 vector of undeformed tetrahedron volumes
//  C,D - material parameters for the Neo-Hookean model
//Output:
//  f - the vector 3xn vector of forces acting on each node of the mass-spring system
void assemble_forces(Eigen::VectorXd &f, Eigen::Ref<const Eigen::VectorXd> q, Eigen::Ref<const Eigen::MatrixXd> qdot, 
                     Eigen::Ref<const Eigen::MatrixXd> V, Eigen::Ref<const Eigen::MatrixXi> T, Eigen::Ref<const Eigen::VectorXd> v0,
                     double C, double D) { 
        

//dV_linear_tetrahedron_dq

//dV_linear_tetrahedron_dq(Eigen::Vector12d &dV, Eigen::Ref<const Eigen::VectorXd> q,  Eigen::Ref<const Eigen::MatrixXd> V, Eigen::Ref<const Eigen::RowVectorXi> element, double volume,  double C, double D)

    int n_tet = T.rows();
    f.resize(q.size());
    f.setZero();
    for (unsigned idx = 0; idx < n_tet; idx++)
    {   

        Eigen::Vector12d temp_dV;
        temp_dV.setZero();
        int v1 = T(idx, 0), v2 = T(idx, 1), v3 = T(idx, 2), v4 = T(idx, 3);
        dV_linear_tetrahedron_dq(temp_dV, q, V, Eigen::RowVector4i(v1, v2, v3, v4), v0(idx), C, D);
        //node1
        f(v1 * 3    ) -= temp_dV(0);
        f(v1 * 3 + 1) -= temp_dV(1);
        f(v1 * 3 + 2) -= temp_dV(2);
        //node2
        f(v2 * 3    ) -= temp_dV(3);
        f(v2 * 3 + 1) -= temp_dV(4);
        f(v2 * 3 + 2) -= temp_dV(5);
        //node3
        f(v3 * 3    ) -= temp_dV(6);
        f(v3 * 3 + 1) -= temp_dV(7);
        f(v3 * 3 + 2) -= temp_dV(8);
        //node4
        f(v4 * 3    ) -= temp_dV(9);
        f(v4 * 3 + 1) -= temp_dV(10);
        f(v4 * 3 + 2) -= temp_dV(11);
    }
    

};
