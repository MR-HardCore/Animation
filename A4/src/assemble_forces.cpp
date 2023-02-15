#include <assemble_forces.h>
#include <iostream>

//Assemble the global force vector for the finite element mesh.

//Input:
//  q - generalized coordinates for the FEM system
//  qdot - generalized velocity for the FEM system
//  dX - an mx9 matrix which stores the flattened dphi/dX matrices for each tetrahedron.
//       Convert this values back to 3x3 matrices using the following code (NOTE YOU MUST USE THE TEMPORARY VARIABLE tmp_row):
//       Eigen::Matrix<double, 1,9> tmp_row
//       tmp_row = dX.row(ei); //ei is the triangle index. 
//       Eigen::Map<const Eigen::Matrix3d>(tmp_row.data())
//  V - the nx3 matrix of undeformed vertex positions. Each row is a single undeformed vertex position.
//  F - the mx3 triangle connectivity matrix. Each row contains to indices into V that indicate a spring between those vertices.
//  a0 - the mx1 vector of undeformed triangle areas
//  mu,lambda - material parameters for the cloth material model
//Output:
//  f - the vector 3xn vector of forces acting on each node of the mass-spring system
void assemble_forces(Eigen::VectorXd &f, Eigen::Ref<const Eigen::VectorXd> q, Eigen::Ref<const Eigen::MatrixXd> qdot, Eigen::Ref<const Eigen::MatrixXd> dX,
                     Eigen::Ref<const Eigen::MatrixXd> V, Eigen::Ref<const Eigen::MatrixXi> F, Eigen::Ref<const Eigen::VectorXd> a0,
                     double mu, double lambda) { 
       
    int n_face = F.rows();
    f.resize(q.size());
    f.setZero();
    Eigen::Matrix<double, 1,9> tmp_row;
    tmp_row.setZero();
    for (unsigned idx = 0; idx < n_face; idx++)
    {
        Eigen::Vector9d temp_dV;
        temp_dV.setZero();
        tmp_row = dX.row(idx);
        
        int v1 = F(idx, 0), v2 = F(idx, 1), v3 = F(idx, 2);
     
        dV_membrane_corotational_dq(temp_dV, q, Eigen::Map<const Eigen::Matrix3d>(tmp_row.data()), V, F.row(idx), a0(idx), mu, lambda);
        
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
        
    }
};

