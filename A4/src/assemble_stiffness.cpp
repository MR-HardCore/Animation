#include <assemble_stiffness.h>
#include <iostream>

//Assemble the global stiffness matrix for the finite element mesh.

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
//  K- the 3n by 3n sparse stiffness matrix. 
void assemble_stiffness(Eigen::SparseMatrixd &K, Eigen::Ref<const Eigen::VectorXd> q, Eigen::Ref<const Eigen::VectorXd> qdot, Eigen::Ref<const Eigen::MatrixXd> dX,
                     Eigen::Ref<const Eigen::MatrixXd> V, Eigen::Ref<const Eigen::MatrixXi> F, Eigen::Ref<const Eigen::VectorXd> a0, 
                     double mu, double lambda) { 
      
    std::vector<Eigen::Triplet<double>> trip_list;
    int n_face = F.rows();
    Eigen::Matrix<double, 1, 9> tmp_row;
    tmp_row.setZero();
    for (unsigned idx = 0; idx < n_face; idx++)
    {
        Eigen::Matrix99d temp_H;
        temp_H.setZero();
        int v1 = F(idx, 0), v2 = F(idx, 1), v3 = F(idx, 2);
        tmp_row = dX.row(idx);
        d2V_membrane_corotational_dq2(temp_H, q, Eigen::Map<const Eigen::Matrix3d>(tmp_row.data()), V, F.row(idx), a0(idx), mu, lambda);

        //TODO:set H to corresponding K
        //v1v1
        trip_list.push_back(Eigen::Triplet<double>(v1 * 3    , v1 * 3    , - temp_H(0 , 0)));
        trip_list.push_back(Eigen::Triplet<double>(v1 * 3    , v1 * 3 + 1, - temp_H(0 , 1)));
        trip_list.push_back(Eigen::Triplet<double>(v1 * 3    , v1 * 3 + 2, - temp_H(0 , 2)));
        trip_list.push_back(Eigen::Triplet<double>(v1 * 3 + 1, v1 * 3    , - temp_H(1 , 0)));
        trip_list.push_back(Eigen::Triplet<double>(v1 * 3 + 1, v1 * 3 + 1, - temp_H(1 , 1)));
        trip_list.push_back(Eigen::Triplet<double>(v1 * 3 + 1, v1 * 3 + 2, - temp_H(1 , 2)));
        trip_list.push_back(Eigen::Triplet<double>(v1 * 3 + 2, v1 * 3    , - temp_H(2 , 0)));
        trip_list.push_back(Eigen::Triplet<double>(v1 * 3 + 2, v1 * 3 + 1, - temp_H(2 , 1)));
        trip_list.push_back(Eigen::Triplet<double>(v1 * 3 + 2, v1 * 3 + 2, - temp_H(2 , 2)));

        //v1v2
        trip_list.push_back(Eigen::Triplet<double>(v1 * 3    , v2 * 3    , - temp_H(0 , 3)));
        trip_list.push_back(Eigen::Triplet<double>(v1 * 3    , v2 * 3 + 1, - temp_H(0 , 4)));
        trip_list.push_back(Eigen::Triplet<double>(v1 * 3    , v2 * 3 + 2, - temp_H(0 , 5)));
        trip_list.push_back(Eigen::Triplet<double>(v1 * 3 + 1, v2 * 3    , - temp_H(1 , 3)));
        trip_list.push_back(Eigen::Triplet<double>(v1 * 3 + 1, v2 * 3 + 1, - temp_H(1 , 4)));
        trip_list.push_back(Eigen::Triplet<double>(v1 * 3 + 1, v2 * 3 + 2, - temp_H(1 , 5)));
        trip_list.push_back(Eigen::Triplet<double>(v1 * 3 + 2, v2 * 3    , - temp_H(2 , 3)));
        trip_list.push_back(Eigen::Triplet<double>(v1 * 3 + 2, v2 * 3 + 1, - temp_H(2 , 4)));
        trip_list.push_back(Eigen::Triplet<double>(v1 * 3 + 2, v2 * 3 + 2, - temp_H(2 , 5)));

        //v1v3
        trip_list.push_back(Eigen::Triplet<double>(v1 * 3    , v3 * 3    , - temp_H(0 , 6)));
        trip_list.push_back(Eigen::Triplet<double>(v1 * 3    , v3 * 3 + 1, - temp_H(0 , 7)));
        trip_list.push_back(Eigen::Triplet<double>(v1 * 3    , v3 * 3 + 2, - temp_H(0 , 8)));
        trip_list.push_back(Eigen::Triplet<double>(v1 * 3 + 1, v3 * 3    , - temp_H(1 , 6)));
        trip_list.push_back(Eigen::Triplet<double>(v1 * 3 + 1, v3 * 3 + 1, - temp_H(1 , 7)));
        trip_list.push_back(Eigen::Triplet<double>(v1 * 3 + 1, v3 * 3 + 2, - temp_H(1 , 8)));
        trip_list.push_back(Eigen::Triplet<double>(v1 * 3 + 2, v3 * 3    , - temp_H(2 , 6)));
        trip_list.push_back(Eigen::Triplet<double>(v1 * 3 + 2, v3 * 3 + 1, - temp_H(2 , 7)));
        trip_list.push_back(Eigen::Triplet<double>(v1 * 3 + 2, v3 * 3 + 2, - temp_H(2 , 8)));

        //v2v1
        trip_list.push_back(Eigen::Triplet<double>(v2 * 3    , v1 * 3    , - temp_H(3 , 0)));
        trip_list.push_back(Eigen::Triplet<double>(v2 * 3    , v1 * 3 + 1, - temp_H(3 , 1)));
        trip_list.push_back(Eigen::Triplet<double>(v2 * 3    , v1 * 3 + 2, - temp_H(3 , 2)));
        trip_list.push_back(Eigen::Triplet<double>(v2 * 3 + 1, v1 * 3    , - temp_H(4 , 0)));
        trip_list.push_back(Eigen::Triplet<double>(v2 * 3 + 1, v1 * 3 + 1, - temp_H(4 , 1)));
        trip_list.push_back(Eigen::Triplet<double>(v2 * 3 + 1, v1 * 3 + 2, - temp_H(4 , 2)));
        trip_list.push_back(Eigen::Triplet<double>(v2 * 3 + 2, v1 * 3    , - temp_H(5 , 0)));
        trip_list.push_back(Eigen::Triplet<double>(v2 * 3 + 2, v1 * 3 + 1, - temp_H(5 , 1)));
        trip_list.push_back(Eigen::Triplet<double>(v2 * 3 + 2, v1 * 3 + 2, - temp_H(5 , 2)));

        //v2v2
        trip_list.push_back(Eigen::Triplet<double>(v2 * 3    , v2 * 3    , - temp_H(3 , 3)));
        trip_list.push_back(Eigen::Triplet<double>(v2 * 3    , v2 * 3 + 1, - temp_H(3 , 4)));
        trip_list.push_back(Eigen::Triplet<double>(v2 * 3    , v2 * 3 + 2, - temp_H(3 , 5)));
        trip_list.push_back(Eigen::Triplet<double>(v2 * 3 + 1, v2 * 3    , - temp_H(4 , 3)));
        trip_list.push_back(Eigen::Triplet<double>(v2 * 3 + 1, v2 * 3 + 1, - temp_H(4 , 4)));
        trip_list.push_back(Eigen::Triplet<double>(v2 * 3 + 1, v2 * 3 + 2, - temp_H(4 , 5)));
        trip_list.push_back(Eigen::Triplet<double>(v2 * 3 + 2, v2 * 3    , - temp_H(5 , 3)));
        trip_list.push_back(Eigen::Triplet<double>(v2 * 3 + 2, v2 * 3 + 1, - temp_H(5 , 4)));
        trip_list.push_back(Eigen::Triplet<double>(v2 * 3 + 2, v2 * 3 + 2, - temp_H(5 , 5)));

        //v2v3
        trip_list.push_back(Eigen::Triplet<double>(v2 * 3    , v3 * 3    , - temp_H(3 , 6)));
        trip_list.push_back(Eigen::Triplet<double>(v2 * 3    , v3 * 3 + 1, - temp_H(3 , 7)));
        trip_list.push_back(Eigen::Triplet<double>(v2 * 3    , v3 * 3 + 2, - temp_H(3 , 8)));
        trip_list.push_back(Eigen::Triplet<double>(v2 * 3 + 1, v3 * 3    , - temp_H(4 , 6)));
        trip_list.push_back(Eigen::Triplet<double>(v2 * 3 + 1, v3 * 3 + 1, - temp_H(4 , 7)));
        trip_list.push_back(Eigen::Triplet<double>(v2 * 3 + 1, v3 * 3 + 2, - temp_H(4 , 8)));
        trip_list.push_back(Eigen::Triplet<double>(v2 * 3 + 2, v3 * 3    , - temp_H(5 , 6)));
        trip_list.push_back(Eigen::Triplet<double>(v2 * 3 + 2, v3 * 3 + 1, - temp_H(5 , 7)));
        trip_list.push_back(Eigen::Triplet<double>(v2 * 3 + 2, v3 * 3 + 2, - temp_H(5 , 8)));


        //v3v1
        trip_list.push_back(Eigen::Triplet<double>(v3 * 3    , v1 * 3    , - temp_H(6 , 0)));
        trip_list.push_back(Eigen::Triplet<double>(v3 * 3    , v1 * 3 + 1, - temp_H(6 , 1)));
        trip_list.push_back(Eigen::Triplet<double>(v3 * 3    , v1 * 3 + 2, - temp_H(6 , 2)));
        trip_list.push_back(Eigen::Triplet<double>(v3 * 3 + 1, v1 * 3    , - temp_H(7 , 0)));
        trip_list.push_back(Eigen::Triplet<double>(v3 * 3 + 1, v1 * 3 + 1, - temp_H(7 , 1)));
        trip_list.push_back(Eigen::Triplet<double>(v3 * 3 + 1, v1 * 3 + 2, - temp_H(7 , 2)));
        trip_list.push_back(Eigen::Triplet<double>(v3 * 3 + 2, v1 * 3    , - temp_H(8 , 0)));
        trip_list.push_back(Eigen::Triplet<double>(v3 * 3 + 2, v1 * 3 + 1, - temp_H(8 , 1)));
        trip_list.push_back(Eigen::Triplet<double>(v3 * 3 + 2, v1 * 3 + 2, - temp_H(8 , 2)));

        //v3v2
        trip_list.push_back(Eigen::Triplet<double>(v3 * 3    , v2 * 3    , - temp_H(6 , 3)));
        trip_list.push_back(Eigen::Triplet<double>(v3 * 3    , v2 * 3 + 1, - temp_H(6 , 4)));
        trip_list.push_back(Eigen::Triplet<double>(v3 * 3    , v2 * 3 + 2, - temp_H(6 , 5)));
        trip_list.push_back(Eigen::Triplet<double>(v3 * 3 + 1, v2 * 3    , - temp_H(7 , 3)));
        trip_list.push_back(Eigen::Triplet<double>(v3 * 3 + 1, v2 * 3 + 1, - temp_H(7 , 4)));
        trip_list.push_back(Eigen::Triplet<double>(v3 * 3 + 1, v2 * 3 + 2, - temp_H(7 , 5)));
        trip_list.push_back(Eigen::Triplet<double>(v3 * 3 + 2, v2 * 3    , - temp_H(8 , 3)));
        trip_list.push_back(Eigen::Triplet<double>(v3 * 3 + 2, v2 * 3 + 1, - temp_H(8 , 4)));
        trip_list.push_back(Eigen::Triplet<double>(v3 * 3 + 2, v2 * 3 + 2, - temp_H(8 , 5)));

        //v3v3
        trip_list.push_back(Eigen::Triplet<double>(v3 * 3    , v3 * 3    , - temp_H(6 , 6)));
        trip_list.push_back(Eigen::Triplet<double>(v3 * 3    , v3 * 3 + 1, - temp_H(6 , 7)));
        trip_list.push_back(Eigen::Triplet<double>(v3 * 3    , v3 * 3 + 2, - temp_H(6 , 8)));
        trip_list.push_back(Eigen::Triplet<double>(v3 * 3 + 1, v3 * 3    , - temp_H(7 , 6)));
        trip_list.push_back(Eigen::Triplet<double>(v3 * 3 + 1, v3 * 3 + 1, - temp_H(7 , 7)));
        trip_list.push_back(Eigen::Triplet<double>(v3 * 3 + 1, v3 * 3 + 2, - temp_H(7 , 8)));
        trip_list.push_back(Eigen::Triplet<double>(v3 * 3 + 2, v3 * 3    , - temp_H(8 , 6)));
        trip_list.push_back(Eigen::Triplet<double>(v3 * 3 + 2, v3 * 3 + 1, - temp_H(8 , 7)));
        trip_list.push_back(Eigen::Triplet<double>(v3 * 3 + 2, v3 * 3 + 2, - temp_H(8 , 8)));
    }
    
    K.resize(q.size(), q.size());
    K.setFromTriplets(trip_list.begin(), trip_list.end());
    
};
