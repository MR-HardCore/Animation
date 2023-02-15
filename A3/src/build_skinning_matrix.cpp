#include <build_skinning_matrix.h>
#include <phi_linear_tetrahedron.h>
#include <vector>
#include <iostream>

//Build the skinning matrix that maps position from the coarse simulation mesh to the high resolution rendering mesh.

//Input:
//  V - the nx3 matrix of undeformed vertex positions. Each row is a single undeformed vertex position.
//  T - the mx4 tetrahedron connectivity matrix. Each row contains to indices into V that indicate a spring between those vertices.
//  V_skin - lx3 matrix of vertices of the display mesh
//Output:
//  N - the lxn sparse skinning matrix 
void build_skinning_matrix(Eigen::SparseMatrixd &N, Eigen::Ref<const Eigen::MatrixXd> V, Eigen::Ref<const Eigen::MatrixXi> T, 
                                                   Eigen::Ref<const Eigen::MatrixXd> V_skin) {
    int l = V_skin.rows(), n = V.rows(), m = T.rows();
    std::vector<Eigen::Triplet<double>> trip_list;
    for (unsigned li = 0; li < l; li++)
    {
	for (unsigned mi = 0; mi < m; mi++)
	{
	    int v1 = T(mi, 0), v2 = T(mi, 1), v3 = T(mi, 2), v4 = T(mi, 3);
	    Eigen::Vector4d phi;
	    phi_linear_tetrahedron(phi, V, Eigen::RowVector4i(v1, v2, v3, v4), Eigen::Vector3d(V_skin(li, 0), V_skin(li, 1), V_skin(li, 2)));
	    if (phi(0) >= 0 && phi(1) >= 0 && phi(2) >= 0 && phi(3) >= 0) // inside current V tet
	    {
		trip_list.push_back(Eigen::Triplet<double>(li , v1, phi(0)));
		trip_list.push_back(Eigen::Triplet<double>(li , v2, phi(1)));
		trip_list.push_back(Eigen::Triplet<double>(li , v3, phi(2)));
		trip_list.push_back(Eigen::Triplet<double>(li , v4, phi(3)));
		//std::cout<<"Check the phi value: phi1: "<<phi(0)<<" phi2: "<<phi(1)<<" phi3: "<<phi(2)<<" phi4: "<<phi(3)<<" the right phi4 should be: "<<(1-phi(0)-phi(1)-phi(2))<<std::endl;
		break;
	    }
	}
    }
    N.resize(l, n);
    N.setFromTriplets(trip_list.begin(), trip_list.end());
}
