#include <V_linear_tetrahedron.h>

#include <dphi_linear_tetrahedron_dX.h>
#include <psi_neo_hookean.h>
#include <quadrature_single_point.h>

//Compute the potential energy of a single tetrahedron. Note: you will need both psi_neo_hookean.cpp and quadrature_single_point.h to do this.

//Input:
// q - generalized coordinates of FEM system
// V - vertex matrix for the mesh (undef vertex n x 3)
// element - vertex indices of the element
// volume - volume of tetrahedron
// C,D - material parameters
//Output:
//  energy - potential energy of tetrahedron

void V_linear_tetrahedron(double &energy, Eigen::Ref<const Eigen::VectorXd> q, 
                          Eigen::Ref<const Eigen::MatrixXd> V, Eigen::Ref<const Eigen::RowVectorXi> element, double volume,
                          double C, double D) {

    auto neohookean_linear_tet = [&](double &e, Eigen::Ref<const Eigen::VectorXd>q, Eigen::Ref<const Eigen::RowVectorXi> element, Eigen::Ref<const Eigen::Vector3d> X) {

    	unsigned  i1 = element(0), i2 = element(1), i3 = element(2), i4 = element(3);
   	 double X_11 = V(i1, 0), X_12 = V(i1, 1), X_13 = V(i1, 2), // X1 - undef vertex1
           	X_21 = V(i2, 0), X_22 = V(i2, 1), X_23 = V(i2, 2), // X2 - undef vertex2
           	X_31 = V(i3, 0), X_32 = V(i3, 1), X_33 = V(i3, 2), // X3 - undef vertex3
           	X_41 = V(i4, 0), X_42 = V(i4, 1), X_43 = V(i4, 2), // X4 - undef vertex4
           	  X1 = X(0), X2 = X(1), X3 = X(2); // X - undef vertex

	double x_11 = q(3 * i1), x_12 = q(3 * i1 + 1), x_13 = q(3 * i1 + 2),
               x_21 = q(3 * i2), x_22 = q(3 * i2 + 1), x_23 = q(3 * i2 + 2),
               x_31 = q(3 * i3), x_32 = q(3 * i3 + 1), x_33 = q(3 * i3 + 2),
               x_41 = q(3 * i4), x_42 = q(3 * i4 + 1), x_43 = q(3 * i4 + 2);
	
	double psi = 0;
	
	Eigen::Matrix43d dphi;
	dphi.setZero();
	dphi_linear_tetrahedron_dX(dphi, V, element, X);
	
	Eigen::Matrix3d F;
	F.setZero();

	F(0, 0) = dphi(0,0) * x_11 + dphi(1,0) * x_21 + dphi(2, 0) * x_31 + dphi(3, 0) * x_41;
    
    	F(0, 1) = dphi(0,1) * x_11 + dphi(1,1) * x_21 + dphi(2, 1) * x_31 + dphi(3, 1) * x_41;
    
    	F(0, 2) = dphi(0,2) * x_11 + dphi(1,2) * x_21 + dphi(2, 2) * x_31 + dphi(3, 2) * x_41;
    
    	F(1, 0) = dphi(0,0) * x_12 + dphi(1,0) * x_22 + dphi(2, 0) * x_32 + dphi(3, 0) * x_42;
                             
    	F(1, 1) = dphi(0,1) * x_12 + dphi(1,1) * x_22 + dphi(2, 1) * x_32 + dphi(3, 1) * x_42;
                             
    	F(1, 2) = dphi(0,2) * x_12 + dphi(1,2) * x_22 + dphi(2, 2) * x_32 + dphi(3, 2) * x_42;
    
    	F(2, 0) = dphi(0,0) * x_13 + dphi(1,0) * x_23 + dphi(2, 0) * x_33 + dphi(3, 0) * x_43;
                             
    	F(2, 1) = dphi(0,1) * x_13 + dphi(1,1) * x_23 + dphi(2, 1) * x_33 + dphi(3, 1) * x_43;
                             
    	F(2, 2) = dphi(0,2) * x_13 + dphi(1,2) * x_23 + dphi(2, 2) * x_33 + dphi(3, 2) * x_43;



		
        psi_neo_hookean(psi, F, C, D);
	e = psi;
       
    };

    quadrature_single_point(energy, q, element, volume, neohookean_linear_tet);  
    
}
