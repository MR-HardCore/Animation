#include <Eigen/Dense>
#include <EigenTypes.h>
#include <rodrigues.h>
#include <iostream>
#include <rigid_body_jacobian.h>
#include <inverse_rigid_body.h>
#include <iostream>

//Implement velocity level collision resolution using progressive Gauss-Seidel and exponential Euler time integration.


//Input:
//  q - 12n vector where n is the number of rigid bodies. Each rigid body is stored as 12 doubles. 
//      The first 9 doubles are the columns of a 3x3 rotation matrix and the final 3 doubles are the world space position of the object's center of mass.
//  qdot - 6n vector of generalied velocities. The first 3 doubles of each body are the world space angular velocity and 
//         the second 3 are the world space linear velocity.
//  dt - the integration time step
//  masses - a vector to mass matrices for each rigid body
//  forces - a 6n vector of generalized forces for n rigid bodies. The first 3 doubles of each rigid body are the torques acting on the object
//           while the second 3 doubles are the linear forces.
//  n - list of collision normals
//  x - list of world space collision points
//  obj - list of collision object ids 
//Output:
//  q - updated generalized coordinates 
//  qdot - updated generalized velocities 
inline void exponential_euler_lcp_contact(Eigen::VectorXd &q, Eigen::VectorXd &qdot, double dt, 
                            std::vector<Eigen::Matrix66d> &masses, Eigen::Ref<const Eigen::VectorXd> forces,
                            std::vector<Eigen::Vector3d> &n, std::vector<Eigen::Vector3d> &x, std::vector<std::pair<int,int> > &obj) {
    
    
    //G is 3x6 rigid_body_jacobian
    int n_rb = masses.size();
    for (unsigned idx = 0; idx < n_rb; idx++)
    {
        double mass = masses[idx](3,3);
        Eigen::Vector3d torque = Eigen::Map<const Eigen::Vector3d>(forces.segment<3>(6*idx).data());
        Eigen::Vector3d force = Eigen::Map<const Eigen::Vector3d>(forces.segment<3>(6*idx + 3).data());
        Eigen::Matrix3d R_t = Eigen::Map<const Eigen::Matrix3d>(q.segment<9>(12*idx).data());
        Eigen::Vector3d p_t = Eigen::Map<const Eigen::Vector3d>(q.segment<3>(12*idx + 9).data());
        Eigen::Vector3d w_t = Eigen::Map<const Eigen::Vector3d>(qdot.segment<3>(6*idx).data());
        Eigen::Vector3d pdot_t = Eigen::Map<const Eigen::Vector3d>(qdot.segment<3>(6*idx + 3).data());
        Eigen::Matrix3d I = masses[idx].block(0, 0, 3, 3);
        
        Eigen::Matrix3d RIRT = R_t * I * R_t.transpose();
        
        Eigen::Matrix3d c_w_t;
        c_w_t.setZero();
        c_w_t(1, 0) =  w_t(2);
        c_w_t(2, 0) = -w_t(1);
        c_w_t(0, 1) = -w_t(2);
        c_w_t(2, 1) =  w_t(0);
        c_w_t(0, 2) =  w_t(1);
        c_w_t(1, 2) = -w_t(0);
        
        Eigen::Matrix66d M = masses[idx];
        M.block(0,0,3,3) = RIRT;
        
        Eigen::Vector6d qdot_t = Eigen::Map<const Eigen::Vector6d>(qdot.segment<6>(6*idx).data());
        
        Eigen::Vector3d tor_t = c_w_t * RIRT * w_t + torque;
        Eigen::Vector6d for_t;
        for_t(0) = tor_t(0);
        for_t(1) = tor_t(1);
        for_t(2) = tor_t(2);
        for_t(3) = force(0);
        for_t(4) = force(1);
        for_t(5) = force(2);
        
        Eigen::Vector6d F;
        F = M * qdot_t + dt * for_t;

        double  I11 = RIRT(0, 0), I12 = RIRT(0, 1), I13 = RIRT(0, 2),
                I21 = RIRT(1, 0), I22 = RIRT(1, 1), I23 = RIRT(1, 2),
                I31 = RIRT(2, 0), I32 = RIRT(2, 1), I33 = RIRT(2, 2);
        Eigen::Matrix66d M_inv;
        M_inv.setZero();
        M_inv(0, 0) = (I22*I33-I23*I32)/(I11*I22*I33-I11*I23*I32-I12*I21*I33+I12*I23*I31+I13*I21*I32-I13*I22*I31);
        M_inv(0, 1) = -(I12*I33-I13*I32)/(I11*I22*I33-I11*I23*I32-I12*I21*I33+I12*I23*I31+I13*I21*I32-I13*I22*I31);
        M_inv(0, 2) = (I12*I23-I13*I22)/(I11*I22*I33-I11*I23*I32-I12*I21*I33+I12*I23*I31+I13*I21*I32-I13*I22*I31);
        M_inv(1, 0) = -(I21*I33-I23*I31)/(I11*I22*I33-I11*I23*I32-I12*I21*I33+I12*I23*I31+I13*I21*I32-I13*I22*I31);
        M_inv(1, 1) = (I11*I33-I13*I31)/(I11*I22*I33-I11*I23*I32-I12*I21*I33+I12*I23*I31+I13*I21*I32-I13*I22*I31);
        M_inv(1, 2) = -(I11*I23-I13*I21)/(I11*I22*I33-I11*I23*I32-I12*I21*I33+I12*I23*I31+I13*I21*I32-I13*I22*I31);
        M_inv(2, 0) = (I21*I32-I22*I31)/(I11*I22*I33-I11*I23*I32-I12*I21*I33+I12*I23*I31+I13*I21*I32-I13*I22*I31);
        M_inv(2, 1) = -(I11*I32-I12*I31)/(I11*I22*I33-I11*I23*I32-I12*I21*I33+I12*I23*I31+I13*I21*I32-I13*I22*I31);
        M_inv(2, 2) = (I11*I22-I12*I21)/(I11*I22*I33-I11*I23*I32-I12*I21*I33+I12*I23*I31+I13*I21*I32-I13*I22*I31);
        M_inv(3, 3) = 1.0/mass;
        M_inv(4, 4) = 1.0/mass;
        M_inv(5, 5) = 1.0/mass;
        
        
        Eigen::Vector6d qdot_unc = M_inv * F;
        
        //Gauss Sidel Method
        //maximum 10 loops for alpha
        int n_c = n.size();
        //initilize alpha to 0
        std::vector<double> alpha(n_c, 0);
        for (unsigned iter = 0; iter < 10; iter++)
        {
            for(unsigned c = 0; c < n_c; c++)
            {
                Eigen::Vector3d c_normal = n[c];
                Eigen::Vector3d c_x = x[c];
                Eigen::Vector3d inv_x;
                inv_x.setZero();
                inverse_rigid_body(inv_x, c_x, R_t, p_t);
                Eigen::Matrix36d G;
                G.setZero();
                rigid_body_jacobian(G, R_t, p_t, inv_x);
                
                //double delta_c = dt * (c_normal.transpose() * G * M_inv * G.transpose() * c_normal);
                double delta_c = dt * (G * M_inv * G.transpose() * c_normal).transpose() * c_normal;
                
                Eigen::Vector6d b_n1, b_p1;
                b_n1.setZero();
                b_p1.setZero();
                
                //b i - 1
                for (unsigned b = 0; b < c; b++)
                {
                    Eigen::Vector3d b_normal = n[b];
                    Eigen::Vector3d b_x = x[b];
                    Eigen::Vector3d inv_b;
                    inv_b.setZero();
                    inverse_rigid_body(inv_b, b_x, R_t, p_t);
                    Eigen::Matrix36d b_G;
                    b_G.setZero();
                    rigid_body_jacobian(b_G, R_t, p_t, inv_x);
                    
                    b_n1 += b_G.transpose() * b_normal * alpha[b];
                }
                
                //b i + 1
                for (unsigned b = c + 1; b < n_c; b++)
                {
                    Eigen::Vector3d b_normal = n[b];
                    Eigen::Vector3d b_x = x[b];
                    Eigen::Vector3d inv_b;
                    inv_b.setZero();
                    inverse_rigid_body(inv_b, b_x, R_t, p_t);
                    Eigen::Matrix36d b_G;
                    b_G.setZero();
                    rigid_body_jacobian(b_G, R_t, p_t, inv_x);
                    
                    b_p1 += b_G.transpose() * b_normal * alpha[b];
                }
                
                b_n1 = dt * M_inv * b_n1;
                b_p1 = dt * M_inv * b_p1;
                
                //double gamma_c = c_normal.transpose() * G * (qdot_unc + b_n1 + b_p1);
                double gamma_c = (G * (qdot_unc + b_n1 + b_p1)).transpose() * c_normal;
                //std::cout<<delta_c<<std::endl;
                //std::cout<<gamma_c<<std::endl;
                
                //alpha[c] = max(0, -gamma_c/delta_c);
                alpha[c] = (-gamma_c / delta_c) < 0? 0:(-gamma_c/delta_c);
            }
        }
        
        //This is an optimized update by storing lots of reused term
        //But Variable B will overflow once iteration time larger than 3
        //Can function well too
        /*
        std::vector<double> alpha(n_c, 0);
        Eigen::Vector6d B;
        B.setZero();
        std::vector<Eigen::Matrix36d> J;
        for (unsigned c = 0; c < n_c; c++)
        {
            Eigen::Matrix36d G;
            G.setZero();
            Eigen::Vector3d X_c;
            inverse_rigid_body(X_c, x[c], R_t, p_t);
            rigid_body_jacobian(G, R_t, p_t, X_c);
            J.push_back(G);
        }
        for(unsigned iter = 0; iter < 2; iter++ )
        {
            for (unsigned c = 0; c < n_c; c++)
            {
                Eigen::Vector6d b_term = B - J[c].transpose() * n[c] * alpha[c];
                //double gamma_c = n[c].transpose() * J[c] * (qdot_unc + dt * M_inv * b_term);
                double gamma_c = (J[c] * (qdot_unc + dt * M_inv * b_term)).transpose() * n[c];
                //double delta_c = dt * n[c].transpose() * J[c] * M_inv * J[c].transpose() *n[c];
                double delta_c = dt * (J[c] * M_inv * J[c].transpose() *n[c]).transpose() * n[c];
                alpha[c] = (-gamma_c / delta_c) < 0? 0:(-gamma_c/delta_c);
                B += J[c].transpose() * n[c] * alpha[c];
            }
        }
        */
        
        Eigen::Vector6d ext;
        ext.setZero();
        for(unsigned c = 0; c < n_c; c++)
        {
            Eigen::Vector3d c_normal = n[c];
            Eigen::Vector3d c_x = x[c];
            Eigen::Vector3d inv_x;
            inv_x.setZero();
            inverse_rigid_body(inv_x, c_x, R_t, p_t);
            Eigen::Matrix36d G;
            G.setZero();
            rigid_body_jacobian(G, R_t, p_t, inv_x);
            
            ext += G.transpose() * c_normal * alpha[c];
        }
        
        Eigen::Vector6d qdot_t1 = qdot_unc + dt * M_inv * ext;
        
        Eigen::Vector3d w_t1;
        w_t1(0) = qdot_t1(0);
        w_t1(1) = qdot_t1(1);
        w_t1(2) = qdot_t1(2);
        
        Eigen::Vector3d pdot_t1;
        pdot_t1(0) = qdot_t1(3);
        pdot_t1(1) = qdot_t1(4);
        pdot_t1(2) = qdot_t1(5);
        
        Eigen::Vector3d p_t1 = p_t + dt * pdot_t1;
        
        Eigen::Matrix3d d_R;
        d_R.setZero();
        rodrigues(d_R, w_t1 * dt);
        
        Eigen::Matrix3d R_t1 = d_R * R_t;
        
        q(idx * 12    )  = R_t1(0, 0);
        q(idx * 12 + 1)  = R_t1(1, 0);
        q(idx * 12 + 2)  = R_t1(2, 0);
        q(idx * 12 + 3)  = R_t1(0, 1);
        q(idx * 12 + 4)  = R_t1(1, 1);
        q(idx * 12 + 5)  = R_t1(2, 1);
        q(idx * 12 + 6)  = R_t1(0, 2);
        q(idx * 12 + 7)  = R_t1(1, 2);
        q(idx * 12 + 8)  = R_t1(2, 2);
        
        q(idx * 12 + 9)  = p_t1(0);
        q(idx * 12 + 10) = p_t1(1);
        q(idx * 12 + 11) = p_t1(2);
        
        qdot(idx * 6 + 0) = w_t1(0);
        qdot(idx * 6 + 1) = w_t1(1);
        qdot(idx * 6 + 2) = w_t1(2);
        qdot(idx * 6 + 3) = pdot_t1(0);
        qdot(idx * 6 + 4) = pdot_t1(1);
        qdot(idx * 6 + 5) = pdot_t1(2);
    }
    
}
