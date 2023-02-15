#include <Eigen/Dense>
#include <EigenTypes.h>
#include <rodrigues.h>
#include <iostream>

//Input:
//  q - 12n vector where n is the number of rigid bodies. Each rigid body is stored as 12 doubles. 
//      The first 9 doubles are the columns of a 3x3 rotation matrix and the final 3 doubles are the world space position of the object's center of mass.
//  qdot - 6n vector of generalied velocities. The first 3 doubles of each body are the world space angular velocity and 
//         the second 3 are the world space linear velocity.
//  dt - the integration time step
//  masses - a vector to mass matrices for each rigid body
//  forces - a 6n vector of generalized forces for n rigid bodies. The first 3 doubles of each rigid body are the torques acting on the object
//           while the second 3 doubles are the linear forces.
//Output:
//  q - updated generalized coordinates 
//  qdot - updated generalized velocities

inline void exponential_euler(Eigen::VectorXd &q, Eigen::VectorXd &qdot, double dt, 
                            std::vector<Eigen::Matrix66d> &masses, Eigen::Ref<const Eigen::VectorXd> forces) {
    
    //q is R and p
    //qdot is w and pdot
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
        
        
        //std::cout << "List of value: q - qdot - I - mass"<<std::endl;
        //std::cout<<q<<std::endl;
        //std::cout<<qdot<<std::endl;
        //std::cout << I << std::endl;
        //std::cout << mass <<std::endl;
        
        //symletic euler to p and pdot
        //m * pdot_t1 = m * pdot_t0 + dt * force
        //p_t1 = p_t0 + dt * pdot_t1
        Eigen::Vector3d pdot_t1 = pdot_t + dt * force / mass;
        Eigen::Vector3d p_t1 = p_t + dt * pdot_t1;
        
        //explicit exponential euler to R and w
        // (R*I*R.T)_t0 * w_t1 = (R*I*R.T)_t0 * w_t0 + dt * (w_t0 x ((R*I*R.T)_t0 * w_t0) + torque)
        // R_t1 = expm([w_t1 * dt]) * R_t0
        Eigen::Matrix3d RIRT = R_t * I * R_t.transpose();
        double  I11 = RIRT(0, 0), I12 = RIRT(0, 1), I13 = RIRT(0, 2),
                I21 = RIRT(1, 0), I22 = RIRT(1, 1), I23 = RIRT(1, 2),
                I31 = RIRT(2, 0), I32 = RIRT(2, 1), I33 = RIRT(2, 2);
        Eigen::Matrix3d inv;
        inv.setZero();
        inv(0, 0) = (I22*I33-I23*I32)/(I11*I22*I33-I11*I23*I32-I12*I21*I33+I12*I23*I31+I13*I21*I32-I13*I22*I31);
        inv(0, 1) = -(I12*I33-I13*I32)/(I11*I22*I33-I11*I23*I32-I12*I21*I33+I12*I23*I31+I13*I21*I32-I13*I22*I31);
        inv(0, 2) = (I12*I23-I13*I22)/(I11*I22*I33-I11*I23*I32-I12*I21*I33+I12*I23*I31+I13*I21*I32-I13*I22*I31);
        inv(1, 0) = -(I21*I33-I23*I31)/(I11*I22*I33-I11*I23*I32-I12*I21*I33+I12*I23*I31+I13*I21*I32-I13*I22*I31);
        inv(1, 1) = (I11*I33-I13*I31)/(I11*I22*I33-I11*I23*I32-I12*I21*I33+I12*I23*I31+I13*I21*I32-I13*I22*I31);
        inv(1, 2) = -(I11*I23-I13*I21)/(I11*I22*I33-I11*I23*I32-I12*I21*I33+I12*I23*I31+I13*I21*I32-I13*I22*I31);
        inv(2, 0) = (I21*I32-I22*I31)/(I11*I22*I33-I11*I23*I32-I12*I21*I33+I12*I23*I31+I13*I21*I32-I13*I22*I31);
        inv(2, 1) = -(I11*I32-I12*I31)/(I11*I22*I33-I11*I23*I32-I12*I21*I33+I12*I23*I31+I13*I21*I32-I13*I22*I31);
        inv(2, 2) = (I11*I22-I12*I21)/(I11*I22*I33-I11*I23*I32-I12*I21*I33+I12*I23*I31+I13*I21*I32-I13*I22*I31);
        
        Eigen::Matrix3d c_w_t;
        c_w_t.setZero();
        c_w_t(1, 0) =  w_t(2);
        c_w_t(2, 0) = -w_t(1);
        c_w_t(0, 1) = -w_t(2);
        c_w_t(2, 1) =  w_t(0);
        c_w_t(0, 2) =  w_t(1);
        c_w_t(1, 2) = -w_t(0);
        
        Eigen::Vector3d b = RIRT * w_t + dt * (c_w_t * RIRT * w_t + torque);

        Eigen::Vector3d w_t1 = inv * b;
        
        Eigen::Matrix3d d_R;
        d_R.setZero();
        rodrigues(d_R, w_t1 * dt);
        
        Eigen::Matrix3d R_t1 = d_R * R_t;
        
        //update q and qdot
        //pdot_t1
        //p_t1
        //w_t1
        //R_t1
        //col major
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
        //std::cout<<q<<std::endl;
        //std::cout<<qdot<<std::endl;
        //exit(0);
    }
}
