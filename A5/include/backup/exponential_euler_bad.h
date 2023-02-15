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
        
        Eigen::Matrix3d RIRT = R_t*I*R_t.transpose();
        Eigen::Matrix3d mI;
        mI.setZero();
        mI(0, 0) = mass;
        mI(1, 1) = mass;
        mI(2, 2) = mass;
        
        Eigen::Vector6d qdot_t;
        qdot_t(0) = w_t(0);
        qdot_t(1) = w_t(1);
        qdot_t(2) = w_t(2);
        qdot_t(3) = pdot_t(0);
        qdot_t(4) = pdot_t(1);
        qdot_t(5) = pdot_t(2);
        
        
        Eigen::SparseMatrixd A_M;
        A_M.resize(6, 6);
        A_M.setZero();
        //A_M.block(0, 0, 3, 3) = RIRT;
        //A_M.block(3, 3, 6, 6) = mI;
        A_M.coeffRef(0, 0) = RIRT(0, 0);
        A_M.coeffRef(0, 1) = RIRT(0, 1);
        A_M.coeffRef(0, 2) = RIRT(0, 2);
        A_M.coeffRef(1, 0) = RIRT(1, 0);
        A_M.coeffRef(1, 1) = RIRT(1, 1);
        A_M.coeffRef(1, 2) = RIRT(1, 2);
        A_M.coeffRef(2, 0) = RIRT(2, 0);
        A_M.coeffRef(2, 1) = RIRT(2, 1);
        A_M.coeffRef(2, 2) = RIRT(2, 2);
        
        A_M.coeffRef(3, 3) = mass;
        A_M.coeffRef(4, 4) = mass;
        A_M.coeffRef(5, 5) = mass;
        
        
        Eigen::Matrix3d c_w_t;
        c_w_t.setZero();
        c_w_t(1, 0) =  w_t(2);
        c_w_t(2, 0) = -w_t(1);
        c_w_t(0, 1) = -w_t(2);
        c_w_t(2, 1) =  w_t(0);
        c_w_t(0, 2) =  w_t(1);
        c_w_t(1, 2) = -w_t(0);
        
        Eigen::Vector3d tor_t = c_w_t * RIRT * w_t + torque;
        
        Eigen::Vector6d b;
        b(0) = tor_t(0);
        b(1) = tor_t(1);
        b(2) = tor_t(2);
        b(3) = force(0);
        b(4) = force(1);
        b(5) = force(2);
        
        b = b * dt + A_M * qdot_t;

        Eigen::SimplicialLDLT<Eigen::SparseMatrixd> solverLDLT;
        solverLDLT.compute(A_M);
        Eigen::Vector6d qdot_t1 = solverLDLT.solve(b);
        
        Eigen::Vector3d w_t1;
        Eigen::Vector3d pdot_t1;
        w_t1(0) = qdot_t1(0);
        w_t1(1) = qdot_t1(1);
        w_t1(2) = qdot_t1(2);
        pdot_t1(0) = qdot_t1(3);
        pdot_t1(1) = qdot_t1(4);
        pdot_t1(2) = qdot_t1(5);
        
        //symletic euler to p and pdot
        //m * pdot_t1 = m * pdot_t0 + dt * force
        //p_t1 = p_t0 + dt * pdot_t1
        Eigen::Vector3d p_t1 = p_t + dt * pdot_t1;
        
        Eigen::Matrix3d d_R;
        d_R.setZero();
        rodrigues(d_R, w_t1 * dt);
        
        Eigen::Matrix3d R_t1 = d_R * R_t;
        
        
        std::cout << "List of value: q - qdot - I - mass"<<std::endl;
        std::cout<<q<<std::endl;
        std::cout<<qdot<<std::endl;
        std::cout << I << std::endl;
        std::cout << mass <<std::endl;
        
        

        //explicit exponential euler to R and w
        // (R*I*R.T)_t0 * w_t1 = (R*I*R.T)_t0 * w_t0 + dt * (w_t0 x ((R*I*R.T)_t0 * w_t0) + torque)
        // R_t1 = expm([w_t1 * dt]) * R_t0
    

        
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
        
        std::cout<<"After iter: "<<std::endl;
        std::cout<<q<<std::endl;
        std::cout<<qdot<<std::endl;
        exit(0);
    }
}
