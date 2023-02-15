#include <velocity_filter_cloth_sphere.h>

//Project out components of the per-vertex velocities which are in the positive direction of the contact normal

//Input:
//  qdot - the 3nx1 generalized velocities of the cloth mesh
//  index - a list of collision vertex indices from the collision detector
//  normals - a list of collision normals from the collision detector
//Output:
//  qdot- the filtered 3nx1 generalized velocities
void velocity_filter_cloth_sphere(Eigen::VectorXd &qdot, const std::vector<unsigned int> &indices, 
                                  const std::vector<Eigen::Vector3d> &normals) {
    
    int n_collision = indices.size();

    for (unsigned idx = 0; idx < n_collision; idx++)
    {
        Eigen::Vector3d qdot_i(qdot(indices[idx] * 3 + 0) , qdot(indices[idx] * 3 + 1) , qdot(indices[idx] * 3 + 2)); 
        double alpha = normals[idx].transpose() * qdot_i;
        if (alpha < 0)
        {
            qdot(indices[idx] * 3 + 0) -= alpha * normals[idx](0);
            qdot(indices[idx] * 3 + 1) -= alpha * normals[idx](1);
            qdot(indices[idx] * 3 + 2) -= alpha * normals[idx](2);
        }
    }
    

}
