#include <collision_detection_cloth_sphere.h>
#include <iostream>

//Detect if any mesh vertex falls inside a sphere centered at (0,0,0.4) with radius 0.22

//  q - generalized coordinates for the FEM system
//  center - the position of the sphere center in the world space
//  radius - the radius of the sphere in the world space 
//Output:
//  cloth_index - the indices of all vertices currently in contact with the sphere
//  normals - the outward facing contact normals for each contacting vertex. 
void collision_detection_cloth_sphere(std::vector<unsigned int> &cloth_index, std::vector<Eigen::Vector3d> &normals, Eigen::Ref<const Eigen::VectorXd> q, Eigen::Ref<const Eigen::Vector3d> center, double radius) {

    cloth_index.clear();
    normals.clear();
    
    unsigned n_vertex = q.size()/3;
    for(unsigned idx = 0; idx < n_vertex; idx++)
    {
        Eigen::Vector3d v(q(idx * 3 + 0), q(idx * 3 + 1), q(idx * 3 + 2));
        double dist = sqrt((v(0) - center(0)) * (v(0) - center(0)) + (v(1) - center(1)) * (v(1) - center(1)) + (v(2) - center(2)) * (v(2) - center(2)));
        if (dist <= radius) 
        {
            cloth_index.push_back(idx);
            normals.push_back(Eigen::Vector3d( (v(0) - center(0)) / dist, (v(1) - center(1)) / dist, (v(2) - center(2)) / dist));
        }
    }     
}
