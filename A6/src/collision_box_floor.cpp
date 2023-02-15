#include <collision_box_floor.h>

//Detect contact between a triangle mesh and an arbitrarily positioned plane.

//Input:
//  R - rotation matrix for rigid body
//  p - world space position of center-of-mass
//  obj_id - id of object being checked for collision with the floor
//  V - the nx3 matrix of undeformed vertex positions
//  dir - the outward facing normal for the floor
//  pos - the world space position of a point on the floor plane
//Output:
//  x - world space, per-vertex collision points. Computed as any vertex that is on the "wrong side" of the floor plane
//  n - collision normals, one for each collision point. These point away from the floor.
//  objs - Pairs of ids for objects involved in collisions. The first id is for the object, away from which the normal points. The second id
//  is for the object towards which the normal points. The floor has an id of -1.
void collision_box_floor(std::vector<Eigen::Vector3d> &x, std::vector<Eigen::Vector3d> &n, std::vector<std::pair<int,int>> &objs,
                         Eigen::Ref<const Eigen::Matrix3d> R, Eigen::Ref<const Eigen::Vector3d> p, unsigned int obj_id,
                         Eigen::Ref<Eigen::MatrixXd> V, 
                         Eigen::Ref<const Eigen::Vector3d> dir, Eigen::Ref<const Eigen::Vector3d> pos) {

    //rigid body obj
    //R: rotation matrix
    //p: position vector com
    //obj_id: id of obj
    //V: undef space of obj
    
    //floor obj
    //dir: the outward normal
    //pos is point
    
    //output:
    //x: world space postion for collision points
    //n: normal for all collision points
    //objs: first - away from normal points | second - towards the normal | floor is -1
    
    // signed distance = normal.t * (p_object - p_floor)
    unsigned n_vertex = V.rows();
    //for each of points on current obj
    for (unsigned idx = 0; idx < n_vertex; idx ++)
    {
        Eigen::Vector3d cur_X = V.row(idx);
        Eigen::Vector3d cur_x = R * cur_X + p;
        double temp_sd = dir.transpose() * (cur_x - pos);
        if (temp_sd <= 0) //collision happened //zero boundary might be significant
        {
            x.push_back(cur_x);
            n.push_back(dir);
            objs.push_back(std::make_pair(obj_id, -1)); // might be flipped
        }
    }

}
