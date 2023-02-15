#include <pick_nearest_vertices.h>
#include <iostream>
#include <igl/Hit.h>
#include <igl/ray_mesh_intersect.h>
#include <igl/unproject.h>
#include <algorithm>

bool pick_nearest_vertices(std::vector<unsigned int> &verts, Eigen::Ref<const Eigen::Vector3d> win, 
                           Eigen::Ref<const Eigen::Matrix44f> view, Eigen::Ref<const Eigen::Matrix44f> proj, Eigen::Vector4f viewport,
                           Eigen::Ref<const Eigen::MatrixXd> V, Eigen::Ref<const Eigen::MatrixXi> F, double radius) {
//Input:
//  win - window coordinate of mouse click (x_window, y_window, 0)
//  view - view transformation matrix
//  proj - projection matrix
//  viewport - viewport coordinates .
//  V - 3xn dense matrix of mesh vertices, each row of the matrix is a single vertex.
//  radius - selection radius for vertex picking
//Output:
//  verts - vertex ids (rows in V) of selected vertices

    verts.clear();

    // Source, destination and direction in world
    Eigen::Vector3f start,dir;
    igl::Hit hit;

    //compute start and direction in the world to check for picked vertex
    //YOUR CODE HERE
    //TODO:make unproject work
    // unproj(win     ,model    ,proj     ,viewport,scene)
    // para:  Vector3f Matrix44f Matrix44f Vector4f Vector3f 

    start(0) = float(win(0));
    start(1) = float(win(1));
    start(2) = float(win(2));

    Eigen::Vector3f unpro_win, unpro_src;
    //unproject the windows coordinate
    igl::unproject(start, view, proj, viewport, unpro_win);
    
    //unproject the rt src coordinate
    start(2) = 1.0f;
    igl::unproject(start, view, proj, viewport, unpro_src);

    //dir = win - src
    dir(0) = unpro_win(0) - unpro_src(0);
    dir(1) = unpro_win(1) - unpro_src(1);
    dir(2) = unpro_win(2) - unpro_src(2);
    
    /*
    //for debug use only: 
    std::cout<<"RT mouse coordinates: X-" <<win(0)<<" Y-"<<win(1) <<" Z-"<<win(2)<<std::endl;
    std::cout<<"RT dir coordinates: X: " <<dir(0)<<" Y: "<<dir(1) <<" Z: "<<dir(2)<<std::endl;
    std::cout<<"RT unpro_win: X: " <<unpro_win(0)<<" Y: "<<unpro_win(1) <<" Z: "<<unpro_win(2)<<std::endl;
    std::cout<<"RT unpro_src: X: " <<unpro_src(0)<<" Y: "<<unpro_src(1) <<" Z: "<<unpro_src(2)<<std::endl;
    */
    
    //set the src as the world coordinate
    start(0) = float(unpro_src(0));
    start(1) = float(unpro_src(1));
    start(2) = float(unpro_src(2));

    const auto & shoot_ray = [&V,&F](const Eigen::Vector3f& s, const Eigen::Vector3f& dir, igl::Hit & hit)->bool
    {
        std::vector<igl::Hit> hits;
        
        if(!igl::ray_mesh_intersect(s,dir,V,F,hits))
        {
            return false;
        }
        hit = hits[0];
        return true;
    };

    if(!shoot_ray(start,dir,hit))
    {
        return false;
    }
    //check if any of the hit vertices are within the selection radius 
    //YOUR CODE HERE

    //TODO:
    //id - primitive  uv - barycentric coordinates
    float lambda1 = hit.u, lambda2 = hit.v, lambda3 = 1 - hit.u - hit.v; 
    int target_id = -1;

    if (lambda1 > lambda2) 
    { 
	if (lambda1 > lambda3) //lambda1 largest -> choose v0
	{ 
	    target_id = F(hit.id, 0); 
	} 
	else //lambda3 largest -> choose v2
	{
	    target_id = F(hit.id, 2); 
	} 
    }
    else 
    { 
	if (lambda2 > lambda3) //lambda2 largest -> choose v1
	{
	    target_id = F(hit.id, 1); 
	} 
	else //lambda3 largest -> choose v2
	{
	    target_id = F(hit.id, 2); 
	} 
    }
    if (target_id == -1 ) std::cout<<"Something wrong with target_id"<<std::endl;

    //check all vectices within the radius in world space
    //NOT Conical space
    double l2_constraint = radius * radius;
    Eigen::Vector3f w_target, v_target(V(target_id, 0), V(target_id, 1), V(target_id, 2));
    igl::unproject(v_target, view, proj, viewport, w_target);
    /*
    //for debug use:
    std::cout<<"Size of F: "<<F.rows() <<" "<<F.size()<< " target_id is: "<<target_id<<" Size of V is: "<<V.rows() << " "<< V.size()<<std::endl;
    std::cout<<"GID: "<<hit.gid << " ID? "<<hit.id<<" l1: "<<lambda1<<" l2: "<<lambda2<<" l3: "<<lambda3 << " R: "<< radius <<std::endl;  
    std::cout<<"target vertex v-space: "<<" x: "<<V(target_id, 0)<<" y: "<<V(target_id, 1)<<" z: "<<V(target_id, 2) <<std::endl;
    std::cout<<"target vertex v-space: "<<" x: "<<v_target(0)<<" y: "<<v_target(1)<<" z: "<<v_target(2) <<std::endl;
    std::cout<<"target vertex world-space: "<<" x: "<<w_target(0)<<" y: "<<w_target(1)<<" z: "<<w_target(2) <<std::endl;
    */
    for (unsigned idx = 0; idx < V.rows(); idx ++) 
    {   
        Eigen::Vector3f w_cur, v_cur(V(idx, 0), V(idx, 1), V(idx, 2));
	igl::unproject(v_cur, view, proj, viewport, w_cur);
	double l2 = (w_target(0) - w_cur(0)) * (w_target(0) - w_cur(0)) + 
		    (w_target(1) - w_cur(1)) * (w_target(1) - w_cur(1)) + 
		    (w_target(2) - w_cur(2)) * (w_target(2) - w_cur(2));
	if (l2 <= l2_constraint) {verts.push_back(idx);} 
	// for debug use: std::cout<<"near vertex world-space: "<<" x: "<<w_cur(0)<<" y: "<<w_cur(1)<<" z: "<<w_cur(2) <<std::endl;
    }
    /*
    // for debug use:
    std::cout<<"size of affected vertex: "<<verts.size()  <<std::endl;
    for (unsigned i = 0; i < verts.size(); i++){std::cout<<verts[i]<<std::endl;}*/
    return (verts.size() == 0 ? false : true);
}
