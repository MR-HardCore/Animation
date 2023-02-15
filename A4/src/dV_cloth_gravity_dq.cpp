#include <dV_cloth_gravity_dq.h>

//Gradient of potential energy due to gravity

//  M - sparse mass matrix for the entire mesh
//  g - the acceleration due to gravity
//Output:
//  fg - the gradient of the gravitational potential for the entire mesh

void dV_cloth_gravity_dq(Eigen::VectorXd &fg, Eigen::SparseMatrixd &M, Eigen::Ref<const Eigen::Vector3d> g) {
    
    int l_fg = M.rows();
    fg.resize(l_fg);
    fg.setZero();
    for (unsigned idx = 0; idx < l_fg; idx++)
    {
        if ( idx % 3 == 0 ) fg(idx) =  - g(0);
        if ( idx % 3 == 1 ) fg(idx) =  - g(1);
        if ( idx % 3 == 2 ) fg(idx) =  - g(2);
    }
    fg = M * fg;
    //fg *= 0;
}
