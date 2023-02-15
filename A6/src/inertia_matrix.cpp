#include <inertia_matrix.h>
#include <cassert>
#include <iostream>

//compute inertia matrix and volume by integrating on surfaces

//Compute the rigid inertia matrix of a 3d object, represented by a surface mesh, via surface only integration.

//Input:
//  V - the nx3 matrix of vertices.
//  F - the mx3 matrix of triangle vertex indices.
//  density - the material density.
//Output:
//  I - the 3x3 angular inertia matrix
//  center - the center of mass of the object
//  mass - the total mass of the object

void integrate(const double w0, const double w1, const double w2, double &f1, double &f2, double &f3, double &g0, double &g1, double &g2)
{
    double temp0 = w0 + w1;
    f1 = temp0 + w2;
    double temp1 = w0 * w0;
    double temp2 = temp1 + w1 * temp0;
    f2 = temp2 + w2 * f1;
    f3 = w0 * temp1 + w1 * temp2 + w2 * f2;
    g0 = f2 + w0 * (f1 + w0);
    g1 = f2 + w1 * (f1 + w1);
    g2 = f2 + w2 * (f1 + w2);
    //std::cout<<f1<<" "<<f2<<" "<<f3<<" "<<g0<<" "<<g1<<" "<<g2<<std::endl;
}


void inertia_matrix(Eigen::Matrix3d &I, Eigen::Vector3d & center, double &mass, Eigen::Ref<const Eigen::MatrixXd> V, Eigen::Ref<const Eigen::MatrixXi> F, double density) {
    
    //mass stands for right bottom entry of M0
    
    //I stands for the left top entry of M0
    
    //center is the com
    
    //density is constant
    center.setZero();
    mass = 0;
    I.setZero();
    unsigned n_triangle = F.rows();
    /*
    for (unsigned idx = 0; idx < n_triangle; idx++)
    {

        unsigned v1 = F(idx, 0), v2 = F(idx, 1), v3 = F(idx, 2);
        double x0 = V(v1, 0), x1 = V(v2, 0), x2 = V(v3, 0),
               y0 = V(v1, 1), y1 = V(v2, 1), y2 = V(v3, 1),
               z0 = V(v1, 2), z1 = V(v2, 2), z2 = V(v3, 2);
        
        //center of mass for this triangle
        Eigen::Vector3d com (x0 + x1 + x2, y0 + y1 + y2, z0 + z1 + z2);
        com = com / 3.0;
        
        //useless information
        //AB: v2 - v1
        //Eigen::Vector3d AB(x1 - x0, y1 - y0, z1 - z0);
        //AC: v3 - v1
        //Eigen::Vector3d AC(x2 - x0, y2 - y0, z2 - z0);
        //Eigen::Matrix3d m_AB <<0,  AB(2), -AB(1),
        //                  -AB(2),      0,  AB(0),
        //                   AB(1), -AB(0),      0;
        //cross product: AB X AC
        //Eigen::Vector3d res = m_AB * AC;
        //double area = 0.5 * sqrt(res(0) * res(0) + res(1) * res(1) + res(2) * res(2));
        
        //addition in center & mass
        //center += area * com;
        //t_area += area;
        //mass += area * density;
        
        double X = x0/6.0+x1/6.0+x2/6.0;
        
        double XXX  = ((x0*x0*x0)*y0)/3.0E+1+((x0*x0*x0)*y1)/1.2E+2+((x1*x1*x1)*y0)/1.2E+2+((x0*x0*x0)*y2)/1.2E+2+((x1*x1*x1)*y1)/3.0E+1+((x2*x2*x2)*y0)/1.2E+2+((x1*x1*x1)*y2)/1.2E+2+((x2*x2*x2)*y1)/1.2E+2+((x2*x2*x2)*y2)/3.0E+1+(x0*(x1*x1)*y0)/6.0E+1+((x0*x0)*x1*y0)/4.0E+1+(x0*(x1*x1)*y1)/4.0E+1+(x0*(x2*x2)*y0)/6.0E+1+((x0*x0)*x1*y1)/6.0E+1+((x0*x0)*x2*y0)/4.0E+1+(x0*(x1*x1)*y2)/1.2E+2+(x0*(x2*x2)*y1)/1.2E+2+(x1*(x2*x2)*y0)/1.2E+2+((x0*x0)*x1*y2)/1.2E+2+((x0*x0)*x2*y1)/1.2E+2+((x1*x1)*x2*y0)/1.2E+2+(x0*(x2*x2)*y2)/4.0E+1+(x1*(x2*x2)*y1)/6.0E+1+((x0*x0)*x2*y2)/6.0E+1+((x1*x1)*x2*y1)/4.0E+1+(x1*(x2*x2)*y2)/4.0E+1+((x1*x1)*x2*y2)/6.0E+1+(x0*x1*x2*y0)/6.0E+1+(x0*x1*x2*y1)/6.0E+1+(x0*x1*x2*y2)/6.0E+1;
        
        double YYY = ((y0*y0*y0)*z0)/3.0E+1+((y0*y0*y0)*z1)/1.2E+2+((y1*y1*y1)*z0)/1.2E+2+((y0*y0*y0)*z2)/1.2E+2+((y1*y1*y1)*z1)/3.0E+1+((y2*y2*y2)*z0)/1.2E+2+((y1*y1*y1)*z2)/1.2E+2+((y2*y2*y2)*z1)/1.2E+2+((y2*y2*y2)*z2)/3.0E+1+(y0*(y1*y1)*z0)/6.0E+1+((y0*y0)*y1*z0)/4.0E+1+(y0*(y1*y1)*z1)/4.0E+1+(y0*(y2*y2)*z0)/6.0E+1+((y0*y0)*y1*z1)/6.0E+1+((y0*y0)*y2*z0)/4.0E+1+(y0*(y1*y1)*z2)/1.2E+2+(y0*(y2*y2)*z1)/1.2E+2+(y1*(y2*y2)*z0)/1.2E+2+((y0*y0)*y1*z2)/1.2E+2+((y0*y0)*y2*z1)/1.2E+2+((y1*y1)*y2*z0)/1.2E+2+(y0*(y2*y2)*z2)/4.0E+1+(y1*(y2*y2)*z1)/6.0E+1+((y0*y0)*y2*z2)/6.0E+1+((y1*y1)*y2*z1)/4.0E+1+(y1*(y2*y2)*z2)/4.0E+1+((y1*y1)*y2*z2)/6.0E+1+(y0*y1*y2*z0)/6.0E+1+(y0*y1*y2*z1)/6.0E+1+(y0*y1*y2*z2)/6.0E+1;
        
        double ZZZ = (x0*(z0*z0*z0))/3.0E+1+(x0*(z1*z1*z1))/1.2E+2+(x1*(z0*z0*z0))/1.2E+2+(x0*(z2*z2*z2))/1.2E+2+(x1*(z1*z1*z1))/3.0E+1+(x2*(z0*z0*z0))/1.2E+2+(x1*(z2*z2*z2))/1.2E+2+(x2*(z1*z1*z1))/1.2E+2+(x2*(z2*z2*z2))/3.0E+1+(x0*z0*(z1*z1))/6.0E+1+(x0*(z0*z0)*z1)/4.0E+1+(x0*z0*(z2*z2))/6.0E+1+(x0*(z0*z0)*z2)/4.0E+1+(x1*z0*(z1*z1))/4.0E+1+(x1*(z0*z0)*z1)/6.0E+1+(x0*z1*(z2*z2))/1.2E+2+(x0*(z1*z1)*z2)/1.2E+2+(x1*z0*(z2*z2))/1.2E+2+(x1*(z0*z0)*z2)/1.2E+2+(x2*z0*(z1*z1))/1.2E+2+(x2*(z0*z0)*z1)/1.2E+2+(x1*z1*(z2*z2))/6.0E+1+(x1*(z1*z1)*z2)/4.0E+1+(x2*z0*(z2*z2))/4.0E+1+(x2*(z0*z0)*z2)/6.0E+1+(x2*z1*(z2*z2))/4.0E+1+(x2*(z1*z1)*z2)/6.0E+1+(x0*z0*z1*z2)/6.0E+1+(x1*z0*z1*z2)/6.0E+1+(x2*z0*z1*z2)/6.0E+1;

        double XXY = ((x0*x0)*y0)/2.0E+1+((x0*x0)*y1)/6.0E+1+((x1*x1)*y0)/6.0E+1+((x0*x0)*y2)/6.0E+1+((x1*x1)*y1)/2.0E+1+((x2*x2)*y0)/6.0E+1+((x1*x1)*y2)/6.0E+1+((x2*x2)*y1)/6.0E+1+((x2*x2)*y2)/2.0E+1+(x0*x1*y0)/3.0E+1+(x0*x1*y1)/3.0E+1+(x0*x2*y0)/3.0E+1+(x0*x1*y2)/6.0E+1+(x0*x2*y1)/6.0E+1+(x1*x2*y0)/6.0E+1+(x0*x2*y2)/3.0E+1+(x1*x2*y1)/3.0E+1+(x1*x2*y2)/3.0E+1;
        
        double YYZ = ((y0*y0)*z0)/2.0E+1+((y0*y0)*z1)/6.0E+1+((y1*y1)*z0)/6.0E+1+((y0*y0)*z2)/6.0E+1+((y1*y1)*z1)/2.0E+1+((y2*y2)*z0)/6.0E+1+((y1*y1)*z2)/6.0E+1+((y2*y2)*z1)/6.0E+1+((y2*y2)*z2)/2.0E+1+(y0*y1*z0)/3.0E+1+(y0*y1*z1)/3.0E+1+(y0*y2*z0)/3.0E+1+(y0*y1*z2)/6.0E+1+(y0*y2*z1)/6.0E+1+(y1*y2*z0)/6.0E+1+(y0*y2*z2)/3.0E+1+(y1*y2*z1)/3.0E+1+(y1*y2*z2)/3.0E+1;
        
        double ZZX = (x0*(z0*z0))/2.0E+1+(x0*(z1*z1))/6.0E+1+(x1*(z0*z0))/6.0E+1+(x0*(z2*z2))/6.0E+1+(x1*(z1*z1))/2.0E+1+(x2*(z0*z0))/6.0E+1+(x1*(z2*z2))/6.0E+1+(x2*(z1*z1))/6.0E+1+(x2*(z2*z2))/2.0E+1+(x0*z0*z1)/3.0E+1+(x0*z0*z2)/3.0E+1+(x1*z0*z1)/3.0E+1+(x0*z1*z2)/6.0E+1+(x1*z0*z2)/6.0E+1+(x2*z0*z1)/6.0E+1+(x1*z1*z2)/3.0E+1+(x2*z0*z2)/3.0E+1+(x2*z1*z2)/3.0E+1;


        // 1 in volume
        // X in surface
        mass += X * density;
        
        // mass * com
        center += X * density * com;
        
        // y^2 + z^2 in volume
        // 1/3 YYY ZZZ in surface
        I(0, 0) += 1/3 * (YYY + ZZZ);
        // -xy
        // - 1/2 XXY in surface
        I(0, 1) += -0.5 * XXY;
        // -xz
        // - 1/2 XZZ in surface
        I(0, 2) += -0.5 * ZZX;
        // -xy
        // - 1/2 XXY in surface
        I(1, 0) += -0.5 * XXY;
        // x^2 + z^2
        // 1/3 XXX ZZZ in surface
        I(1, 1) += 1/3 * (XXX + ZZZ);
        // -yz
        // - 1/2 YYZ in surface
        I(1, 2) += -0.5 * YYZ;
        // -xz
        // - 1/2 XZZ in surface
        I(2, 0) += -0.5 * ZZX;
        // -yz
        // - 1/2 YYZ in surface
        I(2, 1) += -0.5 * YYZ;
        // x^2 + y^2
        // 1/3 XXX YYY in surface
        I(2, 2) += 1/3 * (XXX + YYY);
        
    }
    // sum(mass_i * com_i) / total_mass
    center = center / mass;
    */

    
    

    std::vector<double> intg (10, 0);
    //std::vector<double> mult{1/6.0, 1/24.0, 1/24.0, 1/24.0, 1/60.0, 1/60.0, 1/60.0, 1/120.0, 1/120.0, 1/120.0};
    for (unsigned idx = 0; idx < n_triangle; idx++)
    {

        unsigned v1 = F(idx, 0), v2 = F(idx, 1), v3 = F(idx, 2);
        double x0 = V(v1, 0), x1 = V(v2, 0), x2 = V(v3, 0),
               y0 = V(v1, 1), y1 = V(v2, 1), y2 = V(v3, 1),
               z0 = V(v1, 2), z1 = V(v2, 2), z2 = V(v3, 2);

        double a1 = x1 - x0, b1 = y1 - y0, c1 = z1 - z0, a2 = x2 - x0, b2 = y2 - y0, c2 = z2 -z0, d0 = b1*c2 - b2*c1, d1 = a2*c1 - a1*c2,d2 = a1 * b2 - a2 *b1;
        
        double f1x, f2x, f3x, f1y, f2y, f3y, f1z, f2z, f3z, g0x, g0y, g0z, g1x, g1y, g1z, g2x, g2y, g2z;
        integrate(x0, x1, x2, f1x, f2x, f3x, g0x, g1x, g2x);
        integrate(y0, y1, y2, f1y, f2y, f3y, g0y, g1y, g2y);
        integrate(z0, z1, z2, f1z, f2z, f3z, g0z, g1z, g2z);
        //std::cout<<f1x<<" "<<f2x<<" "<<f3x<<" "<<g0x<<" "<<g1x<<" "<<g2x<<std::endl;
        //exit(0);
        intg[0] +=d0 * f1x;
        intg[1] +=d0 * f2x;
        intg[2] +=d1 * f2y;
        intg[3] +=d2 * f2z;
        intg[4] +=d0 * f3x;
        intg[5] +=d1 * f3y;
        intg[6] +=d2 * f3z;
        intg[7] +=d0 * (y0 * g0x + y1 * g1x + y2 * g2x);
        intg[8] +=d1 * (z0 * g0y + z1 * g1y + z2 * g2y);
        intg[9] +=d2 * (x0 * g0z + x1 * g1z + x2 * g2z);
    }
    //for (unsigned i = 0; i < intg.size(); i++ ) std::cout<<intg[i]<<std::endl;
    //for (unsigned i = 0; i < mult.size(); i++ ) std::cout<<mult[i]<<std::endl;
    //exit(0);
    intg[0] *= 1/6.0;
    intg[1] *=1/24.0;
    intg[2] *=1/24.0;
    intg[3] *=1/24.0;
    intg[4] *=1/60.0;
    intg[5] *=1/60.0;
    intg[6] *=1/60.0;
    intg[7] *=1/120.0;
    intg[8] *=1/120.0;
    intg[9] *=1/120.0;
    
    mass = intg[0] * density;
    double c_x = intg[1]/intg[0], c_y = intg[2]/intg[0], c_z = intg[3]/intg[0];
    center(0) = c_x;
    center(1) = c_y;
    center(2) = c_z;
    
    double XX = intg[5] * density + intg[6] * density - mass * (c_y * c_y + c_z * c_z),
           YY = intg[4] * density + intg[6] * density - mass * (c_z * c_z + c_x * c_x),
           ZZ = intg[4] * density + intg[5] * density - mass * (c_x * c_x + c_y * c_y),
           XY = -(intg[7] * density - mass * c_x * c_y),
           YZ = -(intg[8] * density - mass * c_y * c_z),
           XZ = -(intg[9] * density - mass * c_z * c_x);
    
    for (unsigned i = 0; i < intg.size(); i++ ) std::cout<<intg[i]<<std::endl;
    
    I(0, 0) = YY + ZZ;
    I(0, 1) = XY;
    I(0, 2) = XZ;
    I(1, 0) = XY;
    I(1, 1) = XX + ZZ;
    I(1, 2) = YZ;
    I(2, 0) = XZ;
    I(2, 1) = YZ;
    I(2, 2) = XX + YY;
    
    
    //std::cout<<"COM: "<<std::endl;
    //std::cout<<center<<std::endl;
    //std::cout<<"mass: "<<std::endl;
    //std::cout<<mass<<std::endl;
    //std::cout<<"I: "<<std::endl;
    //std::cout<<I<<std::endl;
    //exit(0);
}
