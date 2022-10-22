#include "heat_magnet.h"

Vector3d HeatB(Vector3d &position){
    double dconvert = 1.e2;
    double bconvert = 1.e-4;
    double x,y,z,rho;
    double bx,by,bz,br,b;

    x = position(0)*dconvert;
    y = position(1)*dconvert;
    z = position(2)*dconvert;

    heat_mag_(&x,&y,&z,&br,&bz,&b);
    br *= bconvert;
    bz *= bconvert;
    b *= bconvert;

    rho = sqrt(x*x+y*y);
    if(rho == 0.){
        bx = br;
        by = br;
    }
    else{
        bx = br*x/rho;
        by = br*y/rho;
    }

    Vector3d bout(bx,by,bz);
    return bout;
}
//Vector3d HeatSubcoil(
//        double rout,
//        double rin,
//        double width,
//        double current,
//        double turns,
//        double rsubdiv,
//        double zsubdiv,
//        Vector3d &position
//        );
//Vector3d HeatIdealCoil(
//        double current,
//        double radius,
//        double turns,
//        Vector3d &position
//        );
