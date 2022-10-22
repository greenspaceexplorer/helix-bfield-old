#ifndef HEAT_MAGNET_H
#define HEAT_MAGNET_H

#include "common.h"
#include <cmath>
#include <Eigen/Geometry> 

using namespace Eigen;
Vector3d HeatB(Vector3d &position);
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
//
#endif
