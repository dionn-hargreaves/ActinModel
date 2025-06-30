#ifndef ELLIPSECONFIDENCE_H
#define ELLIPSECONFIDENCE_H

#include <Eigen/Dense>


std::array<double,3>  ellipsoidfit(const std::vector< std::array<double,3> > &data,
                                   int dimensions, Eigen::Matrix3d &rotationMatrix,
                                   std::array<double,3> &centrePoint);


#endif
