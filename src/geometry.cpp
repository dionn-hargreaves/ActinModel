/*
 *  geometry.cpp
 *
 *  C++ file containing the definition of the geometry functions.

 *  The declaration of this class (including the definitions of the trivial
 *  member functions) is found in the associated header file "Actin.h"
 *
 *  James Bradford
 *  University of Sheffield
 *  May 2017
 */

#include "configHeader.h"
#include "geometry.h"

double distanceBetPoints2DSQR(std::array<double,2> point1, std::array<double,2> point2)
{
    return ( (point1[0] - point2[0]) * (point1[0] - point2[0]) +
             (point1[1] - point2[1]) * (point1[1] - point2[1]) );
}
