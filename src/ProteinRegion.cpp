/*
 *  ProteinRegion.cpp
 *
 *  C++ file containing the definition of the ProteinRegion class.
 *  This contains member variable initialisations and non-trivial member
 *  functions.
 *  The declaration of this class (including the definitions of the trivial
 *  member functions) is found in the associated header file "ProteinRegion.h"
 *
 *  James Bradford
 *  University of Sheffield
 *  Apr 2019
 */

// Constructor for non axis-aligned rectangle

// More general than ProteinRegion.cpp and this will eventually replace it after
// some refinement

#include "configHeader.h"
#include "globals.h"
#include "Actin.h"
#include "ProteinRegion.h"
#include "RNG.h"
#include "geometry.h"

extern RNG rng;

// Rectangular region
ProteinRegion::ProteinRegion(std::array<double,2> xUVec, double width,
                             double height, std::array<double,2> pointBL, double arpConc,
                             double capCoeff, double nucCoeff,
                             double sevCoeff, bool coupledToMem, int coupledMemSub)
    : m_uVec { xUVec },
    m_pointBL { pointBL },
    m_width { width },
    m_height { height },
    m_arpConc { arpConc },
    m_capCoeff { capCoeff },
    m_nucCoeff { nucCoeff },
    m_sevCoeff { sevCoeff },
    m_coupledBool { coupledToMem },
    m_coupledMemSub { coupledMemSub },
    m_isCirc { false },
    m_isRing { false }
{
    // Thoughts: Have below two as generating in the frame of the rectangle
    // so 0 to m_width and 0 to m_height. Then convert to global frame
    m_rectangleWidth = std::uniform_real_distribution<double> (0, m_width);
    m_rectangleHeight = std::uniform_real_distribution<double> (0, m_height);
    m_area = m_width * m_height;
}

// Circular region
ProteinRegion::ProteinRegion(double radius, std::array<double,2> pointCentre, double arpConc,
                             double capCoeff, double nucCoeff,
                             double sevCoeff, bool coupledToMem, int coupledMemSub)
    : m_arpConc { arpConc },
    m_capCoeff { capCoeff },
    m_nucCoeff { nucCoeff },
    m_sevCoeff { sevCoeff },
    m_coupledBool { coupledToMem },
    m_coupledMemSub { coupledMemSub },
    m_isCirc { true },
    m_isRing { false },
    m_pointC { pointCentre },
    m_radius { radius }
{
    // Thoughts: Have below two as generating in the frame of the rectangle
    // so 0 to m_width and 0 to m_height. Then convert to global frame
    m_RDistro = std::uniform_real_distribution<double> (0, m_radius);
    m_ThetaDistro = std::uniform_real_distribution<double> (0, 2*M_PI);
    m_area = M_PI*m_radius*m_radius;
}

// Ring region (nuc only)
ProteinRegion::ProteinRegion(std::array<double,2> pointCentre, double radius,
                             double innerRadius, double nucCoeff,
                             bool coupledToMem, int coupledMemSub)
    : m_arpConc { 0 },
    m_capCoeff { 0 },
    m_nucCoeff { nucCoeff },
    m_sevCoeff { 0 },
    m_coupledBool { coupledToMem },
    m_coupledMemSub { coupledMemSub },
    m_isCirc { false },
    m_isRing { true },
    m_pointC { pointCentre },
    m_radius { radius },
    m_innerRadius { innerRadius }
{
    assert(m_radius > m_innerRadius);
    // Thoughts: Have below two as generating in the frame of the rectangle
    // so 0 to m_width and 0 to m_height. Then convert to global frame
    m_RDistro = std::uniform_real_distribution<double> (m_innerRadius, m_radius);
    m_ThetaDistro = std::uniform_real_distribution<double> (0, 2*M_PI);
    m_area = M_PI*(m_radius*m_radius - m_innerRadius*m_innerRadius);
}

bool ProteinRegion::checkPointwithin(std::array<double,2> Point) const
{
    if (m_isCirc)
    {
        return checkPointwithinCirc(Point);
    }
    else if (m_isRing)
    {
        return checkPointwithinRing(Point);
    }
    else
    {
        return checkPointwithinRect(Point);
    }
}


bool ProteinRegion::checkPointwithinRect(std::array<double,2> Point) const
{
    // Simple function that returns true if the given 2D point lies inside the
    // rectangle

    // Step 1 - move the point so that the bottom left of the rectangle lies on the origin

    Point[0] -= m_pointBL[0];
    Point[1] -= m_pointBL[1];

    // Step 2 - rotate the point from the origin, clockwise by the angle the rectangle
    // makes with the x axis

    std::array<double,2> PointPrime;
    PointPrime[0] = m_uVec[0]*Point[0] + m_uVec[1]*Point[1];
    PointPrime[1] = -m_uVec[1]*Point[0] + m_uVec[0]*Point[1];

    // Step 3 - Compare with the axis aligned rectangle that has its bottom left
    // point at the origin
    if ( (PointPrime[0] < m_width) && (PointPrime[0] > 0) && (PointPrime[1] < m_height) && (PointPrime[1] > 0) )
    {
        // Point lies inside the Rectangle
        return true;
    }


    return false;
}

bool ProteinRegion::checkPointwithinCirc(std::array<double,2> Point) const
{
    // Simple function that returns true if the given 2D point lies inside the
    // circle


    if (sqrt(distanceBetPoints2DSQR(m_pointC, Point)) < m_radius)
    {
        // in the circle
        return true;
    }
    else
    {
        return false;
    }

}

bool ProteinRegion::checkPointwithinRing(std::array<double,2> Point) const
{
    // Simple function that returns true if the given 2D point lies inside the
    // circle

    double dist = sqrt(distanceBetPoints2DSQR(m_pointC, Point));
    if (dist < m_radius && dist > m_innerRadius)
    {
        // in the ring
        return true;
    }
    else
    {
        return false;
    }

}

void ProteinRegion::splitRegion(std::vector<ProteinRegion> &regionVec,
                                std::vector<MembraneWall> &memWalls)
{
    // Function that is called when one region needs to be split into two equal
    // regions
    assert(!m_isCirc);
    assert(!m_isRing);
    // 1. Decide which axis to split by, x or y
    if (m_width > m_height)
    {
        // split along the x axis
        // adjust new bottom left point
        std::array<double,2> newBLPoint = { m_pointBL[0] + m_uVec[0]*(m_width/2), m_pointBL[1] + m_uVec[1]*(m_width/2) };

        // 2. create new region
        ProteinRegion newRegion(m_uVec, m_width/2, m_height, newBLPoint,
                                m_arpConc, m_capCoeff, m_nucCoeff, m_sevCoeff,
                                m_coupledBool, m_coupledMemSub);
        if (memWalls.size() == 1 && m_coupledBool)
        {
          memWalls[0].appendToNucRegions(regionVec.size());
        }
        // 3. Adjust pre-existing region
        m_width /= 2;
        m_rectangleWidth = std::uniform_real_distribution<double> (0, m_width);
        m_area = m_width * m_height;
        regionVec.push_back(newRegion);
    }
    else
    {
        // split along the y axis

        // Need unit vector representing y direction
        std::array<double,2> uVECY;
        uVECY[0] = -m_uVec[1];
        uVECY[1] = m_uVec[0];
        std::array<double,2> newBLPoint = { m_pointBL[0] + uVECY[0]*(m_height/2), m_pointBL[1] + uVECY[1]*(m_height/2) };

        ProteinRegion newRegion(m_uVec, m_width, m_height/2, newBLPoint,
                                m_arpConc, m_capCoeff, m_nucCoeff, m_sevCoeff,
                                m_coupledBool, m_coupledMemSub);

        if (memWalls.size() == 1 && m_coupledBool)
        {
          memWalls[0].appendToNucRegions(regionVec.size());
        }

        m_height /= 2;
        m_rectangleHeight = std::uniform_real_distribution<double> (0, m_height);

        m_area = m_width * m_height;
        regionVec.push_back(newRegion);

    }
}

int ProteinRegion::s_chooseRegion(std::vector<ProteinRegion> regions, int nRegions)
{
    // Given multiple nucleation regions, this function chooses where to
    // nucleate a new filament, based on area
    double rnd = rng.m_probDist(rng.m_mersenne); // generate number between 0 and 1
    for (int i = 0; i < nRegions; ++i)
    {
        if (rnd < regions[i].getNormalisedArea())
            return i;

        rnd -= regions[i].getNormalisedArea();
    }
    std::cout << "Error in ProteinRegion.cpp line 138" << std::endl;
    exit(1);
    return 0;
}

void ProteinRegion::moveRegionTo(std::array<double,2> newPoint,
                                 std::array<double,2> newUVec)
{
    if (m_isCirc || m_isRing)
    {
        moveRegionToCirc(newPoint);
    }
    else
    {
        moveRegionToRect(newPoint, newUVec);
    }
}

void ProteinRegion::moveRegionToRect(std::array<double,2> newBLPoint,
                                    std::array<double,2> newUVec)
{
    m_uVec = newUVec;
    m_pointBL = newBLPoint;
}

void ProteinRegion::moveRegionToCirc(std::array<double,2> newCPoint)
{
    m_pointC = newCPoint;
}

void ProteinRegion::moveRegion(double dx, double dy)
{
    if (m_isCirc || m_isRing)
    {
        m_pointC[0] += dx;
        m_pointC[1] += dy;
    }
    else
    {
        m_pointBL[0] += dx;
        m_pointBL[1] += dy;
    }
}

void ProteinRegion::stickToMem(Membrane &membrane)
{
    if (!m_isCirc && !m_isRing)
    {
        // Change m_pointBL to be stuck to the membrane equal to point coupledMemSub
        m_pointBL[0] = membrane.getPoints()[m_coupledMemSub][0];
        m_pointBL[1] = membrane.getPoints()[m_coupledMemSub][1];

        m_uVec[0] = membrane.getUnitVecs()[m_coupledMemSub][0];
        m_uVec[1] = membrane.getUnitVecs()[m_coupledMemSub][1];
    }

}

std::array<double,2> ProteinRegion::getCentre() const
{

    std::array<double,2> centre;
    // Need to build vector to get from BL to centre: cVec
    // Need unit vector representing y direction

    if (m_isCirc || m_isRing)
    {
        centre = m_pointC;
    }
    else
    {
        std::array<double,2> cVec;
        cVec[0] = m_uVec[0]*(m_width/2) - m_uVec[1]*(m_height/2);
        cVec[1] = m_uVec[1]*(m_height/2) + m_uVec[0]*(m_height/2);


        centre[0] = m_pointBL[0] + cVec[0];
        centre[1] = m_pointBL[1] + cVec[1];
    }

    return centre;
}

std::array<double,8> ProteinRegion::getCorners() const
{
    // Function that returns 4 corners of the rectangle indexed as...
    // 0, 1 : bottom left
    // 2, 3 : top left
    // 4, 5 : top right
    // 6, 7 : bottom right

    std::array<double,8> corners;
    if (!m_isCirc && !m_isRing)
    {
        // Bottom left
        corners[0] = m_pointBL[0];
        corners[1] = m_pointBL[1];

        // Top left
        corners[2] = m_pointBL[0] - m_uVec[1]*m_height;
        corners[3] = m_pointBL[1] + m_uVec[0]*m_height;

        // Top right
        corners[4] = corners[2] + m_uVec[0]*m_width;
        corners[5] = corners[3] + m_uVec[1]*m_width;

        // Bottom right
        corners[6] = m_pointBL[0] + m_uVec[0]*m_width;
        corners[7] = m_pointBL[1] + m_uVec[1]*m_width;
    }
    else
    {
        std::cout << "Circle no corners!" << std::endl;
        exit(1);
    }

    return corners;

}

void ProteinRegion::adjustWidth(double newWidth)
{
    m_width = newWidth;
    m_rectangleWidth = std::uniform_real_distribution<double> (0, m_width);
    m_area = m_width * m_height;
}
