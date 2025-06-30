/*
 *  ExcZone.cpp
 *
 *  C++ file containing the definition of the Exclusion zone class.
 *  This contains member variable initialisations and non-trivial member
 *  functions.
 *  The declaration of this class (including the definitions of the trivial
 *  member functions) is found in the associated header file "ExcZone.h"
 *
 *  James Bradford
 *  University of Sheffield
 *  March 2021
 */

#include "configHeader.h"
#include "geometry.h"
#include "globals.h"
#include "ExcZone.h"
#include "GactinGrid.h"
#include "Actin.h"
#include "Membrane.h"
#include "Cortex.h"

extern RNG rng; // this means RNG rng is defined elsewhere (global variable) - random number defined in main.cpp
typedef gte::DCPQuery<double, gte::Vector<2,double>, gte::Segment<2, double>> DistQuery; // distance between point and segment
typedef gte::DCPQuery<double, gte::Segment<2, double>, gte::Segment<2, double>> SegDistQuery; // distance between two segments

// Constructor for a circular exclusion zone
ExcZone::ExcZone(double x, double y, double radius, int ID, bool motion, bool directedmotion)
    : m_id { ID },
    m_x { x },
    m_y { y },
    m_radius { radius },
    m_isCircle { 1 },
    m_isRectangle { 0 },
    m_motion { motion },
    m_directedmotion { directedmotion },
    m_totalMoveDown { 0 }
{
}

// Constructor for a rectangular exclusion zone
ExcZone::ExcZone(double x_min, double y_min, double x_max, double y_max, int ID, bool dummy1, bool dummy2)
    : m_id { ID },
    m_x_min { x_min },
    m_y_min { y_min },
    m_x_max { x_max },
    m_y_max { y_max },
    m_isCircle { 0 },
    m_isRectangle { 1 },
    m_totalMoveDown { 0 }
{
}

// Constructor for a triangular exclusion zone
ExcZone::ExcZone(double p1X, double p1Y, double len1, double len2, int ID, bool directedmotion, bool dummy1, bool dummy2)
    : m_id { ID },
    m_len1 { len1 },
    m_len2 { len2 },
    m_isCircle { 0 },
    m_isRectangle { 0 },
    m_directedmotion { directedmotion },
    m_totalMoveDown { 0 }
{
    m_p1[0] = p1X;
    m_p1[1] = p1Y;

    m_p2[0] = m_p1[0] - m_len2/2;
    m_p2[1] = m_p1[1] + m_len1;

    m_p3[0] = m_p2[0] + m_len2;
    m_p3[1] = m_p2[1];

}


bool ExcZone::check_ex_vol(const Actin &actin) const
{
    // Called upon nucleation or branching of filament
    if (m_isCircle)
    {
        //return false;
        return circ_check_ex_vol(actin);
    }
    else if (m_isRectangle)
    {
        return rect_check_ex_vol(actin);
    }
    else
    {
      // triangle
      return false;

    }
    // Error
    std::cout << "Error, ExcZone.cpp line 58" << std::endl;
    exit(1);
    return 0;
}


bool ExcZone::s_checkExVol(int subid, const Actin &actin,
                          const std::vector<ExcZone> excZones)
{
    // Called upon bending or diffusion of actin
    for (unsigned int i = 0; i < excZones.size(); ++i)
    {
        if (excZones[i].getCircBool())
        {
            if (excZones[i].circ_check_ex_vol(subid, actin))
            {
                return true;
            }
        }
        else
        {
            if (excZones[i].rect_check_ex_vol(subid, actin))
            {
                return true;
            }
        }

    }

    return false;
}

bool ExcZone::circ_check_ex_vol(int subid, const Actin &actin) const
{
    // Does a distance check between the actin filament and the centre point
    // of the circle, if this is less than the combined radii of both filament
    // and circle then the filament is inside the circle.

    // This should work for 3D and 2D
    // Returns true if violated

    // Changed to work with bent filaments
    // Define centre point of circle
    gte::Vector2<double> point = {m_x, m_y};

    gte::Segment<2,double> subFilament;

    DistQuery min_d_GTE;

    subFilament.p[0][0] = actin.getPoints()[subid][0];
    subFilament.p[1][0] = actin.getPoints()[subid+1][0];
    subFilament.p[0][1] = actin.getPoints()[subid][1];
    subFilament.p[1][1] = actin.getPoints()[subid+1][1];

    auto result = min_d_GTE(point, subFilament);
    if (result.distance < (actin.getStericRadius()+m_radius))
        return true;

    return false;
}

bool ExcZone::circ_check_ex_vol(const Actin &actin) const
{
    // Does a distance check between the actin filament and the centre point
    // of the circle, if this is less than the combined radii of both filament
    // and circle then the filament is inside the circle.

    // This should work for 3D and 2D
    // Returns true if violated

    // Changed to work with bent filaments
    // Define centre point of circle
    gte::Vector2<double> point = {m_x, m_y};

    gte::Segment<2,double> subFilament;

    DistQuery min_d_GTE;

    for (int i = 0; i < actin.getNumSubs(); ++i)
    {
        subFilament.p[0][0] = actin.getPoints()[i][0];
        subFilament.p[1][0] = actin.getPoints()[i+1][0];
        subFilament.p[0][1] = actin.getPoints()[i][1];
        subFilament.p[1][1] = actin.getPoints()[i+1][1];

        auto result = min_d_GTE(point, subFilament);
        if (result.distance < (actin.getStericRadius()+m_radius))
        {
            return true;
        }
    }
    return false;
}


bool ExcZone::rect_check_ex_vol(int subid, const Actin &actin) const
{
    // Does a distance check between the actin filament all sides of the
    // rectangle/square

    // Modified for bending

    // This should work for 3D and 2D
    // Returns true if violated

    // Define sides of rectangle

    gte::Segment<2,double> side1;
    side1.p[0][0] = m_x_min;
    side1.p[1][0] = m_x_min;
    side1.p[0][1] = m_y_min;
    side1.p[1][1] = m_y_max;

    gte::Segment<2,double> side2;
    side2.p[0][0] = m_x_min;
    side2.p[1][0] = m_x_max;
    side2.p[0][1] = m_y_max;
    side2.p[1][1] = m_y_max;

    gte::Segment<2,double> side3;
    side3.p[0][0] = m_x_max;
    side3.p[1][0] = m_x_max;
    side3.p[0][1] = m_y_max;
    side3.p[1][1] = m_y_min;

    gte::Segment<2,double> side4;
    side4.p[0][0] = m_x_max;
    side4.p[1][0] = m_x_min;
    side4.p[0][1] = m_y_min;
    side4.p[1][1] = m_y_min;

    SegDistQuery min_d_GTE;
    gte::Segment<2,double> subFilament;



    // Check the subfila to see if intersect with any sides

    subFilament.p[0][0] = actin.getPoints()[subid][0];
    subFilament.p[1][0] = actin.getPoints()[subid+1][0];
    subFilament.p[0][1] = actin.getPoints()[subid][1];
    subFilament.p[1][1] = actin.getPoints()[subid+1][1];

    auto result = min_d_GTE(side1, subFilament);
    if (result.distance < (actin.getStericRadius()))
        return true;

    result = min_d_GTE(side2, subFilament);
    if (result.distance < (actin.getStericRadius()))
        return true;

    result = min_d_GTE(side3, subFilament);
    if (result.distance < (actin.getStericRadius()))
        return true;

    result = min_d_GTE(side4, subFilament);
    if (result.distance < (actin.getStericRadius()))
        return true;


    // the filament is either completely inside the rectangle or outside, so
    // we need to check if a point on the filament is inside

    std::array<double,2> point { actin.getPoints()[subid][0], actin.getPoints()[subid][1] };

    return rect_checkPointwithin(point);

}

bool ExcZone::rect_check_ex_vol(const Actin &actin) const
{
    // Does a distance check between the actin filament all sides of the
    // rectangle/square

    // Modified for bending

    // This should work for 3D and 2D
    // Returns true if violated

    // Define sides of rectangle

    gte::Segment<2,double> side1;
    side1.p[0][0] = m_x_min;
    side1.p[1][0] = m_x_min;
    side1.p[0][1] = m_y_min;
    side1.p[1][1] = m_y_max;

    gte::Segment<2,double> side2;
    side2.p[0][0] = m_x_min;
    side2.p[1][0] = m_x_max;
    side2.p[0][1] = m_y_max;
    side2.p[1][1] = m_y_max;

    gte::Segment<2,double> side3;
    side3.p[0][0] = m_x_max;
    side3.p[1][0] = m_x_max;
    side3.p[0][1] = m_y_max;
    side3.p[1][1] = m_y_min;

    gte::Segment<2,double> side4;
    side4.p[0][0] = m_x_max;
    side4.p[1][0] = m_x_min;
    side4.p[0][1] = m_y_min;
    side4.p[1][1] = m_y_min;

    SegDistQuery min_d_GTE;
    gte::Segment<2,double> subFilament;



    // Loop over all subunits to check if filament intersects the sides of
    // the rectangle

    for (int i = 0; i < actin.getNumSubs(); ++i)
    {
        subFilament.p[0][0] = actin.getPoints()[i][0];
        subFilament.p[1][0] = actin.getPoints()[i+1][0];
        subFilament.p[0][1] = actin.getPoints()[i][1];
        subFilament.p[1][1] = actin.getPoints()[i+1][1];

        auto result = min_d_GTE(side1, subFilament);
        if (result.distance < (actin.getStericRadius()))
            return true;

        result = min_d_GTE(side2, subFilament);
        if (result.distance < (actin.getStericRadius()))
            return true;

        result = min_d_GTE(side3, subFilament);
        if (result.distance < (actin.getStericRadius()))
            return true;

        result = min_d_GTE(side4, subFilament);
        if (result.distance < (actin.getStericRadius()))
            return true;


    }

    // the filament is either completely inside the rectangle or outside, so
    // we need to check if a point on the filament is inside

    std::array<double,2> barbedEnd { actin.getBarbedEnd()[0], actin.getBarbedEnd()[1] };

    return rect_checkPointwithin(barbedEnd);

}

bool ExcZone::rect_checkPointwithin(std::array<double,2> Point) const
{
    // Simple function that returns true if the given 2D point lies inside the
    // rectangle
    if ( (Point[0] < m_x_max) && (Point[0] > m_x_min) && (Point[1] < m_y_max) && (Point[1] > m_y_min) )
    {
        // Point lies inside the axis aligned Rectangle
        return true;
    }

    return false;
}

bool ExcZone::s_check_polymerisation_barbed(const Actin &actin, const std::vector<ExcZone> &excZones)
{
    if (excZones.empty())
        return false;

    for (ExcZone excZone : excZones)
    {
        if (excZone.getCircBool())
        {
            if (excZone.circ_check_polymerisation_barbed(actin))
                return true;
        }
        else if (excZone.getRectBool())
        {
            if (excZone.rect_check_polymerisation_barbed(actin))
                return true;
        }
    }
    return false;
}

bool ExcZone::circ_check_polymerisation_barbed(const Actin &actin)
{
    double barbx = actin.getBarbedEnd()[0];
    double barby = actin.getBarbedEnd()[1];

    double d = sqrt( ( (barbx-m_x)*(barbx-m_x) + (barby-m_y)*(barby-m_y) ) );

    return (d < (actin.getStericRadius() + m_radius));
}

bool ExcZone::rect_check_polymerisation_barbed(const Actin &actin)
{

    gte::Segment<2,double> endSub;
    endSub.p[0][0] = actin.getBarbedEnd()[0];
    endSub.p[1][0] = actin.getPreBarbedEnd()[0];
    endSub.p[0][1] = actin.getBarbedEnd()[1];
    endSub.p[1][1] = actin.getPreBarbedEnd()[1];

    gte::Segment<2,double> side1;
    side1.p[0][0] = m_x_min;
    side1.p[1][0] = m_x_min;
    side1.p[0][1] = m_y_min;
    side1.p[1][1] = m_y_max;

    SegDistQuery min_d_GTE;
    auto result = min_d_GTE(endSub, side1);

    if (result.distance < actin.getStericRadius())
        return true;

    gte::Segment<2,double> side2;
    side2.p[0][0] = m_x_min;
    side2.p[1][0] = m_x_max;
    side2.p[0][1] = m_y_max;
    side2.p[1][1] = m_y_max;

    result = min_d_GTE(endSub, side2);
    if (result.distance < actin.getStericRadius())
        return true;

    gte::Segment<2,double> side3;
    side3.p[0][0] = m_x_max;
    side3.p[1][0] = m_x_max;
    side3.p[0][1] = m_y_max;
    side3.p[1][1] = m_y_min;

    result = min_d_GTE(endSub, side3);
    if (result.distance < actin.getStericRadius())
        return true;

    gte::Segment<2,double> side4;
    side4.p[0][0] = m_x_max;
    side4.p[1][0] = m_x_min;
    side4.p[0][1] = m_y_min;
    side4.p[1][1] = m_y_min;

    result = min_d_GTE(endSub, side4);
    if (result.distance < actin.getStericRadius())
        return true;

    std::array<double,2> barbedEnd { actin.getBarbedEnd()[0], actin.getBarbedEnd()[1] };

    return rect_checkPointwithin(barbedEnd);
}

bool ExcZone::s_check_polymerisation_pointed(const Actin &actin, const std::vector<ExcZone> &excZones)
{
    if (excZones.empty())
        return false;

    for (ExcZone excZone : excZones)
    {
        if (excZone.getCircBool())
        {
            if (excZone.circ_check_polymerisation_pointed(actin))
                return true;
        }
        else if (excZone.getRectBool())
        {
            if (excZone.rect_check_polymerisation_pointed(actin))
                return true;
        }
    }
    return false;
}

bool ExcZone::circ_check_polymerisation_pointed(const Actin &actin)
{
    double pointedx = actin.getPointedEnd()[0];
    double pointedy = actin.getPointedEnd()[1];

    double d = sqrt( ( (pointedx-m_x)*(pointedx-m_x) + (pointedy-m_y)*(pointedy-m_y) ) );

    return (d < (actin.getStericRadius() + m_radius));
}

bool ExcZone::rect_check_polymerisation_pointed(const Actin &actin)
{

    gte::Segment<2,double> firstSub;
    firstSub.p[0][0] = actin.getPointedEnd()[0];
    firstSub.p[1][0] = actin.getPrePointedEnd()[0];
    firstSub.p[0][1] = actin.getPointedEnd()[1];
    firstSub.p[1][1] = actin.getPrePointedEnd()[1];

    gte::Segment<2,double> side1;
    side1.p[0][0] = m_x_min;
    side1.p[1][0] = m_x_min;
    side1.p[0][1] = m_y_min;
    side1.p[1][1] = m_y_max;

    SegDistQuery min_d_GTE;
    auto result = min_d_GTE(firstSub, side1);

    if (result.distance < actin.getStericRadius())
        return true;

    gte::Segment<2,double> side2;
    side2.p[0][0] = m_x_min;
    side2.p[1][0] = m_x_max;
    side2.p[0][1] = m_y_max;
    side2.p[1][1] = m_y_max;

    result = min_d_GTE(firstSub, side2);
    if (result.distance < actin.getStericRadius())
        return true;

    gte::Segment<2,double> side3;
    side3.p[0][0] = m_x_max;
    side3.p[1][0] = m_x_max;
    side3.p[0][1] = m_y_max;
    side3.p[1][1] = m_y_min;

    result = min_d_GTE(firstSub, side3);
    if (result.distance < actin.getStericRadius())
        return true;

    gte::Segment<2,double> side4;
    side4.p[0][0] = m_x_max;
    side4.p[1][0] = m_x_min;
    side4.p[0][1] = m_y_min;
    side4.p[1][1] = m_y_min;

    result = min_d_GTE(firstSub, side4);
    if (result.distance < actin.getStericRadius())
        return true;

    std::array<double,2> pointedEnd { actin.getPointedEnd()[0], actin.getPointedEnd()[1] };

    return rect_checkPointwithin(pointedEnd);
}

bool ExcZone::s_checkExVolMemFila(const Membrane &membrane,
                                  const std::vector<ExcZone> &excZones,
                                  int pointID)
{
    if (excZones.empty())
        return false;

    for (ExcZone excZone : excZones)
    {
        if (excZone.getCircBool())
        {
            if (excZone.checkExVolMemFilaCircSUB(membrane, pointID))
                return true;
        }
        else if (excZone.getRectBool())
        {
            if (excZone.checkExVolMemFilaRectSUB(membrane, pointID))
                return true;

        }
        else
        {
            // Triangle
            if (excZone.checkExVolMemFilaTriSUB(membrane, pointID))
                return true;
        }
    }
    return false;
}

bool ExcZone::checkExVolMemFilaCirc(const Membrane &membrane)
{
    /* Function to check excluded volume of a membrane filament against a
    circular exclusion zone

    Returns true if violated
    */

    // Define centre point of circle
    gte::Vector2<double> point = {m_x, m_y};

    // Define our membrane wall
    gte::Segment<2,double> membraneSeg;
    for (int i = 0; i < membrane.getNumPoints(); ++i)
    {
        int nxtSub = i+1;
        if (nxtSub == membrane.getNumPoints())
        {
            nxtSub = 0;
        }
        // The membrane
        membraneSeg.p[0][0] = membrane.getPoints()[i][0];
        membraneSeg.p[1][0] = membrane.getPoints()[nxtSub][0];
        membraneSeg.p[0][1] = membrane.getPoints()[i][1];
        membraneSeg.p[1][1] = membrane.getPoints()[nxtSub][1];

        DistQuery min_d_GTE;
        auto result = min_d_GTE(point, membraneSeg);
        double distance = result.distance;

        if (distance < (membrane.getThickness()/2 + m_radius))
        {
            return true;
        }
    }

    return false;
}

bool ExcZone::checkExVolMemFilaCirc(const std::vector<Membrane> &membranes)
{
    /* Function to check excluded volume of a membrane filament against a
    circular exclusion zone

    Returns true if violated
    */

    // Define centre point of circle
    gte::Vector2<double> point = {m_x, m_y};

    for (unsigned int j = 0; j < membranes.size(); ++j)
    {
        // Define our membrane
        gte::Segment<2,double> membraneSeg;
        for (int i = 0; i < membranes[j].getNumPoints(); ++i)
        {
            int nxtSub = i+1;
            if (nxtSub == membranes[j].getNumPoints())
            {
                nxtSub = 0;
            }
            // The membrane
            membraneSeg.p[0][0] = membranes[j].getPoints()[i][0];
            membraneSeg.p[1][0] = membranes[j].getPoints()[nxtSub][0];
            membraneSeg.p[0][1] = membranes[j].getPoints()[i][1];
            membraneSeg.p[1][1] = membranes[j].getPoints()[nxtSub][1];

            DistQuery min_d_GTE;
            auto result = min_d_GTE(point, membraneSeg);
            double distance = result.distance;

            if (distance < (membranes[j].getThickness()/2 + m_radius))
            {
                return true;
            }
        }
    }

    return false;
}

bool ExcZone::checkExVolMemFilaCircSUB(const Membrane &membrane, int pointID)
{
    /* Function to check excluded volume of a membrane filament against a
    circular exclusion zone

    Returns true if violated
    */

    // Define centre point of circle
    gte::Vector2<double> point = {m_x, m_y};

    // Define our membrane wall
    gte::Segment<2,double> membraneSeg;

    std::array<int,3> pointsToCheck;

    int start = pointID - 1;
    int end = pointID + 1;

    if (pointID == 0)
    {
        start = membrane.getNumPoints()-1;
    }
    else if (pointID == membrane.getNumPoints()-1)
    {
        end = 0;
    }

    pointsToCheck = {start, pointID, end};

    for (int i = 0; i < 2; ++i)
    {
        int p = pointsToCheck[i]; // point
        int np = pointsToCheck[i+1]; // next point
        // The membrane
        membraneSeg.p[0][0] = membrane.getPoints()[p][0];
        membraneSeg.p[1][0] = membrane.getPoints()[np][0];
        membraneSeg.p[0][1] = membrane.getPoints()[p][1];
        membraneSeg.p[1][1] = membrane.getPoints()[np][1];

        DistQuery min_d_GTE;
        auto result = min_d_GTE(point, membraneSeg);
        double distance = result.distance;
        if (distance < (membrane.getThickness()/2 + m_radius))
        {
            return true;
        }
    }

    return false;
}

bool ExcZone::checkExVolMemFilaRectSUB(const Membrane &membrane, int pointID)
{

    // Define membrane subunit
    gte::Segment<2,double> membraneSeg;

    std::array<int,3> pointsToCheck;

    int start = pointID - 1;
    int end = pointID + 1;

    if (pointID == 0)
    {
        start = membrane.getNumPoints()-1;
    }
    else if (pointID == membrane.getNumPoints()-1)
    {
        end = 0;
    }
    pointsToCheck = {start, pointID, end};
    for (int i = 0; i < 2; ++i)
    {
        // Membrane
        int p = pointsToCheck[i]; // point
        int np = pointsToCheck[i+1]; // next point
        // The membrane
        membraneSeg.p[0][0] = membrane.getPoints()[p][0];
        membraneSeg.p[1][0] = membrane.getPoints()[np][0];
        membraneSeg.p[0][1] = membrane.getPoints()[p][1];
        membraneSeg.p[1][1] = membrane.getPoints()[np][1];

        // Rectangle
        gte::Segment<2,double> side1;
        side1.p[0][0] = m_x_min;
        side1.p[1][0] = m_x_min;
        side1.p[0][1] = m_y_min;
        side1.p[1][1] = m_y_max;

        SegDistQuery min_d_GTE;
        auto result = min_d_GTE(membraneSeg, side1);

        if (result.distance < membrane.getThickness()/2)
            return true;

        gte::Segment<2,double> side2;
        side2.p[0][0] = m_x_min;
        side2.p[1][0] = m_x_max;
        side2.p[0][1] = m_y_max;
        side2.p[1][1] = m_y_max;

        result = min_d_GTE(membraneSeg, side2);
        if (result.distance < membrane.getThickness()/2)
            return true;

        gte::Segment<2,double> side3;
        side3.p[0][0] = m_x_max;
        side3.p[1][0] = m_x_max;
        side3.p[0][1] = m_y_max;
        side3.p[1][1] = m_y_min;

        result = min_d_GTE(membraneSeg, side3);
        if (result.distance < membrane.getThickness()/2)
            return true;

        gte::Segment<2,double> side4;
        side4.p[0][0] = m_x_max;
        side4.p[1][0] = m_x_min;
        side4.p[0][1] = m_y_min;
        side4.p[1][1] = m_y_min;

        result = min_d_GTE(membraneSeg, side4);
        if (result.distance < membrane.getThickness()/2)
            return true;


        // Check that it is not completely within
        if (rect_checkPointwithin(membrane.getPoints()[p]))
           return true;
    }

    return false;
}

bool ExcZone::checkExVolMemFilaTri(const std::vector<Membrane> &membranes)
{

    // Define triangle
    // Triangle
    gte::Segment<2,double> side1;
    side1.p[0][0] = m_p1[0];
    side1.p[1][0] = m_p2[0];
    side1.p[0][1] = m_p1[1];
    side1.p[1][1] = m_p2[1];

    gte::Segment<2,double> side2;
    side2.p[0][0] = m_p2[0];
    side2.p[1][0] = m_p3[0];
    side2.p[0][1] = m_p2[1];
    side2.p[1][1] = m_p3[1];

    gte::Segment<2,double> side3;
    side3.p[0][0] = m_p3[0];
    side3.p[1][0] = m_p1[0];
    side3.p[0][1] = m_p3[1];
    side3.p[1][1] = m_p3[1];

    // Define membrane subunit
    gte::Segment<2,double> membraneSeg;


    for (unsigned int j = 0; j < membranes.size(); ++j)
    {
        for (int i = 0; i < membranes[j].getNumPoints(); ++i)
        {
            // Membrane
            int nxtSub = i+1;
            if (nxtSub == membranes[j].getNumPoints())
            {
                nxtSub = 0;
            }

            // The membrane
            membraneSeg.p[0][0] = membranes[j].getPoints()[i][0];
            membraneSeg.p[1][0] = membranes[j].getPoints()[nxtSub][0];
            membraneSeg.p[0][1] = membranes[j].getPoints()[i][1];
            membraneSeg.p[1][1] = membranes[j].getPoints()[nxtSub][1];

            // Triangle
            SegDistQuery min_d_GTE;
            auto result = min_d_GTE(membraneSeg, side1);

            if (result.distance < membranes[j].getThickness()/2)
                return true;


            result = min_d_GTE(membraneSeg, side2);
            if (result.distance < membranes[j].getThickness()/2)
                return true;


            result = min_d_GTE(membraneSeg, side3);
            if (result.distance < membranes[j].getThickness()/2)
                return true;

        }
    }
    return false;
}

bool ExcZone::checkExVolMemFilaTriSUB(const Membrane &membrane, int pointID)
{

    // Define membrane subunit
    gte::Segment<2,double> membraneSeg;

    std::array<int,3> pointsToCheck;

    int start = pointID - 1;
    int end = pointID + 1;

    if (pointID == 0)
    {
        start = membrane.getNumPoints()-1;
    }
    else if (pointID == membrane.getNumPoints()-1)
    {
        end = 0;
    }
    pointsToCheck = {start, pointID, end};
    for (int i = 0; i < 2; ++i)
    {
        // Membrane
        int p = pointsToCheck[i]; // point
        int np = pointsToCheck[i+1]; // next point
        // The membrane
        membraneSeg.p[0][0] = membrane.getPoints()[p][0];
        membraneSeg.p[1][0] = membrane.getPoints()[np][0];
        membraneSeg.p[0][1] = membrane.getPoints()[p][1];
        membraneSeg.p[1][1] = membrane.getPoints()[np][1];

        // Triangle
        gte::Segment<2,double> side1;
        side1.p[0][0] = m_p1[0];
        side1.p[1][0] = m_p2[0];
        side1.p[0][1] = m_p1[1];
        side1.p[1][1] = m_p2[1];

        SegDistQuery min_d_GTE;
        auto result = min_d_GTE(membraneSeg, side1);

        if (result.distance < membrane.getThickness()/2)
            return true;

        gte::Segment<2,double> side2;
        side2.p[0][0] = m_p2[0];
        side2.p[1][0] = m_p3[0];
        side2.p[0][1] = m_p2[1];
        side2.p[1][1] = m_p3[1];

        result = min_d_GTE(membraneSeg, side2);
        if (result.distance < membrane.getThickness()/2)
            return true;

        gte::Segment<2,double> side3;
        side3.p[0][0] = m_p3[0];
        side3.p[1][0] = m_p1[0];
        side3.p[0][1] = m_p3[1];
        side3.p[1][1] = m_p3[1];

        result = min_d_GTE(membraneSeg, side3);
        if (result.distance < membrane.getThickness()/2)
            return true;

    }

    return false;
}


bool ExcZone::s_checkExVolCortex(const Cortex &cortex,
                                  const std::vector<ExcZone> &excZones,
                                  int pointID)
{
    if (excZones.empty())
        return false;

    for (ExcZone excZone : excZones)
    {
        if (excZone.getCircBool())
        {
            if (excZone.checkExVolCortexCircSUB(cortex, pointID))
                return true;
        }
        else if (excZone.getRectBool())
        {
            if (excZone.checkExVolCortexRectSUB(cortex, pointID))
                return true;
        }
        else
        {
            // Triangle
            if (excZone.checkExVolCortexTriSUB(cortex, pointID))
                return true;
        }
    }
    return false;
}

bool ExcZone::checkExVolCortexCircSUB(const Cortex &cortex, int pointID)
{
    /* Function to check excluded volume of a cortex against a
    circular exclusion zone

    Returns true if violated
    */

    // Define centre point of circle
    gte::Vector2<double> point = {m_x, m_y};

    // Define our cortex segment
    gte::Segment<2,double> cortexSeg;

    std::array<int,3> pointsToCheck;

    int start = pointID - 1;
    int end = pointID + 1;

    if (pointID == 0)
    {
        start = cortex.getNumPoints()-1;
    }
    else if (pointID == cortex.getNumPoints()-1)
    {
        end = 0;
    }

    pointsToCheck = {start, pointID, end};

    for (int i = 0; i < 2; ++i)
    {
        int p = pointsToCheck[i]; // point
        int np = pointsToCheck[i+1]; // next point
        // The cortex
        cortexSeg.p[0][0] = cortex.getPoints()[p][0];
        cortexSeg.p[1][0] = cortex.getPoints()[np][0];
        cortexSeg.p[0][1] = cortex.getPoints()[p][1];
        cortexSeg.p[1][1] = cortex.getPoints()[np][1];

        DistQuery min_d_GTE;
        auto result = min_d_GTE(point, cortexSeg);
        double distance = result.distance;
        if (distance < (cortex.getThickness()/2 + m_radius))
        {
            return true;
        }
    }

    return false;
}

bool ExcZone::checkExVolCortexRectSUB(const Cortex &cortex, int pointID)
{

    // Define membrane subunit
    gte::Segment<2,double> cortexSeg;

    std::array<int,3> pointsToCheck;

    int start = pointID - 1;
    int end = pointID + 1;

    if (pointID == 0)
    {
        start = cortex.getNumPoints()-1;
    }
    else if (pointID == cortex.getNumPoints()-1)
    {
        end = 0;
    }
    pointsToCheck = {start, pointID, end};
    for (int i = 0; i < 2; ++i)
    {
        // Membrane
        int p = pointsToCheck[i]; // point
        int np = pointsToCheck[i+1]; // next point
        // The cortex
        cortexSeg.p[0][0] = cortex.getPoints()[p][0];
        cortexSeg.p[1][0] = cortex.getPoints()[np][0];
        cortexSeg.p[0][1] = cortex.getPoints()[p][1];
        cortexSeg.p[1][1] = cortex.getPoints()[np][1];

        // Rectangle
        gte::Segment<2,double> side1;
        side1.p[0][0] = m_x_min;
        side1.p[1][0] = m_x_min;
        side1.p[0][1] = m_y_min;
        side1.p[1][1] = m_y_max;

        SegDistQuery min_d_GTE;
        auto result = min_d_GTE(cortexSeg, side1);

        if (result.distance < cortex.getThickness()/2)
            return true;

        gte::Segment<2,double> side2;
        side2.p[0][0] = m_x_min;
        side2.p[1][0] = m_x_max;
        side2.p[0][1] = m_y_max;
        side2.p[1][1] = m_y_max;

        result = min_d_GTE(cortexSeg, side2);
        if (result.distance < cortex.getThickness()/2)
            return true;

        gte::Segment<2,double> side3;
        side3.p[0][0] = m_x_max;
        side3.p[1][0] = m_x_max;
        side3.p[0][1] = m_y_max;
        side3.p[1][1] = m_y_min;

        result = min_d_GTE(cortexSeg, side3);
        if (result.distance < cortex.getThickness()/2)
            return true;

        gte::Segment<2,double> side4;
        side4.p[0][0] = m_x_max;
        side4.p[1][0] = m_x_min;
        side4.p[0][1] = m_y_min;
        side4.p[1][1] = m_y_min;

        result = min_d_GTE(cortexSeg, side4);
        if (result.distance < cortex.getThickness()/2)
            return true;


        // Check that it is not completely within
        if (rect_checkPointwithin(cortex.getPoints()[p]))
           return true;
    }

    return false;
}

bool ExcZone::checkExVolCortexTriSUB(const Cortex &cortex, int pointID)
{

    // Define membrane subunit
    gte::Segment<2,double> cortexSeg;

    std::array<int,3> pointsToCheck;

    int start = pointID - 1;
    int end = pointID + 1;

    if (pointID == 0)
    {
        start = cortex.getNumPoints()-1;
    }
    else if (pointID == cortex.getNumPoints()-1)
    {
        end = 0;
    }
    pointsToCheck = {start, pointID, end};
    for (int i = 0; i < 2; ++i)
    {
        // Membrane
        int p = pointsToCheck[i]; // point
        int np = pointsToCheck[i+1]; // next point
        // The cortex
        cortexSeg.p[0][0] = cortex.getPoints()[p][0];
        cortexSeg.p[1][0] = cortex.getPoints()[np][0];
        cortexSeg.p[0][1] = cortex.getPoints()[p][1];
        cortexSeg.p[1][1] = cortex.getPoints()[np][1];

        // Triangle
        gte::Segment<2,double> side1;
        side1.p[0][0] = m_p1[0];
        side1.p[1][0] = m_p2[0];
        side1.p[0][1] = m_p1[1];
        side1.p[1][1] = m_p2[1];

        SegDistQuery min_d_GTE;
        auto result = min_d_GTE(cortexSeg, side1);

        if (result.distance < cortex.getThickness()/2)
            return true;

        gte::Segment<2,double> side2;
        side2.p[0][0] = m_p2[0];
        side2.p[1][0] = m_p3[0];
        side2.p[0][1] = m_p2[1];
        side2.p[1][1] = m_p3[1];

        result = min_d_GTE(cortexSeg, side2);
        if (result.distance < cortex.getThickness()/2)
            return true;

        gte::Segment<2,double> side3;
        side3.p[0][0] = m_p3[0];
        side3.p[1][0] = m_p1[0];
        side3.p[0][1] = m_p3[1];
        side3.p[1][1] = m_p3[1];

        result = min_d_GTE(cortexSeg, side3);
        if (result.distance < cortex.getThickness()/2)
            return true;

    }

    return false;
}

bool ExcZone::checkExVolCircTar(const std::vector<ExcZone> excZones)
{
    // Check against all other targets (if any exist)

    assert(m_isCircle);
    int numExcZones = excZones.size();
    for (int i = 0; i < numExcZones; ++i)
    {
        if (i == m_id)
        {
            continue;
        }
        std::array<double, 2> p1 = {m_x, m_y};
        std::array<double, 2> p2 = {excZones[i].getX(), excZones[i].getY()};

        double dist = sqrt(distanceBetPoints2DSQR(p1, p2));
        if (dist < (m_radius+excZones[i].getR()))
        {
            return true;
        }
    }

    return false;
}

void ExcZone::targetBM(double temp, double viscosity, double dt,
                       const std::vector<Membrane> &membranes,
                       std::vector<ExcZone> &excZones,
                       std::vector<ProteinRegion> &branchRegions,
                       std::vector<ProteinRegion> &nucRegions,
                       std::vector<ProteinRegion> &aCapRegions,
                       std::vector<Actin> &actinVec)
{
    /*
    Function to move a circular target modelling Brownian motion
    Use Stokes Law D = (K*T)/(6*pi*eta*r)
    and then sqrt(<R^2>) = sqrt(2*d*D*dt)
    where little d is the dimensions, so do this separetely for X and Y
    Define a Normal distribution and then randomly choose from it
    */

    assert(m_isCircle);

    double D = (g_Kb*temp) / (6*M_PI*viscosity*m_radius);
    double stdDev = sqrt(2*dt*D);

    std::normal_distribution<double> tarBM(0,stdDev);

    double dx = tarBM(rng.m_mersenne);
    double dy = tarBM(rng.m_mersenne);

    dx = 0;

    m_x += dx;
    m_y += dy;

    // Check steric hindrance
    if (checkExVolCircTar(excZones) || checkExVolMemFilaCirc(membranes))
    {
        m_x -= dx;
        m_y -= dy;
    }
    else
    {
        // Extra check against actin
        for (unsigned int i = 0; i < actinVec.size(); ++i)
        {
          if (check_ex_vol(actinVec[i]))
          {
              m_x -= dx;
              m_y -= dy;
          }
        }
    }
}

void ExcZone::targetdirectedMovement(double dt,
                                     const std::vector<Membrane> &membranes,
                                     std::vector<ExcZone> &excZones,
                                     std::vector<ProteinRegion> &branchRegions,
                                     std::vector<ProteinRegion> &nucRegions,
                                     std::vector<ProteinRegion> &aCapRegions)
{
    /*
    Directed movement, deterministic with a velocity as if there is an external
    force pushing it in a direction
    */

    // Velocity, to start with lets have 400 nm/s
    std::array<double,2> velVector; // m/s
    velVector[0] = 0; // x component
    velVector[1] = -4E-7; // y component


    double dx = velVector[0]*dt;
    double dy = velVector[1]*dt;

    if (m_isCircle)
    {
        m_x += dx;
        m_y += dy;
        // Check steric hindrance
        if (checkExVolMemFilaCirc(membranes))
        {
            m_x -= dx;
            m_y -= dy;
        }
        else
        {
            for (unsigned int i = 0; i < branchRegions.size(); ++i)
            {
                if (!branchRegions[i].getCoupledBool())
                {
                    branchRegions[i].moveRegion(dx, dy);
                }
            }

            for (unsigned int i = 0; i < nucRegions.size(); ++i)
            {
                if (!nucRegions[i].getCoupledBool())
                {
                    nucRegions[i].moveRegion(dx, dy);
                }
            }

            for (unsigned int i = 0; i < aCapRegions.size(); ++i)
            {
                if (!aCapRegions[i].getCoupledBool())
                {
                    aCapRegions[i].moveRegion(dx, dy);
                }
            }
        }
    }
    else if (!m_isRectangle)
    {
        // Triangle
        m_p1[0] += dx;
        m_p2[0] += dx;
        m_p3[0] += dx;

        m_p1[1] += dy;
        m_p2[1] += dy;
        m_p3[1] += dy;

        // Check steric hindrance
        if (checkExVolMemFilaTri(membranes))
        {
            m_p1[0] -= dx;
            m_p2[0] -= dx;
            m_p3[0] -= dx;

            m_p1[1] -= dy;
            m_p2[1] -= dy;
            m_p3[1] -= dy;
        }
        else
        {
            for (unsigned int i = 0; i < branchRegions.size(); ++i)
            {
                if (!branchRegions[i].getCoupledBool())
                {
                    branchRegions[i].moveRegion(dx, dy);
                }
            }

            for (unsigned int i = 0; i < nucRegions.size(); ++i)
            {
                if (!nucRegions[i].getCoupledBool())
                {
                    nucRegions[i].moveRegion(dx, dy);
                }
            }

            for (unsigned int i = 0; i < aCapRegions.size(); ++i)
            {
                if (!aCapRegions[i].getCoupledBool())
                {
                    aCapRegions[i].moveRegion(dx, dy);
                }
            }
        }
    }

}

void ExcZone::s_targetdirectedMovementCircANDTri(double dt, bool &firstTouch,
                                                 std::array<double,2> velVector,
                                                 double indentDepth,
                                                 const std::vector<Membrane> &membranes,
                                                 std::vector<ExcZone> &excZones)
{
    /*
    Directed movement, deterministic with a velocity as if there is an external
    force pushing it in a direction

    MOVES CIRCLE AND TRIANGLE TOGETHER! BASED ON STERIC OF BOTH
    */

    double dx = velVector[0]*dt;
    double dy = velVector[1]*dt;
    // Move down
    for (unsigned int i = 0; i < excZones.size(); ++i)
    {

        if (excZones[i].getRectBool() || !excZones[i].getDirMove())
        {
            continue;
        }
        else if (excZones[i].getCircBool())
        {
            // Move it
            if (fabs(excZones[i].getTotalMoveDown()) >= indentDepth)
            {
                // Moved down far enough already!
                return;
            }
            excZones[i].addX(dx);
            excZones[i].addY(dy);
            if (firstTouch)
            {
                  // If the probe has touched the cell already
                excZones[i].addTotalMoveDown(dy);
            }
        }
        else
        {
            // Triangle
            if (fabs(excZones[i].getTotalMoveDown()) >= indentDepth)
            {
                // Moved down far enough already!
                return;
            }
            excZones[i].addTriX(dx);
            excZones[i].addTriY(dy);
            if (firstTouch)
            {
                  // If the probe has touched the cell already
                excZones[i].addTotalMoveDown(dy);
            }

        }
    }

    // Check to see if should reject the move
    for (unsigned int i = 0; i < excZones.size(); ++i)
    {
        if (excZones[i].getRectBool() || !excZones[i].getDirMove())
        {
            continue;
        }
        else if (excZones[i].getCircBool())
        {
            if (excZones[i].checkExVolMemFilaCirc(membranes))
            {
                for (unsigned int j = 0; j < excZones.size(); ++j)
                {
                    // Move all back
                    if (excZones[j].getRectBool() || !excZones[j].getDirMove())
                    {
                        continue;
                    }
                    else if (excZones[j].getCircBool())
                    {
                        // Move it
                        excZones[j].addX(-dx);
                        excZones[j].addY(-dy);
                        if (firstTouch)
                        {
                            excZones[j].addTotalMoveDown(-dy);
                        }

                    }
                    else
                    {
                        // Triangle
                        excZones[j].addTriX(-dx);
                        excZones[j].addTriY(-dy);
                        if (firstTouch)
                        {
                            excZones[j].addTotalMoveDown(-dy);
                        }
                    }
                }
                // If this is the first rejection then can now start adding total
                // move down
                firstTouch = true;
                break;
            }

        }
        else
        {
            // Triangle
            if (excZones[i].checkExVolMemFilaTri(membranes))
            {

                for (unsigned int j = 0; j < excZones.size(); ++j)
                {
                    // Move all back
                    if (excZones[j].getRectBool() || !excZones[j].getDirMove())
                    {
                        continue;
                    }
                    else if (excZones[j].getCircBool())
                    {
                        // Move it
                        excZones[j].addX(-dx);
                        excZones[j].addY(-dy);
                        if (firstTouch)
                        {
                            excZones[j].addTotalMoveDown(-dy);
                        }
                    }
                    else
                    {
                        // Triangle
                        excZones[j].addTriX(-dx);
                        excZones[j].addTriY(-dy);
                        if (firstTouch)
                        {
                            excZones[j].addTotalMoveDown(-dy);
                        }
                    }
                }
                // If this is the first rejection then can now start adding total
                // move down
                firstTouch = true;
                break;
            }
        }
    }

}

void ExcZone::move(double dx, double dy)
{
    assert(m_isCircle);
    m_x += dx;
    m_y += dy;
}
