/*
 *  phagoAnalysis.cpp
 *
 *
 *
 *  James Bradford
 *  University of Sheffield
 *  Dec 2018
 */

#include "configHeader.h"
#include "ExcZone.h"
#include "Membrane.h"
#include "GactinGrid.h"
#include "Actin.h"
#include "phagoAnalysis.h"
#include "geometry.h"

typedef gte::DCPQuery<double, gte::Vector<2,double>, gte::Segment<2, double>> DistQuery; // distance between point and segment

double calcEngulfment(std::vector<Membrane> &membranes,
                      std::vector<ExcZone> &targets)
{
    /*
     * Code that calculates the length of membrane in contact with the target,
     * Therefore its the wrapping at any given timepoint
     */


    // Define centre point of circle
    gte::Vector2<double> point = {targets[0].getX(), targets[0].getY()};

    // Define our membrane
    gte::Segment<2,double> membraneSeg;

    int contactSegs = 0;

    for (int j = 0; j < membranes[0].getNumPoints()-1; ++j)
    {
        // The membrane
        membraneSeg.p[0][0] = membranes[0].getPoints()[j][0];
        membraneSeg.p[1][0] = membranes[0].getPoints()[j+1][0];
        membraneSeg.p[0][1] = membranes[0].getPoints()[j][1];
        membraneSeg.p[1][1] = membranes[0].getPoints()[j+1][1];

        DistQuery min_d_GTE;
        auto result = min_d_GTE(point, membraneSeg);
        double distance = result.distance;

        if (distance <= (membranes[0].getThickness()/2 + targets[0].getR() + 1E-9))
        {
            // it is in contact, 1nm precision
            // 12 nm? bc of size of receptor?
            contactSegs += 1;
        }

    }


    // engulfLen in units of the circumference
    double engulfLen = (contactSegs * membranes[0].getSegLength())/(M_PI*2*targets[0].getR());
    return engulfLen;
}

int calcNumContactTips(std::vector<Membrane> &membranes,
                       std::vector<ExcZone> &targets, std::vector<Actin> &actins)
{
    /*
     * Code that calculates the number of tips (either barbed end or pointed end)
     * That are in contact with the target
     */


    // Define centre point of circle
    std::array<double,2> point = {targets[0].getX(), targets[0].getY()};


    int contactTips = 0;

    for (unsigned int i = 0; i < actins.size(); ++i)
    {

        double distToBeat = targets[0].getR() + actins[i].getStericRadius() + membranes[0].getThickness() + Actin::s_monomerLength;

        // for each actin we check both the barbed and the pointed end
        std::array<double,2> barbedEnd = {actins[i].getBarbedEnd()[0], actins[i].getBarbedEnd()[1]};
        if (sqrt(distanceBetPoints2DSQR(barbedEnd, point)) <= distToBeat)
        {
            contactTips += 1;
        }

        // Also check pointed end
        std::array<double,2> pointedEnd = {actins[i].getPointedEnd()[0], actins[i].getPointedEnd()[1]};
        if (sqrt(distanceBetPoints2DSQR(pointedEnd, point)) <= distToBeat)
        {
            contactTips += 1;
        }



    }

    return contactTips;
}


int calcNumContactTipsInRegion(std::vector<Membrane> &membranes,
                               std::vector<Actin> &actins, double regXmin,
                               double regXmax, double regYmin, double regYmax)
{
    /*
     * Code that calculates the number of tips (either barbed end or pointed end)
     * That are in contact with the membrane in a particular square region
     * Region will be related to the cup
     */




    int contactTips = 0;
    // Define our membrane
    gte::Segment<2,double> membraneSeg;

    for (unsigned int i = 0; i < actins.size(); ++i)
    {

        double distToBeat = actins[i].getStericRadius() + membranes[0].getThickness()/2 + Actin::s_monomerLength;

        // for each actin we check both the barbed and the pointed end
        // Define barbed end
        gte::Vector2<double> barbedEnd = {actins[i].getBarbedEnd()[0], actins[i].getBarbedEnd()[1]};

        // Define pointed end
        gte::Vector2<double> pointedEnd = {actins[i].getPointedEnd()[0], actins[i].getPointedEnd()[1]};

        if (barbedEnd[0] <= regXmax && barbedEnd[0] >= regXmin && barbedEnd[1] <= regYmax && barbedEnd[1] >= regYmin)
        {
            // it is in the region we have defined
            // is it in contact with the membrane?

            for (int j = 0; j < membranes[0].getNumPoints()-1; ++j)
            {
                // The membrane
                membraneSeg.p[0][0] = membranes[0].getPoints()[j][0];
                membraneSeg.p[1][0] = membranes[0].getPoints()[j+1][0];
                membraneSeg.p[0][1] = membranes[0].getPoints()[j][1];
                membraneSeg.p[1][1] = membranes[0].getPoints()[j+1][1];

                DistQuery min_d_GTE;
                auto result = min_d_GTE(barbedEnd, membraneSeg);
                double distance = result.distance;
                if (distance <= distToBeat)
                {
                    // it is in contact,
                    contactTips += 1;
                    break;
                }

            }

        }

        if (pointedEnd[0] <= regXmax && pointedEnd[0] >= regXmin && pointedEnd[1] <= regYmax && pointedEnd[1] >= regYmin)
        {
            for (int j = 0; j < membranes[0].getNumPoints()-1; ++j)
            {
                // The membrane
                membraneSeg.p[0][0] = membranes[0].getPoints()[j][0];
                membraneSeg.p[1][0] = membranes[0].getPoints()[j+1][0];
                membraneSeg.p[0][1] = membranes[0].getPoints()[j][1];
                membraneSeg.p[1][1] = membranes[0].getPoints()[j+1][1];

                DistQuery min_d_GTE;
                auto result = min_d_GTE(pointedEnd, membraneSeg);
                double distance = result.distance;

                if (distance <= distToBeat)
                {
                    // it is in contact,
                    contactTips += 1;
                    break;
                }

            }
        }
    }
    return contactTips;
}
