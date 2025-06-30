/*
 *  MembraneWall.cpp
 *
 *  C++ file containing the definition of the Membrane-Wall class.
 *  This contains member variable initialisations and non-trivial member
 *  functions.
 *  The declaration of this class (including the definitions of the trivial
 *  member functions) is found in the associated header file "MembraneWall.h"
 *
 *  James Bradford
 *  University of Sheffield
 *  Nov 2017
 */

#include "configHeader.h"
#include "Actin.h"
#include "MembraneWall.h"

typedef gte::DCPQuery<double, gte::Segment<2, double>, gte::Segment<2, double>> SegDistQuery; // distance between two segments

MembraneWall::MembraneWall(double ypos, double length, double extForce)
    : m_id { 0 },
    m_exist { true },
    m_ypos { ypos },
    m_length { length }, // for now centre at 0
    m_extForce { extForce },
    m_origYPos { ypos }
{
}

double MembraneWall::distToMoveWall_barbed(Actin &actin)
{
    // Returns the distance the membrane has to move, d_i in the Indian paper
    std::array<double,2> unitVec = actin.getUnitVec(actin.getNumSubs()-1);
    double newBarby = actin.getBarbedEnd()[1] + (Actin::s_monomerLength * unitVec[1]);


    if (newBarby + Actin::s_stericRadius > m_ypos)
    {
        // Membrane has to move
        // Include the steric radius here so the membrane touches the end of the filament
        return (newBarby + Actin::s_stericRadius - m_ypos);
    }
    else
    {
        // membrane doesn't move
        return 0;
    }
}

double MembraneWall::distToMoveWall_pointed(Actin &actin)
{
    // Returns the distance the membrane has to move, d_i in the Indian paper
    std::array<double,2> unitVec = actin.getUnitVec(0);
    double newPointedy = actin.getPointedEnd()[1] - (Actin::s_monomerLength * unitVec[1]);


    if (newPointedy + Actin::s_stericRadius > m_ypos)
    {
        // Membrane has to move
        // Include the steric radius here so the membrane touches the end of the filament
        return (newPointedy + Actin::s_stericRadius - m_ypos);
    }
    else
    {
        // membrane doesn't move
        return 0;
    }
}


void MembraneWall::moveCoupledRegions(double dx, double dy,
                                     std::vector<ProteinRegion> &branchRegions,
                                     std::vector<ProteinRegion> &nucRegions,
                                     std::vector<ProteinRegion> &capRegions,
                                     std::vector<ProteinRegion> &sevRegions)
{
    for (unsigned int i = 0; i < m_coupledBranchRegionsIDs.size(); ++i)
    {
        branchRegions[m_coupledBranchRegionsIDs[i]].moveRegion(dx,dy);
    }
    for (unsigned int i = 0; i < m_coupledNucRegionsIDs.size(); ++i)
    {
        nucRegions[m_coupledNucRegionsIDs[i]].moveRegion(dx,dy);
    }
    for (unsigned int i = 0; i < m_coupledCapRegionsIDs.size(); ++i)
    {
        capRegions[m_coupledCapRegionsIDs[i]].moveRegion(dx,dy);
    }
    for (unsigned int i = 0; i < m_coupledSevRegionsIDs.size(); ++i)
    {
        sevRegions[m_coupledSevRegionsIDs[i]].moveRegion(dx,dy);
    }


}

bool MembraneWall::checkExVolNuc(const Actin &actin)
{
    // Function that carries out a distance check between the new actin filament
    // either a branch or free trimer, and the membrane wall

    // Returns true if violated

    // Define our wall
    gte::Segment<2,double> wall;
    wall.p[0][0] = -m_length/2;
    wall.p[1][0] = m_length/2;
    wall.p[0][1] = m_ypos;
    wall.p[1][1] = m_ypos;

    SegDistQuery min_d_GTE;

    // Define our filament
    gte::Segment<2,double> filament;

    filament.p[0][0] = actin.getPoints()[0][0];
    filament.p[1][0] = actin.getPoints()[4][0];
    filament.p[0][1] = actin.getPoints()[0][1];
    filament.p[1][1] = actin.getPoints()[4][1];

    // Do the distance query
    auto result = min_d_GTE(wall,filament);
    if (result.distance < actin.getStericRadius())
    {
        return true;
    }
    else
    {
        return false;
    }
}

bool MembraneWall::checkExVolBM(const Actin &actin) const
{
    // Function that checks whether any point along the actin filament exists
    // beyond the membrane wall

    // add a small error
    double err = 1E-14;
    for (int i = 0; i <= actin.getNumSubs(); ++i)
    {
        // just check y dimension right now
        if (actin.getPoints()[i][1] + Actin::s_stericRadius > m_ypos + err)
        {
            return true;
        }
    }
    return false;
}

double MembraneWall::distToMoveWallBack(const std::vector<Actin> &actinvec)
{
    // Function that calculates the distances between the wall and finds the min
    // if any of the min is 0 it can break out

    // Not currently parralised but maybe this can be?

    // Define our wall

    if (actinvec.size() == 0)
    {
        // Special case: there's no filaments, do not move the wall
        return 0;
    }

    gte::Segment<2,double> wall;
    wall.p[0][0] = -m_length/2;
    wall.p[1][0] = m_length/2;
    wall.p[0][1] = m_ypos;
    wall.p[1][1] = m_ypos;

    SegDistQuery min_d_GTE;
    // Define our subfilament
    gte::Segment<2,double> subFilament;

    double minDist = 1E-6; // maximum the membrane can move in one timestep

    for (unsigned int i = 0; i < actinvec.size(); ++i)
    {
        for (int j = 0; j < actinvec[i].getNumSubs(); ++j)
        {
            subFilament.p[0][0] = actinvec[i].getPoints()[j][0];
            subFilament.p[1][0] = actinvec[i].getPoints()[j+1][0];
            subFilament.p[0][1] = actinvec[i].getPoints()[j][1];
            subFilament.p[1][1] = actinvec[i].getPoints()[j+1][1];

            auto result = min_d_GTE(wall, subFilament);

            if (result.distance  == Actin::s_stericRadius) // Double comparision - Never going to be true, I must have been on something when i wrote this
            {
                return 0;
            }

            minDist = (result.distance - Actin::s_stericRadius < minDist) ? (result.distance - Actin::s_stericRadius) : minDist;


        }
    }
    return minDist;
}

double MembraneWall::distToMoveWallBack2(const std::vector<Actin> &actinvec,
                                         StericGrid &stericGrid)
{

    std::vector<int> cellsToCheck;
    cellsToCheck = stericGrid.getCellsToCheckMemWall(m_stericCells);

    return calc_min_dist_Grid(actinvec, cellsToCheck, stericGrid);

}

double MembraneWall::calc_min_dist_Grid(const std::vector<Actin> &actinvec,
                                        std::vector<int> &cellsToCheck,
                                        StericGrid &stericGrid)
{

    gte::Segment<2,double> wall;
    wall.p[0][0] = -m_length/2;
    wall.p[1][0] = m_length/2;
    wall.p[0][1] = m_ypos;
    wall.p[1][1] = m_ypos;


    // To avoid repeated checks against the same subunits that may appear in
    // more than one cell, have a recorded of which ones we have checked
    std::vector<std::array<int,2>> checkedSubs;
    double minDist = 1E-6;

    for (unsigned int i = 0; i < cellsToCheck.size(); ++i)
    {

        std::vector<std::array<int,2>> cellContents = stericGrid.getCellContents(cellsToCheck[i]);

        for (unsigned int j = 0; j < cellContents.size(); ++j)
        {
            int actinID = cellContents[j][0];

            if (std::find(checkedSubs.begin(), checkedSubs.end(), cellContents[j]) != checkedSubs.end())
            {
                // We have checked this filament before!
                continue;
            }
            checkedSubs.push_back(cellContents[j]);

            int actinSubID = cellContents[j][1];
            gte::Segment<2,double> filament;

            filament.p[0][0] = actinvec[actinID].getPoints()[actinSubID][0];
            filament.p[1][0] = actinvec[actinID].getPoints()[actinSubID+1][0];
            filament.p[0][1] = actinvec[actinID].getPoints()[actinSubID][1];
            filament.p[1][1] = actinvec[actinID].getPoints()[actinSubID+1][1];

            SegDistQuery min_d_GTE;
            auto result = min_d_GTE(wall, filament);

            minDist = (result.distance - Actin::s_stericRadius < minDist) ? (result.distance - Actin::s_stericRadius) : minDist;

        }
    }

    return minDist;
}

void MembraneWall::moveWall(double deltaD, StericGrid &stericGrid)
{
    m_ypos += deltaD;
    stericGrid.resetMemWallandUpdate(*this);
}
