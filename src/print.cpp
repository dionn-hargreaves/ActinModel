/*
 *  print.cpp
 *
 *  C++ file containing the printing functions for printing out to file.
 *
 *  The declaration of these functions are found in the associated header file
 *  "print.h".
 *
 *  James Bradford
 *  University of Sheffield
 *  Jan 2017
 */

#include "Actin.h"
#include "RNG.h"
#include "ProteinRegion.h"
#include "ExcZone.h"
#include "MembraneWall.h"
#include "GactinGrid.h"
#include "ArpGrid.h"
#include "print.h"
#include "Cortex.h"

extern RNG rng;

void printActinframes(std::ofstream &outf, int nActin,
                      const std::vector<Actin> &actinVec, double currTime,
                      double dt, const std::vector<ExcZone> &excZones,
                      const std::vector<ProteinRegion> &nucRegions,
                      const std::vector<ProteinRegion> &branchRegions,
                      const std::vector<ProteinRegion> &capRegions,
                      const std::vector<ProteinRegion> &antiCapRegions,
                      const std::vector<ProteinRegion> &sevRegions,
                      const std::vector<MembraneWall> &memWalls,
                      const std::vector<Membrane> &membranes,
                      bool tether, bool crossLinking, const Cortex &cortex)
{
    //outf.precision(17);
    outf.precision(10);
    outf << "Time (s): " << currTime << "\n";
    if (excZones.size() > 0)
    {
        outf << "# Exclusion regions: " << "\n";
        outf << "# x pos, y pos, radius" << "\n";
    }
    for (unsigned int i = 0; i < excZones.size(); ++i)
    {
        if (excZones[i].getCircBool())
        {
            // Circle
            outf << "!C, " << excZones[i].getX() << ", " << excZones[i].getY() << ", ";
            outf << excZones[i].getR() << std::endl;
        }
        else if (excZones[i].getRectBool())
        {
            // Rectangle
            outf << "!R, " << excZones[i].getXmin() << ", " << excZones[i].getYmin() << ", ";
            outf << excZones[i].getXmax() << ", " << excZones[i].getYmax() << std::endl;
        }
        else
        {
            //Triangle
            outf << "!TR, " << excZones[i].getP1()[0] << ", " << excZones[i].getP1()[1] << ", ";
            outf << excZones[i].getP2()[0] << ", " << excZones[i].getP2()[1] << ", ";
            outf << excZones[i].getP3()[0] << ", " << excZones[i].getP3()[1] << std::endl;
        }
    }
    if (nucRegions.size() > 0)
    {
        outf << "# NucRegions: " << "\n";
        outf << "# xmin, ymin, xmax, ymax" << "\n";
    }
    for (unsigned int i = 0; i < nucRegions.size(); ++i)
    {
        if (nucRegions[i].getCircBool())
        {
            outf << "!NC, " << nucRegions[i].getCentre()[0] << ", ";
            outf << nucRegions[i].getCentre()[1] << ", ";
            outf << nucRegions[i].getR() << std::endl;
        }
        else if (nucRegions[i].getRingBool())
        {
            outf << "!NR, " << nucRegions[i].getCentre()[0] << ", ";
            outf << nucRegions[i].getCentre()[1] << ", ";
            outf << nucRegions[i].getInR() << ", ";
            outf << nucRegions[i].getR() << std::endl;
        }
        else
        {
            std::array<double,8> cornersCoords = nucRegions[i].getCorners();
            outf << "!N, " << cornersCoords[0] << ", " << cornersCoords[1] << ", ";
            outf << cornersCoords[2] << ", " << cornersCoords[3] << ", ";
            outf << cornersCoords[4] << ", " << cornersCoords[5] << ", ";
            outf << cornersCoords[6] << ", " << cornersCoords[7] << std::endl;
        }
    }
    if (branchRegions.size() > 0)
    {
        outf << "# ArpRegions: " << "\n";
        outf << "# xmin, ymin, xmax, ymax" << "\n";
    }
    for (unsigned int i = 0; i < branchRegions.size(); ++i)
    {
        if (branchRegions[i].getCircBool())
        {
            outf << "!BC, " << branchRegions[i].getCentre()[0] << ", ";
            outf << branchRegions[i].getCentre()[1] << ", ";
            outf << branchRegions[i].getR() << std::endl;
        }
        else
        {
            std::array<double,8> cornersCoords = branchRegions[i].getCorners();
            outf << "!B, " << cornersCoords[0] << ", " << cornersCoords[1] << ", ";
            outf << cornersCoords[2] << ", " << cornersCoords[3] << ", ";
            outf << cornersCoords[4] << ", " << cornersCoords[5] << ", ";
            outf << cornersCoords[6] << ", " << cornersCoords[7] << std::endl;
        }
    }
    if (capRegions.size() > 0)
    {
        outf << "# CapRegions: " << "\n";
        outf << "# xmin, ymin, xmax, ymax" << "\n";
    }
    for (unsigned int i = 0; i < capRegions.size(); ++i)
    {
        if (capRegions[i].getCircBool())
        {
            outf << "!CapC, " << capRegions[i].getCentre()[0] << ", ";
            outf << capRegions[i].getCentre()[1] << ", ";
            outf << capRegions[i].getR() << std::endl;
        }
        else
        {
            std::array<double,8> cornersCoords = capRegions[i].getCorners();
            outf << "!Cap, " << cornersCoords[0] << ", " << cornersCoords[1] << ", ";
            outf << cornersCoords[2] << ", " << cornersCoords[3] << ", ";
            outf << cornersCoords[4] << ", " << cornersCoords[5] << ", ";
            outf << cornersCoords[6] << ", " << cornersCoords[7] << std::endl;
        }
    }
    if (antiCapRegions.size() > 0)
    {
        outf << "# AntiCapRegions: " << "\n";
        outf << "# xmin, ymin, xmax, ymax" << "\n";
    }
    for (unsigned int i = 0; i < antiCapRegions.size(); ++i)
    {
        if (antiCapRegions[i].getCircBool())
        {
            outf << "!ACapC, " << antiCapRegions[i].getCentre()[0] << ", ";
            outf << antiCapRegions[i].getCentre()[1] << ", ";
            outf << antiCapRegions[i].getR() << std::endl;
        }
        else
        {
            std::array<double,8> cornersCoords = antiCapRegions[i].getCorners();
            outf << "!ACap, " << cornersCoords[0] << ", " << cornersCoords[1] << ", ";
            outf << cornersCoords[2] << ", " << cornersCoords[3] << ", ";
            outf << cornersCoords[4] << ", " << cornersCoords[5] << ", ";
            outf << cornersCoords[6] << ", " << cornersCoords[7] << std::endl;
        }
    }

    for (unsigned int i = 0; i < sevRegions.size(); ++i)
    {
        if (sevRegions[i].getCircBool())
        {
            std::cout << "Circular sev region error" << std::endl;
            exit(1);
            outf << "!S, " << sevRegions[i].getCentre()[0] << ", ";
            outf << sevRegions[i].getCentre()[1] << ", ";
            outf << sevRegions[i].getR() << std::endl;
        }
        else
        {
            std::array<double,8> cornersCoords = sevRegions[i].getCorners();
            outf << "!S, " << cornersCoords[0] << ", " << cornersCoords[1] << ", ";
            outf << cornersCoords[2] << ", " << cornersCoords[3] << ", ";
            outf << cornersCoords[4] << ", " << cornersCoords[5] << ", ";
            outf << cornersCoords[6] << ", " << cornersCoords[7] << std::endl;
        }
    }
    if (memWalls.size() > 0)
    {
        outf << "# Membrane walls: " << "\n";
        outf<< "# ypos, xlength, extForce" << "\n";
    }
    for (unsigned int i = 0; i < memWalls.size(); ++i)
    {
        outf << "!M, " << memWalls[i].getYpos() << ", " << memWalls[i].getXlength();
        outf << ", " << memWalls[i].getForce() << std::endl;
    }

    if (membranes[0].getExist())
    {
        outf << "# Membrane filament points: " << "\n";
        outf<< "# x pos, ypos" << "\n";
    }

    for (unsigned int i = 0; i < membranes.size(); ++i)
    {
        for (unsigned int j = 0; j < membranes[i].getNumPoints(); ++j)
        {
            if (j == membranes[i].getNumPoints()-1)
            {
                outf << "!MF, " << membranes[i].getPoints()[j][0] << ", " << membranes[i].getPoints()[j][1];
                outf << ", " << membranes[i].getPoints()[0][0] << ", " << membranes[i].getPoints()[0][1];
            }
            else
            {
                outf << "!MF, " << membranes[i].getPoints()[j][0] << ", " << membranes[i].getPoints()[j][1];
                outf << ", " << membranes[i].getPoints()[j+1][0] << ", " << membranes[i].getPoints()[j+1][1];
            }

            outf << "\n";
        }
    }

    if (cortex.getExist())
    {
        outf << "# Cortex points: " << "\n";
        outf<< "# x pos, ypos" << "\n";

        for (unsigned int j = 0; j < cortex.getNumPoints(); ++j)
        {
            if (j == cortex.getNumPoints()-1)
            {
                outf << "!CX, " << cortex.getPoints()[j][0] << ", " << cortex.getPoints()[j][1];
                outf << ", " << cortex.getPoints()[0][0] << ", " << cortex.getPoints()[0][1];
            }
            else
            {
                outf << "!CX, " << cortex.getPoints()[j][0] << ", " << cortex.getPoints()[j][1];
                outf << ", " << cortex.getPoints()[j+1][0] << ", " << cortex.getPoints()[j+1][1];
            }

            outf << "\n";
        }
    }



    for (unsigned int i = 0; i < actinVec.size(); ++i)
    {
        for (unsigned int j = 0; j < actinVec[i].getNumTethers(); ++j)
        {
            outf << "!TP, " << actinVec[i].getTetherPoints(j)[0] << ", ";
            outf << actinVec[i].getTetherPoints(j)[1] << ", ";
            outf << actinVec[i].getTetherPoints(j)[2] << ", ";
            outf << actinVec[i].getTetherPoints(j)[3] << std::endl;
        }
    }


    for (unsigned int i = 0; i < actinVec.size(); ++i)
    {
        for (int j = 0; j < actinVec[i].getNumCLinks(); ++j)
        {
            outf << "!CL, " << actinVec[i].getCLinksPoint(j)[0] << ", ";
            outf << actinVec[i].getCLinksPoint(j)[1] << ", ";
            outf << actinVec[i].getCLinksPoint(j)[2] << ", ";
            outf << actinVec[i].getCLinksPoint(j)[3] << std::endl;
        }
    }


    if (cortex.getExist())
    {
        int numTethers = membranes[0].getTetherDists().size();

        for (int i = 0; i < numTethers; ++i)
        {
            int j = membranes[0].getTetherSubIDs()[i];
            double a = membranes[0].getTetherDists()[i];
            assert (a < 1);
            double tetherA_x = membranes[0].getPoints()[j][0] + membranes[0].getUnitVecs()[j][0]*a*membranes[0].getPresSubLengths()[j];
            double tetherA_y = membranes[0].getPoints()[j][1] + membranes[0].getUnitVecs()[j][1]*a*membranes[0].getPresSubLengths()[j];


            int k = cortex.getTetherSubIDs()[i];
            double b = cortex.getTetherDists()[i];
            assert (b < 1);
            double tetherB_x = cortex.getPoints()[k][0] + cortex.getUnitVecs()[k][0]*b*cortex.getPresSubLengths()[k];
            double tetherB_y = cortex.getPoints()[k][1] + cortex.getUnitVecs()[k][1]*b*cortex.getPresSubLengths()[k];


            outf << "!TP, " << tetherA_x << ", ";
            outf << tetherA_y << ", " << tetherB_x << ", " << tetherB_y << std::endl;

        }
    }


    outf << "nActin: " << nActin << " nSubunits: " << Actin::s_total_subunits << "\n";
    outf << "# actinID, SubunitID, birthtime (s), xStart (m), yStart (m), ";
    outf << "xEnd (m), yEnd (m), pointed capped, ";
    outf << "barbed capped, branch" << std::endl;

    int subunitid { 0 };
    for (Actin actinobj : actinVec)
    {
        for (int i = 0; i < actinobj.getNumSubs(); ++i)
        {
            outf << actinobj.getID() << ", " << subunitid << ", ";
            outf << actinobj.getbirthTime() << ", " << actinobj.getPoints()[i][0];
            outf << ", " << actinobj.getPoints()[i][1] << ", ";
            outf << actinobj.getPoints()[i+1][0] << ", ";
            outf << actinobj.getPoints()[i+1][1] << ", ";

            if (i == 0)
            {
                // We are on the first subunit
                outf << actinobj.getPointedCapped() << ", ";
            }
            else
            {
                // Not on the first subunit
                outf << 0 << ", ";
            }

            if (i == (actinobj.getNumSubs()-1))
            {
                // We are on the last subunit
                outf << actinobj.getBarbedCapped() << ", ";
            }
            else
            {
                // Not on the last subunit
                outf << 0 << ", ";
            }

            if ( (i==0) && (actinobj.getParentID() != -1) )
            {
                // First subunit and it is a branched filament
                outf << 1 << std::endl;
            }
            else
            {
                outf << 0 << std::endl;
            }

            ++subunitid;
        }
    }
    outf << "#--------------------------------------------------" << std::endl;
}

void printActinHeaders(std::ofstream &outf, double dt_bw_f)
{
    outf << "Time between frames (s): " << dt_bw_f << std::endl;
}

void printLog(std::ofstream &outf, double dt, std::uint32_t origSeed, double fracErr)
{
    auto t = std::time(nullptr);
    //auto tm = *std::localtime(&t);
    char datetmp[100];
    std::strftime(datetmp, 100, "%H:%M:%S %d/%m/%Y", std::localtime(&t));
    std::string date (datetmp);

    outf << "#### Simulation run details ####" << "\n";
    outf << "Date and time of run: " << date << "\n";
    outf << "Random number seed: " << origSeed << "\n";
    outf << "Num of OMP threads: " << omp_get_max_threads() << "\n";
    outf << "\n";
    outf << "Fractional error used: " << fracErr << "\n";
    outf << "Time step (s): " << dt << std::endl;;

}

void printRunTime(std::ofstream &outf, double runTime)
{
    outf << "Time taken to finish sim (s): " << runTime << std::endl;
}

void printGActinGridHeader(std::ofstream &outf, double dt_bw_f, const GactinGrid &gActinGrid)
{
    outf << "Time between frames (s): " << dt_bw_f << std::endl;
    outf << "Dimensions of a cell (x,y,z): ";
    outf << gActinGrid.getCellWidth() << ", " << gActinGrid.getCellHeight();
    outf << ", " << gActinGrid.getCellDepth() << std::endl;
}

void printGActinGrid(std::ofstream &outf, const GactinGrid &gActinGrid,
                     double currTime)
{
    outf << "#--------------------------------------------------" << "\n";
    outf << "Time (s): " << currTime << "\n";
    for (int j = gActinGrid.getNumCellsY()-1; j >=0; --j) // height, start at the top
    {
        for (int i = 0; i < gActinGrid.getNumCellsX(); ++i) // width, start at the left
        {
            // Conversion for grid size of 250nm
            //m_num_vec[i][j] = round(m_conc_vec[i][j] * 9.409);
            if (i == gActinGrid.getNumCellsX()-1)
                outf << gActinGrid.getNumber(i,j);
            else
                outf << gActinGrid.getNumber(i,j) << ',';
        }
        outf << std::endl;
    }
}

void printLatAGrid(std::ofstream &outf, const GactinGrid &gActinGrid,
                     double currTime)
{
    outf << "#--------------------------------------------------" << "\n";
    outf << "Time (s): " << currTime << "\n";
    for (int j = gActinGrid.getNumCellsY()-1; j >=0; --j) // height, start at the top
    {
        for (int i = 0; i < gActinGrid.getNumCellsX(); ++i) // width, start at the left
        {
            // Conversion for grid size of 250nm
            if (i == gActinGrid.getNumCellsX()-1)
                outf << gActinGrid.getLatANumber(i,j);
            else
                outf << gActinGrid.getLatANumber(i,j) << ',';
        }
        outf << std::endl;
    }
}

void printLatActinGrid(std::ofstream &outf, const GactinGrid &gActinGrid,
                     double currTime)
{
    outf << "#--------------------------------------------------" << "\n";
    outf << "Time (s): " << currTime << "\n";
    for (int j = gActinGrid.getNumCellsY()-1; j >=0; --j) // height, start at the top
    {
        for (int i = 0; i < gActinGrid.getNumCellsX(); ++i) // width, start at the left
        {
            // Conversion for grid size of 250nm
            if (i == gActinGrid.getNumCellsX()-1)
                outf << gActinGrid.getLatActinNumber(i,j);
            else
                outf << gActinGrid.getLatActinNumber(i,j) << ',';
        }
        outf << std::endl;
    }
}

void printArpGridHeader(std::ofstream &outf, double dt_bw_f, const ArpGrid &arpGrid)
{
    outf << "Time between frames (s): " << dt_bw_f << std::endl;
    outf << "Dimensions of a cell (x,y,z): ";
    outf << arpGrid.getCellWidth() << ", " << arpGrid.getCellHeight();
    outf << ", " << arpGrid.getCellDepth() << std::endl;
}

void printArpGrid(std::ofstream &outf, const ArpGrid &arpGrid,
                  double currTime)
{
    outf << "#--------------------------------------------------" << "\n";
    outf << "Time (s): " << currTime << "\n";
    for (int j = arpGrid.getNumCellsY()-1; j >=0; --j) // height, start at the top
    {
        for (int i = 0; i < arpGrid.getNumCellsX(); ++i) // width, start at the left
        {
            // Conversion for grid size of 250nm
            if (i == arpGrid.getNumCellsX()-1)
                outf << arpGrid.getNumber(i,j);
            else
                outf << arpGrid.getNumber(i,j) << ',';
        }
        outf << std::endl;
    }
}

void printMemLenHeader(std::ofstream &outf, double dt_bw_f)
{
    outf << "Time between frames (s): " << dt_bw_f << std::endl;
    outf << "#--------------------------------------------------" << std::endl;
}

void printMemLen(std::ofstream &outf, const std::vector<Membrane> &membranes, double currTime)
{
    outf << "Time (s): " << currTime << "\n";
    for (unsigned int i = 0; i < membranes.size(); ++i)
    {
        outf << "Mem Id: " << i << "\n";
        for (int j = 0; j < membranes[i].getNumPoints(); ++j)
        {
            outf << membranes[i].getActualSubLengths()[j] << std::endl;
        }
    }
    outf << "#--------------------------------------------------" << std::endl;
}

void printAnglesHeader(std::ofstream &outf, double dt_bw_f)
{
    outf << "Time between frames (s): " << dt_bw_f << std::endl;
}

void printAngles(std::ofstream &outf, std::vector<Actin> &actinvec)
{
    for (unsigned int i = 0; i < actinvec.size(); ++i)
    {

        int centreSub;
        actinvec[i].findCentrePoint(centreSub);
        // ignore the ends?
        for (int j = 1; j < actinvec[i].getNumSubs()-1; ++j)
        {

            if (j < centreSub)
            {
                double globAngle = actinvec[i].getSumAngle(j, centreSub) + atan2(-actinvec[i].getUnitVec(centreSub)[1], -actinvec[i].getUnitVec(centreSub)[0]);
                double preAngle = actinvec[i].getSumAngle(j+1, centreSub) + atan2(-actinvec[i].getUnitVec(centreSub)[1], -actinvec[i].getUnitVec(centreSub)[0]);
                double relAngle = preAngle - globAngle;
                relAngle = relAngle - 2*M_PI*floor(relAngle/(2*M_PI)+0.5);
                outf << relAngle << ", ";

            }
            else if (j > centreSub)
            {
                double globAngle = actinvec[i].getSumAngle(j, centreSub) + atan2(actinvec[i].getUnitVec(centreSub)[1], actinvec[i].getUnitVec(centreSub)[0]);
                double preAngle = actinvec[i].getSumAngle(j-1, centreSub) + atan2(actinvec[i].getUnitVec(centreSub)[1], actinvec[i].getUnitVec(centreSub)[0]);
                double relAngle = globAngle - preAngle;
                relAngle = relAngle - 2*M_PI*floor(relAngle/(2*M_PI)+0.5);
                if (j == actinvec[i].getNumSubs()-2)
                {
                    outf << relAngle << "\n";
                }
                else
                {
                    outf << relAngle << ", ";
                }
            }

        }
    }
    outf << "#--------------------------------------------------" << std::endl;;

}

void printEngulfment(std::ofstream &outf, double currTime, double engulfment)
{
    outf << currTime << ", " << engulfment << std::endl;
}
void printContact(std::ofstream &outf, double currTime, int numContact)
{
    outf << currTime << ", " << numContact << std::endl;
}

void printLength(std::ofstream &outf, std::vector<Actin> &actinvec, double currTime)
{
    // totals up and prints the length of all the F-actin in system
    double totalLen = 0;
    for (unsigned int i = 0; i < actinvec.size(); ++i)
    {
        totalLen += actinvec[i].getLength();
    }

    outf << currTime << ", " << totalLen << std::endl;
}

void printMonoTimeFrames(std::ofstream &outf, int nActin,
                         const std::vector<Actin> &actinVec, double currTime,
                         double dt)
{
    //outf.precision(17);
    outf.precision(10);
    outf << "Time (s): " << currTime << "\n";

    for (unsigned int i = 0; i < actinVec.size(); ++i)
    {
        for (unsigned int j = 0; j < actinVec[i].getNumTethers(); ++j)
        {
            outf << "!TP, " << actinVec[i].getTetherPoints(j)[0] << ", ";
            outf << actinVec[i].getTetherPoints(j)[1] << ", ";
            outf << actinVec[i].getTetherPoints(j)[2] << ", ";
            outf << actinVec[i].getTetherPoints(j)[3] << std::endl;
        }
    }


    for (unsigned int i = 0; i < actinVec.size(); ++i)
    {
        for (int j = 0; j < actinVec[i].getNumCLinks(); ++j)
        {
            outf << "!CL, " << actinVec[i].getCLinksPoint(j)[0] << ", ";
            outf << actinVec[i].getCLinksPoint(j)[1] << ", ";
            outf << actinVec[i].getCLinksPoint(j)[2] << ", ";
            outf << actinVec[i].getCLinksPoint(j)[3] << std::endl;
        }
    }
    outf << "nActin: " << nActin << "\n";
    outf << "# actinID, MonoID, birthtime (s), xStart (m), yStart (m), ";
    outf << "xEnd (m), yEnd (m)" << std::endl;

    int monoid { 0 };
    for (Actin actinobj : actinVec)
    {
        for (int i = 0; i < actinobj.getNumMonomers(); ++i)
        {
            std::array<double,4> monoStartEnd = actinobj.findMonoStartEndCoord(i);
            outf << actinobj.getID() << ", " << monoid << ", ";
            outf << actinobj.getMonoBirthTime(i) << ", " << monoStartEnd[0];
            outf << ", " << monoStartEnd[1] << ", ";
            outf << monoStartEnd[2] << ", ";
            outf << monoStartEnd[3] << std::endl;

            ++monoid;
        }
    }
    outf << "#--------------------------------------------------" << std::endl;
}
