/*
 *  nucleation.cpp
 *
 *
 *  Functions that nucleate filaments, both to initialise the system and during
 *  the simulation.
 *
 *  James Bradford
 *  University of Sheffield
 *  May 2017
 */

#include "configHeader.h"
#include "RNG.h"
#include "Actin.h"
#include "ExcZone.h"
#include "GactinGrid.h"
#include "nucleation.h"

extern RNG rng;

void initNuc(std::vector<Actin> &actinvec, const int nActin,
             const std::vector<ProteinRegion> &nucRegions,
             const std::vector<ExcZone> &excZones,
             const std::vector<Membrane> &membranes,
             const std::vector<MembraneWall> &memWalls,
             Cortex &cortex,
             const bool steric, StericGrid &stericGrid,
             const bool tether,
             std::normal_distribution<double> tetherDistri,
             bool nucDir, double mean, double stdev)
{
    /*
     Initial nucleation, setups up and starting filaments in the nucRegion
    */

    if (nucRegions.size() == 0 && nActin != 0)
    {
        std::cout << "Can't have initial nucleation if no nuc regions" << std::endl;
        assert(nucRegions.size() > 0);
    }

    std::normal_distribution<double> angleDist;
    if (nucDir)
    {
        angleDist = std::normal_distribution<double> (mean, stdev);
    }

    for (int i = 0; i < nActin; ++i)
    {
        int j = ProteinRegion::s_chooseRegion(nucRegions, nucRegions.size());
        if (nucDir)
        {
            if (nucRegions[j].getCircBool() || nucRegions[j].getRingBool())
            {
                actinvec.push_back(Actin(i, nucRegions[j].getRDist(),
                                   nucRegions[j].getThetaDist(), j, nucRegions, angleDist));
            }
            else
            {
                actinvec.push_back(Actin(i, nucRegions[j].getWidthDist(),
                                   nucRegions[j].getHeightDist(), j, nucRegions, angleDist));
            }
        }
        else
        {
            if (nucRegions[j].getCircBool() || nucRegions[j].getRingBool())
            {
                actinvec.push_back(Actin(i, nucRegions[j].getRDist(), nucRegions[j].getThetaDist(), j, nucRegions));
            }
            else
            {
                actinvec.push_back(Actin(i, nucRegions[j].getWidthDist(), nucRegions[j].getHeightDist(), j, nucRegions));
            }
        }

        if (steric)
        {
            for (int k = 0; k < actinvec[i].getNumSubs(); ++k)
            {
                // update the steric grid for the test seed
                stericGrid.checkSubandUpdate(actinvec[i], k);
            }
        }

        const int cutOffAttempts = 10; // maximum number of times to place each one
        int nAttempts = 1;
        while (actinvec[i].check_nucleation_excVol_Grid(actinvec, excZones,
                                                        membranes, memWalls,
                                                        cortex, steric,
                                                        stericGrid))
        {
            nAttempts++;
            if (nAttempts > cutOffAttempts)
            {
                std::cout << "Nucleation region too dense, either try fewer init";
                std::cout << " filaments or increase size of area. Aborting" << std::endl;
                exit(0);
            }


            actinvec[i].regenerate_pos(nucRegions[j]);

            if (steric)
            {
                for (int k = 0; k < actinvec[i].getNumSubs(); ++k)
                {
                    // update the steric grid for the test seed
                    stericGrid.resetSubandUpdate(actinvec[i], k);
                }
            }
        }
        if (tether)
        {
            actinvec[i].addTetherPointNuc();
            int predetBarb = round(tetherDistri(rng.m_mersenne)/Actin::s_monomerLength);
            int predetPoint = round(tetherDistri(rng.m_mersenne)/Actin::s_monomerLength);

            actinvec[i].setPreDetTethDistB(predetBarb);
            actinvec[i].setPreDetTethDistP(predetPoint);

        }

    }
}

void nuc(std::vector<Actin> &actinvec, int &nActin, const double k_NucDens,
         const double dt,
         const double gActinConc, const double Time,
         const int nActinLimit,
         const std::vector<ProteinRegion> &nucRegions,
         const std::vector<ExcZone> &excZones,
         const std::vector<Membrane> &membranes,
         const std::vector<MembraneWall> &memWalls,
         Cortex &cortex,
         const bool steric, StericGrid &stericGrid,
         const bool tether, std::normal_distribution<double> tetherDistri,
         bool nucDir, double mean, double stdev)
{
    std::normal_distribution<double> angleDist = std::normal_distribution<double> (mean, stdev);

    // Nucleation rate should be a density,
    if (nActin >= nActinLimit)
    {
        // Reached the limit, no more nucleation now
        return;
    }
    // loop over our regions

    double gActindep = gActinConc*gActinConc*gActinConc; // cubic

    for (unsigned int i = 0; i < nucRegions.size(); ++i)
    {
        // work out nucleation prob for the particular regionID
        double nucProb = k_NucDens * dt * nucRegions[i].getNucCoeff() * gActindep
                         * (nucRegions[i].getArea() * 1E12); // units of kNuc are per um^2 per s
        assert(nucProb < 0.1); // for good statistics

        double randNuc { rng.m_probDist(rng.m_mersenne) };
        if (randNuc < nucProb)
        {
            // Nucleate
            Actin actin_tmp;
            if (nucDir)
            {
                if (nucRegions[i].getCircBool() || nucRegions[i].getRingBool())
                {
                    actin_tmp = Actin(nActin, nucRegions[i].getRDist(), nucRegions[i].getThetaDist(), i, nucRegions, angleDist, Time);
                }
                else
                {
                    actin_tmp = Actin(nActin, nucRegions[i].getWidthDist(), nucRegions[i].getHeightDist(), i, nucRegions, angleDist, Time);
                }
            }
            else
            {
                if (nucRegions[i].getCircBool() || nucRegions[i].getRingBool())
                {
                    actin_tmp = Actin(nActin, nucRegions[i].getRDist(), nucRegions[i].getThetaDist(), i, nucRegions, Time);
                }
                else
                {
                    actin_tmp = Actin(nActin, nucRegions[i].getWidthDist(), nucRegions[i].getHeightDist(), i, nucRegions, Time);
                }
            }

            if (steric)
            {
                for (int j = 0; j < actin_tmp.getNumSubs(); ++j)
                {
                    // update the steric grid for the test seed
                    stericGrid.checkSubandUpdate(actin_tmp, j);
                }
            }

            if (actin_tmp.check_nucleation_excVol_Grid(actinvec, excZones,
                                                       membranes, memWalls,
                                                       cortex, steric,
                                                       stericGrid))
            {
                // Nucleation failed due to steric hindrance
                Actin::s_total_subunits -= 3;
                if (steric)
                {
                    for (int j = 0; j < actin_tmp.getNumSubs(); ++j)
                    {
                        // Remove from steric Grid
                        stericGrid.resetCells(actin_tmp, j);
                    }
                }
            }
            else
            {
                // Success
                if (tether)
                {

                    actin_tmp.addTetherPointNuc();
                    int predetBarb = round(tetherDistri(rng.m_mersenne)/Actin::s_monomerLength);
                    int predetPoint = round(tetherDistri(rng.m_mersenne)/Actin::s_monomerLength);
                    actin_tmp.setPreDetTethDistB(predetBarb);
                    actin_tmp.setPreDetTethDistP(predetPoint);
                }

                actinvec.push_back(actin_tmp);

                ++nActin;
            }
        }
    }
}

void nucGrid(std::vector<Actin> &actinvec, int &nActin,
             const double k_NucDens, const double dt, const double Time,
             const int nActinLimit,
             const std::vector<ProteinRegion> &nucRegions,
             const std::vector<ExcZone> &excZones,
             const std::vector<Membrane> &membranes,
             const std::vector<MembraneWall> &memWalls,
             Cortex &cortex,
             const bool steric, StericGrid &stericGrid,
             const bool tether, std::normal_distribution<double> tetherDistri,
             GactinGrid &gActinGrid, bool nucDir, double mean, double stdev)
{
    std::normal_distribution<double> angleDist = std::normal_distribution<double> (mean, stdev);
    // Nucleation rate should be a density,
    if (nActin >= nActinLimit)
    {
        // Reached the limit, no more nucleation now
        return;
    }
    // loop over our regions

    for (unsigned int i = 0; i < nucRegions.size(); ++i)
    {
        // work out nucleation prob for the particular regionID
        std::array<double,2> nucCentre = nucRegions[i].getCentre();
        std::array<int, 2> gridCoords;
        gridCoords = gActinGrid.findGrid(nucCentre[0],
                                         nucCentre[1]);
        // Check to see if nuc site is outside the grid
        if (gridCoords[0] == -1)
        {
            std::cout << "Error: Nucleation region is outside the gActin grid!" << std::endl;
            exit(1);
        }

        double gActindep = gActinGrid.getConcentration(gridCoords[0],gridCoords[1])*gActinGrid.getConcentration(gridCoords[0],gridCoords[1])*gActinGrid.getConcentration(gridCoords[0],gridCoords[1]); // cubic

        double nucProb = k_NucDens * dt * nucRegions[i].getNucCoeff() * gActindep
                         * (nucRegions[i].getArea() * 1E12); // units of pNuc are per um^2
        assert(nucProb < 0.1); // for good statistics

        double randNuc { rng.m_probDist(rng.m_mersenne) };
        if (randNuc < nucProb)
        {
            // Nucleate
            Actin actin_tmp;
            if (nucDir)
            {
                if (nucRegions[i].getCircBool() || nucRegions[i].getRingBool())
                {
                    actin_tmp = Actin(nActin, nucRegions[i].getRDist(), nucRegions[i].getThetaDist(), i, nucRegions, angleDist, Time);
                }
                else
                {
                    actin_tmp = Actin(nActin, nucRegions[i].getWidthDist(), nucRegions[i].getHeightDist(), i, nucRegions, angleDist, Time);
                }
            }
            else
            {
                if (nucRegions[i].getCircBool() || nucRegions[i].getRingBool())
                {
                    actin_tmp = Actin(nActin, nucRegions[i].getRDist(), nucRegions[i].getThetaDist(), i, nucRegions, Time);
                }
                else
                {
                    actin_tmp = Actin(nActin, nucRegions[i].getWidthDist(), nucRegions[i].getHeightDist(), i, nucRegions, Time);
                }
            }

            if (steric)
            {
                for (int j = 0; j < actin_tmp.getNumSubs(); ++j)
                {
                    // update the steric grid for the test seed
                    stericGrid.checkSubandUpdate(actin_tmp, j);
                }
            }
            gridCoords = gActinGrid.findGrid(actin_tmp.getPoints()[2][0],
                                             actin_tmp.getPoints()[2][1]);

            if ((actin_tmp.check_nucleation_excVol_Grid(actinvec, excZones,
                                                       membranes, memWalls,
                                                      cortex, steric,
                                                      stericGrid)) ||
                (gActinGrid.getNumber(gridCoords[0],gridCoords[1]) < 3))
            {
                // Nucleation failed either sterically or not enough monomers in the space
                Actin::s_total_subunits -= 3;
                if (steric)
                {
                    for (int j = 0; j < actin_tmp.getNumSubs(); ++j)
                    {
                        // Remove from steric Grid
                        stericGrid.resetCells(actin_tmp, j);
                    }
                }
            }
            else
            {
                if (tether)
                {

                    actin_tmp.addTetherPointNuc();
                    int predetBarb = round(tetherDistri(rng.m_mersenne)/Actin::s_monomerLength);
                    int predetPoint = round(tetherDistri(rng.m_mersenne)/Actin::s_monomerLength);

                    actin_tmp.setPreDetTethDistB(predetBarb);
                    actin_tmp.setPreDetTethDistP(predetPoint);

                }


                gActinGrid.nucleation(actin_tmp.getPoints()[2][0],actin_tmp.getPoints()[2][1]);
                actinvec.push_back(actin_tmp);
                ++nActin;
            }
        }
    }

}

void setSteric(Actin &filament, const double quasi3dProb)
{
    // Given the bool for quasi-3d is true this function is called upon nucleation
    // quasi-3d prob is the probability to IGNORE stericity, therefore a higher
    // prob means lower constraints

    double randSteric { rng.m_probDist(rng.m_mersenne) };
    if (randSteric < quasi3dProb)
    {
        filament.turnStericOff();
    }
}
