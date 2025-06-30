/*
 *  branching.cpp
 *
 *
 *  Functions that create branches and dislocate branches
 *
 *  James Bradford
 *  University of Sheffield
 *  May 2017
 */

#include "configHeader.h"
#include "RNG.h"
#include "Actin.h"
#include "ProteinRegion.h"
#include "GactinGrid.h"
#include "ArpGrid.h"
#include "branching.h"
#include "polymerisation.h"
#include <iterator>
extern RNG rng;


void deBranch(std::vector<Actin> &actinVec, const double p_detach,
              int &nActin, bool arpPool, ArpGrid &arpGrid,
              const bool steric, StericGrid &stericGrid)

{
    // Dislocation of branches.
    std::vector<int> dissociateIDs;


    std::uniform_real_distribution<double> probDist(0, 1);
    for (int i = 0; i < nActin; ++i)
    {
        if (actinVec[i].getParentID() == -1)
            continue;

        int threadID = omp_get_thread_num();
        // The actin filament is a branched filament
        double randdetach { probDist(rng.m_PRNGs[threadID]) };
        if (randdetach < p_detach)
        {
            // Dislocate the particular branch
            int parentid = actinVec[i].getParentID();
            int branchMonoID = actinVec[i].getMotherMonoID();

            std::array<double,3> branchCoords = actinVec[i].getPointedEnd();
            std::array<double,2> branchCoords2D;
            branchCoords2D[0] = branchCoords[0];
            branchCoords2D[1] = branchCoords[1];

            actinVec[i].detach(actinVec, nActin);

            // return arp 2/3 to correct cell
            if (arpPool)
                arpGrid.increment(branchCoords2D[0], branchCoords2D[1]);


            if (actinVec[i].getNumMonomers() < Actin::s_seedSize)
            {
                // it should dissociate
                dissociateIDs.push_back(actinVec[i].getID());
            }

            // recalc distance to leading branch, although this might not change!
            if (branchMonoID == actinVec[parentid].getNumMonomers()-1-actinVec[parentid].getDistToLeadBR())
            {
                int olddist = actinVec[parentid].getDistToLeadBR();
                actinVec[parentid].recalcDistToBR2(actinVec);

                assert(olddist !=  actinVec[parentid].getDistToLeadBR());
            }
            actinVec[parentid].recalcMonoCompat(actinVec, branchMonoID);
        }
    }
    dissociationRout(dissociateIDs, actinVec, nActin, steric, stericGrid);
}


void deBranch2(Actin &actin, std::vector<Actin> &actinVec, int &nActin, bool arpPool,
               ArpGrid &arpGrid, const bool steric, StericGrid &stericGrid)

{
    // as above but it will definitely debranch
    // called when necessary, in severing

    std::vector<int> dissociateIDs;



    // Dislocate the particular branch
    int parentid = actin.getParentID();
    int branchMonoID = actin.getMotherMonoID();

    std::array<double,3> branchCoords = actin.getPointedEnd();
    std::array<double,2> branchCoords2D;
    branchCoords2D[0] = branchCoords[0];
    branchCoords2D[1] = branchCoords[1];

    actin.detach(actinVec, nActin);

    // return arp 2/3 to correct cell
    if (arpPool)
        arpGrid.increment(branchCoords2D[0], branchCoords2D[1]);


    if (actin.getNumMonomers() < Actin::s_seedSize)
    {
        // it should dissociate
        dissociateIDs.push_back(actin.getID());
    }

    // recalc distance to leading branch, although this might not change!
    if (branchMonoID == actinVec[parentid].getNumMonomers()-1-actinVec[parentid].getDistToLeadBR())
    {
        int olddist = actinVec[parentid].getDistToLeadBR();
        actinVec[parentid].recalcDistToBR2(actinVec);

        assert(olddist !=  actinVec[parentid].getDistToLeadBR());
    }

    actinVec[parentid].recalcMonoCompat(actinVec, branchMonoID);

    dissociationRout(dissociateIDs, actinVec, nActin, steric, stericGrid);

}


void branchingRoutine(std::vector<Actin> &actinVec, const double k_branch,
                      const double tStep,
                      const double arpConc, const double gActinConc,
                      int &nActin,
                      const int nActinLimit, const double currTime,
                      const std::vector<ExcZone> &excZones,
                      const std::vector<Membrane> &membranes,
                      const std::vector<MembraneWall> &memWalls,
                      Cortex &cortex,
                      const bool steric, StericGrid &stericGrid,
                      const bool tether,
                      std::normal_distribution<double> tetherDistri)
{
    /*
    Simplest branching subroutine, no depletion of arp2/3 or gActin, not confined
    to regions

    */

    // Now we need to loop over all filaments and then loop over all monomers
    double p_branch_base { k_branch * arpConc * gActinConc * gActinConc * tStep };

    const int imax = nActin;
    for (int i = 0; i < imax; ++i)
    {
        if (nActin >= nActinLimit)
        {
            // Reached the limit, no more branching now
            return;
        }

        if (actinVec[i].getNumMonomers() <= Actin::s_branchSpacing)
        {
            // Not long enough to branch, ignore
            continue;
        }
        // Calculate the probability of branching per monomer

        // Rather than loop over all the monomers we can course grain this to
        // calc the probability of a single branch being created PER FILAMENT
        // Then if the dice roll passes we randomly choose a monomer

        int numAvailMono = actinVec[i].getAvailMonoVec().size();
        double p_branch = p_branch_base * numAvailMono;

        assert(p_branch < 0.1); // statistics check if this is flagging up then dt needs to be decreased

        double randBranch { rng.m_probDist(rng.m_mersenne)};

        if (randBranch < p_branch)
        {
            // So A MONOMER is going to branch, but we need to select which one
            //std::cout << "Branching line 218" << std::endl;

            std::uniform_int_distribution<int> availMonos(0,numAvailMono-1);
            int j = actinVec[i].getAvailMonoID(availMonos(rng.m_mersenne));
            assert(actinVec[i].getMonoCompat(j));
            assert(j >= Actin::s_branchSpacing);

            // Ok this monomer is going to branch
            // We need to add the branch and adjust the compatibility vector
            // of actinVec[i] to exclude any monomers within a branch spacing
            // from this new branch

            // We have to check excluded volume I'm afraid

            /*
            So
            1. Create the branch as a tmp actin object
            2. Check excluded volume, the new 2 monomer length branch against
               everything else! (but not the parent or itself)
            3a. If this fails then end the subroutine, the tmp object gets
                deleted. Will have to do some cleanup.
            3b If this passes then adjust the compatibility vector of the
               mother accordingly, increment nActin and push_back to the main
               actin vector

            */

            Actin actinTMP = Actin(nActin, currTime, actinVec[i], j);

            if (steric)
            {
                for (int k = 0; k < actinTMP.getNumSubs(); ++k)
                {
                    // update the steric grid for the test branch
                    stericGrid.checkSubandUpdate(actinTMP, k);
                }
            }

            if (!actinTMP.checkBranchExcVolGrid(actinVec, excZones, membranes,
                                                memWalls, cortex, steric,
                                                stericGrid))
            {
                // It passed the check
                assert(actinTMP.getMotherMonoID() >= Actin::s_maxSpacing);
                actinTMP.setChangedBoolTrue();
                actinVec.push_back(actinTMP);
                actinVec[i].incrementDaughternum();
                ++nActin;

                // Need to adjust compatibility vector of mother
                actinVec[i].adjustBranchCompat(j);
                actinVec[i].recalcDistToBR(actinVec, nActin);

                if (tether)
                {
                    int predetBarb = round(tetherDistri(rng.m_mersenne)/Actin::s_monomerLength);
                    actinVec[nActin-1].setPreDetTethDistB(predetBarb);
                    int tetLength = actinVec[nActin-1].getNumMonomers();

                    actinVec[nActin-1].setTetherDistB(tetLength);
                }

            }
            else
            {
                // Clean up
                if (steric)
                {
                    for (int k = 0; k < actinTMP.getNumSubs(); ++k)
                    {
                        // reset cells
                        stericGrid.resetCells(actinTMP, k);
                    }
                }
                actinVec[i].purgeBranchFromDaughterVec(nActin);
                Actin::s_total_subunits -= 3;
                // Should we try a different branch?
            }
        }
    }
}




void branchRegionsRoutine(std::vector<Actin> &actinVec, const double k_branch,
                          const double tStep,
                          const double gActinConc,
                          int &nActin,
                          const int nActinLimit, const double currTime,
                          std::vector<ProteinRegion> &branchingRegions,
                          const std::vector<ExcZone> &excZones,
                          const std::vector<Membrane> &membranes,
                          const std::vector<MembraneWall> &memWalls,
                          Cortex &cortex,
                          const bool steric, StericGrid &stericGrid)
{
    if (nActin >= nActinLimit)
    {
        // Reached the limit, no more branching now
        return;
    }

    // Now we need to loop over all filaments and then loop over all monomers

    const int imax = nActin;
    const int nRegions = branchingRegions.size();

    for (int i = 0; i < imax; ++i)
    {
        // Loop over our regions

        if(actinVec[i].getNumMonomers() <= Actin::s_branchSpacing)
        {
            // Not long enough to branch, ignore
            continue;
        }

        for (int j = 0; j < nRegions; ++j)
        {
            double p_branch { k_branch * branchingRegions[j].getArpConc() * gActinConc * gActinConc * tStep };

            std::vector<int> availMonoVec = actinVec[i].getAvailMonoVec();
            std::vector<int> availAndInReg;

            // Need to record the available monomers that are inside the region,
            // no way around this, brute force, but parallise to help!

            #pragma omp parallel
            {
                std::vector<int> availAndInReg_priv;
                #pragma omp for nowait
                for (unsigned int k = 0; k < availMonoVec.size(); ++k)
                {
                    std::array<double,2> monoPoint = actinVec[i].findMonoCoord(availMonoVec[k]);
                    if (branchingRegions[j].checkPointwithin(monoPoint))
                    {
                        availAndInReg_priv.push_back(availMonoVec[k]);
                    }

                }
                #pragma omp critical
                availAndInReg.insert(availAndInReg.end(), availAndInReg_priv.begin(), availAndInReg_priv.end());
            }
            std::sort(availAndInReg.begin(), availAndInReg.end());
            // sort the vector to maintain determinism if same rng seed!

            p_branch *= availAndInReg.size();

            assert(p_branch < 0.1); // statistics check if this is flagging up then dt needs to be decreased

            double randBranch { rng.m_probDist(rng.m_mersenne)};

            if (randBranch < p_branch)
            {
                // Ok we know A MONOMER is going to branch but we have to choose
                std::uniform_int_distribution<int> availMonos(0,availAndInReg.size()-1);
                int k = availAndInReg[availMonos(rng.m_mersenne)];

                // Ok this monomer is going to branch
                // We need to add the branch and adjust the compatibility vector
                // of actinVec[i] to exclude any monomers within a branch spacing
                // from this new branch

                // We have to check excluded volume I'm afraid

                /*
                So
                1. Create the branch as a tmp actin object,
                2. Check excluded volume, the new 2 monomer length branch against
                   everything else! (but not the parent or itself)
                3a. If this fails then end the subroutine, the tmp object gets
                    deleted. Will have to do some cleanup.
                3b If this passes then adjust the compatibility vector of the
                   mother accordingly, increment nActin and push_back to the main
                   actin vector

                */

                Actin actinTMP = Actin(nActin, currTime, actinVec[i], k);
                if (steric)
                {
                    for (int k = 0; k < actinTMP.getNumSubs(); ++k)
                    {
                        // update the steric grid for the test branch
                        stericGrid.checkSubandUpdate(actinTMP, k);
                    }
                }

                if (!actinTMP.checkBranchExcVolGrid(actinVec, excZones, membranes,
                                                    memWalls, cortex, steric,
                                                    stericGrid))
                {
                    // It passed the check
                    assert(actinTMP.getMotherMonoID() >= Actin::s_branchSpacing);
                    actinTMP.setChangedBoolTrue();

                    actinVec.push_back(actinTMP);
                    actinVec[i].incrementDaughternum();
                    ++nActin;

                    // Need to adjust compatibility vector of mother
                    actinVec[i].adjustBranchCompat(k);
                    actinVec[i].recalcDistToBR(actinVec, nActin);
                }
                else
                {
                    // Clean up
                    if (steric)
                    {
                        for (int k = 0; k < actinTMP.getNumSubs(); ++k)
                        {
                            // reset cells
                            stericGrid.resetCells(actinTMP, k);
                        }
                    }
                    actinVec[i].purgeBranchFromDaughterVec(nActin);
                    Actin::s_total_subunits -= 3;
                }
            }
        }
    }
}


void branchArpGridRoutine(std::vector<Actin> &actinVec, const double k_branch,
                          const double tStep,
                          const double gActinConc,
                          int &nActin,
                          const int nActinLimit, const double currTime,
                          ArpGrid &arpGrid,
                          const std::vector<ExcZone> &excZones,
                          const std::vector<Membrane> &membranes,
                          const std::vector<MembraneWall> &memWalls,
                          Cortex &cortex,
                          const bool steric, StericGrid &stericGrid)
{
    if (nActin >= nActinLimit)
    {
        // Reached the limit, no more branching now
        return;
    }

    // Now we need to loop over all filaments and then loop over all monomers

    const int imax = nActin;

    for (int i = 0; i < imax; ++i)
    {
        // Loop over our regions

        if(actinVec[i].getNumMonomers() <= Actin::s_branchSpacing)
        {
            // Not long enough to branch, ignore
            continue;
        }


        std::vector<int> availMonoVec = actinVec[i].getAvailMonoVec();
        std::vector<int> availAndInReg;

        // Need to record the available monomers that are inside the region,
        // no way around this, brute force, but parallise to help!

        // following looks grim but its a 3d list, for each cell [x][y] there's a vector which gives all the available monomers in it
        int NumXCells = arpGrid.getNumCellsX();
        int NumYCells = arpGrid.getNumCellsY();
        std::vector< std::vector< std::vector<int> > > cellAndMonomerRec;
        cellAndMonomerRec.resize(NumXCells);
        for (int x = 0; x < NumXCells; ++x)
        {
            cellAndMonomerRec[x].resize(NumYCells);
        }

        #pragma omp parallel
        {
            std::vector< std::vector< std::vector<int> > > cellAndMonomerRec_priv;
            cellAndMonomerRec_priv.resize(NumXCells);
            for (int x = 0; x < NumXCells; ++x)
            {
                cellAndMonomerRec_priv[x].resize(NumYCells);
            }

            #pragma omp for nowait
            for (unsigned int k = 0; k < availMonoVec.size(); ++k)
            {
                std::array<double,2> monoPoint = actinVec[i].findMonoCoord(availMonoVec[k]);
                std::array<int,2> cellID = arpGrid.findGrid(monoPoint[0], monoPoint[1]);
                std::array<int,2> outside = {-1, -1};
                if (cellID != outside)
                {
                    cellAndMonomerRec_priv[cellID[0]][cellID[1]].push_back(availMonoVec[k]);

                }

            }

            #pragma omp critical
            for (int x = 0; x < NumXCells; ++x)
            {
                for (int y = 0; y < NumYCells; ++y)
                {
                    cellAndMonomerRec[x][y].insert(cellAndMonomerRec[x][y].end(), cellAndMonomerRec_priv[x][y].begin(), cellAndMonomerRec_priv[x][y].end());

                }
            }
        }

        for (int x = 0; x < NumXCells; ++x)
        {
            for (int y = 0; y < NumYCells; ++y)
            {
                if (cellAndMonomerRec[x][y].size() == 0)
                {
                    // No available monomers in this cell, skip
                    continue;
                }

                std::sort(cellAndMonomerRec[x][y].begin(), cellAndMonomerRec[x][y].end());
                double p_branch { k_branch * arpGrid.getConcentration(x,y) * gActinConc * gActinConc * tStep };
                p_branch *= cellAndMonomerRec[x][y].size();
                assert(p_branch < 0.1);
                double randBranch { rng.m_probDist(rng.m_mersenne)};

                if (randBranch < p_branch)
                {
                    // Ok we know A MONOMER is going to branch but we have to choose
                    std::uniform_int_distribution<int> availMonos(0,cellAndMonomerRec[x][y].size()-1);
                    int k = cellAndMonomerRec[x][y][availMonos(rng.m_mersenne)];

                    // Ok this monomer is going to branch
                    // We need to add the branch and adjust the compatibility vector
                    // of actinVec[i] to exclude any monomers within a branch spacing
                    // from this new branch

                    // We have to check excluded volume I'm afraid

                    /*
                    So
                    1. Create the branch as a tmp actin object,
                    2. Check excluded volume, the new 2 monomer length branch against
                       everything else! (but not the parent or itself)
                    3a. If this fails then end the subroutine, the tmp object gets
                        deleted. Will have to do some cleanup.
                    3b If this passes then adjust the compatibility vector of the
                       mother accordingly, increment nActin and push_back to the main
                       actin vector

                    */

                    Actin actinTMP = Actin(nActin, currTime, actinVec[i], k);

                    if (steric)
                    {
                        for (int k = 0; k < actinTMP.getNumSubs(); ++k)
                        {
                            // update the steric grid for the test branch
                            stericGrid.checkSubandUpdate(actinTMP, k);
                        }
                    }

                    if (!actinTMP.checkBranchExcVolGrid(actinVec, excZones,
                                                        membranes, memWalls,
                                                        cortex, steric,
                                                        stericGrid))
                    {
                        // It passed the check
                        assert(actinTMP.getMotherMonoID() >= Actin::s_branchSpacing);
                        actinTMP.setChangedBoolTrue();
                        arpGrid.decrement(x,y);

                        actinVec.push_back(actinTMP);
                        actinVec[i].incrementDaughternum();
                        ++nActin;

                        // Need to adjust compatibility vector of mother
                        actinVec[i].adjustBranchCompat(k);
                        actinVec[i].recalcDistToBR(actinVec, nActin);
                    }
                    else
                    {
                        // Clean up
                        if (steric)
                        {
                            for (int k = 0; k < actinTMP.getNumSubs(); ++k)
                            {
                                // reset cells
                                stericGrid.resetCells(actinTMP, k);
                            }
                        }
                        actinVec[i].purgeBranchFromDaughterVec(nActin);
                        Actin::s_total_subunits -= 3;
                    }
                }
            }
        }
    }
}

void branchArpAndActinGridRoutine(std::vector<Actin> &actinVec, const double k_branch,
                                  const double tStep,
                                  int &nActin,
                                  const int nActinLimit, const double currTime,
                                  ArpGrid &arpGrid,GactinGrid &gActinGrid,
                                  const std::vector<ExcZone> &excZones,
                                  const std::vector<Membrane> &membranes,
                                  const std::vector<MembraneWall> &memWalls,
                                  Cortex &cortex,
                                  const bool steric, StericGrid &stericGrid)
{
    if (nActin >= nActinLimit)
    {
        // Reached the limit, no more branching now
        return;
    }

    // Now we need to loop over all filaments and then loop over all monomers

    const int imax = nActin;

    for (int i = 0; i < imax; ++i)
    {
        // Loop over our regions

        if(actinVec[i].getNumMonomers() <= Actin::s_branchSpacing)
        {
            // Not long enough to branch, ignore
            continue;
        }


        std::vector<int> availMonoVec = actinVec[i].getAvailMonoVec();
        std::vector<int> availAndInReg;

        // Need to record the available monomers that are inside the region,
        // no way around this, brute force, but parallise to help!

        // following looks grim but its a 3d list, for each cell [x][y] there's a vector which gives all the available monomers in it
        int NumXCells = arpGrid.getNumCellsX();
        int NumYCells = arpGrid.getNumCellsY();
        std::vector< std::vector< std::vector<int> > > cellAndMonomerRec;
        cellAndMonomerRec.resize(NumXCells);

        int ActinXCells = gActinGrid.getNumCellsX();
        int ActinYCells = gActinGrid.getNumCellsY();
        assert(ActinXCells == NumXCells && ActinYCells==NumYCells);

        for (int x = 0; x < NumXCells; ++x)
        {
            cellAndMonomerRec[x].resize(NumYCells);
        }

        #pragma omp parallel
        {
            std::vector< std::vector< std::vector<int> > > cellAndMonomerRec_priv;
            cellAndMonomerRec_priv.resize(NumXCells);
            for (int x = 0; x < NumXCells; ++x)
            {
                cellAndMonomerRec_priv[x].resize(NumYCells);
            }

            #pragma omp for nowait
            for (unsigned int k = 0; k < availMonoVec.size(); ++k)
            {
                std::array<double,2> monoPoint = actinVec[i].findMonoCoord(availMonoVec[k]);
                std::array<int,2> cellID = arpGrid.findGrid(monoPoint[0], monoPoint[1]);
                std::array<int,2> outside = {-1, -1};
                if (cellID != outside)
                {
                    cellAndMonomerRec_priv[cellID[0]][cellID[1]].push_back(availMonoVec[k]);

                }

            }

            #pragma omp critical
            for (int x = 0; x < NumXCells; ++x)
            {
                for (int y = 0; y < NumYCells; ++y)
                {
                    cellAndMonomerRec[x][y].insert(cellAndMonomerRec[x][y].end(), cellAndMonomerRec_priv[x][y].begin(), cellAndMonomerRec_priv[x][y].end());

                }
            }
        }

        // sort the vector to maintain determinism if same rng seed!

        for (int x = 0; x < NumXCells; ++x)
        {
            for (int y = 0; y < NumYCells; ++y)
            {
                if (cellAndMonomerRec[x][y].size() == 0)
                {
                    // No available monomers in this cell, skip
                    continue;
                }

                std::sort(cellAndMonomerRec[x][y].begin(), cellAndMonomerRec[x][y].end());
                double p_branch { k_branch * arpGrid.getConcentration(x,y) * gActinGrid.getConcentration(x,y) * gActinGrid.getConcentration(x,y) * tStep };
                p_branch *= cellAndMonomerRec[x][y].size();
                assert(p_branch < 0.1);

                double randBranch { rng.m_probDist(rng.m_mersenne)};

                if (randBranch < p_branch)
                {
                    // Ok we know A MONOMER is going to branch but we have to choose
                    std::uniform_int_distribution<int> availMonos(0,cellAndMonomerRec[x][y].size()-1);
                    int k = cellAndMonomerRec[x][y][availMonos(rng.m_mersenne)];

                    // Ok this monomer is going to branch
                    // We need to add the branch and adjust the compatibility vector
                    // of actinVec[i] to exclude any monomers within a branch spacing
                    // from this new branch

                    // We have to check excluded volume I'm afraid

                    /*
                    So
                    1. Create the branch as a tmp actin object,
                    2. Check excluded volume, the new 2 monomer length branch against
                       everything else! (but not the parent or itself)
                    3a. If this fails then end the subroutine, the tmp object gets
                        deleted. Will have to do some cleanup.
                    3b If this passes then adjust the compatibility vector of the
                       mother accordingly, increment nActin and push_back to the main
                       actin vector

                    */

                    Actin actinTMP = Actin(nActin, currTime, actinVec[i], k);
                    if (steric)
                    {
                        for (int k = 0; k < actinTMP.getNumSubs(); ++k)
                        {
                            // update the steric grid for the test branch
                            stericGrid.checkSubandUpdate(actinTMP, k);
                        }
                    }

                    if ((gActinGrid.getNumber(x,y) > 1) &&
                        (!actinTMP.checkBranchExcVolGrid(actinVec, excZones,
                                                        membranes, memWalls,
                                                        cortex, steric,
                                                        stericGrid)))
                    {
                        // It passed the check
                        assert(actinTMP.getMotherMonoID() >= Actin::s_branchSpacing);
                        gActinGrid.branchingCell(x,y);
                        actinTMP.setChangedBoolTrue();
                        arpGrid.decrement(x,y);

                        actinVec.push_back(actinTMP);
                        actinVec[i].incrementDaughternum();
                        ++nActin;

                        // Need to adjust compatibility vector of mother
                        actinVec[i].adjustBranchCompat(k);
                        actinVec[i].recalcDistToBR(actinVec, nActin);
                    }
                    else
                    {
                        // Clean up
                        if (steric)
                        {
                            for (int k = 0; k < actinTMP.getNumSubs(); ++k)
                            {
                                // reset cells
                                stericGrid.resetCells(actinTMP, k);
                            }
                        }
                        actinVec[i].purgeBranchFromDaughterVec(nActin);
                        Actin::s_total_subunits -= 3;
                    }
                }
            }
        }
    }
}
