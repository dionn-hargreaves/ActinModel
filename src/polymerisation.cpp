/*
 *  polymerisation.cpp
 *
 *
 *  Functions that polymerise or depolymerise
 *
 *  James Bradford
 *  University of Sheffield
 *  May 2017
 */

#include "configHeader.h"
#include "RNG.h"
#include "Actin.h"
#include "ExcZone.h"
#include "MembraneWall.h"
#include "Membrane.h"
#include "StericGrid.h"
#include "polymerisation.h"
#include "GactinGrid.h"
#include "ArpGrid.h"
#include "globals.h"
#include "Cortex.h"
#include "crosslinking.h"

extern RNG rng;

void polymerisationBarbed(std::vector<Actin> &actinVec, const double gActinConc,
                          const double k_on, const double dt, const double currTime,
                          const std::vector<ExcZone> &excZones, const bool steric,
                          std::vector<MembraneWall> &memWalls,
                          std::vector<Membrane> &membranes,
                          const double temp,
                          const bool tether,
                          std::normal_distribution<double> tetherDistri,
                          std::vector<ProteinRegion> &branchRegions,
                          std::vector<ProteinRegion> &nucRegions,
                          std::vector<ProteinRegion> &capRegions,
                          std::vector<ProteinRegion> &sevRegions,
                          StericGrid &stericGrid,
                          const bool bDynamics, Cortex &cortex,
                          int &nActin, bool arpPool, ArpGrid &arpGrid)
{
    for (Actin &actinobj : actinVec)
    {

        if (!actinobj.getBarbedCapped())
        {
            double randPoly { rng.m_probDist(rng.m_mersenne) };

            double p_poly = k_on * gActinConc * dt;

            bool moveWall = false;
            double d_i = 0;
            if (memWalls.size() > 0)
            {
                d_i = memWalls[0].distToMoveWall_barbed(actinobj);
                if (d_i != 0)
                {
                    double exponent = -(memWalls[0].getForce() * d_i) / (g_Kb*temp);
                    p_poly *= exp(exponent);
                    moveWall = true;
                }

            }

            if (randPoly < p_poly)
            {

                // Polymerise and then check

                actinobj.polymerise_barbed(actinVec, steric, stericGrid,
                                           bDynamics, tether, currTime, tetherDistri,
                                           memWalls, moveWall, d_i,
                                           branchRegions, nucRegions,
                                           capRegions, sevRegions);


                if (steric)
                {
                    stericGrid.updateCellsBarbPoly(actinobj);
                }

                if ( (ExcZone::s_check_polymerisation_barbed(actinobj, excZones)) ||
                     (Membrane::s_checkExVolBarbPoly(actinobj, membranes, stericGrid)) ||
                     (Cortex::s_checkExVolBarbPoly(actinobj, cortex, stericGrid)) ||
                     (actinobj.check_barbed_polymerisation_Grid(actinVec, steric, stericGrid)) )
                {
                    //reverse the pol

                    actinobj.depolymerise_barbed(actinVec, nActin,
                                                 arpPool, arpGrid, steric,
                                                 stericGrid, excZones, memWalls,
                                                 membranes, cortex, bDynamics,
                                                 tether, moveWall, -d_i,
                                                 branchRegions, nucRegions,
                                                 capRegions, sevRegions);
                    // Steric check failed, update back
                    if (steric)
                        stericGrid.updateCellsBarbPoly(actinobj);
                }
            }
        }
    }
}


void polymerisationBarbedGrid(std::vector<Actin> &actinVec, GactinGrid &gActinGrid,
                              const double k_on, const double dt, const double currTime,
                              const std::vector<ExcZone> &excZones, const bool steric,
                              std::vector<MembraneWall> &memWalls,
                              std::vector<Membrane> &membranes,
                              const double temp,
                              const bool tether,
                              std::normal_distribution<double> tetherDistri,
                              std::vector<ProteinRegion> &branchRegions,
                              std::vector<ProteinRegion> &nucRegions,
                              std::vector<ProteinRegion> &capRegions,
                              std::vector<ProteinRegion> &sevRegions,
                              StericGrid &stericGrid,
                              const bool bDynamics, Cortex &cortex,
                              int &nActin, bool arpPool, ArpGrid &arpGrid)
{
    for (Actin &actinobj : actinVec)
    {
        if (!actinobj.getBarbedCapped())
        {

            std::array<int, 2> gridCoords;
            gridCoords = gActinGrid.findGrid(actinobj.getBarbedEnd()[0],
                                             actinobj.getBarbedEnd()[1]);
            // Check to see if barbed end is outside the grid
            if (gridCoords[0] == -1)
                continue;


            double gActinConc = gActinGrid.getConcentration(gridCoords[0],gridCoords[1]);
            double randPoly { rng.m_probDist(rng.m_mersenne) };

            double p_poly = k_on * gActinConc * dt;


            bool moveWall = false;
            double d_i = 0;
            if (memWalls.size() > 0)
            {
                d_i = memWalls[0].distToMoveWall_barbed(actinobj);
                if (d_i != 0)
                {
                    double exponent = -(memWalls[0].getForce() * d_i) / (g_Kb*temp);
                    p_poly *= exp(exponent);
                    moveWall = true;
                }
            }
            if (randPoly < p_poly)
            {
                actinobj.polymerise_barbed(actinVec, steric, stericGrid,
                                           bDynamics, tether, currTime, tetherDistri,
                                           memWalls, moveWall, d_i,
                                           branchRegions, nucRegions,
                                           capRegions, sevRegions);

                if (steric)
                    stericGrid.updateCellsBarbPoly(actinobj);


                if ( (!ExcZone::s_check_polymerisation_barbed(actinobj, excZones)) &&
                     (!actinobj.check_barbed_polymerisation_Grid(actinVec, steric, stericGrid)) )
                {
                    gActinGrid.decrement(actinobj.getBarbedEnd()[0], actinobj.getBarbedEnd()[1]);
                }
                else
                {
                    actinobj.depolymerise_barbed(actinVec, nActin,
                                                 arpPool, arpGrid, steric,
                                                 stericGrid, excZones, memWalls,
                                                 membranes, cortex, bDynamics,
                                                 tether, moveWall, -d_i,
                                                 branchRegions, nucRegions,
                                                 capRegions, sevRegions);
                    // Steric check failed, update back
                    if (steric)
                        stericGrid.updateCellsBarbPoly(actinobj);
                }
            }
        }
    }
}


void polymerisationPointed(std::vector<Actin> &actinVec, const double gActinConc,
                           const double k_on, const double dt, const double currTime,
                           const std::vector<ExcZone> &excZones, const bool steric,
                           std::vector<MembraneWall> &memWalls,
                           std::vector<Membrane> &membranes,
                           const double temp,
                           const bool tether,
                           std::normal_distribution<double> tetherDistri,
                           std::vector<ProteinRegion> &branchRegions,
                           std::vector<ProteinRegion> &nucRegions,
                           std::vector<ProteinRegion> &capRegions,
                           std::vector<ProteinRegion> &sevRegions,
                           StericGrid &stericGrid, const bool bDynamics,
                           Cortex &cortex,
                           int &nActin, bool arpPool, ArpGrid &arpGrid)
{
    for (Actin &actinobj : actinVec)
    {
        if (!actinobj.getPointedCapped())
        {
            double randPoly { rng.m_probDist(rng.m_mersenne) };

            double p_poly = k_on * gActinConc * dt;

            bool moveWall = false;
            double d_i = 0;

            if (memWalls.size() > 0)
            {
                d_i = memWalls[0].distToMoveWall_pointed(actinobj);
                if (d_i != 0)
                {
                    double exponent = -(memWalls[0].getForce() * d_i) / (g_Kb*temp);
                    p_poly *= exp(exponent);
                    moveWall = true;
                }
            }
            if (randPoly < p_poly)
            {
                actinobj.polymerise_pointed(actinVec, steric, stericGrid,
                                            bDynamics, tether, currTime, tetherDistri,
                                            memWalls, moveWall, d_i,
                                            branchRegions, nucRegions,
                                            capRegions, sevRegions);
                if (steric)
                    stericGrid.updateCellsPointPoly(actinobj);

                if ( (ExcZone::s_check_polymerisation_pointed(actinobj, excZones)) ||
                     (Membrane::s_checkExVolPointPoly(actinobj, membranes, stericGrid)) ||
                     (Cortex::s_checkExVolPointPoly(actinobj, cortex, stericGrid)) ||
                     (actinobj.check_pointed_polymerisation_excVol_Grid(actinVec, steric, stericGrid)) )
                {

                    actinobj.depolymerise_pointed(actinVec, nActin,
                                                     arpPool,
                                                     arpGrid,
                                                     steric, stericGrid,
                                                     excZones, memWalls,
                                                     membranes, cortex,
                                                     bDynamics, tether, moveWall,
                                                     -d_i,
                                                     branchRegions, nucRegions,
                                                     capRegions, sevRegions);
                    // Steric check failed, update back
                    if (steric)
                        stericGrid.updateCellsPointPoly(actinobj);
                }
            }
        }
    }
}
void polymerisationPointedGrid(std::vector<Actin> &actinVec, GactinGrid &gActinGrid,
                               const double k_on, const double dt, const double currTime,
                               const std::vector<ExcZone> &excZones, const bool steric,
                               std::vector<MembraneWall> &memWalls,
                               std::vector<Membrane> &membranes,
                               const double temp,
                               const bool tether,
                               std::normal_distribution<double> tetherDistri,
                               std::vector<ProteinRegion> &branchRegions,
                               std::vector<ProteinRegion> &nucRegions,
                               std::vector<ProteinRegion> &capRegions,
                               std::vector<ProteinRegion> &sevRegions,
                               StericGrid &stericGrid, const bool bDynamics,
                               Cortex &cortex,
                               int &nActin, bool arpPool, ArpGrid &arpGrid)
{
    for (Actin &actinobj : actinVec)
    {
        if (!actinobj.getPointedCapped())
        {
            std::array<int, 2> gridCoords;
            gridCoords = gActinGrid.findGrid(actinobj.getPointedEnd()[0],
                                             actinobj.getPointedEnd()[1]);
            // Check to see if barbed end is outside the grid
            if (gridCoords[0] == -1)
                continue;

            double gActinConc = gActinGrid.getConcentration(gridCoords[0],gridCoords[1]);
            double randPoly { rng.m_probDist(rng.m_mersenne) };

            double p_poly = k_on * gActinConc * dt;

            bool moveWall = false;
            double d_i = 0;

            if (memWalls.size() > 0)
            {
                d_i = memWalls[0].distToMoveWall_pointed(actinobj);
                if (d_i != 0)
                {
                    double exponent = -(memWalls[0].getForce() * d_i) / (g_Kb*temp);
                    p_poly *= exp(exponent);
                    moveWall = true;
                }
            }
            if (randPoly < p_poly)
            {
                actinobj.polymerise_pointed(actinVec, steric, stericGrid,
                                            bDynamics, tether, currTime, tetherDistri,
                                            memWalls, moveWall, d_i,
                                            branchRegions, nucRegions,
                                            capRegions, sevRegions);

                if (steric)
                    stericGrid.updateCellsPointPoly(actinobj);

                if ( (!ExcZone::s_check_polymerisation_pointed(actinobj, excZones)) &&
                     (!actinobj.check_pointed_polymerisation_excVol_Grid(actinVec, steric, stericGrid)) )
                {

                    gActinGrid.decrement(actinobj.getPointedEnd()[0], actinobj.getPointedEnd()[1]);
                }
                else
                {
                    actinobj.depolymerise_pointed(actinVec, nActin,
                                                     arpPool,
                                                     arpGrid,
                                                     steric, stericGrid,
                                                     excZones, memWalls,
                                                     membranes, cortex,
                                                     bDynamics, tether, moveWall,
                                                     -d_i,
                                                     branchRegions, nucRegions,
                                                     capRegions, sevRegions);
                    if (steric)
                        stericGrid.updateCellsPointPoly(actinobj);
                }
            }
        }
    }
}

void depolymerisationBarbed(std::vector<Actin> &actinVec, const double k_off,
                            const double dt, const bool tether, int &nActin,
                            const bool dissociate,
                            bool arpPool, ArpGrid &arpGrid, const bool steric,
                            StericGrid &stericGrid,
                            const std::vector<ExcZone> &excZones,
                            std::vector<MembraneWall> &memWalls,
                            std::vector<Membrane> &membranes,
                            Cortex &cortex,
                            const bool bDynamics,
                            std::vector<ProteinRegion> &branchRegions,
                            std::vector<ProteinRegion> &nucRegions,
                            std::vector<ProteinRegion> &capRegions,
                            std::vector<ProteinRegion> &sevRegions)
{
    std::uniform_real_distribution<double> probDist(0, 1);

    std::vector<int> dissociateIDs;

    //#pragma omp parallel for shared(nActin) private(probDist)
    for (int i = 0; i < nActin; ++i)
    {
        int threadID = omp_get_thread_num();
        if (!actinVec[i].getBarbedCapped())
        {
            double randDepoly { probDist(rng.m_PRNGs[threadID]) };
            double p_depoly = k_off * dt;
            if (randDepoly < p_depoly)
            {
                if (dissociate && actinVec[i].getNumMonomers() == Actin::s_seedSize && actinVec[i].getParentID() == -1)
                {
                    /*
                    So this filament is going to break apart into monomers
                    It needs to be deleted
                    This is tricky because the id system we use is based on
                    filament's position in the actinVec. If we delete a fila
                    we alter this so that id != index in the array
                    we would have to subtract 1 from everyone's id in the array
                    It would also affect parent IDs - any daughters of filas
                    that come after filament i in the actinVec need to have their
                    parent id reduced by 1

                    Structure ids are fine

                    Would also have to deal with branches, its unlikely but possible
                    that it would still have a branch attached

                    Be careful because altering actinVec will affect this loop we
                    are in, especially when running in parallel

                    */

                    dissociateIDs.push_back(actinVec[i].getID());
                    continue;
                }
                else if (actinVec[i].getNumMonomers() == Actin::s_seedSize && actinVec[i].getParentID() == -1)
                {
                    // dissociation is turned off, but we dont want to depoly it
                    // just leave it

                    continue;
                }
                else if (dissociate && actinVec[i].getNumMonomers() == Actin::s_branchSeedSize && actinVec[i].getParentID() != -1)
                {
                    // It is a branch and dissociation is turned on, so it is removed

                    dissociateIDs.push_back(actinVec[i].getID());
                    if (arpPool)
                    {
                        // we have an arp grid and so must increment by 1
                        arpGrid.increment(actinVec[i].getBarbedEnd()[0], actinVec[i].getBarbedEnd()[1]);
                    }
                    continue;
                }
                else if (actinVec[i].getNumMonomers() == Actin::s_branchSeedSize && actinVec[i].getParentID() != -1)
                {
                    // It is a branch, we allow it to go down to 2 but branches
                    continue;
                }


                actinVec[i].depolymerise_barbed(actinVec, nActin,
                                                arpPool, arpGrid, steric,
                                                stericGrid, excZones, memWalls,
                                                membranes, cortex, bDynamics,
                                                tether, 0, 0,
                                                branchRegions, nucRegions,
                                                capRegions, sevRegions);

                if (steric)
                    stericGrid.updateCellsBarbPoly(actinVec[i]);

            }
        }
    }

    // After the loop, loop over the dissociateIDs vector
    dissociationRout(dissociateIDs, actinVec, nActin, steric, stericGrid);

}

void depolymerisationBarbedForce(int i, std::vector<Actin> &actinVec,
                                 const bool tether, int &nActin,
                                 bool arpPool, ArpGrid &arpGrid, const bool steric,
                                 StericGrid &stericGrid,
                                 std::vector<ExcZone> &excZones,
                                 std::vector<MembraneWall> &memWalls,
                                 std::vector<Membrane> &membranes,
                                 Cortex &cortex,
                                 const bool bDynamics,
                                 std::vector<ProteinRegion> &branchRegions,
                                 std::vector<ProteinRegion> &nucRegions,
                                 std::vector<ProteinRegion> &capRegions,
                                 std::vector<ProteinRegion> &sevRegions)
{
    /*
     * Definitely will happen with a single filament passed to it
    */
    std::vector<int> dissociateIDs;

    if (actinVec[i].getNumMonomers() == Actin::s_seedSize && actinVec[i].getParentID() == -1)
    {
        /*
        So this filament is going to break apart into monomers
        It needs to be deleted
        This is tricky because the id system we use is based on
        filament's position in the actinVec. If we delete a fila
        we alter this so that id != index in the array
        we would have to subtract 1 from everyone's id in the array
        It would also affect parent IDs - any daughters of filas
        that come after filament i in the actinVec need to have their
        parent id reduced by 1

        Structure ids are fine

        Would also have to deal with branches, its unlikely but possible
        that it would still have a branch attached

        Be careful because altering actinVec will affect this loop we
        are in, especially when running in parallel

        */

        dissociateIDs.push_back(actinVec[i].getID());
        dissociationRout(dissociateIDs, actinVec, nActin, steric, stericGrid);
        return;
    }

    else if (actinVec[i].getNumMonomers() == Actin::s_branchSeedSize && actinVec[i].getParentID() != -1)
    {
        // It is a branch and dissociation is turned on, so it is removed
        dissociateIDs.push_back(actinVec[i].getID());
        if (arpPool)
        {
            // we have an arp grid and so must increment by 1
            arpGrid.increment(actinVec[i].getBarbedEnd()[0], actinVec[i].getBarbedEnd()[1]);
        }
        dissociationRout(dissociateIDs, actinVec, nActin, steric, stericGrid);
        return;
    }
    actinVec[i].depolymerise_barbed(actinVec, nActin,
                                    arpPool, arpGrid, steric,
                                    stericGrid, excZones, memWalls,
                                    membranes, cortex, bDynamics,
                                    tether, 0, 0,
                                    branchRegions, nucRegions,
                                    capRegions, sevRegions);





    if (steric)
        stericGrid.updateCellsBarbPoly(actinVec[i]);

}

void depolymerisationBarbedGrid(std::vector<Actin> &actinVec, const double k_off,
                                const double dt, const bool tether, int &nActin,
                                const bool dissociate, GactinGrid &gActinGrid,
                                bool arpPool, ArpGrid &arpGrid,
                                const bool steric, StericGrid &stericGrid,
                                std::vector<ExcZone> &excZones,
                                std::vector<MembraneWall> &memWalls,
                                std::vector<Membrane> &membranes,
                                Cortex &cortex,
                                const bool bDynamics,
                                std::vector<ProteinRegion> &branchRegions,
                                std::vector<ProteinRegion> &nucRegions,
                                std::vector<ProteinRegion> &capRegions,
                                std::vector<ProteinRegion> &sevRegions)
{
    // For now this in serial code but will parralise when figure out best way
    // to increment the gActingrid inside a parallel thread

    std::uniform_real_distribution<double> probDist(0, 1);
    std::vector<int> dissociateIDs;

    for (int i = 0; i < nActin; ++i)
    {
        if (!actinVec[i].getBarbedCapped())
        {
            double randDepoly { probDist(rng.m_mersenne) };
            double p_depoly = k_off * dt;
            if (randDepoly < p_depoly)
            {
                // Is the following line thread safe?

                if (dissociate && actinVec[i].getNumMonomers() == Actin::s_seedSize && actinVec[i].getParentID() == -1)
                {
                    /*
                    So this filament is going to break apart into monomers
                    It needs to be deleted
                    This is tricky because the id system we use is based on
                    filament's position in the actinVec. If we delete a fila
                    we alter this so that id != index in the array
                    we would have to subtract 1 from everyone's id in the array
                    It would also affect parent IDs - any daughters of filas
                    that come after filament i in the actinVec need to have their
                    parent id reduced by 1

                    Structure ids are fine

                    Would also have to deal with branches, its unlikely but possible
                    that it would still have a branch attached

                    Be careful because altering actinVec will affect this loop we
                    are in, especially when running in parallel

                    */

                    dissociateIDs.push_back(actinVec[i].getID());
                    gActinGrid.dissociate(actinVec[i].getBarbedEnd()[0], actinVec[i].getBarbedEnd()[1]);

                    continue;
                }
                else if (actinVec[i].getNumMonomers() == Actin::s_seedSize && actinVec[i].getParentID() == -1)
                {
                    // dissociation is turned off, but we dont want to depoly it
                    // just leave it
                    continue;
                }
                else if (dissociate && actinVec[i].getNumMonomers() == Actin::s_branchSeedSize && actinVec[i].getParentID() != -1)
                {
                    // It is a branch and dissociation is turned on, so it is removed
                    dissociateIDs.push_back(actinVec[i].getID());
                    gActinGrid.dissociateBR(actinVec[i].getPointedEnd()[0], actinVec[i].getPointedEnd()[1]);
                    if (arpPool)
                    {
                        // we have an arp grid and so must increment by 1
                        arpGrid.increment(actinVec[i].getBarbedEnd()[0], actinVec[i].getBarbedEnd()[1]);
                    }
                    continue;
                }
                else if (actinVec[i].getNumMonomers() == Actin::s_branchSeedSize && actinVec[i].getParentID() != -1)
                {
                    // It is a branch, we allow it to go down to 2 but branches
                    // cannot dissociate in our model
                    continue;
                }

                gActinGrid.increment(actinVec[i].getBarbedEnd()[0], actinVec[i].getBarbedEnd()[1]);


                actinVec[i].depolymerise_barbed(actinVec, nActin,
                                                arpPool,
                                                arpGrid, steric, stericGrid,
                                                excZones, memWalls, membranes,
                                                cortex,
                                                bDynamics,
                                                tether, 0, 0,
                                                branchRegions, nucRegions,
                                                capRegions, sevRegions);

                if (steric)
                    stericGrid.updateCellsBarbPoly(actinVec[i]);

            }
        }
    }

    dissociationRout(dissociateIDs, actinVec, nActin, steric, stericGrid);

}


void depolymerisationPointed(std::vector<Actin> &actinVec, const double k_off,
                            const double dt, const bool tether, int &nActin,
                            const bool dissociate,
                            bool arpPool, ArpGrid &arpGrid, const bool steric,
                            StericGrid &stericGrid,
                            std::vector<ExcZone> &excZones,
                            std::vector<MembraneWall> &memWalls,
                            std::vector<Membrane> &membranes,
                            Cortex &cortex,
                            const bool bDynamics,
                            std::vector<ProteinRegion> &branchRegions,
                            std::vector<ProteinRegion> &nucRegions,
                            std::vector<ProteinRegion> &capRegions,
                            std::vector<ProteinRegion> &sevRegions)
{

    std::uniform_real_distribution<double> probDist(0, 1);

    std::vector<int> dissociateIDs;

    //#pragma omp parallel for private(probDist)
    for (int i = 0; i < nActin; ++i)
    {
        int threadID = omp_get_thread_num();
        if (!actinVec[i].getPointedCapped())
        {
            double randDepoly { probDist(rng.m_PRNGs[threadID]) };
            double p_depoly = k_off * dt;
            if (randDepoly < p_depoly)
            {
                if (dissociate && actinVec[i].getNumMonomers() == Actin::s_seedSize)
                {
                    /*
                    So this filament is going to break apart into monomers
                    It needs to be deleted
                    This is tricky because the id system we use is based on
                    filament's position in the actinVec. If we delete a fila
                    we alter this so that id != index in the array
                    we would have to subtract 1 from everyone's id in the array
                    It would also affect parent IDs - any daughters of filas
                    that come after filament i in the actinVec need to have their
                    parent id reduced by 1

                    Structure ids are fine

                    Would also have to deal with branches, its unlikely but possible
                    that it would still have a branch attached

                    Be careful because altering actinVec will affect this loop we
                    are in, especially when running in parallel

                    */

                    dissociateIDs.push_back(actinVec[i].getID());
                    continue;
                }
                else if (actinVec[i].getNumMonomers() == Actin::s_seedSize)
                {
                    // dissociation is turned off, but we dont want to depoly it
                    // just leave it
                    continue;
                }


                actinVec[i].depolymerise_pointed(actinVec, nActin,
                                                 arpPool,
                                                 arpGrid,
                                                 steric, stericGrid,
                                                 excZones, memWalls,
                                                 membranes, cortex, bDynamics,
                                                 tether, 0, 0,
                                                 branchRegions, nucRegions,
                                                 capRegions, sevRegions);

                if (steric)
                    stericGrid.updateCellsPointPoly(actinVec[i]);

            }
        }
    }

    // After the loop, loop over the dissociateIDs vector
    dissociationRout(dissociateIDs, actinVec, nActin, steric, stericGrid);

}

void depolymerisationPointedForce(int i, std::vector<Actin> &actinVec,
                                  const bool tether, int &nActin,
                                  bool arpPool, ArpGrid &arpGrid,
                                  const bool steric,
                                  StericGrid &stericGrid,
                                  std::vector<ExcZone> &excZones,
                                  std::vector<MembraneWall> &memWalls,
                                  std::vector<Membrane> &membranes,
                                  Cortex &cortex,
                                  const bool bDynamics,
                                  std::vector<ProteinRegion> &branchRegions,
                                  std::vector<ProteinRegion> &nucRegions,
                                  std::vector<ProteinRegion> &capRegions,
                                  std::vector<ProteinRegion> &sevRegions)
{


    std::vector<int> dissociateIDs;

    if (actinVec[i].getNumMonomers() == Actin::s_seedSize)
    {
        /*
        So this filament is going to break apart into monomers
        It needs to be deleted
        This is tricky because the id system we use is based on
        filament's position in the actinVec. If we delete a fila
        we alter this so that id != index in the array
        we would have to subtract 1 from everyone's id in the array
        It would also affect parent IDs - any daughters of filas
        that come after filament i in the actinVec need to have their
        parent id reduced by 1

        Structure ids are fine

        Would also have to deal with branches, its unlikely but possible
        that it would still have a branch attached

        Be careful because altering actinVec will affect this loop we
        are in, especially when running in parallel

        */

        dissociateIDs.push_back(actinVec[i].getID());
        dissociationRout(dissociateIDs, actinVec, nActin, steric, stericGrid);
        return;
    }


    actinVec[i].depolymerise_pointed(actinVec, nActin,
                                     arpPool,
                                     arpGrid,
                                     steric, stericGrid,
                                     excZones, memWalls,
                                     membranes, cortex, bDynamics, tether, 0, 0,
                                     branchRegions, nucRegions,
                                     capRegions, sevRegions);

    if (steric)
        stericGrid.updateCellsPointPoly(actinVec[i]);



}

void depolymerisationPointedGrid(std::vector<Actin> &actinVec, const double k_off,
                                 const double dt, const bool tether, int &nActin,
                                 const bool dissociate, GactinGrid &gActinGrid,
                                 bool arpPool, ArpGrid &arpGrid,
                                 const bool steric, StericGrid &stericGrid,
                                 std::vector<ExcZone> &excZones,
                                 std::vector<MembraneWall> &memWalls,
                                 std::vector<Membrane> &membranes,
                                 Cortex &cortex,
                                 const bool bDynamics,
                                 std::vector<ProteinRegion> &branchRegions,
                                 std::vector<ProteinRegion> &nucRegions,
                                 std::vector<ProteinRegion> &capRegions,
                                 std::vector<ProteinRegion> &sevRegions)
{
    // Again, this is serial code but will parallelise

    std::uniform_real_distribution<double> probDist(0, 1);
    std::vector<int> dissociateIDs;

    for (int i = 0; i < nActin; ++i)
    {
        if (!actinVec[i].getPointedCapped())
        {
            double randDepoly { probDist(rng.m_mersenne) };
            double p_depoly = k_off * dt;
            if (randDepoly < p_depoly)
            {

                if (dissociate && actinVec[i].getNumMonomers() == Actin::s_seedSize)
                {
                    /*
                    So this filament is going to break apart into monomers
                    It needs to be deleted
                    This is tricky because the id system we use is based on
                    filament's position in the actinVec. If we delete a fila
                    we alter this so that id != index in the array
                    we would have to subtract 1 from everyone's id in the array
                    It would also affect parent IDs - any daughters of filas
                    that come after filament i in the actinVec need to have their
                    parent id reduced by 1

                    Structure ids are fine

                    Would also have to deal with branches, its unlikely but possible
                    that it would still have a branch attached

                    Be careful because altering actinVec will affect this loop we
                    are in, especially when running in parallel

                    */

                    dissociateIDs.push_back(actinVec[i].getID());
                    gActinGrid.dissociate(actinVec[i].getPointedEnd()[0], actinVec[i].getPointedEnd()[1]);
                    continue;
                }
                else if (actinVec[i].getNumMonomers() == Actin::s_seedSize)
                {
                    // dissociation is turned off, but we dont want to depoly it
                    // just leave it
                    continue;
                }

                gActinGrid.increment(actinVec[i].getPointedEnd()[0], actinVec[i].getPointedEnd()[1]);


                actinVec[i].depolymerise_pointed(actinVec, nActin,
                                                 arpPool,
                                                 arpGrid, steric, stericGrid,
                                                 excZones, memWalls,
                                                 membranes, cortex, bDynamics,
                                                 tether, 0, 0,
                                                 branchRegions, nucRegions,
                                                 capRegions, sevRegions);

                if (steric)
                    stericGrid.updateCellsPointPoly(actinVec[i]);

            }
        }
    }
    // After the loop, loop over the dissociateIDs vector
    dissociationRout(dissociateIDs, actinVec, nActin, steric, stericGrid);

}


void autoDebranchBarb(Actin &parent, int &nActin, std::vector<Actin> &actinVec,
                      bool arpPool, ArpGrid &arpGrid, const bool steric,
                      StericGrid &stericGrid)
{
    // The filament has branches
    std::vector<int> dissociateIDs;

    for (int i = 0; i < parent.getNumSubs(); ++i)
    {
        for (unsigned int j = 0; j < parent.getBranchIDSubVector(i).size(); ++j)
        {
            int k = parent.getBranchIDSubVector(i)[j];
            if (actinVec[k].getMotherMonoID() >= parent.getNumMonomers())
            {
                // Branch lies beyond barbed end, needs to detach
                std::array<double,3> branchCoords = actinVec[k].getPointedEnd();
                std::array<double,2> branchCoords2D;
                branchCoords2D[0] = branchCoords[0];
                branchCoords2D[1] = branchCoords[1];

                actinVec[k].detach(actinVec, nActin);

                if (arpPool)
                    arpGrid.increment(branchCoords2D[0], branchCoords2D[1]);

                if (actinVec[k].getNumMonomers() < Actin::s_seedSize)
                {
                    // it should dissociate
                    dissociateIDs.push_back(actinVec[k].getID());
                }

                parent.recalcDistToBR2(actinVec);
                parent.recalcMonoCompat(actinVec, parent.getNumMonomers());

            }
        }
    }

    dissociationRout(dissociateIDs, actinVec, nActin, steric, stericGrid);

}


void autoDebranchPoint(Actin &parent, int &nActin, std::vector<Actin> &actinVec,
                       bool arpPool, ArpGrid &arpGrid, const bool steric,
                       StericGrid &stericGrid)
{
    // The filament has branches
    std::vector<int> dissociateIDs;

    for (int i = 0; i < parent.getNumSubs(); ++i)
    {
        for (unsigned int j = 0; j < parent.getBranchIDSubVector(i).size(); ++j)
        {
            int k = parent.getBranchIDSubVector(i)[j];
            if (actinVec[k].getMotherMonoID() < 0)
            {
                std::array<double,3> branchCoords = actinVec[k].getPointedEnd();
                std::array<double,2> branchCoords2D;
                branchCoords2D[0] = branchCoords[0];
                branchCoords2D[1] = branchCoords[1];

                actinVec[k].detach(actinVec, nActin);


                if (arpPool)
                    arpGrid.increment(branchCoords2D[0], branchCoords2D[1]);

                if (actinVec[k].getNumMonomers() < Actin::s_seedSize)
                {
                    // it should dissociate
                    dissociateIDs.push_back(actinVec[k].getID());
                }
            }

        }
    }
    dissociationRout(dissociateIDs, actinVec, nActin, steric, stericGrid);

}

void dissociationRout(std::vector<int> dissociateIDs,
                      std::vector<Actin> &actinVec, int &nActin,
                      const bool steric, StericGrid &stericGrid)
{
    for (unsigned int i = 0; i < dissociateIDs.size(); ++i)
    {
        int id = dissociateIDs[i];

        int parentID = actinVec[id].getParentID();

        // Step 1, if it is a branch, remove it from its mother
        if (parentID != -1)
        {
            // it is a branch, have to take extra steps
            // Easier if we detach it before deleting it to avoid floating pointers
            int branchMonoID = actinVec[id].getMotherMonoID();


            std::array<double,3> branchCoords = actinVec[id].getPointedEnd();
            std::array<double,2> branchCoords2D;
            branchCoords2D[0] = branchCoords[0];
            branchCoords2D[1] = branchCoords[1];

            actinVec[id].detach(actinVec, nActin);


            // recalc distance to leading branch, although this might not change!
            // if it "was" the leading branch...
            if (branchMonoID == actinVec[parentID].getNumMonomers()-1-actinVec[parentID].getDistToLeadBR())
            {
                int olddist = actinVec[parentID].getDistToLeadBR();
                actinVec[parentID].recalcDistToBR2(actinVec);
                assert(olddist != actinVec[parentID].getDistToLeadBR());
            }
            actinVec[parentID].recalcMonoCompat(actinVec, branchMonoID);
        }


        if (steric)
            for (int j = 0; j < actinVec[id].getNumSubs(); ++j)
            {
                // reset cells
                stericGrid.resetCells(actinVec[id], j);
            }

        // Check if any crosslinks exist on this fila?
        while (actinVec[id].getNumCLinks() > 0)
        {
            unLinkBase(id, 0, actinVec);
        }
        assert(actinVec[id].getNumCLinks() == 0);

        // Check if any branches exist on this fila
        bool removeDaughter = false;
        if (actinVec[id].getDaughternum() != 0)
        {
            // Need to remove branches from the filament
            // Should only be one
            for (int j = 0; j < actinVec[id].getNumSubs(); ++j)
            {
                for (unsigned int k = 0; k < actinVec[id].getBranchIDSubVector(j).size(); ++k)
                {
                    int daughterID = actinVec[id].getBranchIDSubVector(j)[k];
                    // What if the daughter has only 2 monomers? This too needs to dissociate
                    if (actinVec[daughterID].getNumMonomers() < Actin::s_seedSize)
                    {
                        std::vector<int> dissociateIDs2;
                        dissociateIDs2.push_back(daughterID);
                        dissociationRout(dissociateIDs2, actinVec, nActin, steric,
                                         stericGrid);
                        removeDaughter = true;
                    }
                    else
                    {
                        actinVec[daughterID].detach(actinVec, nActin);
                        removeDaughter = true;
                    }
                }

                if (removeDaughter)
                    break;
            }
        }
        assert(actinVec[id].getDaughternum() == 0);

        // Have to be careful here also
        for (unsigned int j = id+1; j < actinVec.size(); ++j)
        {
            // Step 1, reduce the id numbers by 1
            // this also reduces parent ids of any daughters, have to be careful
            actinVec[j].reduceID(actinVec);
            if (steric)
                stericGrid.dissociate(actinVec[j]);
        }

        for (unsigned int j = 0; j < actinVec.size(); ++j)
        {
            // Step 2, reduce any daughterIDs
            actinVec[j].checkandReduceDaughtID(id);
        }



        // Reduce the crosslink pointers
        // Also reduce the crosslinked filaments pointers
        for (unsigned int j = 0; j < actinVec.size(); ++j)
        {
            for (int k = 0; k < actinVec[j].getNumCLinks(); ++k)
            {
                // have to do both filaments at once
                int otherFila = actinVec[j].getCLinkActinAndSite(k)[2];

                if (otherFila > id)
                {
                    // the "other filament" will have just had its id lowered,
                    // this needs to be changed in m_cLinkActinAndSites
                    actinVec[j].setClOther(k, otherFila-1);
                }
            }
        }

        // Step 3, erase the element in the vector, deleting the extra filament
        actinVec.erase(actinVec.begin()+id);

        nActin--;

        // Adjust the higher dissociateIDs entries (if there are any)
        for (unsigned int j = i+1; j < dissociateIDs.size(); ++j)
        {
            dissociateIDs[j] -= 1;
        }

        Actin::s_total_subunits -= 3;
    }
}
