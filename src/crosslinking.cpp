/*
 *  crosslinking.cpp
 *
 *
 *  Functions that create crosslinks between filaments
 *
 *  James Bradford
 *  University of Sheffield
 *  Feb 2020
 */

#include "configHeader.h"
#include "RNG.h"
#include "Actin.h"
#include "StericGrid.h"
#include "crosslinking.h"
#include "geometry.h"
extern RNG rng;

void crosslinkingRoutine(std::vector<Actin> &actinVec, const double k_cLink,
                         const double tStep, const double currTime,
                         const bool steric, StericGrid &stericGrid,
                         bool &crossLinking)
{
   /* First attempt at crosslinking subroutine
    * Crosslinking discrete sites are per monomer, similar to branches
    * This makes it easier when doing things like severing
    *
    */

   double p_cLink_base = k_cLink*tStep; // probability of forming crosslink per site

   for (unsigned int i = 0; i < actinVec.size(); ++i)
   {

       if (actinVec[i].getNumMonomers() < 4)
       {
            // AVOID due to error this causes when identifying which subunit the
            // crosslink is on. CSpacing will regularly be more than 3 anyway!
            continue;
       }


       // First get the number of available sites on the i'th filament
       std::vector< std::array<int,5> > possCrossLinks;
       int numAvailSites = actinVec[i].getAvailMonoVec().size();

       // For each of these sites calculate the number of nearby available sites
       // on other filaments. Nearby meaning less than the "Actin::s_cLDist" threshold

       for (int j = 0; j < numAvailSites; ++j)
       {
           // Use the steric grid to find nearby actin subunits

           // what subunit is j'th site on?
           int monoID = actinVec[i].getAvailMonoID(j);

           assert(monoID < actinVec[i].getNumMonomers());

           assert(actinVec[i].getMonoCompat(monoID));


           double lengthAlongSub;
           int subID = actinVec[i].findSubunit(monoID, lengthAlongSub);

           // What is close by?
           std::array<double,2> jSitePoint = actinVec[i].findMonoCoord(monoID);

           // Find steric grid cell that monomer is in
           int cellID = stericGrid.getCellFromCoord(jSitePoint);

           std::vector<int> cellIDs;
           cellIDs.push_back(cellID); // quick and dirty way
           std::vector<int> cellsToCheck = stericGrid.getCellsToCheckBIG(cellIDs, Actin::s_cLDist); // gets 8 surrounding cells plus cellID

           std::vector<std::array<int,2> > checkedSubs;

           for (unsigned int k = 0; k < cellsToCheck.size(); ++k)
           {
                   std::vector<std::array<int,2> > cellContents = stericGrid.getCellContents(cellsToCheck[k]);

                   for (unsigned int l = 0; l < cellContents.size(); ++l)
                   {

                           int actinID = cellContents[l][0];
                           int actinSubID = cellContents[l][1];

                           if (actinID == actinVec[i].getID())
                           {
                               // Same filament, ignore
                               continue;
                           }

                           if (actinVec[actinID].getNumMonomers() < 4)
                           {
                                // AVOID
                                continue;
                           }

                           if (actinID < i)
                           {
                                // We have checked this filament's potential crosslinks already
                                // Avoid to avoid double counting each potential clink
                                continue;
                           }

                           if ((actinVec[i].getParentID() == actinID && subID <= 1) || (actinVec[actinID].getParentID() == i && actinSubID <= 1))
                           {
                                // no point connecting branched actin to mothers, already a branch here!
                                continue;
                           }
                           if (std::find(checkedSubs.begin(), checkedSubs.end(), cellContents[l]) != checkedSubs.end())
                           {
                                   // We have checked this subfilament before!
                                   continue;
                           }
                           checkedSubs.push_back(cellContents[l]);

                           std::array<int,2> neighbourMonoLims = actinVec[actinID].findMonoIDLimitsOfSub(actinSubID);
                           std::vector<int> availMonos;
                           for (int m = neighbourMonoLims[0]; m <= neighbourMonoLims[1]; ++m)
                           {

                              double junk;
                              assert(actinVec[actinID].findSubunit(m, junk) == actinSubID);

                              if (actinVec[actinID].getMonoCompat(m))
                              {
                                  // If this monomer is available to link
                                  availMonos.push_back(m);
                              }
                           }
                           // Find unoccupied sites on filament "actinID"
                           // Would be better to do this on subunit level but for now do this
                           int numAvailSitesNeighbourSub = availMonos.size();


                           for (int m = 0; m < numAvailSitesNeighbourSub; ++m) // should restrict to just the sites on actinSubID
                           {
                               // Check distance
                               int idx = availMonos[m];
                               assert(actinVec[actinID].getMonoCompat(idx));

                               std::array<double,2> mSitePoint = actinVec[actinID].findMonoCoord(idx);

                               double distance = sqrt(distanceBetPoints2DSQR(jSitePoint, mSitePoint));
                               if (distance < Actin::s_cLDist + actinVec[i].getStericRadius() + actinVec[actinID].getStericRadius())
                               {
                                   // Add to possible crosslinks
                                   std::array<int,5> possCLink;
                                   possCLink[0] = monoID;
                                   possCLink[1] = subID;
                                   possCLink[2] = actinID;
                                   possCLink[3] = idx;
                                   possCLink[4] = actinSubID;
                                   possCrossLinks.push_back(possCLink);
                               }


                           }


                   }


           }

       }

       // Now have a vector of possible crosslinks
       // MC Step to determine if one of these is created
       int numPossLinks = possCrossLinks.size();
       if (numPossLinks == 0)
         return;

       double p_cLink = p_cLink_base * numPossLinks;
       assert(p_cLink < 0.1); // statistics check if this is flagging up then dt needs to be decreased

       double randCLink { rng.m_probDist(rng.m_mersenne)};

       if (randCLink < p_cLink)
       {
           // Randomly choose one of the possible links for this filament
           // and create it

           std::uniform_int_distribution<int> possLinks(0,numPossLinks-1);
           std::array<int,5> chosenLink = possCrossLinks[possLinks(rng.m_mersenne)];

           actinVec[i].createCrossLink(actinVec[chosenLink[2]], chosenLink[0],
                                       chosenLink[1], chosenLink[3],
                                       chosenLink[4], actinVec);
       }
   }
}

void unLinkBase(int i, int j, std::vector<Actin> &actinVec)
{
    // i is the filament id, j is the crosslink id

    // Remove this crosslink
    int monoID = actinVec[i].getCLinkActinAndSite(j)[0];
    int otherFila = actinVec[i].getCLinkActinAndSite(j)[2];
    int otherFilaIDX = actinVec[i].getCLinkActinAndSite(j)[3];
    int otherMonoID = actinVec[i].getCLinkActinAndSite(j)[4];



    actinVec[i].eraseCLink(j);
    for (int k = 0; k < actinVec[i].getNumCLinks(); ++k)
    {
      // Make sure any linked actin still point to the right cl idx

      int otherFila_k = actinVec[i].getCLinkActinAndSite(k)[2];
      int otherFilaCLIDX_k = actinVec[i].getCLinkActinAndSite(k)[3];
      actinVec[otherFila_k].setClOtherIDX(otherFilaCLIDX_k, k);
    }

    actinVec[i].setDistToLeadCL(actinVec[i].getNumMonomers());
    for (int k = 0; k < actinVec[i].getNumCLinks(); ++k)
    {
        int monomerID = actinVec[i].getCLinkActinAndSite(k)[0];
        actinVec[i].recalcDistToLeadCL(monomerID);
    }



    actinVec[otherFila].eraseCLink(otherFilaIDX);
    for (int k = 0; k < actinVec[otherFila].getNumCLinks(); ++k)
    {
      // Make sure any linked actin still point to the right cl idx

      int otherFila_k = actinVec[otherFila].getCLinkActinAndSite(k)[2];
      int otherFilaCLIDX_k = actinVec[otherFila].getCLinkActinAndSite(k)[3];
      actinVec[otherFila_k].setClOtherIDX(otherFilaCLIDX_k, k);
    }

    actinVec[otherFila].setDistToLeadCL(actinVec[otherFila].getNumMonomers());
    for (int k = 0; k < actinVec[otherFila].getNumCLinks(); ++k)
    {
        int monomerID = actinVec[otherFila].getCLinkActinAndSite(k)[0];
        actinVec[otherFila].recalcDistToLeadCL(monomerID);
    }


    actinVec[i].recalcMonoCompat(actinVec, monoID);
    actinVec[otherFila].recalcMonoCompat(actinVec, otherMonoID);
    if (actinVec[i].getCLStructureID() != actinVec[otherFila].getCLStructureID())
    {
        // different branched structures

        // Do the structure ids need to changed, are the two filaments not
        // associated with one another anymore?

        if (!actinVec[i].checkForSameStructure(actinVec, otherFila))
        {
            // They are no longer connected
            // Where is the current master filament?
            actinVec[i].setStructureID(Actin::s_structure_idGenerator++);
            actinVec[i].changeStructure(actinVec, actinVec.size());

            if (!actinVec[i].checkForMasterInStructure(actinVec))
            {
                // No master in i now
                actinVec[i].findNewMaster(actinVec);
            }
            else if (!actinVec[otherFila].checkForMasterInStructure(actinVec))
            {
                // No master in other now
                actinVec[otherFila].findNewMaster(actinVec);
            }

        }
        else if (actinVec[i].getCLinkMaster() && actinVec[i].getNumCLinks() == 0)
        {
            actinVec[i].setSelfMasterBoolFalse();
            actinVec[i].findNewMaster(actinVec);
        }
        else if (actinVec[otherFila].getCLinkMaster() && actinVec[otherFila].getNumCLinks() == 0)
        {
            actinVec[otherFila].setSelfMasterBoolFalse();
            actinVec[otherFila].findNewMaster(actinVec);
        }

    }
    else if (actinVec[i].getCLinkMaster() && actinVec[i].getNumCLinks() == 0)
    {
        actinVec[i].setSelfMasterBoolFalse();
        actinVec[i].findNewMaster(actinVec);
    }
    else if (actinVec[otherFila].getCLinkMaster() && actinVec[otherFila].getNumCLinks() == 0)
    {
        actinVec[otherFila].setSelfMasterBoolFalse();
        actinVec[otherFila].findNewMaster(actinVec);
    }

}

void unLink(std::vector<Actin> &actinVec, const double k_unLink,
            const double tStep, const double currTime,
            const bool steric, StericGrid &stericGrid,
            bool &crossLinking)
{
  /*
   * Function that removes a crossLink based on a flat off rate
   *
   */

   double p_unLink = k_unLink*tStep;

   for (unsigned int i = 0; i < actinVec.size(); ++i)
   {
      int numClinks = actinVec[i].getNumCLinks();

      for (int j = 0; j < numClinks; ++j)
      {
          // each crosslink has p_unLink chance of being removed
          double randUnLink { rng.m_probDist(rng.m_mersenne)};

          if (randUnLink < p_unLink)
          {
              unLinkBase(i, j, actinVec);
              break;
          }
      }
   }
}

void autoUnlinkBarb(Actin &actin, std::vector<Actin> &actinVec)
{
    // The filament has crosslinks
    for (int i = 0; i < actin.getNumCLinks(); ++i)
    {
        if (actin.getCLinkActinAndSite(i)[0] >= actin.getNumMonomers())
        {
            // Needs to unlink
            // Remove this crosslink
            unLinkBase(actin.getID(), i, actinVec);
            break;
        }
    }
}

void autoUnlinkPoint(Actin &actin, std::vector<Actin> &actinVec)
{
  // The filament has crosslinks

  for (int i = 0; i < actin.getNumCLinks(); ++i)
  {
      if (actin.getCLinkActinAndSite(i)[0] < 0)
      {
          // Needs to unlink
          // Remove this crosslink
          unLinkBase(actin.getID(), i, actinVec);
          break;
      }
  }
}
