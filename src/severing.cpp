/*
 *  severing.cpp
 *
 *
 *  Functions that sever filaments
 *
 *  James Bradford
 *  University of Sheffield
 *  June 2017
 */

#include "configHeader.h"
#include "GactinGrid.h"
#include "StericGrid.h"
#include "ArpGrid.h"
#include "severing.h"
#include "Actin.h"
#include "RNG.h"

extern RNG rng;

void severfilament(std::vector<Actin> &actinvec, const double k_sever,
                   const double t_step, int &nActin, const int nActinLimit,
                   bool &severing, const double currTime, const bool steric,
                   StericGrid &stericGrid,
                   GactinGrid &gActinGrid, const bool arpPool, ArpGrid &arpGrid)
{
    // Serial version
    if (nActin == nActinLimit)
    {
        // Cannot sever anymore reached the nActinLimit
        return;
    }


    int imax = nActin;
    for (int i = 0; i < imax; ++i)
    {
        double severProb = actinvec[i].getLength() * 1E6 * k_sever * t_step;
        double randSever { rng.m_probDist(rng.m_mersenne) };

        if (randSever < severProb)
        {
            // So this filament is going to sever, need to randomly select where
            std::uniform_int_distribution<> severDist(0, actinvec[i].getNumMonomers()-2);
            int severMonomer = severDist(rng.m_mersenne);

            actinvec[i].sever(actinvec, currTime, severMonomer, nActin, steric,
                              stericGrid, gActinGrid, arpPool, arpGrid);
        }
    }
}


void severfilamentRegions(std::vector<Actin> &actinvec, const double k_sever,
                          const double t_step, int &nActin, const int nActinLimit,
                          bool &severing, const double currTime, const bool steric,
                          StericGrid &stericGrid,
                          GactinGrid &gActinGrid, const bool arpPool, ArpGrid &arpGrid,
                          std::vector<ProteinRegion> &sevRegions)
{
    if (nActin == nActinLimit)
    {
        // Cannot sever anymore reached the nActinLimit
        std::cout << "Reached the limit of number of filaments: ";
        std::cout << nActinLimit << std::endl;
        severing = false;
        return;
    }

    int imax = nActin;
    int nRegions = sevRegions.size();
    for (int i = 0; i < imax; ++i)
    {
        for (int j = 0; j < nRegions; ++j)
        {
            std::vector<int> inReg;

            #pragma omp parallel
            {
                std::vector<int> inReg_priv;
                #pragma omp for nowait
                for (unsigned int k = 0; k < actinvec[i].getNumMonomers(); ++k)
                {
                    std::array<double,2> monoPoint = actinvec[i].findMonoCoord(k);
                    if (sevRegions[j].checkPointwithin(monoPoint))
                    {
                        inReg_priv.push_back(k);
                    }

                }
                #pragma omp critical
                inReg.insert(inReg.end(), inReg_priv.begin(), inReg_priv.end());
            }
            std::sort(inReg.begin(), inReg.end());

            double severProb = inReg.size() * Actin::s_monomerLength * 1E6 * k_sever * sevRegions[j].getSevCoeff() * t_step;
            double randSever { rng.m_probDist(rng.m_mersenne) };

            if (randSever < severProb)
            {
                std::uniform_int_distribution<int> availMonos(0,inReg.size()-1);
                int severMonomer = inReg[availMonos(rng.m_mersenne)];

                actinvec[i].sever(actinvec, currTime, severMonomer, nActin, steric,
                                        stericGrid, gActinGrid, arpPool, arpGrid);
            }
        }
    }
}
