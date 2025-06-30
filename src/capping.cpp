/*
 *  capping.cpp
 *
 *
 *  Functions that caps and uncaps filaments
 *
 *  James Bradford
 *  University of Sheffield
 *  May 2017
 */

#include "configHeader.h"
#include "capping.h"
#include "RNG.h"
#include "Actin.h"

extern RNG rng;


void capBarbed(std::vector<Actin> &actinvec, const double pCap)
{
    int nActin = actinvec.size();
    std::uniform_real_distribution<double> probDist(0, 1);

    #pragma omp parallel for private(probDist)
    for (int i = 0; i < nActin; ++i)
    {
        int threadID = omp_get_thread_num();
        if (actinvec[i].getBarbedCapped())
            continue;

        double randCap { probDist(rng.m_PRNGs[threadID]) };
        if (randCap < pCap)
            actinvec[i].setBarbedCapped();
    }
}

void capBarbedRegion(std::vector<Actin> &actinvec, double pCap,
                     const std::vector<ProteinRegion> &capRegions)
{
    for (Actin &actinobj : actinvec)
    {
        if (actinobj.getBarbedCapped())
            continue;

        std::array<double,2> barbedPoint { actinobj.getBarbedEnd()[0], actinobj.getBarbedEnd()[1] };

        for (unsigned int i = 0; i < capRegions.size(); ++i)
        {
            if (capRegions[i].checkPointwithin(barbedPoint))
            {
                double randCap { rng.m_probDist(rng.m_mersenne) };
                pCap *= capRegions[i].getCapCoeff();
                if (randCap < pCap)
                    actinobj.setBarbedCapped();
            }
        }
    }
}
void capBarbedAntiRegion(std::vector<Actin> &actinvec, double pCap,
                     const std::vector<ProteinRegion> &antiCapRegions)
{
    for (Actin &actinobj : actinvec)
    {
        if (actinobj.getBarbedCapped())
            continue;

        std::array<double,2> barbedPoint { actinobj.getBarbedEnd()[0], actinobj.getBarbedEnd()[1] };
        bool inRegion = false;
        // Need to make sure the barbed point is not in ANY antiCapRegion
        for (unsigned int i = 0; i < antiCapRegions.size(); ++i)
        {
            if (antiCapRegions[i].checkPointwithin(barbedPoint))
            {
                // it is in a anti cap region, end this
                inRegion = true;
                break;
            }
        }


        if (!inRegion)
        {
            double randCap { rng.m_probDist(rng.m_mersenne) };
            if (randCap < pCap)
            {
                actinobj.setBarbedCapped();
            }
        }

    }
}

void capPointed(std::vector<Actin> &actinvec, const double pCap)
{
    int nActin = actinvec.size();
    std::uniform_real_distribution<double> probDist(0, 1);

    #pragma omp parallel for private(probDist)
    for (int i = 0; i < nActin; ++i)
    {
        int threadID = omp_get_thread_num();
        if (actinvec[i].getPointedCapped())
            continue;

        double randCap { probDist(rng.m_PRNGs[threadID]) };
        if (randCap < pCap)
            actinvec[i].setPointedCapped();
    }
}

void capPointedRegion(std::vector<Actin> &actinvec, double pCap,
                     const std::vector<ProteinRegion> &capRegions)
{
    for (Actin &actinobj : actinvec)
    {
        if (actinobj.getPointedCapped())
            continue;

        std::array<double,2> pointedPoint { actinobj.getPointedEnd()[0], actinobj.getPointedEnd()[1] };

        for (unsigned int i = 0; i < capRegions.size(); ++i)
        {
            if (capRegions[i].checkPointwithin(pointedPoint))
            {
                pCap *= capRegions[i].getCapCoeff();
                double randCap { rng.m_probDist(rng.m_mersenne) };
                if (randCap < pCap)
                    actinobj.setPointedCapped();
            }
        }
    }
}
void capPointedAntiRegion(std::vector<Actin> &actinvec, double pCap,
                     const std::vector<ProteinRegion> &antiCapRegions)
{
    for (Actin &actinobj : actinvec)
    {
        if (actinobj.getPointedCapped())
            continue;

        std::array<double,2> pointedPoint { actinobj.getPointedEnd()[0], actinobj.getPointedEnd()[1] };
        bool inRegion = false;
        for (unsigned int i = 0; i < antiCapRegions.size(); ++i)
        {
            if (antiCapRegions[i].checkPointwithin(pointedPoint))
            {
                inRegion = true;
                break;
            }

        }

        if (!inRegion)
        {
            double randCap { rng.m_probDist(rng.m_mersenne) };
            if (randCap < pCap)
                actinobj.setPointedCapped();
        }
    }
}

void uncapBarbed(std::vector<Actin> &actinvec, const double pUncap)
{
    int nActin = actinvec.size();
    std::uniform_real_distribution<double> probDist(0, 1);

    #pragma omp parallel for private(probDist)
    for (int i = 0; i < nActin; ++i)
    {
        int threadID = omp_get_thread_num();
        if (!actinvec[i].getBarbedCapped())
            continue;

        double randUncap { probDist(rng.m_PRNGs[threadID]) };
        if (randUncap < pUncap)
            actinvec[i].setBarbedUncapped();
    }
}

void uncapPointed(std::vector<Actin> &actinvec, const double pUncap)
{
    int nActin = actinvec.size();
    std::uniform_real_distribution<double> probDist(0, 1);

    #pragma omp parallel for private(probDist)
    for (int i = 0; i < nActin; ++i)
    {
        int threadID = omp_get_thread_num();
        if (!actinvec[i].getPointedCapped())
            continue;

        double randUncap { probDist(rng.m_PRNGs[threadID]) };
        if ((actinvec[i].getParentID() == -1) && (randUncap < pUncap))
            actinvec[i].setPointedUncapped();
    }
}
