/*
 *  Actin.cpp
 *
 *  C++ file containing the definition of the Actin class.
 *  This contains member variable initialisations and non-trivial member
 *  functions.
 *  The declaration of this class (including the definitions of the trivial
 *  member functions) is found in the associated header file "Actin.h"
 *
 *  James Bradford
 *  University of Sheffield
 *  March 2020
 */

#include <cassert> // for assert()
#include "Actin.h"
#include "branching.h"
#include "ellipticsolver.h"
#include "ellipseconfidence.h"
#include "RNG.h"
#include "ProteinRegion.h"
#include "ExcZone.h"
#include "geometry.h"
#include "globals.h"
#include "MembraneWall.h"
#include "Membrane.h"
#include "polymerisation.h"
#include "crosslinking.h"
#include "phagoAnalysis.h"
#include "Cortex.h"
#include <fstream>

extern RNG rng;

typedef gte::DCPQuery<double, gte::Segment<2, double>, gte::Segment<2, double> > RobustQuery;

// Default constructor for single filaments (includes branches but not bundles!)
Actin::Actin() //Constructor 0
        : m_radius { s_trueRadius },
        m_steric { true },
        m_stericRadius { s_stericRadius },
        m_daughter_num { 0 },
        m_capped_b { false },
        m_num_subunits { 3 },
        m_des_subunit_len { s_segmentationLength },
        m_changed { false },
        m_flexible { false },
        m_cLinkMaster { false }
{
        // start with an empty vector
        m_daughterIDs.push_back(std::vector<int>());
        m_daughterIDs.push_back(std::vector<int>());
        m_daughterIDs.push_back(std::vector<int>());

        // For steric grid
        m_stericCells.push_back(std::vector<int>());
        m_stericCells.push_back(std::vector<int>());
        m_stericCells.push_back(std::vector<int>());

        m_Darray = {-1, -1, -1};

}

// filaments nucleated within a region
Actin::Actin(int ID, std::uniform_real_distribution<double> dist1,
             std::uniform_real_distribution<double> dist2, int regionID,
             const std::vector<ProteinRegion> &nucRegions,
             double birthtime) // Constructor 1
        : Actin() // Calls Constructor 0
{
        m_id =  ID;
        m_length = s_seedSize * s_monomerLength;
        m_num_monomers = s_seedSize;
        m_distToLeadBR = m_num_monomers;
        m_distToLeadCL = m_num_monomers;

        for (int i = 0; i < m_num_monomers; ++i)
        {
            if (i < s_maxSpacing)
            {
                m_monosCompat.push_back(false);
            }
            else
            {
                m_monosCompat.push_back(true);
                m_availMonos.push_back(i);
            }
        }


        // Initialise with 4 points, the generated point is the second one

        m_points = {{0,0,0},
                    {0,0,0},
                    {0,0,0},
                    {0,0,0}};

        if (nucRegions[regionID].getCircBool() || nucRegions[regionID].getRingBool())
        {
                double r = dist1(rng.m_mersenne);
                double theta = dist2(rng.m_mersenne);

                m_points[1][0] = nucRegions[regionID].getCentre()[0] + r*cos(theta);
                m_points[1][1] = nucRegions[regionID].getCentre()[1] + r*sin(theta);

        }
        else
        {
                double x_loc = dist1(rng.m_mersenne);
                double y_loc = dist2(rng.m_mersenne);

                m_points[1][0] = nucRegions[regionID].getBLPoint()[0] + ((x_loc*nucRegions[regionID].getUVec()[0]) - (y_loc*nucRegions[regionID].getUVec()[1]));
                m_points[1][1] = nucRegions[regionID].getBLPoint()[1] + ((x_loc*nucRegions[regionID].getUVec()[1]) + (y_loc*nucRegions[regionID].getUVec()[0]));
        }



        double xyangle = rng.m_angledist(rng.m_mersenne);
        m_subunit_unitVecs.push_back({cos(xyangle),sin(xyangle)});
        m_subunit_unitVecs.push_back({cos(xyangle),sin(xyangle)});
        m_subunit_unitVecs.push_back({cos(xyangle),sin(xyangle)});

        m_birthtime = birthtime;
        for (int i = 0; i < m_num_monomers; ++i)
        {
            m_monoBirthTime.push_back(birthtime);
        }

        m_parent_id = -1; // invalid parent id means no parent
        m_structure_id = s_structure_idGenerator++;

        m_subStructure_id = s_subStructure_idGenerator++;

        m_crossStructure_id = s_crossStructure_idGenerator++;

        m_capped_p = false;

        m_points[0] = { m_points[1][0] - ((m_length/4) * m_subunit_unitVecs[0][0]), m_points[1][1] - ((m_length/4) * m_subunit_unitVecs[0][1]), 0 };

        m_points[2] = { m_points[1][0] + ((m_length/2) * m_subunit_unitVecs[1][0]), m_points[1][1] + ((m_length/2) * m_subunit_unitVecs[1][1]), 0 };
        m_points[3] = { m_points[2][0] + ((m_length/4) * m_subunit_unitVecs[2][0]), m_points[2][1] + ((m_length/4) * m_subunit_unitVecs[2][1]), 0 };


        m_actualSubLengths.push_back( m_length / 4);
        m_actualSubLengths.push_back( m_length / 2);
        m_actualSubLengths.push_back( m_length / 4);

        m_prescribedSubLengths.push_back( m_length / 4);
        m_prescribedSubLengths.push_back( m_length / 2);
        m_prescribedSubLengths.push_back( m_length / 4);

        s_total_subunits += 3;

        m_regionID = regionID;

}

// Regional but with directionality determine by gaussian
Actin::Actin(int ID, std::uniform_real_distribution<double> xpos_dist,
             std::uniform_real_distribution<double> ypos_dist, int regionID,
             const std::vector<ProteinRegion> &nucRegions,
             std::normal_distribution<double> angle_dist, double birthtime) // Constructor 2
        : Actin(ID, xpos_dist, ypos_dist, regionID, nucRegions, birthtime) // Calls Constructor 1
{
        // Change the angle, which obv changes the end points

        double xyangle = angle_dist(rng.m_mersenne);
        m_subunit_unitVecs[0] = {cos(xyangle),sin(xyangle)};
        m_subunit_unitVecs[1] = {cos(xyangle),sin(xyangle)};
        m_subunit_unitVecs[2] = {cos(xyangle),sin(xyangle)};

        m_points[0] = { m_points[1][0] - ((m_length/4) * m_subunit_unitVecs[0][0]), m_points[1][1] - ((m_length/4) * m_subunit_unitVecs[0][1]), 0 };

        m_points[2] = { m_points[1][0] + ((m_length/2) * m_subunit_unitVecs[1][0]), m_points[1][1] + ((m_length/2) * m_subunit_unitVecs[1][1]), 0 };
        m_points[3] = { m_points[2][0] + ((m_length/4) * m_subunit_unitVecs[2][0]), m_points[2][1] + ((m_length/4) * m_subunit_unitVecs[2][1]), 0 };

}

// Constructor for branched filaments
Actin::Actin(int ID, double birthtime, Actin &parent, int motherMonomer) // Constructor 3
        : Actin() // Calls Constructor 0
{

        m_id = ID;
        m_length = s_branchSeedSize * s_monomerLength;
        m_num_monomers = s_branchSeedSize;
        m_distToLeadBR = m_num_monomers;
        m_distToLeadCL = m_num_monomers;

        for (int i = 0; i < m_num_monomers; ++i)
        {
            if (i < s_maxSpacing)
            {
                m_monosCompat.push_back(false);
            }
            else
            {
                m_monosCompat.push_back(true);
                m_availMonos.push_back(i);
            }
        }


        m_birthtime = birthtime;
        for (int i = 0; i < m_num_monomers; ++i)
        {
            m_monoBirthTime.push_back(birthtime);
        }

        m_parent_id = parent.getID();
        m_structure_id = parent.getStructureID();
        m_subStructure_id = parent.getSubStructureID();
        m_crossStructure_id = parent.getCLStructureID();
        m_capped_p = true;
        m_motherMonomer = motherMonomer;

        m_branchSubUnit = parent.findSubunit(motherMonomer, m_branchSubLen);

        parent.addBranchtoParentVector(m_branchSubUnit, ID);

        m_branchdir = branchDir();
        Eigen::Matrix2d rotMatrix;
        rotMatrix << cos(m_branchdir*s_branchAngle), -sin(m_branchdir*s_branchAngle),
                sin(m_branchdir*s_branchAngle), cos(m_branchdir*s_branchAngle);

        std::array<double,2> parentVec = parent.getUnitVec(m_branchSubUnit);
        Eigen::Vector2d parentVecEig;
        parentVecEig << parentVec[0], parentVec[1];
        Eigen::Vector2d branchVec = rotMatrix*parentVecEig;
        m_subunit_unitVecs.push_back({branchVec(0),branchVec(1)});
        m_subunit_unitVecs.push_back({branchVec(0),branchVec(1)});
        m_subunit_unitVecs.push_back({branchVec(0),branchVec(1)});

        double pointedEndX = parent.getPoint(m_branchSubUnit)[0]
                             + (m_branchSubLen * parentVec[0]);

        double pointedEndY = parent.getPoint(m_branchSubUnit)[1]
                             + (m_branchSubLen * parentVec[1]);

        m_points = {{pointedEndX, pointedEndY, 0}};

        m_points.push_back({ m_points[0][0] + ((m_length/4) * m_subunit_unitVecs[0][0]), m_points[0][1] + ((m_length/4) * m_subunit_unitVecs[0][1]), 0 });
        m_points.push_back({ m_points[1][0] + ((m_length/2) * m_subunit_unitVecs[1][0]), m_points[1][1] + ((m_length/2) * m_subunit_unitVecs[1][1]), 0 });
        m_points.push_back({ m_points[2][0] + ((m_length/4) * m_subunit_unitVecs[2][0]), m_points[2][1] + ((m_length/4) * m_subunit_unitVecs[2][1]), 0 });

        m_actualSubLengths.push_back( m_length / 4);
        m_actualSubLengths.push_back( m_length / 2);
        m_actualSubLengths.push_back( m_length / 4);

        m_prescribedSubLengths.push_back( m_length / 4);
        m_prescribedSubLengths.push_back( m_length / 2);
        m_prescribedSubLengths.push_back( m_length / 4);

        s_total_subunits += 3;

}

//Constructor for new severed filament
Actin::Actin(int ID, const std::array<double,2> &severPoint, int severMonomer,
             int severSub, double lenAlongSub, double severTime, const Actin &original)
        : Actin() // Calls constructor 0
{
        m_id =  ID;
        m_birthtime = severTime;
        m_structure_id = s_structure_idGenerator++;
        m_subStructure_id = s_subStructure_idGenerator++;
        m_crossStructure_id = s_crossStructure_idGenerator++;
        m_points = {{severPoint[0], severPoint[1], 0}};
        m_num_monomers = original.getNumMonomers() - severMonomer - 1;
        m_length = m_num_monomers*s_monomerLength;

        m_distToLeadCL = m_num_monomers;
        if (m_length < 3*s_segmentationLength)
        {
                // needs 4 points
                // needs to be straight
                // end point same as original
                m_num_subunits = 3;

                for (int i = 0; i < m_num_subunits; ++i)
                {
                        m_points.push_back({0,0,0});
                        m_subunit_unitVecs.push_back({0,0});
                        m_actualSubLengths.push_back(0);
                        m_prescribedSubLengths.push_back(0);
                }

                m_actualSubLengths[0] = m_length/4;
                m_actualSubLengths[1] = m_length/2;
                m_actualSubLengths[2] = m_length/4;

                m_prescribedSubLengths[0] = m_length/4;
                m_prescribedSubLengths[1] = m_length/2;
                m_prescribedSubLengths[2] = m_length/4;

                m_points[3] = original.getBarbedEnd();

                straightenFila();
        }
        else
        {

                // Can bend, needs 5 points or more
                m_flexible = true;
                for (int i = severSub+1; i <= original.getNumSubs(); ++i)
                {
                        m_points.push_back({original.getPoint(i)[0], original.getPoint(i)[1], 0});
                        m_subunit_unitVecs.push_back({0,0});
                        m_actualSubLengths.push_back(0);
                        m_prescribedSubLengths.push_back(original.getPresSubLengths()[i-1]);
                }
                m_prescribedSubLengths[0] -= lenAlongSub;
                m_num_subunits = original.getNumSubs() - severSub;
                if (m_num_subunits == 3)
                {
                    // Only has four points!
                    // Must have one more point
                    // Add the barbed end again

                    //

                    m_points.push_back(original.getBarbedEnd());
                    m_subunit_unitVecs.push_back({0,0});
                    m_actualSubLengths.push_back(0);
                    m_prescribedSubLengths.push_back(0);

                    m_num_subunits = 4;
                    m_actualSubLengths[0] = s_segmentationLength/2;
                    m_actualSubLengths[1] = s_segmentationLength;
                    double end2Subs = m_length - (3*s_segmentationLength)/2;
                    double endSub = (end2Subs-s_segmentationLength/2)/2;
                    double penultSub = end2Subs-endSub;

                    m_actualSubLengths[2] = penultSub;
                    m_actualSubLengths[3] = endSub;

                    m_prescribedSubLengths[0] = s_segmentationLength/2;
                    m_prescribedSubLengths[1] = s_segmentationLength;
                    m_prescribedSubLengths[2] = penultSub;
                    m_prescribedSubLengths[3] = endSub;

                    straightenFila();
                    m_flexible = true;
                    // Best way around this is just to straighten it
                }

                updateUnitVecsNoCheck();
                updateSubLengths();

                for (int i = 0; i < m_num_subunits-3; ++i)
                {
                        m_daughterIDs.push_back(std::vector<int>());
                        m_stericCells.push_back(std::vector<int>());
                }
                assert(m_num_subunits > 3);
        }


        for (int i = severMonomer+1; i < original.getNumMonomers(); ++i)
        {
                if (i-severMonomer-1 >= s_maxSpacing)
                {
                        m_monosCompat.push_back(original.getMonoCompat(i));
                }
                else
                {
                        m_monosCompat.push_back(false);
                }

                m_monoBirthTime.push_back(severTime);
        }

        for (unsigned int i = 0; i < original.getAvailMonoVec().size(); ++i)
        {
                if (original.getAvailMonoID(i) > severMonomer && original.getAvailMonoID(i)-severMonomer-1 >= s_maxSpacing)
                {
                        // Now belongs to new filament
                        m_availMonos.push_back(original.getAvailMonoID(i)-severMonomer-1);
                }
        }
        m_parent_id = -1; // invalid parent id means no parent
        m_regionID = original.getRegionID();
        m_capped_p = false;
        m_capped_b = original.getBarbedCapped();

        assert(m_num_subunits >= 3);


}

Actin::~Actin(){
        // Destructor
        //std::cout << "Actin filament correctly deleted" << std::endl;
}

int Actin::branchDir()
{
        // For 2D branching, function that determines whether branch is clockwise
        // returning -1 or anti-clockwise returning +1

        return (rng.m_probDist(rng.m_mersenne) < 0.5) ? 1 : -1;
}

void Actin::polymerise_barbed(std::vector<Actin> &actinVec, const bool steric,
                              StericGrid &stericGrid, const bool bDynamics,
                              const bool tether, const double currTime,
                              std::normal_distribution<double> tetherDistri,
                              std::vector<MembraneWall> &memWalls,
                              bool moveWall, double d_i,
                              std::vector<ProteinRegion> &branchRegions,
                              std::vector<ProteinRegion> &nucRegions,
                              std::vector<ProteinRegion> &capRegions,
                              std::vector<ProteinRegion> &sevRegions)
{
        // Function that polymerises the filament, increasing the length by the
        // monomer length AND adjusting the end position accordingly.
        m_length += s_monomerLength;


        if (bDynamics && (m_length - s_monomerLength) <= 3*s_segmentationLength && m_length > 3*s_segmentationLength)
        {
          // if bdynamics is on, its now flexible

            m_flexible = true;
        }

        if ((m_length - s_monomerLength) <= 2*s_segmentationLength && m_length > 2*s_segmentationLength)
        {
            if (m_parent_id != -1)
            {
              // if its a branch it now belongs to a new substructure
              m_subStructure_id = ++s_subStructure_idGenerator;
              changeSubStructure(actinVec);
            }

            if (m_cLinkMaster)
            {
                // An Other filament should become master, if all other filaments
                // are flexible then there is no master
                //std::cout << "switch master?" << std::endl;
                m_cLinkMaster = false;
                findNewMaster(actinVec);
            }
        }

        // New scheme

        m_points[m_num_subunits][0] += (s_monomerLength * m_subunit_unitVecs[m_num_subunits-1][0]);
        m_points[m_num_subunits][1] += (s_monomerLength * m_subunit_unitVecs[m_num_subunits-1][1]);

        if (!m_flexible)
        {
                m_points[m_num_subunits-1][0] += ((3./4)*s_monomerLength) * m_subunit_unitVecs[m_num_subunits-2][0];
                m_points[m_num_subunits-1][1] += ((3./4)*s_monomerLength) * m_subunit_unitVecs[m_num_subunits-2][1];

                m_points[m_num_subunits-2][0] += ((1./4)*s_monomerLength) * m_subunit_unitVecs[m_num_subunits-3][0];
                m_points[m_num_subunits-2][1] += ((1./4)*s_monomerLength) * m_subunit_unitVecs[m_num_subunits-3][1];

                m_actualSubLengths[m_num_subunits-1] += s_monomerLength/4;
                m_prescribedSubLengths[m_num_subunits-1] += s_monomerLength/4;

                m_actualSubLengths[m_num_subunits-2] += s_monomerLength/2;
                m_prescribedSubLengths[m_num_subunits-2] += s_monomerLength/2;

                m_actualSubLengths[m_num_subunits-3] += s_monomerLength/4;
                m_prescribedSubLengths[m_num_subunits-3] += s_monomerLength/4;

        }
        else
        {
                m_points[m_num_subunits-1][0] += ((s_monomerLength/2) * m_subunit_unitVecs[m_num_subunits-2][0]);
                m_points[m_num_subunits-1][1] += ((s_monomerLength/2) * m_subunit_unitVecs[m_num_subunits-2][1]);

                m_actualSubLengths[m_num_subunits-1] += s_monomerLength/2;
                m_prescribedSubLengths[m_num_subunits-1] += s_monomerLength/2;
                m_actualSubLengths[m_num_subunits-2] += s_monomerLength/2;
                m_prescribedSubLengths[m_num_subunits-2] += s_monomerLength/2;
        }


        ++m_num_monomers;

        ++m_distToLeadBR;
        ++m_distToLeadCL;
        if (m_distToLeadBR > m_distToLeadCL)
        {
            // CL is closest to barbed end, not branch
            if (m_distToLeadCL > s_cLSpacing)
            {
                m_availMonos.push_back(m_num_monomers-1);
                m_monosCompat.push_back(true);
            }
            else
            {
                m_monosCompat.push_back(false);

            }
        }
        else if (m_distToLeadCL > m_distToLeadBR)
        {
            // Br is closest
            if (m_distToLeadBR > s_branchSpacing)
            {
                m_availMonos.push_back(m_num_monomers-1);
                m_monosCompat.push_back(true);
            }
            else
            {
                m_monosCompat.push_back(false);
            }
        }
        else
        {
            //filament has no branches and no crosslinks
            if (m_num_monomers <= s_maxSpacing)
            {
                m_monosCompat.push_back(false);
            }
            else
            {
                m_monosCompat.push_back(true);
                m_availMonos.push_back(m_num_monomers-1);
            }

        }

        // Need to update points vector
        m_changed = true;
        addPointsBarbed(actinVec, steric, stericGrid);

        if (!m_flexible)
        {

                moveBranchSub_down(actinVec);
                moveTetherSub_down();
                moveClSub_down(actinVec);
        }
        else
        {

                moveBranchSub_poly_barb(actinVec);
                moveTetherPoint_poly_barb();
                moveClPoint_poly_barb(actinVec);
        }

        if (moveWall)
        {

            memWalls[0].moveWall(d_i, stericGrid);
            memWalls[0].moveCoupledRegions(0, d_i, branchRegions,
                                           nucRegions, capRegions,
                                           sevRegions);
        }

        m_monoBirthTime.push_back(currTime);

        if (tether)
        {
            m_distToTetBarb += 1;

            if (getNumTethers() == 0)
                m_distToTetPoint += 1;


            if (m_distToTetBarb >= m_preDetDistToTetBarb)
            {
                if (m_parent_id == -1 || getLength() >= 2*s_segmentationLength)
                {
                    // is not a "short" branch, and so can tether
                    addTetherPointBarb();
                    // set the next length
                    int predet = round(tetherDistri(rng.m_mersenne)/Actin::s_monomerLength);

                    m_preDetDistToTetBarb = predet;
                }
            }
        }

}

void Actin::polymerise_pointed(std::vector<Actin> &actinVec, const bool steric,
                               StericGrid &stericGrid, const bool bDynamics,
                               const bool tether, const double currTime,
                               std::normal_distribution<double> tetherDistri,
                               std::vector<MembraneWall> &memWalls,
                               bool moveWall, double d_i,
                               std::vector<ProteinRegion> &branchRegions,
                               std::vector<ProteinRegion> &nucRegions,
                               std::vector<ProteinRegion> &capRegions,
                               std::vector<ProteinRegion> &sevRegions)
{
        // Function that polymerises the filament AT THE POINTED END,
        //increasing the length by the
        // monomer length AND adjusting the end position accordingly.
        m_length += s_monomerLength;
        if (bDynamics && (m_length - s_monomerLength) <= 3*s_segmentationLength && m_length > 3*s_segmentationLength)
        {
            // if bdynamics is on, its now flexible
            m_flexible = true;
        }

        if ((m_length - s_monomerLength) <= 2*s_segmentationLength && m_length > 2*s_segmentationLength)
        {
            if (m_parent_id != -1)
            {
              // if its a branch it now belongs to a new substructure
              m_subStructure_id = ++s_subStructure_idGenerator;
              changeSubStructure(actinVec);
            }

            if (m_cLinkMaster)
            {
                // An Other filament should become master, if all other filaments
                // are flexible then there is no master
                m_cLinkMaster = false;
                findNewMaster(actinVec);
            }
        }


        m_points[0][0] -= (s_monomerLength * m_subunit_unitVecs[0][0]);
        m_points[0][1] -= (s_monomerLength * m_subunit_unitVecs[0][1]);

        if (!m_flexible)
        {
                m_points[1][0] -= ((3./4)*s_monomerLength) * m_subunit_unitVecs[1][0];
                m_points[1][1] -= ((3./4)*s_monomerLength) * m_subunit_unitVecs[1][1];

                m_points[2][0] -= ((1./4)*s_monomerLength) * m_subunit_unitVecs[2][0];
                m_points[2][1] -= ((1./4)*s_monomerLength) * m_subunit_unitVecs[2][1];

                m_actualSubLengths[0] += s_monomerLength/4;
                m_prescribedSubLengths[0] += s_monomerLength/4;

                m_actualSubLengths[1] += s_monomerLength/2;
                m_prescribedSubLengths[1] += s_monomerLength/2;

                m_actualSubLengths[2] += s_monomerLength/4;
                m_prescribedSubLengths[2] += s_monomerLength/4;
        }
        else
        {
                m_points[1][0] -= ((s_monomerLength/2) * m_subunit_unitVecs[1][0]);
                m_points[1][1] -= ((s_monomerLength/2) * m_subunit_unitVecs[1][1]);

                m_actualSubLengths[0] += s_monomerLength/2;
                m_prescribedSubLengths[0] += s_monomerLength/2;
                m_actualSubLengths[1] += s_monomerLength/2;
                m_prescribedSubLengths[1] += s_monomerLength/2;
        }


        ++m_num_monomers;
        shiftMonosCompatUp();
        shiftAvailMonosUp();

        shiftMotherMonosUp(actinVec);
        shiftTetherMonosUp();
        shiftClMonosUp(actinVec);

        // If there are no daughters currently then the distToBR needs to be incremented
        if (m_daughter_num == 0)
        {
            ++m_distToLeadBR;
        }

        if (getNumCLinks() == 0)
        {
            ++m_distToLeadCL;
        }

        // Do we need to add an available monomer?
        // Need to find out if m_monosCompat[s_maxSpacing] should be true or false
        if (m_num_monomers > s_maxSpacing)
        {
            if (m_daughter_num == 0 && getNumCLinks() == 0)
            {
                // Nothing on the fila so its available (true)
                m_monosCompat[s_maxSpacing] = true;
                addAvailMono();
            }
            else
            {
                std::vector<int> unAvailMonos = getUnavailMonos(actinVec);
                std::sort(unAvailMonos.begin(), unAvailMonos.end());

                if (std::find(unAvailMonos.begin(), unAvailMonos.end(), s_maxSpacing) == unAvailMonos.end())
                {
                    // Should be available
                    m_monosCompat[s_maxSpacing] = true;
                    addAvailMono();
                }
            }
        }

        if (m_availMonos.size() > 0)
        {
                assert(m_availMonos[m_availMonos.size()-1] < m_num_monomers);
        }

        shiftMonosBirthTimeUp();
        m_monoBirthTime[0] = currTime;

        // Need to update points vector
        m_changed = true;
        addPointsPointed(actinVec, steric, stericGrid);
        if (!m_flexible)
        {
                moveBranchSub_up(actinVec);
                moveTetherSub_up();
                moveClSub_up(actinVec);
        }
        else
        {
                moveBranchSub_poly_point(actinVec);
                moveTetherPoint_poly_point();
                moveClPoint_poly_point(actinVec);
        }

        if (moveWall)
        {
            memWalls[0].moveWall(d_i, stericGrid);
            memWalls[0].moveCoupledRegions(0,d_i, branchRegions,
                                            nucRegions, capRegions,
                                            sevRegions);
        }

        if (tether)
        {
            m_distToTetPoint += 1;

            if (getNumTethers() == 0)
                m_distToTetBarb += 1;

            if (m_distToTetPoint >= m_preDetDistToTetPoint)
            {
                addTetherPointPointed();
                // set the next length
                int predet = round(tetherDistri(rng.m_mersenne)/Actin::s_monomerLength);
                m_preDetDistToTetPoint = predet;
            }
        }

}

void Actin::depolymerise_barbed(std::vector<Actin> &actinVec, int &nActin,
                                bool arpPool, ArpGrid &arpGrid,
                                const bool steric, StericGrid &stericGrid,
                                const std::vector<ExcZone> &excZones,
                                std::vector<MembraneWall> &memWalls,
                                std::vector<Membrane> &membranes,
                                Cortex &cortex,
                                const bool bDynamics, const bool tether,
                                bool moveWall, double d_i,
                                std::vector<ProteinRegion> &branchRegions,
                                std::vector<ProteinRegion> &nucRegions,
                                std::vector<ProteinRegion> &capRegions,
                                std::vector<ProteinRegion> &sevRegions)
{
        // Function that removes one monomer from the barbed end
        m_length -= s_monomerLength;

        m_points[m_num_subunits][0] -= (s_monomerLength * m_subunit_unitVecs[m_num_subunits-1][0]);
        m_points[m_num_subunits][1] -= (s_monomerLength * m_subunit_unitVecs[m_num_subunits-1][1]);

        if (!m_flexible)
        {
                m_points[m_num_subunits-1][0] -= ((3./4)*s_monomerLength) * m_subunit_unitVecs[m_num_subunits-2][0];
                m_points[m_num_subunits-1][1] -= ((3./4)*s_monomerLength) * m_subunit_unitVecs[m_num_subunits-2][1];

                m_points[m_num_subunits-2][0] -= ((1./4)*s_monomerLength) * m_subunit_unitVecs[m_num_subunits-3][0];
                m_points[m_num_subunits-2][1] -= ((1./4)*s_monomerLength) * m_subunit_unitVecs[m_num_subunits-3][1];

                m_actualSubLengths[m_num_subunits-1] -= s_monomerLength/4;
                m_prescribedSubLengths[m_num_subunits-1] -= s_monomerLength/4;

                m_actualSubLengths[m_num_subunits-2] -= s_monomerLength/2;
                m_prescribedSubLengths[m_num_subunits-2] -= s_monomerLength/2;

                m_actualSubLengths[m_num_subunits-3] -= s_monomerLength/4;
                m_prescribedSubLengths[m_num_subunits-3] -= s_monomerLength/4;
        }
        else
        {
                m_points[m_num_subunits-1][0] -= ((s_monomerLength/2) * m_subunit_unitVecs[m_num_subunits-2][0]);
                m_points[m_num_subunits-1][1] -= ((s_monomerLength/2) * m_subunit_unitVecs[m_num_subunits-2][1]);

                m_actualSubLengths[m_num_subunits-1] -= s_monomerLength/2;
                m_prescribedSubLengths[m_num_subunits-1] -= s_monomerLength/2;
                m_actualSubLengths[m_num_subunits-2] -= s_monomerLength/2;
                m_prescribedSubLengths[m_num_subunits-2] -= s_monomerLength/2;
        }


        if ((m_length+s_monomerLength >= 3*m_des_subunit_len) && (m_length < 3*m_des_subunit_len))
        {
                bool straightFail = straightenFila(actinVec, excZones, memWalls,
                                                   membranes, cortex, steric,
                                                   stericGrid, 1);
                if (straightFail)
                {
                        return;
                }
        }

        if ((m_length+s_monomerLength >= 2*m_des_subunit_len) && (m_length < 2*m_des_subunit_len))
        {
            if (bDynamics)
            {
                if (m_parent_id != -1)
                {
                  // if its a branch it now belongs to a the substructure of it's parent
                  m_subStructure_id = actinVec[m_parent_id].getSubStructureID();
                  changeSubStructure(actinVec);
                }


                if (getNumCLinks() > 0)
                {
                    // Check to see if at least one in the crosslink STRUCTURE is still a master
                    // if not set this one to be master
                    if (!checkForMasterInStructure(actinVec))
                    {
                        // No master in structure, set it to this one
                        m_cLinkMaster = true;
                    }
                }
            }
        }

        --m_num_monomers;
        --m_distToLeadBR;
        --m_distToLeadCL;

        m_monosCompat.pop_back();


        if (m_availMonos.size() > 0 && m_availMonos[m_availMonos.size()-1] == m_num_monomers)
        {
                m_availMonos.pop_back();
        }

        // Need to check if remove branches here!
        if (m_daughter_num != 0)
        {
                autoDebranchBarb(*this, nActin, actinVec, arpPool,
                                 arpGrid, steric, stericGrid);
        }

        if (getNumCLinks() != 0)
        {
            // auto unLinkBarb
            autoUnlinkBarb(*this, actinVec);
        }

        m_changed = true;
        removePointsBarbed(actinVec, steric, stericGrid);

        if (!m_flexible)
        {
            moveBranchSub_up(actinVec);
            moveTetherSub_up();
            moveClSub_up(actinVec);
        }
        else
        {
            moveBranchSub_depoly_barb(actinVec);
            moveTetherPoint_depoly_barb();
            moveClPoint_depoly_barb(actinVec);
        }

        m_monoBirthTime.pop_back();

        if (moveWall)
        {
            memWalls[0].moveWall(d_i, stericGrid);
            memWalls[0].moveCoupledRegions(0, d_i, branchRegions,
                                           nucRegions, capRegions,
                                           sevRegions);
        }

        if (tether)
        {
            m_distToTetBarb -= 1;
            // If this has gone negative need to remove the tether

            if (getNumTethers() > 0)
            {
                if (m_distToTetBarb < 0)
                    remTetherB(actinVec);
            }
            else
            {
                m_distToTetPoint -= 1;
            }

        }


}

void Actin::cancelBarbDepol(std::vector<Actin> &actinVec, const bool steric,
                            StericGrid &stericGrid)
{
        // Dealing with the rare event that if a filament straightening after depol,
        // results in steric hindrance breaking, avoid depoly

        m_length += s_monomerLength;


        m_points[m_num_subunits][0] += (s_monomerLength * m_subunit_unitVecs[m_num_subunits-1][0]);
        m_points[m_num_subunits][1] += (s_monomerLength * m_subunit_unitVecs[m_num_subunits-1][1]);
        m_points[m_num_subunits-1][0] += ((s_monomerLength/2) * m_subunit_unitVecs[m_num_subunits-2][0]);
        m_points[m_num_subunits-1][1] += ((s_monomerLength/2) * m_subunit_unitVecs[m_num_subunits-2][1]);

        m_actualSubLengths[m_num_subunits-1] += s_monomerLength/2;
        m_prescribedSubLengths[m_num_subunits-1] += s_monomerLength/2;
        m_actualSubLengths[m_num_subunits-2] += s_monomerLength/2;
        m_prescribedSubLengths[m_num_subunits-2] += s_monomerLength/2;


}

void Actin::cancelPointDepol(std::vector<Actin> &actinVec, const bool steric,
                             StericGrid &stericGrid)
{
        // Dealing with the rare event that if a filament straightening after depol,
        // results in steric hindrance breaking, avoid depoly

        m_length += s_monomerLength;


        m_points[0][0] -= (s_monomerLength * m_subunit_unitVecs[0][0]);
        m_points[0][1] -= (s_monomerLength * m_subunit_unitVecs[0][1]);
        m_points[1][0] -= ((s_monomerLength/2) * m_subunit_unitVecs[1][0]);
        m_points[1][1] -= ((s_monomerLength/2) * m_subunit_unitVecs[1][1]);

        m_actualSubLengths[0] += s_monomerLength/2;
        m_prescribedSubLengths[0] += s_monomerLength/2;
        m_actualSubLengths[1] += s_monomerLength/2;
        m_prescribedSubLengths[1] += s_monomerLength/2;


}

void Actin::depolymerise_pointed(std::vector<Actin> &actinVec, int &nActin,
                                 bool arpPool, ArpGrid &arpGrid,
                                 const bool steric, StericGrid &stericGrid,
                                 const std::vector<ExcZone> &excZones,
                                 std::vector<MembraneWall> &memWalls,
                                 std::vector<Membrane> &membranes,
                                 Cortex &cortex,
                                 const bool bDynamics,
                                 const bool tether,
                                 bool moveWall, double d_i,
                                 std::vector<ProteinRegion> &branchRegions,
                                 std::vector<ProteinRegion> &nucRegions,
                                 std::vector<ProteinRegion> &capRegions,
                                 std::vector<ProteinRegion> &sevRegions)
{
        // Function that removes one monomer from the pointed end

        m_length -= s_monomerLength;


        m_points[0][0] += (s_monomerLength * m_subunit_unitVecs[0][0]);
        m_points[0][1] += (s_monomerLength * m_subunit_unitVecs[0][1]);

        if (!m_flexible)
        {
                m_points[1][0] += ((3./4)*s_monomerLength) * m_subunit_unitVecs[1][0];
                m_points[1][1] += ((3./4)*s_monomerLength) * m_subunit_unitVecs[1][1];

                m_points[2][0] += ((1./4)*s_monomerLength) * m_subunit_unitVecs[2][0];
                m_points[2][1] += ((1./4)*s_monomerLength) * m_subunit_unitVecs[2][1];


                m_actualSubLengths[0] -= s_monomerLength/4;
                m_prescribedSubLengths[0] -= s_monomerLength/4;

                m_actualSubLengths[1] -= s_monomerLength/2;
                m_prescribedSubLengths[1] -= s_monomerLength/2;

                m_actualSubLengths[2] -= s_monomerLength/4;
                m_prescribedSubLengths[2] -= s_monomerLength/4;
        }
        else
        {
                m_points[1][0] += ((s_monomerLength/2) * m_subunit_unitVecs[1][0]);
                m_points[1][1] += ((s_monomerLength/2) * m_subunit_unitVecs[1][1]);

                m_actualSubLengths[0] -= s_monomerLength/2;
                m_prescribedSubLengths[0] -= s_monomerLength/2;
                m_actualSubLengths[1] -= s_monomerLength/2;
                m_prescribedSubLengths[1] -= s_monomerLength/2;
        }

        if ((m_length+s_monomerLength >= 3*m_des_subunit_len) && (m_length < 3*m_des_subunit_len))
        {
                std::cout << "Straightening filament from pointed end" << std::endl;
                bool straightFail = straightenFila(actinVec, excZones, memWalls,
                                                   membranes, cortex, steric,
                                                   stericGrid, 0);
                if (straightFail)
                {
                        return;
                }
        }

        if ((m_length+s_monomerLength >= 2*m_des_subunit_len) && (m_length < 2*m_des_subunit_len))
        {
            if (bDynamics)
            {
              // if bdynamics is on, its now rigid
                if (m_parent_id != -1)
                {
                  // if its a branch it now belongs to a the substructure of its parent
                  m_subStructure_id = actinVec[m_parent_id].getSubStructureID();
                  changeSubStructure(actinVec);
                }

                if (getNumCLinks() > 0)
                {
                    // Check to see if at least one in the crosslink STRUCTURE is still a master
                    // if not set this one to be master
                    if (!checkForMasterInStructure(actinVec))
                    {
                        // No master in structure, set it to this one
                        m_cLinkMaster = true;
                    }
                }
            }
        }

        --m_num_monomers;

        shiftMotherMonosDown(actinVec);
        shiftTetherMonosDown();
        shiftClMonosDown(actinVec);

        shiftMonosCompatDown();
        if (m_num_monomers >= s_maxSpacing)
        {
            if (s_maxSpacing > 0)
            {
                m_monosCompat[s_maxSpacing-1] = false; // must always be false
            }
        }

        shiftAvailMonosDown();
        if (m_daughter_num == 0)
        {
            --m_distToLeadBR;
        }
        if (getNumCLinks() == 0)
        {
            --m_distToLeadCL;
        }
        // NEed to check if remove branches here!
        if (m_daughter_num != 0)
        {
                autoDebranchPoint(*this, nActin, actinVec, arpPool,
                                  arpGrid, steric, stericGrid);
        }

        if (getNumCLinks() != 0)
        {
            autoUnlinkPoint(*this, actinVec);
        }

        m_changed = true;

        removePointsPointed(actinVec, steric, stericGrid);

        if (!m_flexible)
        {
                moveBranchSub_down(actinVec);
                moveTetherSub_down();
                moveClSub_down(actinVec);
        }
        else
        {
                moveBranchSub_depoly_point(actinVec);
                moveTetherPoint_depoly_point();
                moveClPoint_depoly_point(actinVec);
        }

        shiftMonosBirthTimeDown();

        if (moveWall)
        {
            memWalls[0].moveWall(d_i, stericGrid);
            memWalls[0].moveCoupledRegions(0, d_i, branchRegions,
                                           nucRegions, capRegions,
                                           sevRegions);
        }

        if (tether)
        {
            m_distToTetPoint -= 1;
            if (getNumTethers() > 0)
            {
                // If this has gone negative need to remove the tether
                if (m_distToTetPoint < 0)
                    remTetherP(actinVec);
            }
            else
            {
                m_distToTetBarb -= 1;
            }
        }


}

void Actin::addPointsBarbed(std::vector<Actin> &actinVec, const bool steric,
                            StericGrid &stericGrid)
{
        // Function that checks whether a new subunit needs to be added at the
        // barbed end

        // New point is added to the middle of the barbed end subunit!
        if (m_length < 3*m_des_subunit_len)
        {
            return;
        }

        if (m_num_subunits == 3)
        {

              // Need an additional point
              ++s_total_subunits;
              ++m_num_subunits;

              // "Add" the endpoint to the vector
              std::array<double,3> tmp { };
              tmp[0] = m_points[m_num_subunits-1][0];
              tmp[1] = m_points[m_num_subunits-1][1];

              m_points.push_back(tmp);
              // the subunit lengths will be adjusted in the function call at the end
              // of this function

              m_actualSubLengths.push_back(0);
              m_prescribedSubLengths.push_back(0);
              m_subunit_unitVecs.push_back({0,0});

              m_daughterIDs.push_back(std::vector<int>());

              if (steric)
                      m_stericCells.push_back(std::vector<int>());



              // Readjust all the points
              m_points[1][0] = m_points[0][0] + ((m_des_subunit_len/2)*m_subunit_unitVecs[0][0]); // x
              m_points[1][1] = m_points[0][1] + ((m_des_subunit_len/2)*m_subunit_unitVecs[0][1]); // y

              m_points[2][0] = m_points[1][0] + (m_des_subunit_len*m_subunit_unitVecs[1][0]); // x
              m_points[2][1] = m_points[1][1] + (m_des_subunit_len*m_subunit_unitVecs[1][1]); // y

              double lengthToEnd = m_length - m_des_subunit_len*(3./2);
              lengthToEnd -= m_des_subunit_len/2;
              lengthToEnd /= 2;
              double endsubL = lengthToEnd;
              lengthToEnd += m_des_subunit_len/2;

              m_points[3][0] = m_points[2][0] + (lengthToEnd*m_subunit_unitVecs[2][0]); // x
              m_points[3][1] = m_points[2][1] + (lengthToEnd*m_subunit_unitVecs[2][1]); // y

              m_points[4][0] = m_points[3][0] + (endsubL*m_subunit_unitVecs[2][0]); // x
              m_points[4][1] = m_points[3][1] + (endsubL*m_subunit_unitVecs[2][1]); // y

              m_actualSubLengths[0] = m_des_subunit_len/2;
              m_actualSubLengths[1] = m_des_subunit_len;
              m_actualSubLengths[2] = lengthToEnd;
              m_actualSubLengths[3] = endsubL;

              m_prescribedSubLengths[0] = m_des_subunit_len/2;
              m_prescribedSubLengths[1] = m_des_subunit_len;
              m_prescribedSubLengths[2] = lengthToEnd;
              m_prescribedSubLengths[3] = endsubL;


              updateUnitVecsNoCheck();


              for (int i = 0; i < m_num_subunits-1; ++i)
              {
                  checkAndMoveBranchUp(actinVec, i);
                  checkAndMoveClUp(actinVec, i);
                  checkAndMoveTetherUp(0, 0, true);
              }

              if (steric)
                      stericGrid.updateCellsAll(*this);
        }
        else if (m_prescribedSubLengths[m_num_subunits-2] + m_prescribedSubLengths[m_num_subunits-1] >= 2.5*m_des_subunit_len)
        {
                ++s_total_subunits;
                ++m_num_subunits;

                // "Add" the endpoint to the vector
                std::array<double,3> tmp { };
                tmp[0] = m_points[m_num_subunits-1][0];
                tmp[1] = m_points[m_num_subunits-1][1];

                m_points.push_back(tmp);

                // the subunit lengths will be adjusted in the function call at the end
                // of this function
                m_actualSubLengths.push_back(0);
                m_prescribedSubLengths.push_back(0);
                m_subunit_unitVecs.push_back({0,0});

                m_daughterIDs.push_back(std::vector<int>());
                if (steric)
                        m_stericCells.push_back(std::vector<int>());
                // however branches from end 2 subs may have moved up 1, this needs to be
                // checked

                // readjust the now-p-pemultimate point,
                m_points[m_num_subunits-2][0] = m_points[m_num_subunits-3][0] + (m_des_subunit_len*m_subunit_unitVecs[m_num_subunits-3][0]); // x
                m_points[m_num_subunits-2][1] = m_points[m_num_subunits-3][1] + (m_des_subunit_len*m_subunit_unitVecs[m_num_subunits-3][1]); // y

                std::array<double,2> endpoint { m_points[m_num_subunits][0], m_points[m_num_subunits][1] };
                std::array<double,2> ppenultimate { m_points[m_num_subunits-2][0], m_points[m_num_subunits-2][1] };

                double lengthToEnd = sqrt(distanceBetPoints2DSQR(endpoint, ppenultimate));

                lengthToEnd -= m_des_subunit_len/2;
                lengthToEnd /= 2;
                lengthToEnd += m_des_subunit_len/2;



                m_points[m_num_subunits-1][0] = m_points[m_num_subunits-2][0] + (lengthToEnd*m_subunit_unitVecs[m_num_subunits-2][0]); // x
                m_points[m_num_subunits-1][1] = m_points[m_num_subunits-2][1] + (lengthToEnd*m_subunit_unitVecs[m_num_subunits-2][1]); // y



                updateSubLengths();
                // Now do the same for prescribed subunit lengths

                lengthToEnd = m_prescribedSubLengths[m_num_subunits-2] + m_prescribedSubLengths[m_num_subunits-3] - m_des_subunit_len;

                double endsubL;

                lengthToEnd -= m_des_subunit_len/2;
                lengthToEnd /= 2;
                endsubL = lengthToEnd;
                if (endsubL < 0.5*m_des_subunit_len)
                {
                        std::cout << "Error in adding subunit, line 920 of Actin.cpp" << std::endl;
                        exit(1);
                }
                lengthToEnd += m_des_subunit_len/2;

                m_prescribedSubLengths[m_num_subunits-3] = m_des_subunit_len;
                m_prescribedSubLengths[m_num_subunits-2] = lengthToEnd;
                m_prescribedSubLengths[m_num_subunits-1] = endsubL;

                // Check for moving branches

                updateUnitVecsNoCheck();


                moveBranchSub_Add_barb(actinVec);
                moveTetherPoint_Add_barb();
                moveClPoint_Add_barb(actinVec);

                if (steric)
                        stericGrid.updateCellsAddPointBarb(*this);
        }
}

void Actin::addPointsPointed(std::vector<Actin> &actinVec, const bool steric,
                             StericGrid &stericGrid)
{
        // As above but for pointed end, need to shift some vectors around here
        if (m_length < 3*m_des_subunit_len)
        {
            return;
        }

        if (m_num_subunits == 3)
        {
              // Need an additional point
              ++s_total_subunits;
              ++m_num_subunits;

              // "Add" the endpoint to the vector
              m_points.push_back( m_points[m_num_subunits-1]);
              for (int i = m_num_subunits-1; i > 1; --i)
              {
                      // moves everything but point 0 up
                      m_points[i] = m_points[i-1];
              }

              m_actualSubLengths.push_back(0);
              m_prescribedSubLengths.push_back(0);
              m_subunit_unitVecs.push_back({0,0});
              for (int i = m_subunit_unitVecs.size()-1; i > 0; --i)
              {
                      // moves everything up one
                      m_subunit_unitVecs[i] = m_subunit_unitVecs[i-1];
              }
              m_subunit_unitVecs[0] = {0,0};

              m_daughterIDs.push_back(std::vector<int>());
              if (steric)
              {
                  m_stericCells.push_back(std::vector<int>());
                  stericGrid.moveCellsUp(*this);
              }

              for (int i = m_num_subunits-1; i > 0; --i)
              {
                      for (unsigned int j = 0; j < m_daughterIDs[i-1].size(); ++j)
                      {
                              int k = m_daughterIDs[i-1][j];
                              actinVec[k].incrementBRsub();
                      }
                      m_daughterIDs[i] = m_daughterIDs[i-1];
                      if (steric)
                              m_stericCells[i] = m_stericCells[i-1];
              }
              m_daughterIDs[0].clear();
              if (steric)
                      m_stericCells[0].clear();

              // Move any tether points up
              int numTethers = m_tetherMono.size();
              for (int i = 0; i < numTethers; ++i)
              {
                      m_tetherSubunits[i] += 1;
              }

              // Move any crosslink points up too
              int numCLinks = getNumCLinks();
              for (int i = 0; i < numCLinks; ++i)
              {
                  int oldCLSub = m_cLinkActinAndSites[i][1];
                  setClSub(i, oldCLSub+1, actinVec);
              }

              // Readjust all the points
              m_points[3][0] = m_points[4][0] - ((m_des_subunit_len/2)*m_subunit_unitVecs[3][0]); // x
              m_points[3][1] = m_points[4][1] - ((m_des_subunit_len/2)*m_subunit_unitVecs[3][1]); // y

              m_points[2][0] = m_points[3][0] - (m_des_subunit_len*m_subunit_unitVecs[2][0]); // x
              m_points[2][1] = m_points[3][1] - (m_des_subunit_len*m_subunit_unitVecs[2][1]); // y

              double lengthToEnd = m_length - m_des_subunit_len*(3./2);
              lengthToEnd -= m_des_subunit_len/2;
              lengthToEnd /= 2;
              double endsubL = lengthToEnd;
              lengthToEnd += m_des_subunit_len/2;


              m_points[1][0] = m_points[2][0] - (lengthToEnd*m_subunit_unitVecs[1][0]); // x
              m_points[1][1] = m_points[2][1] - (lengthToEnd*m_subunit_unitVecs[1][1]); // y

              m_points[0][0] = m_points[1][0] - (endsubL*m_subunit_unitVecs[1][0]); // x
              m_points[0][1] = m_points[1][1] - (endsubL*m_subunit_unitVecs[1][1]); // y


              m_actualSubLengths[0] = endsubL;
              m_actualSubLengths[1] = lengthToEnd;
              m_actualSubLengths[2] = m_des_subunit_len;
              m_actualSubLengths[3] = m_des_subunit_len/2;

              m_prescribedSubLengths[0] = endsubL;
              m_prescribedSubLengths[1] = lengthToEnd;
              m_prescribedSubLengths[2] = m_des_subunit_len;
              m_prescribedSubLengths[3] = m_des_subunit_len/2;

              updateUnitVecsNoCheck();

              for (int i = 1; i < m_num_subunits; ++i)
              {
                  checkAndMoveBranchDown(actinVec, i);
                  checkAndMoveClDown(actinVec, i);
                  checkAndMoveTetherDown(0, 0, true);
              }

              if (steric)
                      stericGrid.updateCellsAll(*this);


        }
        else if (m_prescribedSubLengths[1]+m_prescribedSubLengths[0] >= 2.5*m_des_subunit_len)
        {

                ++s_total_subunits;
                ++m_num_subunits;

                // "Add" the endpoint to the vector

                m_points.push_back( m_points[m_num_subunits-1]);
                for (int i = m_num_subunits-1; i > 1; --i)
                {
                        // moves everything but point 0 up
                        m_points[i] = m_points[i-1];
                }

                // the subunit lengths will be adjusted in the function call at the end
                // of this function
                m_actualSubLengths.push_back(0);
                m_prescribedSubLengths.push_back(0);
                m_subunit_unitVecs.push_back({0,0});
                for (int i = m_subunit_unitVecs.size()-1; i > 0; --i)
                {
                        // moves everything up one
                        m_subunit_unitVecs[i] = m_subunit_unitVecs[i-1];
                }
                m_subunit_unitVecs[0] = {0,0};

                m_daughterIDs.push_back(std::vector<int>());
                if (steric)
                        m_stericCells.push_back(std::vector<int>());

                if (steric)
                        stericGrid.moveCellsUp(*this);

                for (int i = m_num_subunits-1; i > 0; --i)
                {
                        for (unsigned int j = 0; j < m_daughterIDs[i-1].size(); ++j)
                        {
                                int k = m_daughterIDs[i-1][j];
                                actinVec[k].incrementBRsub();
                        }
                        m_daughterIDs[i] = m_daughterIDs[i-1];
                        if (steric)
                                m_stericCells[i] = m_stericCells[i-1];
                }
                m_daughterIDs[0].clear();
                if (steric)
                        m_stericCells[0].clear();

                // Move any tether points up
                int numTethers = m_tetherMono.size();
                for (int i = 0; i < numTethers; ++i)
                {
                        m_tetherSubunits[i] += 1;
                }

                // Move any crosslink points up too
                int numCLinks = getNumCLinks();
                for (int i = 0; i < numCLinks; ++i)
                {
                    int oldCLSub = m_cLinkActinAndSites[i][1];
                    setClSub(i, oldCLSub+1, actinVec);
                }

                // readjust the now-p-pemultimate point,

                m_points[2][0] = m_points[3][0] - (m_des_subunit_len*m_subunit_unitVecs[2][0]); // x
                m_points[2][1] = m_points[3][1] - (m_des_subunit_len*m_subunit_unitVecs[2][1]); // y

                std::array<double,2> endpoint { m_points[0][0], m_points[0][1] };
                std::array<double,2> ppenultimate { m_points[2][0], m_points[2][1] };

                double lengthToEnd = sqrt(distanceBetPoints2DSQR(endpoint, ppenultimate));

                lengthToEnd -= m_des_subunit_len/2;
                lengthToEnd /= 2;
                lengthToEnd += m_des_subunit_len/2;


                m_points[1][0] = m_points[2][0] - (lengthToEnd*m_subunit_unitVecs[1][0]); // x
                m_points[1][1] = m_points[2][1] - (lengthToEnd*m_subunit_unitVecs[1][1]); // y

                for (int i = m_num_subunits-1; i > 0; --i)
                {
                        m_prescribedSubLengths[i] = m_prescribedSubLengths[i-1];
                }

                updateSubLengths();
                lengthToEnd = m_prescribedSubLengths[1] + m_prescribedSubLengths[2] - m_des_subunit_len;

                double endsubL;

                lengthToEnd -= m_des_subunit_len/2;
                lengthToEnd /= 2;
                endsubL = lengthToEnd;
                lengthToEnd += m_des_subunit_len/2;



                m_prescribedSubLengths[0] = endsubL;
                m_prescribedSubLengths[1] = lengthToEnd;
                m_prescribedSubLengths[2] = m_des_subunit_len;

                moveBranchSub_Add_point(actinVec);
                moveTetherPoint_Add_point();
                moveClPoint_Add_point(actinVec);
                updateUnitVecsNoCheck();

                if (steric)
                        stericGrid.updateCellsAddPointPoint(*this);

        }

}

void Actin::removePointsBarbed(std::vector<Actin> &actinVec, const bool steric,
                               StericGrid &stericGrid, bool severing)
{
        // Function that is called after barbed end depolymerisation to check if
        // we need to remove the barbed end point
        if (m_num_subunits == 3)
        {
                return;
        }

        if (m_num_subunits == 4 && m_prescribedSubLengths[m_num_subunits-2] + m_prescribedSubLengths[m_num_subunits-1] < 1.5*m_des_subunit_len)
        {
            // Special case
            if (m_length >= 3*m_des_subunit_len)
            {
                // DONT REMOVE A POINT!
                // Just have to rejig
                double lengthofEnds = m_length - m_des_subunit_len;
                lengthofEnds /= 4;
                double endsubL = lengthofEnds;
                lengthofEnds += m_des_subunit_len/2;

                m_prescribedSubLengths[0] = endsubL;
                m_prescribedSubLengths[1] = lengthofEnds;
                m_prescribedSubLengths[2] = lengthofEnds;
                m_prescribedSubLengths[3] = endsubL;

                m_actualSubLengths[0] = endsubL;
                m_actualSubLengths[1] = lengthofEnds;
                m_actualSubLengths[2] = lengthofEnds;
                m_actualSubLengths[3] = endsubL;

                updatePointsUp(2);
                updatePointsDown(2);

                // Moving anything attached
                for (int i = 0; i < m_num_subunits; ++i)
                {
                    checkAndMoveBranchUp(actinVec, i);
                    checkAndMoveClUp(actinVec, i);
                    checkAndMoveTetherUp(0, 0, true);
                }

                if (steric)
                        stericGrid.updateCellsAll(*this);
            }
        }

        // Change the centre first
        if (m_prescribedSubLengths[m_num_subunits-2] + m_prescribedSubLengths[m_num_subunits-1] < 1.5*m_des_subunit_len)
        {
                // change the number of subunits here
                --s_total_subunits;
                --m_num_subunits;

                if (severing)
                {
                    ++s_total_subunits;
                    // Need to reverse the above, if this is called during severing
                }

                m_points[m_num_subunits] = m_points[m_num_subunits+1]; // new endpont equals old one
                m_points.pop_back();
                m_actualSubLengths.pop_back();

                if (steric)
                {
                        stericGrid.resetCells(*this,m_num_subunits); // delete from steric grid too
                        m_stericCells.pop_back();
                }

                // Any daughters on this final subunit needs to move onto the new end sub
                std::vector<int> daughterIDs_copy = m_daughterIDs[m_num_subunits];

                for (unsigned int i = 0; i < daughterIDs_copy.size(); ++i)
                {
                        int j = daughterIDs_copy[i];
                        actinVec[j].decrementBRsub();
                        purgeBranchFromDaughterVec(j);
                        addBranchtoParentVector(m_num_subunits-1,j);
                }

                // Same for tether points
                for (unsigned int i = 0; i < m_tetherMono.size(); ++i)
                {
                        if (m_tetherSubunits[i]==m_num_subunits)
                        {
                                m_tetherSubunits[i] -= 1;
                        }
                }

                // and for crosslinkers
                for (int i = 0; i < getNumCLinks(); ++i)
                {
                      if (m_cLinkActinAndSites[i][1] == m_num_subunits)
                      {
                          setClSub(i, m_cLinkActinAndSites[i][1]-1, actinVec);
                      }
                }

                m_daughterIDs.pop_back();
                m_subunit_unitVecs.pop_back();
                // Need to recalculate the end angle

                double dy = m_points[m_num_subunits][1] - m_points[m_num_subunits-2][1];
                double dx = m_points[m_num_subunits][0] - m_points[m_num_subunits-2][0];
                double mag = sqrt(dx*dx + dy*dy);

                // Should these be indices of -2 and -1 instead of -3 and -2???
                m_subunit_unitVecs[m_num_subunits-2][0] = dx/mag;
                m_subunit_unitVecs[m_num_subunits-2][1] = dy/mag;

                m_subunit_unitVecs[m_num_subunits-1][0] = dx/mag;
                m_subunit_unitVecs[m_num_subunits-1][1] = dy/mag;


                // now we need to rejig penultimate point
                std::array<double,2> endpoint { m_points[m_num_subunits][0], m_points[m_num_subunits][1] };
                std::array<double,2> ppenultimate { m_points[m_num_subunits-2][0], m_points[m_num_subunits-2][1] };
                std::array<double,2> penultimate { m_points[m_num_subunits-1][0], m_points[m_num_subunits-1][1] };
                double lengthToEnd = sqrt(distanceBetPoints2DSQR(endpoint, penultimate)) + sqrt(distanceBetPoints2DSQR(ppenultimate, penultimate));

                double endsubL;

                lengthToEnd -= m_actualSubLengths[m_num_subunits-2]/2;
                lengthToEnd /= 2;
                endsubL = lengthToEnd;
                lengthToEnd += m_actualSubLengths[m_num_subunits-2]/2;

                //lengthToEnd is now the length of the penultimate subunit

                m_points[m_num_subunits-1][0] = m_points[m_num_subunits-2][0] + (lengthToEnd*m_subunit_unitVecs[m_num_subunits-2][0]); // x
                m_points[m_num_subunits-1][1] = m_points[m_num_subunits-2][1] + (lengthToEnd*m_subunit_unitVecs[m_num_subunits-2][1]); // y

                m_points[m_num_subunits][0] = m_points[m_num_subunits-1][0] + (endsubL*m_subunit_unitVecs[m_num_subunits-1][0]); // x
                m_points[m_num_subunits][1] = m_points[m_num_subunits-1][1] + (endsubL*m_subunit_unitVecs[m_num_subunits-1][1]); // y



                // shorten the last subunit by the difference

                // recalculate barbed end based on this

                updateSubLengths();

                lengthToEnd = m_prescribedSubLengths[m_num_subunits-1] + m_prescribedSubLengths[m_num_subunits-2] + m_prescribedSubLengths[m_num_subunits];
                m_prescribedSubLengths.pop_back();


                lengthToEnd -= m_prescribedSubLengths[m_num_subunits-2]/2;
                lengthToEnd /= 2;
                endsubL = lengthToEnd;
                lengthToEnd += m_prescribedSubLengths[m_num_subunits-2]/2;


                m_prescribedSubLengths[m_num_subunits-2] = lengthToEnd;

                m_prescribedSubLengths[m_num_subunits-1] = endsubL;

                moveBranchSub_Rem_barb(actinVec);
                moveTetherPoint_Rem_barb();
                moveClPoint_Rem_barb(actinVec);

                updateTetherLoc();
                updateCrossLinkLoc(actinVec);
                updateUnitVecs();
                if (steric)
                        stericGrid.updateCellsRemPointBarb(*this);
        }

}

void Actin::removePointsPointed(std::vector<Actin> &actinVec, const bool steric,
                                StericGrid &stericGrid, bool severing)
{
        // Function that is called after pointed end depolymerisation to check if
        // we need to remove the pointed end point

        // Don't remove if just 4 points
        if (m_num_subunits == 3)
        {
                return;
        }


        if (m_num_subunits == 4 && m_prescribedSubLengths[0] + m_prescribedSubLengths[1] < 1.5*m_des_subunit_len)
        {
            // Special case going down to three perhaps?
            if (m_length >= 3*m_des_subunit_len)
            {
                // DONT REMOVE A POINT!
                // Just have to rejig
                double lengthofEnds = m_length - m_des_subunit_len;
                lengthofEnds /= 4;
                double endsubL = lengthofEnds;
                lengthofEnds += m_des_subunit_len/2;

                m_prescribedSubLengths[0] = endsubL;
                m_prescribedSubLengths[1] = lengthofEnds;
                m_prescribedSubLengths[2] = lengthofEnds;
                m_prescribedSubLengths[3] = endsubL;

                m_actualSubLengths[0] = endsubL;
                m_actualSubLengths[1] = lengthofEnds;
                m_actualSubLengths[2] = lengthofEnds;
                m_actualSubLengths[3] = endsubL;

                updatePointsUp(2);
                updatePointsDown(2);

                // Moving anything attached
                for (int i = 0; i < m_num_subunits; ++i)
                {
                    checkAndMoveBranchDown(actinVec, i);
                    checkAndMoveClDown(actinVec, i);
                    checkAndMoveTetherDown(0, 0, true);
                }

                if (steric)
                        stericGrid.updateCellsAll(*this);
            }
        }

        if (m_prescribedSubLengths[0] + m_prescribedSubLengths[1] < 1.5*m_des_subunit_len)
        {
                // First we change the centre subunit
                --s_total_subunits;
                --m_num_subunits;
                if (severing)
                {
                    ++s_total_subunits;
                    // Need to reverse the above, if this is called during severing
                }


                for (unsigned int i = 0; i < m_subunit_unitVecs.size()-1; ++i)
                {
                        m_subunit_unitVecs[i] = m_subunit_unitVecs[i+1];
                }

                // Dealing with daughters on the old subunit, they need to move up to
                // subunit 1 before they ALL move back down again
                std::vector<int> daughterIDs_copy = m_daughterIDs[0];
                for (unsigned int i = 0; i < daughterIDs_copy.size(); ++i)
                {
                        int j = daughterIDs_copy[i];
                        actinVec[j].incrementBRsub();
                        purgeBranchFromDaughterVec(j);
                        addBranchtoParentVector(1,j);
                }

                if (steric)
                        stericGrid.resetCells(*this, 0); // delete from steric grid too

                for (unsigned int i = 0; i < m_actualSubLengths.size()-1; ++i)
                {
                        m_actualSubLengths[i] = m_actualSubLengths[i+1];
                        if (steric)
                                m_stericCells[i] = m_stericCells[i+1];

                        for (unsigned int j = 0; j < m_daughterIDs[i+1].size(); ++j)
                        {
                                int k = m_daughterIDs[i+1][j];
                                actinVec[k].decrementBRsub();
                        }

                        m_daughterIDs[i] = m_daughterIDs[i+1];
                }

                for (unsigned int i = 1; i < m_points.size()-1; ++i)
                {
                        // start at 1 to not overwrite point 0 here,
                        // same as first line after centre loop in rem Barb points func
                        m_points[i] = m_points[i+1];
                }

                // Move any tether points down but if any are on subunit 0 they need to
                // stay on sub 0
                int numTethers = m_tetherMono.size();
                for (int i = 0; i < numTethers; ++i)
                {
                    if (m_tetherSubunits[i] > 0)
                        m_tetherSubunits[i] -= 1;
                }

                // Move any crosslink points down too
                int numCLinks = getNumCLinks();
                for (int i = 0; i < numCLinks; ++i)
                {
                    int oldCLSub = m_cLinkActinAndSites[i][1];
                    if (m_cLinkActinAndSites[i][1] > 0)
                        setClSub(i, oldCLSub-1, actinVec);
                }

                // Then pop
                m_points.pop_back();
                m_actualSubLengths.pop_back();
                if (steric)
                {
                        m_stericCells.pop_back();
                        stericGrid.moveCellsDown(*this);
                }
                m_daughterIDs.pop_back();
                m_subunit_unitVecs.pop_back();
                // recalculate end angle

                double dy = m_points[2][1] - m_points[0][1];
                double dx = m_points[2][0] - m_points[0][0];
                double mag = sqrt(dx*dx + dy*dy);
                m_subunit_unitVecs[1][0] = dx/mag;
                m_subunit_unitVecs[1][1] = dy/mag;
                m_subunit_unitVecs[0][0] = dx/mag;
                m_subunit_unitVecs[0][1] = dy/mag;



                // Now we need to rejig the penultimate point (the new point 1)
                std::array<double,2> endpoint { m_points[0][0], m_points[0][1] };
                std::array<double,2> ppenultimate { m_points[2][0], m_points[2][1] };
                std::array<double,2> penultimate { m_points[1][0], m_points[1][1] };
                double lengthToEnd = sqrt(distanceBetPoints2DSQR(endpoint, penultimate)) + sqrt(distanceBetPoints2DSQR(ppenultimate, penultimate));

                double endsubL;

                lengthToEnd -= m_actualSubLengths[1]/2;
                lengthToEnd /= 2;
                endsubL = lengthToEnd;
                lengthToEnd += m_actualSubLengths[1]/2;

                // lengthToEnd is now the length of the penultimate subunit

                m_points[1][0] = m_points[2][0] - (lengthToEnd*m_subunit_unitVecs[1][0]); // x
                m_points[1][1] = m_points[2][1] - (lengthToEnd*m_subunit_unitVecs[1][1]); // y

                m_points[0][0] = m_points[1][0] - (endsubL*m_subunit_unitVecs[0][0]); // x
                m_points[0][1] = m_points[1][1] - (endsubL*m_subunit_unitVecs[0][1]); // y


                updateSubLengths();

                lengthToEnd = m_prescribedSubLengths[0] + m_prescribedSubLengths[1] + m_prescribedSubLengths[2];
                for (unsigned int i = 0; i < m_actualSubLengths.size(); ++i)
                {
                        m_prescribedSubLengths[i] = m_prescribedSubLengths[i+1];
                }
                m_prescribedSubLengths.pop_back();


                lengthToEnd -= m_prescribedSubLengths[1]/2;
                lengthToEnd /= 2;
                endsubL = lengthToEnd;
                lengthToEnd += m_prescribedSubLengths[1]/2;

                m_prescribedSubLengths[1] = lengthToEnd;
                m_prescribedSubLengths[0] = endsubL;

                moveBranchSub_Rem_point(actinVec);
                moveTetherPoint_Rem_point();
                moveClPoint_Rem_point(actinVec);
                updateTetherLoc();
                updateCrossLinkLoc(actinVec);
                updateUnitVecs();
                if (steric)
                        stericGrid.updateCellsRemPointPoint(*this);
        }

}

bool Actin::straightenFila(std::vector<Actin> &actinVec,
                           const std::vector<ExcZone> &excZones,
                           const std::vector<MembraneWall> &memWalls,
                           const std::vector<Membrane> &membranes,
                           Cortex &cortex,
                           const bool steric, StericGrid &stericGrid,
                           bool barbedEnd)
{
        /*
           Function that is called in the case of depolymerisation from a bent filament
           (L > 3*m_des_subunit_len) to a straight one
         */
        // Filament is now straight, the global angle needs to be recalculated
        // we recalculate the new global angle (the angle of the centre subunit)

        // should end with 3 subs!
        saveOldPos(actinVec); // saves any branches any crosslinks

        if (m_num_subunits != 4)
        {
            std::cout << "Straighten num subs: " << m_num_subunits << std::endl;
        }
        assert(m_num_subunits == 4);
        std::vector< std::array<double, 2> > oldSubunitVecs;
        std::vector<double> oldPresLens;
        std::vector<double> oldActLens;

        double dyG = m_points[m_num_subunits][1] - m_points[0][1];
        double dxG = m_points[m_num_subunits][0] - m_points[0][0];
        double mag = sqrt(dxG*dxG + dyG*dyG);

        for (int i = 0; i < m_num_subunits; ++i)
        {
                oldSubunitVecs.push_back(m_subunit_unitVecs[i]);
                m_subunit_unitVecs[i][0] = dxG/mag;
                m_subunit_unitVecs[i][1] = dyG/mag;

                oldPresLens.push_back(m_prescribedSubLengths[i]);
                oldActLens.push_back(m_actualSubLengths[i]);
        }



        if (barbedEnd)
        {
            m_prescribedSubLengths[0] = m_length/4;
            m_actualSubLengths[0] = m_length/4;

            m_prescribedSubLengths[1] = m_length/2;
            m_actualSubLengths[1] = m_length/2;

            m_prescribedSubLengths[2] = m_length/4;
            m_actualSubLengths[2] = m_length/4;

            for (int i = 3; i < m_num_subunits; ++i)
            {
                // Extra points are (temporarily) on top of each other
                m_prescribedSubLengths[i] = 0;
                m_actualSubLengths[i] = 0;
            }
        }
        else
        {
            m_prescribedSubLengths[m_num_subunits-1] = m_length/4;
            m_actualSubLengths[m_num_subunits-1] = m_length/4;

            m_prescribedSubLengths[m_num_subunits-2] = m_length/2;
            m_actualSubLengths[m_num_subunits-2] = m_length/2;

            m_prescribedSubLengths[m_num_subunits-3] = m_length/4;
            m_actualSubLengths[m_num_subunits-3] = m_length/4;

            for (int i = m_num_subunits-4; i >= 0; --i)
            {
                m_prescribedSubLengths[i] = 0;
                m_actualSubLengths[i] = 0;
            }
        }

        for (int i = 0; i < m_num_subunits; ++i)
        {
            checkAndMoveBranch(actinVec, i);
        }
        checkAndMoveCl(actinVec);

        updatePointsUp(0);
        m_flexible = false;
        rotbranches(actinVec, 0);
        rotRigidCLinked(actinVec, 0);

        updateTetherLoc();
        updateCrossLinkLocs(actinVec);

        if (steric)
        {
            stericGrid.resetAndUpdateAllCLinksAndDaughters(*this, actinVec);
        }

        if (check_bend(actinVec, excZones, memWalls,
                       membranes, cortex, steric, stericGrid))
        {
                // rebend filament
                for (int i = 0; i < m_num_subunits; ++i)
                {
                        m_subunit_unitVecs[i][0] = oldSubunitVecs[i][0];
                        m_subunit_unitVecs[i][1] = oldSubunitVecs[i][1];

                        m_prescribedSubLengths[i] = oldPresLens[i];
                        m_actualSubLengths[i] = oldActLens[i];
                }

                for (int i = 0; i < m_num_subunits; ++i)
                {
                    checkAndMoveBranch(actinVec, i);
                }
                checkAndMoveCl(actinVec);

                updatePointsUp(0);

                m_flexible = true;

                if (m_cLinkMaster)
                    m_cLinkMaster = false;

                resetToOldPos(actinVec);

                updateTetherLoc();

                updateCrossLinkLocs(actinVec);

                if (barbedEnd)
                {
                        cancelBarbDepol(actinVec, steric, stericGrid);
                }
                else
                {
                        cancelPointDepol(actinVec, steric, stericGrid);
                }

                // Have to do this again
                for (int i = 0; i < m_num_subunits; ++i)
                {
                    checkAndMoveBranch(actinVec, i);
                }
                checkAndMoveCl(actinVec);


                // Update back
                if (steric)
                {
                    stericGrid.resetAndUpdateAllCLinksAndDaughters(*this, actinVec);
                }



                return true;
        }
        else if (m_num_subunits > 3)
        {
            if (barbedEnd)
            {
                removePointsBarbed(actinVec, steric, stericGrid);
            }
            else
            {
                removePointsPointed(actinVec, steric, stericGrid);
            }
        }

        assert(m_num_subunits == 3);

        return false;
}

bool Actin::straightenFila()
{
        /*
           Reduced version
           Function that is called in the case of severing
           Subunit lengths already been changed in severing function
         */


        double dyG = m_points[m_num_subunits][1] - m_points[0][1];
        double dxG = m_points[m_num_subunits][0] - m_points[0][0];
        double mag = sqrt(dxG*dxG + dyG*dyG);
        for (int i = 0; i < m_num_subunits; ++i)
        {
                m_subunit_unitVecs[i][0] = dxG/mag;
                m_subunit_unitVecs[i][1] = dyG/mag;
        }

        updatePointsUp(0);
        m_flexible = false;

        return false; // uncomment this to turn off the rebending feature
}

void Actin::straightenEnd(bool barb)
{
        /*
           After severing need to straighten the end, so the last two subunits are the same
           direction
           if barb == true then its the barbed end
           else: pointed end
         */

        int end = 2;
        if (barb)
        {
                end = m_num_subunits;
        }

        double dyG = m_points[end][1] - m_points[end-2][1];
        double dxG = m_points[end][0] - m_points[end-2][0];
        double mag = sqrt(dxG*dxG + dyG*dyG);

        m_subunit_unitVecs[end-1][0] = dxG/mag;
        m_subunit_unitVecs[end-1][1] = dyG/mag;

        m_subunit_unitVecs[end-2][0] = dxG/mag;
        m_subunit_unitVecs[end-2][1] = dyG/mag;

        if (barb)
        {
                updatePointsUp(end-2);
        }
        else
        {
                updatePointsDown(end);
        }

}


void Actin::checkAndMoveBranchDown(std::vector<Actin> &actinVec, int i)
{
        /*
           General function that has input of i which is the subunit of the mother
           Checks all branches on the sub and determines if they need to move down
         */

        int vecsize = actinVec.size();
        int numDaughts = m_daughterIDs[i].size();
        std::vector<int> daughterIDs_copy = m_daughterIDs[i];
        for (int j = 0; j < numDaughts; ++j)
        {
                int k = daughterIDs_copy[j];
                if (k >= vecsize)
                {
                        // filament doesn't exist in the vector yet
                        continue;
                }

                // j is now the id of the branched filament
                // m_id is the id of the mother filament
                // i is the mother's subunit ID


                int oldbranchSub = actinVec[k].getBranchSubunit();
                assert(oldbranchSub == i);
                double junk;
                int newbranchSub = findSubunit(actinVec[k].getMotherMonoID(), junk);


                if (newbranchSub !=  oldbranchSub)
                {
                        // remove this daughter from this vector and add to
                        // the next sub down
                        assert(newbranchSub == i-1);
                        purgeBranchFromDaughterVec(k);
                        addBranchtoParentVector(i-1,k);
                        actinVec[k].setBranchSubunit(newbranchSub);

                }


        }
}

void Actin::checkAndMoveBranchUp(std::vector<Actin> &actinVec, int i)
{
        /*
           General function that has input of i which is the subunit of the mother
           Checks all branches on the sub and determines if they need to move up
         */

        int vecsize = actinVec.size();
        int numDaughts = m_daughterIDs[i].size();
        std::vector<int> daughterIDs_copy = m_daughterIDs[i];
        for (int j = 0; j < numDaughts; ++j)
        {
                int k = daughterIDs_copy[j];
                if (k >= vecsize)
                {
                        // filament doesn't exist in the vector yet
                        continue;
                }

                // j is now the id of the branched filament
                // m_id is the id of the mother filament
                // i is the mother's subunit ID


                int oldbranchSub = actinVec[k].getBranchSubunit();
                assert(oldbranchSub == i);
                double junk;
                int newbranchSub = findSubunit(actinVec[k].getMotherMonoID(), junk);

                if (newbranchSub !=  oldbranchSub)
                {
                        // remove this daughter from this vector and add to
                        // the next sub down
                        assert(newbranchSub == i+1);
                        purgeBranchFromDaughterVec(k);
                        addBranchtoParentVector(i+1,k);
                        actinVec[k].setBranchSubunit(newbranchSub);

                }

        }
}

void Actin::checkAndMoveBranch(std::vector<Actin> &actinVec, int i)
{
        /*
           General function that has input of i which is the subunit of the mother
           Checks all branches on the sub and determines if they need to move
           Can move up or down, we dont know!
         */

        int vecsize = actinVec.size();
        int numDaughts = m_daughterIDs[i].size();
        std::vector<int> daughterIDs_copy = m_daughterIDs[i];
        for (int j = 0; j < numDaughts; ++j)
        {
                int k = daughterIDs_copy[j];
                if (k >= vecsize)
                {
                        // filament doesn't exist in the vector yet
                        continue;
                }

                // j is now the id of the branched filament
                // m_id is the id of the mother filament
                // i is the mother's subunit ID


                int oldbranchSub = actinVec[k].getBranchSubunit();
                assert(oldbranchSub == i);
                double junk;
                int newbranchSub = findSubunit(actinVec[k].getMotherMonoID(), junk);

                if (newbranchSub !=  oldbranchSub)
                {
                        // remove this daughter from this vector and add to
                        // the next sub down
                        purgeBranchFromDaughterVec(k);
                        addBranchtoParentVector(newbranchSub, k);
                        actinVec[k].setBranchSubunit(newbranchSub);

                }

        }
}

void Actin::moveBranchSub_down(std::vector<Actin> &actinVec)
{
        // Check for moving branches
        /*
           This is for either barbed end poly or pointed end depoly, when the branch
           can move down a subunit.
           Called for when we are in initial straight regime
         */

        for (int i = 0; i < m_num_subunits; ++i)
        {
                checkAndMoveBranchDown(actinVec, i);
        }

}

void Actin::moveBranchSub_up(std::vector<Actin> &actinVec)
{
        // Check for moving branches
        /*
           This is for either barbed end depoly or pointed end poly, when the branch
           can move up a subunit.
           Called for when we are in initial straight regime
         */

        for (int i = 0; i < m_num_subunits; ++i)
        {
                checkAndMoveBranchUp(actinVec, i);
        }

}

void Actin::moveBranchSub_poly_barb(std::vector<Actin> &actinVec)
{
        // Check for moving branches
        // After a polymerisation event, check all branches on the final subunit
        // if any have "moved" we know to move the branches down to the penultimate
        // subunit

        // check daughter vecs
        // For each subunit we need to find all the branches
        int i = m_num_subunits - 1;
        checkAndMoveBranchDown(actinVec, i);
        i = m_num_subunits - 2;
        checkAndMoveBranchDown(actinVec, i);
        i = m_num_subunits - 3;
        checkAndMoveBranchDown(actinVec, i);

}

void Actin::moveBranchSub_poly_point(std::vector<Actin> &actinVec)
{
        // Check for moving branches
        // After a polymerisation event, check all branches on the final subunit
        // if any have "moved" we know to move the branches down to the penultimate
        // subunit

        // check daughter vecs
        // For each subunit we need to find all the branches
        int i = 0;
        checkAndMoveBranchUp(actinVec, i);
        i = 1;
        checkAndMoveBranchUp(actinVec, i);
        i = 2;
        checkAndMoveBranchUp(actinVec, i);

}

void Actin::moveBranchSub_depoly_barb(std::vector<Actin> &actinVec)
{
        // Check for moving branches
        // After a polymerisation event, check all branches on the final subunit
        // if any have "moved" we know to move the branches down to the penultimate
        // subunit
        // check daughter vecs
        // For each subunit we need to find all the branches
        int i = m_num_subunits - 2;
        checkAndMoveBranchUp(actinVec, i);
        i = m_num_subunits - 3;
        checkAndMoveBranchUp(actinVec, i);
        i = m_num_subunits - 4;
        checkAndMoveBranchUp(actinVec, i);

}

void Actin::moveBranchSub_depoly_point(std::vector<Actin> &actinVec)
{
        // Check for moving branches
        // After a polymerisation event, check all branches on the final subunit
        // if any have "moved" we know to move the branches down to the penultimate
        // subunit
        // check daughter vecs
        // For each subunit we need to find all the branches
        int i = 1;
        checkAndMoveBranchDown(actinVec, i);
        i = 2;
        checkAndMoveBranchDown(actinVec, i);
        i = 3;
        checkAndMoveBranchDown(actinVec, i);

}

void Actin::moveBranchSub_Add_barb(std::vector<Actin> &actinVec)
{
        // Check for moving branches
        // After a polymerisation event, check all branches on the final subunit
        // if any have "moved" we know to move the branches down to the penultimate
        // subunit

        // check daughter vecs
        // For each subunit we need to find all the branches
        for (int i = m_num_subunits-3; i < m_num_subunits-1; ++i)
        {
                checkAndMoveBranchUp(actinVec, i);
        }
}

void Actin::moveBranchSub_Add_point(std::vector<Actin> &actinVec)
{
        // Check for moving branches
        // After a polymerisation event, check all branches on the final subunit
        // if any have "moved" we know to move the branches down to the penultimate
        // subunit

        // check daughter vecs
        // For each subunit we need to find all the branches
        for (int i = 1; i < 3; ++i)
        {
                checkAndMoveBranchDown(actinVec, i);
        }
}

void Actin::moveBranchSub_Rem_barb(std::vector<Actin> &actinVec)
{
        // Check for moving branches
        // After a polymerisation event, check all branches on the final subunit
        // if any have "moved" we know to move the branches down to the penultimate
        // subunit

        // check daughter vecs
        // For each subunit we need to find all the branches
        // do we need to check both m_num_subunits-2 and m_num_subunits-1?

        for (int i = m_num_subunits-2; i < m_num_subunits; ++i)
        {
                checkAndMoveBranchDown(actinVec, i);
        }
}

void Actin::moveBranchSub_Rem_point(std::vector<Actin> &actinVec)
{
        // Check for moving branches
        // After a polymerisation event, check all branches on the final subunit
        // if any have "moved" we know to move the branches down to the penultimate
        // subunit

        // check daughter vecs
        // For each subunit we need to find all the branches
        for(int i = 0; i < 2; ++i)
        {
                checkAndMoveBranchUp(actinVec, i);
        }
}

void Actin::updateSubLengths()
{
        for (unsigned int i = 0; i < m_actualSubLengths.size(); ++i)
        {
                std::array<double,2> point1 { m_points[i][0], m_points[i][1] };
                std::array<double,2> point2 { m_points[i+1][0], m_points[i+1][1] };
                m_actualSubLengths[i] = sqrt(distanceBetPoints2DSQR(point1, point2));
        }

}

void Actin::updateSubLength(int p)
{
    // Update the sublength of sub p

    std::array<double,2> point1 { m_points[p][0], m_points[p][1] };
    std::array<double,2> point2 { m_points[p+1][0], m_points[p+1][1] };
    m_actualSubLengths[p] = sqrt(distanceBetPoints2DSQR(point1, point2));
}


void Actin::updatePointsUp(int startpoint)
{
        // Called after bending

        for (int i = startpoint+1; i <= m_num_subunits; ++i)
        {
                m_points[i][0] = m_points[i-1][0] + (m_actualSubLengths[i-1]*m_subunit_unitVecs[i-1][0]); // x dimension
                m_points[i][1] = m_points[i-1][1] + (m_actualSubLengths[i-1]*m_subunit_unitVecs[i-1][1]); // y dimension
                m_points[i][2] = 0.0; // z dimension
        }
}

void Actin::updatePointsDown(int startpoint)
{
        // Called after bending

        for (int i = startpoint-1; i >= 0; --i)
        {
                m_points[i][0] = m_points[i+1][0] - (m_actualSubLengths[i]*m_subunit_unitVecs[i][0]); // x dimension
                m_points[i][1] = m_points[i+1][1] - (m_actualSubLengths[i]*m_subunit_unitVecs[i][1]); // y dimension
                m_points[i][2] = 0.0; // z dimension

        }
}

void Actin::regenerate_pos(const ProteinRegion &region)
{
        // Function that regenerates trimers and is called if the check initial
        // positions function returns true


        // define points local to axis aligned, origined rectangle
        if (region.getCircBool() || region.getRingBool())
        {
                double r = region.getRDist()(rng.m_mersenne);
                double theta = region.getThetaDist()(rng.m_mersenne);

                m_points[1][0] = region.getCentre()[0] + r*cos(theta);
                m_points[1][1] = region.getCentre()[1] + r*sin(theta);
        }
        else
        {
                double x_loc = region.getWidthDist()(rng.m_mersenne);
                double y_loc = region.getHeightDist()(rng.m_mersenne);

                m_points[1][0] = region.getBLPoint()[0] + ((x_loc*region.getUVec()[0]) - (y_loc*region.getUVec()[1]));
                m_points[1][1] = region.getBLPoint()[1] + ((x_loc*region.getUVec()[1]) + (y_loc*region.getUVec()[0]));
        }

        m_points[0] = { m_points[1][0] - ((m_length/4) * m_subunit_unitVecs[0][0]), m_points[1][1] - ((m_length/4) * m_subunit_unitVecs[0][1]), 0 };

        m_points[2] = { m_points[1][0] + ((m_length/2) * m_subunit_unitVecs[1][0]), m_points[2][1] + ((m_length/2) * m_subunit_unitVecs[1][1]), 0 };
        m_points[3] = { m_points[2][0] + ((m_length/4) * m_subunit_unitVecs[2][0]), m_points[3][1] + ((m_length/4) * m_subunit_unitVecs[2][1]), 0 };
}

bool Actin::check_barbed_polymerisation_Grid(const std::vector<Actin> &actinVec,
                                             const bool steric,
                                             StericGrid &stericGrid)
{
        /* Function that does a distance check between one filament that is
           polymerised, and every other filament in the simulation using GeoTools
           Parallised using OpenMP
         */

        if (!steric || !m_steric)
                return false;

        int startSub = 0;
        int endSub = m_num_subunits;
        if (getNumCLinks() == 0)
        {
            // If no crosslinks then only need to check the very end sub
            startSub = m_num_subunits-1;
        }
        else if (m_flexible)
        {
            // if flexible: Only need to check the last two subs because they have changed
            startSub = m_num_subunits-2;
        }
        // If none of the two above are true then have to check everything!

        for (int h = startSub; h < endSub; ++h)
        {
            gte::Segment<2,double> filament_1;
            filament_1.p[1][0] = m_points[h][0];
            filament_1.p[0][0] = m_points[h+1][0]; // check dist of end sub
            filament_1.p[1][1] = m_points[h][1];
            filament_1.p[0][1] = m_points[h+1][1];


            std::vector<int> barbCellIDs = m_stericCells[h];
            std::vector<int> cellsToCheck = stericGrid.getCellsToCheck(barbCellIDs);
            std::vector<std::array<int,2> > checkedSubs;

            for (unsigned int i = 0; i < cellsToCheck.size(); ++i)
            {
                    std::vector<std::array<int,2> > cellContents = stericGrid.getCellContents(cellsToCheck[i]);

                    for (unsigned int j = 0; j < cellContents.size(); ++j)
                    {

                            int actinID = cellContents[j][0];
                            int actinSubID = cellContents[j][1];

                            if (std::find(checkedSubs.begin(), checkedSubs.end(), cellContents[j]) != checkedSubs.end())
                            {
                                    // We have checked this filament before!
                                    continue;
                            }
                            checkedSubs.push_back(cellContents[j]);

                            if (actinID == m_id)
                            {
                              // Same filament, ignore same sub, and adjacent subs
                                if (!m_flexible)
                                {
                                    // Rigid filament, cant hinder self
                                    continue;
                                }

                                if (actinSubID == m_num_subunits-1 || actinSubID == m_num_subunits-2 || actinSubID == m_num_subunits-3)
                                {
                                    continue;
                                }

                            }

                            if (actinID == m_parent_id)
                            {
                                // other filament is parent, ignore sub that the branch sits on AND ADJACENT SUBS
                                if (m_length < 2*s_segmentationLength && actinVec[actinID].getLength() < 2*s_segmentationLength)
                                {
                                    // both mother and daughter rigid, cant hinder
                                    continue;
                                }

                                if (actinSubID == m_branchSubUnit || actinSubID == m_branchSubUnit-1 || actinSubID == m_branchSubUnit+1)
                                {
                                    continue;
                                }

                            }

                            if (actinVec[actinID].getParentID() == m_id)
                            {
                                if (m_length < 2*s_segmentationLength && actinVec[actinID].getLength() < 2*s_segmentationLength)
                                {
                                    // both mother and daughter rigid, cant hinder
                                    continue;
                                }

                                if (m_num_subunits-1 == actinVec[actinID].getBranchSubunit() || m_num_subunits-2 == actinVec[actinID].getBranchSubunit() || m_num_subunits-3 == actinVec[actinID].getBranchSubunit())
                                {
                                    continue;
                                }

                            }

                            bool cont = false; // continue
                            for (unsigned int k = 0; k < m_cLinkActinAndSites.size(); ++k)
                            {
                                // if other filament is linked
                                if (m_cLinkActinAndSites[k][2] == actinID)
                                {
                                    // ignore if link is on last two subs
                                    if (h-1 == m_cLinkActinAndSites[k][1] || h == m_cLinkActinAndSites[k][1] || h+1 == m_cLinkActinAndSites[k][1])
                                    {
                                          // if the other sub has the link on too!
                                          if (m_cLinkActinAndSites[k][5] == actinSubID || m_cLinkActinAndSites[k][5] == actinSubID-1 || m_cLinkActinAndSites[k][5] == actinSubID+1)
                                          {
                                              cont = true;
                                              break;
                                          }
                                    }
                                }
                            }
                            if (cont)
                            {
                                continue;
                            }

                            gte::Segment<2,double> filament_2;

                            filament_2.p[0][0] = actinVec[actinID].getPoint(actinSubID)[0];
                            filament_2.p[1][0] = actinVec[actinID].getPoint(actinSubID+1)[0];
                            filament_2.p[0][1] = actinVec[actinID].getPoint(actinSubID)[1];
                            filament_2.p[1][1] = actinVec[actinID].getPoint(actinSubID+1)[1];

                            RobustQuery min_d_GTE;
                            auto result = min_d_GTE(filament_1, filament_2);
                            double distance = result.distance;

                            if (distance < (m_stericRadius + actinVec[actinID].getStericRadius()))
                            {
                                    return true;
                            }
                    }
            }
        }
        return false;
}

bool Actin::check_pointed_polymerisation_excVol_Grid(const std::vector<Actin> &actinVec,
                                                     const bool steric,
                                                     StericGrid &stericGrid)
{

        if (!steric)
                return false;


        int startSub = 0;
        int endSub = m_num_subunits;
        if (getNumCLinks() == 0)
        {
            // If no crosslinks then only need to check the very end sub
            endSub = 1;
        }
        else if (m_flexible)
        {
            // if flexible: Only need to check the last two subs because they have changed
            startSub = 2;
        }
        // If none of the two above are true then have to check everything!

        for (int h = startSub; h < endSub; ++h)
        {

            gte::Segment<2,double> filament_1; // new monomer
            filament_1.p[0][0] = m_points[h+1][0];
            filament_1.p[1][0] = m_points[h][0];
            filament_1.p[0][1] = m_points[h+1][1];
            filament_1.p[1][1] = m_points[h][1];



            std::vector<int> pointCellIDs = m_stericCells[h];
            std::vector<int> cellsToCheck = stericGrid.getCellsToCheck(pointCellIDs);

            std::vector<std::array<int,2> > checkedSubs;

            for (unsigned int i = 0; i < cellsToCheck.size(); ++i)
            {

                    std::vector<std::array<int,2> > cellContents = stericGrid.getCellContents(cellsToCheck[i]);

                    for (unsigned int j = 0; j < cellContents.size(); ++j)
                    {
                            int actinID = cellContents[j][0];
                            int actinSubID = cellContents[j][1];

                            if (std::find(checkedSubs.begin(), checkedSubs.end(), cellContents[j]) != checkedSubs.end())
                            {
                                    // We have checked this filament before!
                                    continue;
                            }
                            checkedSubs.push_back(cellContents[j]);


                            if (actinID == m_id)
                            {
                              // Same filament, ignore same sub, and adjacent subs
                                if (!m_flexible)
                                {
                                    // Rigid filament, cant hinder self
                                    continue;
                                }

                                if (actinSubID == 0 || actinSubID == 1 || actinSubID == 2)
                                {
                                    continue;
                                }

                            }

                            if (actinID == m_parent_id)
                            {
                                  //Other filament is parent,
                                  // branches cannot poitned end poly anyway!
                                  continue;
                            }

                            if (actinVec[actinID].getParentID() == m_id)
                            {
                                if (m_length < 2*s_segmentationLength && actinVec[actinID].getLength() < 2*s_segmentationLength)
                                {
                                    // both mother and daughter rigid, cant hinder
                                    continue;
                                }

                                if (0 == actinVec[actinID].getBranchSubunit() || 1 == actinVec[actinID].getBranchSubunit() || 2 == actinVec[actinID].getBranchSubunit())
                                {
                                    continue;
                                }

                            }

                            bool cont = false;
                            for (unsigned int k = 0; k < m_cLinkActinAndSites.size(); ++k)
                            {
                                // if the other filament is linked to this one
                                if (m_cLinkActinAndSites[k][2] == actinID)
                                {   // if this link is on the first two subs
                                    if (h-1 == m_cLinkActinAndSites[k][1] || h == m_cLinkActinAndSites[k][1] || h+1 == m_cLinkActinAndSites[k][1])
                                    {
                                        // if the other sub has the link on too!
                                        if (m_cLinkActinAndSites[k][5] == actinSubID || m_cLinkActinAndSites[k][5] == actinSubID-1 || m_cLinkActinAndSites[k][5] == actinSubID+1)
                                        {
                                            cont = true;
                                            break;
                                        }
                                    }
                                }
                            }
                            if (cont)
                                continue;


                            gte::Segment<2,double> filament_2;

                            filament_2.p[0][0] = actinVec[actinID].getPoint(actinSubID)[0];
                            filament_2.p[1][0] = actinVec[actinID].getPoint(actinSubID+1)[0];
                            filament_2.p[0][1] = actinVec[actinID].getPoint(actinSubID)[1];
                            filament_2.p[1][1] = actinVec[actinID].getPoint(actinSubID+1)[1];


                            RobustQuery min_d_GTE;
                            auto result = min_d_GTE(filament_1, filament_2);
                            double distance = result.distance;

                            if (distance < (m_stericRadius + actinVec[actinID].getStericRadius()))
                            {
                                    return true;
                            }

                    }

            }
        }
        return false;
}

bool Actin::check_nucleation_excVol_Grid(const std::vector<Actin> &actinVec,
                                         const std::vector<ExcZone> &excZones,
                                         const std::vector<Membrane> &membranes,
                                         const std::vector<MembraneWall> &memWalls,
                                         Cortex &cortex,
                                         const bool steric, StericGrid &stericGrid)
{

        for (unsigned int i = 0; i < excZones.size(); ++i)
        {
                if (excZones[i].check_ex_vol(*this))
                {
                        return true;
                }
        }

        for (unsigned int j = 0; j < memWalls.size(); ++j)
        {
                if (memWalls[j].checkExVolBM(*this))
                        return true;
        }

        if (!steric)
                return false;


        std::vector<int> occCells; // cells that are occupied

        for (int i = 0; i < m_num_subunits; ++i)
        {
                for (unsigned int j = 0; j < m_stericCells[i].size(); ++j)
                {
                        if (std::find(occCells.begin(), occCells.end(),m_stericCells[i][j])==occCells.end()) // if it is not in occCells already
                        {
                                occCells.push_back(m_stericCells[i][j]);
                        }
                }
        }


        std::vector<int> cellsToCheck = stericGrid.getCellsToCheck(occCells);

        // Need to avoid having repeats in cellsToCheck
        if (check_min_dist_Grid(actinVec, cellsToCheck, stericGrid)
            || Membrane::s_checkExVolNuc(*this, membranes, cellsToCheck, stericGrid)
            || Cortex::s_checkExVolNuc(*this, cortex, cellsToCheck, stericGrid))
        {
                return true;
        }
        else
        {
                return false;
        }
}

bool Actin::checkDiffusionStructure(std::vector<Actin> &actinVec,
                                       const std::vector<ExcZone> &excZones,
                                       const std::vector<MembraneWall> &memWalls,
                                       const std::vector<Membrane> &membranes,
                                       Cortex &cortex, const bool steric,
                                       StericGrid &stericGrid,
                                       const bool diffBool)
{
    // called for checking all structure, avoids double counting
    std::vector<int> checkedFilas;
    for (unsigned int i = 0; i < actinVec.size(); ++i)
    {
        if (actinVec[i].getStructureID() == m_structure_id)
        {
            if (actinVec[i].check_diffusion_excVolGrid2(actinVec,
                                                       excZones,
                                                       memWalls,
                                                       membranes,
                                                       cortex,
                                                       steric,
                                                       stericGrid,
                                                       diffBool, checkedFilas))
            {
                    // If we have violated exc area we move back
                    return false;
            }
        }
    }

    return true;

}

bool Actin::checkDiffusionClStructure(std::vector<Actin> &actinVec,
                                       const std::vector<ExcZone> &excZones,
                                       const std::vector<MembraneWall> &memWalls,
                                       const std::vector<Membrane> &membranes,
                                       Cortex &cortex, const bool steric,
                                       StericGrid &stericGrid,
                                       const bool diffBool)
{
    // called for checking just the crosslinked structure, avoids double counting
    std::vector<int> checkedFilas;
    for (unsigned int i = 0; i < actinVec.size(); ++i)
    {
        if (actinVec[i].getCLStructureID() == m_crossStructure_id)
        {
            if (actinVec[i].check_diffusion_excVolGrid2(actinVec,
                                                       excZones,
                                                       memWalls,
                                                       membranes,
                                                       cortex,
                                                       steric,
                                                       stericGrid,
                                                       diffBool, checkedFilas))
            {
                    // If we have violated exc area we move back
                    return false;
            }
        }
    }

    return true;

}

bool Actin::check_diffusion_excVolGrid(const std::vector<Actin> &actinVec,
                                       const std::vector<ExcZone> &excZones,
                                       const std::vector<MembraneWall> &memWalls,
                                       const std::vector<Membrane> &membranes,
                                       Cortex &cortex, const bool steric,
                                       StericGrid &stericGrid,
                                       const bool diffBool)
{

        for (unsigned int i = 0; i < memWalls.size(); ++i)
        {
                if (memWalls[i].checkExVolBM(*this))
                        return true;
        }

        if (!steric)
                return false;

        for (int i = 0; i < m_num_subunits; ++i)
        {
                std::vector<int> occCells; // cells that are occupied
                for (unsigned int j = 0; j < m_stericCells[i].size(); ++j)
                {

                        occCells.push_back(m_stericCells[i][j]);
                }

                std::vector<int> cellsToCheck = stericGrid.getCellsToCheck(occCells);

                if (check_min_dist_Diff_Grid(i, actinVec, cellsToCheck, stericGrid, diffBool)
                    || Membrane::s_checkExVolBM(i, *this, membranes, cellsToCheck, stericGrid)
                    || ExcZone::s_checkExVol(i, *this, excZones)
                    || Cortex::s_checkExVolBM(i, *this, cortex, cellsToCheck, stericGrid))
                {
                        return true;
                }

        }
        return false;

}

bool Actin::check_diffusion_excVolGrid2(const std::vector<Actin> &actinVec,
                                       const std::vector<ExcZone> &excZones,
                                       const std::vector<MembraneWall> &memWalls,
                                       const std::vector<Membrane> &membranes,
                                       Cortex &cortex, const bool steric,
                                       StericGrid &stericGrid,
                                       const bool diffBool,
                                       std::vector<int> &checkedFilas)
{

        for (unsigned int i = 0; i < memWalls.size(); ++i)
        {
                if (memWalls[i].checkExVolBM(*this))
                        return true;
        }

        if (!steric)
                return false;

        for (int i = 0; i < m_num_subunits; ++i)
        {
                //std::cout << m_id << " ; " << i << std::endl;
                std::vector<int> occCells; // cells that are occupied
                for (unsigned int j = 0; j < m_stericCells[i].size(); ++j)
                {
                        occCells.push_back(m_stericCells[i][j]);
                }

                std::vector<int> cellsToCheck = stericGrid.getCellsToCheck(occCells);

                if (check_min_dist_Diff_Grid2(i, actinVec, cellsToCheck, stericGrid, diffBool, checkedFilas)
                    || Membrane::s_checkExVolBM(i, *this, membranes, cellsToCheck, stericGrid)
                    || ExcZone::s_checkExVol(i, *this, excZones)
                    || Cortex::s_checkExVolBM(i, *this, cortex, cellsToCheck, stericGrid))
                {
                        return true;
                }

        }
        checkedFilas.push_back(m_id);
        return false;

}

bool Actin::check_bend(std::vector<Actin> &actinVec,
                       const std::vector<ExcZone> &excZones,
                       const std::vector<MembraneWall> &memWalls,
                       const std::vector<Membrane> &membranes,
                       Cortex &cortex,
                       const bool steric, StericGrid &stericGrid)
{
        std::vector<int> checkedFilas;
        checkedFilas.push_back(m_id);

        std::vector<std::array<int,2> > checkedFilaSubs;
        // Function that calls check bending excVol
        // recursive so checks all the daughters, grandaughters etc

        if (check_bending_excVolGrid(actinVec, excZones, memWalls, membranes, cortex, steric, stericGrid, checkedFilaSubs))
        {

                // Checks the filament against all others
                return true;

        }
        else
        {
                // Must also consider all branches etc
                for (int i = 0; i < m_num_subunits; ++i)
                {
                        for (int j : m_daughterIDs[i])
                        {
                                // only consider short branches as longer branches will have not moved
                                if (std::find(checkedFilas.begin(), checkedFilas.end(), j) != checkedFilas.end())
                                {
                                        // We have checked this filament before!
                                        continue;
                                }
                                checkedFilas.push_back(j);

                                if ((actinVec[j].getLength() < 2*s_segmentationLength) && (actinVec[j].check_bend(actinVec,excZones, memWalls, membranes, cortex, steric, stericGrid, checkedFilas, checkedFilaSubs)))
                                {
                                        return true;
                                }
                        }
                }

                // Must also consider all short crosslinks too
                for (int i = 0; i < getNumCLinks(); ++i)
                {
                    int otherID = m_cLinkActinAndSites[i][2];

                    if (std::find(checkedFilas.begin(), checkedFilas.end(), otherID) != checkedFilas.end())
                    {
                            // We have checked this filament before!
                            continue;
                    }
                    checkedFilas.push_back(otherID);

                    if (actinVec[otherID].getCLStructureID() == m_crossStructure_id)
                    {
                        continue; // already accounted for with branches above
                    }

                    if (actinVec[otherID].getLength() < 2*s_segmentationLength && actinVec[otherID].check_bend(actinVec,excZones, memWalls, membranes, cortex, steric, stericGrid, checkedFilas, checkedFilaSubs))
                    {
                          return true;
                    }

                }

                return false;
        }
}

bool Actin::check_bend(std::vector<Actin> &actinVec,
                       const std::vector<ExcZone> &excZones,
                       const std::vector<MembraneWall> &memWalls,
                       const std::vector<Membrane> &membranes,
                       Cortex &cortex,
                       const bool steric, StericGrid &stericGrid,
                       std::vector<int> &checkedFilas,
                       std::vector<std::array<int,2> > &checkedFilaSubs)
{
        // Downstream of original
        // recursive so checks all the daughters, grandaughters etc
        if (check_bending_excVolGrid(actinVec, excZones, memWalls, membranes, cortex, steric, stericGrid, checkedFilaSubs))
        {

                // Checks the filament against all others
                return true;

        }
        else
        {
                // Must also consider all branches etc
                for (int i = 0; i < m_num_subunits; ++i)
                {
                        for (int j : m_daughterIDs[i])
                        {
                                // only consider short branches as longer branches will have not moved

                                if (std::find(checkedFilas.begin(), checkedFilas.end(), j) != checkedFilas.end())
                                {
                                        // We have checked this filament before!
                                        continue;
                                }
                                checkedFilas.push_back(j);

                                if ((actinVec[j].getLength() < 2*s_segmentationLength) && (actinVec[j].check_bend(actinVec,excZones, memWalls, membranes, cortex, steric, stericGrid, checkedFilas, checkedFilaSubs)))
                                {
                                        return true;
                                }
                        }
                }

                // Must also consider all short crosslinks too
                for (int i = 0; i < getNumCLinks(); ++i)
                {
                    int otherID = m_cLinkActinAndSites[i][2];

                    if (actinVec[otherID].getCLStructureID() == m_crossStructure_id)
                    {
                        continue; // already accounted for with branches above
                        // but hasnt been checked yet
                    }

                    if (std::find(checkedFilas.begin(), checkedFilas.end(), otherID) != checkedFilas.end())
                    {
                            // We have checked this filament before!
                            continue;
                    }
                    checkedFilas.push_back(otherID);


                    if (actinVec[otherID].getLength() < 2*s_segmentationLength && actinVec[otherID].check_bend(actinVec,excZones, memWalls, membranes, cortex, steric, stericGrid, checkedFilas, checkedFilaSubs))
                    {
                          return true;
                    }

                }

                return false;
        }
}

bool Actin::check_bending_excVolGrid(const std::vector<Actin> &actinVec,
                                     const std::vector<ExcZone> &excZones,
                                     const std::vector<MembraneWall> &memWalls,
                                     const std::vector<Membrane> &membranes,
                                     Cortex &cortex,
                                     const bool steric, StericGrid &stericGrid,
                                     std::vector<std::array<int,2> > &checkedFilaSubs)
{

        for (unsigned int i = 0; i < memWalls.size(); ++i)
        {
                if (memWalls[i].checkExVolBM(*this))
                        return true;
        }

        if (!steric)
                return false;

        for (int i = 0; i < m_num_subunits; ++i)
        {
                std::vector<int> occCells; // cells that are occupied
                for (unsigned int j = 0; j < m_stericCells[i].size(); ++j)
                {
                        occCells.push_back(m_stericCells[i][j]);
                }

                std::vector<int> cellsToCheck = stericGrid.getCellsToCheck(occCells);

                if (check_min_dist_BD_Grid(i, actinVec, cellsToCheck, stericGrid, checkedFilaSubs)
                    || Membrane::s_checkExVolBM(i, *this, membranes, cellsToCheck, stericGrid)
                    || ExcZone::s_checkExVol(i, *this, excZones)
                    || Cortex::s_checkExVolBM(i, *this, cortex, cellsToCheck, stericGrid))
                {
                        return true;
                }

        }
        return false;

}

bool Actin::checkBranchExcVolGrid(std::vector<Actin> &actinVec,
                                  const std::vector<ExcZone> &excZones,
                                  const std::vector<Membrane> &membranes,
                                  const std::vector<MembraneWall> &memWalls,
                                  Cortex &cortex,
                                  const bool steric, StericGrid &stericGrid)
{
        for (unsigned int i = 0; i < excZones.size(); ++i)
        {
                if (excZones[i].check_ex_vol(*this))
                        return true;
        }

        for (unsigned int j = 0; j < memWalls.size(); ++j)
        {
                if (memWalls[j].checkExVolBM(*this))
                        return true;
        }

        if (!steric)
                return false;

        std::vector<int> occCells; // cells that are occupied

        for (int i = 0; i < m_num_subunits; ++i)
        {
                for (unsigned int j = 0; j < m_stericCells[i].size(); ++j)
                {
                        if (std::find(occCells.begin(), occCells.end(),m_stericCells[i][j])==occCells.end()) // if it is not in occCells already
                        {
                                occCells.push_back(m_stericCells[i][j]);
                        }
                }
        }


        std::vector<int> cellsToCheck = stericGrid.getCellsToCheck(occCells);

        // Need to avoid having repeats in cellsToCheck
        if (check_min_dist_branch_Grid(actinVec, cellsToCheck, stericGrid)
            || Membrane::s_checkExVolNuc(*this, membranes, cellsToCheck, stericGrid)
            || Cortex::s_checkExVolNuc(*this, cortex, cellsToCheck, stericGrid))
        {
                return true;
        }
        else
        {
                return false;
        }
}

bool Actin::check_min_dist_Grid(const std::vector<Actin> &actinVec,
                                std::vector<int> &cellsToCheck,
                                StericGrid &stericGrid)
{
        // new trimer vs everything else?
        gte::Segment<2,double> filament_1;
        filament_1.p[0][0] = m_points[0][0];
        filament_1.p[1][0] = m_points[4][0];
        filament_1.p[0][1] = m_points[0][1];
        filament_1.p[1][1] = m_points[4][1];


        // To avoid repeated checks against the same subunits that may appear in
        // more than one cell, have a recorded of which ones we have checked
        std::vector<std::array<int,2> > checkedSubs;

        for (unsigned int i = 0; i < cellsToCheck.size(); ++i)
        {

                std::vector<std::array<int,2> > cellContents = stericGrid.getCellContents(cellsToCheck[i]);

                for (unsigned int j = 0; j < cellContents.size(); ++j)
                {
                        int actinID = cellContents[j][0];

                        if (std::find(checkedSubs.begin(), checkedSubs.end(), cellContents[j]) != checkedSubs.end())
                        {
                                // We have checked this filament before!
                                continue;
                        }
                        checkedSubs.push_back(cellContents[j]);

                        if (actinID == m_id)
                        {
                                // Ignore the same filament
                                continue;
                        }
                        int actinSubID = cellContents[j][1];

                        gte::Segment<2,double> filament_2;

                        filament_2.p[0][0] = actinVec[actinID].getPoint(actinSubID)[0];
                        filament_2.p[1][0] = actinVec[actinID].getPoint(actinSubID+1)[0];
                        filament_2.p[0][1] = actinVec[actinID].getPoint(actinSubID)[1];
                        filament_2.p[1][1] = actinVec[actinID].getPoint(actinSubID+1)[1];


                        RobustQuery min_d_GTE;
                        auto result = min_d_GTE(filament_1, filament_2);
                        double distance = result.distance;

                        if (distance < (m_stericRadius + actinVec[actinID].getStericRadius()))
                        {
                                return true;
                        }
                }
        }
        return false;
}

bool Actin::check_min_dist_branch_Grid(const std::vector<Actin> &actinVec,
                                       std::vector<int> &cellsToCheck,
                                       StericGrid &stericGrid)
{

        gte::Segment<2,double> filament_1;
        filament_1.p[0][0] = m_points[0][0];
        filament_1.p[1][0] = m_points[4][0];
        filament_1.p[0][1] = m_points[0][1];
        filament_1.p[1][1] = m_points[4][1];


        // To avoid repeated checks against the same subunits that may appear in
        // more than one cell, have a recorded of which ones we have checked
        std::vector<std::array<int,2> > checkedSubs;

        for (unsigned int i = 0; i < cellsToCheck.size(); ++i)
        {

                std::vector<std::array<int,2> > cellContents = stericGrid.getCellContents(cellsToCheck[i]);

                for (unsigned int j = 0; j < cellContents.size(); ++j)
                {
                        int actinID = cellContents[j][0];

                        if (std::find(checkedSubs.begin(), checkedSubs.end(), cellContents[j]) != checkedSubs.end())
                        {
                                // We have checked this filament before!
                                continue;
                        }
                        checkedSubs.push_back(cellContents[j]);

                        if (actinID == m_id || m_parent_id == actinID)
                        {
                                // Ignore the same filament
                                continue;
                        }
                        int actinSubID = cellContents[j][1];


                        gte::Segment<2,double> filament_2;

                        filament_2.p[0][0] = actinVec[actinID].getPoint(actinSubID)[0];
                        filament_2.p[1][0] = actinVec[actinID].getPoint(actinSubID+1)[0];
                        filament_2.p[0][1] = actinVec[actinID].getPoint(actinSubID)[1];
                        filament_2.p[1][1] = actinVec[actinID].getPoint(actinSubID+1)[1];


                        RobustQuery min_d_GTE;
                        auto result = min_d_GTE(filament_1, filament_2);
                        double distance = result.distance;

                        if (distance < (m_stericRadius + actinVec[actinID].getStericRadius()))
                        {
                                return true;
                        }
                }
        }
        return false;
}
bool Actin::check_min_dist_Diff_Grid(int subid, const std::vector<Actin> &actinVec,
                                     std::vector<int> &cellsToCheck,
                                     StericGrid &stericGrid, const bool diffBool)
{
        gte::Segment<2,double> filament_1;
        filament_1.p[0][0] = m_points[subid][0];
        filament_1.p[1][0] = m_points[subid+1][0];
        filament_1.p[0][1] = m_points[subid][1];
        filament_1.p[1][1] = m_points[subid+1][1];


        // To avoid repeated checks against the same subunits that may appear in
        // more than one cell, have a recorded of which ones we have checked
        std::vector<std::array<int,2> > checkedSubs;


        for (unsigned int i = 0; i < cellsToCheck.size(); ++i)
        {

                std::vector<std::array<int,2> > cellContents = stericGrid.getCellContents(cellsToCheck[i]);

                for (unsigned int j = 0; j < cellContents.size(); ++j)
                {
                        int actinID = cellContents[j][0];
                        int actinSubID = cellContents[j][1];

                        if (std::find(checkedSubs.begin(), checkedSubs.end(), cellContents[j]) != checkedSubs.end())
                        {
                                // We have checked this filament before!
                                continue;
                        }
                        checkedSubs.push_back(cellContents[j]);

                          // Brownian dynamics is on
                        if (actinVec[actinID].getSubStructureID() == m_subStructure_id || actinVec[actinID].getID() == m_parent_id || actinVec[actinID].getParentID() == m_id)
                        {
                                // Ignore the same fila, and any branches and the mother
                                continue;
                        }

                        // For rotation of a crosslink
                        bool cont = false;
                        for (unsigned int k = 0; k < m_cLinkActinAndSites.size(); ++k)
                        {
                            if (m_cLinkActinAndSites[k][2] == actinID && (subid >= m_cLinkActinAndSites[k][1]-1 && subid <= m_cLinkActinAndSites[k][1]+1))
                            {
                                // if the other sub has the link on too!
                                if (m_cLinkActinAndSites[k][5] == actinSubID || m_cLinkActinAndSites[k][5] == actinSubID-1 || m_cLinkActinAndSites[k][5] == actinSubID+1)
                                {
                                    cont = true;
                                    break;
                                }
                            }
                        }
                        if (cont)
                            continue;

                        gte::Segment<2,double> filament_2;

                        filament_2.p[0][0] = actinVec[actinID].getPoint(actinSubID)[0];
                        filament_2.p[1][0] = actinVec[actinID].getPoint(actinSubID+1)[0];
                        filament_2.p[0][1] = actinVec[actinID].getPoint(actinSubID)[1];
                        filament_2.p[1][1] = actinVec[actinID].getPoint(actinSubID+1)[1];


                        RobustQuery min_d_GTE;
                        auto result = min_d_GTE(filament_1, filament_2);
                        double distance = result.distance;
                        if (distance < (m_stericRadius + actinVec[actinID].getStericRadius()))
                        {
                                return true;
                        }
                }
        }
        return false;
}

bool Actin::check_min_dist_Diff_Grid2(int subid, const std::vector<Actin> &actinVec,
                                     std::vector<int> &cellsToCheck,
                                     StericGrid &stericGrid, const bool diffBool,
                                     std::vector<int> &checkedFilas)
{
        gte::Segment<2,double> filament_1;
        filament_1.p[0][0] = m_points[subid][0];
        filament_1.p[1][0] = m_points[subid+1][0];
        filament_1.p[0][1] = m_points[subid][1];
        filament_1.p[1][1] = m_points[subid+1][1];


        // To avoid repeated checks against the same subunits that may appear in
        // more than one cell, have a recorded of which ones we have checked
        std::vector<std::array<int,2> > checkedSubs;


        for (unsigned int i = 0; i < cellsToCheck.size(); ++i)
        {

                std::vector<std::array<int,2> > cellContents = stericGrid.getCellContents(cellsToCheck[i]);

                for (unsigned int j = 0; j < cellContents.size(); ++j)
                {

                        if (std::find(checkedSubs.begin(), checkedSubs.end(), cellContents[j]) != checkedSubs.end())
                        {
                                // We have checked this filament before!
                                continue;
                        }
                        checkedSubs.push_back(cellContents[j]);

                        if (std::find(checkedFilas.begin(), checkedFilas.end(), cellContents[j][0]) != checkedFilas.end())
                        {
                            // check this fila before in the structure!
                              continue;
                        }

                        int actinID = cellContents[j][0];
                        int actinSubID = cellContents[j][1];


                          // Brownian dynamics is on
                        if (actinVec[actinID].getSubStructureID() == m_subStructure_id || actinVec[actinID].getID() == m_parent_id || actinVec[actinID].getParentID() == m_id)
                        {
                                // Ignore the same fila, and any branches and the mother
                                continue;
                        }

                        // For rotation of a crosslink
                        bool cont = false;
                        for (unsigned int k = 0; k < m_cLinkActinAndSites.size(); ++k)
                        {
                            if (m_cLinkActinAndSites[k][2] == actinID && (subid >= m_cLinkActinAndSites[k][1]-1 && subid <= m_cLinkActinAndSites[k][1]+1))
                            {
                                // if the other sub has the link on too!
                                if (m_cLinkActinAndSites[k][5] == actinSubID || m_cLinkActinAndSites[k][5] == actinSubID-1 || m_cLinkActinAndSites[k][5] == actinSubID+1)
                                {
                                    cont = true;
                                    break;
                                }
                            }
                        }
                        if (cont)
                            continue;

                        gte::Segment<2,double> filament_2;

                        filament_2.p[0][0] = actinVec[actinID].getPoint(actinSubID)[0];
                        filament_2.p[1][0] = actinVec[actinID].getPoint(actinSubID+1)[0];
                        filament_2.p[0][1] = actinVec[actinID].getPoint(actinSubID)[1];
                        filament_2.p[1][1] = actinVec[actinID].getPoint(actinSubID+1)[1];


                        RobustQuery min_d_GTE;
                        auto result = min_d_GTE(filament_1, filament_2);
                        double distance = result.distance;
                        if (distance < (m_stericRadius + actinVec[actinID].getStericRadius()))
                        {
                                return true;
                        }
                }
        }
        return false;
}

bool Actin::check_min_dist_BD_Grid(int subid, const std::vector<Actin> &actinVec,
                                   std::vector<int> &cellsToCheck,
                                   StericGrid &stericGrid,
                                   std::vector<std::array<int,2> > &checkedFilaSubs)
{

        gte::Segment<2,double> filament_1;
        filament_1.p[0][0] = m_points[subid][0];
        filament_1.p[1][0] = m_points[subid+1][0];
        filament_1.p[0][1] = m_points[subid][1];
        filament_1.p[1][1] = m_points[subid+1][1];

        // To avoid repeated checks against the same subunits that may appear in
        // more than one cell, have a recorded of which ones we have checked
        std::vector<std::array<int,2> > checkedSubs;

        for (unsigned int i = 0; i < cellsToCheck.size(); ++i)
        {

                std::vector<std::array<int,2> > cellContents = stericGrid.getCellContents(cellsToCheck[i]);
                int contentsSize = cellContents.size();
                for (int j = 0; j < contentsSize; ++j)
                {

                        if (std::find(checkedSubs.begin(), checkedSubs.end(), cellContents[j]) != checkedSubs.end())
                        {
                                // We have checked this filament before!
                                continue;
                        }
                        checkedSubs.push_back(cellContents[j]);

                        if (std::find(checkedFilaSubs.begin(), checkedFilaSubs.end(), cellContents[j]) != checkedFilaSubs.end())
                        {
                                // This sub filament has previously been checked against all neighbours
                                // checking against this filament would be a duplicate check
                                continue;
                        }

                        int actinID = cellContents[j][0];
                        int actinSubID = cellContents[j][1];

                        if (actinID == m_id)
                        {
                          // Same filament, ignore same sub, and adjacent subs
                            if (actinSubID == subid || actinSubID == subid + 1 || actinSubID == subid-1)
                            {
                                continue;
                            }

                            if (!m_flexible)
                            {
                                // Rigid filament, cant hinder self
                                continue;
                            }
                        }

                        if (actinID == m_parent_id)
                        {
                            // other filament is parent, ignore sub that the branch sits on AND ADJACENT SUBS
                            if (actinSubID == m_branchSubUnit || actinSubID == m_branchSubUnit-1 || actinSubID == m_branchSubUnit+1)
                            {
                                continue;
                            }


                            if (m_length < 2*s_segmentationLength && actinVec[actinID].getLength() < 2*s_segmentationLength)
                            {
                                // both mother and daughter rigid, cant hinder
                                continue;
                            }

                        }

                        if (actinVec[actinID].getParentID() == m_id)
                        {

                            // other filament is a daughter, ignore first two segements on branch OR ALL OF IT IF RIGID!
                            if (m_length < 2*s_segmentationLength && actinVec[actinID].getLength() < 2*s_segmentationLength)
                            {
                                // both mother and daughter rigid, cant hinder
                                continue;
                            }

                            if (subid == actinVec[actinID].getBranchSubunit() || subid == actinVec[actinID].getBranchSubunit()-1 || subid == actinVec[actinID].getBranchSubunit()+1)
                            {
                                continue;
                            }

                        }
                        bool cont = false;
                        for (unsigned int k = 0; k < m_cLinkActinAndSites.size(); ++k)
                        {
                          // if other filament is linked to this one
                            if (m_cLinkActinAndSites[k][2] == actinID)
                            {
                                // if link is on this sub or adjacent sub
                                if (m_cLinkActinAndSites[k][1] == subid || m_cLinkActinAndSites[k][1] == subid-1 || m_cLinkActinAndSites[k][1] == subid+1)
                                {
                                    // if the other sub has the link on too!
                                    if (m_cLinkActinAndSites[k][5] == actinSubID || m_cLinkActinAndSites[k][5] == actinSubID-1 || m_cLinkActinAndSites[k][5] == actinSubID+1)
                                    {
                                        cont = true;
                                        break;
                                    }
                                }
                            }
                        }
                        if (cont)
                            continue;

                        gte::Segment<2,double> filament_2;

                        filament_2.p[0][0] = actinVec[actinID].getPoint(actinSubID)[0];
                        filament_2.p[1][0] = actinVec[actinID].getPoint(actinSubID+1)[0];
                        filament_2.p[0][1] = actinVec[actinID].getPoint(actinSubID)[1];
                        filament_2.p[1][1] = actinVec[actinID].getPoint(actinSubID+1)[1];


                        RobustQuery min_d_GTE;
                        auto result = min_d_GTE(filament_1, filament_2);
                        double distance = result.distance;

                        if (distance < (m_stericRadius + actinVec[actinID].getStericRadius()))
                        {
                                  return true;
                        }
                }
        }

        // add {m_id, subid} to checkedFilas here!
        std::array<int,2> checkedFilaSub = {m_id, subid};
        checkedFilaSubs.push_back(checkedFilaSub);

        return false;
}

void Actin::diffuse2D(double dxprime, double dyprime, std::vector<Actin> &actinVec,
                      Eigen::Matrix3d &rotationMatrix,
                      const std::vector<ExcZone> &excZones,
                      const std::vector<MembraneWall> &memWalls,
                      const std::vector<Membrane> &membranes,
                      Cortex &cortex,
                      const bool steric, StericGrid &stericGrid,
                      const bool translation,
                      bool crossLinking)
{
        // Moves the filament in 2D
        // Need to check whether we can move here
        // This is called to a parent filament only

        if (m_daughter_num == 0 && getNumCLinks() == 0)
        {
                // One filament, no branches so we create a rotation matrix here
                double dx { (dxprime * m_subunit_unitVecs[0][0])
                            - (dyprime * m_subunit_unitVecs[0][1]) };

                double dy { (dxprime * m_subunit_unitVecs[0][1])
                            + (dyprime * m_subunit_unitVecs[0][0]) };

                moveactin(dx, dy);

                if (steric)
                {
                        for (int i = 0; i < m_num_subunits; ++i)
                        {
                                // update the steric grid for the test seed
                                stericGrid.resetSubandUpdate(*this, i);
                        }
                }
                std::vector<int> dummy;
                if (check_diffusion_excVolGrid(actinVec, excZones, memWalls,
                                               membranes, cortex, steric,
                                               stericGrid, translation))
                {
                        // If we have violated exc area we move back
                        moveactin(-dx,-dy);

                        // Rejected, update back
                        for (int i = 0; i < m_num_subunits; ++i)
                        {
                                // update the steric grid for the test seed
                                stericGrid.resetSubandUpdate(*this, i);
                        }
                }

        }
        else
        {
                // Part of a structure

                Eigen::Vector3d primes;
                Eigen::Vector3d standards;

                primes (0) = dxprime;
                primes (1) = dyprime;
                primes (2) = 0.0; // Z prime for now lets make it 0

                standards = rotationMatrix * primes;

                double dx { standards(0) };
                double dy { standards(1) };

                movebranches(dx, dy, actinVec, translation);
                for (Actin &actinobj : actinVec)
                {

                    if (actinobj.getStructureID() == m_structure_id)
                    {
                            // Update the grid BEFORE CHECKING!?
                            if (steric)
                            {
                                    for (int i = 0; i < actinobj.getNumSubs(); ++i)
                                    {
                                            // update the steric grid for the test seed
                                            stericGrid.resetSubandUpdate(actinobj, i);
                                    }
                            }
                    }
                }

                if (!checkDiffusionStructure(actinVec, excZones, memWalls,
                                            membranes, cortex, steric,
                                            stericGrid, translation))
                {
                        movebranches(-dx, -dy, actinVec, translation);
                        if (steric)
                        {
                            for (Actin &actinobj : actinVec)
                            {
                                if (actinobj.getStructureID() == m_structure_id)
                                {
                                    // The move is rejected, so update the stericGrid back
                                    // have to reset and update for all subs - maybe able to optimise this
                                    for (int i = 0; i < actinobj.getNumSubs(); ++i)
                                    {
                                            // update the steric grid for the test seed
                                            stericGrid.resetSubandUpdate(actinobj, i);
                                    }
                                }
                            }
                        }
                }
                else
                {
                    // Update centre point
                    m_EllipCentrePoint[0] += dx;
                    m_EllipCentrePoint[1] += dy;
                }

        }


}

void Actin::movebranches(double dx, double dy, std::vector<Actin> &actinVec,
                         const bool translation)
{
        // Moves all the points in the structure by the same amount
        for (Actin &actinobj : actinVec)
        {
            // Brownian dynamics, only sub-structure moves
              if (m_structure_id == actinobj.getStructureID())
              {
                  actinobj.moveactin(dx, dy);
              }
        }
}

void Actin::moveactin(double dx, double dy)
{
        for (int i = 0; i <= m_num_subunits; ++i)
        {
                m_points[i][0] += dx;
                m_points[i][1] += dy;
        }

}

void Actin::detach(std::vector<Actin> &actinVec, int nActin)
{
        // Function that dislocates the filament from its parent
        // The actin it is acting on is a branch filament
        // Instead of loop DO THIS actinVec[m_parent_id]

        // find the subunit the branch is on

        actinVec[m_parent_id].purgeBranchFromDaughterVec(m_id);
        actinVec[m_parent_id].decrementDaughternum();
        actinVec[m_parent_id].setChangedBoolTrue(); // we changed the old mothers bool

        m_changed = true; // and we changed the freed branch too
        m_parent_id = -1;
        setPointedUncapped(); // Uncaps the pointed end


        // Need to change the centre now that it is no longer anchored

        // Calculate new global angle
        // angle
        updateUnitVecs();


        m_structure_id = s_structure_idGenerator++;
        m_subStructure_id = s_subStructure_idGenerator++;
        m_crossStructure_id = s_crossStructure_idGenerator++;



        // Also must change any attached daughter filaments to this structure id

        changeStructure(actinVec, nActin);
        changeSubStructure(actinVec);
        changeCLStructure(actinVec); // only change ones that are branches

}

void Actin::purgeBranchFromDaughterVec(int branchid)
{
        // called on the parent
        // this removes the branch from the parent's vector of vectors
        // there is probably a more efficient way to do this
        for (int i = 0; i < m_num_subunits; ++i)
        {
                int vecsize = m_daughterIDs[i].size();
                for (int j = 0; j < vecsize; ++j)
                {
                        if (m_daughterIDs[i][j] == branchid)
                        {

                                // found the branch, its located at i,j
                                // need to remove this from the vector
                                for (int k = j; k < vecsize - 1; ++k)
                                {
                                        // overwrite downwards
                                        m_daughterIDs[i][k] = m_daughterIDs[i][k+1];
                                }
                                m_daughterIDs[i].pop_back();
                                break;
                        }
                }
        }
}

void Actin::changeStructure(std::vector<Actin> &actinVec, int nActin)
{
        // function that finds the all the daughters AND PARENTS of the filament that has
        // debranched and sets their structureID to that of the parent
        std::vector<int> changedFilas;
        changedFilas.push_back(m_id); // this call!

        for (int i = 0; i < nActin; ++i)
        {
                if (m_id == actinVec[i].getParentID() || actinVec[i].getID() == m_parent_id)
                {
                        actinVec[i].setStructureID(m_structure_id);
                        changedFilas.push_back(i);
                        // Need to use recursion here for grandaughters etc
                        actinVec[i].changeStructure(actinVec, nActin,
                                                    changedFilas);
                }
        }

        for (int i = 0; i < getNumCLinks(); ++i)
        {
            int otherFila = m_cLinkActinAndSites[i][2];
            actinVec[otherFila].setStructureID(m_structure_id);
            changedFilas.push_back(otherFila);
            actinVec[otherFila].changeStructure(actinVec, nActin, changedFilas);
        }

}

void Actin::changeStructure(std::vector<Actin> &actinVec, int nActin,
                            std::vector<int> &changedFilas)
{
        // function that finds the all the daughters of the filament that has
        // debranched and sets their structureID to that of the parent

        for (int i = 0; i < nActin; ++i)
        {
                if (m_id == actinVec[i].getParentID() || actinVec[i].getID() == m_parent_id)
                {
                      if (std::find(changedFilas.begin(), changedFilas.end(), i) != changedFilas.end())
                      {
                              // We have changed this filament before!
                              continue;
                      }
                      actinVec[i].setStructureID(m_structure_id);
                      changedFilas.push_back(i);
                      // Need to use recursion here for grandaughters etc
                      actinVec[i].changeStructure(actinVec, nActin, changedFilas);
                }
        }

        for (int i = 0; i < getNumCLinks(); ++i)
        {
            int otherFila = m_cLinkActinAndSites[i][2];
            if (std::find(changedFilas.begin(), changedFilas.end(), otherFila) != changedFilas.end())
            {
                    // We have changed this filament before!
                    continue;
            }

            actinVec[otherFila].setStructureID(m_structure_id);
            changedFilas.push_back(otherFila);
            actinVec[otherFila].changeStructure(actinVec, nActin, changedFilas);
        }

}

void Actin::changeStructureALL(std::vector<Actin> &actinVec, int newStructureID)
{
        // function that finds all other filaments that belong to the
        // same structure and change them

        int oldStructureID = m_structure_id;

        for (unsigned int i = 0; i < actinVec.size(); ++i)
        {
            if (actinVec[i].getStructureID() == oldStructureID)
            {
                actinVec[i].setStructureID(newStructureID);
            }
        }

}

void Actin::changeSubStructure(std::vector<Actin> &actinVec)
{
        // after debranching, change sub structure accordingly
        // called after changeStructure!

        for (int i = 0; i < m_num_subunits; ++i)
        {
            for (unsigned int j = 0; j < m_daughterIDs[i].size(); ++j)
            {
                int daughtID = m_daughterIDs[i][j];
                if (actinVec[daughtID].getLength() < 2*s_segmentationLength)
                {
                    actinVec[daughtID].setSubStructureID(m_subStructure_id);
                    actinVec[daughtID].changeSubStructure(actinVec);
                }
            }
        }
}

void Actin::changeCLStructure(std::vector<Actin> &actinVec)
{
        // function that finds the all the daughters of the filament that has
        // debranched and sets their structureID to that of the parent


        for (unsigned int i = 0; i < actinVec.size(); ++i)
        {
                if (m_id == actinVec[i].getParentID())
                {
                        actinVec[i].setCLStructureID(m_crossStructure_id);
                        // Need to use recursion here for grandaughters etc
                        actinVec[i].changeCLStructure(actinVec);
                }
        }


}


std::array<double,3> Actin::fitEllipsoids(const std::vector<Actin> &actinVec,
                                          const double viscosity, const double KT,
                                          Eigen::Matrix3d &rotationMatrix,
                                          std::array<double,3> &centrePoint,
                                          const bool translation,
                                          const bool rotation)

{
        // Passed to this function is a `mother` filament, not an unbranched one!

        // efficiency: this should only be called once every few timesteps
        // OR at least only called when the structure has changed!

        std::array<double, 3> D_array {};
        // D_array will contain the three translational diffusion coefficients,
        // followed by the three rotational diffusion coefficients

        std::vector< std::array<double,3> > actin_points {  };

        for (unsigned int i = 0; i < actinVec.size(); ++i)
        {
                if (actinVec[i].getStructureID() == m_structure_id)
                {
                        int nMonomers = actinVec[i].getNumMonomers();
            #pragma omp parallel
                        {
                                std::vector< std::array<double,3> > actin_points_private { };
                #pragma omp for
                                for (int j = 0; j <= nMonomers; j += 2)
                                {

                                        double x_loc { actinVec[i].getPointedEnd()[0]
                                                       + actinVec[i].getUnitVec(0)[0]*s_monomerLength*j };


                                        double y_loc { actinVec[i].getPointedEnd()[1]
                                                       + actinVec[i].getUnitVec(0)[1]*s_monomerLength*j };
                                        std::array<double,3> tmp { x_loc, y_loc, 0 };
                                        actin_points_private.push_back(tmp);
                                }
                #pragma omp critical
                                {

                                        actin_points.insert(actin_points.end(), actin_points_private.begin(), actin_points_private.end());
                                }
                        }
                }
        }
        std::sort(actin_points.begin(),actin_points.end());

        std::array<double,3> radii = ellipsoidfit(actin_points, 2, rotationMatrix, centrePoint);

        double a = radii[0];
        double b = radii[1];
        double c = radii[2];
        //std::cout << b << std::endl;
        // Step 2: Find P, Q, R, S
        a = (a < m_radius) ? m_radius : a;
        b = (b < m_radius) ? m_radius : b;
        //c = (c < m_radius) ? m_radius : c;

        // 2D adjustment of third semi-axis
        c = m_radius;

        double S = PerrinSolverS(a, b, c);
        double P = PerrinSolverPQR(a, b, c);
        double Q = PerrinSolverPQR(b, a, c);
        //double R = PerrinSolverPQR(c, a, b);


        if (rotation)
        {
                D_array[2] = KT*rotEllipMobility(viscosity, a, b, c, P, Q);

                //D_array[3] = rotMobility[0] * KT; // around x?
                //D_array[4] = rotMobility[1] * KT; // around y?
                //D_array[5] = rotMobility[2] * KT; // around z?
        }

        if (translation)
        {
                std::array<double,2> transMobility {};
                transMobility = transEllipMobility(viscosity, a, b, c, P, Q, S);
                D_array[0] = transMobility[0] * KT; // x
                D_array[1] = transMobility[1] * KT; // y
                //D_array[2] = transMobility[2] * KT; // z
        }

        return D_array; // also returns by reference the rotationMatrix and centrePoint

}

std::array<double,3> Actin::fitEllipsoidsCLink(const std::vector<Actin> &actinVec,
                                          const double viscosity, const double KT,
                                          Eigen::Matrix3d &rotationMatrix,
                                          std::array<double,3> &centrePoint,
                                          const bool translation,
                                          const bool rotation)

{
        // Passed to this function is a `mother` filament, not an unbranched one!

        // efficiency: this should only be called once every few timesteps
        // OR at least only called when the structure has changed!

        std::array<double, 3> D_array {};
        // D_array will contain the three translational diffusion coefficients,
        // followed by the three rotational diffusion coefficients

        std::vector< std::array<double,3> > actin_points {  };
        for (unsigned int i = 0; i < actinVec.size(); ++i)
        {
                if (actinVec[i].getCLStructureID() == m_crossStructure_id)
                {
                        int nMonomers = actinVec[i].getNumMonomers();
            #pragma omp parallel
                        {
                                std::vector< std::array<double,3> > actin_points_private { };
                #pragma omp for
                                for (int j = 0; j <= nMonomers; j += 2)
                                {

                                        double x_loc { actinVec[i].getPointedEnd()[0]
                                                       + actinVec[i].getUnitVec(0)[0]*s_monomerLength*j };


                                        double y_loc { actinVec[i].getPointedEnd()[1]
                                                       + actinVec[i].getUnitVec(0)[1]*s_monomerLength*j };
                                        std::array<double,3> tmp { x_loc, y_loc, 0 };
                                        actin_points_private.push_back(tmp);
                                }
                #pragma omp critical
                                {

                                        actin_points.insert(actin_points.end(), actin_points_private.begin(), actin_points_private.end());
                                }
                        }
                }
        }
        std::sort(actin_points.begin(),actin_points.end());

        std::array<double,3> radii = ellipsoidfit(actin_points, 2, rotationMatrix, centrePoint);

        double a = radii[0];
        double b = radii[1];
        double c = radii[2];
        //std::cout << b << std::endl;
        // Step 2: Find P, Q, R, S
        a = (a < m_radius) ? m_radius : a;
        b = (b < m_radius) ? m_radius : b;
        c = (c < m_radius) ? m_radius : c;

        // 2D adjustment of third semi-axis
        c = m_radius;

        double S = PerrinSolverS(a, b, c);
        double P = PerrinSolverPQR(a, b, c);
        double Q = PerrinSolverPQR(b, a, c);
        //double R = PerrinSolverPQR(c, a, b);


        if (rotation)
        {
                D_array[2] = KT*rotEllipMobility(viscosity, a, b, c, P, Q);

                //D_array[3] = rotMobility[0] * KT; // around x?
                //D_array[4] = rotMobility[1] * KT; // around y?
                //D_array[5] = rotMobility[2] * KT; // around z?
        }

        if (translation)
        {
                std::array<double,2> transMobility {};
                transMobility = transEllipMobility(viscosity, a, b, c, P, Q, S);
                D_array[0] = transMobility[0] * KT; // x
                D_array[1] = transMobility[1] * KT; // y
                //D_array[2] = transMobility[2] * KT; // z
        }

        return D_array; // also returns by reference the rotationMatrix and centrePoint

}

std::array<double,2> Actin::transEllipMobility(const double viscosity, double a,
                                               double b, double c, double P,
                                               double Q, double S)
{

        std::array<double,2> mobility_array {};
        mobility_array[0] = ( (S + ((a * a) * P))
                              / (16 * M_PI * viscosity) );

        mobility_array[1] = ( (S + ((b * b) * Q))
                              / (16 * M_PI * viscosity) );

        //mobility_array[2] = ( (S + ((c * c) * R))
        //                      / (16 * M_PI * viscosity) );

        return mobility_array;
}

double Actin::rotEllipMobility(const double viscosity, double a,
                                              double b, double c, double P,
                                              double Q)
{
        //std::array<double,3> mobility_array {};

        //mobility_array[0] = ( 3*((b*b)*Q + (c*c)*R)
        //                      / ((16*M_PI*viscosity) * ((b*b)+(c*c))) );

        //mobility_array[1] = ( 3*((c*c)*R + (a*a)*P)
        //                      / ((16*M_PI*viscosity) * ((c*c)+(a*a))) );

        double mobility_z = ( 3*((a*a)*P + (b*b)*Q)
                              / ((16*M_PI*viscosity) * ((a*a)+(b*b))) );

        return mobility_z;

}

std::array<double,3> Actin::kirkwoodDiff(const double viscosity,
                                         const double KT, const bool translation,
                                         const bool rotation)
{
        // called on a single rod
        std::array<double, 3> D_array {};
        // again, this will contain diffusion coefficients for both translation
        // and rotation

        if (translation)
        {
                D_array[0] = KT * x_filamentmobility(viscosity);
                D_array[1] = KT* y_filamentmobility(viscosity);
                //D_array[2] = 0;
        }

        if (rotation)
        {
                //D_array[3] = 0;
                //D_array[4] = 0; // 2D so only one angle
                D_array[2] = KT * rodDRot(viscosity);
        }

        return D_array;
}

double Actin::x_filamentmobility(const double viscosity)
{

        double Yama_corr_X = 0.044;
        double mobility = ( (log(m_length / (2*m_radius)) + Yama_corr_X)
                            / (2 * M_PI * viscosity * m_length) );

        assert(mobility > 0);
        return mobility;
}

double Actin::y_filamentmobility(const double viscosity)
{

        // This is currently the perpendicular mobility from Kirkwood shish kebab model

        double Yama_corr_Y = 1.111;
        double mobility = ( (log(m_length / (2*m_radius)) + Yama_corr_Y)
                            / (4 * M_PI * viscosity * m_length) );
        assert(mobility > 0);
        return mobility;
}

double Actin::rodDRot(const double viscosity)
{

        // This is the mobility to rotate about the z axis
        double Yama_corr_rot = -0.8;

        double rotMobility = ( (3 * (log(m_length / (2*m_radius)) + Yama_corr_rot))
                               / (M_PI * viscosity * m_length * m_length * m_length) );
        if (rotMobility < 0)
        {
                rotMobility = ( (3 * (log(m_length / (2*m_radius))))
                                / (M_PI * viscosity * m_length * m_length * m_length) );
        }

        assert (rotMobility > 0);
        return rotMobility;
}

double Actin::rodDRotTether(const double viscosity)
{

        // This is the mobility to rotate about the z axis

        // For this we take the length to be twice the longer of the two free ends,
        // as the rod will rotate about the tether point, not the centre of mass

        double effL = (m_distToTetBarb > m_distToTetPoint) ? 2*m_distToTetBarb*s_monomerLength : 2*m_distToTetPoint*s_monomerLength;
        effL = (effL < 2*m_radius) ? s_seedSize*s_monomerLength : effL;

        double rotMobility = ( (3 * (log(effL / (2*m_radius)) - 0.8))
                               / (M_PI * viscosity * effL * effL * effL) );

        if (rotMobility < 0)
        {
                rotMobility = ( (3 * (log(effL / (2*m_radius))))
                                / (M_PI * viscosity * effL * effL * effL) );
        }

        assert (rotMobility > 0);
        return rotMobility;
}

void Actin::rotate2D(double d_xyangle, std::vector<Actin> &actinVec,
                     const std::vector<ExcZone> &excZones,
                     const std::vector<MembraneWall> &memWalls,
                     const std::vector<Membrane> &membranes,
                     Cortex &cortex,
                     const bool steric, StericGrid &stericGrid,
                     const bool rotation)
{
        // for an unbranched single filament

        // Rotation from the centre point

        std::array<double,2> centrePoint;
        int centreSub = 0; // gets calculate below and returned by reference
        centrePoint = findCentrePoint(centreSub);

        Eigen::Matrix2d rotMatrix;
        rotMatrix << cos(d_xyangle), -sin(d_xyangle),
                sin(d_xyangle), cos(d_xyangle);

        Eigen::Vector2d oldPoints;
        oldPoints << m_points[centreSub][0] - centrePoint[0], m_points[centreSub][1] - centrePoint[1];

        Eigen::Vector2d newPoints;
        newPoints = rotMatrix*oldPoints;
        m_points[centreSub][0] = newPoints(0) + centrePoint[0];
        m_points[centreSub][1] = newPoints(1) + centrePoint[1];

        for (unsigned int i = 0; i < m_subunit_unitVecs.size(); ++i)
        {
                Eigen::Vector2d unitVec;
                unitVec << m_subunit_unitVecs[i][0], m_subunit_unitVecs[i][1];
                unitVec = rotMatrix*unitVec;

                m_subunit_unitVecs[i][0] = unitVec(0);
                m_subunit_unitVecs[i][1] = unitVec(1);
        }

        m_points[centreSub+1][0] = m_points[centreSub][0] + (m_actualSubLengths[centreSub] * m_subunit_unitVecs[centreSub][0]);
        m_points[centreSub+1][1] = m_points[centreSub][1] + (m_actualSubLengths[centreSub] * m_subunit_unitVecs[centreSub][1]);

        updatePointsUp(centreSub+1);
        updatePointsDown(centreSub);
        // Check excluded volume
        // Update the grid BEFORE CHECKING!?
        if (steric)
        {
                for (int i = 0; i < m_num_subunits; ++i)
                {
                        // update the steric grid for the test seed
                        stericGrid.resetSubandUpdate(*this, i);
                }
        }

        if (check_diffusion_excVolGrid(actinVec, excZones, memWalls, membranes,
                                       cortex, steric, stericGrid, rotation))
        {
                // If we have violated exc area we move back

                rotMatrix << cos(-d_xyangle), -sin(-d_xyangle),
                        sin(-d_xyangle), cos(-d_xyangle);


                for (unsigned int i = 0; i < m_subunit_unitVecs.size(); ++i)
                {

                        Eigen::Vector2d unitVec;
                        unitVec << m_subunit_unitVecs[i][0], m_subunit_unitVecs[i][1];
                        unitVec = rotMatrix*unitVec;
                        m_subunit_unitVecs[i][0] = unitVec(0);
                        m_subunit_unitVecs[i][1] = unitVec(1);
                }

                m_points[centreSub][0] = oldPoints(0) + centrePoint[0];
                m_points[centreSub][1] = oldPoints(1) + centrePoint[1];
                m_points[centreSub+1][0] = m_points[centreSub][0] + (m_actualSubLengths[centreSub] * m_subunit_unitVecs[centreSub][0]);
                m_points[centreSub+1][1] = m_points[centreSub][1] + (m_actualSubLengths[centreSub] * m_subunit_unitVecs[centreSub][1]);

                updatePointsUp(centreSub+1);
                updatePointsDown(centreSub);

                // Rejected, update back
                for (int i = 0; i < m_num_subunits; ++i)
                {
                        // update the steric grid for the test seed
                        stericGrid.resetSubandUpdate(*this, i);
                }
        }


}

void Actin::rotate2DTether(double d_xyangle, std::vector<Actin> &actinVec,
                           const std::vector<ExcZone> &excZones,
                           const std::vector<MembraneWall> &memWalls,
                           const std::vector<Membrane> &membranes,
                           Cortex &cortex,
                           const bool steric, StericGrid &stericGrid,
                           const bool rotation)
{
        /*
           Given a "short" filament and ONE tether point, the filament will rotate around
           the tether point
         */

        // NEW: Rotation from the Tether point
        assert(getNumTethers() == 1);
        int tetherSub = m_tetherSubunits[0];
        std::array<double,2> centrePoint;
        centrePoint[0] = m_tetherPoints[0][0];
        centrePoint[1] = m_tetherPoints[0][1];

        Eigen::Matrix2d rotMatrix;
        rotMatrix << cos(d_xyangle), -sin(d_xyangle),
                sin(d_xyangle), cos(d_xyangle);

        Eigen::Vector2d oldPoints;
        oldPoints << m_points[tetherSub][0] - centrePoint[0], m_points[tetherSub][1] - centrePoint[1];
        Eigen::Vector2d newPoints;
        newPoints = rotMatrix*oldPoints;
        m_points[tetherSub][0] = newPoints(0) + centrePoint[0];
        m_points[tetherSub][1] = newPoints(1) + centrePoint[1];

        for (unsigned int i = 0; i < m_subunit_unitVecs.size(); ++i)
        {
                Eigen::Vector2d unitVec;
                unitVec << m_subunit_unitVecs[i][0], m_subunit_unitVecs[i][1];
                unitVec = rotMatrix*unitVec;
                m_subunit_unitVecs[i][0] = unitVec(0);
                m_subunit_unitVecs[i][1] = unitVec(1);
        }

        m_points[tetherSub+1][0] = m_points[tetherSub][0] + (m_actualSubLengths[tetherSub] * m_subunit_unitVecs[tetherSub][0]);
        m_points[tetherSub+1][1] = m_points[tetherSub][1] + (m_actualSubLengths[tetherSub] * m_subunit_unitVecs[tetherSub][1]);

        updatePointsUp(tetherSub+1);
        updatePointsDown(tetherSub);

        // Check excluded volume
        // Update the grid BEFORE CHECKING!?
        if (steric)
        {
                for (int i = 0; i < m_num_subunits; ++i)
                {
                        // update the steric grid for the test seed
                        stericGrid.resetSubandUpdate(*this, i);
                }
        }

        if (check_diffusion_excVolGrid(actinVec, excZones, memWalls, membranes,
                                       cortex, steric, stericGrid, rotation))
        {
                // If we have violated exc area we move back

                rotMatrix << cos(-d_xyangle), -sin(-d_xyangle),
                        sin(-d_xyangle), cos(-d_xyangle);


                for (unsigned int i = 0; i < m_subunit_unitVecs.size(); ++i)
                {
                        Eigen::Vector2d unitVec;
                        unitVec << m_subunit_unitVecs[i][0], m_subunit_unitVecs[i][1];
                        unitVec = rotMatrix*unitVec;
                        m_subunit_unitVecs[i][0] = unitVec(0);
                        m_subunit_unitVecs[i][1] = unitVec(1);
                }

                m_points[tetherSub][0] = oldPoints(0) + centrePoint[0];
                m_points[tetherSub][1] = oldPoints(1) + centrePoint[1];
                m_points[tetherSub+1][0] = m_points[tetherSub][0] + (m_actualSubLengths[tetherSub] * m_subunit_unitVecs[tetherSub][0]);
                m_points[tetherSub+1][1] = m_points[tetherSub][1] + (m_actualSubLengths[tetherSub] * m_subunit_unitVecs[tetherSub][1]);

                updatePointsUp(tetherSub+1);
                updatePointsDown(tetherSub);

                // Rejected, update back
                for (int i = 0; i < m_num_subunits; ++i)
                {
                        // update the steric grid for the test seed
                        stericGrid.resetSubandUpdate(*this, i);
                }
        }

}

void Actin::rotate2DCLink(double d_xyangle, std::vector<Actin> &actinVec,
                           const std::vector<ExcZone> &excZones,
                           const std::vector<MembraneWall> &memWalls,
                           const std::vector<Membrane> &membranes,
                           Cortex &cortex,
                           const bool steric, StericGrid &stericGrid,
                           const bool rotation)
{
        /*
           Given a "short" filament and ONE crosslink, the filament will rotate around
           the crosslink
         */

        // NEW: Rotation from the crosslink
        assert(getNumCLinks() == 1);
        int cLinkSub = m_cLinkActinAndSites[0][1];
        std::array<double,2> centrePoint; // centre of rotation
        centrePoint[0] = m_cLinkPoints[0][0]; // rotates around crosslink
        centrePoint[1] = m_cLinkPoints[0][1];

        Eigen::Matrix2d rotMatrix;
        rotMatrix << cos(d_xyangle), -sin(d_xyangle),
                sin(d_xyangle), cos(d_xyangle);

        Eigen::Vector2d oldPoints;
        oldPoints << m_points[cLinkSub][0] - centrePoint[0], m_points[cLinkSub][1] - centrePoint[1];
        Eigen::Vector2d newPoints;
        newPoints = rotMatrix*oldPoints;
        m_points[cLinkSub][0] = newPoints(0) + centrePoint[0];
        m_points[cLinkSub][1] = newPoints(1) + centrePoint[1];

        for (unsigned int i = 0; i < m_subunit_unitVecs.size(); ++i)
        {
                Eigen::Vector2d unitVec;
                unitVec << m_subunit_unitVecs[i][0], m_subunit_unitVecs[i][1];
                unitVec = rotMatrix*unitVec;
                m_subunit_unitVecs[i][0] = unitVec(0);
                m_subunit_unitVecs[i][1] = unitVec(1);
        }

        m_points[cLinkSub+1][0] = m_points[cLinkSub][0] + (m_actualSubLengths[cLinkSub] * m_subunit_unitVecs[cLinkSub][0]);
        m_points[cLinkSub+1][1] = m_points[cLinkSub][1] + (m_actualSubLengths[cLinkSub] * m_subunit_unitVecs[cLinkSub][1]);

        updatePointsUp(cLinkSub+1);
        updatePointsDown(cLinkSub);

        // Check excluded volume
        // Update the grid BEFORE CHECKING!?
        if (steric)
        {
                for (int i = 0; i < m_num_subunits; ++i)
                {
                        // update the steric grid for the test seed
                        stericGrid.resetSubandUpdate(*this, i);
                }
        }

        if (check_diffusion_excVolGrid(actinVec, excZones, memWalls, membranes,
                                       cortex, steric, stericGrid, rotation))
        {
                // If we have violated exc area we move back
                rotMatrix << cos(-d_xyangle), -sin(-d_xyangle),
                        sin(-d_xyangle), cos(-d_xyangle);


                for (unsigned int i = 0; i < m_subunit_unitVecs.size(); ++i)
                {
                        Eigen::Vector2d unitVec;
                        unitVec << m_subunit_unitVecs[i][0], m_subunit_unitVecs[i][1];
                        unitVec = rotMatrix*unitVec;
                        m_subunit_unitVecs[i][0] = unitVec(0);
                        m_subunit_unitVecs[i][1] = unitVec(1);
                }

                m_points[cLinkSub][0] = oldPoints(0) + centrePoint[0];
                m_points[cLinkSub][1] = oldPoints(1) + centrePoint[1];
                m_points[cLinkSub+1][0] = m_points[cLinkSub][0] + (m_actualSubLengths[cLinkSub] * m_subunit_unitVecs[cLinkSub][0]);
                m_points[cLinkSub+1][1] = m_points[cLinkSub][1] + (m_actualSubLengths[cLinkSub] * m_subunit_unitVecs[cLinkSub][1]);

                updatePointsUp(cLinkSub+1);
                updatePointsDown(cLinkSub);

                // Rejected, update back
                for (int i = 0; i < m_num_subunits; ++i)
                {
                        // update the steric grid for the test seed
                        stericGrid.resetSubandUpdate(*this, i);
                }
        }
        else
        {
            m_changed = true; // rotation will change the shape of the ellipse
        }

}


std::array<double,2> Actin::findCentrePoint(int &centreSubUnit)
{
        // Function that returns the centre point of a single filament, whether bent
        // or straight

        double lengthToCentre = m_length / 2;
        //int centreSubUnit; // not necessarily m_centre_subunit!
        for (unsigned int i = 0; i < m_actualSubLengths.size(); ++i)
        {
                lengthToCentre -= m_actualSubLengths[i];
                if (lengthToCentre <= 0)
                {
                        lengthToCentre += m_actualSubLengths[i];
                        centreSubUnit = i;
                        // length to Centre is now positive and is the length from the
                        // point previous to the centre, to the centre.
                        break;
                }

        }

        // Calculate new centre subunit
        // dist is the distance from the subunits "pointed end" to the sever point
        std::array<double,2> centrePoint;
        centrePoint[0] = m_points[centreSubUnit][0] + (m_subunit_unitVecs[centreSubUnit][0]*lengthToCentre);
        centrePoint[1] = m_points[centreSubUnit][1] + (m_subunit_unitVecs[centreSubUnit][1]*lengthToCentre);

        return centrePoint;

}

std::array<double,2> Actin::findCentrePoint(int &centreSubUnit) const
{
        // Function that returns the centre point of a single filament, whether bent
        // or straight

        double lengthToCentre = m_length / 2;

        //int centreSubUnit; // not necessarily m_centre_subunit!
        for (unsigned int i = 0; i < m_actualSubLengths.size(); ++i)
        {
                lengthToCentre -= m_actualSubLengths[i];
                if (lengthToCentre <= 0)
                {
                        lengthToCentre += m_actualSubLengths[i];
                        centreSubUnit = i;
                        // length to Centre is now positive and is the length from the
                        // point previous to the centre, to the centre.
                        break;
                }
        }

        // Calculate new centre subunit
        // dist is the distance from the subunits "pointed end" to the sever point
        std::array<double,2> centrePoint;
        centrePoint[0] = m_points[centreSubUnit][0] + (m_subunit_unitVecs[centreSubUnit][0]*lengthToCentre);
        centrePoint[1] = m_points[centreSubUnit][1] + (m_subunit_unitVecs[centreSubUnit][1]*lengthToCentre);

        return centrePoint;

}

void Actin::rotateEllipsoid2D(double d_xyangle, std::vector<Actin> &actinVec,
                              Eigen::Matrix3d &rotationMatrix,
                              std::array<double,3> &centrePoint,
                              const std::vector<ExcZone> &excZones,
                              const std::vector<MembraneWall> &memWalls,
                              const std::vector<Membrane> &membranes,
                              Cortex &cortex,
                              const bool steric, StericGrid &stericGrid,
                              const bool rotation)
{
        // the mother filament of a branched structure is passed to the function
        // this will change the rotationmatrix
        // actually it may be simpler than that, dont need to worry about rotating
        // the ellipsoid, just rotate the mother and all daughters

        // DO WE NEED TO CHANGE THE ROTATION MATRIX IF DOING efficiency???
        // YES!!!!

        // rotate the mother from the centre point of the ellipsoid!
        // centrePoint is passed to this function
        int centreSub = 0; // centre sub of mother gets calculate below and returned by reference
        findCentrePoint(centreSub);


        Eigen::Matrix2d rotMatrix;
        rotMatrix << cos(d_xyangle), -sin(d_xyangle),
                sin(d_xyangle), cos(d_xyangle);

        Eigen::Vector2d oldPoints;
        oldPoints << m_points[centreSub][0] - centrePoint[0], m_points[centreSub][1] - centrePoint[1];
        Eigen::Vector2d newPoints;
        newPoints = rotMatrix*oldPoints;
        m_points[centreSub][0] = newPoints(0) + centrePoint[0];
        m_points[centreSub][1] = newPoints(1) + centrePoint[1];

        for (unsigned int i = 0; i < m_subunit_unitVecs.size(); ++i)
        {
                Eigen::Vector2d unitVec;
                unitVec << m_subunit_unitVecs[i][0], m_subunit_unitVecs[i][1];
                unitVec = rotMatrix*unitVec;
                m_subunit_unitVecs[i][0] = unitVec(0);
                m_subunit_unitVecs[i][1] = unitVec(1);
        }

        m_points[centreSub+1][0] = m_points[centreSub][0] + (m_actualSubLengths[centreSub] * m_subunit_unitVecs[centreSub][0]);
        m_points[centreSub+1][1] = m_points[centreSub][1] + (m_actualSubLengths[centreSub] * m_subunit_unitVecs[centreSub][1]);

        updatePointsDown(centreSub);
        updatePointsUp(centreSub+1);

        saveBranchPosInclFlex(actinVec, rotation);
        rotbranchesInclFlex(actinVec, rotation);

        assert(!checkForCrossLinksInStructure(actinVec));
        // Check no crosslinks

        // Check if violated excluded volume?

        for (Actin &actinobj : actinVec)
        {
                if (actinobj.getStructureID() == m_structure_id)
                {
                        // Lets try something here
                        // Update the grid BEFORE CHECKING!?
                        if (steric)
                        {
                                for (int i = 0; i < actinobj.getNumSubs(); ++i)
                                {
                                        // update the steric grid for the test seed
                                        stericGrid.resetSubandUpdate(actinobj, i);
                                }
                        }
                }
        }

        if (!checkDiffusionStructure(actinVec, excZones, memWalls,
                                          membranes, cortex, steric, stericGrid,
                                          rotation))
        {
            // Need to rotate back

            rotMatrix << cos(-d_xyangle), -sin(-d_xyangle),
                    sin(-d_xyangle), cos(-d_xyangle);


            for (unsigned int i = 0; i < m_subunit_unitVecs.size(); ++i)
            {
                    Eigen::Vector2d unitVec;
                    unitVec << m_subunit_unitVecs[i][0], m_subunit_unitVecs[i][1];
                    unitVec = rotMatrix*unitVec;
                    m_subunit_unitVecs[i][0] = unitVec(0);
                    m_subunit_unitVecs[i][1] = unitVec(1);
            }

            // If we have violated exc area we move back
            m_points[centreSub][0] = oldPoints(0) + centrePoint[0];
            m_points[centreSub][1] = oldPoints(1) + centrePoint[1];
            m_points[centreSub+1][0] = m_points[centreSub][0] + (m_actualSubLengths[centreSub] * m_subunit_unitVecs[centreSub][0]);
            m_points[centreSub+1][1] = m_points[centreSub][1] + (m_actualSubLengths[centreSub] * m_subunit_unitVecs[centreSub][1]);
            updatePointsDown(centreSub);
            updatePointsUp(centreSub+1);

            //rotbranches(actinVec, rotation);
            resetBranchPosInclFlex(actinVec, rotation);

            // Update grid back to what it was
            if (steric)
            {
                    for (Actin &actinobj : actinVec)
                    {
                            if (actinobj.getStructureID() == m_structure_id)
                            {
                                    // The move is accepted, so update the stericGrid
                                    // have to reset and update for all subs - maybe able to optimise this
                                    for (int i = 0; i < actinobj.getNumSubs(); ++i)
                                    {
                                            // update the steric grid for the test seed
                                            stericGrid.resetSubandUpdate(actinobj, i);
                                    }
                            }
                    }
            }
        }
        else
        {
            // change rotation MATRIX

                Eigen::Matrix3d rotMatrix3D;
                rotMatrix3D << cos(d_xyangle), -sin(d_xyangle), 0,
                               sin(d_xyangle), cos(d_xyangle), 0,
                               0, 0, 1;
                rotationMatrix = rotationMatrix * rotMatrix3D;

        }


}

void Actin::rotateEllipsoid2DTether(double d_xyangle, std::vector<Actin> &actinVec,
                                    Eigen::Matrix3d &rotationMatrix,
                                    const std::vector<ExcZone> &excZones,
                                    const std::vector<MembraneWall> &memWalls,
                                    const std::vector<Membrane> &membranes,
                                    Cortex &cortex,
                                    const bool steric, StericGrid &stericGrid,
                                    const bool rotation)
{
        // the mother filament of a branched structure is passed to the function
        // this will change the rotationmatrix
        // actually it may be simpler than that, dont need to worry about rotating
        // the ellipsoid, just rotate the mother and all daughters

        // DO WE NEED TO CHANGE THE ROTATION MATRIX IF DOING efficiency???
        // YES!!!!

        // rotate the mother from the centre point of the ellipsoid!
        std::array<double,2> centrePoint;
        centrePoint[0] = m_tetherPoints[0][0];
        centrePoint[1] = m_tetherPoints[0][1];
        int tetherSub = m_tetherSubunits[0];

        Eigen::Matrix2d rotMatrix;
        rotMatrix << cos(d_xyangle), -sin(d_xyangle),
                sin(d_xyangle), cos(d_xyangle);

        Eigen::Vector2d oldPoints;
        oldPoints << m_points[tetherSub][0] - centrePoint[0], m_points[tetherSub][1] - centrePoint[1];
        Eigen::Vector2d newPoints;
        newPoints = rotMatrix*oldPoints;
        m_points[tetherSub][0] = newPoints(0) + centrePoint[0];
        m_points[tetherSub][1] = newPoints(1) + centrePoint[1];

        for (unsigned int i = 0; i < m_subunit_unitVecs.size(); ++i)
        {
                Eigen::Vector2d unitVec;
                unitVec << m_subunit_unitVecs[i][0], m_subunit_unitVecs[i][1];
                unitVec = rotMatrix*unitVec;
                m_subunit_unitVecs[i][0] = unitVec(0);
                m_subunit_unitVecs[i][1] = unitVec(1);
        }

        m_points[tetherSub+1][0] = m_points[tetherSub][0] + (m_actualSubLengths[tetherSub] * m_subunit_unitVecs[tetherSub][0]);
        m_points[tetherSub+1][1] = m_points[tetherSub][1] + (m_actualSubLengths[tetherSub] * m_subunit_unitVecs[tetherSub][1]);

        updatePointsDown(tetherSub);
        updatePointsUp(tetherSub+1);

        saveBranchPosInclFlex(actinVec, rotation);
        rotbranchesInclFlex(actinVec, rotation);

        // Check if violated excluded volume?
        for (Actin &actinobj : actinVec)
        {
                if (actinobj.getStructureID() == m_structure_id)
                {
                        // Update the grid BEFORE CHECKING!?
                        if (steric)
                        {
                                for (int i = 0; i < actinobj.getNumSubs(); ++i)
                                {
                                        // update the steric grid for the test seed
                                        stericGrid.resetSubandUpdate(actinobj, i);
                                }
                        }
                }
      }

      if (!checkDiffusionStructure(actinVec, excZones,  memWalls, membranes,
                                  cortex, steric, stericGrid, rotation))
      {
          // If we have violated exc area we move back

          rotMatrix << cos(-d_xyangle), -sin(-d_xyangle),
                  sin(-d_xyangle), cos(-d_xyangle);


          for (unsigned int i = 0; i < m_subunit_unitVecs.size(); ++i)
          {
                  Eigen::Vector2d unitVec;
                  unitVec << m_subunit_unitVecs[i][0], m_subunit_unitVecs[i][1];
                  unitVec = rotMatrix*unitVec;
                  m_subunit_unitVecs[i][0] = unitVec(0);
                  m_subunit_unitVecs[i][1] = unitVec(1);
          }

          m_points[tetherSub][0] = oldPoints(0) + centrePoint[0];
          m_points[tetherSub][1] = oldPoints(1) + centrePoint[1];
          m_points[tetherSub+1][0] = m_points[tetherSub][0] + (m_actualSubLengths[tetherSub] * m_subunit_unitVecs[tetherSub][0]);
          m_points[tetherSub+1][1] = m_points[tetherSub][1] + (m_actualSubLengths[tetherSub] * m_subunit_unitVecs[tetherSub][1]);

          updatePointsDown(tetherSub);
          updatePointsUp(tetherSub+1);

          resetBranchPosInclFlex(actinVec, rotation);

          if (steric)
          {
                  for (Actin &actinobj : actinVec)
                  {
                          if (actinobj.getStructureID() == m_structure_id)
                          {
                                  // The move is accepted, so update the stericGrid
                                  // have to reset and update for all subs - maybe able to optimise this
                                  for (int i = 0; i < actinobj.getNumSubs(); ++i)
                                  {
                                          // update the steric grid for the test seed
                                          stericGrid.resetSubandUpdate(actinobj, i);
                                  }
                          }
                  }
          }
        }
        else
        {
                // change rotation MATRIX
                rotationMatrix (0, 0) = m_subunit_unitVecs[tetherSub][0];
                rotationMatrix (0, 1) = -m_subunit_unitVecs[tetherSub][1];
                rotationMatrix (1, 0) = m_subunit_unitVecs[tetherSub][1];
                rotationMatrix (1, 1) = m_subunit_unitVecs[tetherSub][0];

                rotationMatrix (0, 2) = 0;
                rotationMatrix (1, 2) = 0;
                rotationMatrix (2, 0) = 0;
                rotationMatrix (2, 1) = 0;
                rotationMatrix (2, 2) = 1;

        }

}

void Actin::rotateEllipsoidALL(double d_xyangle, std::vector<Actin> &actinVec,
                                   Eigen::Matrix3d &rotationMatrix,
                                   std::array<double,3> &centrePoint,
                                   const std::vector<ExcZone> &excZones,
                                   const std::vector<MembraneWall> &memWalls,
                                   const std::vector<Membrane> &membranes,
                                   Cortex &cortex,
                                   const bool steric, StericGrid &stericGrid,
                                   const bool rotation)
{
        // the mother filament of a branched structure is passed to the function
        // this will change the rotationmatrix
        // actually it may be simpler than that, dont need to worry about rotating
        // the ellipsoid, just rotate the mother and all daughters

        // DO WE NEED TO CHANGE THE ROTATION MATRIX IF DOING efficiency???
        // YES!!!!


        rotateAllStructure(actinVec, d_xyangle, centrePoint);

        // Check if violated excluded volume?

        for (Actin &actinobj : actinVec)
        {
            if (actinobj.getStructureID() == m_structure_id)
            {
                // Lets try something here
                // Update the grid BEFORE CHECKING!?
                if (steric)
                {
                    for (int i = 0; i < actinobj.getNumSubs(); ++i)
                    {
                            // update the steric grid for the test seed
                            stericGrid.resetSubandUpdate(actinobj, i);
                    }
                }
            }
        }

        if (!checkDiffusionStructure(actinVec, excZones, memWalls, membranes,
                                    cortex, steric, stericGrid, rotation))
        {
            rotateAllStructure(actinVec, -d_xyangle, centrePoint);
            // Update grid back to what it was
            if (steric)
            {
                for (Actin &actinobj : actinVec)
                {
                    if (actinobj.getStructureID() == m_structure_id)
                    {
                        // The move is rejected, so update the stericGrid
                        // have to reset and update for all subs - maybe able to optimise this
                        for (int i = 0; i < actinobj.getNumSubs(); ++i)
                        {
                                // update the steric grid for the test seed
                                stericGrid.resetSubandUpdate(actinobj, i);
                        }
                    }
                }
            }
        }
        else
        {
            // change rotation MATRIX
            Eigen::Matrix3d rotMatrix3D;
            rotMatrix3D << cos(d_xyangle), -sin(d_xyangle), 0,
                           sin(d_xyangle), cos(d_xyangle), 0,
                           0, 0, 1;
            rotationMatrix = rotationMatrix * rotMatrix3D;

        }

}

void Actin::rotateMiniEllipsoid(double d_xyangle, std::vector<Actin> &actinVec,
                                   Eigen::Matrix3d &rotationMatrix,
                                   std::array<double,3> &centrePoint,
                                   const std::vector<ExcZone> &excZones,
                                   const std::vector<MembraneWall> &memWalls,
                                   const std::vector<Membrane> &membranes,
                                   Cortex &cortex,
                                   const bool steric, StericGrid &stericGrid,
                                   const bool rotation)
{
        // the mother filament of a branched structure is passed to the function
        // this will change the rotationmatrix
        // actually it may be simpler than that, dont need to worry about rotating
        // the ellipsoid, just rotate the mother and all daughters

        // DO WE NEED TO CHANGE THE ROTATION MATRIX IF DOING efficiency???
        // YES!!!!


        std::array<double,3> rotPoint = getCLCoordInCLStruct(actinVec);
        // rot point not necessary on m_id filament, could be on a branch


        rotateCLStructure(actinVec, d_xyangle, rotPoint);
        // Check if violated excluded volume?

        for (Actin &actinobj : actinVec)
        {
            if (actinobj.getCLStructureID() == m_crossStructure_id)
            {
                // Lets try something here
                // Update the grid BEFORE CHECKING!?
                if (steric)
                {
                    for (int i = 0; i < actinobj.getNumSubs(); ++i)
                    {
                            // update the steric grid for the test seed
                            stericGrid.resetSubandUpdate(actinobj, i);
                    }
                }
            }
        }

        if (!checkDiffusionClStructure(actinVec, excZones, memWalls, membranes,
                                       cortex, steric, stericGrid, rotation))
        {
            rotateCLStructure(actinVec, -d_xyangle, rotPoint);
            if (steric)
            {
                for (Actin &actinobj : actinVec)
                {
                    if (actinobj.getCLStructureID() == m_crossStructure_id)
                    {
                        // The move is accepted, so update the stericGrid
                        // have to reset and update for all subs - maybe able to optimise this
                        for (int i = 0; i < actinobj.getNumSubs(); ++i)
                        {
                              // update the steric grid for the test seed
                              stericGrid.resetSubandUpdate(actinobj, i);
                        }
                    }
                }
            }
        }

}

void Actin::rotateAllStructure(std::vector<Actin> &actinVec, double d_xyangle,
                               std::array<double,3> &centrePoint)
{
    // Function that rotates the entire structure

    Eigen::Matrix2d rotMatrix;
    rotMatrix << cos(d_xyangle), -sin(d_xyangle),
            sin(d_xyangle), cos(d_xyangle);

    int centreSub = 0;
    Eigen::Vector2d oldPoints;
    for (unsigned int i = 0; i < actinVec.size(); ++i)
    {
        if (actinVec[i].getStructureID() == m_structure_id)
        {
            centreSub = 0; // centre sub of mother gets calculate below and returned by reference
            actinVec[i].findCentrePoint(centreSub);
            //Eigen::Vector2d oldPoints;
            oldPoints << actinVec[i].getPoint(centreSub)[0] - centrePoint[0], actinVec[i].getPoint(centreSub)[1] - centrePoint[1];
            Eigen::Vector2d newPoints;
            newPoints = rotMatrix*oldPoints;
            // NEED TO SET POINTS!
            std::array<double,2> newPoint = {newPoints(0) + centrePoint[0], newPoints(1) + centrePoint[1]};
            actinVec[i].setPoints(centreSub, newPoint);

            for (unsigned int j = 0; j < actinVec[i].getUnitVecs().size(); ++j)
            {
                    Eigen::Vector2d unitVec;
                    unitVec << actinVec[i].getUnitVec(j)[0], actinVec[i].getUnitVec(j)[1];
                    unitVec = rotMatrix*unitVec;
                    std::array<double,2> unitVec2 = {unitVec(0), unitVec(1)};
                    actinVec[i].setUnitVec(j, unitVec2);
            }


            double newPointp1x = actinVec[i].getPoint(centreSub)[0] + (actinVec[i].getActualSubLengths()[centreSub]*actinVec[i].getUnitVec(centreSub)[0]);
            double newPointp1y = actinVec[i].getPoint(centreSub)[1] + (actinVec[i].getActualSubLengths()[centreSub]*actinVec[i].getUnitVec(centreSub)[1]);

            std::array<double,2> newPointp1 = {newPointp1x, newPointp1y};
            actinVec[i].setPoints(centreSub+1, newPointp1);
            actinVec[i].updatePointsDown(centreSub);
            actinVec[i].updatePointsUp(centreSub+1);
        }
    }
}

void Actin::rotateCLStructure(std::vector<Actin> &actinVec, double d_xyangle,
                               std::array<double,3> &centrePoint)
{
    // Function that rotates the entire structure

    Eigen::Matrix2d rotMatrix;
    rotMatrix << cos(d_xyangle), -sin(d_xyangle),
            sin(d_xyangle), cos(d_xyangle);

    int centreSub = 0;
    Eigen::Vector2d oldPoints;
    for (unsigned int i = 0; i < actinVec.size(); ++i)
    {
        if (actinVec[i].getCLStructureID() == m_crossStructure_id)
        {
            centreSub = 0; // centre sub of mother gets calculate below and returned by reference
            actinVec[i].findCentrePoint(centreSub);
            oldPoints << actinVec[i].getPoint(centreSub)[0] - centrePoint[0], actinVec[i].getPoint(centreSub)[1] - centrePoint[1];
            Eigen::Vector2d newPoints;
            newPoints = rotMatrix*oldPoints;
            // NEED TO SET POINTS!
            std::array<double,2> newPoint = {newPoints(0) + centrePoint[0], newPoints(1) + centrePoint[1]};
            actinVec[i].setPoints(centreSub, newPoint);

            for (unsigned int j = 0; j < actinVec[i].getUnitVecs().size(); ++j)
            {
                    Eigen::Vector2d unitVec;
                    unitVec << actinVec[i].getUnitVec(j)[0], actinVec[i].getUnitVec(j)[1];
                    unitVec = rotMatrix*unitVec;
                    std::array<double,2> unitVec2 = {unitVec(0), unitVec(1)};
                    actinVec[i].setUnitVec(j, unitVec2);
            }


            double newPointp1x = actinVec[i].getPoint(centreSub)[0] + (actinVec[i].getActualSubLengths()[centreSub]*actinVec[i].getUnitVec(centreSub)[0]);
            double newPointp1y = actinVec[i].getPoint(centreSub)[1] + (actinVec[i].getActualSubLengths()[centreSub]*actinVec[i].getUnitVec(centreSub)[1]);

            std::array<double,2> newPointp1 = {newPointp1x, newPointp1y};
            actinVec[i].setPoints(centreSub+1, newPointp1);
            actinVec[i].updatePointsDown(centreSub);
            actinVec[i].updatePointsUp(centreSub+1);
        }
    }
}

void Actin::rotbranches(std::vector<Actin> &actinVec, const bool rotation)
{

        // function that rotates all branches called during
        // bending

        // called on the parent which has already had its angle and position updated

        for (int i = 0; i < m_num_subunits; ++i)
        {
                // For each subunit we need to find all the branches
                for (int j : m_daughterIDs[i])
                {
                        // Calculate new pointed end of the branch and new global angle
                        // j is now the id of the branched filament
                        // m_id is the id of the mother filament
                        // i is the mother's subunit ID

                        if (actinVec[j].getLength() >= 2*s_segmentationLength)
                        {
                                // avoid doing this as it will use bending alg
                                continue;
                        }


                        double lenAlongSub;
                        int subid = findSubunit(actinVec[j].getMotherMonoID(), lenAlongSub);

                        assert(subid==i);

                        double specialBranchSubPos = lenAlongSub / m_prescribedSubLengths[i];

                        assert(specialBranchSubPos < 1);
                        double newxpos = (m_points[i][0]
                                          + (specialBranchSubPos
                                             * (m_points[i+1][0]
                                                - m_points[i][0])));

                        double newypos = (m_points[i][1]
                                          + (specialBranchSubPos
                                             * (m_points[i+1][1]
                                                - m_points[i][1])));

                        actinVec[j].setPointedEndPoints(newxpos, newypos);

                        // need to update the global angle

                        // NEED TO CHANGE UNIT VECS HERE
                        // GET UNIT VEC OF PARENT AND APPLY ROT MATRIX
                        std::array<double,2> parUnitVec = m_subunit_unitVecs[subid];
                        Eigen::Matrix2d rotMatrix;
                        int branchDir = actinVec[j].getBranchDir();
                        rotMatrix << cos(branchDir*s_branchAngle), -sin(branchDir*s_branchAngle),
                                sin(branchDir*s_branchAngle), cos(branchDir*s_branchAngle);

                        Eigen::Vector2d parentVecEig;
                        parentVecEig << parUnitVec[0], parUnitVec[1];
                        Eigen::Vector2d branchVecEig = rotMatrix*parentVecEig;
                        std::array<double,2> branchVec;
                        branchVec[0] = branchVecEig(0);
                        branchVec[1] = branchVecEig(1);
                        for (int k = 0; k < actinVec[j].getNumSubs(); ++k)
                        {
                                actinVec[j].setUnitVec(k, branchVec);
                        }

                        // update branch filament up from pointed end
                        actinVec[j].updatePointsUp(0);

                        actinVec[j].updateSubLengths();
                        actinVec[j].updateUnitVecs();

                        // recursion, call this on the branch
                        actinVec[j].rotbranches(actinVec, rotation);

                        actinVec[j].updateTetherLoc();
                }
        }

}

void Actin::rotRigidCLinked(std::vector<Actin> &actinVec, const bool rotation)
{
        // original
        std::vector<int> checkedFilas;
        checkedFilas.push_back(m_id);

        // function that rotates all short rigid crosslinks called during
        // brownian dynamics

        for (int i = 0; i < m_num_subunits; ++i)
        {
                // For each subunit we need to find all the rigid crosslinks

                for (int j = 0; j < getNumCLinks(); ++j)
                {
                    int otherID = m_cLinkActinAndSites[j][2];
                    if (actinVec[otherID].getCLStructureID() == m_crossStructure_id)
                    {
                        continue; // already accounted for with branches?
                    }

                    if (actinVec[otherID].getLength() >= 2*s_segmentationLength)
                    {
                        // avoid as will move on its own
                        continue;
                    }


                    if (i == m_cLinkActinAndSites[j][1])
                    {
                        if (std::find(checkedFilas.begin(), checkedFilas.end(), otherID) != checkedFilas.end())
                        {
                                // We have checked this filament before!
                                continue;
                        }
                        checkedFilas.push_back(otherID);

                        // Correct subunit

                        double oldxpos = m_cLinkPoints[j][0];
                        double oldypos = m_cLinkPoints[j][1];

                        updateCrossLinkLoc(actinVec);

                        double newxpos = m_cLinkPoints[j][0];
                        double newypos = m_cLinkPoints[j][1];

                        double dx = newxpos - oldxpos;
                        double dy = newypos - oldypos;

                        actinVec[otherID].moveactin(dx, dy);

                        // Check any short crosslinkins off this crosslink
                        // Dont need to do this?
                        actinVec[otherID].rotRigidCLinked(actinVec, rotation,
                                                          checkedFilas);

                        actinVec[otherID].rotbranches(actinVec, rotation);

                        actinVec[otherID].updateTetherLoc();
                        actinVec[otherID].updateCrossLinkLoc(actinVec);
                    }
                }

                // Also check rigid daughters for crosslinks
                for (int j : m_daughterIDs[i])
                {
                        if (actinVec[j].getLength() >= 2*s_segmentationLength)
                        {
                                // avoid doing this as it will use bending alg
                                continue;
                        }

                        // recursion, call this on the branch
                        actinVec[j].rotRigidCLinked(actinVec, rotation,
                                                    checkedFilas);
                }
        }

}

void Actin::rotRigidCLinked(std::vector<Actin> &actinVec, const bool rotation,
                            std::vector<int> &checkedFilas)
{
        // overloaded, carries on with checkedFilas
        // function that rotates all short rigid crosslinks called during
        // brownian dynamics

        for (int i = 0; i < m_num_subunits; ++i)
        {
                // For each subunit we need to find all the rigid crosslinks

                for (int j = 0; j < getNumCLinks(); ++j)
                {
                    int otherID = m_cLinkActinAndSites[j][2];

                    if (actinVec[otherID].getCLStructureID() == m_crossStructure_id)
                    {
                        continue; // already accounted for with branches?
                    }

                    if (actinVec[otherID].getLength() >= 2*s_segmentationLength)
                    {
                        // avoid as will move on its own
                        continue;
                    }

                    if (i == m_cLinkActinAndSites[j][1])
                    {
                        if (std::find(checkedFilas.begin(), checkedFilas.end(), otherID) != checkedFilas.end())
                        {
                                // We have checked this filament before!
                                continue;
                        }
                        checkedFilas.push_back(otherID);
                        // Correct subunit


                        double oldxpos = m_cLinkPoints[j][0];
                        double oldypos = m_cLinkPoints[j][1];

                        updateCrossLinkLoc(actinVec);

                        double newxpos = m_cLinkPoints[j][0];
                        double newypos = m_cLinkPoints[j][1];

                        double dx = newxpos - oldxpos;
                        double dy = newypos - oldypos;
                        actinVec[otherID].moveactin(dx, dy);

                        // Check any short crosslinkins off this crosslink
                        // Dont need to do this?
                        actinVec[otherID].rotRigidCLinked(actinVec, rotation,
                                                          checkedFilas);

                        actinVec[otherID].rotbranches(actinVec, rotation);

                        actinVec[otherID].updateTetherLoc();
                        actinVec[otherID].updateCrossLinkLoc(actinVec);

                    }

                }

                // Also check rigid daughters for crosslinks
                for (int j : m_daughterIDs[i])
                {
                        if (actinVec[j].getLength() >= 2*s_segmentationLength)
                        {
                            // avoid doing this as it will use bending alg
                            continue;
                        }
                        // recursion, call this on the branch
                        actinVec[j].rotRigidCLinked(actinVec, rotation,
                                                    checkedFilas);
                }
        }

}

void Actin::rotbranchesInclFlex(std::vector<Actin> &actinVec, const bool rotation)
{

        // function that rotates all branches called during
        // bending#
        // Includes flexible branches so for ellipsoid rotation

        // called on the parent which has already had its angle and position updated

        for (int i = 0; i < m_num_subunits; ++i)
        {
                // For each subunit we need to find all the branches
                for (int j : m_daughterIDs[i])
                {
                        // Calculate new pointed end of the branch and new global angle
                        // j is now the id of the branched filament
                        // m_id is the id of the mother filament
                        // i is the mother's subunit ID

                        double lenAlongSub;
                        int subid = findSubunit(actinVec[j].getMotherMonoID(), lenAlongSub);

                        assert(subid==i);

                        double specialBranchSubPos = lenAlongSub / m_prescribedSubLengths[i];

                        assert(specialBranchSubPos < 1);
                        double newxpos = (m_points[i][0]
                                          + (specialBranchSubPos
                                             * (m_points[i+1][0]
                                                - m_points[i][0])));

                        double newypos = (m_points[i][1]
                                          + (specialBranchSubPos
                                             * (m_points[i+1][1]
                                                - m_points[i][1])));

                        actinVec[j].setPointedEndPoints(newxpos, newypos);

                        // need to update the global angle

                        // NEED TO CHANGE UNIT VECS HERE
                        // GET UNIT VEC OF PARENT AND APPLY ROT MATRIX
                        std::array<double,2> parUnitVec = m_subunit_unitVecs[subid];
                        Eigen::Matrix2d rotMatrix;
                        int branchDir = actinVec[j].getBranchDir();
                        rotMatrix << cos(branchDir*s_branchAngle), -sin(branchDir*s_branchAngle),
                                sin(branchDir*s_branchAngle), cos(branchDir*s_branchAngle);

                        Eigen::Vector2d parentVecEig;
                        parentVecEig << parUnitVec[0], parUnitVec[1];
                        Eigen::Vector2d branchVecEig = rotMatrix*parentVecEig;
                        std::array<double,2> branchVec;
                        branchVec[0] = branchVecEig(0);
                        branchVec[1] = branchVecEig(1);
                        for (int k = 0; k < actinVec[j].getNumSubs(); ++k)
                        {
                                actinVec[j].setUnitVec(k, branchVec);
                        }

                        // update branch filament up from pointed end
                        actinVec[j].updatePointsUp(0);

                        actinVec[j].updateSubLengths();
                        actinVec[j].updateUnitVecs();

                        // recursion, call this on the branch
                        actinVec[j].rotbranchesInclFlex(actinVec, rotation);

                        // Update tether positions on the branch?
                        actinVec[j].updateTetherLoc();
                        actinVec[j].updateCrossLinkLoc(actinVec);
                }
        }

}


void Actin::rotbranchesSub(std::vector<Actin> &actinVec, int p, const bool rotation)
{

        // function that rotates branches on a subunit level, called during
        // bending

        // called on the parent which has already had its angle and position updated

        for (int j : m_daughterIDs[p])
        {
                // Calculate new pointed end of the branch and new global angle
                // j is now the id of the branched filament
                // m_id is the id of the mother filament
                // i is the mother's subunit ID

                if (actinVec[j].getLength() >= 2*s_segmentationLength)
                {
                        // avoid doing this as it will use bending alg
                        continue;
                }

                double lenAlongSub;
                int subid = findSubunit(actinVec[j].getMotherMonoID(), lenAlongSub);

                assert(subid==p);

                double specialBranchSubPos = lenAlongSub / m_prescribedSubLengths[p];

                assert(specialBranchSubPos < 1);
                double newxpos = (m_points[p][0]
                                  + (specialBranchSubPos
                                     * (m_points[p+1][0]
                                        - m_points[p][0])));

                double newypos = (m_points[p][1]
                                  + (specialBranchSubPos
                                     * (m_points[p+1][1]
                                        - m_points[p][1])));

                actinVec[j].setPointedEndPoints(newxpos, newypos);

                // need to update the global angle

                // NEED TO CHANGE UNIT VECS HERE
                // GET UNIT VEC OF PARENT AND APPLY ROT MATRIX
                std::array<double,2> parUnitVec = m_subunit_unitVecs[subid];
                Eigen::Matrix2d rotMatrix;
                int branchDir = actinVec[j].getBranchDir();
                rotMatrix << cos(branchDir*s_branchAngle), -sin(branchDir*s_branchAngle),
                        sin(branchDir*s_branchAngle), cos(branchDir*s_branchAngle);

                Eigen::Vector2d parentVecEig;
                parentVecEig << parUnitVec[0], parUnitVec[1];
                Eigen::Vector2d branchVecEig = rotMatrix*parentVecEig;
                std::array<double,2> branchVec;
                branchVec[0] = branchVecEig(0);
                branchVec[1] = branchVecEig(1);
                for (int k = 0; k < actinVec[j].getNumSubs(); ++k)
                {
                        actinVec[j].setUnitVec(k, branchVec);
                }

                // update branch filament up from pointed end
                actinVec[j].updatePointsUp(0);

                actinVec[j].updateSubLengths();
                actinVec[j].updateUnitVecs();

                // recursion, call this on the branch
                actinVec[j].rotbranches(actinVec, rotation);
                actinVec[j].updateTetherLoc();
                actinVec[j].updateCrossLinkLoc(actinVec);
        }


}

void Actin::saveOldPos(std::vector<Actin> &actinVec)
{
    // Going to pass through a vector containing ids that the function has already been called on
    std::vector<int> calledFilas;
    saveOldPos(actinVec, calledFilas);
}

void Actin::saveOldPos(std::vector<Actin> &actinVec,
                       std::vector<int> &calledFilas)
{

    // Called after movement, saves old positions of any rigid branches
    // and crosslinks, so can reset if violated steric hindrance
    // first check if this function has been called on filament before
    if (std::find(calledFilas.begin(), calledFilas.end(), m_id) != calledFilas.end())
    {
        // We have done this filament already, so return
        return;
    }
    // Then add "filament" to calledFilas
    calledFilas.push_back(m_id);


    for (int i = 0; i < m_num_subunits; ++i)
    {
        for (int j : m_daughterIDs[i])
        {
            if (actinVec[j].getLength() < 2*s_segmentationLength)
            {
                actinVec[j].updateUnitVecs();
                actinVec[j].savePointsToPrev();
                actinVec[j].saveOldPos(actinVec, calledFilas);
            }
        }
    }

    for (int i = 0; i < getNumCLinks(); ++i)
    {
        // update the steric grid for the test seed
        int cLinkID = getCLinkActinAndSite(i)[2];
        if (actinVec[cLinkID].getLength() < 2*s_segmentationLength)
        {
          actinVec[cLinkID].updateUnitVecs();
          actinVec[cLinkID].savePointsToPrev();
          actinVec[cLinkID].saveOldPos(actinVec, calledFilas);
        }
    }

}

void Actin::resetToOldPos(std::vector<Actin> &actinVec)
{
    // Going to pass through a vector containing ids that the function has already been called on
    std::vector<int> calledFilas;
    resetToOldPos(actinVec, calledFilas);
}

void Actin::resetToOldPos(std::vector<Actin> &actinVec,
                          std::vector<int> &calledFilas)
{
    // Function that resets the positions of any rigid branches or crosslinks,
    // done after steric hindrance has been violated

    // first check if this function has been called on filament before
    if (std::find(calledFilas.begin(), calledFilas.end(), m_id) != calledFilas.end())
    {
        // We have done this filament already, so return
        return;
    }
    // Then add "filament" to calledFilas
    calledFilas.push_back(m_id);


    for (int i = 0; i < m_num_subunits; ++i)
    {
        for (int j : m_daughterIDs[i])
        {
            if (actinVec[j].getLength() < 2*s_segmentationLength)
            {
                actinVec[j].resetPointsToPrev();
                actinVec[j].updateSubLengths();
                actinVec[j].updateUnitVecs();
                actinVec[j].updateTetherLoc();
                actinVec[j].resetToOldPos(actinVec, calledFilas);
            }
        }
    }

    for (int i = 0; i < getNumCLinks(); ++i)
    {
        int cLinkID = getCLinkActinAndSite(i)[2];
        if (actinVec[cLinkID].getLength() < 2*s_segmentationLength)
        {
          actinVec[cLinkID].resetPointsToPrev();
          actinVec[cLinkID].updateSubLengths();
          actinVec[cLinkID].updateUnitVecs();
          actinVec[cLinkID].updateTetherLoc();
          actinVec[cLinkID].resetToOldPos(actinVec, calledFilas);
        }
    }
}

void Actin::saveBranchPos(std::vector<Actin> &actinVec, const bool rotation)
{

        // function that rotates branches on a subunit level, called during
        // bending

        // called on the parent which has already had its angle and position updated


        for (int i = 0; i < m_num_subunits; ++i)
        {
                // For each subunit we need to find all the branches
                for (int j : m_daughterIDs[i])
                {
                        // Calculate new pointed end of the branch and new global angle
                        // j is now the id of the branched filament
                        // m_id is the id of the mother filament
                        // i is the mother's subunit ID

                        if (actinVec[j].getLength() >= 2*s_segmentationLength)
                        {
                                // avoid doing this as it will use bending alg
                                continue;
                        }

                        actinVec[j].updateUnitVecs();
                        actinVec[j].savePointsToPrev();
                        actinVec[j].saveBranchPos(actinVec, rotation);
                }
        }

}

void Actin::saveBranchPosInclFlex(std::vector<Actin> &actinVec, const bool rotation)
{

        // function that rotates branches on a subunit level, called during
        // bending

        // called on the parent which has already had its angle and position updated

        for (int i = 0; i < m_num_subunits; ++i)
        {
                // For each subunit we need to find all the branches
                for (int j : m_daughterIDs[i])
                {
                        // Calculate new pointed end of the branch and new global angle
                        // j is now the id of the branched filament
                        // m_id is the id of the mother filament
                        // i is the mother's subunit ID

                        actinVec[j].updateUnitVecs();
                        actinVec[j].savePointsToPrev();
                        actinVec[j].saveBranchPosInclFlex(actinVec, rotation);
                }
        }

}

void Actin::saveBranchPosSub(std::vector<Actin> &actinVec, int p, const bool rotation)
{

        // Saves branch pos on subunit p


        for (int j : m_daughterIDs[p])
        {

                if (actinVec[j].getLength() >= 2*s_segmentationLength)
                {
                        // avoid doing this as it will use bending alg
                        continue;
                }

                actinVec[j].savePointsToPrev();
                actinVec[j].saveBranchPos(actinVec, rotation);
        }


}

void Actin::resetBranchPos(std::vector<Actin> &actinVec, const bool rotation)
{

        // function that rotates branches on a subunit level, called during
        // bending

        // called on the parent which has already had its angle and position updated

        for (int i = 0; i < m_num_subunits; ++i)
        {
                // For each subunit we need to find all the branches
                for (int j : m_daughterIDs[i])
                {
                        // Calculate new pointed end of the branch and new global angle
                        // j is now the id of the branched filament
                        // m_id is the id of the mother filament
                        // i is the mother's subunit ID

                        if (actinVec[j].getLength() >= 2*s_segmentationLength)
                        {
                                // avoid doing this as it will use bending alg
                                continue;
                        }


                        actinVec[j].resetPointsToPrev();

                        actinVec[j].updateSubLengths();

                        actinVec[j].updateUnitVecs();

                        actinVec[j].updateTetherLoc();
                        actinVec[j].resetBranchPos(actinVec, rotation);
                }
        }

}

void Actin::resetBranchPosInclFlex(std::vector<Actin> &actinVec, const bool rotation)
{

        // function that rotates branches on a subunit level, called during
        // bending

        // called on the parent which has already had its angle and position updated

        for (int i = 0; i < m_num_subunits; ++i)
        {
                // For each subunit we need to find all the branches
                for (int j : m_daughterIDs[i])
                {
                        // Calculate new pointed end of the branch and new global angle
                        // j is now the id of the branched filament
                        // m_id is the id of the mother filament
                        // i is the mother's subunit ID


                        actinVec[j].resetPointsToPrev();

                        actinVec[j].updateSubLengths();

                        actinVec[j].updateUnitVecs();

                        actinVec[j].updateTetherLoc();
                        actinVec[j].updateCrossLinkLoc(actinVec);
                        actinVec[j].resetBranchPosInclFlex(actinVec, rotation);
                }
        }

}

void Actin::resetBranchPosSub(std::vector<Actin> &actinVec, int p, const bool rotation)
{
        // Resets branch for subunit p

        for (int j : m_daughterIDs[p])
        {

                if (actinVec[j].getLength() >= 2*s_segmentationLength)
                {
                        // avoid doing this as it will use bending alg
                        continue;
                }
                actinVec[j].resetPointsToPrev();
                actinVec[j].updateSubLengths();
                actinVec[j].updateUnitVecs();
                actinVec[j].updateTetherLoc();
                actinVec[j].updateCrossLinkLoc(actinVec);
                actinVec[j].resetBranchPos(actinVec, rotation);
        }
}

Eigen::VectorXd Actin::calcBranchForces(std::vector<Actin> &actinVec, int N)
{

        // function that gets the hookean forces from daughters
        // called on the mother

        Eigen::VectorXd branchF (2*N);
        branchF.setZero();
        if (m_daughter_num == 0)
        {
            return branchF;
        }

        for (int i = 0; i < m_num_subunits; ++i)
        {
                // For each subunit we need to find all the branches
                for (int j : m_daughterIDs[i])
                {
                        // Hookean spring with non-zero rest length

                        // one end is the branchpoint on the mother
                        // other end is the first midpoint of the branch (second point)
                        assert(actinVec[j].getParentID() != -1);

                        if (actinVec[j].getLength() < 2*s_segmentationLength)
                        {
                                // Small branch, ignore the effect
                                continue;
                        }
                        double lenAlongSub;
                        int subid = findSubunit(actinVec[j].getMotherMonoID(), lenAlongSub);


                        assert(subid==i);

                        double specialBranchSubPos = lenAlongSub / m_prescribedSubLengths[i];
                        assert(specialBranchSubPos < 1);
                        double branchx = (m_points[i][0]
                                          + (specialBranchSubPos
                                             * (m_points[i+1][0]
                                                - m_points[i][0])));

                        double branchy = (m_points[i][1]
                                          + (specialBranchSubPos
                                             * (m_points[i+1][1]
                                                - m_points[i][1])));


                        if (i == 0)
                        {
                                //--------------------------------------------------------------

                                // Forces applied to the first midpoint

                                branchF(0) = branchF(0) - s_branchStiff*(branchx - actinVec[j].getPointedEnd()[0]);
                                branchF(1) = branchF(1) - s_branchStiff*(branchy - actinVec[j].getPointedEnd()[1]);

                                //--------------------------------------------------------------

                                // units are Nm/rad
                                // resting angle is +/- 1.22
                                // Need an  angular component which acts on the two points
                                // either end of the mothers subunit with the branch on
                                // of the branch so (branchF((i-1)*2), branchF(((i-1)*2)+1))
                                // and (branchF((i)*2), branchF(((i)*2)+1))


                                std::array<double,2> parentVec = m_subunit_unitVecs[1];
                                Eigen::Vector2d parentVecEig;
                                parentVecEig (0) = parentVec[0];
                                parentVecEig (1) = parentVec[1];
                                std::array<double,2> branchVec = actinVec[j].getUnitVec(1);
                                Eigen::Vector2d branchVecEig;
                                branchVecEig (0) = branchVec[0];
                                branchVecEig (1) = branchVec[1];

                                double dot = branchVecEig.dot(parentVecEig); // dot product between [x1, y1] and [x2, y2]
                                double det = branchVecEig(1)*parentVecEig(0) - branchVecEig(0)*parentVecEig(1); // determinant
                                double trueAngle = atan2(det, dot);

                                double restAngle = actinVec[j].getBranchDir()*s_branchAngle;

                                double dAngle = trueAngle - restAngle;
                                dAngle = dAngle - 2*M_PI*floor(dAngle/(2*M_PI)+0.5);

                                // Now we need to get the direction, which is the normal vector of the branch
                                // pointing away or towards the mother?
                                double xcomp = (m_points[i+1][1] - m_points[i][1]);
                                double ycomp = (m_points[i][0] - m_points[i+1][0]);

                                // Normalise this
                                double mag = sqrt(xcomp*xcomp + ycomp*ycomp);
                                xcomp /= mag;
                                ycomp /= mag;


                                double forceX1 = s_torsionCoeff*dAngle*(xcomp/(m_actualSubLengths[i+1]));
                                double forceY1 = s_torsionCoeff*dAngle*(ycomp/(m_actualSubLengths[i+1]));

                                branchF(0) = branchF(0) + forceX1;
                                branchF(1) = branchF(1) + forceY1;

                                branchF(2) = branchF(2) - forceX1;
                                branchF(3) = branchF(3) - forceY1;


                        }
                        else if (i == m_num_subunits-1)
                        {
                                //--------------------------------------------------------------

                                // Force is put on the end mid-point instead

                                branchF((i-1)*2) = branchF((i-1)*2) - s_branchStiff*(branchx - actinVec[j].getPointedEnd()[0]);
                                branchF(((i-1)*2)+1) = branchF(((i-1)*2)+1) - s_branchStiff*(branchy - actinVec[j].getPointedEnd()[1]);

                                //--------------------------------------------------------------
                                std::array<double,2> parentVec = m_subunit_unitVecs[subid];
                                Eigen::Vector2d parentVecEig;
                                parentVecEig (0) = parentVec[0];
                                parentVecEig (1) = parentVec[1];
                                std::array<double,2> branchVec = actinVec[j].getUnitVec(1);
                                Eigen::Vector2d branchVecEig;
                                branchVecEig (0) = branchVec[0];
                                branchVecEig (1) = branchVec[1];

                                double dot = branchVecEig.dot(parentVecEig); // dot product between [x1, y1] and [x2, y2]
                                double det = branchVecEig(1)*parentVecEig(0) - branchVecEig(0)*parentVecEig(1); // determinant
                                double trueAngle = atan2(det, dot);




                                double restAngle = actinVec[j].getBranchDir()*s_branchAngle;

                                double dAngle = trueAngle - restAngle;
                                dAngle = dAngle - 2*M_PI*floor(dAngle/(2*M_PI)+0.5);

                                double xcomp = (m_points[i+1][1] - m_points[i][1]);
                                double ycomp = (m_points[i][0] - m_points[i+1][0]);

                                // Normalise this
                                double mag = sqrt(xcomp*xcomp + ycomp*ycomp);
                                xcomp /= mag;
                                ycomp /= mag;


                                double forceX1 = s_torsionCoeff*dAngle*(xcomp/(m_actualSubLengths[i-1]));
                                double forceY1 = s_torsionCoeff*dAngle*(ycomp/(m_actualSubLengths[i-1]));


                                branchF(((i-1)*2)-2) = branchF(((i-1)*2)-2) + forceX1;
                                branchF(((i-1)*2)-1) = branchF(((i-1)*2)-1) + forceY1;


                                branchF(((i-1)*2)) = branchF(((i-1)*2)) - forceX1;
                                branchF(((i-1)*2)+1) = branchF(((i-1)*2)+1) - forceY1;
                        }
                        else
                        {
                                //--------------------------------------------------------------
                                double a = lenAlongSub / m_prescribedSubLengths[i];
                                assert (a < 1);
                                // pointed end forces
                                branchF((i-1)*2) = branchF((i-1)*2) + (1-a)*(s_branchStiff*(actinVec[j].getPointedEnd()[0] - branchx));
                                branchF(((i-1)*2)+1) = branchF(((i-1)*2)+1) + (1-a)*(s_branchStiff*(actinVec[j].getPointedEnd()[1] - branchy));

                                // barbed end forces
                                branchF((i)*2) = branchF((i)*2) + a*(s_branchStiff*(actinVec[j].getPointedEnd()[0] - branchx));
                                branchF(((i)*2)+1) = branchF(((i)*2)+1) + a*(s_branchStiff*(actinVec[j].getPointedEnd()[1]- branchy));


                                //--------------------------------------------------------------
                                std::array<double,2> parentVec = m_subunit_unitVecs[subid];
                                Eigen::Vector2d parentVecEig;
                                parentVecEig (0) = parentVec[0];
                                parentVecEig (1) = parentVec[1];
                                std::array<double,2> branchVec = actinVec[j].getUnitVec(1);
                                Eigen::Vector2d branchVecEig;
                                branchVecEig (0) = branchVec[0];
                                branchVecEig (1) = branchVec[1];

                                double dot = branchVecEig.dot(parentVecEig); // dot product between [x1, y1] and [x2, y2]
                                double det = branchVecEig(1)*parentVecEig(0) - branchVecEig(0)*parentVecEig(1); // determinant
                                double trueAngle = atan2(det, dot);

                                double restAngle = actinVec[j].getBranchDir()*s_branchAngle;

                                double dAngle = trueAngle - restAngle;
                                dAngle = dAngle - 2*M_PI*floor(dAngle/(2*M_PI)+0.5);

                                double xcomp = (m_points[i+1][1] - m_points[i][1]);
                                double ycomp = (m_points[i][0] - m_points[i+1][0]);

                                // Normalise this
                                double mag = sqrt(xcomp*xcomp + ycomp*ycomp);
                                xcomp /= mag;
                                ycomp /= mag;


                                double forceX1 = s_torsionCoeff*dAngle*(xcomp/(m_actualSubLengths[i]));
                                double forceY1 = s_torsionCoeff*dAngle*(ycomp/(m_actualSubLengths[i]));
                                // force on the second point is apparently equal and opposite

                                branchF((i-1)*2) = branchF((i-1)*2) + forceX1;
                                branchF(((i-1)*2)+1) = branchF(((i-1)*2)+1) + forceY1;
                                branchF((i)*2) = branchF((i)*2) - forceX1;
                                branchF(((i)*2)+1) = branchF(((i)*2)+1) - forceY1;

                        }
                }
        }

        return branchF;
}

Eigen::VectorXd Actin::calcMotherForces(std::vector<Actin> &actinVec, int N)
{

        // function that gets the hookean forces from the mother

        // called on the branch
        Eigen::VectorXd motherF (2*N);
        motherF.setZero();
        // if this doesnt have a mother then return 0 force
        if (m_parent_id == -1)
        {
                return motherF;
        }

        double lenAlongSub;
        int subid = actinVec[m_parent_id].findSubunit(m_motherMonomer, lenAlongSub);
        double specialBranchSubPos = lenAlongSub / actinVec[m_parent_id].getPresSubLengths()[subid];
        assert(specialBranchSubPos < 1);
        double branchx = (actinVec[m_parent_id].getPoint(subid)[0]
                          + (specialBranchSubPos
                             * (actinVec[m_parent_id].getPoint(subid+1)[0]
                                - actinVec[m_parent_id].getPoint(subid)[0])));

        double branchy = (actinVec[m_parent_id].getPoint(subid)[1]
                          + (specialBranchSubPos
                             * (actinVec[m_parent_id].getPoint(subid+1)[1]
                                - actinVec[m_parent_id].getPoint(subid)[1])));

        //--------------------------------------------------------------------------
        // Linear component

        motherF(0) = motherF(0) + (-s_branchStiff*(m_points[0][0] - branchx));
        motherF(1) = motherF(1) + (-s_branchStiff*(m_points[0][1] - branchy));


        //--------------------------------------------------------------------------
        // Torsional component

        std::array<double,2> parentVec = actinVec[m_parent_id].getUnitVec(subid);
        Eigen::Vector2d parentVecEig;
        parentVecEig (0) = parentVec[0];
        parentVecEig (1) = parentVec[1];
        std::array<double,2> branchVec = m_subunit_unitVecs[1];
        Eigen::Vector2d branchVecEig;
        branchVecEig (0) = branchVec[0];
        branchVecEig (1) = branchVec[1];

        double dot = branchVecEig.dot(parentVecEig);  // dot product between [x1, y1] and [x2, y2]
        double det = branchVecEig(1)*parentVecEig(0) - branchVecEig(0)*parentVecEig(1); // determinant
        double trueAngle = atan2(det, dot);

        double restAngle = m_branchdir*s_branchAngle;

        double dAngle = trueAngle - restAngle;
        dAngle = dAngle - 2*M_PI*floor(dAngle/(2*M_PI)+0.5);

        double xcomp = (m_points[2][1] - m_points[1][1]);
        double ycomp = (m_points[1][0] - m_points[2][0]);
        // Normalise this
        double mag = sqrt(xcomp*xcomp + ycomp*ycomp);
        xcomp /= mag;
        ycomp /= mag;


        double forceX1 = s_torsionCoeff*dAngle*(xcomp/(m_actualSubLengths[1]));
        double forceY1 = s_torsionCoeff*dAngle*(ycomp/(m_actualSubLengths[1]));

        motherF(0) = motherF(0) - forceX1;
        motherF(1) = motherF(1) - forceY1;
        motherF(2) = motherF(2) + forceX1;
        motherF(3) = motherF(3) + forceY1;

        return motherF;

}

void Actin::sever(std::vector<Actin> &actinVec, double currTime, int severMonomer,
                  int &nActin, const bool steric, StericGrid &stericGrid,
                  GactinGrid &gActinGrid,
                  const bool arpPool, ArpGrid &arpGrid)
{
        // Given the filament and after having passed the probability check
        // Need to sever this filament

        // Create a new filament with start point (pointed) equal to the sever point
        // and end point equal to the end point of the original


        double lenAlongSub = 0;
        int severSub = findSubunit(severMonomer, lenAlongSub);
        std::array<double,2> severPoint = findMonoCoord(severMonomer);

        // Decrease this here, and then add to it as we make new filament
        s_total_subunits -= m_num_subunits;

        if (steric)
        {
            for (int i = 0; i < m_num_subunits; ++i)
            {
                    stericGrid.resetCells(*this, i);
            }
        }




        // Create new filament - copy everything from the severPoint

        if (severMonomer > m_num_monomers-s_seedSize-1)
        {
            // new filament is shorter than a trimer so don't bother making it,
            // it will dissociate. Dissociation point is the severPoint
            if (gActinGrid.getExist())
                    gActinGrid.dissociate(severPoint[0], severPoint[1]);
            //WHAT ABOUT ANY BRANCHES ON THIS
            if (m_distToLeadBR <= 1)
            {
                    // Need to debranch this branch that is on the end
                    int leadingBR = 0;
                    int id = 0;
                    for (unsigned int i = 0; i < m_daughterIDs[m_num_subunits-1].size(); ++i)
                    {
                            int k = m_daughterIDs[m_num_subunits-1][i]; // id of daughter
                            if (actinVec[k].getMotherMonoID() >= leadingBR)
                            {
                                    leadingBR = actinVec[k].getMotherMonoID();
                                    id = k;
                            }
                    }
                    deBranch2(actinVec[id], actinVec, nActin, arpPool, arpGrid, steric,
                              stericGrid);

            }
            assert(m_distToLeadBR > 1);

            // Same for crosslinks
            if (m_distToLeadCL <= 1)
            {
                int initNumCLinks = getNumCLinks();
                int numClRemoved = 0;
                for (int i = 0; i < initNumCLinks; ++i)
                {
                    int oldMonoID = m_cLinkActinAndSites[i-numClRemoved][0];

                    if (oldMonoID > severMonomer)
                    {
                        // Crosslink will be removed
                        unLinkBase(m_id, i-numClRemoved, actinVec);
                        numClRemoved++;
                    }
                }
            }
        }
        else
        {
            // Create new filament

            Actin newActin(nActin, severPoint, severMonomer, severSub, lenAlongSub, currTime, *this);

            // Need to move branches aswell

            // First need to check if branches belong to "new filament"
            // This is only done on the branches that belong to the subunit that has the
            // severing point

            int numDaughtersToCheck = m_daughterIDs[severSub].size();
            int numMoved = 0;
            for (int i = 0; i < numDaughtersToCheck; ++i)
            {
                int ID = m_daughterIDs[severSub][i-numMoved];

                if (actinVec[ID].getMotherMonoID() > severMonomer)
                {
                    // Branch needs to be moved to "new" filament

                    // NEed to adjust old filament's branch compat vectors here if necessary

                    actinVec[ID].setParentID(newActin.getID());

                    --m_daughter_num; // reduce num of daughters on "old" fila
                    newActin.incrementDaughternum(); // increase number of daughters on "new" fila
                    // make changes to m_daughterIDs

                    newActin.addBranchtoParentVector(0, ID); // branch is added to 0 index subunit
                    actinVec[ID].setBranchSubunit(0); // this branch is added to the 0th sub
                    actinVec[ID].setMotherMonoID(actinVec[ID].getMotherMonoID()-severMonomer-1);

                    // remove from "old" filament
                    for (int j = i-numMoved; j < numDaughtersToCheck-1-numMoved; ++j)
                    {
                            m_daughterIDs[severSub][j] = m_daughterIDs[severSub][j+1];
                    }
                    m_daughterIDs[severSub].pop_back();
                    ++numMoved;
                }
            }

            // Then we need to change all the branches that come after the severSubUnit
            // as these definitely move!

            for (int i = severSub+1; i < m_num_subunits; ++i)
            {
                while (m_daughterIDs[i].size() > 0)
                {
                    int ID = m_daughterIDs[i][0];

                    // branch needs to be moved to "new" filament
                    actinVec[ID].setParentID(newActin.getID());
                    // What about subsequent branches? They need their structureID changing too
                    --m_daughter_num; // reduce num of daughters on "old" fila
                    newActin.incrementDaughternum(); // increase number of daughters on "new" fila
                    // make changes to m_daughterIDs
                    newActin.addBranchtoParentVector(i-severSub, ID);
                    actinVec[ID].setBranchSubunit(i-severSub); // this branch is added to the 0th sub
                    actinVec[ID].setMotherMonoID(actinVec[ID].getMotherMonoID()-severMonomer-1);
                    // remove from "old" filament

                    for (unsigned int k = 0; k < m_daughterIDs[i].size()-1; ++k)
                    {
                            m_daughterIDs[i][k] = m_daughterIDs[i][k+1];
                    }
                    m_daughterIDs[i].pop_back();

                }
            }

            newActin.recalcDistToBR2(actinVec);
            // change the changed status of both filaments
            newActin.setChangedBoolTrue();

            // NOW DO THE SAME FOR CROSSLINKS!
            int initNumCLinks = getNumCLinks();
            int numClMoved = 0;
            for (int i = 0; i < initNumCLinks; ++i)
            {
                //int oldSubID = m_cLinkActinAndSites[i-numClMoved][1];
                int oldMonoID = m_cLinkActinAndSites[i-numClMoved][0];

                if (oldMonoID > severMonomer)
                {
                    // Crosslink will certainly move!
                    // add crosslink to new actin
                    int monoID = m_cLinkActinAndSites[i-numClMoved][0] - severMonomer - 1;
                    int subID = m_cLinkActinAndSites[i-numClMoved][1] - severSub;
                    int otherFila = m_cLinkActinAndSites[i-numClMoved][2];
                    int otherFilaCLIDX = m_cLinkActinAndSites[i-numClMoved][3];
                    int otherMonoID = m_cLinkActinAndSites[i-numClMoved][4];
                    int otherFilaSubID = m_cLinkActinAndSites[i-numClMoved][5];
                    newActin.addToLinkVect(monoID, subID, otherFila,
                                           otherFilaCLIDX, otherMonoID,
                                           otherFilaSubID);
                    newActin.remMonoFromAvail(monoID);
                    newActin.addEmptyToCLPoints();
                    newActin.recalcDistToLeadCL(monoID);


                    // Have to adjust the linked filament too
                    actinVec[otherFila].setClOther(otherFilaCLIDX, newActin.getID());
                    actinVec[otherFila].setClOtherIDX(otherFilaCLIDX, newActin.getNumCLinks()-1);
                    actinVec[otherFila].setClMonoOther(otherFilaCLIDX, monoID);
                    actinVec[otherFila].setClSubOther(otherFilaCLIDX, subID);

                    // Now remove the side of the crosslink from this filament
                    eraseCLink(i-numClMoved);
                    numClMoved++;

                }

            }

            newActin.changeStructure(actinVec, nActin);
            newActin.changeCLStructure(actinVec); // changes all daughters to have same cl id
            newActin.changeSubStructure(actinVec);



            // reset master bools
            setMasterBoolFalse(actinVec);
            newActin.setMasterBoolFalse(actinVec);
            findNewMaster(actinVec);

            // need the if because could be on same structure still
            if (!newActin.checkForMasterInStructure(actinVec))
            {
                // no master in new
                newActin.findNewMaster(actinVec);
            }

            // Rejig pointed end points if necessary

            if (newActin.getLength() >= 3*s_segmentationLength)
            {
                    newActin.removePointsPointed(actinVec, steric, stericGrid, true);
                    double first2SubLens = newActin.getActualSubLengths()[0] + newActin.getActualSubLengths()[1]; // this needs to be kept the same
                    newActin.setActualSubLength(1, (first2SubLens + s_segmentationLength/2)/2);
                    newActin.setActualSubLength(0, first2SubLens - newActin.getActualSubLengths()[1]);

                    first2SubLens = newActin.getPresSubLengths()[0] + newActin.getPresSubLengths()[1]; // this needs to be kept the same
                    newActin.setPresSubLength(1, (first2SubLens + s_segmentationLength/2)/2);
                    newActin.setPresSubLength(0, first2SubLens - newActin.getPresSubLengths()[1]);

                    newActin.updatePointsDown(newActin.getNumSubs()-1);
            }


            newActin.straightenEnd(false);
            for (int i = 0; i < newActin.getNumSubs(); ++i)
            {
                    newActin.checkAndMoveBranch(actinVec, i);
            }
            newActin.rotbranches(actinVec, true);

            newActin.checkAndMoveCl(actinVec);


            newActin.updateCrossLinkLocs(actinVec);

            assert(newActin.getNumSubs() >= 3);

            actinVec.push_back(newActin);
            ++nActin;
            s_total_subunits += newActin.getNumSubs();

            if (steric)
            {
                    for (int i = 0; i < newActin.getNumSubs(); ++i)
                    {
                            stericGrid.checkSubandUpdate(actinVec[nActin-1], i);
                    }
            }
        }


        // ------------------ Dealing with the "old" filament ----------------------
        // Remove everything after severPoint

        if ( (severMonomer < s_seedSize-1 && m_parent_id == -1) || (severMonomer < s_branchSeedSize-1 && m_parent_id != -1) )
        {
            // The old filament is dissociated
            // There could be branches here, still attached!
            if (m_distToLeadBR >= m_num_monomers-2 && m_distToLeadBR != m_num_monomers) // branch is on pointed end
            {
                    assert(m_daughter_num == 1);
                    // Need to debranch this branch that is on the end
                    int leadingBR = m_num_monomers;
                    int id = 0;
                    for (unsigned int i = 0; i < m_daughterIDs[0].size(); ++i)
                    {
                            int k = m_daughterIDs[0][i]; // id of daughter
                            if (actinVec[k].getMotherMonoID() <= leadingBR)
                            {
                                    leadingBR = actinVec[k].getMotherMonoID();
                                    id = k;
                            }
                    }
                    deBranch2(actinVec[id], actinVec, nActin, arpPool, arpGrid, steric,
                              stericGrid);

            }
            assert(m_daughter_num == 0);

            int initNumCLinks = getNumCLinks();
            int numClRemoved = 0;
            for (int i = 0; i < initNumCLinks; ++i)
            {
                int oldMonoID = m_cLinkActinAndSites[i-numClRemoved][0];

                if (oldMonoID <= severMonomer)
                {
                    // Crosslink will be removed
                    unLinkBase(m_id, i-numClRemoved, actinVec);
                    numClRemoved++;

                }

            }


            std::vector<int> dissociateIDs;
            dissociateIDs.push_back(m_id);
            dissociationRout(dissociateIDs, actinVec, nActin, steric, stericGrid);
            // Following line is needed as above function decreases this by 3,
            // We have already taken this into account!
            s_total_subunits += 3;

            if (gActinGrid.getExist())
                    gActinGrid.dissociate(severPoint[0], severPoint[1]);

        }
        else
        {
            m_num_monomers = severMonomer + 1;
            m_length = m_num_monomers*s_monomerLength;

            if (m_length < 3*s_segmentationLength)
            {
                // still needs 4 points and 3 subs
                for (int i = 3; i < m_num_subunits; ++i)
                {
                    m_actualSubLengths.pop_back();
                    m_prescribedSubLengths.pop_back();
                    m_daughterIDs.pop_back();
                    m_subunit_unitVecs.pop_back();

                    m_points.pop_back();

                }

                m_num_subunits = 3;
                m_points[m_num_subunits][0] = severPoint[0];
                m_points[m_num_subunits][1] = severPoint[1];


                m_actualSubLengths[0] = m_length/4;
                m_actualSubLengths[1] = m_length/2;
                m_actualSubLengths[2] = m_length/4;

                m_prescribedSubLengths[0] = m_length/4;
                m_prescribedSubLengths[1] = m_length/2;
                m_prescribedSubLengths[2] = m_length/4;

                straightenFila(); // this assumes the 0th point is correct

            }
            else
            {

                int subToStart = severSub + 1;
                if (subToStart < 3)
                        subToStart = 3;

                for (int i = subToStart; i < m_num_subunits; ++i)
                {
                    m_actualSubLengths.pop_back();
                    m_prescribedSubLengths.pop_back();
                    m_daughterIDs.pop_back();
                    m_subunit_unitVecs.pop_back();

                    m_points.pop_back();

                }

                m_num_subunits = subToStart;
                m_points[m_num_subunits][0] = severPoint[0];
                m_points[m_num_subunits][1] = severPoint[1];
                m_prescribedSubLengths[m_num_subunits-1] = lenAlongSub + (0.5*s_monomerLength);
                updateSubLengths();
                updateUnitVecsNoCheck();

                // ----------

                removePointsBarbed(actinVec, steric, stericGrid, true);

                double a = m_actualSubLengths[m_num_subunits-2] + m_actualSubLengths[m_num_subunits-1]; // this needs to be kept the same
                m_actualSubLengths[m_num_subunits-2] = (a + s_segmentationLength/2)/2;
                m_actualSubLengths[m_num_subunits-1] = a - m_actualSubLengths[m_num_subunits-2];

                a = m_prescribedSubLengths[m_num_subunits-2] + m_prescribedSubLengths[m_num_subunits-1]; // this needs to be kept the same
                m_prescribedSubLengths[m_num_subunits-2] = (a + s_segmentationLength/2)/2;
                m_prescribedSubLengths[m_num_subunits-1] = a - m_prescribedSubLengths[m_num_subunits-2];

                updatePointsUp(0);
            }

            straightenEnd(true);
            for (int i = 0; i < m_num_subunits; ++i)
            {
                    checkAndMoveBranch(actinVec, i);
            }
            rotbranches(actinVec, true);

            // We may have lost some crosslinks to the new filament, need to
            // check everything points to the correct crossLink

            for (int i = 0; i < getNumCLinks(); ++i)
            {
              // Make sure any linked actin still point to the right cl idx

              int otherFila = m_cLinkActinAndSites[i][2];
              int otherFilaCLIDX = m_cLinkActinAndSites[i][3];
              actinVec[otherFila].setClOtherIDX(otherFilaCLIDX, i);
            }

            checkAndMoveCl(actinVec);

            // If a branch was on the barbed-end side of the sever point, but within
            // s_branchSpacing of the sever point then need to update the compatibility
            // vectors for the "original" filament, as those end monomers are now
            // available to branch


            int numRemoved = 0;
            int numToCheck = m_availMonos.size();
            for (int i = 0; i < numToCheck; ++i)
            {
                if (m_availMonos[i-numRemoved] > severMonomer)
                {
                        for (int j = i-numRemoved; j < numToCheck-1-numRemoved; ++j)
                        {
                                m_availMonos[j] = m_availMonos[j+1];
                        }
                        m_availMonos.pop_back();
                        ++numRemoved;
                }
            }



            int monosCompatSize = m_monosCompat.size();
            while (monosCompatSize-1 > severMonomer)
            {
                m_monosCompat.pop_back();
                m_monoBirthTime.pop_back();
                monosCompatSize -= 1;
            }

            recalcDistToBR2(actinVec);
            m_changed = true;

            m_distToLeadCL = m_num_monomers;
            for (int i = 0; i < getNumCLinks(); ++i)
            {
                int monomerID = m_cLinkActinAndSites[i][0];
                recalcDistToLeadCL(monomerID);
            }

            recalcMonoCompatALL(actinVec);


            s_total_subunits += m_num_subunits;

            // If the filament was a long branch, it will have had a different
            // substructure id to its parent.
            // We need to check to see if it is no longer a long branch, so we
            // can reassign the parent's substructure
            if (m_parent_id != -1)
            {
                if (m_subStructure_id != actinVec[m_parent_id].getSubStructureID() && m_length < 2*Actin::s_segmentationLength)
                {
                    m_subStructure_id = actinVec[m_parent_id].getSubStructureID();
                    // Need to repeat for any branches off this fila
                    changeSubStructure(actinVec);
                }
            }


            if (steric)
            {
                    for (int i = 0; i < m_num_subunits; ++i)
                    {
                            stericGrid.checkSubandUpdate(*this, i);
                    }
            }
            assert(m_num_subunits >= 3);
        }
}

void Actin::eraseCLink(int idx)
{
    // Only removes from this filament not the other!
    m_cLinkActinAndSites.erase(m_cLinkActinAndSites.begin()+idx);
    m_cLinkPoints.erase(m_cLinkPoints.begin()+idx);
}

double Actin::getSumAngle(int subUnit, int centreSub) const
{
        // Function that retrieves the sumLocalAngles of the particular subunit

        double sumLocalAngles = 0;
        if (centreSub > subUnit)
        {
                // Need to move down the filament
                double subAngle = atan2(-m_subunit_unitVecs[subUnit][1], -m_subunit_unitVecs[subUnit][0]);
                double centAngle = atan2(-m_subunit_unitVecs[centreSub][1], -m_subunit_unitVecs[centreSub][0]);
                sumLocalAngles = subAngle - centAngle;



        }
        else if (centreSub < subUnit)
        {
                // move up the filament
                double subAngle = atan2(m_subunit_unitVecs[subUnit][1], m_subunit_unitVecs[subUnit][0]);
                double centAngle = atan2(m_subunit_unitVecs[centreSub][1], m_subunit_unitVecs[centreSub][0]);
                sumLocalAngles = subAngle - centAngle;

        }
        else
        {
                // the sever subunit IS the centre
                sumLocalAngles = 0;
        }

        return sumLocalAngles;

}

void Actin::changeToNewAngle(int index, int centreSub, double newAngle)
{
        // As before but we need to do this using our unit vectors now, so it's
        // a little more complicated

        // global angle of previous sub
        double globAngle;
        // Update along the filament
        if (index >= centreSub)
        {
                globAngle = getSumAngle(index-1, centreSub) + atan2(m_subunit_unitVecs[centreSub][1], m_subunit_unitVecs[centreSub][0]);
                m_subunit_unitVecs[index][0] = cos(globAngle+newAngle);
                m_subunit_unitVecs[index][1] = sin(globAngle+newAngle);
                if (index+1 < m_num_subunits)
                {
                        changeToNewAngle(index+1, centreSub, 0);

                }

        }
        else
        {
                globAngle = getSumAngle(index+1, centreSub) + atan2(-m_subunit_unitVecs[centreSub][1], -m_subunit_unitVecs[centreSub][0]);
                m_subunit_unitVecs[index][0] = -cos(globAngle-newAngle);
                m_subunit_unitVecs[index][1] = -sin(globAngle-newAngle);
                if (index-1 >= 0)
                {
                        changeToNewAngle(index-1, centreSub, 0);

                }
        }
}

void Actin::initialiseLength(int numMonomers, std::vector<Actin> &actinVec,
                             const bool steric, StericGrid &stericGrid,
                             const std::vector<Membrane> &membranes,
                             Cortex &cortex,
                             const std::vector<ExcZone> &excZones,
                             const std::vector<MembraneWall> &memWalls,
                             const bool tethering)
{
        // Called when want to change the length manually for some reason in testing
        m_num_monomers = numMonomers;
        m_length = numMonomers * s_monomerLength;
        m_distToLeadBR = m_num_monomers;
        m_distToLeadCL = m_num_monomers;

        for (int i = s_seedSize; i < m_num_monomers; ++i)
        {
            if (i < s_maxSpacing)
            {
                m_monosCompat.push_back(false);
            }
            else
            {
                m_monosCompat.push_back(true);
                m_availMonos.push_back(i);
            }

            m_monoBirthTime.push_back(0);
        }

        if (m_length > 3*s_segmentationLength)
        {
            m_flexible = true;
        }

        // Set old endpoint to be end of filament
        m_points[3][0] = m_points[1][0] + ((3*m_length/4)*m_subunit_unitVecs[0][0]); // x
        m_points[3][1] = m_points[1][1] + ((3*m_length/4)*m_subunit_unitVecs[0][1]); // y
        m_points[0][0] = m_points[1][0] - ((m_length/4)*m_subunit_unitVecs[0][0]); // x
        m_points[0][1] = m_points[1][1] - ((m_length/4)*m_subunit_unitVecs[0][1]); // y

        if (m_length >= 3*m_des_subunit_len)
        {
                m_points[1][0] = m_points[0][0] + ((m_des_subunit_len/2)*m_subunit_unitVecs[0][0]); // x
                m_points[1][1] = m_points[0][1] + ((m_des_subunit_len/2)*m_subunit_unitVecs[0][1]); // y
                m_actualSubLengths[0] = m_des_subunit_len/2;
                m_prescribedSubLengths[0] = m_des_subunit_len/2;

                m_points[2][0] = m_points[1][0] + (m_des_subunit_len*m_subunit_unitVecs[0][0]); // x
                m_points[2][1] = m_points[1][1] + (m_des_subunit_len*m_subunit_unitVecs[0][1]); // y
                m_actualSubLengths[1] = m_des_subunit_len;
                m_prescribedSubLengths[1] = m_des_subunit_len;

                std::array<double,3> tmp { };
                tmp[0] = m_points[3][0];
                tmp[1] = m_points[3][1];

                m_points[3][0] = m_points[2][0] + (m_des_subunit_len*m_subunit_unitVecs[0][0]); // x
                m_points[3][1] = m_points[2][1] + (m_des_subunit_len*m_subunit_unitVecs[0][1]); // y
                m_actualSubLengths[2] = m_des_subunit_len;
                m_prescribedSubLengths[2] = m_des_subunit_len;

                m_points.push_back(tmp);
                m_actualSubLengths.push_back(m_length - (5./2)*m_des_subunit_len);
                m_prescribedSubLengths.push_back(m_length - (5./2)*m_des_subunit_len);
                ++s_total_subunits;
                ++m_num_subunits;
                m_subunit_unitVecs.push_back({0,0});
                updateUnitVecsNoCheck();

                m_daughterIDs.push_back(std::vector<int>());

                if (steric)
                        m_stericCells.push_back(std::vector<int>());

        }
        else
        {
                // Still only 4 points
                m_points[1][0] = m_points[0][0] + ((m_length/4)*m_subunit_unitVecs[0][0]); // x
                m_points[1][1] = m_points[0][1] + ((m_length/4)*m_subunit_unitVecs[0][1]); // y

                m_points[2][0] = m_points[1][0] + ((m_length/2)*m_subunit_unitVecs[0][0]); // x
                m_points[2][1] = m_points[1][1] + ((m_length/2)*m_subunit_unitVecs[0][1]); // y

                m_points[3][0] = m_points[2][0] + ((m_length/4)*m_subunit_unitVecs[0][0]); // x
                m_points[3][1] = m_points[2][1] + ((m_length/4)*m_subunit_unitVecs[0][1]); // y

                m_actualSubLengths[0] = m_length/4;
                m_actualSubLengths[1] = m_length/2;
                m_actualSubLengths[2] = m_length/4;

                m_prescribedSubLengths[0] = m_length/4;
                m_prescribedSubLengths[1] = m_length/2;
                m_prescribedSubLengths[2] = m_length/4;

                // Have to reset and update steric grid cells
                if (steric)
                {
                        for (int i = 0; i < m_num_subunits; ++i)
                        {
                                stericGrid.resetSubandUpdate(*this, i);
                        }
                }
                if (steric)
                {
                    if (check_bend(actinVec, excZones, memWalls, membranes, cortex, steric, stericGrid))
                    {
                        // initialisation failed, end simulation
                        std::cout << "Init length failed due to steric hindrance, either have shorter filaments or a larger nucleation region.\n";
                        std::cout << "Or considering adding code to reposition filaments that are hindered, although this might be inefficient!\n";
                        std::cout << "(Relevant code is in Actin.cpp, lines 5753ish)\n";
                        std::cout << "Ending simulation" << std::endl;
                        exit(1);
                    }
                }

                return;
        }

        while (m_actualSubLengths[m_num_subunits-1] > (m_des_subunit_len/2))
        {
                // Do something, either...
                // Add a new point,
                // or move the penultimate point
                if ((m_actualSubLengths[m_num_subunits-1]+m_actualSubLengths[m_num_subunits-2]) < (2.5*m_des_subunit_len))
                {
                        // move penultimate point
                        // dont need to add another point
                        std::array<double,2> endpoint { m_points[m_num_subunits][0], m_points[m_num_subunits][1] };
                        std::array<double,2> ppenultimate { m_points[m_num_subunits-2][0], m_points[m_num_subunits-2][1] };

                        double lengthToEnd = sqrt(distanceBetPoints2DSQR(endpoint, ppenultimate));
                        double endsubL;
                        if (m_num_subunits==4)
                        {

                                lengthToEnd -= m_actualSubLengths[0];
                                lengthToEnd /= 2;
                                endsubL = lengthToEnd;
                                lengthToEnd += m_actualSubLengths[0];
                        }
                        else
                        {
                                lengthToEnd -= m_actualSubLengths[m_num_subunits-4]/2;
                                lengthToEnd /= 2;
                                endsubL = lengthToEnd;
                                lengthToEnd += m_actualSubLengths[m_num_subunits-4]/2;

                        }
                        // now lengthToEnd is the subunit length of the penultimate sub

                        m_points[m_num_subunits-1][0] =  m_points[m_num_subunits-2][0] + (lengthToEnd*m_subunit_unitVecs[0][0]);
                        m_points[m_num_subunits-1][1] =  m_points[m_num_subunits-2][1] + (lengthToEnd*m_subunit_unitVecs[0][1]);

                        m_actualSubLengths[m_num_subunits-1] = endsubL;
                        m_actualSubLengths[m_num_subunits-2] = lengthToEnd;
                        // its not had chance to go through the bend algorithm so actual and prescribed will be equal
                        m_prescribedSubLengths[m_num_subunits-1] = endsubL;
                        m_prescribedSubLengths[m_num_subunits-2] = lengthToEnd;
                        break;
                }
                else
                {
                        // add new point
                        ++s_total_subunits;
                        ++m_num_subunits;
                        // Change the new end point to be the old end point
                        // the old end point
                        std::array<double,3> tmp { };
                        tmp[0] = m_points[m_num_subunits-1][0];
                        tmp[1] = m_points[m_num_subunits-1][1];
                        m_points.push_back(tmp);
                        m_actualSubLengths[m_num_subunits-2] = m_des_subunit_len;
                        m_prescribedSubLengths[m_num_subunits-2] = m_des_subunit_len;



                        // Determine new point
                        m_points[m_num_subunits-1][0] = m_points[m_num_subunits-2][0] + (m_des_subunit_len*m_subunit_unitVecs[0][0]); // x
                        m_points[m_num_subunits-1][1] = m_points[m_num_subunits-2][1] + (m_des_subunit_len*m_subunit_unitVecs[0][1]); // y

                        std::array<double,2> point1 { m_points[m_num_subunits][0], m_points[m_num_subunits][1] }; // end point
                        std::array<double,2> point2 { m_points[m_num_subunits-1][0], m_points[m_num_subunits-1][1] }; // new point

                        m_actualSubLengths.push_back(sqrt(distanceBetPoints2DSQR(point1, point2)));
                        m_prescribedSubLengths.push_back(sqrt(distanceBetPoints2DSQR(point1, point2)));
                        m_subunit_unitVecs.push_back({0,0});
                        updateUnitVecs();

                        m_daughterIDs.push_back(std::vector<int>());

                        if (steric)
                                m_stericCells.push_back(std::vector<int>());


                }
        }

        // Have to reset and update steric grid cells
        if (steric)
        {
                for (int i = 0; i < m_num_subunits; ++i)
                {
                        stericGrid.resetSubandUpdate(*this, i);
                }
        }
        if (steric)
        {
            if (check_bend(actinVec, excZones, memWalls, membranes, cortex, steric, stericGrid))
            {
                // initialisation failed, end simulation
                std::cout << "Init length failed due to steric hindrance, either have shorter filaments or a larger nucleation region.\n";
                std::cout << "Or considering adding code to reposition filaments that are hindered, although this might be inefficient!\n";
                std::cout << "(Relevant code is in Actin.cpp, lines 5753ish)\n";
                std::cout << "Ending simulation" << std::endl;
                exit(1);
            }
        }

        if (tethering)
        {
          // Move tetherPoint to pointed end
          m_tetherMono[0] = 0;
          m_tetherSubunits[0] = 0;

          std::array<double,4> newTether;

          // Tether point half a monomer length from end
          newTether[0] = m_points[0][0] + m_subunit_unitVecs[0][0]*(s_monomerLength/2);
          newTether[1] = m_points[0][1] + m_subunit_unitVecs[0][1]*(s_monomerLength/2);
          newTether[2] = newTether[0];
          newTether[3] = newTether[1];

          m_tetherPoints[0] = newTether;


          m_distToTetBarb = m_num_monomers-1;
          m_distToTetPoint = 0;
          m_indexLeadingBarbTet = 0;
          m_indexLeadingPointTet = 0;
        }

}



void Actin::updateUnitVecs()
{
        // Function that updates the unit vectors of the filament
        for (int i = 0; i < m_num_subunits; ++i)
        {
                // x term
                double dx = m_points[i+1][0] - m_points[i][0];
                double dy = m_points[i+1][1] - m_points[i][1];
                double mag = sqrt(dx*dx + dy*dy);

                assert(mag > 0);
                m_subunit_unitVecs[i][0] = dx/mag;
                m_subunit_unitVecs[i][1] = dy/mag;
                assert(m_subunit_unitVecs[i][0] != 0 || m_subunit_unitVecs[i][1] != 0);
        }

        assert(fabs(m_subunit_unitVecs[0][0] - m_subunit_unitVecs[1][0]) < 1e-7 && fabs(m_subunit_unitVecs[0][1] - m_subunit_unitVecs[1][1]) < 1e-7);
        m_subunit_unitVecs[0] = m_subunit_unitVecs[1];
        assert(fabs(m_subunit_unitVecs[m_num_subunits-1][0] - m_subunit_unitVecs[m_num_subunits-2][0]) < 1e-7 && fabs(m_subunit_unitVecs[m_num_subunits-1][1] - m_subunit_unitVecs[m_num_subunits-2][1]) < 1e-7);
        m_subunit_unitVecs[m_num_subunits-1] = m_subunit_unitVecs[m_num_subunits-2];

}

void Actin::updateUnitVecsNoCheck()
{
        // Function that updates the unit vectors of the filament
        //std::cout << m_id << std::endl;
        for (int i = 0; i < m_num_subunits; ++i)
        {
                // x term
                double dx = m_points[i+1][0] - m_points[i][0];
                double dy = m_points[i+1][1] - m_points[i][1];
                double mag = sqrt(dx*dx + dy*dy);
                m_subunit_unitVecs[i][0] = dx/mag;
                m_subunit_unitVecs[i][1] = dy/mag;
                assert(m_subunit_unitVecs[i][0] != 0 || m_subunit_unitVecs[i][1] != 0);
        }
}

void Actin::updateUnitVec(int p)
{
        // Function that updates the unit vector of subunit p

        double dx = m_points[p+1][0] - m_points[p][0];
        double dy = m_points[p+1][1] - m_points[p][1];
        double mag = sqrt(dx*dx + dy*dy);
        m_subunit_unitVecs[p][0] = dx/mag;
        m_subunit_unitVecs[p][1] = dy/mag;
        assert(m_subunit_unitVecs[p][0] != 0 || m_subunit_unitVecs[p][1] != 0);

}

Eigen::VectorXd Actin::calcBendingForces(double temp, int N)
{
        double k = (g_Kb*temp*s_persistenceLength) / m_des_subunit_len; // will need to divide this by the relevant segment length

        Eigen::VectorXd bendingForces (2*N);
        if (N == 2)
        {
                //Cant bend
                bendingForces = Eigen::VectorXd::Zero(N*2);
                return bendingForces;
        }

        for (int n = 0; n < N; ++n)
        {
                // n being the bead id
                if (n == 0)
                {
                        bendingForces[0] = -k*(m_subunit_unitVecs[2][0] - m_subunit_unitVecs[1][0])/m_actualSubLengths[1];
                        bendingForces[1] = -k*(m_subunit_unitVecs[2][1] - m_subunit_unitVecs[1][1])/m_actualSubLengths[1];
                }
                else if (n == 1)
                {
                        if (N == 3)
                        {
                                //ignore negative term?
                                bendingForces[2] = k*(m_subunit_unitVecs[2][0] - m_subunit_unitVecs[1][0])*(1/m_actualSubLengths[1] + 1/m_actualSubLengths[2]);

                                bendingForces[3] = k*(m_subunit_unitVecs[2][1] - m_subunit_unitVecs[1][1])*(1/m_actualSubLengths[1] + 1/m_actualSubLengths[2]);
                        }
                        else
                        {
                                bendingForces[2] = k*(m_subunit_unitVecs[2][0] - m_subunit_unitVecs[1][0])*(1/m_actualSubLengths[1] + 1/m_actualSubLengths[2])
                                                   - k*(m_subunit_unitVecs[3][0] - m_subunit_unitVecs[2][0])/m_actualSubLengths[2];

                                bendingForces[3] = k*(m_subunit_unitVecs[2][1] - m_subunit_unitVecs[1][1])*(1/m_actualSubLengths[1] + 1/m_actualSubLengths[2])
                                                   - k*(m_subunit_unitVecs[3][1] - m_subunit_unitVecs[2][1])/m_actualSubLengths[2];
                        }
                }
                else if (n == N-2)
                {
                        // N-1 case
                        bendingForces[(2*N)-4] = -k*(m_subunit_unitVecs[N-2][0] - m_subunit_unitVecs[N-3][0])/m_actualSubLengths[N-2]
                                                 +k*(m_subunit_unitVecs[N-1][0] - m_subunit_unitVecs[N-2][0])*(1/m_actualSubLengths[N-2]+1/m_actualSubLengths[N-1]);

                        bendingForces[(2*N)-3] = -k*(m_subunit_unitVecs[N-2][1] - m_subunit_unitVecs[N-3][1])/m_actualSubLengths[N-2]
                                                 +k*(m_subunit_unitVecs[N-1][1] - m_subunit_unitVecs[N-2][1])*(1/m_actualSubLengths[N-2]+1/m_actualSubLengths[N-1]);

                }
                else if (n == N-1)
                {
                        bendingForces[(2*N)-2] = -k*(m_subunit_unitVecs[N-1][0] - m_subunit_unitVecs[N-2][0])/m_actualSubLengths[N-1];
                        bendingForces[(2*N)-1] = -k*(m_subunit_unitVecs[N-1][1] - m_subunit_unitVecs[N-2][1])/m_actualSubLengths[N-1];
                }
                else
                {
                        // 2 <= n <= N-3 where N is the num of points
                        bendingForces[2*n] = -k*(m_subunit_unitVecs[n][0] - m_subunit_unitVecs[n-1][0])/m_actualSubLengths[n]
                                             + k*(m_subunit_unitVecs[n+1][0] - m_subunit_unitVecs[n][0])*(1/m_actualSubLengths[n] + 1/m_actualSubLengths[n+1])
                                             - k*(m_subunit_unitVecs[n+2][0] - m_subunit_unitVecs[n+1][0])/m_actualSubLengths[n+1];




                        bendingForces[(2*n)+1] = -k*(m_subunit_unitVecs[n][1] - m_subunit_unitVecs[n-1][1])/m_actualSubLengths[n]
                                                 + k*(m_subunit_unitVecs[n+1][1] - m_subunit_unitVecs[n][1])*(1/m_actualSubLengths[n] + 1/m_actualSubLengths[n+1])
                                                 - k*(m_subunit_unitVecs[n+2][1] - m_subunit_unitVecs[n+1][1])/m_actualSubLengths[n+1];

                }

        }

        return bendingForces;
}


Eigen::MatrixXd Actin::buildDMatrixAniso(double temp, double viscosity, int N)
{
        double KT = (g_Kb*temp);
        Eigen::Matrix2d I;
        I (0,0) = 1;
        I (0,1) = 0;
        I (1,0) = 0;
        I (1,1) = 1;
        Eigen::MatrixXd transDiffMat (N*2,N*2);
        transDiffMat = Eigen::MatrixXd::Zero(N*2,N*2);

        double Yama_corr_X = 0.044;
        double Yama_corr_Y = 1.111;
        double drag_para;
        double drag_perp;

        // j and k are the indices of the 2x2 subblocks
        // they represent the beads
        for (int j = 0; j < N; ++j)
        {
                // Calculate the local tangent vector
                Eigen::Vector2d u_tilde;
                Eigen::MatrixXd u_tilde_T (1, 2);
                if (j > 0 && j < N-1)
                {
                        Eigen::Vector2d u_j;
                        u_j << m_subunit_unitVecs[j+1][0], m_subunit_unitVecs[j+1][1];
                        Eigen::Vector2d u_jm1;
                        u_jm1 << m_subunit_unitVecs[j][0], m_subunit_unitVecs[j][1];

                        u_tilde = (u_j + u_jm1) / (u_j + u_jm1).norm();

                        u_tilde_T << u_tilde(0), u_tilde(1);

                        drag_para = (2*M_PI*viscosity*m_des_subunit_len)/(log(m_length/(2*s_trueRadius))+Yama_corr_X);
                        drag_perp = (4*M_PI*viscosity*m_des_subunit_len)/(log(m_length/(2*s_trueRadius))+Yama_corr_Y);

                }
                else if (j == 0)
                {
                        // end rod, need to adjust length
                        double rodL = 2*m_prescribedSubLengths[0];

                        drag_para = (2*M_PI*viscosity*rodL)/(log(m_length/(2*s_trueRadius))+Yama_corr_X);
                        drag_perp = (4*M_PI*viscosity*rodL)/(log(m_length/(2*s_trueRadius))+Yama_corr_Y);

                        u_tilde <<  m_subunit_unitVecs[j+1][0], m_subunit_unitVecs[j+1][1];
                        u_tilde_T << u_tilde(0), u_tilde(1);

                }
                else
                {
                        // end rod, need to adjust length
                        double rodL = 2*m_prescribedSubLengths[m_num_subunits-1];

                        drag_para = (2*M_PI*viscosity*rodL)/(log(m_length/(2*s_trueRadius))+Yama_corr_X);
                        drag_perp = (4*M_PI*viscosity*rodL)/(log(m_length/(2*s_trueRadius))+Yama_corr_Y);

                        u_tilde <<  m_subunit_unitVecs[N-1][0], m_subunit_unitVecs[N-1][1];
                        u_tilde_T << u_tilde(0), u_tilde(1);
                }
                assert(std::isfinite(m_subunit_unitVecs[j][0]));
                assert(std::isfinite(m_subunit_unitVecs[j][1]));
                assert(std::isfinite(m_subunit_unitVecs[j+1][0]));
                assert(std::isfinite(m_subunit_unitVecs[j+1][1]));

                assert(std::isfinite(u_tilde(0)));
                assert(std::isfinite(u_tilde(1)));

                Eigen::Matrix2d mobKT = KT*((1/drag_para)*(u_tilde*u_tilde.transpose()) + (1/drag_perp)*(I-(u_tilde*u_tilde.transpose()))); // equation 25
                transDiffMat.block<2,2>(j*2,j*2) = mobKT;
        }

        return transDiffMat;
}

Eigen::VectorXd Actin::getRand_dR(double dt, double temp, double viscosity, int N)
{
        // xi 2N vector

        Eigen::VectorXd randDR (2*N);
        double KT = (g_Kb*temp);

        double Yama_corr_X = 0.044;
        double Yama_corr_Y = 1.111;
        double mu_para;
        double mu_perp;

        std::normal_distribution<double> randNumber(0,1);

        for (int j = 0; j < N; ++j)
        {
                // Calculate the local tangent vector
                Eigen::Vector2d u_tilde;
                if (j > 0 && j < N-1)
                {
                        Eigen::Vector2d u_j;
                        u_j << m_subunit_unitVecs[j+1][0], m_subunit_unitVecs[j+1][1];
                        Eigen::Vector2d u_jm1;
                        u_jm1 << m_subunit_unitVecs[j][0], m_subunit_unitVecs[j][1];

                        u_tilde = (u_j + u_jm1) / (u_j + u_jm1).norm();

                        mu_para = (log(m_length/(2*s_trueRadius))+Yama_corr_X)/(2*M_PI*viscosity*m_des_subunit_len);
                        mu_perp = (log(m_length/(2*s_trueRadius))+Yama_corr_Y)/(4*M_PI*viscosity*m_des_subunit_len);

                }
                else if (j == 0)
                {
                        // end rod, need to adjust length
                        double rodL = 2*m_prescribedSubLengths[0];
                        mu_para = (log(m_length/(2*s_trueRadius))+Yama_corr_X)/(2*M_PI*viscosity*rodL);
                        mu_perp = (log(m_length/(2*s_trueRadius))+Yama_corr_Y)/(4*M_PI*viscosity*rodL);
                        u_tilde <<  m_subunit_unitVecs[j+1][0], m_subunit_unitVecs[j+1][1];

                }
                else
                {
                        // end rod, need to adjust length
                        double rodL = 2*m_prescribedSubLengths[m_num_subunits-1];
                        mu_para = (log(m_length/(2*s_trueRadius))+Yama_corr_X)/(2*M_PI*viscosity*rodL);
                        mu_perp = (log(m_length/(2*s_trueRadius))+Yama_corr_Y)/(4*M_PI*viscosity*rodL);
                        u_tilde <<  m_subunit_unitVecs[N-1][0], m_subunit_unitVecs[N-1][1];

                }

                assert(std::isfinite(u_tilde(0)));
                assert(std::isfinite(u_tilde(1)));

                double dxprime = sqrt(2*KT*mu_para*dt)*randNumber(rng.m_mersenne);
                double dyprime = sqrt(2*KT*mu_perp*dt)*randNumber(rng.m_mersenne);


                double dx = dxprime*u_tilde(0) - dyprime*u_tilde(1);
                double dy = dxprime*u_tilde(1) + dyprime*u_tilde(0);

                randDR(j*2) = dx;
                randDR(j*2 + 1) = dy;

        }

        return randDR;
}

Eigen::MatrixXd Actin::getBMatrix(int N)
{
        Eigen::MatrixXd B (N-1, 2*N);
        // fill with zeroes
        B = Eigen::MatrixXd::Zero(N-1,2*N);
        double mag = 0;

        for (int row = 0; row < N-1; ++row)
        {

                int startcol = 2*row;
                int ptI = row + 1; // point index, to avoid end points
                //First pair
                B (row, startcol) = m_points[ptI][0] - m_points[ptI+1][0];
                B (row, startcol+1) = m_points[ptI][1] - m_points[ptI+1][1];

                mag = sqrt( B(row,startcol)*B(row,startcol) + B(row,startcol+1)*B(row,startcol+1) );
                B (row, startcol) = B (row, startcol) / mag;
                B (row, startcol+1) = B (row, startcol+1) / mag;


                // Second pair
                B (row, startcol+2) = m_points[ptI+1][0] - m_points[ptI][0];
                B (row, startcol+3) = m_points[ptI+1][1] - m_points[ptI][1];

                mag = sqrt( B(row,startcol+2)*B(row,startcol+2) + B(row,startcol+3)*B(row,startcol+3) );
                B (row, startcol+2) = B (row, startcol+2) / mag;
                B (row, startcol+3) = B (row, startcol+3) / mag;

        }
        return B;
}

void Actin::addTetherPointNuc()
{
        /*
           Function that is called when required to add a new tethering spring. Give the
           coordiates as inputs
           Should put on middle monomer
           Called during nucleation
         */
        m_tetherMono.push_back(1); // put on 2nd mono (0-based index)
        m_tetherSubunits.push_back(2);

        std::array<double,4> newTether;

        // Tether point half a monomer length from end
        newTether[0] = m_points[2][0];
        newTether[1] = m_points[2][1];
        newTether[2] = newTether[0];
        newTether[3] = newTether[1];

        m_tetherPoints.push_back(newTether);


        m_distToTetBarb = 1;
        m_distToTetPoint = 1;
        m_indexLeadingBarbTet = 0;
        m_indexLeadingPointTet = 0;

}

void Actin::addTetherPointBarb()
{
        /*
           Function that is called when required to add a new tethering spring. Give the
           distance from the barbed end as input
           puts on barbed end
           SIMPLER, PUTS TETHER POINT AT BARBED END
           End causes problems, for now have 0.5 monomer length from end?
           DISCRETE VERSION
         */

        m_distToTetBarb = 0;

        m_tetherSubunits.push_back(m_num_subunits-1);
        m_tetherMono.push_back(m_num_monomers-1);

        std::array<double,4> newTether;


        newTether[0] = m_points[m_num_subunits][0] - (s_monomerLength*0.5)*m_subunit_unitVecs[m_num_subunits-1][0];
        newTether[1] = m_points[m_num_subunits][1] - (s_monomerLength*0.5)*m_subunit_unitVecs[m_num_subunits-1][1];
        newTether[2] = newTether[0];
        newTether[3] = newTether[1];
        m_tetherPoints.push_back(newTether);

        int numTethers = m_tetherMono.size();
        m_indexLeadingBarbTet = numTethers - 1;
        if (numTethers == 1)
        {
                // This is the first tether being put in, so need to calculate properties
                // from the opposite end too
                m_indexLeadingPointTet = 0;
                m_distToTetPoint = m_num_monomers-1;
        }


}

void Actin::addTetherPointPointed()
{
        /*
           Function that is called when required to add a new tethering spring.
           puts on pointed end
           SIMPLER, PUTS TETHER POINT AT pointed END
           End causes problems, for now have 1 monomer length from end?
         */


        m_distToTetPoint = 0;

        m_tetherSubunits.push_back(0);
        m_tetherMono.push_back(0);
        std::array<double,4> newTether;

        newTether[0] = m_points[0][0] + (s_monomerLength*0.5)*m_subunit_unitVecs[0][0];
        newTether[1] = m_points[0][1] + (s_monomerLength*0.5)*m_subunit_unitVecs[0][1];
        newTether[2] = newTether[0];
        newTether[3] = newTether[1];
        m_tetherPoints.push_back(newTether);

        int numTethers = m_tetherMono.size();
        m_indexLeadingPointTet = numTethers - 1;
        if (numTethers == 1)
        {
                // This is the first tether being put in, so need to calculate properties
                // from the opposite end too
                m_indexLeadingBarbTet = 0;
                m_distToTetBarb = m_num_monomers-1;
        }


}


void Actin::remTetherB(std::vector<Actin> &actinVec)
{
        // Removes the end tether (barbed end) upon barbed end depolymerisation

        // Need to remove it from the 3 vectors

        int numTethers = m_tetherMono.size();

        int j = m_indexLeadingBarbTet;


        for (int i = j; i < numTethers-1; ++i)
        {
                m_tetherMono[i] = m_tetherMono[i+1];
                m_tetherSubunits[i] = m_tetherSubunits[i+1];
                m_tetherPoints[i] = m_tetherPoints[i+1];
        }
        m_tetherMono.pop_back();
        m_tetherSubunits.pop_back();
        m_tetherPoints.pop_back();

        // Then need to recalc distance to leading tether points,
        // (if there is at least one left)

        numTethers -= 1;
        if (numTethers > 0)
        {
                // Need to find the nearest tether point to the barbed end

                int tetherMono = 0;

                for (int i = 0; i < numTethers; ++i)
                {
                        if (m_tetherMono[i] >= tetherMono)
                        {
                                m_indexLeadingBarbTet = i;
                        }
                }

                m_distToTetBarb = m_num_monomers - m_tetherMono[m_indexLeadingBarbTet] - 1;
                // it should be equal to 0 if its on the end
        }
        else
        {
                // No tethers left, what do we do with distance to tether point?
                m_distToTetBarb = m_num_monomers;
                m_distToTetPoint = m_num_monomers;

        }

}

void Actin::remTetherP(std::vector<Actin> &actinVec)
{
        // Removes the end tether (pointed end) upon pointed end depolymerisation
        // Need to remove it from the 3 vectors

        int numTethers = m_tetherMono.size();

        int j = m_indexLeadingPointTet;


        for (int i = j; i < numTethers-1; ++i)
        {
                m_tetherMono[i] = m_tetherMono[i+1];
                m_tetherSubunits[i] = m_tetherSubunits[i+1];
                m_tetherPoints[i] = m_tetherPoints[i+1];
        }
        m_tetherMono.pop_back();
        m_tetherSubunits.pop_back();
        m_tetherPoints.pop_back();

        // Then need to recalc distance to leading tether points,
        // (if there is at least one left)

        numTethers -= 1;
        if (numTethers > 0)
        {
                // Need to find the nearest tether point to the pointed end

                int tetherMono = m_num_monomers-1;

                for (int i = 0; i < numTethers; ++i)
                {
                        //std::cout << m_tetherSubunits[i] << std::endl;
                        if (m_tetherMono[i] <= tetherMono)
                        {
                                m_indexLeadingPointTet = i;
                        }
                }



                m_distToTetPoint = m_tetherMono[m_indexLeadingPointTet];
                // equal to 0 if on the end
        }
        else
        {
                // No tethers left, what do we do with distance to tether point?

                m_distToTetBarb = m_num_monomers;
                m_distToTetPoint = m_num_monomers;

        }

}

void Actin::updateTetherLoc()
{
        // Function that uses the tetherMono to update the coordinates of the tether
        // point
        // Called after the filament has moved

        int numTethers = m_tetherMono.size();

        for (int i = 0; i < numTethers; ++i)
        {
                int j = m_tetherSubunits[i];
                double lenAlongSub;
                int subid = findSubunit(m_tetherMono[i], lenAlongSub);
                assert(subid == j);
                m_tetherPoints[i][2] = m_points[j][0] + lenAlongSub*m_subunit_unitVecs[j][0];
                m_tetherPoints[i][3] = m_points[j][1] + lenAlongSub*m_subunit_unitVecs[j][1];
                assert(std::isfinite(m_tetherPoints[i][2]));
                assert(std::isfinite(m_tetherPoints[i][3]));

        }
}

void Actin::shiftTetherMonosUp()
{
        // Called on the filament, after pointed polymerisation
        // Any tethers will have thier mono id moved
        // up by 1
        int numTethers = getNumTethers();
        for (int i = 0; i < numTethers; ++i)
        {
                m_tetherMono[i]++;
        }
}

void Actin::shiftTetherMonosDown()
{
        int numTethers = getNumTethers();
        for (int i = 0; i < numTethers; ++i)
        {
                m_tetherMono[i]--;
        }
}

void Actin::checkAndMoveTetherDown(int subOne, int subTwo, bool checkAll)
{
        /*
           General function that has input of two subids (subOne and subTwo)
           which are the subunit of the fila where the tethers are that will be checked
           Checks all tethers on the subs and determines if they need to move DOWN.
           subOne and subTwo are defaulted to -1 and -2 respectively (which are invalid)
           Therefore if no or just one argument is given the other is invalid and ignored
           checkAll is defaulted to true
         */

        int numTethers = m_tetherMono.size();

        for (int i = 0; i < numTethers; ++i)
        {
                int j = m_tetherSubunits[i];
                if (j != subOne && j!= subTwo && !checkAll)
                        continue;

                // j is the old sub


                double junk;
                int newTetherSub = findSubunit(m_tetherMono[i], junk);

                if (j != newTetherSub)
                {

                        assert(newTetherSub == j-1);
                        // The tether subunit has moved down
                        setTetherSub(i, j-1);
                }
        }
}

void Actin::checkAndMoveTetherUp(int subOne, int subTwo, bool checkAll)
{
        /*
           General function that has input of two subids (subOne and subTwo)
           which are the subunit of the fila where the tethers are that will be checked
           Checks all tethers on the subs and determines if they need to move UP.
           subOne and subTwo are defaulted to -1 and -2 respectively (which are invalid)
           Therefore if just one argument is given the other is invalid and ignored
         */

        int numTethers = m_tetherMono.size();

        for (int i = 0; i < numTethers; ++i)
        {
                int j = m_tetherSubunits[i];
                if (j != subOne && j!= subTwo && !checkAll)
                        continue;

                // j is the old sub


                double junk;
                int newTetherSub = findSubunit(m_tetherMono[i], junk);

                if (j != newTetherSub)
                {
                        assert(newTetherSub == j+1);
                        // The tether subunit has moved up
                        setTetherSub(i, j+1);
                }
        }
}


void Actin::moveTetherSub_down()
{
        // Check for moving tethers
        /*
           This is for either barbed end poly or pointed end depoly, when the tether
           can move down a subunit.
           Called for when we are in initial straight regime
           Check ALL tethers along filament
         */
        checkAndMoveTetherDown(); // use default arguments

}

void Actin::moveTetherSub_up()
{
        // Check for moving tethers
        /*
           This is for either barbed end depoly or pointed end poly, when the branch
           can move up a subunit.
           Called for when we are in initial straight regime
           Check ALL tethers along filament
         */

        checkAndMoveTetherUp(); // use default arguments

}

void Actin::moveTetherPoint_poly_barb()
{
        // Check for moving tether points
        // After a polymerisation event, check any tether points on the final subunit
        // If any have "moved", they should be moved down to the penultimate subunit

        checkAndMoveTetherDown(m_num_subunits-1,-2,false);

}

void Actin::moveTetherPoint_poly_point()
{
        // Check for moving tether points
        // After a polymerisation event, check any tether points on the first subunit
        // If any have "moved", they should be moved up to the first subunit

        checkAndMoveTetherUp(0,-2,false);

}

void Actin::moveTetherPoint_depoly_barb()
{
        // Check for moving tether points
        // After a depolymerisation event, check any tether points on the penultimate
        // subunit. If any have moved we know to move them up to the final sub

        checkAndMoveTetherUp(m_num_subunits-2,-2,false);

}

void Actin::moveTetherPoint_depoly_point()
{
        // Check for moving tether points
        // After a depolymerisation event, check any tether points on the second
        // subunit. If any have moved we know to move them down to the first sub
        checkAndMoveTetherDown(1,-2,false);

}

void Actin::moveTetherPoint_Add_barb()
{
        // Check for moving tether points
        // After the addition of a barbed end point,
        // check any tether points on the penultimate and pen-penultimate subunit
        // If any have "moved", they should be moved up one

        checkAndMoveTetherUp(m_num_subunits-3,m_num_subunits-2,false);

}

void Actin::moveTetherPoint_Add_point()
{
        // Check for moving tether points
        // After the addition of a pointed end point,
        // check any tether points on the second and third subunit
        // If any have "moved", they should be moved down one
        checkAndMoveTetherDown(1,2,false);

}

void Actin::moveTetherPoint_Rem_barb()
{
        // Check for moving tether points
        // After the removal of a barbed end point,
        // check any tether points on the penultimate sub
        // If any have "moved", they should be moved down one
        checkAndMoveTetherDown(m_num_subunits-2,m_num_subunits-1,false);

}

void Actin::moveTetherPoint_Rem_point()
{
        // Check for moving tether points
        // After the removal of a pointed end point,
        // check any tether points on the first sub
        // If any have "moved", they should be moved up one
        checkAndMoveTetherUp(0,1,false);
}

Eigen::VectorXd Actin::calcTetherForces(int N, bool tethering)
{
        // Forces that will be applied to each end point of our membraneFilament
        // Just uses Hookes Law F = -k(r_spring - r_point)
        // Returns 2 doubles per spring, (Fx, Fy)
        // These should then be added to the bending force

        Eigen::VectorXd tetherF (2*N);
        tetherF.setZero();
        if (!tethering)
            return tetherF;

        int numTethers = m_tetherMono.size();

        for (int i = 0; i < numTethers; ++i)
        {
                int j = m_tetherSubunits[i];
                if (j == 0)
                {
                        tetherF(0) = tetherF(0) + s_tetherStiff*(m_tetherPoints[i][0] - m_tetherPoints[i][2]);
                        tetherF(1) = tetherF(1) + s_tetherStiff*(m_tetherPoints[i][1] - m_tetherPoints[i][3]);

                }
                else if (j == m_num_subunits-1)
                {
                        tetherF((2*N)-2) = tetherF((2*N)-2) + s_tetherStiff*(m_tetherPoints[i][0] - m_tetherPoints[i][2]);
                        tetherF((2*N)-1) = tetherF((2*N)-1) + s_tetherStiff*(m_tetherPoints[i][1] - m_tetherPoints[i][3]);
                }
                else
                {
                        double lenAlongSub;
                        int subid = findSubunit(m_tetherMono[i], lenAlongSub);
                        assert (subid == j);
                        double a = lenAlongSub / m_prescribedSubLengths[j];
                        assert (a < 1);
                        // pointed end forces
                        tetherF((j-1)*2) = tetherF((j-1)*2) + (1-a)*(s_tetherStiff*(m_tetherPoints[i][0] - m_tetherPoints[i][2]));
                        tetherF(((j-1)*2)+1) = tetherF(((j-1)*2)+1) + (1-a)*(s_tetherStiff*(m_tetherPoints[i][1] - m_tetherPoints[i][3]));

                        // barbed end forces
                        tetherF((j)*2) = tetherF((j)*2) + a*(s_tetherStiff*(m_tetherPoints[i][0] - m_tetherPoints[i][2]));
                        tetherF(((j)*2)+1) = tetherF(((j)*2)+1) + a*(s_tetherStiff*(m_tetherPoints[i][1] - m_tetherPoints[i][3]));
                }
        }

        return tetherF;
}

int Actin::findSubunit(int monoID, double &lenAlongSub)
{
        /* Function that returns the id of the subunit containing the monomer with
           the id of monoID
           returns by reference the length along the sub too
         */

        double lengthAlongParent = monoID*s_monomerLength + (0.5*s_monomerLength);

        for (int i = 0; i < m_num_subunits; ++i)
        {
                // go along the filament and minus the subunit lengths to find which
                // subunit has the branch point
                // Use prescibed lengths for this because we are comparing to the
                // prescribed total length
                lengthAlongParent -= m_prescribedSubLengths[i];

                if (lengthAlongParent < -1E-14) // avoid floating point error , should be 0
                {
                        lengthAlongParent += m_prescribedSubLengths[i];
                        lenAlongSub = lengthAlongParent;
                        if (lenAlongSub < 0)
                        {
                                lenAlongSub = 0;
                        }
                        return i;
                }
        }
        std::cout << "Code should not get here! Actin.cpp line 7997" << std::endl;
        std::cout << lengthAlongParent << std::endl;
        std::cout << monoID << " : " << m_num_monomers << std::endl;
        assert(monoID < m_num_monomers);
        exit(1);
}

std::array<double,2> Actin::findMonoCoord(int monoID)
{
        /* Function that returns the 2D spatial coordinate corresponding to the
           middle of the monomer monoID
         */
        //assert(monoID < m_num_monomers);
        double lengthAlongParent = monoID*s_monomerLength + (0.5*s_monomerLength);
        for (int i = 0; i < m_num_subunits; ++i)
        {
                // go along the filament and minus the subunit lengths to find which
                // subunit has the branch point
                // Use prescibed lengths for this because we are comparing to the
                // prescribed total length
                lengthAlongParent -= m_prescribedSubLengths[i];
                if (lengthAlongParent < -1E-14) // avoid floating point error , should be 0
                {
                        lengthAlongParent += m_prescribedSubLengths[i];
                        std::array<double, 2> coords;
                        coords[0] = m_points[i][0] + lengthAlongParent*m_subunit_unitVecs[i][0];
                        coords[1] = m_points[i][1] + lengthAlongParent*m_subunit_unitVecs[i][1];
                        return coords;
                }
        }
        std::cout << "Code should not get here! Actin.cpp line 6156 ";
        std::cout << lengthAlongParent << std::endl;
        std::cout << m_num_monomers << ", " << monoID << std::endl;
        std::cout << m_id << std::endl;
        exit(1);
}

std::array<double,4> Actin::findMonoStartEndCoord(int monoID)
{
    std::array<double,2> monoCentre = findMonoCoord(monoID);
    double junk;
    int subid = findSubunit(monoID, junk);

    double monoX1 = monoCentre[0] - (m_subunit_unitVecs[subid][0]*(s_monomerLength/2));
    double monoY1 = monoCentre[1] - (m_subunit_unitVecs[subid][1]*(s_monomerLength/2));
    double monoX2 = monoCentre[0] + (m_subunit_unitVecs[subid][0]*(s_monomerLength/2));
    double monoY2 = monoCentre[1] + (m_subunit_unitVecs[subid][1]*(s_monomerLength/2));

    std::array<double,4> monoStartEnd = {monoX1, monoY1, monoX2, monoY2};
    return monoStartEnd;
}

void Actin::shiftAvailMonosUp()
{
        // function that is called on pointed end polymerisation that adds 1 to each
        // value in m_availMonos.
        for (unsigned int i = 0; i < m_availMonos.size(); ++i)
        {
                m_availMonos[i] += 1;
        }

}

void Actin::addAvailMono()
{
        // Add available monomer if it is needed
        m_availMonos.push_back(0); // push back garbage
        for (unsigned int i = m_availMonos.size()-1; i > 0; --i)
        {
                m_availMonos[i] = m_availMonos[i-1];
        }
        m_availMonos[0] = s_maxSpacing;
}

void Actin::shiftAvailMonosDown()
{
        // function that is called on pointed end DEpolymerisation that subtracts 1
        // from each value in m_availMonos.
        for (unsigned int i = 0; i < m_availMonos.size(); ++i)
        {
                m_availMonos[i] -= 1;
        }

        if (m_availMonos.size() > 0 && m_availMonos[0] == s_maxSpacing-1)
        {
                // need to remove from vector
                m_availMonos.erase(m_availMonos.begin());
        }
}

void Actin::adjustBranchCompat(int monoID)
{
        // Function that is called on branching, changes all neighboring monomers to
        // be incompatible for branching

        int upperLim;
        if (monoID+s_branchSpacing >= m_num_monomers)
        {
                upperLim = m_num_monomers;
        }
        else
        {
                upperLim = monoID+s_branchSpacing+1;
        }

        for (int i = monoID-s_branchSpacing; i < upperLim; ++i)
        {
                if (i < 0)
                {
                        continue;
                }
                assert(i < m_num_monomers);

                if (m_monosCompat[i] == true)
                {
                        // it needs to change to false
                        m_monosCompat[i] = false;

                        // find monomer in m_availMonos and remove it
                        for (unsigned int j = 0; j < m_availMonos.size(); ++j)
                        {
                                if (m_availMonos[j] == i)
                                {

                                        // found the monomer, its located at i,j
                                        // need to remove this from the vector
                                        m_availMonos.erase(m_availMonos.begin()+j);
                                        break;
                                }
                        }
                }


        }

}

void Actin::recalcDistToBR(std::vector<Actin> &actinVec, const int nActin)
{
        // When creating a new branch the distance from the barbed end of the mother
        // to the branch point *may* need to be recalculated (it will if the new
        // branch is the leading branch)

        // Called on the mother
        // The new branch will be the last one in the actinVec

        int branchMono = actinVec[nActin-1].getMotherMonoID();
        int barbedMono = m_num_monomers-1;

        if (barbedMono-branchMono < m_distToLeadBR)
        {
                m_distToLeadBR = barbedMono-branchMono;
        }


}

void Actin::recalcDistToBR2(std::vector<Actin> &actinVec)
{
        // Version of the above function but called on autodebranch barb

        // Need to find the NEW leading branch


        if (m_daughter_num == 0)
        {
                m_distToLeadBR = m_num_monomers;
        }
        else
        {
                int leadingBR = 0;
                for (int i = 0; i < m_num_subunits; ++i)
                {
                        for (unsigned int j = 0; j < m_daughterIDs[i].size(); ++j)
                        {
                                int k = m_daughterIDs[i][j]; // id of daughter
                                if (actinVec[k].getMotherMonoID() > leadingBR)
                                {
                                        leadingBR = actinVec[k].getMotherMonoID();
                                }
                        }
                }
                m_distToLeadBR = m_num_monomers - leadingBR -1;
        }

}

std::vector<int> Actin::getUnavailMonos(std::vector<Actin> &actinVec)
{
  /*
   * Calculates a list of all the monomers that are unavailable to either branch
   * Or crosslink
   * i is the id of the filament
   */
    std::vector<int> unAvailMonos;

    int end = (s_maxSpacing < m_num_monomers) ? s_maxSpacing : m_num_monomers;
    for (int i = 0; i < end; ++i)
    {
        unAvailMonos.push_back(i);
    }

    for (int j = 0; j < m_num_subunits; ++j)
    {
        for (unsigned int k = 0; k < getBranchIDSubVector(j).size(); ++k)
        {
            int id = getBranchIDSubVector(j)[k];

            int monoID = actinVec[id].getMotherMonoID();
            int lowLim = (monoID-Actin::s_branchSpacing < Actin::s_maxSpacing) ? Actin::s_maxSpacing : monoID-Actin::s_branchSpacing;
            int upLim = (monoID+Actin::s_branchSpacing >= getNumMonomers()) ? getNumMonomers() : monoID+Actin::s_branchSpacing+1;

            for (int l = lowLim; l < upLim; ++l)
            {
                if (std::find(unAvailMonos.begin(), unAvailMonos.end(), l) == unAvailMonos.end())
                {
                    // if not already here add it
                    unAvailMonos.push_back(l);
                }
            }

        }
    }

    // Crosslinking
    for (int j = 0; j < getNumCLinks(); ++j)
    {
        // Make sure any linked actin still point to the right cl idx
        int otherFila = getCLinkActinAndSite(j)[2];
        int otherFilaCLIDX = getCLinkActinAndSite(j)[3];
        assert(otherFilaCLIDX < actinVec[otherFila].getNumCLinks());
        int otherFilaotherIDX = actinVec[otherFila].getCLinkActinAndSite(otherFilaCLIDX)[3];

        assert(j == otherFilaotherIDX);

        // Check availability vector - that it is unavailable
        int monoID = getCLinkActinAndSite(j)[0];
        int lowLim = (monoID-Actin::s_cLSpacing < Actin::s_maxSpacing) ? Actin::s_maxSpacing : monoID-Actin::s_cLSpacing;
        int upLim = (monoID+Actin::s_cLSpacing >= m_num_monomers) ? m_num_monomers : monoID+Actin::s_cLSpacing+1;

        for (int k = lowLim; k < upLim; ++k)
        {
            if (std::find(unAvailMonos.begin(), unAvailMonos.end(), k) == unAvailMonos.end())
            {
                // if not already here add it
                unAvailMonos.push_back(k);
            }
        }

    }

    return unAvailMonos;
}

void Actin::recalcMonoCompat(std::vector<Actin> &actinVec, int monoID)
{
    // Called after either debranching or unlinking

    int checkStart = (monoID - s_maxSpacing*2 > 0) ? monoID - s_maxSpacing*2 : 0;
    int checkEnd = (monoID + s_maxSpacing*2 < m_num_monomers) ? monoID + s_maxSpacing*2 +1 : m_num_monomers;
    std::vector<int> unAvailMonos = getUnavailMonos(actinVec);
    std::sort(unAvailMonos.begin(), unAvailMonos.end());

    for (int i = checkStart; i < checkEnd; ++i)
    {
        if (m_monosCompat[i] == false && (std::find(unAvailMonos.begin(), unAvailMonos.end(), i) == unAvailMonos.end()))
        {
            // the monomer is designated unavailable but it is now available!
            m_monosCompat[i] = true;
            m_availMonos.push_back(i);
        }
    }

    sortAvailMonos(); // get these in correct order!
}

void Actin::recalcMonoCompatALL(std::vector<Actin> &actinVec)
{
    // Called during severing
    // whole filament
    std::vector<int> unAvailMonos = getUnavailMonos(actinVec);
    std::sort(unAvailMonos.begin(), unAvailMonos.end());

    for (int i = 0; i < m_num_monomers; ++i)
    {
        if (m_monosCompat[i] == false && (std::find(unAvailMonos.begin(), unAvailMonos.end(), i) == unAvailMonos.end()))
        {
            // the monomer is designated unavailable but it is now available!
            m_monosCompat[i] = true;
            m_availMonos.push_back(i);
        }

    }

    sortAvailMonos(); // get these in correct order!
}

void Actin::sortAvailMonos()
{
        std::sort(m_availMonos.begin(), m_availMonos.end());
}




void Actin::shiftMotherMonosUp(std::vector<Actin> &actinVec)
{
        // Called on the mother, after pointed polymerisation
        // Any branches attached to the mother will have their motherMonomer id moved
        // up by 1
        for (int i = 0; i < m_num_subunits; ++i)
        {
            for (unsigned int j = 0; j < m_daughterIDs[i].size(); ++j)
            {
                    int id = m_daughterIDs[i][j];
                    actinVec[id].incrementMotMono();
            }
        }

}

void Actin::shiftMotherMonosDown(std::vector<Actin> &actinVec)
{
        for (int i = 0; i < m_num_subunits; ++i)
        {
            for (unsigned int j = 0; j < m_daughterIDs[i].size(); ++j)
            {
                    int id = m_daughterIDs[i][j];
                    actinVec[id].decrementMotMono();
            }
        }

}

void Actin::reduceID(std::vector<Actin> &actinVec)
{
        // Upon deletion or something the id of this filament needs to be reduced by
        // 1. But also the parent id of any daughters of this filament has also has to be
        // reduced by 1.

        m_id -= 1;
        for (int i = 0; i < m_num_subunits; ++i)
        {
                for (unsigned int j = 0; j < m_daughterIDs[i].size(); ++j)
                {
                        int daughtID = m_daughterIDs[i][j];
                        actinVec[daughtID].reduceParentID();
                }
        }


}

void Actin::checkandReduceDaughtID(int doomedID)
{
        for (int i = 0; i < m_num_subunits; ++i)
        {
                for (unsigned int j = 0; j < m_daughterIDs[i].size(); ++j)
                {
                        int daughtID = m_daughterIDs[i][j];
                        if (daughtID > doomedID)
                        {
                                // reduce
                                m_daughterIDs[i][j]--;
                        }
                }
        }
}

void Actin::createCrossLink(Actin &otherFila, int monoID, int subID,
                            int otherMonoID, int otherSubID,
                            std::vector<Actin> &actinVec)
{
    /*
     * Function to create a link between *this and otherFila via monoID and
     * other monoID
     */

     addToLinkVect(monoID, subID, otherFila.getID(), otherFila.getNumCLinks(), otherMonoID, otherSubID);
     otherFila.addToLinkVect(otherMonoID, otherSubID, m_id, getNumCLinks()-1, monoID, subID);


     // Remove from available vector
     remMonoFromAvail(monoID);
     otherFila.remMonoFromAvail(otherMonoID);


     addEmptyToCLPoints();
     otherFila.addEmptyToCLPoints();


     recalcDistToLeadCL(monoID);
     otherFila.recalcDistToLeadCL(otherMonoID);


     if (otherFila.getStructureID() != m_structure_id || !checkForMasterInStructure(actinVec))
     {
        // set "otherFila" structure id to m_structure_id, but will need to do any branches
        otherFila.changeStructureALL(actinVec, m_structure_id);
        otherFila.setMasterBoolFalse(actinVec);
        findNewMaster(actinVec);
        // have to make sure any other filas connected to otherFila also change
     }

     m_changed = true;
     otherFila.setChangedBoolTrue();

     updateCrossLinkLoc(actinVec);
     otherFila.updateCrossLinkLoc(actinVec);

}

void Actin::setMasterBoolFalse(std::vector<Actin> &actinVec)
{
        // function that changes all filaments in structure to non-masters

        for (unsigned int i = 0; i < actinVec.size(); ++i)
        {
            if (actinVec[i].getStructureID() == m_structure_id)
            {
                actinVec[i].turnOffCLinkMaster();
            }
        }

}

bool Actin::checkForMasterInStructure(std::vector<Actin> &actinVec)
{
        // function that scans the whole structure and checks to see if there is
        // a master

        int numMasters = 0;
        for (unsigned int i = 0; i < actinVec.size(); ++i)
        {
            if (actinVec[i].getStructureID() == m_structure_id)
            {
                if (actinVec[i].getCLinkMaster())
                {
                    //return true;
                    assert(actinVec[i].getNumCLinks() > 0);

                    assert(actinVec[i].getLength() < 2*s_segmentationLength);
                    numMasters += 1;
                }
            }
        }

        if (numMasters == 1)
        {
            return true;
        }
        else if (numMasters == 0)
        {
            return false;
        }
        else
        {
            std::cout << "Error, more than one master in structure!" << std::endl;
            std::cout << numMasters << std::endl;
            exit(1);
        }

}

bool Actin::checkForCrossLinksInStructure(std::vector<Actin> &actinVec)
{
    // function that scans the whole structure and checks to see if there is
    // a single crosslink
    if (getNumCLinks() > 0)
    {
        return true;
    }

    std::vector<int> checkedFilas;
    checkedFilas.push_back(m_id); // this call!

    for (unsigned int i = 0; i < actinVec.size(); ++i)
    {
        if (m_id == actinVec[i].getParentID() || actinVec[i].getID() == m_parent_id)
        {
            if (actinVec[i].getNumCLinks() > 0)
            {
                return true;
            }
            else
            {
                // Need to use recursion here for grandaughters etc
                if (actinVec[i].checkForCrossLinksInStructure(actinVec, checkedFilas))
                {
                    return true;
                }
            }
        }
    }

    return false;
}

bool Actin::checkForCrossLinksInStructure(std::vector<Actin> &actinVec,
                                          std::vector<int> &checkedFilas)
{
    // function that scans the whole structure and checks to see if there is
    // a single crosslink

    if (std::find(checkedFilas.begin(), checkedFilas.end(), m_id) != checkedFilas.end())
    {
        // We have done this filament already, so return
        return false;
    }
    checkedFilas.push_back(m_id); // this call!


    if (getNumCLinks() > 0)
    {
        return true;
    }


    for (unsigned int i = 0; i < actinVec.size(); ++i)
    {
        if (m_id == actinVec[i].getParentID() || actinVec[i].getID() == m_parent_id)
        {
            checkedFilas.push_back(i);
            if (actinVec[i].getNumCLinks() > 0)
            {
                return true;
            }
            else
            {
                // Need to use recursion here for grandaughters etc
                if (actinVec[i].checkForCrossLinksInStructure(actinVec, checkedFilas))
                {
                    return true;
                }
            }
        }
    }

    return false;
}


bool Actin::findNewMaster(std::vector<Actin> &actinVec)
{
    // function that scans the whole structure and checks to see if there is
    // a single crosslink and that filament is rigid,
    // first one it hits is set to be the master
    // BRANCHES CANT BE MASTERS!

    if (getNumCLinks() > 0 && m_length < 2*s_segmentationLength && m_parent_id == -1)
    {
        m_cLinkMaster = true;
        return true;
    }

    std::vector<int> checkedFilas;
    checkedFilas.push_back(m_id); // this call!

    for (unsigned int i = 0; i < actinVec.size(); ++i)
    {
        if (actinVec[i].getID() == m_parent_id || m_id == actinVec[i].getParentID())
        {
            if (actinVec[i].getNumCLinks() > 0 && actinVec[i].getLength() < 2*s_segmentationLength && actinVec[i].getParentID() == -1)
            {
                actinVec[i].setSelfMasterBoolTrue();
                return true;
            }
            else
            {
                // Need to use recursion here for grandaughters etc
                if (actinVec[i].findNewMaster(actinVec, checkedFilas))
                {
                    return true;
                }
            }
        }
    }

    for (int i = 0; i < getNumCLinks(); ++i)
    {
        // if has crosslinks but flexible
        int otherFila = m_cLinkActinAndSites[i][2];
        if (actinVec[otherFila].getLength() < 2*s_segmentationLength && actinVec[otherFila].getParentID() == -1)
        {
            actinVec[otherFila].setSelfMasterBoolTrue();
            return true;
        }
        else
        {
            if (actinVec[otherFila].findNewMaster(actinVec, checkedFilas))
            {
                return true;
            }
        }
    }

    return false;
}

bool Actin::findNewMaster(std::vector<Actin> &actinVec,
                          std::vector<int> &checkedFilas)
{
    // function that scans the whole structure and checks to see if there is
    // a single crosslink and that filament is rigid,
    // first one it hits is set to be the master

    checkedFilas.push_back(m_id); // this call!

    for (unsigned int i = 0; i < actinVec.size(); ++i)
    {
        if (actinVec[i].getID() == m_parent_id || m_id == actinVec[i].getParentID())
        {
            if (std::find(checkedFilas.begin(), checkedFilas.end(), i) != checkedFilas.end())
            {
                // We have checked this filament before!
                continue;
            }
            checkedFilas.push_back(i);

            if (actinVec[i].getNumCLinks() > 0 && actinVec[i].getLength() < 2*s_segmentationLength && actinVec[i].getParentID() == -1)
            {
                actinVec[i].setSelfMasterBoolTrue();
                return true;
            }
            else
            {
                // Need to use recursion here for grandaughters etc
                if (actinVec[i].findNewMaster(actinVec, checkedFilas))
                {
                    return true;
                }
            }
        }
    }

    for (int i = 0; i < getNumCLinks(); ++i)
    {
        // if has crosslinks but flexible
        int otherFila = m_cLinkActinAndSites[i][2];
        if (std::find(checkedFilas.begin(), checkedFilas.end(), otherFila) != checkedFilas.end())
        {
            // We have checked this filament before!
            continue;
        }

        if (actinVec[otherFila].getLength() < 2*s_segmentationLength && actinVec[otherFila].getParentID() == -1)
        {
            actinVec[otherFila].setSelfMasterBoolTrue();
            return true;
        }
        else
        {
            if (actinVec[otherFila].findNewMaster(actinVec, checkedFilas))
            {
                return true;
            }
        }
    }

    return false;
}

bool Actin::checkForSameStructure(std::vector<Actin> &actinVec, int otherFilaID)
{
        // function that finds the all the daughters AND PARENTS of the filament that has
        // debranched and sets their structureID to that of the parent
        std::vector<int> checkedFilas;
        checkedFilas.push_back(m_id); // this call!

        for (unsigned int i = 0; i < actinVec.size(); ++i)
        {
            if (i == otherFilaID)
            {
                // Found it!
                return true;
            }

            if (m_id == actinVec[i].getParentID() || actinVec[i].getID() == m_parent_id)
            {
                    checkedFilas.push_back(i);
                    // Need to use recursion here for grandaughters etc
                    actinVec[i].checkForSameStructure(actinVec, otherFilaID,
                                                      checkedFilas);
            }
        }

        for (int i = 0; i < getNumCLinks(); ++i)
        {
            int otherFila = m_cLinkActinAndSites[i][2];
            if (otherFila == otherFilaID)
            {
                // Found it!
                return true;
            }

            checkedFilas.push_back(otherFila);
            actinVec[otherFila].checkForSameStructure(actinVec, otherFilaID, checkedFilas);
        }
        return false;
}



bool Actin::checkForSameStructure(std::vector<Actin> &actinVec, int otherFilaID,
                                  std::vector<int> &checkedFilas)
{
        // function that finds the all the daughters of the filament that has
        // debranched and sets their structureID to that of the parent

        for (unsigned int i = 0; i < actinVec.size(); ++i)
        {
            if (m_id == actinVec[i].getParentID() || actinVec[i].getID() == m_parent_id)
            {
                  if (i == otherFilaID)
                  {
                      // Found it!
                      return true;
                  }

                  if (std::find(checkedFilas.begin(), checkedFilas.end(), i) != checkedFilas.end())
                  {
                          // We have changed this filament before!
                          continue;
                  }
                  checkedFilas.push_back(i);
                  // Need to use recursion here for grandaughters etc
                  actinVec[i].checkForSameStructure(actinVec, otherFilaID, checkedFilas);
            }
        }

        for (int i = 0; i < getNumCLinks(); ++i)
        {
            int otherFila = m_cLinkActinAndSites[i][2];
            if (otherFila == otherFilaID)
            {
                // Found it!
                return true;
            }

            if (std::find(checkedFilas.begin(), checkedFilas.end(), otherFila) != checkedFilas.end())
            {
                    // We have changed this filament before!
                    continue;
            }

            checkedFilas.push_back(otherFila);
            actinVec[otherFila].checkForSameStructure(actinVec, otherFilaID, checkedFilas);
        }

        return false;

}

int Actin::countCrossLinksInCLStruct(std::vector<Actin> &actinVec)
{
  /*
   * Scans the "crosslink" structure for crosslinks BUT only crosslinks that
   * connect to other "crosslink" structures
   * Used for rotation of the "mini-ellipse"
   */

    int numCLinks = 0;

    for (unsigned int i = 0; i < actinVec.size(); ++i)
    {
        if (actinVec[i].getCLStructureID() == m_crossStructure_id)
        {
            // Same "crosslink" structure
            for (int j = 0; j < actinVec[i].getNumCLinks(); ++j)
            {
                // Scan all crosslinks on this fila
                int otherFilaID = actinVec[i].getCLinkActinAndSite(j)[2];
                if (actinVec[otherFilaID].getCLStructureID() != m_crossStructure_id)
                {
                    // The crosslink connects to another structure!
                    numCLinks += 1;
                }
            }
        }
    }

    return numCLinks;
}

std::array<double,3> Actin::getCLCoordInCLStruct(std::vector<Actin> &actinVec)
{
  /*
   * Similar to above but returns the coordinate
   */

   for (unsigned int i = 0; i < actinVec.size(); ++i)
   {
       if (actinVec[i].getCLStructureID() == m_crossStructure_id)
       {
           // Same "crosslink" structure
           for (int j = 0; j < actinVec[i].getNumCLinks(); ++j)
           {
               // Scan all crosslinks on this fila
               int otherFilaID = actinVec[i].getCLinkActinAndSite(j)[2];
               if (actinVec[otherFilaID].getCLStructureID() != m_crossStructure_id)
               {
                   // The crosslink connects to another structure!
                   std::array<double,3> cLCoord = {actinVec[i].getCLinksPoint(0)[0], actinVec[i].getCLinksPoint(0)[1], 0};
                   return cLCoord;
               }
           }
       }
   }

   std::cout << "Error code should not get here in function getCLCoordInCLStruct" << std::endl;
   exit(1);

}

void Actin::addToLinkVect(int mono, int subid, int otherID, int otherCLIDX,
                          int otherMono, int otherSubID)
{

    std::array<int,6> tmp;
    tmp[0] = mono;
    tmp[1] = subid;
    tmp[2] = otherID;
    tmp[3] = otherCLIDX; // the index of this crosslink on the other filament
    tmp[4] = otherMono;
    tmp[5] = otherSubID;
    m_cLinkActinAndSites.push_back(tmp);
}

void Actin::setMonoToUnavailable(int mono)
{
    // Go below by s_cLSpacing and above by s_cLSpacing
    int start = mono - s_cLSpacing;
    if (start < 0)
        start = 0;

    int end = mono + s_cLSpacing;
    if (end >= m_num_monomers)
        end = m_num_monomers-1;

    for (int i = start; i <= end; ++i)
    {
        m_monosCompat[i] = false;
    }

}

void Actin::remMonoFromAvail(int mono)
{

    int upperLim;
    if (mono+s_cLSpacing >= m_num_monomers)
    {
            upperLim = m_num_monomers;
    }
    else
    {
            upperLim = mono+s_cLSpacing+1;
    }

    for (int i = mono-s_cLSpacing; i < upperLim; ++i)
    {
            if (i < 0)
            {
                    continue;
            }
            assert(i < m_num_monomers);
            if (m_monosCompat[i] == true)
            {
                // it needs to change to false
                m_monosCompat[i] = false;

                // find monomer in m_availMonoBranch and remove it
                for (unsigned int j = 0; j < m_availMonos.size(); ++j)
                {
                    if (m_availMonos[j] == i)
                    {
                        // found the monomer, its located at i,j
                        // need to remove this from the vector
                        m_availMonos.erase(m_availMonos.begin()+j);
                        break;
                    }
                }
            }


    }

}

void Actin::shiftClMonosUp(std::vector<Actin> &actinVec)
{
        // Called on the filament, after pointed polymerisation
        // Any tethers will have thier mono id moved
        // up by 1
        int numCrossLinks = m_cLinkActinAndSites.size();
        for (int i = 0; i < numCrossLinks; ++i)
        {
                m_cLinkActinAndSites[i][0]++;
                int otherID = m_cLinkActinAndSites[i][2];
                int otherIDX = m_cLinkActinAndSites[i][3];
                actinVec[otherID].shiftClMonoUp(otherIDX);
        }
}

void Actin::shiftClMonoUp(int clIDX)
{
  m_cLinkActinAndSites[clIDX][4]++;
}


void Actin::shiftClMonosDown(std::vector<Actin> &actinVec)
{
        // Called on the filament, after pointed polymerisation
        // Any tethers will have thier mono id moved
        // up by 1
        int numCrossLinks = m_cLinkActinAndSites.size();
        for (int i = 0; i < numCrossLinks; ++i)
        {
                m_cLinkActinAndSites[i][0]--;
                int otherID = m_cLinkActinAndSites[i][2];
                int otherIDX = m_cLinkActinAndSites[i][3];
                actinVec[otherID].shiftClMonoDown(otherIDX);
        }
}

void Actin::shiftClMonoDown(int clIDX)
{
  m_cLinkActinAndSites[clIDX][4]--;
}

void Actin::checkAndMoveClDown(std::vector<Actin> &actinVec, int subOne,
                               int subTwo, bool checkAll)
{
        /*
           General function that has input of two subids (subOne and subTwo)
           which are the subunit of the fila where the tethers are that will be checked
           Checks all tethers on the subs and determines if they need to move DOWN.
           subOne and subTwo are defaulted to -1 and -2 respectively (which are invalid)
           Therefore if no or just one argument is given the other is invalid and ignored
           checkAll is defaulted to true
         */

        int numCrossLinks = m_cLinkActinAndSites.size();

        for (int i = 0; i < numCrossLinks; ++i)
        {
                int j = m_cLinkActinAndSites[i][1]; // subid
                if (j != subOne && j!= subTwo && !checkAll)
                        continue;

                // j is the old sub


                double junk;
                int newClSub = findSubunit(m_cLinkActinAndSites[i][0], junk);

                if (j != newClSub)
                {

                        assert(newClSub == j-1);
                        // The tether subunit has moved down
                        setClSub(i, j-1, actinVec);
                }
        }
}

void Actin::checkAndMoveClUp(std::vector<Actin> &actinVec,
                             int subOne, int subTwo, bool checkAll)
{
        /*
           General function that has input of two subids (subOne and subTwo)
           which are the subunit of the fila where the tethers are that will be checked
           Checks all tethers on the subs and determines if they need to move UP.
           subOne and subTwo are defaulted to -1 and -2 respectively (which are invalid)
           Therefore if just one argument is given the other is invalid and ignored
         */

        int numCrossLinks = m_cLinkActinAndSites.size();

        for (int i = 0; i < numCrossLinks; ++i)
        {
                int j = m_cLinkActinAndSites[i][1]; // subid
                if (j != subOne && j!= subTwo && !checkAll)
                        continue;

                // j is the old sub


                double junk;
                int newClSub = findSubunit(m_cLinkActinAndSites[i][0], junk);

                if (j != newClSub)
                {
                        assert(newClSub == j+1);
                        // The tether subunit has moved up
                        setClSub(i, j+1, actinVec);
                }
        }
}

void Actin::checkAndMoveCl(std::vector<Actin> &actinVec)
{
        /*
           During severing moves the crosslink points up or down
         */

        int numCrossLinks = m_cLinkActinAndSites.size();

        for (int i = 0; i < numCrossLinks; ++i)
        {
                int j = m_cLinkActinAndSites[i][1]; // subid

                double junk;
                int newClSub = findSubunit(m_cLinkActinAndSites[i][0], junk);

                if (j != newClSub)
                {
                        setClSub(i, newClSub, actinVec);
                }
        }
}

void Actin::setClSub(int i, int subid, std::vector<Actin> &actinVec)
{
    m_cLinkActinAndSites[i][1] = subid;
    int otherFilaID = m_cLinkActinAndSites[i][2];
    int otherFilaCLIDX = m_cLinkActinAndSites[i][3];
    actinVec[otherFilaID].setClSubOther(otherFilaCLIDX, subid);

}


void Actin::moveClSub_down(std::vector<Actin> &actinVec)
{
        // Check for moving tethers
        /*
           This is for either barbed end poly or pointed end depoly, when the tether
           can move down a subunit.
           Called for when we are in initial straight regime
           Check ALL tethers along filament
         */
        checkAndMoveClDown(actinVec); // use default arguments

}

void Actin::moveClSub_up(std::vector<Actin> &actinVec)
{
        // Check for moving tethers
        /*
           This is for either barbed end depoly or pointed end poly, when the branch
           can move up a subunit.
           Called for when we are in initial straight regime
           Check ALL tethers along filament
         */

        checkAndMoveClUp(actinVec); // use default arguments

}

void Actin::moveClPoint_poly_barb(std::vector<Actin> &actinVec)
{
        // Check for moving tether points
        // After a polymerisation event, check any tether points on the final subunit
        // If any have "moved", they should be moved down to the penultimate subunit

        checkAndMoveClDown(actinVec, m_num_subunits-1,-2,false);

}

void Actin::moveClPoint_poly_point(std::vector<Actin> &actinVec)
{
        // Check for moving tether points
        // After a polymerisation event, check any tether points on the first subunit
        // If any have "moved", they should be moved up to the first subunit

        checkAndMoveClUp(actinVec, 0,-2,false);

}

void Actin::moveClPoint_depoly_barb(std::vector<Actin> &actinVec)
{
        // Check for moving tether points
        // After a depolymerisation event, check any tether points on the penultimate
        // subunit. If any have moved we know to move them up to the final sub

        checkAndMoveClUp(actinVec, m_num_subunits-2,-2,false);

}

void Actin::moveClPoint_depoly_point(std::vector<Actin> &actinVec)
{
        // Check for moving tether points
        // After a depolymerisation event, check any tether points on the second
        // subunit. If any have moved we know to move them down to the first sub
        checkAndMoveClDown(actinVec, 1,-2,false);

}

void Actin::moveClPoint_Add_barb(std::vector<Actin> &actinVec)
{
        // Check for moving tether points
        // After the addition of a barbed end point,
        // check any tether points on the penultimate and pen-penultimate subunit
        // If any have "moved", they should be moved up one

        checkAndMoveClUp(actinVec, m_num_subunits-3,m_num_subunits-2,false);

}

void Actin::moveClPoint_Add_point(std::vector<Actin> &actinVec)
{
        // Check for moving tether points
        // After the addition of a pointed end point,
        // check any tether points on the second and third subunit
        // If any have "moved", they should be moved down one
        checkAndMoveClDown(actinVec, 1,2,false);

}

void Actin::moveClPoint_Rem_barb(std::vector<Actin> &actinVec)
{
        // Check for moving tether points
        // After the removal of a barbed end point,
        // check any tether points on the penultimate sub
        // If any have "moved", they should be moved down one
        checkAndMoveClDown(actinVec, m_num_subunits-2,m_num_subunits-1,false);

}

void Actin::moveClPoint_Rem_point(std::vector<Actin> &actinVec)
{
        // Check for moving tether points
        // After the removal of a pointed end point,
        // check any tether points on the first sub
        // If any have "moved", they should be moved up one
        checkAndMoveClUp(actinVec, 0,1,false);
}

void Actin::shiftMonosCompatUp()
{
    // Upon pointed end polymerisation
    m_monosCompat.push_back(false);
    for (unsigned int i = m_monosCompat.size()-1; i > 0; --i)
    {
        m_monosCompat[i] = m_monosCompat[i-1];
    }
    m_monosCompat[0] = false;
}

void Actin::shiftMonosCompatDown()
{
    // Upon pointed end DEpolymerisation
    for (unsigned int i = 0; i < m_monosCompat.size()-1; ++i)
    {
        m_monosCompat[i] = m_monosCompat[i+1];
    }
    m_monosCompat.pop_back();
}

void Actin::shiftMonosBirthTimeUp()
{
    // Upon pointed end polymerisation
    m_monoBirthTime.push_back(0);
    for (unsigned int i = m_monoBirthTime.size()-1; i > 0; --i)
    {
        m_monoBirthTime[i] = m_monoBirthTime[i-1];
    }
    m_monoBirthTime[0] = 0;
}

void Actin::shiftMonosBirthTimeDown()
{
    // Upon pointed end DEpolymerisation
    for (unsigned int i = 0; i < m_monoBirthTime.size()-1; ++i)
    {
        m_monoBirthTime[i] = m_monoBirthTime[i+1];
    }
    m_monoBirthTime.pop_back();
}

void Actin::updateCrossLinkLoc(std::vector<Actin> &actinVec)
{
        // Function that uses the crosslink monomer to update the coordinates
        // of the crosslink point
        // point
        // Called after the filament has moved

        int numCrossLinks = m_cLinkActinAndSites.size();

        for (int i = 0; i < numCrossLinks; ++i)
        {
                int j = m_cLinkActinAndSites[i][1];
                double lenAlongSub;
                int subid = findSubunit(m_cLinkActinAndSites[i][0], lenAlongSub);
                //std::cout << "subid: " << subid << " j: " << j << std::endl;
                assert(subid == j);
                m_cLinkPoints[i][0] = m_points[j][0] + lenAlongSub*m_subunit_unitVecs[j][0];
                m_cLinkPoints[i][1] = m_points[j][1] + lenAlongSub*m_subunit_unitVecs[j][1];
                assert(std::isfinite(m_cLinkPoints[i][0]));
                assert(std::isfinite(m_cLinkPoints[i][1]));

                // Update other filament!
                int otherFilaID = m_cLinkActinAndSites[i][2];
                int clIDX = m_cLinkActinAndSites[i][3];
                actinVec[otherFilaID].setOtherCrossLinkPoint(clIDX, m_cLinkPoints[i][0], m_cLinkPoints[i][1]);

        }
}

void Actin::updateCrossLinkLocs(std::vector<Actin> &actinVec)
{
        // Function that uses the crosslink monomer to update the coordinates
        // of the crosslink point
        // Done for all filaments in same structure

        for (unsigned int i = 0; i < actinVec.size(); ++i)
        {
            if (m_structure_id == actinVec[i].getStructureID())
            {
                int numCrossLinks = actinVec[i].getNumCLinks();

                for (int j = 0; j < numCrossLinks; ++j)
                {
                        int k = actinVec[i].getCLinkActinAndSites()[j][1];
                        double lenAlongSub;
                        int subid = actinVec[i].findSubunit(actinVec[i].getCLinkActinAndSites()[j][0], lenAlongSub);
                        //std::cout << "subid: " << subid << " j: " << j << std::endl;
                        assert(subid == k);
                        double cLinkPx = actinVec[i].getPoint(k)[0] + lenAlongSub*actinVec[i].getUnitVec(k)[0];
                        double cLinkPy = actinVec[i].getPoint(k)[1] + lenAlongSub*actinVec[i].getUnitVec(k)[1];
                        actinVec[i].setCrossLinkPoint(j, cLinkPx, cLinkPy);
                        assert(std::isfinite(cLinkPx));
                        assert(std::isfinite(cLinkPy));

                        // Update other filament!
                        int otherFilaID = actinVec[i].getCLinkActinAndSites()[j][2];
                        int clIDX = actinVec[i].getCLinkActinAndSites()[j][3];
                        actinVec[otherFilaID].setOtherCrossLinkPoint(clIDX, cLinkPx, cLinkPy);

                }
            }
        }
}

void Actin::setCrossLinkPoint(int cLID, double pX, double pY)
{
    m_cLinkPoints[cLID][0] = pX;
    m_cLinkPoints[cLID][1] = pY;
}

void Actin::setOtherCrossLinkPoint(int cLID, double pX, double pY)
{
    m_cLinkPoints[cLID][2] = pX;
    m_cLinkPoints[cLID][3] = pY;
}

void Actin::recalcDistToLeadCL(int cLMono)
{
        // When creating a new crosslink the distance from the barbed end of the mother
        // to the crosslink point *may* need to be recalculated (it will if the new
        // cl is the leading branch)

        int barbedMono = m_num_monomers-1;

        if (barbedMono-cLMono < m_distToLeadCL)
        {
            m_distToLeadCL = barbedMono-cLMono;
        }


}


Eigen::VectorXd Actin::calcCrossLinkForces(int N, bool crosslinking,
                                            std::vector<Actin> &actinVec)
{
        // Forces that will be applied to each end point of our membraneFilament
        // Just uses Hookes Law F = -k(r_spring - r_point)
        // Returns 2 doubles per spring, (Fx, Fy)
        // These should then be added to the bending force

        Eigen::VectorXd crossLF (2*N);
        crossLF.setZero();
        int numCrossLinks = m_cLinkActinAndSites.size();
        if (!crosslinking || numCrossLinks == 0)
            return crossLF;

        for (int i = 0; i < numCrossLinks; ++i)
        {
                int j = m_cLinkActinAndSites[i][1];
                int otherID = m_cLinkActinAndSites[i][2];
                if (actinVec[otherID].getLength() < 2*s_segmentationLength)
                {
                    // ignore "one-sided" springs
                    continue;
                }

                double restLen = s_cLDist + m_radius + actinVec[otherID].getRadius();
                if (j == 0)
                {
                        double actualDist = sqrt((m_cLinkPoints[i][0] - m_cLinkPoints[i][2])*(m_cLinkPoints[i][0] - m_cLinkPoints[i][2]) + (m_cLinkPoints[i][1] - m_cLinkPoints[i][3])*(m_cLinkPoints[i][1] - m_cLinkPoints[i][3]));

                        crossLF(0) = crossLF(0) - s_cLStiff*(1-(restLen/actualDist))*(m_cLinkPoints[i][0] - m_cLinkPoints[i][2]);
                        crossLF(1) = crossLF(1) - s_cLStiff*(1-(restLen/actualDist))*(m_cLinkPoints[i][1] - m_cLinkPoints[i][3]);

                }
                else if (j == m_num_subunits-1)
                {
                        double actualDist = sqrt((m_cLinkPoints[i][0] - m_cLinkPoints[i][2])*(m_cLinkPoints[i][0] - m_cLinkPoints[i][2]) + (m_cLinkPoints[i][1] - m_cLinkPoints[i][3])*(m_cLinkPoints[i][1] - m_cLinkPoints[i][3]));

                        crossLF((2*N)-2) = crossLF((2*N)-2) - s_cLStiff*(1-(restLen/actualDist))*(m_cLinkPoints[i][0] - m_cLinkPoints[i][2]);
                        crossLF((2*N)-1) = crossLF((2*N)-1) - s_cLStiff*(1-(restLen/actualDist))*(m_cLinkPoints[i][1] - m_cLinkPoints[i][3]);
                }
                else
                {
                        double lenAlongSub;
                        int subid = findSubunit(m_cLinkActinAndSites[i][0], lenAlongSub);

                        assert (subid == j);
                        double a = lenAlongSub / m_prescribedSubLengths[j];
                        assert (a < 1);
                        // pointed end forces

                        double actualDist = sqrt((m_cLinkPoints[i][0] - m_cLinkPoints[i][2])*(m_cLinkPoints[i][0] - m_cLinkPoints[i][2]) + (m_cLinkPoints[i][1] - m_cLinkPoints[i][3])*(m_cLinkPoints[i][1] - m_cLinkPoints[i][3]));

                        crossLF((j-1)*2) = crossLF((j-1)*2) - (1-a)*(s_cLStiff*(1-(restLen/actualDist))*(m_cLinkPoints[i][0] - m_cLinkPoints[i][2]));
                        crossLF(((j-1)*2)+1) = crossLF(((j-1)*2)+1) - (1-a)*(s_cLStiff*(1-(restLen/actualDist))*(m_cLinkPoints[i][1] - m_cLinkPoints[i][3]));

                        // barbed end forces
                        crossLF((j)*2) = crossLF((j)*2) - a*(s_cLStiff*(1-(restLen/actualDist))*(m_cLinkPoints[i][0] - m_cLinkPoints[i][2]));
                        crossLF(((j)*2)+1) = crossLF(((j)*2)+1) - a*(s_cLStiff*(1-(restLen/actualDist))*(m_cLinkPoints[i][1] - m_cLinkPoints[i][3]));
                }
        }

        return crossLF;
}


std::array<int,2> Actin::findMonoIDLimitsOfSub(int subid)
{
  /*
   * Function that finds the monomer ids of the edges of the given subunit
   * inclusive
   */



    int lowLim, upLim;

    double startLen = 0;
    for (int i = 0; i < subid; ++i)
    {
        startLen += m_prescribedSubLengths[i];
    }
    lowLim = ((startLen-1E-14) / s_monomerLength) + 0.5;



    double endLen = 0;
    for (int i = 0; i <= subid; ++i)
    {
        endLen += m_prescribedSubLengths[i];
    }

    upLim = ((endLen-1E-14) / s_monomerLength) - 0.5;

    std::array<int,2> monoLimits;
    monoLimits[0] = lowLim;
    monoLimits[1] = upLim;

    return monoLimits;
}



// Static constants always defined at the bottom, here we have the radius of
// the actin filaments and the length of a single G-actin monomer.
const double Actin::s_stericRadius = 7.5E-9; // Excluded volume radius
const double Actin::s_trueRadius = 3.5E-9; // true radius (for diffusion)
const double Actin::s_monomerLength = 2.7E-9;

const int Actin::s_seedSize = 3; // Size of seed in monomers
const int Actin::s_branchSeedSize = 2; // Size of branch seed in monomers

int Actin::s_branchSpacing; // number of monomers between branches
int Actin::s_structure_idGenerator = 1; // start our ID generator with value 1
int Actin::s_subStructure_idGenerator = 1;
int Actin::s_crossStructure_idGenerator = 1;
// We are going to post-increment this when assigning

const double Actin::s_persistenceLength = 17E-6; // 17 microns

const double Actin::s_branchAngle = 1.2; // 70 deg in radians
double Actin::s_segmentationLength; // Have this here for now

// Spring constants have to be tuned carefully, they will depend on timestep and viscosity
double Actin::s_tetherStiff;
double Actin::s_branchStiff;
double Actin::s_torsionCoeff; // angular stiffness

double Actin::s_cLStiff; // crosslink spring stiffness
int Actin::s_cLSpacing;
double Actin::s_cLDist; // crosslink distance threshold and restLength

int Actin::s_maxSpacing;

int Actin::s_total_subunits = 0;
