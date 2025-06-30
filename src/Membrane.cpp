/*
 *  Membrane.cpp
 *
 *  C++ file containing the definition of the Membrane class.
 *  This is an inextensible membrane, made up of bead-rod components
 *  This contains member variable initialisations and non-trivial member
 *  functions.
 *  The declaration of this class (including the definitions of the trivial
 *  member functions) is found in the associated header file "Membrane.h"
 *
 *  James Bradford
 *  University of Sheffield
 *  Mar 2019
 */

#include "configHeader.h"
#include "Actin.h"
#include "ArpGrid.h"
#include "MembraneWall.h"
#include "Membrane.h"
#include "globals.h"
#include "geometry.h"
#include "polymerisation.h"
#include "Cortex.h"
#include <Eigen/SparseCore>
#include<Eigen/IterativeLinearSolvers>

extern RNG rng;
typedef gte::DCPQuery<double, gte::Segment<2, double>, gte::Segment<2, double>> SegDistQuery; // distance between two segments
typedef gte::TIQuery<double, gte::Segment<2, double>, gte::Segment<2, double>> Check_I; // check intersection between two lines
typedef gte::DCPQuery<double, gte::Vector<2,double>, gte::Segment<2, double>> DistQuery; // distance between point and segment

Membrane::Membrane()
    : m_exist { false }
{}

Membrane::Membrane(double ypos, double length, double temp, double bendMod,
                   bool fusable, bool active, double activeDist, bool initCortex,
                   double memTension)
    : m_exist { true },
    m_cellMem { true },
    m_id { 0 }, // Should be 0 for the original cell membrane
    m_ypos { ypos },
    m_length { length }, // for now centre at 0
    m_segLength { 2.5E-7 },
    //m_bendingMod { 20*g_Kb*temp }, // membrane + cortex?
    //m_bendingMod { 6*g_Kb*temp}, // pure membrane
    //m_bendingMod { 50*g_Kb*temp}, // membrane+cortex
    m_bendingMod { bendMod*g_Kb*temp },
    m_memTension { memTension },
    m_thickness { 5E-9 },
    m_activate { active },
    m_activeDist { activeDist },
    m_fusable { fusable },
    m_initCortex { initCortex }
{
    // Partition the membrane into points
    m_numPoints = length / m_segLength;

    double radius = length / (2*M_PI);
    m_points.push_back({+m_segLength/2, ypos-radius});
    m_points.push_back({0, 0});

    std::vector<double> tmpAngles = initCircle();
    m_points[1][0] = m_points[0][0] + m_segLength*cos(tmpAngles[0]);
    m_points[1][1] = m_points[0][1] + m_segLength*sin(tmpAngles[0]);

    double sumAngle = tmpAngles[0];
    for (int i = 2; i < m_numPoints; ++i)
    {
        sumAngle += tmpAngles[i-2];
        std::array<double,2> tmp;
        tmp[0] = m_points[i-1][0] + m_segLength*cos(sumAngle);
        tmp[1] = m_points[i-1][1] + m_segLength*sin(sumAngle);
        m_points.push_back(tmp);
    }

    for (int i = 0; i < m_numPoints; ++i)
    {
        std::array<double,2> tmp;
        m_subunit_unitVecs.push_back(tmp);
        m_actualSubLengths.push_back(m_segLength);
        m_prescribedSubLengths.push_back(m_segLength);

        m_stericCells.push_back(std::vector<int>());

        m_coupBranchRegIDs.push_back(std::vector<int>());
        m_coupNucRegIDs.push_back(std::vector<int>());
        m_coupCapRegIDs.push_back(std::vector<int>());
        m_coupAntiCapRegIDs.push_back(std::vector<int>());
        m_coupSevRegIDs.push_back(std::vector<int>());

    }

    calcSubunitVecs();


    m_zWidth =  m_length/M_PI;
    m_effPersLen = m_bendingMod*m_zWidth / (g_Kb*temp);
    m_stiffness = m_memTension*m_numPoints;

    setCOM();
    m_tetherStiff = 1E-7;

}

// Phagosome membrane - using average method
Membrane::Membrane(int numPoints, double temp, std::array<double,2> fusePoint,
                   Membrane &membrane, std::array<int,3> fusePointsToAv, int fuseType,
                   int ID)
    : m_exist { true },
    m_cellMem { false },
    m_id { ID },
    m_numPoints { numPoints },
    m_segLength { 2.5E-7 },
    m_bendingMod { 6*g_Kb*temp }, // bare membrane, no cortex
    m_thickness { 5E-9 },
    m_activate { false }, // phagosome doesn't stimulate actin polymerisation?
    m_activeDist { 0 },
    m_fusable { false },
    m_initCortex { false }
{
    // Partition the membrane into points

    int start, end;
    if (fuseType == 0 || fuseType == 1)
    {
        start = fusePointsToAv[0];
        end = fusePointsToAv[1];
    }
    else
    {
        start = fusePointsToAv[1];
        end = fusePointsToAv[2];
    }

    // Different order
    int newStart = (start+end)/2;
    for (int i = newStart; i < end; ++i)
    {
        m_points.push_back(membrane.getPoints()[i]);
    }

    m_points.push_back({fusePoint[0], fusePoint[1]});

    for (int i = start+1; i < newStart; ++i)
    {
        m_points.push_back(membrane.getPoints()[i]);

    }

    for (int i = 0; i < m_numPoints; ++i)
    {
        std::array<double,2> tmp;
        m_subunit_unitVecs.push_back(tmp); // will update
        m_prescribedSubLengths.push_back(m_segLength); // will change
        m_actualSubLengths.push_back(m_segLength); // will update
        m_stericCells.push_back(std::vector<int>()); // will update in fusion function

        m_coupBranchRegIDs.push_back(std::vector<int>());
        m_coupNucRegIDs.push_back(std::vector<int>());
        m_coupCapRegIDs.push_back(std::vector<int>());
        m_coupAntiCapRegIDs.push_back(std::vector<int>());
        m_coupSevRegIDs.push_back(std::vector<int>());
    }
    calcSubLengths();

    calcSubunitVecs();


    // Correct m_prescribedSubLengths, want them all equal
    m_length = 0;
    for (int i = 0; i < m_numPoints; ++i)
    {
        m_length += m_actualSubLengths[i];
    }

    m_segLength = m_length/m_numPoints;

    for (int i = 0; i < m_numPoints; ++i)
    {
        m_prescribedSubLengths[i] = m_segLength;
    }

    m_zWidth =  (m_length/M_PI); // w in calculations
    m_effPersLen = m_bendingMod*m_zWidth / (g_Kb*temp);
    m_memTension = membrane.getMemTension();
    m_stiffness = m_memTension*m_numPoints;
}

Membrane::~Membrane(){
    //std::cout << "Membrane correctly deleted" << std::endl;
}

std::vector<double> Membrane::initCircle()
{
    /*
    Function that initialise the membrane filament in a (near)perfect circle
    */
    double angleDiff = (2*M_PI) / m_numPoints;

    std::vector<double> tmpAngles;
    for (int i = 0; i < m_numPoints-2; ++i)
    {
        tmpAngles.push_back(angleDiff);

    }
    return tmpAngles;
}

void Membrane::calcSubunitVecs()
{
    for (int i = 0; i < m_numPoints; ++i)
    {
        int p1 = i + 1; // next point
        if (i == m_numPoints -1)
        {
            p1 = 0;
        }

        double dx = m_points[p1][0] - m_points[i][0];
        double dy = m_points[p1][1] - m_points[i][1];
        double mag = sqrt(dx*dx + dy*dy);
        m_subunit_unitVecs[i][0] = dx/mag;
        m_subunit_unitVecs[i][1] = dy/mag;
    }
}

void Membrane::calcSubLengths()
{
    for (int i = 0; i < m_numPoints; ++i)
    {
        int p1 = i + 1; // next point
        if (i == m_numPoints -1)
        {
            p1 = 0;
        }
        std::array<double,2> point1 { m_points[i][0], m_points[i][1] };
        std::array<double,2> point2 { m_points[p1][0], m_points[p1][1] };
        m_actualSubLengths[i] = sqrt(distanceBetPoints2DSQR(point1, point2));
    }
}

Eigen::VectorXd Membrane::calcBendingForces(double temp)
{

    double g = ((m_bendingMod * m_zWidth) / m_segLength);
    int N = m_numPoints; // number of points

    Eigen::VectorXd bendingForces (2*N);
    for (int n = 0; n < N; ++n)
    {
        // Define nearby points
        int m2 = n-2;  // minus 2
        int m1 = n-1; // minus 1
        int p1 = n+1; // plus 1

        if (n == 0)
        {
            m2 = N-2;
            m1 = N-1;
        }
        else if (n == 1)
        {
            m2 = N-1;
        }
        else if (n == N-1)
        {
            p1 = 0;
        }
        assert(m2 >= 0);
        assert(m1 >= 0);
        assert(p1 >= 0);


        bendingForces[2*n] = -g*(m_subunit_unitVecs[m1][0] - m_subunit_unitVecs[m2][0])/m_actualSubLengths[m1]
                               + g*(m_subunit_unitVecs[n][0] - m_subunit_unitVecs[m1][0])*(1/m_actualSubLengths[m1] + 1/m_actualSubLengths[n])
                               - g*(m_subunit_unitVecs[p1][0] - m_subunit_unitVecs[n][0])/m_actualSubLengths[n];




        bendingForces[(2*n)+1] = -g*(m_subunit_unitVecs[m1][1] - m_subunit_unitVecs[m2][1])/m_actualSubLengths[m1]
                               + g*(m_subunit_unitVecs[n][1] - m_subunit_unitVecs[m1][1])*(1/m_actualSubLengths[m1] + 1/m_actualSubLengths[n])
                               - g*(m_subunit_unitVecs[p1][1] - m_subunit_unitVecs[n][1])/m_actualSubLengths[n];



    }

    for (int i = 0; i < bendingForces.size(); ++i)
    {

        if (!std::isfinite(bendingForces[i]))
        {
            std::cout << m_length << "; " << m_effPersLen << "; " << std::endl;
            for (int j = 0; j < m_numPoints; ++j)
            {
                std::cout << m_points[j][0] << ", " << m_points[j][1] << std::endl;
            }

            for (int j = 0; j < m_numPoints; ++j)
            {
                std::cout << m_actualSubLengths[j] << std::endl;
            }
        }
        assert(std::isfinite(bendingForces[i]));
    }

    return bendingForces;
}

Eigen::VectorXd Membrane::calcBSForces()
{

    Eigen::VectorXd bsForces (2*m_numPoints);
    for (int n = 0; n < m_numPoints; ++n)
    {
        // Define nearby points
        int m1 = n-1; // minus 1

        if (n == 0)
        {
            m1 = m_numPoints-1;
        }

        assert(m1 >= 0);

        // Below mag actual have direction, and it differs for each spring
        double mag1 = m_stiffness*(m_prescribedSubLengths[m1] - m_actualSubLengths[m1]);
        double mag2 = -m_stiffness*(m_prescribedSubLengths[n] - m_actualSubLengths[n]);


        bsForces[2*n] = m_subunit_unitVecs[n][0]*mag2 + m_subunit_unitVecs[m1][0]*mag1;

        bsForces[(2*n)+1] = m_subunit_unitVecs[n][1]*mag2 + m_subunit_unitVecs[m1][1]*mag1;


    }

    for (int i = 0; i < bsForces.size(); ++i)
    {

        if (!std::isfinite(bsForces[i]))
        {
            std::cout << m_length << "; " << m_effPersLen << "; " << std::endl;
            for (int j = 0; j < m_numPoints; ++j)
            {
                std::cout << m_points[j][0] << ", " << m_points[j][1] << std::endl;
            }

            for (int j = 0; j < m_numPoints; ++j)
            {
                std::cout << m_actualSubLengths[j] << std::endl;
            }
        }
        assert(std::isfinite(bsForces[i]));
    }

    return bsForces;
}

Eigen::MatrixXd Membrane::buildDMatrix(double temp, double viscosity)
{
    int N = m_numPoints;

    double KT = (g_Kb*temp);
    Eigen::Matrix2d I;
    I (0,0) = 1;
    I (0,1) = 0;
    I (1,0) = 0;
    I (1,1) = 1;
    Eigen::MatrixXd transDiffMat (N*2+2,N*2+2);
    transDiffMat = Eigen::MatrixXd::Zero(N*2+2,N*2+2);

    double Yama_corr_X = 0.044;
    double Yama_corr_Y = 1.111;
    double drag_para = (2*M_PI*viscosity*m_segLength)/(log(m_segLength/(m_thickness))+Yama_corr_X);
    double drag_perp = (4*M_PI*viscosity*m_segLength)/(log(m_segLength/(m_thickness))+Yama_corr_Y);

    // j and k are the indices of the 2x2 subblocks
    // they represent the beads
    for (int j = 0; j < N; ++j)
    {
        // Calculate the local tangent vector
        int m1 = j-1;
        if (j == 0)
        {
            m1 = N - 1;
        }
        Eigen::Vector2d u_tilde;

        Eigen::Vector2d u_j;
        u_j << m_subunit_unitVecs[j][0], m_subunit_unitVecs[j][1];
        Eigen::Vector2d u_jm1;
        u_jm1 << m_subunit_unitVecs[m1][0], m_subunit_unitVecs[m1][1];

        u_tilde = (u_j + u_jm1) / (u_j + u_jm1).norm();


        Eigen::Matrix2d mobKT = KT*((1/drag_para)*(u_tilde*u_tilde.transpose()) + (1/drag_perp)*(I-(u_tilde*u_tilde.transpose()))); // equation 25
        transDiffMat.block<2,2>(j*2,j*2) = mobKT;
    }

    return transDiffMat;
}

Eigen::VectorXd Membrane::getRand_dR(double dt, Eigen::MatrixXd D)
{
    // xi 2N vector

    int N = m_numPoints;
    Eigen::VectorXd rand_dR (2*N);
    // For projection method
    //Eigen
    for (int i = 0; i < 2*N; ++i)
    {
        std::normal_distribution<double> noise(0,sqrt(2*D(i,i)*dt)); // it is a dR
        rand_dR (i) = noise(rng.m_mersenne);
    }

    return rand_dR;
}

Eigen::MatrixXd Membrane::getBMatrix()
{
    int N = m_numPoints;
    Eigen::MatrixXd B (N, 2*N);
    // fill with zeroes
    B = Eigen::MatrixXd::Zero(N,2*N);
    double mag = 0;
    for (int row = 0; row < N; ++row)
    {

        int startcol = 2*row;
        int nextPoint = row + 1;
        if (row == N-1)
        {
            // Constraint connecting end to start point
            nextPoint = 0;
        }

        //First pair
        B (row, startcol) = m_points[row][0] - m_points[nextPoint][0];
        B (row, startcol+1) = m_points[row][1] - m_points[nextPoint][1];

        mag = sqrt( B(row,startcol)*B(row,startcol) + B(row,startcol+1)*B(row,startcol+1) );
        B (row, startcol) = B (row, startcol) / mag;
        B (row, startcol+1) = B (row, startcol+1) / mag;



        // Second pair
        if (row != N-1)
        {
          B (row, startcol+2) = m_points[nextPoint][0] - m_points[row][0];
          B (row, startcol+3) = m_points[nextPoint][1] - m_points[row][1];

          mag = sqrt( B(row,startcol+2)*B(row,startcol+2) + B(row,startcol+3)*B(row,startcol+3) );
          B (row, startcol+2) = B (row, startcol+2) / mag;
          B (row, startcol+3) = B (row, startcol+3) / mag;
        }
        else
        {
          // Periodicity
          B (row, 0) = m_points[nextPoint][0] - m_points[row][0];
          B (row, 1) = m_points[nextPoint][1] - m_points[row][1];

          mag = sqrt( B(row,0)*B(row,0) + B(row,1)*B(row,1) );
          B (row, 0) = B (row, 0) / mag;
          B (row, 1) = B (row, 1) / mag;
        }

    }
    return B;
}



void Membrane::fluctuate(double temp, double viscosity, double dt,
                             std::vector<Actin> &actinvec,
                             std::vector<ExcZone> &excZones,
                             const bool steric, unsigned int &timesFused,
                             StericGrid &stericGrid, double currTime,
                             std::vector<ProteinRegion> &branchRegions,
                             bool &branching, double arpConc,
                             std::vector<ProteinRegion> &capRegions,
                             std::vector<ProteinRegion> &antiCapRegions,
                             std::vector<ProteinRegion> &nucRegions,
                             std::vector<ProteinRegion> &sevRegions,
                             std::vector<Membrane> &membranes,
                             int &nActin, const bool tether,
                             bool arpPool, ArpGrid &arpGrid,
                             std::vector<MembraneWall> &memWalls,
                             std::ofstream &fusionTime, const bool bDynamics,
                             Cortex &cortex, bool COMadjust,
                             const bool memSprings,
                             const bool activeBranch,
                             const double activeBranchHeight,
                             const bool activeNuc, const double activeNucHeight,
                             const bool activeCap, const double activeCapHeight,
                             const bool activeAntiCap,
                             const double activeAntiCapHeight,
                             const bool activeSever,
                             const double activeSeverHeight)
{
    /*
    Function that is to be called every timestep to fluctuate
    the Membrane

    Same algorithm used for dynamics of our actin filaments

    */


    // N is the number of points (midpoints of the 'rods')
    int N = m_numPoints;

    Eigen::MatrixXd B (N, 2*N);
    Eigen::VectorXd F_spr (2*N);
    if (memSprings)
    {
        F_spr = calcBSForces();
    }
    else
    {
        B = getBMatrix();
    }

    Eigen::VectorXd F_bend (2*N);
    F_bend = calcBendingForces(temp);


    Eigen::VectorXd F_tet (2*N);
    F_tet = calcTetherForces(cortex);

    Eigen::MatrixXd D (2*N, 2*N);

    D = buildDMatrix(temp, viscosity);

    Eigen::VectorXd xi (2*N);

    xi = getRand_dR(dt, D);

    Eigen::VectorXd r_old (2*N);
    for (int j = 0; j < (2*N)-1; j+=2)
    {
        // input our points

        int k = j / 2;
        r_old (j) = m_points[k][0]; // x
        r_old (j+1) = m_points[k][1]; // y
    }

    Eigen::VectorXd r_new (2*N);
    if (memSprings)
    {
        r_new = r_old + (dt/(g_Kb*temp))*D*(F_bend+F_tet+F_spr) + xi;
    }
    else
    {
        //Eigen::MatrixXd I (2*N, 2*N);
        //I = Eigen::MatrixXd::Identity(2*N,2*N);
        /*
        Eigen::MatrixXd T (2*N, N);
        // Following line is cpu intensive due to inverse, Wang and Gao suggest
        // approximating T

        Eigen::MatrixXd B_T (2*N, N);
        B_T = B.transpose();

        T = B_T * (B*B_T).inverse();

        // Here is the algorithm
        // d is the vector contained the prescribed bond lengths (between points)


        Eigen::VectorXd d = Eigen::VectorXd::Map(&m_prescribedSubLengths[0], m_prescribedSubLengths.size());

        // Modify to change forces
        //r_new = (I - T*B)*(r_old + (dt/(g_Kb*temp))*D*(F_bend + xi)) + T*d;


        //r_new = (I - T*B)*(r_old + (dt/(g_Kb*temp))*D*(F_bend+F_tet) + xi) + T*d;
        // Avoid matrix-matrix multiplication
        Eigen::VectorXd vector (2*N);
        vector = (r_old + (dt/(g_Kb*temp))*D*(F_bend+F_tet) + xi);
        //r_new = vector - T*(B*vector) + T*d;
        r_new = vector - T*((B*vector) - d);
        */
        Eigen::MatrixXd B_T (2*N, N);
        B_T = B.transpose();

        Eigen::VectorXd d = Eigen::VectorXd::Map(&m_prescribedSubLengths[0], m_prescribedSubLengths.size());
        Eigen::VectorXd vector (2*N);
        vector = (r_old + (dt/(g_Kb*temp))*D*(F_bend+F_tet) + xi);

        Eigen::VectorXd BtimesVec (N);
        BtimesVec = B*vector;

        r_new = vector - B_T*((B*B_T).llt().solve(BtimesVec-d));
    }


    for (int j = 0; j < N; ++j)
    {
        // j is the correct index for r_new
        m_points[j][0] = r_new (j*2);
        m_points[j][1] = r_new (j*2 + 1);


        assert(std::isfinite(m_points[j][0]));
        assert(std::isfinite(m_points[j][1]));


        calcSubLengths();
        calcSubunitVecs();

        if (steric)
        {
            // Update steric grid before check
            if (j == 0)
            {
                stericGrid.resetMemandUpdate(*this, 0);
                stericGrid.resetMemandUpdate(*this, m_numPoints-1);
            }
            else
            {
                stericGrid.resetMemandUpdate(*this, j-1);
                stericGrid.resetMemandUpdate(*this, j);
            }
        }

        if (steric && (checkExVolSUBGrid(actinvec, j, stericGrid) || ExcZone::s_checkExVolMemFila(*this, excZones, j) || checkExVolSUBGridCortex(cortex, j, stericGrid) ))
        {
            // This move has violated excvol and therefore we reject it
            m_points[j][0] = r_old (j*2);
            m_points[j][1] = r_old (j*2 + 1);
            calcSubLengths();
            calcSubunitVecs();
            // Update steric grid back
            if (j == 0)
            {
                stericGrid.resetMemandUpdate(*this, 0);
                stericGrid.resetMemandUpdate(*this, m_numPoints-1);
            }
            else
            {
                stericGrid.resetMemandUpdate(*this, j-1);
                stericGrid.resetMemandUpdate(*this, j);
            }
        }
        else
        {
            if (m_activate)
            {
                // 1. Calc distance from each sub to target(s)

                //2. Add/remove accordingly
                // What to add remove? Nuc? Branch? ACap?


                if (calcDistToTarsForActive(excZones, j))
                {
                    // membrane should be active

                    // if not already a region there, make one
                    if (activeBranch && m_coupBranchRegIDs[j].size() == 0)
                    {
                        std::array<double,2> uVec = { m_subunit_unitVecs[j][0], m_subunit_unitVecs[j][1] };
                        std::array<double,2> pointBL = { m_points[j][0], m_points[j][1] };
                        branchRegions.push_back((ProteinRegion(uVec, m_prescribedSubLengths[j], activeBranchHeight, pointBL, arpConc, 0, 0, 0, 1, j)));

                        branchRegions[branchRegions.size()-1].stickToMem(*this);
                        addToBranchReg(j, branchRegions.size()-1);
                    }

                    if (activeNuc && m_coupNucRegIDs[j].size() == 0)
                    {
                        std::array<double,2> uVec = { m_subunit_unitVecs[j][0], m_subunit_unitVecs[j][1] };
                        std::array<double,2> pointBL = { m_points[j][0], m_points[j][1] };
                        nucRegions.push_back((ProteinRegion(uVec, m_prescribedSubLengths[j], activeNucHeight, pointBL, 0, 0, 1.0, 0, 1, j)));

                        nucRegions[nucRegions.size()-1].stickToMem(*this);
                        addToNucReg(j, nucRegions.size()-1);
                    }

                    if (activeCap && m_coupCapRegIDs[j].size() == 0)
                    {
                        std::array<double,2> uVec = { m_subunit_unitVecs[j][0], m_subunit_unitVecs[j][1] };
                        std::array<double,2> pointBL = { m_points[j][0], m_points[j][1] };
                        capRegions.push_back((ProteinRegion(uVec, m_prescribedSubLengths[j], activeCapHeight, pointBL, 0, 1.0, 0, 0, 1, j)));

                        capRegions[capRegions.size()-1].stickToMem(*this);
                        addToCapReg(j, capRegions.size()-1);
                    }

                    if (activeAntiCap && m_coupAntiCapRegIDs[j].size() == 0)
                    {
                        std::array<double,2> uVec = { m_subunit_unitVecs[j][0], m_subunit_unitVecs[j][1] };
                        std::array<double,2> pointBL = { m_points[j][0], m_points[j][1] };
                        antiCapRegions.push_back((ProteinRegion(uVec, m_prescribedSubLengths[j], activeAntiCapHeight, pointBL, 0, 0, 0, 0, 1, j)));

                        antiCapRegions[antiCapRegions.size()-1].stickToMem(*this);
                        addToAntiCapReg(j, antiCapRegions.size()-1);
                    }

                    if (activeSever && m_coupSevRegIDs[j].size() == 0)
                    {
                        std::array<double,2> uVec = { m_subunit_unitVecs[j][0], m_subunit_unitVecs[j][1] };
                        std::array<double,2> pointBL = { m_points[j][0], m_points[j][1] };
                        sevRegions.push_back((ProteinRegion(uVec, m_prescribedSubLengths[j], activeSeverHeight, pointBL, 0, 0, 0, 1.0, 1, j)));

                        sevRegions[sevRegions.size()-1].stickToMem(*this);
                        addToSevReg(j, sevRegions.size()-1);
                    }

                    // New idea - turn on neighbours too
                    int plus1 = j+1;
                    int minus1 = j-1;
                    if (j == 0)
                    {
                        minus1 = m_numPoints-1;
                    }
                    else if (j == m_numPoints - 1)
                    {
                        plus1 = 0;
                    }

                    // if not already a region there, make one
                    if (activeBranch && m_coupBranchRegIDs[minus1].size() == 0)
                    {
                        std::array<double,2> uVec = { m_subunit_unitVecs[minus1][0], m_subunit_unitVecs[minus1][1] };
                        std::array<double,2> pointBL = { m_points[minus1][0], m_points[minus1][1] };
                        branchRegions.push_back((ProteinRegion(uVec, m_prescribedSubLengths[minus1], activeBranchHeight, pointBL, arpConc, 0, 0, 0, 1, minus1)));

                        branchRegions[branchRegions.size()-1].stickToMem(*this);
                        addToBranchReg(minus1, branchRegions.size()-1);
                    }

                    if (activeNuc && m_coupNucRegIDs[minus1].size() == 0)
                    {
                        std::array<double,2> uVec = { m_subunit_unitVecs[minus1][0], m_subunit_unitVecs[minus1][1] };
                        std::array<double,2> pointBL = { m_points[minus1][0], m_points[minus1][1] };
                        nucRegions.push_back((ProteinRegion(uVec, m_prescribedSubLengths[minus1], activeNucHeight, pointBL, 0, 0, 1.0, 0, 1, minus1)));

                        nucRegions[nucRegions.size()-1].stickToMem(*this);
                        addToNucReg(minus1, nucRegions.size()-1);
                    }

                    if (activeCap && m_coupCapRegIDs[minus1].size() == 0)
                    {
                        std::array<double,2> uVec = { m_subunit_unitVecs[minus1][0], m_subunit_unitVecs[minus1][1] };
                        std::array<double,2> pointBL = { m_points[minus1][0], m_points[minus1][1] };
                        capRegions.push_back((ProteinRegion(uVec, m_prescribedSubLengths[minus1], activeCapHeight, pointBL, 0, 1.0, 0, 0, 1, minus1)));

                        capRegions[capRegions.size()-1].stickToMem(*this);
                        addToCapReg(minus1, capRegions.size()-1);
                    }

                    if (activeAntiCap && m_coupAntiCapRegIDs[minus1].size() == 0)
                    {
                        std::array<double,2> uVec = { m_subunit_unitVecs[minus1][0], m_subunit_unitVecs[minus1][1] };
                        std::array<double,2> pointBL = { m_points[minus1][0], m_points[minus1][1] };
                        antiCapRegions.push_back((ProteinRegion(uVec, m_prescribedSubLengths[minus1], activeAntiCapHeight, pointBL, 0, 0, 0, 0, 1, minus1)));

                        antiCapRegions[antiCapRegions.size()-1].stickToMem(*this);
                        addToAntiCapReg(minus1, antiCapRegions.size()-1);
                    }

                    if (activeSever && m_coupSevRegIDs[minus1].size() == 0)
                    {
                        std::array<double,2> uVec = { m_subunit_unitVecs[minus1][0], m_subunit_unitVecs[minus1][1] };
                        std::array<double,2> pointBL = { m_points[minus1][0], m_points[minus1][1] };
                        sevRegions.push_back((ProteinRegion(uVec, m_prescribedSubLengths[minus1], activeSeverHeight, pointBL, 0, 0, 0, 1.0, 1, minus1)));

                        sevRegions[sevRegions.size()-1].stickToMem(*this);
                        addToSevReg(minus1, sevRegions.size()-1);
                    }

                    // if not already a region there, make one
                    if (activeBranch && m_coupBranchRegIDs[plus1].size() == 0)
                    {
                        std::array<double,2> uVec = { m_subunit_unitVecs[plus1][0], m_subunit_unitVecs[plus1][1] };
                        std::array<double,2> pointBL = { m_points[plus1][0], m_points[plus1][1] };
                        branchRegions.push_back((ProteinRegion(uVec, m_prescribedSubLengths[plus1], activeBranchHeight, pointBL, arpConc, 0, 0, 0, 1, plus1)));

                        branchRegions[branchRegions.size()-1].stickToMem(*this);
                        addToBranchReg(plus1, branchRegions.size()-1);
                    }

                    if (activeNuc && m_coupNucRegIDs[plus1].size() == 0)
                    {
                        std::array<double,2> uVec = { m_subunit_unitVecs[plus1][0], m_subunit_unitVecs[plus1][1] };
                        std::array<double,2> pointBL = { m_points[plus1][0], m_points[plus1][1] };
                        nucRegions.push_back((ProteinRegion(uVec, m_prescribedSubLengths[plus1], activeNucHeight, pointBL, 0, 0, 1.0, 0, 1, plus1)));

                        nucRegions[nucRegions.size()-1].stickToMem(*this);
                        addToNucReg(plus1, nucRegions.size()-1);
                    }

                    if (activeCap && m_coupCapRegIDs[plus1].size() == 0)
                    {
                        std::array<double,2> uVec = { m_subunit_unitVecs[plus1][0], m_subunit_unitVecs[plus1][1] };
                        std::array<double,2> pointBL = { m_points[plus1][0], m_points[plus1][1] };
                        capRegions.push_back((ProteinRegion(uVec, m_prescribedSubLengths[plus1], activeCapHeight, pointBL, 0, 1.0, 0, 0, 1, plus1)));

                        capRegions[capRegions.size()-1].stickToMem(*this);
                        addToCapReg(plus1, capRegions.size()-1);
                    }

                    if (activeAntiCap && m_coupAntiCapRegIDs[plus1].size() == 0)
                    {
                        std::array<double,2> uVec = { m_subunit_unitVecs[plus1][0], m_subunit_unitVecs[plus1][1] };
                        std::array<double,2> pointBL = { m_points[plus1][0], m_points[plus1][1] };
                        antiCapRegions.push_back((ProteinRegion(uVec, m_prescribedSubLengths[plus1], activeAntiCapHeight, pointBL, 0, 0, 0, 0, 1, plus1)));

                        antiCapRegions[antiCapRegions.size()-1].stickToMem(*this);
                        addToAntiCapReg(plus1, antiCapRegions.size()-1);
                    }

                    if (activeSever && m_coupSevRegIDs[plus1].size() == 0)
                    {
                        std::array<double,2> uVec = { m_subunit_unitVecs[plus1][0], m_subunit_unitVecs[plus1][1] };
                        std::array<double,2> pointBL = { m_points[plus1][0], m_points[plus1][1] };
                        sevRegions.push_back((ProteinRegion(uVec, m_prescribedSubLengths[plus1], activeSeverHeight, pointBL, 0, 0, 0, 1.0, 1, plus1)));

                        sevRegions[sevRegions.size()-1].stickToMem(*this);
                        addToSevReg(plus1, sevRegions.size()-1);
                    }

                    // Delete any tethers on sub j
                    int numTethers = m_tetherRelDists.size();
                    int k = 0;
                    while (k < numTethers)
                    {
                        if (m_tetherSubIDs[k] == j || m_tetherSubIDs[k] == plus1 || m_tetherSubIDs[k] == minus1)
                        {
                            // Delete the tether
                            numTethers -= 1;
                            m_tetherSubIDs.erase(m_tetherSubIDs.begin()+k);
                            m_tetherRelDists.erase(m_tetherRelDists.begin()+k);
                            // Delete on cortex
                            cortex.eraseTether(k);
                        }
                        else
                        {
                            ++k;
                        }
                    }

                }

                // Check prev sub
                int prevSub = j-1;
                if (prevSub == -1)
                {
                    prevSub = m_numPoints - 1;
                }

                if (calcDistToTarsForActive(excZones, prevSub))
                {
                    // membrane should be active

                    // if not already a region there, make one
                    if (activeBranch && m_coupBranchRegIDs[prevSub].size() == 0)
                    {
                        std::array<double,2> uVec = { m_subunit_unitVecs[prevSub][0], m_subunit_unitVecs[prevSub][1] };
                        std::array<double,2> pointBL = { m_points[prevSub][0], m_points[prevSub][1] };
                        branchRegions.push_back((ProteinRegion(uVec, m_prescribedSubLengths[prevSub], activeBranchHeight, pointBL, arpConc, 0, 0, 0, 1, prevSub)));

                        branchRegions[branchRegions.size()-1].stickToMem(*this);
                        addToBranchReg(prevSub, branchRegions.size()-1);
                    }

                    if (activeNuc && m_coupNucRegIDs[prevSub].size() == 0)
                    {
                        std::array<double,2> uVec = { m_subunit_unitVecs[prevSub][0], m_subunit_unitVecs[prevSub][1] };
                        std::array<double,2> pointBL = { m_points[prevSub][0], m_points[prevSub][1] };
                        nucRegions.push_back((ProteinRegion(uVec, m_prescribedSubLengths[prevSub], activeNucHeight, pointBL, 0, 0, 1.0, 0, 1, prevSub)));

                        nucRegions[nucRegions.size()-1].stickToMem(*this);
                        addToNucReg(prevSub, nucRegions.size()-1);
                    }

                    if (activeCap && m_coupCapRegIDs[prevSub].size() == 0)
                    {
                        std::array<double,2> uVec = { m_subunit_unitVecs[prevSub][0], m_subunit_unitVecs[prevSub][1] };
                        std::array<double,2> pointBL = { m_points[prevSub][0], m_points[prevSub][1] };
                        capRegions.push_back((ProteinRegion(uVec, m_prescribedSubLengths[prevSub], activeCapHeight, pointBL, 0, 1.0, 0, 0, 1, prevSub)));

                        capRegions[capRegions.size()-1].stickToMem(*this);
                        addToCapReg(prevSub, capRegions.size()-1);
                    }

                    if (activeAntiCap && m_coupAntiCapRegIDs[prevSub].size() == 0)
                    {
                        std::array<double,2> uVec = { m_subunit_unitVecs[prevSub][0], m_subunit_unitVecs[prevSub][1] };
                        std::array<double,2> pointBL = { m_points[prevSub][0], m_points[prevSub][1] };
                        antiCapRegions.push_back((ProteinRegion(uVec, m_prescribedSubLengths[prevSub], activeAntiCapHeight, pointBL, 0, 0, 0, 0, 1, prevSub)));

                        antiCapRegions[antiCapRegions.size()-1].stickToMem(*this);
                        addToAntiCapReg(prevSub, antiCapRegions.size()-1);
                    }

                    if (activeSever && m_coupSevRegIDs[prevSub].size() == 0)
                    {
                        std::array<double,2> uVec = { m_subunit_unitVecs[prevSub][0], m_subunit_unitVecs[prevSub][1] };
                        std::array<double,2> pointBL = { m_points[prevSub][0], m_points[prevSub][1] };
                        sevRegions.push_back((ProteinRegion(uVec, m_prescribedSubLengths[prevSub], activeSeverHeight, pointBL, 0, 0, 0, 1.0, 1, prevSub)));

                        sevRegions[sevRegions.size()-1].stickToMem(*this);
                        addToSevReg(prevSub, sevRegions.size()-1);
                    }

                    // New idea - turn on neighbours too
                    int plus1 = prevSub + 1;
                    int minus1 = prevSub - 1;
                    if (prevSub == 0)
                    {
                        minus1 = m_numPoints-1;
                    }
                    else if (prevSub == m_numPoints - 1)
                    {
                        plus1 = 0;
                    }

                    // if not already a region there, make one
                    if (activeBranch && m_coupBranchRegIDs[minus1].size() == 0)
                    {
                        std::array<double,2> uVec = { m_subunit_unitVecs[minus1][0], m_subunit_unitVecs[minus1][1] };
                        std::array<double,2> pointBL = { m_points[minus1][0], m_points[minus1][1] };
                        branchRegions.push_back((ProteinRegion(uVec, m_prescribedSubLengths[minus1], activeBranchHeight, pointBL, arpConc, 0, 0, 0, 1, minus1)));

                        branchRegions[branchRegions.size()-1].stickToMem(*this);
                        addToBranchReg(minus1, branchRegions.size()-1);
                    }

                    if (activeNuc && m_coupNucRegIDs[minus1].size() == 0)
                    {
                        std::array<double,2> uVec = { m_subunit_unitVecs[minus1][0], m_subunit_unitVecs[minus1][1] };
                        std::array<double,2> pointBL = { m_points[minus1][0], m_points[minus1][1] };
                        nucRegions.push_back((ProteinRegion(uVec, m_prescribedSubLengths[minus1], activeNucHeight, pointBL, 0, 0, 1.0, 0, 1, minus1)));

                        nucRegions[nucRegions.size()-1].stickToMem(*this);
                        addToNucReg(minus1, nucRegions.size()-1);
                    }

                    if (activeCap && m_coupCapRegIDs[minus1].size() == 0)
                    {
                        std::array<double,2> uVec = { m_subunit_unitVecs[minus1][0], m_subunit_unitVecs[minus1][1] };
                        std::array<double,2> pointBL = { m_points[minus1][0], m_points[minus1][1] };
                        capRegions.push_back((ProteinRegion(uVec, m_prescribedSubLengths[minus1], activeCapHeight, pointBL, 0, 1.0, 0, 0, 1, minus1)));

                        capRegions[capRegions.size()-1].stickToMem(*this);
                        addToCapReg(minus1, capRegions.size()-1);
                    }

                    if (activeAntiCap && m_coupAntiCapRegIDs[minus1].size() == 0)
                    {
                        std::array<double,2> uVec = { m_subunit_unitVecs[minus1][0], m_subunit_unitVecs[minus1][1] };
                        std::array<double,2> pointBL = { m_points[minus1][0], m_points[minus1][1] };
                        antiCapRegions.push_back((ProteinRegion(uVec, m_prescribedSubLengths[minus1], activeAntiCapHeight, pointBL, 0, 0, 0, 0, 1, minus1)));

                        antiCapRegions[antiCapRegions.size()-1].stickToMem(*this);
                        addToAntiCapReg(minus1, antiCapRegions.size()-1);
                    }

                    if (activeSever && m_coupSevRegIDs[minus1].size() == 0)
                    {
                        std::array<double,2> uVec = { m_subunit_unitVecs[minus1][0], m_subunit_unitVecs[minus1][1] };
                        std::array<double,2> pointBL = { m_points[minus1][0], m_points[minus1][1] };
                        sevRegions.push_back((ProteinRegion(uVec, m_prescribedSubLengths[minus1], activeSeverHeight, pointBL, 0, 0, 0, 1.0, 1, minus1)));

                        sevRegions[sevRegions.size()-1].stickToMem(*this);
                        addToSevReg(minus1, sevRegions.size()-1);
                    }

                    // if not already a region there, make one
                    if (activeBranch && m_coupBranchRegIDs[plus1].size() == 0)
                    {
                        std::array<double,2> uVec = { m_subunit_unitVecs[plus1][0], m_subunit_unitVecs[plus1][1] };
                        std::array<double,2> pointBL = { m_points[plus1][0], m_points[plus1][1] };
                        branchRegions.push_back((ProteinRegion(uVec, m_prescribedSubLengths[plus1], activeBranchHeight, pointBL, arpConc, 0, 0, 0, 1, plus1)));

                        branchRegions[branchRegions.size()-1].stickToMem(*this);
                        addToBranchReg(plus1, branchRegions.size()-1);
                    }

                    if (activeNuc && m_coupNucRegIDs[plus1].size() == 0)
                    {
                        std::array<double,2> uVec = { m_subunit_unitVecs[plus1][0], m_subunit_unitVecs[plus1][1] };
                        std::array<double,2> pointBL = { m_points[plus1][0], m_points[plus1][1] };
                        nucRegions.push_back((ProteinRegion(uVec, m_prescribedSubLengths[plus1], activeNucHeight, pointBL, 0, 0, 1.0, 0, 1, plus1)));

                        nucRegions[nucRegions.size()-1].stickToMem(*this);
                        addToNucReg(plus1, nucRegions.size()-1);
                    }

                    if (activeCap && m_coupCapRegIDs[plus1].size() == 0)
                    {
                        std::array<double,2> uVec = { m_subunit_unitVecs[plus1][0], m_subunit_unitVecs[plus1][1] };
                        std::array<double,2> pointBL = { m_points[plus1][0], m_points[plus1][1] };
                        capRegions.push_back((ProteinRegion(uVec, m_prescribedSubLengths[plus1], activeCapHeight, pointBL, 0, 1.0, 0, 0, 1, plus1)));

                        capRegions[capRegions.size()-1].stickToMem(*this);
                        addToCapReg(plus1, capRegions.size()-1);
                    }

                    if (activeAntiCap && m_coupAntiCapRegIDs[plus1].size() == 0)
                    {
                        std::array<double,2> uVec = { m_subunit_unitVecs[plus1][0], m_subunit_unitVecs[plus1][1] };
                        std::array<double,2> pointBL = { m_points[plus1][0], m_points[plus1][1] };
                        antiCapRegions.push_back((ProteinRegion(uVec, m_prescribedSubLengths[plus1], activeAntiCapHeight, pointBL, 0, 0, 0, 0, 1, plus1)));

                        antiCapRegions[antiCapRegions.size()-1].stickToMem(*this);
                        addToAntiCapReg(plus1, antiCapRegions.size()-1);
                    }

                    if (activeSever && m_coupSevRegIDs[plus1].size() == 0)
                    {
                        std::array<double,2> uVec = { m_subunit_unitVecs[plus1][0], m_subunit_unitVecs[plus1][1] };
                        std::array<double,2> pointBL = { m_points[plus1][0], m_points[plus1][1] };
                        sevRegions.push_back((ProteinRegion(uVec, m_prescribedSubLengths[plus1], activeSeverHeight, pointBL, 0, 0, 0, 1.0, 1, plus1)));

                        sevRegions[sevRegions.size()-1].stickToMem(*this);
                        addToSevReg(plus1, sevRegions.size()-1);
                    }

                    // Delete any tethers on sub j
                    int numTethers = m_tetherRelDists.size();
                    int k = 0;
                    while (k < numTethers)
                    {
                        if (m_tetherSubIDs[k] == prevSub || m_tetherSubIDs[k] == plus1 || m_tetherSubIDs[k] == minus1)
                        {
                            // Delete the tether
                            numTethers -= 1;
                            m_tetherSubIDs.erase(m_tetherSubIDs.begin()+k);
                            m_tetherRelDists.erase(m_tetherRelDists.begin()+k);
                            // Delete on other membrane
                            cortex.eraseTether(k);
                        }
                        else
                        {
                            ++k;
                        }
                    }
                }

            }

            // Need to move/update any regions that are coupled to subunits that
            // are adjacent to this point

            // Make region uVEc equal to subunit vector
            // Work out dx and dy of point j (membrane), apply this to the
            // bottom left point of region - i think this is wrong but ok
            double dx = r_new(j*2) - r_old(j*2);
            double dy = r_new(j*2 + 1) - r_old(j*2 + 1);

            for (unsigned int k = 0; k < m_coupBranchRegIDs[j].size(); ++k)
            {

                // Move it and rotate
                int branchID = m_coupBranchRegIDs[j][k];
                std::array<double,2> newUVec;
                newUVec[0] = m_subunit_unitVecs[branchRegions[branchID].getCoupledSub()][0];
                newUVec[1] = m_subunit_unitVecs[branchRegions[branchID].getCoupledSub()][1];

                std::array<double,2> newBLPoint = branchRegions[branchID].getBLPoint();
                newBLPoint[0] += dx;
                newBLPoint[1] += dy;

                branchRegions[branchID].moveRegionTo(newBLPoint, newUVec);
            }

            for (unsigned int k = 0; k < m_coupNucRegIDs[j].size(); ++k)
            {
                // Move it and rotate
                int nucID = m_coupNucRegIDs[j][k];
                std::array<double,2> newUVec;
                newUVec[0] = m_subunit_unitVecs[nucRegions[nucID].getCoupledSub()][0];
                newUVec[1] = m_subunit_unitVecs[nucRegions[nucID].getCoupledSub()][1];

                std::array<double,2> newBLPoint = nucRegions[nucID].getBLPoint();
                newBLPoint[0] += dx;
                newBLPoint[1] += dy;

                nucRegions[nucID].moveRegionTo(newBLPoint, newUVec);
            }

            for (unsigned int k = 0; k < m_coupCapRegIDs[j].size(); ++k)
            {
                // Move it and rotate

                int capID = m_coupCapRegIDs[j][k];
                std::array<double,2> newUVec;
                newUVec[0] = m_subunit_unitVecs[capRegions[capID].getCoupledSub()][0];
                newUVec[1] = m_subunit_unitVecs[capRegions[capID].getCoupledSub()][1];

                std::array<double,2> newBLPoint = capRegions[capID].getBLPoint();
                newBLPoint[0] += dx;
                newBLPoint[1] += dy;

                capRegions[capID].moveRegionTo(newBLPoint, newUVec);
            }

            for (unsigned int k = 0; k < m_coupAntiCapRegIDs[j].size(); ++k)
            {
                // Move it and rotate
                int antiCapID = m_coupAntiCapRegIDs[j][k];
                std::array<double,2> newUVec;
                newUVec[0] = m_subunit_unitVecs[antiCapRegions[antiCapID].getCoupledSub()][0];
                newUVec[1] = m_subunit_unitVecs[antiCapRegions[antiCapID].getCoupledSub()][1];

                std::array<double,2> newBLPoint = antiCapRegions[antiCapID].getBLPoint();
                newBLPoint[0] += dx;
                newBLPoint[1] += dy;

                antiCapRegions[antiCapID].moveRegionTo(newBLPoint, newUVec);
            }

            for (unsigned int k = 0; k < m_coupSevRegIDs[j].size(); ++k)
            {
                // Move it and rotate
                int SevID = m_coupSevRegIDs[j][k];
                std::array<double,2> newUVec;
                newUVec[0] = m_subunit_unitVecs[sevRegions[SevID].getCoupledSub()][0];
                newUVec[1] = m_subunit_unitVecs[sevRegions[SevID].getCoupledSub()][1];

                std::array<double,2> newBLPoint = sevRegions[SevID].getBLPoint();
                newBLPoint[0] += dx;
                newBLPoint[1] += dy;

                sevRegions[SevID].moveRegionTo(newBLPoint, newUVec);
            }

            int prevSub = j-1;
            if (prevSub == -1)
            {
                prevSub = m_numPoints - 1;
            }

            for (unsigned int k = 0; k < m_coupBranchRegIDs[prevSub].size(); ++k)
            {

                // Rotate it
                int branchID = m_coupBranchRegIDs[prevSub][k];
                std::array<double,2> newUVec;
                newUVec[0] = m_subunit_unitVecs[branchRegions[branchID].getCoupledSub()][0];
                newUVec[1] = m_subunit_unitVecs[branchRegions[branchID].getCoupledSub()][1];

                std::array<double,2> newBLPoint = branchRegions[branchID].getBLPoint();

                branchRegions[branchID].moveRegionTo(newBLPoint, newUVec);
            }

            for (unsigned int k = 0; k < m_coupNucRegIDs[prevSub].size(); ++k)
            {

                // Rotate it
                int nucID = m_coupNucRegIDs[prevSub][k];
                std::array<double,2> newUVec;
                newUVec[0] = m_subunit_unitVecs[nucRegions[nucID].getCoupledSub()][0];
                newUVec[1] = m_subunit_unitVecs[nucRegions[nucID].getCoupledSub()][1];

                std::array<double,2> newBLPoint = nucRegions[nucID].getBLPoint();

                nucRegions[nucID].moveRegionTo(newBLPoint, newUVec);
            }

            for (unsigned int k = 0; k < m_coupCapRegIDs[prevSub].size(); ++k)
            {

                // Rotate it
                int capID = m_coupCapRegIDs[prevSub][k];
                std::array<double,2> newUVec;
                newUVec[0] = m_subunit_unitVecs[capRegions[capID].getCoupledSub()][0];
                newUVec[1] = m_subunit_unitVecs[capRegions[capID].getCoupledSub()][1];

                std::array<double,2> newBLPoint = capRegions[capID].getBLPoint();

                capRegions[capID].moveRegionTo(newBLPoint, newUVec);
            }

            for (unsigned int k = 0; k < m_coupAntiCapRegIDs[prevSub].size(); ++k)
            {
                // Rotate it
                int antiCapID = m_coupAntiCapRegIDs[prevSub][k];
                std::array<double,2> newUVec;
                newUVec[0] = m_subunit_unitVecs[antiCapRegions[antiCapID].getCoupledSub()][0];
                newUVec[1] = m_subunit_unitVecs[antiCapRegions[antiCapID].getCoupledSub()][1];

                std::array<double,2> newBLPoint = antiCapRegions[antiCapID].getBLPoint();

                antiCapRegions[antiCapID].moveRegionTo(newBLPoint, newUVec);
            }

            for (unsigned int k = 0; k < m_coupSevRegIDs[prevSub].size(); ++k)
            {
                // Rotate it
                int SevID = m_coupSevRegIDs[prevSub][k];
                std::array<double,2> newUVec;
                newUVec[0] = m_subunit_unitVecs[sevRegions[SevID].getCoupledSub()][0];
                newUVec[1] = m_subunit_unitVecs[sevRegions[SevID].getCoupledSub()][1];

                std::array<double,2> newBLPoint = sevRegions[SevID].getBLPoint();

                sevRegions[SevID].moveRegionTo(newBLPoint, newUVec);
            }

        }
    }

    if (m_fusable && excZones.size() > 0)
    {
        int fuseType;
        std::array<double,2> fusePoint;
        std::array<int,3> fusePointsToAv;
        if (checkFusion(timesFused, excZones.size(), fuseType, fusePoint, fusePointsToAv))
        {
            doFusion(temp, membranes, stericGrid, timesFused,
                     excZones.size(), branchRegions, branching, capRegions,
                     nucRegions, antiCapRegions, sevRegions, actinvec, nActin,
                     steric, arpPool, arpGrid, excZones, memWalls, cortex,
                     tether, currTime, fusionTime, bDynamics, fuseType,
                     fusePoint, fusePointsToAv);
        }
    }

    // Option to have reference frame of centre of membrane - must not do for phagosome
    if (COMadjust && m_cellMem)
    {
        adjustCOM(actinvec, excZones, membranes, cortex, steric, stericGrid, branchRegions,
                  capRegions, antiCapRegions, nucRegions, sevRegions);
    }

    setCOM();

}

void Membrane::adjustCOM(std::vector<Actin> &actinvec,
                         std::vector<ExcZone> &excZones,
                         std::vector<Membrane> &membranes, Cortex &cortex,
                         const bool steric, StericGrid &stericGrid,
                         std::vector<ProteinRegion> &branchRegions,
                         std::vector<ProteinRegion> &capRegions,
                         std::vector<ProteinRegion> &antiCapRegions,
                         std::vector<ProteinRegion> &nucRegions,
                         std::vector<ProteinRegion> &sevRegions)
{
    // Move all points back to keep COM constant

    std::array<double,2> newCOM = calcCOM();
    std::array<double,2> dR;
    dR[0] = newCOM[0] - m_COM[0];
    dR[1] = newCOM[1] - m_COM[1];

    for (unsigned int i = 0; i < membranes.size(); ++i)
    {
        membranes[i].moveMem(dR, steric, stericGrid);
    }
    // Move Cortex if exists
    if (cortex.getExist())
    {
        cortex.moveCortex(dR, steric, stericGrid);
    }

    // Move and update steric grid for actin

    for (unsigned int i = 0; i < actinvec.size(); ++i)
    {
        actinvec[i].moveactin(-dR[0], -dR[1]);
        if (steric)
            stericGrid.resetAndUpdateAllDaughters(actinvec[i], actinvec);

        if (actinvec[i].getParentID() == -1 && actinvec[i].getDaughternum() != 0)
        {
          // Trunk of a branch, need ellipse centre moving
          std::array<double,3> centrePoint;
          centrePoint = actinvec[i].getEllipCentre();
          centrePoint[0] -= dR[0];
          centrePoint[1] -= dR[1];
          actinvec[i].setEllipCentre(centrePoint);
        }
    }

    // move all regions with COM
    for (unsigned int k = 0; k < branchRegions.size(); ++k)
    {
        branchRegions[k].moveRegion(-dR[0], -dR[1]);
    }
    for (unsigned int k = 0; k < nucRegions.size(); ++k)
    {
        nucRegions[k].moveRegion(-dR[0], -dR[1]);
    }
    for (unsigned int k = 0; k < antiCapRegions.size(); ++k)
    {
        antiCapRegions[k].moveRegion(-dR[0], -dR[1]);
    }
    for (unsigned int k = 0; k < capRegions.size(); ++k)
    {
        capRegions[k].moveRegion(-dR[0], -dR[1]);
    }
    for (unsigned int k = 0; k < sevRegions.size(); ++k)
    {
        sevRegions[k].moveRegion(-dR[0], -dR[1]);
    }
    // Move target(s) too!
    for (unsigned int i = 0; i < excZones.size(); ++i)
    {
        excZones[i].move(-dR[0], -dR[1]);
    }

}

void Membrane::moveMem(std::array<double,2> dR, const bool steric,
                       StericGrid &stericGrid)
{
    /*
      Function that simply moves the membrane for reference frame adjustment
     */

     for (int i = 0; i < m_numPoints; ++i)
     {
         m_points[i][0] -= dR[0];
         m_points[i][1] -= dR[1];
     }

     // Update steric grid for membrane
     if (steric)
     {
         for (int i = 0; i < m_numPoints; ++i)
         {
            stericGrid.resetMemandUpdate(*this, i);
         }
     }

}

bool Membrane::checkExVolSUB(const std::vector<Actin> &actinvec, int pointID)
{
    /*
    Function that should be called after fluctuating the membrane, to see if
    that move has violated excluded volume

    Check against all actin filaments

    Returns true if exc vol is violated
    */

    // First define our membrane
    gte::Segment<2,double> membrane;
    gte::Segment<2,double> filament;
    int start = pointID - 1;
    int end = pointID + 1;

    if (pointID == 0)
    {
        start = 0;
    }
    else if (pointID == m_numPoints-1)
    {
        end = pointID;
    }

    for (int i = start; i < end; ++i)
    {
        // The membrane
        membrane.p[0][0] = m_points[i][0];
        membrane.p[1][0] = m_points[i+1][0];
        membrane.p[0][1] = m_points[i][1];
        membrane.p[1][1] = m_points[i+1][1];

        // Now define our actin filament
        for (unsigned int j = 0; j < actinvec.size(); ++j)
        {
            for (int k = 0; k < actinvec[j].getNumSubs(); ++k)
            {
                filament.p[0][0] = actinvec[j].getPoints()[k][0];
                filament.p[1][0] = actinvec[j].getPoints()[k+1][0];
                filament.p[0][1] = actinvec[j].getPoints()[k][1];
                filament.p[1][1] = actinvec[j].getPoints()[k+1][1];


                SegDistQuery min_d_GTE;
                auto result = min_d_GTE(membrane, filament);
                double distance = result.distance;

                if (distance < (m_thickness/2 + actinvec[j].getStericRadius()))
                {
                    return true;
                }

            }
        }
    }

    if (pointID == 0 || pointID == m_numPoints-1)
    {
        // Also should check the spring
        // The membrane
        membrane.p[0][0] = m_points[0][0];
        membrane.p[1][0] = m_points[m_numPoints-1][0];
        membrane.p[0][1] = m_points[0][1];
        membrane.p[1][1] = m_points[m_numPoints-1][1];

        // Now define our actin filament
        for (unsigned int j = 0; j < actinvec.size(); ++j)
        {
            for (int k = 0; k < actinvec[j].getNumSubs(); ++k)
            {
                filament.p[0][0] = actinvec[j].getPoints()[k][0];
                filament.p[1][0] = actinvec[j].getPoints()[k+1][0];
                filament.p[0][1] = actinvec[j].getPoints()[k][1];
                filament.p[1][1] = actinvec[j].getPoints()[k+1][1];


                SegDistQuery min_d_GTE;
                auto result = min_d_GTE(membrane, filament);
                double distance = result.distance;

                if (distance < (m_thickness/2 + actinvec[j].getStericRadius()))
                {
                    return true;
                }

            }
        }
    }

    return false;

}

bool Membrane::checkExVolSUBGrid(const std::vector<Actin> &actinvec,
                                         int pointID, StericGrid &stericGrid)
{
    /*
    Function that should be called after fluctuating the membrane, to see if
    that move has violated excluded volume

    Check against all actin filaments

    Returns true if exc vol is violated

    */


    std::array<int,3> pointsToCheck;

    int start = pointID - 1;

    if (pointID == 0)
    {
        start = m_numPoints-1;
    }

    pointsToCheck = {start, pointID};
    // Do we need to check the end sub? surely just start and pointID??
    for (int i = 0; i < 2; ++i)
    {
        int p = pointsToCheck[i];
        std::vector<int> occCells; // cells that are occupied

        for (unsigned int j = 0; j < m_stericCells[p].size(); ++j)
        {
            occCells.push_back(m_stericCells[p][j]);
        }
        std::vector<int> cellsToCheck;
        cellsToCheck = stericGrid.getCellsToCheck(occCells);


        if (check_min_dist_Grid(p, actinvec, cellsToCheck, stericGrid))
            return true;
    }

    return false;

}

bool Membrane::checkExVolSUBGridMembrane(const std::vector<Membrane> &membranes,
                                             int pointID, StericGrid &stericGrid)
{
    /*
    Function that should be called after fluctuating the membrane, to see if
    that move has violated excluded volume by checking against other membranes

    Check against all actin filaments

    Returns true if exc vol is violated
    */


    std::array<int,3> pointsToCheck;

    int start = pointID - 1;

    if (pointID == 0)
    {
        start = m_numPoints-1;
    }

    pointsToCheck = {start, pointID};

    gte::Segment<2,double> membrane;
    for (int i = 0; i < 2; ++i)
    {
        if (pointsToCheck[i] == m_numPoints-1)
        {
            // its the end sub
            membrane.p[0][0] = m_points[0][0];
            membrane.p[1][0] = m_points[m_numPoints-1][0];
            membrane.p[0][1] = m_points[0][1];
            membrane.p[1][1] = m_points[m_numPoints-1][1];
        }
        else
        {
            membrane.p[0][0] = m_points[pointsToCheck[i]][0];
            membrane.p[1][0] = m_points[pointsToCheck[i]+1][0];
            membrane.p[0][1] = m_points[pointsToCheck[i]][1];
            membrane.p[1][1] = m_points[pointsToCheck[i]+1][1];
        }


        // To avoid repeated checks against the same subunits that may appear in
        // more than one cell, have a recorded of which ones we have checked
        std::vector<std::array<int,2>> checkedSubs;
        int numMembranes = membranes.size();

        for (int j = 0; j < numMembranes; ++j)
        {
            if (j == m_id)
            {
              // Dont check against self!
                continue;
            }

            for (int k = 0; k < membranes[j].getNumPoints(); ++k)
            {
                gte::Segment<2,double> other;

                other.p[0][0] = membranes[j].getPoints()[k][0];
                other.p[0][1] = membranes[j].getPoints()[k][1];

                if (k == membranes[j].getNumPoints() -1 )
                {
                    // end sub
                    other.p[1][0] = membranes[j].getPoints()[0][0];
                    other.p[1][1] = membranes[j].getPoints()[0][1];
                }
                else
                {
                    other.p[1][0] = membranes[j].getPoints()[k+1][0];
                    other.p[1][1] = membranes[j].getPoints()[k+1][1];
                }

                SegDistQuery min_d_GTE;
                auto result = min_d_GTE(membrane, other);
                double distance = result.distance;

                if (distance < (m_thickness/2 + membranes[j].getThickness()/2))
                {
                    return true;
                }
            }
        }
    }

    return false;

}

bool Membrane::checkExVolSUBGridCortex(const Cortex &cortex,
                                       int pointID, StericGrid &stericGrid)
{
    /*
    Function that should be called after fluctuating the membrane, to see if
    that move has violated excluded volume by checking against other membranes

    Check against all actin filaments

    Returns true if exc vol is violated
    */
    if (!cortex.getExist())
      return false;

    std::array<int,3> pointsToCheck;

    int start = pointID - 1;

    if (pointID == 0)
    {
        start = m_numPoints-1;
    }

    pointsToCheck = {start, pointID};


    gte::Segment<2,double> membrane;
    for (int i = 0; i < 2; ++i)
    {
        if (pointsToCheck[i] == m_numPoints-1)
        {
            // its the end sub
            membrane.p[0][0] = m_points[0][0];
            membrane.p[1][0] = m_points[m_numPoints-1][0];
            membrane.p[0][1] = m_points[0][1];
            membrane.p[1][1] = m_points[m_numPoints-1][1];
        }
        else
        {
            membrane.p[0][0] = m_points[pointsToCheck[i]][0];
            membrane.p[1][0] = m_points[pointsToCheck[i]+1][0];
            membrane.p[0][1] = m_points[pointsToCheck[i]][1];
            membrane.p[1][1] = m_points[pointsToCheck[i]+1][1];
        }

        // To avoid repeated checks against the same subunits that may appear in
        // more than one cell, have a recorded of which ones we have checked
        std::vector<std::array<int,2>> checkedSubs;

        for (int k = 0; k < cortex.getNumPoints(); ++k)
        {
            gte::Segment<2,double> cortexSub;

            cortexSub.p[0][0] = cortex.getPoints()[k][0];
            cortexSub.p[0][1] = cortex.getPoints()[k][1];

            if (k == cortex.getNumPoints() -1 )
            {
                // end sub
                cortexSub.p[1][0] = cortex.getPoints()[0][0];
                cortexSub.p[1][1] = cortex.getPoints()[0][1];
            }
            else
            {
                cortexSub.p[1][0] = cortex.getPoints()[k+1][0];
                cortexSub.p[1][1] = cortex.getPoints()[k+1][1];
            }

            SegDistQuery min_d_GTE;
            auto result = min_d_GTE(membrane, cortexSub);
            double distance = result.distance;

            if (distance < (m_thickness/2 + cortex.getThickness()/2))
            {
                return true;
            }
        }

    }

    return false;

}

bool Membrane::check_min_dist_Grid(int subid, const std::vector<Actin> &actinvec,
                                              std::vector<int> &cellsToCheck,
                                              StericGrid &stericGrid)
{
    gte::Segment<2,double> membrane;

    if (subid == m_numPoints-1)
    {
        // its the end sub
        membrane.p[0][0] = m_points[0][0];
        membrane.p[1][0] = m_points[m_numPoints-1][0];
        membrane.p[0][1] = m_points[0][1];
        membrane.p[1][1] = m_points[m_numPoints-1][1];
    }
    else
    {
        membrane.p[0][0] = m_points[subid][0];
        membrane.p[1][0] = m_points[subid+1][0];
        membrane.p[0][1] = m_points[subid][1];
        membrane.p[1][1] = m_points[subid+1][1];
    }


    // To avoid repeated checks against the same subunits that may appear in
    // more than one cell, have a recorded of which ones we have checked
    std::vector<std::array<int,2>> checkedSubs;

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
            auto result = min_d_GTE(membrane, filament);
            double distance = result.distance;

            if (distance < (m_thickness/2 + actinvec[actinID].getStericRadius()))
            {
                return true;
            }
        }
    }
    return false;
}

std::array<int,2> Membrane::check_min_dist_Grid_Identify(int subid, const std::vector<Actin> &actinvec,
                                              std::vector<int> &cellsToCheck,
                                              StericGrid &stericGrid)
{
    /*
     * Same as above but returns the actinID and actinSubID of the filament sub
     * that breaks Steric hindrance. If steric hindrance is not broken
     * returns {-1,-1}
     */

    gte::Segment<2,double> membrane;


    if (subid == m_numPoints-1)
    {
        // its the end sub
        membrane.p[0][0] = m_points[0][0];
        membrane.p[1][0] = m_points[m_numPoints-1][0];
        membrane.p[0][1] = m_points[0][1];
        membrane.p[1][1] = m_points[m_numPoints-1][1];
    }
    else
    {
        membrane.p[0][0] = m_points[subid][0];
        membrane.p[1][0] = m_points[subid+1][0];
        membrane.p[0][1] = m_points[subid][1];
        membrane.p[1][1] = m_points[subid+1][1];
    }


    // To avoid repeated checks against the same subunits that may appear in
    // more than one cell, have a recorded of which ones we have checked
    std::vector<std::array<int,2>> checkedSubs;

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
            auto result = min_d_GTE(membrane, filament);
            double distance = result.distance;

            if (distance < (m_thickness/2 + actinvec[actinID].getStericRadius()))
            {
                return {actinID, actinSubID};
            }
        }
    }
    return {-1, -1};
}

bool Membrane::s_checkExVolBarbPoly(const Actin &filament,
                                        std::vector<Membrane> &membranes,
                                        StericGrid &stericGrid)
{
    for (unsigned int i = 0; i < membranes.size(); ++i)
    {
        if (membranes[i].getExist() && membranes[i].checkExVolBarbPolyGrid(filament, stericGrid))
            return true;
    }

    return false;
}

bool Membrane::checkExVolBarbPoly(const Actin &filament) const
{
    /*
    Function that should be called during barbed end polymerisation of actin filaments

    This checks the new monomer against the whole Membrane

    Returns true if exc vol is violated
    */

    // Define our actin monomer
    std::array<double,2> unitVec = filament.getUnitVec(filament.getNumSubs()-1);
    double newBarbx = filament.getBarbedEnd()[0] + (Actin::s_monomerLength * unitVec[0]);
    double newBarby = filament.getBarbedEnd()[1] + (Actin::s_monomerLength * unitVec[1]);
    gte::Segment<2,double> monomer;
    monomer.p[0][0] = filament.getBarbedEnd()[0];
    monomer.p[1][0] = newBarbx;
    monomer.p[0][1] = filament.getBarbedEnd()[1];
    monomer.p[1][1] = newBarby;

    // define our membrane
    gte::Segment<2,double> membrane;
    for (int i = 0; i < m_numPoints; ++i)
    {
        // The membrane
        int np = i+1;
        if (i == m_numPoints-1)
            np = 0;

        membrane.p[0][0] = m_points[i][0];
        membrane.p[1][0] = m_points[np][0];
        membrane.p[0][1] = m_points[i][1];
        membrane.p[1][1] = m_points[np][1];



        SegDistQuery min_d_GTE;
        auto result = min_d_GTE(membrane, monomer);
        double distance = result.distance;

        if (distance < (m_thickness/2 + filament.getStericRadius()))
        {
            return true;
        }
    }

    return false;
}


bool Membrane::checkExVolBarbPolyGrid(const Actin &filament,
                                              StericGrid &stericGrid) const
{
    /*
    Function that should be called during barbed end polymerisation of actin filaments

    This checks the new monomer against the whole Membrane

    Returns true if exc vol is violated
    */

    // Define our actin end sub
    gte::Segment<2,double> endSub;
    endSub.p[0][0] = filament.getBarbedEnd()[0];
    endSub.p[1][0] = filament.getPreBarbedEnd()[0];
    endSub.p[0][1] = filament.getBarbedEnd()[1];
    endSub.p[1][1] = filament.getPreBarbedEnd()[1];

    std::vector<int> barbCellIDs = filament.getStericCells(filament.getNumSubs()-1);
    std::vector<int> cellsToCheck;
    cellsToCheck = stericGrid.getCellsToCheck(barbCellIDs);

    for (unsigned int i = 0; i < cellsToCheck.size(); ++i)
    {
        //std::cout << cellsToCheck[i] << std::endl;
        std::vector<std::array<int,2> > cellContents = stericGrid.getCellContentsMem(cellsToCheck[i]);
        for (unsigned int j = 0; j < cellContents.size(); ++j)
        {

            int memSubID = cellContents[j][1];
            gte::Segment<2,double> membrane;

            if (memSubID == m_numPoints-1)
            {
                // the spring
                membrane.p[0][0] = m_points[memSubID][0];
                membrane.p[1][0] = m_points[0][0];
                membrane.p[0][1] = m_points[memSubID][1];
                membrane.p[1][1] = m_points[0][1];
            }
            else
            {
                membrane.p[0][0] = m_points[memSubID][0];
                membrane.p[1][0] = m_points[memSubID+1][0];
                membrane.p[0][1] = m_points[memSubID][1];
                membrane.p[1][1] = m_points[memSubID+1][1];
            }

            SegDistQuery min_d_GTE;
            auto result = min_d_GTE(membrane, endSub);
            double distance = result.distance;

            if (distance < (m_thickness/2 + filament.getStericRadius()))
            {
                return true;
            }
        }
    }

    return false;
}
bool Membrane::s_checkExVolPointPoly(const Actin &filament,
                                         std::vector<Membrane> &membranes,
                                         StericGrid &stericGrid)
{
    for (unsigned int i = 0; i < membranes.size(); ++i)
    {
        if (membranes[i].getExist() && membranes[i].checkExVolPointPolyGrid(filament, stericGrid))
            return true;
    }

    return false;

}

bool Membrane::checkExVolPointPoly(const Actin &filament) const
{
    /*
    Function that should be called during pointed end polymerisation of actin filaments

    This checks the new monomer against the whole Membrane

    Returns true if exc vol is violated
    */

    // Define our actin monomer

    gte::Segment<2,double> firstSub;
    firstSub.p[0][0] = filament.getPointedEnd()[0];
    firstSub.p[1][0] = filament.getPrePointedEnd()[0];
    firstSub.p[0][1] = filament.getPointedEnd()[1];
    firstSub.p[1][1] = filament.getPrePointedEnd()[1];


    // define our membrane
    gte::Segment<2,double> membrane;
    for (int i = 0; i < m_numPoints; ++i)
    {
        // The membrane
        int np = i+1;
        if (i == m_numPoints-1)
            np = 0;

        membrane.p[0][0] = m_points[i][0];
        membrane.p[1][0] = m_points[np][0];
        membrane.p[0][1] = m_points[i][1];
        membrane.p[1][1] = m_points[np][1];


        SegDistQuery min_d_GTE;
        auto result = min_d_GTE(membrane, firstSub);
        double distance = result.distance;

        if (distance < (m_thickness/2 + filament.getStericRadius()))
        {
            return true;
        }
    }

    return false;
}

bool Membrane::checkExVolPointPolyGrid(const Actin &filament,
                                               StericGrid &stericGrid) const
{
    // Define the first sub
    gte::Segment<2,double> firstSub;
    firstSub.p[0][0] = filament.getPointedEnd()[0];
    firstSub.p[1][0] = filament.getPrePointedEnd()[0];
    firstSub.p[0][1] = filament.getPointedEnd()[1];
    firstSub.p[1][1] = filament.getPrePointedEnd()[1];

    std::vector<int> PointCellIDs = filament.getStericCells(0);
    std::vector<int> cellsToCheck;
    cellsToCheck = stericGrid.getCellsToCheck(PointCellIDs);


    for (unsigned int i = 0; i < cellsToCheck.size(); ++i)
    {
        std::vector<std::array<int,2>> cellContents = stericGrid.getCellContentsMem(cellsToCheck[i]);
        for (unsigned int j = 0; j < cellContents.size(); ++j)
        {
            int memSubID = cellContents[j][1];
            gte::Segment<2,double> membrane;

            if (memSubID == m_numPoints-1)
            {
                // the spring
                membrane.p[0][0] = m_points[memSubID][0];
                membrane.p[1][0] = m_points[0][0];
                membrane.p[0][1] = m_points[memSubID][1];
                membrane.p[1][1] = m_points[0][1];
            }
            else
            {
                membrane.p[0][0] = m_points[memSubID][0];
                membrane.p[1][0] = m_points[memSubID+1][0];
                membrane.p[0][1] = m_points[memSubID][1];
                membrane.p[1][1] = m_points[memSubID+1][1];
            }

            SegDistQuery min_d_GTE;
            auto result = min_d_GTE(membrane, firstSub);
            double distance = result.distance;

            if (distance < (m_thickness/2 + filament.getStericRadius()))
            {
                return true;
            }
        }
    }

    return false;

}

bool Membrane::checkExVolNuc(const Actin &actin) const
{
    /*
    Function that carries out a distance check between the membrane
    and the new actin filament (either a free trimer or a branch)

    Returns true if exc vol is violated
    */

    // Define our new filament
    gte::Segment<2,double> filament;
    filament.p[0][0] = actin.getPoints()[0][0];
    filament.p[1][0] = actin.getPoints()[4][0];
    filament.p[0][1] = actin.getPoints()[0][1];
    filament.p[1][1] = actin.getPoints()[4][1];

    // Define our membrane wall
    gte::Segment<2,double> membrane;
    for (int i = 0; i < m_numPoints; ++i)
    {
        // The membrane
        int np = i+1;
        if (i == m_numPoints-1)
            i = 0;

        membrane.p[0][0] = m_points[i][0];
        membrane.p[1][0] = m_points[np][0];
        membrane.p[0][1] = m_points[i][1];
        membrane.p[1][1] = m_points[np][1];

        SegDistQuery min_d_GTE;
        auto result = min_d_GTE(filament, membrane);
        double distance = result.distance;

        if (distance < (m_thickness/2 + actin.getStericRadius()))
        {
            return true;
        }
    }

    return false;

}

bool Membrane::s_checkExVolNuc(const Actin &actin,
                                   const std::vector<Membrane> &membranes,
                                   std::vector<int> cellsToCheck,
                                   StericGrid &stericGrid)
{
    for (unsigned int i = 0; i < membranes.size(); ++i)
    {
        if (membranes[i].getExist() && membranes[i].checkExVolNuc_Grid(actin, cellsToCheck, stericGrid))
            return true;
    }

    return false;

}

bool Membrane::checkExVolNuc_Grid(const Actin &actin,
                                          std::vector<int> cellsToCheck,
                                          StericGrid &stericGrid) const
{



    gte::Segment<2,double> filament;
    filament.p[0][0] = actin.getPoints()[0][0];
    filament.p[1][0] = actin.getPoints()[4][0];
    filament.p[0][1] = actin.getPoints()[0][1];
    filament.p[1][1] = actin.getPoints()[4][1];


    // To avoid repeated checks against the same subunits that may appear in
    // more than one cell, have a recorded of which ones we have checked
    std::vector<int> checkedSubs;

    // Need to check that it is inside the membrane!
    if (!checkInOut(actin))
    {
        return true;
    }

    for (unsigned int i = 0; i < cellsToCheck.size(); ++i)
    {

        std::vector<std::array<int,2>> cellContents = stericGrid.getCellContentsMem(cellsToCheck[i]);

        for (unsigned int j = 0; j < cellContents.size(); ++j)
        {
            // j is membrane sub

            if (std::find(checkedSubs.begin(), checkedSubs.end(), cellContents[j][1]) != checkedSubs.end())
            {
                // We have checked this mem sub before!
                continue;
            }
            checkedSubs.push_back(cellContents[j][1]);

            int memSubID = cellContents[j][1];
            gte::Segment<2,double> membrane;

            if (memSubID == m_numPoints-1)
            {
                // the spring
                membrane.p[0][0] = m_points[memSubID][0];
                membrane.p[1][0] = m_points[0][0];
                membrane.p[0][1] = m_points[memSubID][1];
                membrane.p[1][1] = m_points[0][1];
            }
            else
            {
                membrane.p[0][0] = m_points[memSubID][0];
                membrane.p[1][0] = m_points[memSubID+1][0];
                membrane.p[0][1] = m_points[memSubID][1];
                membrane.p[1][1] = m_points[memSubID+1][1];
            }

            SegDistQuery min_d_GTE;
            auto result = min_d_GTE(filament, membrane);
            double distance = result.distance;

            if (distance < (m_thickness/2 + actin.getStericRadius()))
            {
                return true;
            }
        }
    }
    return false;
}

bool Membrane::checkExVolBM(const Actin &actin) const
{
    /*
    Function that calculates the distance between all of the filament and the bent
    Membrane. For now let's just have this as checking for an intersection
    This means that filaments will be able to jump across though, so need to
    address this
    */


    if (!checkInOut(actin))
    {
        return true;
    }

    // First define our membrane
    gte::Segment<2,double> membrane;
    gte::Segment<2,double> filament;
    for (int i = 0; i < m_numPoints; ++i)
    {
        // The membrane
        int np = i+1; // next point
        if (i == m_numPoints-1)
            np = 0;

        membrane.p[0][0] = m_points[i][0];
        membrane.p[1][0] = m_points[np][0];
        membrane.p[0][1] = m_points[i][1];
        membrane.p[1][1] = m_points[np][1];


        // Now define our actin filament
        for (int j = 0; j < actin.getNumSubs(); ++j)
        {
            filament.p[0][0] = actin.getPoints()[j][0];
            filament.p[1][0] = actin.getPoints()[j+1][0];
            filament.p[0][1] = actin.getPoints()[j][1];
            filament.p[1][1] = actin.getPoints()[j+1][1];


            SegDistQuery min_d_GTE;
            auto result = min_d_GTE(filament, membrane);
            double distance = result.distance;

            if (distance < (m_thickness/2 + actin.getStericRadius()))
            {
                return true;
            }

        }
    }

    return false;
}

bool Membrane::s_checkExVolBM(int subid, const Actin &actin,
                                  const std::vector<Membrane> &membranes,
                                  std::vector<int> cellsToCheck,
                                  StericGrid &stericGrid)
{
    for (unsigned int i = 0; i < membranes.size(); ++i)
    {
        if (membranes[i].getExist() && membranes[i].check_min_dist_Actin_Grid(subid, actin, cellsToCheck, stericGrid))
        {
            return true;
        }
    }

    return false;
}


bool Membrane::check_min_dist_Actin_Grid(int subid, const Actin &actin,
                                                 std::vector<int> cellsToCheck,
                                                 StericGrid &stericGrid) const
{
    /*
    Function that calculates the distance between all of the filament and the bent
    Membrane. For now let's just have this as checking for an intersection
    This means that filaments will be able to jump across though, so need to
    address this
    */


    if (!checkInOut(actin))
    {
        return true;
    }


    gte::Segment<2,double> filament;
    filament.p[0][0] = actin.getPoints()[subid][0];
    filament.p[1][0] = actin.getPoints()[subid+1][0];
    filament.p[0][1] = actin.getPoints()[subid][1];
    filament.p[1][1] = actin.getPoints()[subid+1][1];


    // To avoid repeated checks against the same subunits that may appear in
    // more than one cell, have a recorded of which ones we have checked
    std::vector<int> checkedSubs;


    for (unsigned int i = 0; i < cellsToCheck.size(); ++i)
    {

        std::vector<std::array<int,2>> cellContents = stericGrid.getCellContentsMem(cellsToCheck[i]);

        for (unsigned int j = 0; j < cellContents.size(); ++j)
        {
            // j is membrane sub

            if (std::find(checkedSubs.begin(), checkedSubs.end(), cellContents[j][1]) != checkedSubs.end())
            {
                // We have checked this mem sub before!
                continue;
            }
            checkedSubs.push_back(cellContents[j][1]);

            int memSubID = cellContents[j][1];
            gte::Segment<2,double> membrane;

            if (memSubID == m_numPoints-1)
            {
                // the spring
                membrane.p[0][0] = m_points[memSubID][0];
                membrane.p[1][0] = m_points[0][0];
                membrane.p[0][1] = m_points[memSubID][1];
                membrane.p[1][1] = m_points[0][1];
            }
            else
            {
                membrane.p[0][0] = m_points[memSubID][0];
                membrane.p[1][0] = m_points[memSubID+1][0];
                membrane.p[0][1] = m_points[memSubID][1];
                membrane.p[1][1] = m_points[memSubID+1][1];
            }

            SegDistQuery min_d_GTE;
            auto result = min_d_GTE(filament, membrane);
            double distance = result.distance;

            if (distance < (m_thickness/2 + actin.getStericRadius()))
            {
                return true;
            }
        }
    }
    return false;
}

bool Membrane::checkInOut(const Actin &actin) const
{
    /*
    Uses the crossing point algorithm to determine if a given point
    (so the centre) of our actin filament, lies in or out of the filament
    If it is in return true if out return false
    Used as an extra check to ensure no filaments "pop-out" of our "cell"
    */

    // Should only work on cell membrane NOT PHAGOSOME!


    // Define a point we know with certainity will be outside of our cell,
    // here it is at (1m, 0m)
    std::array<double,2> outPoint = {1,0};
    outPoint[0] += m_COM[0];
    outPoint[1] += m_COM[1];
    int junk; // not needed
    std::array<double,2> actinP = actin.findCentrePoint(junk);

    // Use geometry tools to find intersections
    // define our "ray"

    gte::Segment<2,double> ray;

    ray.p[0][0] = outPoint[0];
    ray.p[1][0] = actinP[0];
    ray.p[0][1] = outPoint[1];
    ray.p[1][1] = actinP[1];

    // Now we need to check this against our whole membrane and count intersections

    gte::Segment<2,double> membrane;
    int numIntersections = 0;
    for (int i = 0; i < m_numPoints; ++i)
    {
        // The membrane
        int np = i+1;
        if (i == m_numPoints-1)
            np = 0;
        membrane.p[0][0] = m_points[i][0];
        membrane.p[1][0] = m_points[np][0];
        membrane.p[0][1] = m_points[i][1];
        membrane.p[1][1] = m_points[np][1];

        Check_I IntersectQuery;
        auto result = IntersectQuery(ray, membrane);
        numIntersections += result.numIntersections;
    }


    if (numIntersections % 2 == 0)
    {
        // even, so outside
        if (m_cellMem)
        {
            return false;
        }
        else
        {
            return true;
        }
    }
    else
    {
        if (m_cellMem)
        {
            return true;
        }
        else
        {
            return false;
        }
    }

}

bool Membrane::checkInOutAll(const std::vector<Actin> &actinVec) const
{
    for (unsigned int i = 0; i < actinVec.size(); ++i)
    {
      if (!checkInOut(actinVec[i]))
      {
          std::cout << "FIlament outside, abort" << std::endl;
          return false;
      }
    }
    return true;
}

bool Membrane::checkPointIn(std::array<double,2> point) const
{
    /*
    Uses the crossing point algorithm to determine if a given point
    (so the centre) of our actin filament, lies in or out of the filament
    If it is in return true if out return false
    Used as an extra check to ensure no filaments "pop-out" of our "cell"
    */

    // Should only work on cell membrane NOT PHAGOSOME!


    // Define a point we know with certainity will be outside of our cell,
    // here it is at (1m, 0m)
    std::array<double,2> outPoint = {1,0};
    outPoint[0] += m_COM[0];
    outPoint[1] += m_COM[1];
    // Use geometry tools to find intersections
    // define our "ray"

    gte::Segment<2,double> ray;

    ray.p[0][0] = outPoint[0];
    ray.p[1][0] = point[0];
    ray.p[0][1] = outPoint[1];
    ray.p[1][1] = point[1];

    // Now we need to check this against our whole membrane and count intersections

    gte::Segment<2,double> membrane;
    int numIntersections = 0;
    for (int i = 0; i < m_numPoints; ++i)
    {
        // The membrane
        int np = i+1;
        if (i == m_numPoints-1)
            np = 0;
        membrane.p[0][0] = m_points[i][0];
        membrane.p[1][0] = m_points[np][0];
        membrane.p[0][1] = m_points[i][1];
        membrane.p[1][1] = m_points[np][1];

        Check_I IntersectQuery;
        auto result = IntersectQuery(ray, membrane);
        numIntersections += result.numIntersections;
    }


    if (numIntersections % 2 == 0)
    {
        // even, so outside
        return false;
    }
    else
    {
        return true;
    }

}

bool Membrane::calcDistToTarsForActive(const std::vector<ExcZone> &excZones,
                                           int subID)
{
    for (unsigned int i = 0; i < excZones.size(); ++i)
    {
        if (excZones[i].getCircBool())
        {
            // Define centre point of circle
            gte::Vector2<double> point = {excZones[i].getX(), excZones[i].getY()};

            // Define our membrane wall
            gte::Segment<2,double> membraneSeg;
            int endP = subID + 1;
            if (endP == m_numPoints)
            {
                endP = 0;
            }
            membraneSeg.p[0][0] = m_points[subID][0];
            membraneSeg.p[1][0] = m_points[endP][0];
            membraneSeg.p[0][1] = m_points[subID][1];
            membraneSeg.p[1][1] = m_points[endP][1];


            DistQuery min_d_GTE;
            auto result = min_d_GTE(point, membraneSeg);
            double distance = result.distance;
            if (distance < (m_thickness/2 + excZones[i].getR() + m_activeDist))
            {
                // membrane segment is withing threshold distance and should be activated
                return true;
            }

        }
        else
        {
            std::cout << i << std::endl;
            std::cout << "Membrane.cpp line 3016 Target not circle, exiting" << std::endl;
            exit(0);
        }
    }

    return false;
}

bool Membrane::checkFusion(unsigned int &timesFused,
                           unsigned int numTars, int &fuseType,
                           std::array<double,2> &fusePoint,
                           std::array<int,3> &fusePointsToAv)
{
    // Check whether membrane should fuse
    // Subroutine to check if the membrane should fuse and at which point
    // Uses average position of the triangle to create fusionPoint
    // Removes three points and replaces with shared point

    if (timesFused == numTars)
    {
        // Hard limit on fusion
        return false;
    }

    if (!m_cellMem)
    {
      // Do not allow fusion of phagosome
        return false;
    }

    double distThreshold = m_thickness; // so its touching

    bool fusion = false;
    // Calculate distances between non-adjacent points along the membrane

    for (int i = 0; i < m_numPoints; ++i)
    {
        for (int j = i+2; j < m_numPoints; ++j)
        {
            if (i==0 && j == m_numPoints-1)
            {
                // Adjacent subs
                continue;
            }

            gte::Segment<2,double> membraneSeg1;
            gte::Segment<2,double> membraneSeg2;

            int nextPI = i+1;
            if (nextPI == m_numPoints)
            {
                nextPI = 0;
            }

            int nextPJ = j+1;
            if (nextPJ == m_numPoints)
            {
                nextPJ = 0;
            }

            membraneSeg1.p[0][0] = m_points[i][0];
            membraneSeg1.p[1][0] = m_points[nextPI][0];
            membraneSeg1.p[0][1] = m_points[i][1];
            membraneSeg1.p[1][1] = m_points[nextPI][1];

            membraneSeg2.p[0][0] = m_points[j][0];
            membraneSeg2.p[1][0] = m_points[nextPJ][0];
            membraneSeg2.p[0][1] = m_points[j][1];
            membraneSeg2.p[1][1] = m_points[nextPJ][1];

            SegDistQuery min_d_GTE;
            auto result = min_d_GTE(membraneSeg1, membraneSeg2);
            double distance = result.distance;

            if (distance < distThreshold)
            {
                fusion = true;
                // Should fuse
                // Where is the fusion point?
                // Should be one of four endpoints

                gte::Vector<2,double> pointA;
                pointA = membraneSeg1.p[0]; // i
                gte::Vector<2,double> pointB;
                pointB = membraneSeg1.p[1]; // nextPI

                gte::Vector<2,double> pointC;
                pointC = membraneSeg2.p[0]; // j
                gte::Vector<2,double> pointD;
                pointD = membraneSeg2.p[1]; // nextPJ

                DistQuery min_d; // distance between point and segment

                // point i and j-j+1 sub
                auto result = min_d(pointA, membraneSeg2);
                double minDist = result.distance;

                fusePointsToAv[0] = i;
                fusePointsToAv[1] = j;
                fusePointsToAv[2] = nextPJ;

                fuseType = 0;
                // point nextPI and j-j+1 sub
                result = min_d(pointB, membraneSeg2);
                if (result.distance < minDist)
                {
                    minDist = result.distance;

                    fusePointsToAv[0] = nextPI;
                    fusePointsToAv[1] = j;
                    fusePointsToAv[2] = nextPJ;

                    fuseType = 1;
                }

                // point j and i-i+1 sub
                result = min_d(pointC, membraneSeg1);
                if (result.distance < minDist)
                {
                    minDist = result.distance; // point j

                    fusePointsToAv[0] = i;
                    fusePointsToAv[1] = nextPI;
                    fusePointsToAv[2] = j;

                    fuseType = 2;

                }

                // point nextPJ and i-i+1 sub
                result = min_d(pointD, membraneSeg1);
                if (result.distance < minDist)
                {
                    minDist = result.distance; // point nextPJ

                    fusePointsToAv[0] = i;
                    fusePointsToAv[1] = nextPI;
                    fusePointsToAv[2] = nextPJ;
                    fuseType = 3;
                }

                fusePoint[0] = (m_points[fusePointsToAv[0]][0]
                               + m_points[fusePointsToAv[1]][0]
                               + m_points[fusePointsToAv[2]][0]) / 3;

                fusePoint[1] = (m_points[fusePointsToAv[0]][1]
                              + m_points[fusePointsToAv[1]][1]
                              + m_points[fusePointsToAv[2]][1]) / 3;

                break;
            }
        }

        if (fusion)
            break;
    }

    return fusion;
}

void Membrane::doFusion(double temp, std::vector<Membrane> &membranes,
                           StericGrid &stericGrid, unsigned int &timesFused,
                           unsigned int numTars,
                           std::vector<ProteinRegion> &branchRegions,
                           bool &branching,
                           std::vector<ProteinRegion> &capRegions,
                           std::vector<ProteinRegion> &nucRegions,
                           std::vector<ProteinRegion> &antiCapRegions,
                           std::vector<ProteinRegion> &sevRegions,
                           std::vector<Actin> &actinVec,
                           int &nActin, const bool steric,
                           bool arpPool, ArpGrid &arpGrid,
                           std::vector<ExcZone> &excZones,
                           std::vector<MembraneWall> &memWalls,
                           Cortex &cortex,
                           const bool tether, double currTime,
                           std::ofstream &fusionTime,
                           const bool bDynamics, int &fuseType,
                           std::array<double,2> &fusePoint,
                           std::array<int,3> &fusePointsToAv)
{
    // before
    int numMemPoints = m_numPoints;

    // Create new membrane with points from old membrane
    int numPoints;
    if (fuseType == 0 || fuseType == 1)
    {
        numPoints = fusePointsToAv[1] - fusePointsToAv[0];
    }
    else
    {
        numPoints = fusePointsToAv[2] - fusePointsToAv[1];
    }

    if (numPoints <= 2)
    {
        std::cout << "Fusion with less than 3 points, numPoints: " << numPoints << std::endl;
        return;
    }
    assert(numPoints > 2);
    Membrane phagosome = Membrane(numPoints, temp, fusePoint, *this,
                                          fusePointsToAv, fuseType,
                                          membranes.size());


    for (int i = 0; i < phagosome.getNumPoints(); ++i)
    {
        stericGrid.updateMembrane(phagosome, i);
    }

    membranes.push_back(phagosome);


    // Now need to cut out of the original membrane
    // membrane and phagosome should initially share the fuse point

    // Update steric Grid for all membranes
    for (int i = 0; i < m_numPoints; ++i)
    {
        stericGrid.resetMem(*this, i);
    }

    int start, end;
    if (fuseType == 0 || fuseType == 1)
    {
        start = fusePointsToAv[0];
        end = fusePointsToAv[1];
    }
    else
    {
        start = fusePointsToAv[1];
        end = fusePointsToAv[2];
    }
    int newStart = (start+end)/2;

    //Loop over part of the membrane that will become phagosome
    // Remove any "activated" regions

    // If all targets done then delete all regions
    int regStart = 0;
    int regEnd = m_numPoints;
    if (timesFused+1 != numTars) // timesFused not increased yet!
    {
        // Still fusion events to do
        // loop one either side because of how activated regions have neighbours
        regStart = fusePointsToAv[0]-2;
        regEnd = fusePointsToAv[2]+2;
    }


    for (int i = regStart; i < regEnd; ++i)
    {
        if (m_coupBranchRegIDs[i].size() == 1)
        {
            branchRegions.erase(branchRegions.begin() + m_coupBranchRegIDs[i][0]);
            for (int j = 0; j < m_numPoints; ++j)
            {
                if (m_coupBranchRegIDs[j].size() == 1)
                {
                    if (m_coupBranchRegIDs[j][0] > m_coupBranchRegIDs[i][0])
                    {
                        m_coupBranchRegIDs[j][0] -= 1;
                    }
                }
            }
            m_coupBranchRegIDs[i].clear();
        }

        if (m_coupNucRegIDs[i].size() == 1)
        {
            nucRegions.erase(nucRegions.begin() + m_coupNucRegIDs[i][0]);
            for (int j = 0; j < m_numPoints; ++j)
            {
                if (m_coupNucRegIDs[j].size() == 1)
                {
                    if (m_coupNucRegIDs[j][0] > m_coupNucRegIDs[i][0])
                    {
                        m_coupNucRegIDs[j][0] -= 1;
                    }
                }
            }
            m_coupNucRegIDs[i].clear();
        }

        if (m_coupCapRegIDs[i].size() == 1)
        {
            capRegions.erase(capRegions.begin() + m_coupCapRegIDs[i][0]);
            for (int j = 0; j < m_numPoints; ++j)
            {
                if (m_coupCapRegIDs[j].size() == 1)
                {
                    if (m_coupCapRegIDs[j][0] > m_coupCapRegIDs[i][0])
                    {
                        m_coupCapRegIDs[j][0] -= 1;
                    }
                }
            }
            m_coupCapRegIDs[i].clear();
        }

        if (m_coupAntiCapRegIDs[i].size() == 1)
        {
            antiCapRegions.erase(antiCapRegions.begin() + m_coupAntiCapRegIDs[i][0]);
            for (int j = 0; j < m_numPoints; ++j)
            {
                if (m_coupAntiCapRegIDs[j].size() == 1)
                {
                    if (m_coupAntiCapRegIDs[j][0] > m_coupAntiCapRegIDs[i][0])
                    {
                        m_coupAntiCapRegIDs[j][0] -= 1;
                    }
                }
            }
            m_coupAntiCapRegIDs[i].clear();
        }

        if (m_coupSevRegIDs[i].size() == 1)
        {
            sevRegions.erase(sevRegions.begin() + m_coupSevRegIDs[i][0]);
            for (int j = 0; j < m_numPoints; ++j)
            {
                if (m_coupSevRegIDs[j].size() == 1)
                {
                    if (m_coupSevRegIDs[j][0] > m_coupSevRegIDs[i][0])
                    {
                        m_coupSevRegIDs[j][0] -= 1;
                    }
                }
            }
            m_coupSevRegIDs[i].clear();
        }

    }


    // Move everything down in relevant vectors
    // The membrane loses the number of points in the phagosome plus 2
    // but we only do it numPoints+1 times because we are going to insert fusePoint
    // Remember number of points is the number in the (pre) phagosome
    for (int i = 0; i < numPoints+1; ++i)
    {
        // j starts at fusePointsToAv[0]+1 because we are going to replace
        // fusePointsToAv[0] with the new shared average point (fusePoint)
        for (int j = fusePointsToAv[0]+1; j < m_numPoints-1; ++j)
        {
            m_points[j] = m_points[j+1];
            m_subunit_unitVecs[j] = m_subunit_unitVecs[j+1];
            m_prescribedSubLengths[j] = m_prescribedSubLengths[j+1];
            m_actualSubLengths[j] = m_actualSubLengths[j+1];
            m_stericCells[j] = m_stericCells[j+1];

            m_coupBranchRegIDs[j] = m_coupBranchRegIDs[j+1];
            m_coupNucRegIDs[j] = m_coupNucRegIDs[j+1];
            m_coupCapRegIDs[j] = m_coupCapRegIDs[j+1];
            m_coupAntiCapRegIDs[j] = m_coupAntiCapRegIDs[j+1];
            m_coupSevRegIDs[j] = m_coupSevRegIDs[j+1];
        }
    }

    //Remove numPoints+1
    for (int i = 0; i < numPoints+1; ++i)
    {
        m_points.pop_back();
        m_subunit_unitVecs.pop_back();
        m_prescribedSubLengths.pop_back();
        m_actualSubLengths.pop_back();
        m_stericCells.pop_back();

        m_coupBranchRegIDs.pop_back();
        m_coupNucRegIDs.pop_back();
        m_coupCapRegIDs.pop_back();
        m_coupAntiCapRegIDs.pop_back();
        m_coupSevRegIDs.pop_back();
    }

    // ------------------ add the average point --------------------------------

    m_points[fusePointsToAv[0]] = fusePoint;

    //--------------------------------------------------------------------------
    m_numPoints -= (numPoints+1);

    // Length correction
    m_length -= phagosome.getLength();
    m_segLength = m_length / m_numPoints;

    for (int i = 0; i < m_numPoints; ++i)
    {
        m_prescribedSubLengths[i] = m_segLength;
    }

    assert(m_numPoints + numPoints == numMemPoints - 1);

    calcSubLengths();
    calcSubunitVecs();

    // Reset to membrane

    // Move regions to correspond to new membrane sub

    if (m_coupBranchRegIDs[fusePointsToAv[0]].size()==1)
    {
        int index = m_coupBranchRegIDs[fusePointsToAv[0]][0];
        branchRegions[index].stickToMem(*this);
    }
    if (m_coupNucRegIDs[fusePointsToAv[0]].size()==1)
    {
        int index = m_coupNucRegIDs[fusePointsToAv[0]][0];
        nucRegions[index].stickToMem(*this);
    }
    if (m_coupCapRegIDs[fusePointsToAv[0]].size()==1)
    {
        int index = m_coupCapRegIDs[fusePointsToAv[0]][0];
        capRegions[index].stickToMem(*this);
    }
    if (m_coupAntiCapRegIDs[fusePointsToAv[0]].size()==1)
    {
        int index = m_coupAntiCapRegIDs[fusePointsToAv[0]][0];
        antiCapRegions[index].stickToMem(*this);
    }
    if (m_coupSevRegIDs[fusePointsToAv[0]].size()==1)
    {
        int index = m_coupSevRegIDs[fusePointsToAv[0]][0];
        sevRegions[index].stickToMem(*this);
    }


    // The width of any coupled regions has to be adjusted to fit m_prescribedSubLengths
    // Do for all regions
    for (int i = 0; i < m_numPoints; ++i)
    {

      if (m_coupBranchRegIDs[i].size()==1)
      {
          int index = m_coupBranchRegIDs[i][0];
          branchRegions[index].adjustWidth(m_prescribedSubLengths[i]);
          branchRegions[index].setCoupledSub(i);
      }
      if (m_coupNucRegIDs[i].size()==1)
      {
          int index = m_coupNucRegIDs[i][0];
          nucRegions[index].adjustWidth(m_prescribedSubLengths[i]);
          nucRegions[index].setCoupledSub(i);
      }
      if (m_coupCapRegIDs[i].size()==1)
      {
          int index = m_coupCapRegIDs[i][0];
          capRegions[index].adjustWidth(m_prescribedSubLengths[i]);
          capRegions[index].setCoupledSub(i);
      }
      if (m_coupAntiCapRegIDs[i].size()==1)
      {
          int index = m_coupAntiCapRegIDs[i][0];
          antiCapRegions[index].adjustWidth(m_prescribedSubLengths[i]);
          antiCapRegions[index].setCoupledSub(i);
      }
      if (m_coupSevRegIDs[i].size()==1)
      {
          int index = m_coupSevRegIDs[i][0];
          sevRegions[index].adjustWidth(m_prescribedSubLengths[i]);
          sevRegions[index].setCoupledSub(i);
      }
    }


    m_zWidth =  (m_length/M_PI); // w in calculations
    m_effPersLen = m_bendingMod*m_zWidth / (g_Kb*temp);

    // Update steric Grid for all membranes
    for (int i = 0; i < m_numPoints; ++i)
    {
        stericGrid.updateMembrane(*this, i);
    }

    // Force depol of any filaments which now break steric hindrance
    for (int i = fusePointsToAv[0]-2; i < fusePointsToAv[0]+2; ++i)
    {
        forceDepol(i, actinVec, steric, stericGrid, nActin, tether,
                   arpPool, arpGrid, excZones, memWalls, membranes, cortex,
                   bDynamics, branchRegions, nucRegions, capRegions, sevRegions);
    }


    forceDissociate(actinVec, steric, stericGrid, nActin, tether,
                    arpPool, arpGrid, excZones, memWalls, membranes, cortex,
                    bDynamics, branchRegions, nucRegions, capRegions, sevRegions);

    int phagFuseSub = end - newStart;
    for (int i = phagFuseSub-2; i < phagFuseSub+2; ++i)
    {
        phagosome.forceDepol(i, actinVec, steric, stericGrid, nActin, tether,
                             arpPool, arpGrid, excZones, memWalls, membranes, cortex,
                             bDynamics, branchRegions, nucRegions, capRegions,
                             sevRegions);
    }

    phagosome.forceDissociate(actinVec, steric, stericGrid, nActin, tether,
                              arpPool, arpGrid, excZones, memWalls, membranes,
                              cortex, bDynamics, branchRegions, nucRegions,
                              capRegions, sevRegions);

    std::cout << "Fusion complete" << std::endl;
    timesFused += 1;

    // output the time to file
    fusionTime << currTime << std::endl;


    // Ones associated with phagosome??
    // Turn off branching so it doesn't reduce to global branching routine

    if (timesFused == numTars)
    {
        // All possible fusion events done
        branching = false;
        m_activate = false;
    }

}

void Membrane::forceDepol(int subid, std::vector<Actin> &actinVec,
                          const bool steric, StericGrid &stericGrid,
                          int &nActin, const bool tether,
                          bool arpPool, ArpGrid &arpGrid,
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
     * Given a successful fusion check, there could be filaments that now
     * break steric hindrance, easiest way is to remove the actin causing this
     * Run this function for each new subunit
     */

    // Are any filaments completely outside the cell membrane?
    // or completly inside phagosome (unlikely)
    // do once in separate function

    std::vector<int> occCells; // cells that are occupied
    for (unsigned int j = 0; j < m_stericCells[subid].size(); ++j)
    {
       occCells.push_back(m_stericCells[subid][j]);
    }

    std::vector<int> cellsToCheck = stericGrid.getCellsToCheck(occCells);

    std::array<int,2> checkResult = check_min_dist_Grid_Identify(subid, actinVec, cellsToCheck, stericGrid);
    while (checkResult[0] != -1)
    {
        // A filament has broken steric hindrance, depol it?
        // Which one was it?

        // Is it an end sub?

        if (checkResult[1] == 0 || checkResult[1] == 1)
        {
            // Pointed end
            // Force pointed end depol

            std::array<double,2> barbedEnd;
            barbedEnd[0] = actinVec[checkResult[0]].getBarbedEnd()[0];
            barbedEnd[1] = actinVec[checkResult[0]].getBarbedEnd()[1];

            if (actinVec[checkResult[0]].getParentID() != -1)
            {
                //Filament is a branch, must depol barbed end instead

                depolymerisationBarbedForce(checkResult[0], actinVec, tether, nActin,
                                            arpPool, arpGrid, steric,
                                            stericGrid, excZones, memWalls,
                                            membranes, cortex, bDynamics,
                                            branchRegions, nucRegions,
                                            capRegions, sevRegions);
            }
            else if ((!m_cellMem && checkPointIn(barbedEnd)) || (m_cellMem && !checkPointIn(barbedEnd)))
            {
                // Barbed end inside phagosome or outside cell membrane,
                // must depol barbed end instead

                depolymerisationBarbedForce(checkResult[0], actinVec, tether, nActin,
                                            arpPool, arpGrid, steric,
                                            stericGrid, excZones, memWalls,
                                            membranes, cortex, bDynamics,
                                            branchRegions, nucRegions,
                                            capRegions, sevRegions);
            }
            else
            {
                depolymerisationPointedForce(checkResult[0], actinVec, tether, nActin,
                                            arpPool, arpGrid, steric,
                                            stericGrid, excZones, memWalls,
                                            membranes, cortex, bDynamics,
                                            branchRegions, nucRegions,
                                            capRegions, sevRegions);
            }


        }
        else if (checkResult[1] == actinVec[checkResult[0]].getNumSubs()-1 || checkResult[1] == actinVec[checkResult[0]].getNumSubs()-2)
        {

            // Barbed end
            std::array<double,2> pointedEnd;
            pointedEnd[0] = actinVec[checkResult[0]].getPointedEnd()[0];
            pointedEnd[1] = actinVec[checkResult[0]].getPointedEnd()[1];

            if (actinVec[checkResult[0]].getParentID() == -1 && ((!m_cellMem && checkPointIn(pointedEnd)) || (m_cellMem && !checkPointIn(pointedEnd))))
            {
                // Pointed end inside phagosome, must depol pointed end instead
                // MUST NOT BE A BRANCH THOUGH!
                depolymerisationPointedForce(checkResult[0], actinVec, tether, nActin,
                                            arpPool, arpGrid, steric,
                                            stericGrid, excZones, memWalls,
                                            membranes, cortex, bDynamics,
                                            branchRegions, nucRegions,
                                            capRegions, sevRegions);

            }
            else
            {
              depolymerisationBarbedForce(checkResult[0], actinVec, tether, nActin,
                                          arpPool, arpGrid, steric,
                                          stericGrid, excZones, memWalls,
                                          membranes, cortex, bDynamics,
                                          branchRegions, nucRegions,
                                          capRegions, sevRegions);
            }
        }
        else
        {
            // Middle sub!
            std::cout << "Middle sub! Membrane.cpp line 4709" << std::endl;
            exit(1);
        }

        // Very inefficient!
        checkResult = check_min_dist_Grid_Identify(subid, actinVec, cellsToCheck, stericGrid);
    }

}


void Membrane::forceDissociate(std::vector<Actin> &actinVec,
                              const bool steric, StericGrid &stericGrid,
                              int &nActin, const bool tether,
                              bool arpPool, ArpGrid &arpGrid,
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
  // Called after force depol, if any filaments are completely outside
  // they are deleted (or if they are inside the phagosome)

  std::vector<int> dissociateIDs;
  unsigned int nActinStart = actinVec.size();
  for (unsigned int i = 0; i < nActinStart; ++i)
  {
      std::array<double,2> barbedEnd;
      barbedEnd[0] = actinVec[i].getBarbedEnd()[0];
      barbedEnd[1] = actinVec[i].getBarbedEnd()[1];
      if ((!m_cellMem && checkPointIn(barbedEnd)) || (m_cellMem && !checkPointIn(barbedEnd)))
      {
        // Inside phagosome or outside cell membrane, delete
          int nMonos = actinVec[i].getNumMonomers();
          int smallestNMonos = 3;
          if (actinVec[i].getParentID() != -1)
          {
            // Branch, so is dimer
            smallestNMonos = 2;
          }

          std::cout << "Filament to be deleted: " << i << std::endl;
          for (int j = 0; j < nMonos-smallestNMonos; ++j)
          {
              std::cout << "Depolymerising monomer to trimer/dimer" << j << " Num monos: ";
              std::cout << actinVec[i].getNumMonomers() << std::endl;;
              depolymerisationBarbedForce(i, actinVec, tether, nActin,
                                          arpPool, arpGrid, steric,
                                          stericGrid, excZones, memWalls,
                                          membranes, cortex, bDynamics,
                                          branchRegions, nucRegions,
                                          capRegions, sevRegions);
          }
          dissociateIDs.push_back(actinVec[i].getID());

      }
  }

  dissociationRout(dissociateIDs, actinVec, nActin, steric, stericGrid);

}

void Membrane::initActRegion(std::vector<ProteinRegion> &branchRegions,
                             std::vector<ProteinRegion> &capRegions,
                             std::vector<ProteinRegion> &antiCapRegions,
                             std::vector<ProteinRegion> &nucRegions,
                             std::vector<ProteinRegion> &sevRegions,
                             std::vector<ExcZone> &excZones,
                             const bool activeBranch,
                             const double activeBranchHeight,
                             double arpConc,
                             const bool activeNuc, const double activeNucHeight,
                             const bool activeCap, const double activeCapHeight,
                             const bool activeAntiCap,
                             const double activeAntiCapHeight,
                             const bool activeSever,
                             const double activeSeverHeight)
{
    for (int j = 0; j < m_numPoints; ++j)
    {
        if (calcDistToTarsForActive(excZones, j))
        {
            // membrane should be active

            // if not already a region there, make one

            if (activeBranch && m_coupBranchRegIDs[j].size() == 0)
            {
                std::array<double,2> uVec = { m_subunit_unitVecs[j][0], m_subunit_unitVecs[j][1] };
                std::array<double,2> pointBL = { m_points[j][0], m_points[j][1] };
                branchRegions.push_back((ProteinRegion(uVec, m_prescribedSubLengths[j], activeBranchHeight, pointBL, arpConc, 0, 0, 0, 1, j)));

                branchRegions[branchRegions.size()-1].stickToMem(*this);
                addToBranchReg(j, branchRegions.size()-1);
            }

            if (activeNuc && m_coupNucRegIDs[j].size() == 0)
            {
                std::array<double,2> uVec = { m_subunit_unitVecs[j][0], m_subunit_unitVecs[j][1] };
                std::array<double,2> pointBL = { m_points[j][0], m_points[j][1] };
                nucRegions.push_back((ProteinRegion(uVec, m_prescribedSubLengths[j], activeNucHeight, pointBL, 0, 0, 1.0, 0, 1, j)));

                nucRegions[nucRegions.size()-1].stickToMem(*this);
                addToNucReg(j, nucRegions.size()-1);
            }

            if (activeCap && m_coupCapRegIDs[j].size() == 0)
            {
                std::array<double,2> uVec = { m_subunit_unitVecs[j][0], m_subunit_unitVecs[j][1] };
                std::array<double,2> pointBL = { m_points[j][0], m_points[j][1] };
                capRegions.push_back((ProteinRegion(uVec, m_prescribedSubLengths[j], activeCapHeight, pointBL, 0, 1.0, 0, 0, 1, j)));

                capRegions[capRegions.size()-1].stickToMem(*this);
                addToCapReg(j, capRegions.size()-1);
            }

            if (activeAntiCap && m_coupAntiCapRegIDs[j].size() == 0)
            {
                std::array<double,2> uVec = { m_subunit_unitVecs[j][0], m_subunit_unitVecs[j][1] };
                std::array<double,2> pointBL = { m_points[j][0], m_points[j][1] };
                antiCapRegions.push_back((ProteinRegion(uVec, m_prescribedSubLengths[j], activeAntiCapHeight, pointBL, 0, 0, 0, 0, 1, j)));

                antiCapRegions[antiCapRegions.size()-1].stickToMem(*this);
                addToAntiCapReg(j, antiCapRegions.size()-1);
            }

            if (activeSever && m_coupSevRegIDs[j].size() == 0)
            {
                std::array<double,2> uVec = { m_subunit_unitVecs[j][0], m_subunit_unitVecs[j][1] };
                std::array<double,2> pointBL = { m_points[j][0], m_points[j][1] };
                sevRegions.push_back((ProteinRegion(uVec, m_prescribedSubLengths[j], activeSeverHeight, pointBL, 0, 0, 0, 1.0, 1, j)));

                sevRegions[sevRegions.size()-1].stickToMem(*this);
                addToSevReg(j, sevRegions.size()-1);
            }

        }
    }
}

void Membrane::initCortexRegions(std::vector<ProteinRegion> &branchRegions,
                                 std::vector<ProteinRegion> &nucRegions,
                                 double arpConc)
{
    for (int j = 0; j < m_numPoints; ++j)
    {
        if (m_coupBranchRegIDs[j].size() == 0)
        {
            std::array<double,2> uVec = { m_subunit_unitVecs[j][0], m_subunit_unitVecs[j][1] };
            std::array<double,2> pointBL = { m_points[j][0], m_points[j][1] };
            branchRegions.push_back((ProteinRegion(uVec, m_prescribedSubLengths[j], 200E-9, pointBL, arpConc, 0, 0, 0, 1, j)));

            branchRegions[branchRegions.size()-1].stickToMem(*this);
            addToBranchReg(j, branchRegions.size()-1);
        }

        if (m_coupNucRegIDs[j].size() == 0)
        {
            std::array<double,2> uVec = { m_subunit_unitVecs[j][0], m_subunit_unitVecs[j][1] };
            std::array<double,2> pointBL = { m_points[j][0], m_points[j][1] };
            nucRegions.push_back((ProteinRegion(uVec, m_prescribedSubLengths[j], 200E-9, pointBL, 0, 0, 1.0, 0, 1, j)));

            nucRegions[nucRegions.size()-1].stickToMem(*this);
            addToNucReg(j, nucRegions.size()-1);
        }
    }
}

void Membrane::setCOM()
{
    // Get centre of mass coordinates for membrane

    double xSum = 0;
    double ySum = 0;
    for (int i = 0; i < m_numPoints; ++i)
    {
        xSum += m_points[i][0];
        ySum += m_points[i][1];
    }

    m_COM[0] = xSum/m_numPoints;
    m_COM[1] = ySum / m_numPoints;
}

std::array<double,2> Membrane::calcCOM()
{
    // Get centre of mass coordinates for membrane

    double xSum = 0;
    double ySum = 0;
    for (int i = 0; i < m_numPoints; ++i)
    {
        xSum += m_points[i][0];
        ySum += m_points[i][1];
    }

    std::array<double,2> COM;
    COM[0] = xSum/m_numPoints;
    COM[1] = ySum / m_numPoints;

    return COM;
}

void Membrane::exocytosis(double dt, StericGrid &stericGrid, double temp,
                          double k_exo, double exo_mean, double exo_stdev)
{
  /*
   * Growth of membrane by discrete chunks modelling exocytosis events
   * Stochastic
   * For now let's make WHERE the exocytosis happens random
   */


   if (!m_cellMem)
   {
      // Dont grow phagosome!
      return;
   }


   double randNum { rng.m_probDist(rng.m_mersenne) };
   if (randNum < k_exo*dt)
   {
     // exocytosis event
      auto exoLenDist = std::normal_distribution<double> (exo_mean, exo_stdev);
      double exoLength = exoLenDist(rng.m_mersenne);
      while (exoLength < 0)
      {
          exoLength = exoLenDist(rng.m_mersenne);
      }
      //std::cout << exoLength << std::endl;
      // Where does it happen? Which sub gets bigger
      auto subDist = std::uniform_int_distribution<int> (0, (m_numPoints-1));
      // Choose sub to embiggen
      int sub = subDist(rng.m_mersenne);
      // Allow BrownianDynamics algorithm to handle this by relaxing the constraint
      // on the sub
      m_prescribedSubLengths[sub] += exoLength;
      m_length += exoLength;
      m_zWidth =  (m_length/M_PI); // w in calculations
      m_effPersLen = m_bendingMod*m_zWidth / (g_Kb*temp);
      // Do need to add points at some point
      // Allow a sub to get twice as big and then split it?

      if (m_prescribedSubLengths[sub] >= 2*m_segLength)
      {

          double pointX, pointY;
          if (sub != m_numPoints-1)
          {
             pointX = (m_points[sub][0] + m_points[sub+1][0])/2;
             pointY = (m_points[sub][1] + m_points[sub+1][1])/2;

          }
          else
          {
              pointX = (m_points[sub][0] + m_points[0][0])/2;
              pointY = (m_points[sub][1] + m_points[0][1])/2;
          }

          m_subunit_unitVecs.push_back(std::array<double,2>());
          m_points.push_back(std::array<double,2>());
          m_prescribedSubLengths.push_back(0);
          m_actualSubLengths.push_back(0);

          m_stericCells.push_back(std::vector<int>());

          m_coupBranchRegIDs.push_back(std::vector<int>());
          m_coupNucRegIDs.push_back(std::vector<int>());
          m_coupCapRegIDs.push_back(std::vector<int>());
          m_coupAntiCapRegIDs.push_back(std::vector<int>());
          m_coupSevRegIDs.push_back(std::vector<int>());

          m_numPoints += 1;

          for (int i = m_numPoints-1; i > sub+1; --i)
          {
              m_subunit_unitVecs[i] = m_subunit_unitVecs[i-1];
              m_points[i] = m_points[i-1];
              m_prescribedSubLengths[i] = m_prescribedSubLengths[i-1];
              m_actualSubLengths[i] = m_actualSubLengths[i-1];

              m_coupBranchRegIDs[i] = m_coupBranchRegIDs[i-1];
              m_coupNucRegIDs[i] = m_coupNucRegIDs[i-1];
              m_coupCapRegIDs[i] = m_coupCapRegIDs[i-1];
              m_coupAntiCapRegIDs[i] = m_coupAntiCapRegIDs[i-1];
              m_coupSevRegIDs[i] = m_coupSevRegIDs[i-1];
          }

          m_prescribedSubLengths[sub] /= 2;
          m_prescribedSubLengths[sub+1] = m_prescribedSubLengths[sub];
          m_points[sub+1][0] = pointX;
          m_points[sub+1][1] = pointY;

          calcSubLengths();

          for (int i = 0; i < m_numPoints; ++i)
          {
              stericGrid.resetMemandUpdate(*this, i);
          }

          calcSubunitVecs();

      }

   }
}

Eigen::VectorXd Membrane::calcTetherForces(Cortex &cortex)
{
        // Forces that will be applied to each end point of our membraneFilament
        // Just uses Hookes Law F = -k(r_spring - r_point)
        // Returns 2 doubles per spring, (Fx, Fy)
        // These should then be added to the bending force

        int N = m_numPoints;
        Eigen::VectorXd tetherF (2*N);
        tetherF.setZero();
        if (!cortex.getExist())
          return tetherF;

        int numTethers = m_tetherRelDists.size();
        double restLength = 2E-8 + m_thickness/2 + cortex.getThickness()/2;

        for (int i = 0; i < numTethers; ++i)
        {
                int j = m_tetherSubIDs[i];
                double a = m_tetherRelDists[i];

                assert (a < 1);
                double tetherA_x = m_points[j][0] + m_subunit_unitVecs[j][0]*a*m_prescribedSubLengths[j];
                double tetherA_y = m_points[j][1] + m_subunit_unitVecs[j][1]*a*m_prescribedSubLengths[j];

                int k = cortex.getTetherSubIDs()[i];
                double b = cortex.getTetherDists()[i];

                double tetherB_x = cortex.getPoints()[k][0] + cortex.getUnitVecs()[k][0]*b*cortex.getPresSubLengths()[k];
                double tetherB_y = cortex.getPoints()[k][1] + cortex.getUnitVecs()[k][1]*b*cortex.getPresSubLengths()[k];

                // pointed end forces
                double tetherDist = sqrt((tetherB_x - tetherA_x)*(tetherB_x - tetherA_x) + (tetherB_y - tetherA_y)*(tetherB_y - tetherA_y));
                double forceX = -m_tetherStiff*(1-(restLength/tetherDist))*(tetherA_x - tetherB_x);
                double forceY = -m_tetherStiff*(1-(restLength/tetherDist))*(tetherA_y - tetherB_y);

                tetherF(j*2) = tetherF(j*2) + (1-a)*forceX;
                tetherF((j*2)+1) = tetherF((j*2)+1) + (1-a)*forceY;

                // barbed end forces
                if (j == m_numPoints-1)
                {
                    tetherF(0) = tetherF(0) + a*forceX;
                    tetherF(1) = tetherF(1) + a*forceY;
                }
                else
                {
                    tetherF((j+1)*2) = tetherF((j+1)*2) + a*forceX;
                    tetherF(((j+1)*2)+1) = tetherF(((j+1)*2)+1) + a*forceY;
                }
        }
        return tetherF;
}
