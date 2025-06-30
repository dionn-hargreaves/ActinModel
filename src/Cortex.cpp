/*
 *  Cortex.cpp
 *
 *  C++ file containing the definition of the Cortex class.
 *  This is an inextensible 1d object similar to the membrane, made up of bead-rod components
 *  This contains member variable initialisations and non-trivial member
 *  functions.
 *  The declaration of this class (including the definitions of the trivial
 *  member functions) is found in the associated header file "Cortex.h"
 *
 *  James Bradford
 *  University of Sheffield
 *  Nov 2019
 */

#include "configHeader.h"
#include "Actin.h"
#include "ArpGrid.h"
#include "MembraneWall.h"
#include "Cortex.h"
#include "globals.h"
#include "geometry.h"
#include "polymerisation.h"
#include <Eigen/SparseCore>
#include<Eigen/IterativeLinearSolvers>

extern RNG rng;
typedef gte::DCPQuery<double, gte::Segment<2, double>, gte::Segment<2, double>> SegDistQuery; // distance between two segments
typedef gte::TIQuery<double, gte::Segment<2, double>, gte::Segment<2, double>> Check_I; // check intersection between two lines

Cortex::Cortex()
    : m_exist { false }
{}
// Cortex has the membrane it is attahed to passed to it
Cortex::Cortex(double temp, Membrane &membrane)
    : m_id { 0 },
    m_bendingMod { 50*g_Kb*temp },
    m_segLength { 2.5E-7 },
    m_thickness { 2E-7 }
{
    m_tetherStiff = 1E-7; // stiffness of attachment to membrane
    double radius = membrane.getRadius() - membrane.getThickness()/2 - m_thickness/2 - 2E-8; // have 20nm gap between edge of cortex and membrane (reported in literature)
    m_ypos = membrane.getYpos();

    m_length = radius*2*M_PI;
    m_exist = true;
    // Partition the cortex into points
    m_numPoints = ceil(m_length / m_segLength);
    m_segLength = m_length / m_numPoints;

    m_points.push_back({+m_segLength/2, m_ypos-radius});
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
    }

    calcSubunitVecs();


    m_zWidth =  2*radius; // w in calculations
    m_effPersLen = m_bendingMod*m_zWidth / (g_Kb*temp); // do we need this?

    setCOM();

    createTetherPointsCortex(membrane);
}

Cortex::~Cortex(){
    //std::cout << "Cortex correctly deleted" << std::endl;
}

std::vector<double> Cortex::initCircle()
{
    /*
    Function that initialise the Cortex filament in a (near)perfect circle
    */
    double angleDiff = (2*M_PI) / m_numPoints;

    std::vector<double> tmpAngles;
    for (int i = 0; i < m_numPoints-2; ++i)
    {
        tmpAngles.push_back(angleDiff);

    }
    return tmpAngles;
}

void Cortex::calcSubunitVecs()
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

void Cortex::calcSubLengths()
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

Eigen::VectorXd Cortex::calcBendingForces(double temp)
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

Eigen::MatrixXd Cortex::buildDMatrix(double temp, double viscosity)
{
    int N = m_numPoints;

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

Eigen::VectorXd Cortex::getRand_dR(double dt, Eigen::MatrixXd D)
{
    // xi 2N vector

    int N = m_numPoints;
    Eigen::VectorXd rand_dR (2*N);
    for (int i = 0; i < 2*N; ++i)
    {
        std::normal_distribution<double> noise(0,sqrt(2*D(i,i)*dt)); // it is a dR
        rand_dR (i) = noise(rng.m_mersenne);
    }

    return rand_dR;
}

Eigen::MatrixXd Cortex::getBMatrix()
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

void Cortex::fluctuate(double temp, double viscosity, double dt,
                             std::vector<Actin> &actinvec,
                             std::vector<ExcZone> &excZones,
                             const bool steric,
                             StericGrid &stericGrid, double currTime,
                             std::vector<ProteinRegion> &branchRegions,
                             bool &branching,
                             std::vector<ProteinRegion> &capRegions,
                             std::vector<ProteinRegion> &antiCapRegions,
                             std::vector<ProteinRegion> &nucRegions,
                             std::vector<Membrane> &membranes,
                             int &nActin, const bool tether,
                             bool arpPool, ArpGrid &arpGrid,
                             std::vector<MembraneWall> &memWalls,
                             std::ofstream &fusionTime, const bool bDynamics)
{
    /*
    Function that is to be called every timestep to fluctuate
    the Cortex

    Same algorithm used for dynamics of our actin filaments

    */

    // N is the number of points (midpoints of the 'rods')
    int N = m_numPoints;

    Eigen::MatrixXd B (N, 2*N);

    B = getBMatrix();

    Eigen::VectorXd F_bend (2*N);
    F_bend = calcBendingForces(temp);

    Eigen::VectorXd F_tet (2*N);
    F_tet = calcTetherForces(membranes[0]);

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


    //Eigen::MatrixXd I (2*N, 2*N);
    //I = Eigen::MatrixXd::Identity(2*N,2*N);

    Eigen::MatrixXd T (2*N, N);
    // Following line is cpu intensive due to inverse, Wang and Gao suggest
    // approximating T

    Eigen::MatrixXd B_T (2*N, N);
    B_T = B.transpose();
    /*
    // Avoid inverting matrix!
    T = B_T * (B*B_T).inverse();

    // Here is the algorithm
    // d is the vector contained the prescribed bond lengths (between points)


    Eigen::VectorXd d = Eigen::VectorXd::Map(&m_prescribedSubLengths[0], m_prescribedSubLengths.size());

    Eigen::VectorXd r_new (2*N);

    // Modify to change forces

    r_new = (I - T*B)*(r_old + (dt/(g_Kb*temp))*D*(F_bend+F_tet) + xi) + T*d;
    */

    Eigen::VectorXd d = Eigen::VectorXd::Map(&m_prescribedSubLengths[0], m_prescribedSubLengths.size());
    Eigen::VectorXd vector (2*N);
    vector = (r_old + (dt/(g_Kb*temp))*D*(F_bend+F_tet) + xi);

    Eigen::VectorXd BtimesVec (N);
    BtimesVec = B*vector;
    Eigen::VectorXd r_new (2*N);

    r_new = vector - B_T*((B*B_T).llt().solve(BtimesVec-d));


    for (int j = 0; j < N; ++j)
    {
        // j is the correct index for r_new
        m_points[j][0] = r_new (j*2);
        m_points[j][1] = r_new (j*2 + 1);
        calcSubLengths();
        calcSubunitVecs();

        if (steric)
        {
            // Update steric grid before check
            if (j == 0)
            {
                stericGrid.resetCortexandUpdate(*this, 0);
                stericGrid.resetCortexandUpdate(*this, m_numPoints-1);
            }
            else
            {
                stericGrid.resetCortexandUpdate(*this, j-1);
                stericGrid.resetCortexandUpdate(*this, j);
            }
        }

        if (steric && (checkExVolSUBGrid(actinvec, j, stericGrid) || ExcZone::s_checkExVolCortex(*this, excZones, j) || checkExVolSUBGridMembrane(membranes, j, stericGrid) ))
        {
            // This move has violated excvol and therefore we reject it
            m_points[j][0] = r_old (j*2);
            m_points[j][1] = r_old (j*2 + 1);
            calcSubLengths();
            calcSubunitVecs();
            // Update steric grid back
            if (j == 0)
            {
                stericGrid.resetCortexandUpdate(*this, 0);
                stericGrid.resetCortexandUpdate(*this, m_numPoints-1);
            }
            else
            {
                stericGrid.resetCortexandUpdate(*this, j-1);
                stericGrid.resetCortexandUpdate(*this, j);
            }
        }
    }
}

bool Cortex::checkExVolSUBGrid(const std::vector<Actin> &actinvec,
                                         int pointID, StericGrid &stericGrid)
{
    /*
    Function that should be called after fluctuating the Cortex, to see if
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
    for (int i = 0; i < 2; ++i)
    {
        int p = pointsToCheck[i];
        std::vector<int> occCells; // cells that are occupied

        for (unsigned int j = 0; j < m_stericCells[p].size(); ++j)
        {
            occCells.push_back(m_stericCells[p][j]);
        }
        std::vector<int> cellsToCheck;
        cellsToCheck = stericGrid.getCellsToCheckBIG(occCells, m_thickness);





        if (check_min_dist_Grid(p, actinvec, cellsToCheck, stericGrid))
            return true;
    }

    return false;

}

bool Cortex::checkExVolSUBGridMembrane(const std::vector<Membrane> &membranes,
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

    gte::Segment<2,double> cortex;
    for (int i = 0; i < 2; ++i)
    {
        if (pointsToCheck[i] == m_numPoints-1)
        {
            // its the end sub
            cortex.p[0][0] = m_points[0][0];
            cortex.p[1][0] = m_points[m_numPoints-1][0];
            cortex.p[0][1] = m_points[0][1];
            cortex.p[1][1] = m_points[m_numPoints-1][1];
        }
        else
        {
            cortex.p[0][0] = m_points[pointID][0];
            cortex.p[1][0] = m_points[pointID+1][0];
            cortex.p[0][1] = m_points[pointID][1];
            cortex.p[1][1] = m_points[pointID+1][1];
        }


        // To avoid repeated checks against the same subunits that may appear in
        // more than one cell, have a recorded of which ones we have checked
        std::vector<std::array<int,2>> checkedSubs;

        for (unsigned int j = 0; j < membranes.size(); ++j)
        {
            for (int k = 0; k < membranes[j].getNumPoints(); ++k)
            {
                gte::Segment<2,double> membrane;

                membrane.p[0][0] = membranes[j].getPoints()[k][0];
                membrane.p[0][1] = membranes[j].getPoints()[k][1];

                if (k == membranes[j].getNumPoints() -1 )
                {
                    // end sub
                    membrane.p[1][0] = membranes[j].getPoints()[0][0];
                    membrane.p[1][1] = membranes[j].getPoints()[0][1];
                }
                else
                {
                    membrane.p[1][0] = membranes[j].getPoints()[k+1][0];
                    membrane.p[1][1] = membranes[j].getPoints()[k+1][1];
                }

                SegDistQuery min_d_GTE;
                auto result = min_d_GTE(membrane, cortex);
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

bool Cortex::check_min_dist_Grid(int subid, const std::vector<Actin> &actinvec,
                                              std::vector<int> &cellsToCheck,
                                              StericGrid &stericGrid)
{
    gte::Segment<2,double> cortex;

    if (subid == m_numPoints-1)
    {
        // its the end sub
        cortex.p[0][0] = m_points[0][0];
        cortex.p[1][0] = m_points[m_numPoints-1][0];
        cortex.p[0][1] = m_points[0][1];
        cortex.p[1][1] = m_points[m_numPoints-1][1];
    }
    else
    {
        cortex.p[0][0] = m_points[subid][0];
        cortex.p[1][0] = m_points[subid+1][0];
        cortex.p[0][1] = m_points[subid][1];
        cortex.p[1][1] = m_points[subid+1][1];
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
            auto result = min_d_GTE(cortex, filament);
            double distance = result.distance;

            if (distance < (m_thickness/2 + actinvec[actinID].getStericRadius()))
            {
                return true;
            }
        }
    }
    return false;
}

std::array<int,2> Cortex::check_min_dist_Grid_Identify(int subid, const std::vector<Actin> &actinvec,
                                              std::vector<int> &cellsToCheck,
                                              StericGrid &stericGrid)
{
    /*
     * Same as above but returns the actinID and actinSubID of the filament sub
     * that breaks Steric hindrance. If steric hindrance is not broken
     * returns {-1,-1}
     */

    gte::Segment<2,double> cortex;

    if (subid == m_numPoints-1)
    {
        // its the end sub
        cortex.p[0][0] = m_points[0][0];
        cortex.p[1][0] = m_points[m_numPoints-1][0];
        cortex.p[0][1] = m_points[0][1];
        cortex.p[1][1] = m_points[m_numPoints-1][1];
    }
    else
    {
        cortex.p[0][0] = m_points[subid][0];
        cortex.p[1][0] = m_points[subid+1][0];
        cortex.p[0][1] = m_points[subid][1];
        cortex.p[1][1] = m_points[subid+1][1];
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
            auto result = min_d_GTE(cortex, filament);
            double distance = result.distance;

            if (distance < (m_thickness/2 + actinvec[actinID].getStericRadius()))
            {
                return {actinID, actinSubID};
            }
        }
    }
    return {-1, -1};
}

bool Cortex::s_checkExVolBarbPoly(const Actin &filament, Cortex &cortex,
                                  StericGrid &stericGrid)
{
    if (cortex.getExist() && cortex.checkExVolBarbPolyGrid(filament, stericGrid))
        return true;


    return false;
}

bool Cortex::checkExVolBarbPoly(const Actin &filament) const
{
    /*
    Function that should be called during barbed end polymerisation of actin filaments

    This checks the new monomer against the whole Cortex

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


bool Cortex::checkExVolBarbPolyGrid(const Actin &filament,
                                              StericGrid &stericGrid) const
{
    /*
    Function that should be called during barbed end polymerisation of actin filaments

    This checks the new monomer against the whole Cortex

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
    cellsToCheck = stericGrid.getCellsToCheckBIG(barbCellIDs, m_thickness);


    for (unsigned int i = 0; i < cellsToCheck.size(); ++i)
    {
        std::vector<std::array<int,2> > cellContents = stericGrid.getCellContentsCortex(cellsToCheck[i]);
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
bool Cortex::s_checkExVolPointPoly(const Actin &filament, Cortex &cortex,
                                   StericGrid &stericGrid)
{

    if (cortex.getExist() && cortex.checkExVolPointPolyGrid(filament, stericGrid))
        return true;


    return false;

}

bool Cortex::checkExVolPointPoly(const Actin &filament) const
{
    /*
    Function that should be called during pointed end polymerisation of actin filaments

    This checks the new monomer against the whole Cortex

    Returns true if exc vol is violated
    */

    // Define our actin monomer
    std::array<double,2> unitVec = filament.getUnitVec(0);
    double newPointx = filament.getPointedEnd()[0] - (Actin::s_monomerLength * unitVec[0]);
    double newPointy = filament.getPointedEnd()[1] - (Actin::s_monomerLength * unitVec[1]);
    gte::Segment<2,double> monomer;
    monomer.p[0][0] = filament.getPointedEnd()[0];
    monomer.p[1][0] = newPointx;
    monomer.p[0][1] = filament.getPointedEnd()[1];
    monomer.p[1][1] = newPointy;


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

bool Cortex::checkExVolPointPolyGrid(const Actin &filament,
                                               StericGrid &stericGrid) const
{
    // Define our actin monomer

    gte::Segment<2,double> firstSub;
    firstSub.p[0][0] = filament.getPointedEnd()[0];
    firstSub.p[1][0] = filament.getPrePointedEnd()[0];
    firstSub.p[0][1] = filament.getPointedEnd()[1];
    firstSub.p[1][1] = filament.getPrePointedEnd()[1];

    std::vector<int> PointCellIDs = filament.getStericCells(0);
    std::vector<int> cellsToCheck;
    cellsToCheck = stericGrid.getCellsToCheckBIG(PointCellIDs, m_thickness);


    for (unsigned int i = 0; i < cellsToCheck.size(); ++i)
    {
        std::vector<std::array<int,2>> cellContents = stericGrid.getCellContentsCortex(cellsToCheck[i]);
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

bool Cortex::checkExVolNuc(const Actin &actin) const
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

bool Cortex::s_checkExVolNuc(const Actin &actin, Cortex &cortex,
                                   std::vector<int> cellsToCheck,
                                   StericGrid &stericGrid)
{

    if (cortex.getExist() && cortex.checkExVolNuc_Grid(actin, cellsToCheck, stericGrid))
        return true;


    return false;

}

bool Cortex::checkExVolNuc_Grid(const Actin &actin,
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

    // Need to check that it is outside the cortex!
    if (!checkInOut(actin))
    {
        return true;
    }

    std::vector<int> occCells; // cells that are occupied

    for (int i = 0; i < actin.getNumSubs(); ++i)
    {
            for (unsigned int j = 0; j < actin.getStericCells(i).size(); ++j)
            {
                    if (std::find(occCells.begin(), occCells.end(),actin.getStericCells(i)[j])==occCells.end()) // if it is not in occCells already
                    {
                            occCells.push_back(actin.getStericCells(i)[j]);
                    }
            }
    }
    cellsToCheck = stericGrid.getCellsToCheckBIG(occCells, m_thickness);


    for (unsigned int i = 0; i < cellsToCheck.size(); ++i)
    {

        std::vector<std::array<int,2>> cellContents = stericGrid.getCellContentsCortex(cellsToCheck[i]);

        for (unsigned int j = 0; j < cellContents.size(); ++j)
        {
            // j is membrane sub

            if (std::find(checkedSubs.begin(), checkedSubs.end(), cellContents[j][1]) != checkedSubs.end())
            {
                // We have checked this mem sub before!
                continue;
            }
            checkedSubs.push_back(cellContents[j][1]);

            //std::cout << "Does it ever get here?" << std::endl;
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

bool Cortex::checkExVolBM(const Actin &actin) const
{
    /*
    Function that calculates the distance between all of the filament and the bent
    Cortex. For now let's just have this as checking for an intersection
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

bool Cortex::s_checkExVolBM(int subid, const Actin &actin, Cortex &cortex,
                                  std::vector<int> cellsToCheck,
                                  StericGrid &stericGrid)
{

    if (cortex.getExist() && cortex.check_min_dist_Actin_Grid(subid, actin, cellsToCheck, stericGrid))
    {
        return true;
    }


    return false;
}


bool Cortex::check_min_dist_Actin_Grid(int subid, const Actin &actin,
                                                 std::vector<int> cellsToCheck,
                                                 StericGrid &stericGrid) const
{
    /*
    Function that calculates the distance between all of the filament and the bent
    Cortex. For now let's just have this as checking for an intersection
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

        std::vector<std::array<int,2>> cellContents = stericGrid.getCellContentsCortex(cellsToCheck[i]);

        for (unsigned int j = 0; j < cellContents.size(); ++j)
        {
            // j is membrane sub

            if (std::find(checkedSubs.begin(), checkedSubs.end(), cellContents[j][1]) != checkedSubs.end())
            {
                // We have checked this mem sub before!
                continue;
            }
            checkedSubs.push_back(cellContents[j][1]);

            //std::cout << "Does it ever get here?" << std::endl;
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

bool Cortex::checkInOut(const Actin &actin) const
{
    /*
    Uses the crossing point algorithm to determine if a given point
    (so the centre) of our actin filament, lies in or out of the cortex shell
    If it is in return false if out return true
    */


    // Define a point we know with certainity will be outside of our cell,
    // here it is at (1m, 0m)
    std::array<double,2> outPoint = {1,0};
    int junk; // not needed
    std::array<double,2> actinP = actin.findCentrePoint(junk);

    // Use geometry tools to find intersections
    // define our "ray"

    gte::Segment<2,double> ray;

    ray.p[0][0] = outPoint[0];
    ray.p[1][0] = actinP[0];
    ray.p[0][1] = outPoint[1];
    ray.p[1][1] = actinP[1];

    // Now we need to check this against our whole cortex and count intersections

    gte::Segment<2,double> cortex;
    int numIntersections = 0;
    for (int i = 0; i < m_numPoints; ++i)
    {
        // The cortex
        int np = i+1;
        if (i == m_numPoints-1)
            np = 0;
        cortex.p[0][0] = m_points[i][0];
        cortex.p[1][0] = m_points[np][0];
        cortex.p[0][1] = m_points[i][1];
        cortex.p[1][1] = m_points[np][1];

        Check_I IntersectQuery;
        auto result = IntersectQuery(ray, cortex);
        numIntersections += result.numIntersections;
    }


    if (numIntersections % 2 == 0)
    {
        // even, so outside
        return true;
    }
    else
    {
        return false;
    }

}

bool Cortex::checkPointIn(std::array<double,2> point) const
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

void Cortex::setCOM()
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

std::array<double,2> Cortex::calcCOM()
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

void Cortex::exocytosis(double dt, StericGrid &stericGrid, double temp)
{
  /*
   * Growth of membrane by discrete chunks modelling exocytosis events
   * Stochastic
   * For now let's make WHERE the exocytosis happens random
   */

   double k_exo = 0.01; // have this as input to function. units of per second
   // Should this be a density? so rate per um or per um^2? Bigger cells=more events?
   k_exo = 5;


   double randNum { rng.m_probDist(rng.m_mersenne) };
   if (randNum < k_exo*dt)
   {
     // exocytosis event
      auto exoLenDist = std::normal_distribution<double> (220E-9, 45E-9);
      double exoLength = exoLenDist(rng.m_mersenne);
      while (exoLength < 0)
      {
          exoLength = exoLenDist(rng.m_mersenne);
      }
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
      // Allow a sub to get twice as big and then split it

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

          m_numPoints += 1;

          for (int i = m_numPoints-1; i > sub+1; --i)
          {
              m_subunit_unitVecs[i] = m_subunit_unitVecs[i-1];
              m_points[i] = m_points[i-1];
              m_prescribedSubLengths[i] = m_prescribedSubLengths[i-1];
              m_actualSubLengths[i] = m_actualSubLengths[i-1];

          }

          m_prescribedSubLengths[sub] /= 2;
          m_prescribedSubLengths[sub+1] = m_prescribedSubLengths[sub];
          m_actualSubLengths[sub] /= 2;
          m_actualSubLengths[sub+1] = m_actualSubLengths[sub];

          m_points[sub+1][0] = pointX;
          m_points[sub+1][1] = pointY;

          for (int i = 0; i < m_numPoints; ++i)
          {
              stericGrid.resetCortexandUpdate(*this, i);
          }

          calcSubunitVecs();

      }

   }
}
void Cortex::createTetherPointsCortex(Membrane &membrane)
{
  // Called on the cortex in its constructor

  // For now have a tether point at each point on cortex object

  std::vector<double> tempAnglesMem;
  for (int i = 0; i < membrane.getNumPoints(); ++i)
  {
      double theta = atan2(membrane.getPoints()[i][1]-m_ypos, membrane.getPoints()[i][0]);
      if (theta < 0)
          theta += 2*M_PI;

      tempAnglesMem.push_back(theta);
  }

  std::vector<int> tetherSubIDs;
  std::vector<double> tetherRelDists;
  for (int i = 0; i < m_numPoints; ++i)
  {
      m_tetherSubIDs.push_back(i);
      m_tetherRelDists.push_back(0);

      double theta = atan2(m_points[i][1]-m_ypos, m_points[i][0]);
      if (theta < 0)
          theta += 2*M_PI;

      for (int j = 0; j < membrane.getNumPoints(); ++j)
      {
          int np = j + 1; // next point
          if (j == membrane.getNumPoints()-1)
          {
              np = 0;
          }

          if (theta > tempAnglesMem[j] && theta < tempAnglesMem[np])
          {
              int subID = j;
              double relDistAlongSub = (theta-tempAnglesMem[j]) / (tempAnglesMem[np] - tempAnglesMem[j]);
              tetherSubIDs.push_back(subID);
              tetherRelDists.push_back(relDistAlongSub);
              break;
          }
          else if (tempAnglesMem[j] > tempAnglesMem[np]) // cross positive x axis
          {
              int subID = j;
              double distToX = 2*M_PI - tempAnglesMem[j];
              double distFromX = tempAnglesMem[np];

              if (theta > tempAnglesMem[j])
              {
                  double relDistAlongSub = (theta-tempAnglesMem[j]) / (distToX + distFromX);
                  tetherSubIDs.push_back(subID);
                  tetherRelDists.push_back(relDistAlongSub);
                  break;
              }
              else if (theta < tempAnglesMem[np])
              {
                  double relDistAlongSub = 1 - (tempAnglesMem[np] - theta) / (distToX + distFromX);
                  tetherSubIDs.push_back(subID);
                  tetherRelDists.push_back(relDistAlongSub);
                  break;
              }

            }

          }
      }

      membrane.setTetherSubIDs(tetherSubIDs);
      membrane.setTetherDists(tetherRelDists);
}

Eigen::VectorXd Cortex::calcTetherForces(Membrane &membrane)
{
        // Forces that will be applied to each end point of our membraneFilament
        // Just uses Hookes Law F = -k(r_spring - r_point)
        // Returns 2 doubles per spring, (Fx, Fy)
        // These should then be added to the bending force

        int N = m_numPoints;
        Eigen::VectorXd tetherF (2*N+2);
        tetherF.setZero();

        int numTethers = m_tetherRelDists.size();
        double restLength = 2E-8 + m_thickness/2 + membrane.getThickness()/2;

        for (int i = 0; i < numTethers; ++i)
        {
                int j = m_tetherSubIDs[i];
                double a = m_tetherRelDists[i];

                assert (a < 1);
                double tetherA_x = m_points[j][0] + m_subunit_unitVecs[j][0]*a*m_prescribedSubLengths[j];
                double tetherA_y = m_points[j][1] + m_subunit_unitVecs[j][1]*a*m_prescribedSubLengths[j];

                int k = membrane.getTetherSubIDs()[i];
                double b = membrane.getTetherDists()[i];

                double tetherB_x = membrane.getPoints()[k][0] + membrane.getUnitVecs()[k][0]*b*membrane.getPresSubLengths()[k];
                double tetherB_y = membrane.getPoints()[k][1] + membrane.getUnitVecs()[k][1]*b*membrane.getPresSubLengths()[k];

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

void Cortex::moveCortex(std::array<double,2> dR, const bool steric,
                       StericGrid &stericGrid)
{
    /*
      Function that simply moves the cortex for reference frame adjustment
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
            stericGrid.resetCortexandUpdate(*this, i);
         }
     }
}
