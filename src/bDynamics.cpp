/*
 *  bDynamics.cpp
 *
 *
 *  Function that bends filaments due to thermal fluctuations
 *
 *  James Bradford
 *  University of Sheffield
 *  September 2017
 */


#include "configHeader.h"
#include "globals.h"
#include "RNG.h"
#include "Actin.h"
#include "bDynamics.h"
#include "crosslinking.h"
#include <Eigen/SparseCore>

extern RNG rng;

void BrownianDynamics(std::vector<Actin> &actinVec, const std::vector<ExcZone> &excZones,
                      const std::vector<MembraneWall> &memWalls,
                      const std::vector<Membrane> &membranes,
                      Cortex &cortex, const bool steric,
                      StericGrid &stericGrid,
                      const double dt, double currTime, double temp, double viscosity,
                      bool tethering, bool crossLinking)
{
    // Bead-Rod model - anistrophic drag
    // updated for new points scheme

    std::vector<int> movedStructures;
    for (unsigned int i = 0; i < actinVec.size(); ++i)
    {
        if (actinVec[i].getParentID() == -1 && actinVec[i].getLength() < 2*Actin::s_segmentationLength)
        {
            // straight filament
            if (!tethering || actinVec[i].getNumTethers() < 2)
            {
                brownianMotionRandWalk(actinVec[i], actinVec, viscosity, temp, dt,
                                       excZones, memWalls, membranes, cortex, steric,
                                       stericGrid, tethering, crossLinking);

            }
            continue;
        }
        else if (actinVec[i].getParentID() != -1 && actinVec[i].getLength() < 2*Actin::s_segmentationLength)
        {
            // The filament is a short (straight) branch, that will just move
            // with its mother
            continue;
        }
        else if (actinVec[i].getParentID() == -1 && (actinVec[i].getDaughternum() != 0 || actinVec[i].getNumCLinks() != 0))
        {
            // It the trunk of a branch structure and is flexible
            // Move the whole structure as an ellipse before bending individual branches

            // Check for no master in structure - as it will have already moved as a structure

            if (!crossLinking && (!tethering || actinVec[i].getNumTethers() < 2))
            {
                brownianMotionRandWalk(actinVec[i], actinVec, viscosity, temp, dt,
                                       excZones, memWalls, membranes, cortex, steric,
                                       stericGrid, tethering, crossLinking);
            }
            else if (crossLinking && !actinVec[i].checkForMasterInStructure(actinVec))
            {
                if (std::find(movedStructures.begin(), movedStructures.end(), actinVec[i].getStructureID()) == movedStructures.end())
                {
                        // We have NOT moved this structure before
                        brownianMotionRandWalk(actinVec[i], actinVec, viscosity, temp, dt,
                                               excZones, memWalls, membranes, cortex, steric,
                                               stericGrid, tethering, crossLinking);


                        movedStructures.push_back(actinVec[i].getStructureID());
                        // Have to ensure no master filaments else structure would have moved already
                        // in first if condition
                }
            }
        }

        // N is the number of points (midpoints of the 'rods')
        // N is always the number of subs - 1
        int N = actinVec[i].getNumSubs() - 1;

        Eigen::MatrixXd B (N-1, 2*N);
        B = actinVec[i].getBMatrix(N);
        Eigen::MatrixXd B_T (2*N, N-1);
        B_T = B.transpose();

        Eigen::VectorXd F_bend (2*N);
        F_bend = actinVec[i].calcBendingForces(temp, N);

        Eigen::MatrixXd D (2*N, 2*N);
        D = actinVec[i].buildDMatrixAniso(temp, viscosity, N);

        Eigen::VectorXd F_branch (2*N);
        F_branch = actinVec[i].calcBranchForces(actinVec, N);

        Eigen::VectorXd F_mot (2*N);
        F_mot = actinVec[i].calcMotherForces(actinVec, N);

        Eigen::VectorXd F_tet (2*N);
        F_tet = actinVec[i].calcTetherForces(N, tethering);

        Eigen::VectorXd F_cL (2*N);
        F_cL = actinVec[i].calcCrossLinkForces(N, crossLinking, actinVec);

        Eigen::VectorXd xi (2*N);
        xi = actinVec[i].getRand_dR(dt, temp, viscosity, N);

        Eigen::VectorXd r_old (2*N);
        for (int j = 0; j < (2*N)-1; j+=2)
        {
            // input our points, but we ignore the end points as they are not
            // midpoints of vectors
            int k = j / 2;
            r_old (j) = actinVec[i].getPoints()[k+1][0]; // x
            r_old (j+1) = actinVec[i].getPoints()[k+1][1]; // y
        }


        //Eigen::MatrixXd I (2*N, 2*N);
        //I = Eigen::MatrixXd::Identity(2*N,2*N);

        // Avoid explicitly inverting the matrix B*B_T
        // instead use choleskly factorization
        /*
        Eigen::MatrixXd T (2*N, N-1);
        T = B_T * (B*B_T).inverse();

        // Here is the algorithm
        // d is the vector contained the prescribed bond lengths (between points)
        // we ignore end subunits
        std::vector<double> presLens = actinVec[i].getPresSubLengths();
        Eigen::VectorXd d = Eigen::VectorXd::Map(&presLens[1], presLens.size()-2);

        Eigen::VectorXd r_new (2*N);

        // Equation 17
        //r_new = (I - T*B)*(r_old + (dt/(g_Kb*temp))*D*(F_bend + F_mot + F_branch + F_tet + F_cL) + xi) + T*d;

        // W/O Brownian aspect
        //r_new = (I - T*B)*(r_old + (dt/(g_Kb*temp))*D*(F_bend + F_mot + F_branch + F_tet + F_cL)) + T*d;


        // Avoid matrix-matrix multiplication
        Eigen::VectorXd vector (2*N);
        vector = (r_old + (dt/(g_Kb*temp))*D*(F_bend + F_mot + F_branch + F_tet + F_cL) + xi);
        //r_new = vector - T*(B*vector) + T*d;
        r_new = vector - T*((B*vector) - d);

        */


        // choleskly

        Eigen::VectorXd r_new (2*N);
        std::vector<double> presLens = actinVec[i].getPresSubLengths();
        Eigen::VectorXd d = Eigen::VectorXd::Map(&presLens[1], presLens.size()-2);

        Eigen::VectorXd vector (2*N);
        vector = (r_old + (dt/(g_Kb*temp))*D*(F_bend + F_mot + F_branch + F_tet + F_cL) + xi);

        Eigen::VectorXd BtimesVec (N-1);
        BtimesVec = B*vector;

        //r_new = vector - B_T*((B*B_T).llt().solve(BtimesVec)) + B_T*((B*B_T).llt().solve(d));
        r_new = vector - B_T*((B*B_T).llt().solve(BtimesVec-d));



        //#pragma omp parallel for
        for (int j = 0; j < N; ++j)
        {
            // j is the correct index for r_new
            std::array<double,2> point = { r_new (j*2), r_new (j*2 + 1) };
            // need an index for output
            int k = j + 1;
            actinVec[i].setPoints(k,point);
        }

        // However now we need to recalculate our end points
        actinVec[i].updateUnitVecsNoCheck();

        // save the old ones in case we need

        std::array<double,2> oldPointedEnd, oldBarbedEnd;
        oldPointedEnd[0] = actinVec[i].getPoints()[0][0];
        oldPointedEnd[1] = actinVec[i].getPoints()[0][1];
        oldBarbedEnd[0] = actinVec[i].getPoints()[actinVec[i].getNumSubs()][0];
        oldBarbedEnd[1] = actinVec[i].getPoints()[actinVec[i].getNumSubs()][1];

        // project end points from penultimate points
        std::array<double,2> pointedEnd;
        std::array<double,2> unitVec = actinVec[i].getUnitVec(1);
        double ptEndLength = actinVec[i].getActualSubLengths()[0];
        pointedEnd[0] = actinVec[i].getPoints()[1][0] - ptEndLength*unitVec[0];
        pointedEnd[1] = actinVec[i].getPoints()[1][1] - ptEndLength*unitVec[1];

        std::array<double,2> barbedEnd;
        unitVec = actinVec[i].getUnitVec(actinVec[i].getNumSubs()-2);
        double brEndLength = actinVec[i].getActualSubLengths()[actinVec[i].getNumSubs()-1];

        barbedEnd[0] = actinVec[i].getPoints()[actinVec[i].getNumSubs()-1][0] + brEndLength*unitVec[0];
        barbedEnd[1] = actinVec[i].getPoints()[actinVec[i].getNumSubs()-1][1] + brEndLength*unitVec[1];

        actinVec[i].setPoints(0,pointedEnd);
        actinVec[i].setPoints(actinVec[i].getNumSubs(),barbedEnd);

        actinVec[i].updateSubLengths();
        actinVec[i].updateUnitVecs();

        actinVec[i].saveOldPos(actinVec);

        actinVec[i].rotbranches(actinVec, 0);
        actinVec[i].rotRigidCLinked(actinVec, 0);

        actinVec[i].updateTetherLoc();
        actinVec[i].updateCrossLinkLocs(actinVec);

        // Do steric check here
        // If failed set points with r_old
        if (steric)
        {
            stericGrid.resetAndUpdateAllCLinksAndDaughters(actinVec[i], actinVec);
        }

        // Update grid before the check
        // Change below to a while loop to debug

        if (actinVec[i].check_bend(actinVec, excZones, memWalls,
                                    membranes, cortex, steric, stericGrid))
        {

            // steric check failed
            actinVec[i].setPoints(0, oldPointedEnd);

            actinVec[i].setPoints(N+1, oldBarbedEnd);
            for (int j = 0; j < N; ++j)
            {
                std::array<double,2> point = { r_old (j*2), r_old (j*2 + 1) };
                int k = j+1;
                actinVec[i].setPoints(k,point);
            }

            actinVec[i].updateUnitVecs();

            actinVec[i].updateSubLengths();
            actinVec[i].resetToOldPos(actinVec);


            actinVec[i].updateTetherLoc();

            actinVec[i].updateCrossLinkLocs(actinVec);

            // Update back
            if (steric)
            {
                stericGrid.resetAndUpdateAllCLinksAndDaughters(actinVec[i], actinVec);
            }

        }
        else
        {
            // accepted
            actinVec[i].setChangedBoolTrue();
        }

    }

}

void brownianMotionRandWalk(Actin &actin, std::vector<Actin> &actinVec, const double viscosity,
                            const double temperature, const double dt,
                            const std::vector<ExcZone> &excZones,
                            const std::vector<MembraneWall> &memWalls,
                            const std::vector<Membrane> &membranes,
                            Cortex &cortex,
                            const bool steric, StericGrid &stericGrid, bool tethering,
                            bool crossLinking)
{
    // bundled translational BM and rotational BM into one overarching functions
    // more in line with actin class

    // Always called on a trunk

    if (actin.getDaughternum() == 0 && actin.getNumCLinks() == 0)
    {

        // Unbranched, we use Kirkwood model
        std::array<double, 3> D_array = actin.kirkwoodDiff(viscosity,
                                                           (g_Kb*temperature),
                                                           true,
                                                           true);


        // Rotation
        double xyangle_rms = sqrt(2*dt*D_array[2]);
        std::normal_distribution<double> xy_rotDiff(0, xyangle_rms);
        double d_xyangle { xy_rotDiff(rng.m_mersenne) };

        if (!tethering || (tethering && actin.getNumTethers() == 0))
        {
            actin.rotate2D(d_xyangle, actinVec, excZones, memWalls, membranes,
                           cortex, steric, stericGrid, 0);
        }
        else if (tethering && actin.getNumTethers() != 0)
        {
            actin.rotate2DTether(d_xyangle, actinVec, excZones, memWalls, membranes,
                                 cortex, steric, stericGrid, 0);
        }


        if (!tethering || (tethering && actin.getNumTethers() == 0))
        {
            // Translation
            double x_diffusion_rmsDist = sqrt(2*dt*D_array[0]);
            double y_diffusion_rmsDist = sqrt(2*dt*D_array[1]);
            std::normal_distribution<double> x_diff_gauss(0,x_diffusion_rmsDist);
            std::normal_distribution<double> y_diff_gauss(0,y_diffusion_rmsDist);
            // If the filament is an original parent
            double dxprime { x_diff_gauss(rng.m_mersenne) };
            double dyprime { y_diff_gauss(rng.m_mersenne) };

            Eigen::Matrix3d rotationMatrix;
            actin.diffuse2D(dxprime, dyprime, actinVec, rotationMatrix,
                            excZones, memWalls, membranes, cortex, steric,
                            stericGrid, 0, crossLinking);
        }

    }
    else if (actin.getNumCLinks() > 0 && actin.getDaughternum() == 0 && actin.checkForMasterInStructure(actinVec))
    {
        // Unbranched with at least one crosslink
        // Rotation

        if (actin.getNumCLinks() == 1)
        {
            //std::cout << "CaseB 1 CLink" << std::endl;

            std::array<double, 3> D_array2 = actin.kirkwoodDiff(viscosity,
                                                               (g_Kb*temperature),
                                                               true,
                                                               true);

            double xyangle_rms = sqrt(2*dt*D_array2[2]);
            std::normal_distribution<double> xy_rotDiff(0, xyangle_rms);
            double d_xyangle { xy_rotDiff(rng.m_mersenne) };

            actin.rotate2DCLink(d_xyangle, actinVec, excZones, memWalls, membranes,
                                cortex, steric, stericGrid, 0);

            // rotate around each other
        }

        if (actin.getCLinkMaster()) // only do it once
        {
            // Make sure other filament is not a short branch?

            // Ellipse
            // Check if the structure has changed
            bool structChange = false;
            for (unsigned int i = 0; i < actinVec.size(); ++i)
            {
                if (actin.getStructureID() == actinVec[i].getStructureID())
                {
                    if (actinVec[i].getChanged())
                    {
                        // structure has changed
                        structChange = true;
                        break;
                    }
                }
            }


            Eigen::Matrix3d rotationMatrix;
            std::array<double,3> centrePoint;
            std::array<double,3> D_array;

            if (structChange)
            {
                D_array = actin.fitEllipsoids(actinVec,
                                                 viscosity,
                                                 (g_Kb*temperature),
                                                 rotationMatrix,
                                                 centrePoint,
                                                 true,
                                                 true);


                actin.setDarray(D_array);
                actin.setRotMat(rotationMatrix);
                actin.setEllipCentre(centrePoint);

                // ReSet changed bool to false
                for (unsigned int i = 0; i < actinVec.size(); ++i)
                {
                    if (actin.getStructureID() == actinVec[i].getStructureID())
                    {
                        actinVec[i].setChangedBoolFalse();
                    }
                }
            }
            else
            {
                D_array = actin.getDarray();
                rotationMatrix = actin.getRotMat();
                centrePoint = actin.getEllipCentre();
            }

            assert(D_array[2] > 0);
            // Rotation
            double xyangle_rms = sqrt(2*dt*D_array[2]);

            std::normal_distribution<double> xy_rotDiff(0, xyangle_rms);
            double d_xyangle { xy_rotDiff(rng.m_mersenne) };

            actin.rotateEllipsoidALL(d_xyangle, actinVec, rotationMatrix,
                                           centrePoint, excZones, memWalls,
                                           membranes, cortex,
                                           steric, stericGrid, 0);
            // rotate around each other
            actin.updateCrossLinkLocs(actinVec);
            actin.setRotMat(rotationMatrix);
            //}

            // Translation as ellipse
            double x_diffusion_rmsDist = sqrt(2*dt*D_array[0]);
            double y_diffusion_rmsDist = sqrt(2*dt*D_array[1]);
            assert(D_array[0] > 0);
            assert(D_array[1] > 0);

            std::normal_distribution<double> x_diff_gauss(0,x_diffusion_rmsDist);
            std::normal_distribution<double> y_diff_gauss(0,y_diffusion_rmsDist);
            // If the filament is an original parent
            // Prime is in ref frame of ellipse
            double dxprime { x_diff_gauss(rng.m_mersenne) };
            double dyprime { y_diff_gauss(rng.m_mersenne) };
            actin.diffuse2D(dxprime, dyprime, actinVec, rotationMatrix,
                            excZones, memWalls, membranes, cortex, steric,
                            stericGrid, 0, crossLinking);

            actin.updateCrossLinkLocs(actinVec);
            //std::cout << "CaseB Master translate" << std::endl;

        }

    }
    else if (actin.getDaughternum() != 0 && actin.countCrossLinksInCLStruct(actinVec) > 0 && actin.checkForMasterInStructure(actinVec))
    {
        // As above but branched

        if (actin.countCrossLinksInCLStruct(actinVec) == 1)
        {
          // small ellipse rotates around crosslink point

          Eigen::Matrix3d rotationMatrixCL;
          std::array<double,3> centrePointCL;
          std::array<double,3>  D_arrayCL = actin.fitEllipsoidsCLink(actinVec,
                                                                    viscosity,
                                                           (g_Kb*temperature),
                                                             rotationMatrixCL,
                                                                centrePointCL,
                                                                         true,
                                                                        true);

          double xyangle_rmsCL = sqrt(2*dt*D_arrayCL[2]);
          std::normal_distribution<double> xy_rotDiffCL(0, xyangle_rmsCL);
          double d_xyangleCL { xy_rotDiffCL(rng.m_mersenne) };
          // rotate around the one crosslink in the structure
          actin.rotateMiniEllipsoid(d_xyangleCL, actinVec, rotationMatrixCL,
                                       centrePointCL, excZones, memWalls,
                                       membranes, cortex,
                                       steric, stericGrid, 0);
            // rotate around each other
          actin.updateCrossLinkLocs(actinVec);

        }


        if (actin.getCLinkMaster())
        {
            assert(actin.getParentID() == -1);
            // Check if the structure has changed
            bool structChange = false;
            for (unsigned int i = 0; i < actinVec.size(); ++i)
            {
                if (actin.getStructureID() == actinVec[i].getStructureID())
                {
                    if (actinVec[i].getChanged())
                    {
                        // structure has changed
                        structChange = true;
                        break;
                    }
                }
            }

            Eigen::Matrix3d rotationMatrix;
            std::array<double,3> centrePoint;
            std::array<double,3> D_array;

            if (structChange)
            {

                D_array = actin.fitEllipsoids(actinVec,
                                                 viscosity,
                                                 (g_Kb*temperature),
                                                 rotationMatrix,
                                                 centrePoint,
                                                 true,
                                                 true);


                actin.setDarray(D_array);
                actin.setRotMat(rotationMatrix);
                actin.setEllipCentre(centrePoint);

                // Set changed bool to false
                for (unsigned int i = 0; i < actinVec.size(); ++i)
                {
                    if (actin.getStructureID() == actinVec[i].getStructureID())
                    {
                        actinVec[i].setChangedBoolFalse();
                    }
                }

            }
            else
            {

                D_array = actin.getDarray();
                rotationMatrix = actin.getRotMat();
                centrePoint = actin.getEllipCentre();
            }

            assert(D_array[2] > 0);
            double xyangle_rms = sqrt(2*dt*D_array[2]);
            std::normal_distribution<double> xy_rotDiff(0, xyangle_rms);
            double d_xyangle { xy_rotDiff(rng.m_mersenne) };
            actin.rotateEllipsoidALL(d_xyangle, actinVec, rotationMatrix,
                                    centrePoint, excZones, memWalls,
                                    membranes, cortex, steric, stericGrid, 0);
              // rotate around each other
            actin.updateCrossLinkLocs(actinVec);

            // Translation as ellipse
            assert(D_array[0] > 0);
            assert(D_array[1] > 0);
            double x_diffusion_rmsDist = sqrt(2*dt*D_array[0]);
            double y_diffusion_rmsDist = sqrt(2*dt*D_array[1]);
            std::normal_distribution<double> x_diff_gauss(0,x_diffusion_rmsDist);
            std::normal_distribution<double> y_diff_gauss(0,y_diffusion_rmsDist);
            // If the filament is an original parent
            // Prime is in ref frame of ellipse
            double dxprime { x_diff_gauss(rng.m_mersenne) };
            double dyprime { y_diff_gauss(rng.m_mersenne) };

            actin.diffuse2D(dxprime, dyprime, actinVec, rotationMatrix,
                            excZones, memWalls, membranes, cortex, steric,
                            stericGrid, 0, crossLinking);

            actin.updateCrossLinkLocs(actinVec);


            actin.setRotMat(rotationMatrix);

        }

    }
    else
    {
        // branched
        // either no crosslinks ON THE TRUNK
        // or there is no master in WHOLE STRUCTURE

        // Check if the structure has changed
        bool structChange = false;
        for (unsigned int i = 0; i < actinVec.size(); ++i)
        {
            if (actin.getStructureID() == actinVec[i].getStructureID())
            {
                if (actinVec[i].getChanged())
                {
                    // structure has changed
                    structChange = true;
                    break;
                }
            }
        }


        Eigen::Matrix3d rotationMatrix;
        std::array<double,3> centrePoint;
        std::array<double,3> D_array;

        if (structChange)
        {

            D_array = actin.fitEllipsoids(actinVec,
                                             viscosity,
                                             (g_Kb*temperature),
                                             rotationMatrix,
                                             centrePoint,
                                             true,
                                             true);


            actin.setDarray(D_array);
            actin.setRotMat(rotationMatrix);
            actin.setEllipCentre(centrePoint);

            // Set changed bool to false
            for (unsigned int i = 0; i < actinVec.size(); ++i)
            {
                if (actin.getStructureID() == actinVec[i].getStructureID())
                {
                    actinVec[i].setChangedBoolFalse();
                }
            }
        }
        else
        {

            D_array = actin.getDarray();
            rotationMatrix = actin.getRotMat();
            centrePoint = actin.getEllipCentre();
        }

        // Rotation
        assert(D_array[2] > 0);
        double xyangle_rms = sqrt(2*dt*D_array[2]);
        std::normal_distribution<double> xy_rotDiff(0, xyangle_rms);
        double d_xyangle { xy_rotDiff(rng.m_mersenne) };

        if (actin.getNumTethers() != 0 && actin.getNumCLinks() == 0)
        {
            // ellipse rotates around tether point
            actin.rotateEllipsoid2DTether(d_xyangle, actinVec, rotationMatrix,
                                          excZones, memWalls, membranes, cortex,
                                          steric, stericGrid, 0);


            actin.updateCrossLinkLocs(actinVec);
        }
        else if (actin.getNumTethers() == 0)
        {
            actin.rotateEllipsoidALL(d_xyangle, actinVec, rotationMatrix,
                                     centrePoint, excZones, memWalls,
                                     membranes, cortex, steric, stericGrid, 0);
              // rotate around each other
            actin.updateCrossLinkLocs(actinVec);
        }

        actin.setRotMat(rotationMatrix);

        if (!tethering || actin.getNumTethers() == 0)
        {
            // Translation as ellipse
            assert(D_array[0] > 0);
            assert(D_array[1] > 0);
            double x_diffusion_rmsDist = sqrt(2*dt*D_array[0]);
            double y_diffusion_rmsDist = sqrt(2*dt*D_array[1]);
            std::normal_distribution<double> x_diff_gauss(0,x_diffusion_rmsDist);
            std::normal_distribution<double> y_diff_gauss(0,y_diffusion_rmsDist);
            // If the filament is an original parent
            // Prime is in ref frame of ellipse
            double dxprime { x_diff_gauss(rng.m_mersenne) };
            double dyprime { y_diff_gauss(rng.m_mersenne) };
            actin.diffuse2D(dxprime, dyprime, actinVec, rotationMatrix,
                            excZones, memWalls, membranes, cortex, steric,
                            stericGrid, 0, crossLinking);


            actin.updateCrossLinkLocs(actinVec);

        }
    }
}

void initConfiguration(std::vector<Actin> &actinVec, const std::vector<ExcZone> &excZones,
                       const std::vector<MembraneWall> &memWalls,
                       const std::vector<Membrane> &membranes, Cortex &cortex,
                       const bool steric, StericGrid &stericGrid)
{
    // Function that initialises a bent configuration for all the actin filaments
    // based on the static distribution - see Shiladetya's paper
    // Changed for our new points scheme (original is below)

    for (unsigned int i = 0; i < actinVec.size(); ++i)
    {
        if (actinVec[i].getLength() <= 3*Actin::s_segmentationLength)
        {
            continue;
        }
        else
        {
            actinVec[i].setFlexibleBool(true);
        }

        int centreSub;
        actinVec[i].findCentrePoint(centreSub);


        for (int j = centreSub+1; j < actinVec[i].getNumSubs()-1; ++j)
        {
            // Bending upwards from centre
            double stdDev_static = sqrt(actinVec[i].getActualSubLengths()[j] / Actin::s_persistenceLength);
            if (j == actinVec[i].getNumSubs()-2)
            {
                // Last proper subunit must take into account extra length
                stdDev_static = sqrt((actinVec[i].getActualSubLengths()[j]+actinVec[i].getActualSubLengths()[j+1]) / Actin::s_persistenceLength);
            }

            std::normal_distribution<double> bendDist(0,stdDev_static);
            double angle = bendDist(rng.m_mersenne);
            actinVec[i].changeToNewAngle(j, centreSub, angle);
        }
        for (int j = centreSub-1; j >= 1; --j)
        {
            double stdDev_static = sqrt(actinVec[i].getActualSubLengths()[j] / Actin::s_persistenceLength);
            if (j == 1)
            {
                // Last proper subunit must take into account extra length
                stdDev_static = sqrt((actinVec[i].getActualSubLengths()[1]+actinVec[i].getActualSubLengths()[0]) / Actin::s_persistenceLength);
            }

            std::normal_distribution<double> bendDist(0,stdDev_static);
            double angle = bendDist(rng.m_mersenne);
            actinVec[i].changeToNewAngle(j, centreSub, angle);

        }

        actinVec[i].updatePointsUp(centreSub+1);
        actinVec[i].updatePointsDown(centreSub);
        actinVec[i].updateUnitVecs();

        // Update grid before check
        if (steric)
        {
            for (int j = 0; j < actinVec[i].getNumSubs(); ++j)
            {
                stericGrid.resetSubandUpdate(actinVec[i], j);
            }
        }

        if (actinVec[i].check_bend(actinVec, excZones, memWalls, membranes,
                                    cortex, steric, stericGrid))
        {
            // initialisation failed, end simulation
            std::cout << "Bent config failed, either have shorter filaments or a larger nucleation region.\n";
            std::cout << "Or considering adding code to reposition filaments that are hindered, although this might be inefficient!\n";
            std::cout << "(Relevant code is in bDynamics.cpp, lines 960ish)\n";
            std::cout << "Ending simulation" << std::endl;
            exit(1);
        }

    }
}
