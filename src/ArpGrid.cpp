/*
 *  ArpGrid.cpp
 *
 *  C++ file containing the definition of the ArpGrid class.
 *  This contains member variable initialisations and non-trivial member
 *  functions.
 *  The declaration of this class (including the definitions of the trivial
 *  member functions) is found in the associated header file "ArpGrid.h"
 *
 *  James Bradford
 *  University of Sheffield
 *  Nov 2018
 */

// Constructor for axis-aligned rectangle, could allow rotation in the future
// Let's stick with this

#include "configHeader.h"
#include "geometry.h"
#include "ArpGrid.h"
#include "Actin.h"
#include "globals.h"
#include "RNG.h"

extern RNG rng;

typedef gte::TIQuery<double, gte::Segment<2, double>, gte::Segment<2, double>> Check_I;
typedef gte::FIQuery<double, gte::Segment<2, double>, gte::Segment<2, double>> Find_I;

ArpGrid::ArpGrid()
{}

ArpGrid::ArpGrid(double x_min, double y_min, double x_max, double y_max,
                       double cellWidth, double cellHeight, double cellDepth,
                       double arpConc, bool flow)
    : m_x_min { x_min },
    m_y_min { y_min },
    m_x_max { x_max },
    m_y_max { y_max },
    m_x_width { x_max - x_min },
    m_y_height { y_max - y_min },
    m_cellWidth { cellWidth },
    m_cellHeight { cellHeight },
    m_cellDepth { cellDepth },
    m_flow { flow }
{
    if (m_cellWidth > m_x_width || m_cellHeight > m_y_height)
    {
        std::cout << "Arp grid: The cell size is greater than the whole grid.";
        std::cout << " Check config file. Exiting simulation" << std::endl;
        exit(1);
    }

    if (m_cellWidth != m_cellHeight && m_flow)
    {
        std::cout << "Arp grid: Cells need to be squares for there to be diffusion";
        std::cout << " between them. Check config file" << std::endl;
        exit(1);
    }

    m_grids_x = m_x_width / m_cellWidth;
    m_grids_y = m_y_height / m_cellHeight;

    m_conc_vec.resize(m_grids_x);  // resize top level vector
    m_num_vec.resize(m_grids_x);  // resize top level vector

    // Conversion factor
    // from conc to number
    m_convert = m_cellHeight*m_cellWidth*m_cellDepth*1E-3*g_NA;
    // so N = C*m_convert

    for (int i = 0; i < m_grids_x; ++i)
    {
        m_conc_vec[i].resize(m_grids_y);  // resize each contained vector
        m_num_vec[i].resize(m_grids_y);

        for (int j = 0; j < m_grids_y; ++j)
        {
            m_conc_vec[i][j] = arpConc;
            m_num_vec[i][j] = round(arpConc * m_convert);
        }
    }

}

std::array<int, 2> ArpGrid::findGrid(double xpos, double ypos)
{
    //function that given a point, finds the grid square it is in
    // relative to minimum so must be positive!
    double x_rel = xpos - m_x_min;
    double y_rel = ypos - m_y_min;

    if (x_rel < 0 || y_rel < 0 || x_rel > m_x_width || y_rel > m_y_height)
    {
        std::array<int, 2> outside;
        outside[0] = -1;
        outside[1] = -1;
        return outside;
    }

    std::array<int, 2> gridnum;
    gridnum[0] = (int)(x_rel / m_cellWidth);
    gridnum[1] = (int)(y_rel / m_cellHeight);

    return gridnum;
}

void ArpGrid::decrement(int cellIDx, int cellIDy)
{
    // Due to branching decrement the arp2/3 in the grid
    if (cellIDx == -1)
    {
        std::cout << "Error: line 123 in ArpGrid.cpp. Edge effect, system is";
        std::cout << "gaining monomers from outside the grid" << std::endl;
        exit(1);
    }
    m_num_vec[cellIDx][cellIDy] -= 1;
    m_conc_vec[cellIDx][cellIDy] =  m_num_vec[cellIDx][cellIDy] / m_convert;
}

void ArpGrid::increment(double xpos, double ypos)
{
    // Due to debranching increment the arp2/3 in the grid
    // turn this off
    //return;
    std::array<int,2> gridnum;
    gridnum = findGrid(xpos, ypos);
    if (gridnum[0] == -1)
    {
        gridnum = findNearestCell(xpos,ypos);
    }

    m_num_vec[gridnum[0]][gridnum[1]] += 1;
    m_conc_vec[gridnum[0]][gridnum[1]] =  m_num_vec[gridnum[0]][gridnum[1]] / m_convert;
}


double ArpGrid::getConcentration(int xcoord, int ycoord) const
{
    return m_conc_vec[xcoord][ycoord];
}

int ArpGrid::getNumber(int xcoord, int ycoord) const
{
    return m_num_vec[xcoord][ycoord];
}

std::array<int,2> ArpGrid::findNearestCell(double xpos, double ypos)
{
    // Given that we have depolarismation outside the grid we may want to conserve
    // monomers, so they are not lost. Current plan is to return the monomer to
    // the nearest cell immediately.

    // This function finds the nearest cell to the position (that is outside the
    // grid)

    std::array<double, 2> outsidePoint = { xpos, ypos };
    std::vector<double> sqDistances;
    // Absolute values for corrdinates of cell edges. Start at 0,0 (bottom left)
    double xstart = m_x_min;
    double ystart = m_y_min;

    // First do distance check with cell 0

    std::array<double,2> centrePoint = { xstart+m_cellWidth/2, ystart+m_cellHeight/2 };
    double minSqDist = distanceBetPoints2DSQR(outsidePoint, centrePoint);
    std::array<int, 2> nearestCell = { 0, 0 };

    for (int i = 0; i < m_grids_x; ++i)
    {
        for (int j = 0; j < m_grids_y; ++j)
        {
            if (i == 0 && j == 0)
            {
                // Cell 0, already checked this
                ystart += m_cellHeight;
                continue;
            }

            std::array<double,2> centrePoint = { xstart+m_cellWidth/2, ystart+m_cellHeight/2 };
            double sqDist = distanceBetPoints2DSQR(outsidePoint, centrePoint);
            if (sqDist < minSqDist)
            {
                // Overwrite nearest cell
                minSqDist = sqDist;
                nearestCell = { i, j };
            }

            ystart += m_cellHeight;
        }

        xstart += m_cellWidth;
    }

    return nearestCell;
}

void ArpGrid::calcFlow(double dt)
{
    if (m_grids_y == 1 || m_grids_x == 1)
        return;
    double D = 5E-12; // diffusion constant of g-actin
    std::vector< std::vector<double> > conc_vec_TMP;
    conc_vec_TMP.resize(m_grids_x);

    double FO = (D * dt) / (m_cellHeight*m_cellWidth);

    if (FO > 0.25)
    {
        std::cout << "Unstable monomer diffusion" << std::endl;
        std::cout << "Need to make timestep smaller" << std::endl;
        std::cout << "Or cell sizes larger" << std::endl;
        exit(1);
    }
    for (int i = 0; i < m_grids_x; ++i)
    {
        conc_vec_TMP[i].resize(m_grids_y);

        for (int j = 0; j < m_grids_y; ++j)
        {
            if (i == 0)
            {
                // left surface of grid
                if (j == 0)
                {
                    // bottom left corner of grid
                    conc_vec_TMP[i][j] = FO * (2*m_conc_vec[i+1][j]
                                                 + 2*m_conc_vec[i][j+1])
                                           + (1-4*FO)*m_conc_vec[i][j];

                }
                else if (j == m_grids_y -1)
                {
                    // top left corner of grid
                    conc_vec_TMP[i][j] = FO * (2*m_conc_vec[i+1][j]
                                                 + 2*m_conc_vec[i][j-1])
                                           + (1-4*FO)*m_conc_vec[i][j];
                }
                else
                {
                    // left side
                    conc_vec_TMP[i][j] = FO * (2*m_conc_vec[i+1][j]
                                                 + m_conc_vec[i][j+1]
                                                 + m_conc_vec[i][j-1])
                                           + (1-4*FO)*m_conc_vec[i][j];
                }
            }
            else if (i == m_grids_x -1)
            {
                // right surface of grid
                if (j == 0)
                {
                    // bottom right corner of grid
                    conc_vec_TMP[i][j] = FO * (2*m_conc_vec[i-1][j]
                                                 + 2*m_conc_vec[i][j+1])
                                           + (1-4*FO)*m_conc_vec[i][j];

                }
                else if (j == m_grids_y -1)
                {
                    // top right corner of grid
                    conc_vec_TMP[i][j] = FO * (2*m_conc_vec[i-1][j]
                                                 + 2*m_conc_vec[i][j-1])
                                           + (1-4*FO)*m_conc_vec[i][j];
                }
                else
                {
                    // right side
                    conc_vec_TMP[i][j] = FO * (2*m_conc_vec[i-1][j]
                                                 + m_conc_vec[i][j+1]
                                                 + m_conc_vec[i][j-1])
                                           + (1-4*FO)*m_conc_vec[i][j];
                }
            }
            else if (j == 0)
            {
                // bottom side
                conc_vec_TMP[i][j] = FO * (m_conc_vec[i+1][j]
                                             + m_conc_vec[i-1][j]
                                             + 2*m_conc_vec[i][j+1])
                                       + (1-4*FO)*m_conc_vec[i][j];
            }
            else if (j == m_grids_y -1)
            {
                // top side
                conc_vec_TMP[i][j] = FO * (m_conc_vec[i+1][j]
                                             + m_conc_vec[i-1][j]
                                             + 2*m_conc_vec[i][j-1])
                                       + (1-4*FO)*m_conc_vec[i][j];
            }
            else
            {
                // interior

                conc_vec_TMP[i][j] = FO * (m_conc_vec[i+1][j]
                                             + m_conc_vec[i-1][j]
                                             + m_conc_vec[i][j+1]
                                             + m_conc_vec[i][j-1])
                                       + (1-4*FO)*m_conc_vec[i][j];
            }
        }
    }

    for (int i = 0; i < m_grids_x; ++i)
    {
        for (int j = 0; j < m_grids_y; ++j)
        {
            m_conc_vec[i][j] = conc_vec_TMP[i][j];
            m_num_vec[i][j] = round(m_conc_vec[i][j] * m_convert);
        }
    }
}
