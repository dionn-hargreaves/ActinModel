




#ifndef ARPGRID_H
#define ARPGRID_H

typedef gte::FIQuery<double, gte::Segment<2, double>, gte::Segment<2, double>> Find_I;

class ArpGrid
{
private:
    double m_x_min;
    double m_y_min;
    double m_x_max;
    double m_y_max;
    double m_x_width;
    double m_y_height;
    double m_cellWidth;
    double m_cellHeight;
    double m_cellDepth;
    double m_convert;
    int m_grids_y;
    int m_grids_x;
    bool m_flow;
    std::vector< std::vector<double> > m_conc_vec;
    std::vector< std::vector<int>  > m_num_vec;
public:
    ArpGrid();

    ArpGrid(double x_min, double y_min, double x_max, double y_max,
               double cellWidth, double cellHeight, double cellDepth,
               double arpConc, bool flow = false);


    std::array<int, 2> findGrid(double xpos, double ypos);
    void decrement(int cellIDx, int cellIDy);
    void increment(double xpos, double ypos);
    double getConcentration(int xcoord, int ycoord) const;
    int getNumber(int xcoord, int ycoord) const;
    std::array<int,2> findNearestCell(double xpos, double ypos);
    void calcFlow(double dt);

    bool getFlow() const { return m_flow; }
    int getNumCellsX() const { return m_grids_x; }
    int getNumCellsY() const { return m_grids_y; }

    double getCellWidth() const { return m_cellWidth; }
    double getCellHeight() const { return m_cellHeight; }
    double getCellDepth() const { return m_cellDepth; }









};

#endif
