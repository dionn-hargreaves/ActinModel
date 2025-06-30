




#ifndef PROTEINREGION_H
#define PROTEINREGION_H

#include <GTMathematics.h>
#include <random>
#include "Actin.h"
#include "Membrane.h"

typedef gte::FIQuery<double, gte::Segment<2, double>, gte::Segment<2, double>> Find_I;

class ProteinRegion
{
private:
    std::array<double,2> m_uVec;
    std::array<double,2> m_pointBL; // Bottom left point
    double m_width;
    double m_height;
    double m_area;
    double m_normalisedArea;
    double m_arpConc;
    double m_capCoeff;
    double m_nucCoeff;
    double m_sevCoeff;
    bool m_coupledBool;
    int m_coupledMemSub;

    bool m_isCirc;
    bool m_isRing;
    std::array<double,2> m_pointC; // Centre point of circle
    double m_radius;
    double m_innerRadius;

public:
    // Rectangle
    ProteinRegion(std::array<double,2> xUVec, double width,
                  double height, std::array<double,2> pointBL,
                  double arpConc, double capCoeff = 1.0, double nucCoeff = 1.0,
                  double sevCoeff = 1.0, bool coupledToMem = 0, int coupledMemSub = 0);

    // Circle
    ProteinRegion(double radius, std::array<double,2> pointCentre, double arpConc,
                  double capCoeff = 1.0, double nucCoeff = 1.0,
                  double sevCoeff = 1.0, bool coupledToMem = 0, int coupledMemSub = 0);

    // Ring
    ProteinRegion(std::array<double,2> pointCentre, double radius,
                  double innerRadius, double nucCoeff = 1.0,
                  bool coupledToMem = 0, int coupledMemSub = 0);

    void splitRegion(std::vector<ProteinRegion> &regionVec,
                     std::vector<MembraneWall> &memWalls);
                     
    bool checkPointwithin(std::array<double,2> Point) const;
    bool checkPointwithinRect(std::array<double,2> Point) const;
    bool checkPointwithinCirc(std::array<double,2> Point) const;
    bool checkPointwithinRing(std::array<double,2> Point) const;


    std::array<double,2> FindIntersectside(gte::Segment<2,double> actin, gte::Segment<2,double> side) const;
    static int s_chooseRegion(std::vector<ProteinRegion> regions, int nRegions);
    void moveRegionTo(std::array<double,2> newPoint,
                      std::array<double,2> newUVec);


    void moveRegionToRect(std::array<double,2> newBLPoint,
                          std::array<double,2> newUVec);

    void moveRegionToCirc(std::array<double,2> newCPoint);


    void stickToMem(Membrane &membrane);
    void moveRegion(double dx, double dy);
    void changeRegion(double newXmin, double newYmin, double newXmax, double newYmax);

    double getWidth() const { return m_width; }
    double getHeight() const { return m_height; }
    std::array<double,2> getBLPoint() const { return m_pointBL; }
    std::array<double,2> getUVec() const { return m_uVec; }

    double getArea() const { return m_area; }
    double getNormalisedArea() const { return m_normalisedArea; }
    //double getArpCoeff() const { return m_arpCoeff; }
    double getArpConc() const { return m_arpConc; }
    double getCapCoeff() const { return m_capCoeff; }
    double getNucCoeff() const { return m_nucCoeff; }
    double getSevCoeff() const { return m_sevCoeff; }
    bool getCoupledBool() const { return m_coupledBool; }
    int getCoupledSub() const { return m_coupledMemSub; }
    void setCoupledSub(int coupSub) { m_coupledMemSub = coupSub; }
    std::uniform_real_distribution<double> m_rectangleWidth;
    std::uniform_real_distribution<double> m_rectangleHeight;
    std::uniform_real_distribution<double> getWidthDist() const { return m_rectangleWidth; }
    std::uniform_real_distribution<double> getHeightDist() const { return m_rectangleHeight; }

    std::uniform_real_distribution<double> m_RDistro;
    std::uniform_real_distribution<double> m_ThetaDistro;
    std::uniform_real_distribution<double> getRDist() const { return m_RDistro; }
    std::uniform_real_distribution<double> getThetaDist() const { return m_ThetaDistro; }

    void setNormalisedArea(double totalArea) { m_normalisedArea = m_area / totalArea; }

    std::array<double,2> getCentre() const;
    std::array<double,8> getCorners() const;
    bool getCircBool() const { return m_isCirc; }
    bool getRingBool() const { return m_isRing; }
    double getR() const { return m_radius; }
    double getInR() const { return m_innerRadius; }

    void adjustWidth(double newWidth);



};

#endif
