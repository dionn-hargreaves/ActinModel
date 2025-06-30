#include "configHeader.h"
#include "ellipseconfidence.h"


std::array<double,3>  ellipsoidfit(const std::vector< std::array<double,3> > &data,
                                   int dimensions, Eigen::Matrix3d &rotationMatrix,
                                   std::array<double,3> &centrePoint)
{

    int N_points { static_cast<int>(data.size()) };
    double chisquare_val { };
    if (dimensions == 3)
        chisquare_val = sqrt(4.642); // 80% confidence limit = 4.642
    else if (dimensions == 2)
        chisquare_val = sqrt(3.22); // 80% confidence limit = 3.22, 68% = 2.30


    Eigen::MatrixXd mat (N_points, 3);

    for (int i = 0; i < N_points; ++i)
    {
        mat.row(i) = Eigen::VectorXd::Map(&data[i][0], data[i].size());

    }

    Eigen::MatrixXd centered = mat.rowwise() - mat.colwise().mean();

    Eigen::MatrixXd cov = (centered.adjoint() * centered) / (N_points - 1);

    Eigen::EigenSolver<Eigen::MatrixXd> es(cov);
    // ellipsoid_vector needs to be the largest eigenvector
    Eigen::Vector3d eigenvalues = es.eigenvalues().real();

    std::array< std::pair< double, Eigen::Vector3d >, 3 > eigens;

    for (int i = 0; i < 3; ++i)
        eigens[i] = std::make_pair( eigenvalues(i), es.eigenvectors().real().col(i) );


    std::sort(eigens.begin(), eigens.end(), [](auto &left, auto &right) {
        return left.first > right.first;
    });

    for (int i = 0; i < 3; ++i)
    {
      eigenvalues (i) = eigens[i].first;
    }

    Eigen::Vector3d ellipsoid_vector = eigens[0].second; // maximum eigenvector


    double angle = atan2(ellipsoid_vector(1), ellipsoid_vector(0)); // angle about z axis
    //double angle2 = atan2(ellipsoid_vector(2), ellipsoid_vector(0)); // angle about y axis
    //double angle3 = atan2(ellipsoid_vector(2), ellipsoid_vector(1)); // angle about the x axis

    // When in 2D these angles should always be 0
    double angle2 = 0;
    double angle3 = 0;

    //Center points - only needed for debugging

    // Find average x, y, z
    double cenx { 0 };
    double ceny { 0 };
    //double cenz { 0 };

    for (int i = 0; i < N_points; ++i)
    {
        cenx += mat (i, 0);
        ceny += mat (i, 1);
        //cenz += mat (i, 2);
    }
    cenx /= N_points;
    ceny /= N_points;
    //cenz /= N_points;


    double a = chisquare_val*sqrt(eigenvalues(0));
    double b = chisquare_val*sqrt(eigenvalues(1));
    double c = chisquare_val*sqrt(eigenvalues(2));


    Eigen::AngleAxisd rollAngle(angle, Eigen::Vector3d::UnitZ());
    Eigen::AngleAxisd yawAngle(angle2, Eigen::Vector3d::UnitY());
    Eigen::AngleAxisd pitchAngle(angle3, Eigen::Vector3d::UnitX());

    Eigen::Quaternion<double> q = rollAngle * yawAngle * pitchAngle;

    rotationMatrix = q.matrix();


    std::array<double,3> radii { a, b, c };
    centrePoint[0] = cenx;
    centrePoint[1] = ceny;
    centrePoint[2] = 0; // make 0 for 2d simulation

    return radii;
}
