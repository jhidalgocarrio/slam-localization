#define BOOST_TEST_MODULE template_for_test_test
#include <boost/test/included/unit_test.hpp>


/** Eigen **/
#include <Eigen/Core> /** Core */
#include <Eigen/StdVector> /** For STL container with Eigen types **/


BOOST_AUTO_TEST_CASE( EIGEN_MATRICES )
{
    Eigen::Matrix<double, Eigen::Dynamic, 1> X(4, 1);
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> A(4, 4);
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> B(4, 4);

    X << 3, 5, 4, 8;
    A << 1, 2, 3, 4,
         1, 2, 3, 4,
         1, 2, 3, 4,
         1, 2, 3, 4;

    B << 5, 6, 7, 8,
         5, 6, 7, 8,
         5, 6, 7, 8,
         5, 6, 7, 8;

    std::cout<<"X\n"<<X<<"\n";
    std::cout<<"A\n"<<A<<"\n";
    std::cout<<"B\n"<<B<<"\n";
    std::cout<<"X.transpose() * A * X: "<<X.transpose() * A * X <<"\n";
    const double value = (X.block(0, 0, 2, 1).transpose() * A.block(0, 0, 2, 2) * X.block(0, 0, 2, 1))[0];
    std::cout<<"X.transpose() * A * X: "<<value<<"\n";

}
