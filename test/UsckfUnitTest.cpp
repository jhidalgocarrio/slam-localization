#define BOOST_TEST_MODULE template_for_test_test
#include <boost/test/included/unit_test.hpp>

/** Library **/
#include <rover_localization/filters/Usckf.hpp> /** USCKF class with Manifolds */
#include <rover_localization/filters/MtkWrap.hpp> /** USCKF wrapper for the state vector */
#include <rover_localization/filters/State.hpp> /** Filters State */
#include <rover_localization/Configuration.hpp> /** Constant values of the library */

/** Eigen **/
#include <Eigen/Core> /** Core */
#include <Eigen/StdVector> /** For STL container with Eigen types **/

/** Standard libs **/
#include <iostream>
#include <vector>

std::vector<int> hola(100, 2);
typedef std::vector< MTK::vect<3, double> > MTKintFeature;

MTKintFeature adios(100, Eigen::Vector3d::Zero() );

/** Wrap the VectorState **/
typedef localization::MtkWrap<localization::VectorState> WVectorState;
typedef localization::MtkWrap<localization::SingleState> WSingleState;

//acc is acceleration without any perturbation
//angvelo is angular velocity (gyros) without any perturbation
WSingleState processModel (const WSingleState &serror, const Eigen::Matrix<double, 3, 1> &acc,
                        const Eigen::Matrix<double, 3, 1> &angvelo,
                        const Eigen::Quaternion<double> &orientq, double dt)
{
    WSingleState s2;
    Eigen::Matrix <double,3,3> dFki; /** Discrete form of subsytem models */

    //std::cout<<"[PROCESS_MODEL] acc:\n"<<acc<<"\n";
    //std::cout<<"[PROCESS_MODEL] angevelo:\n"<<angvelo<<"\n";

    /** Propagate the position (discretization approx) **/
    dFki = Eigen::Matrix<double, 3, 3>::Identity() * dt +
        0.5 * Eigen::Matrix<double, 3, 3>::Identity() * Eigen::Matrix<double, 3, 3>::Identity() * pow(dt,2);
    s2.pos = serror.pos + dFki * serror.vel;

    /** Propagate the velocity (discretization approx) **/
    Eigen::Matrix <double,3, 1> deltaVel;
    deltaVel = /*orientq.inverse() */ (acc * dt); /** TO-DO: Increment in velocity in world frame */
    s2.vel = serror.vel +  serror.orient * deltaVel - orientq.inverse() * (serror.abias * dt);

    /** Propagate the error quaternion **/
    Eigen::Matrix <double,3, 1> deltaAngle = (angvelo - 0.5 * serror.gbias) * dt;
    MTK::SO3 <double> rot = MTK::SO3<double>::exp (deltaAngle);
    s2.orient = (serror.orient * rot);

    /** The bias (gyros and acc) update **/
    s2.gbias = serror.gbias;
    s2.abias = serror.abias;

    return s2;
}

localization::Usckf <WVectorState, WSingleState>::SingleStateCovariance processNoiseCov (double dt)
{
    localization::Usckf <WVectorState, WSingleState>::SingleStateCovariance cov =
        localization::Usckf <WVectorState, WSingleState>::SingleStateCovariance::Zero();
    MTK::setDiagonal (cov, &WSingleState::pos, 0.1 * dt);
    MTK::setDiagonal (cov, &WSingleState::vel,  0.1 * dt);
    MTK::setDiagonal (cov, &WSingleState::orient, 0.1 * dt);
    MTK::setDiagonal (cov, &WSingleState::gbias, 0.002 * dt);
    MTK::setDiagonal (cov, &WSingleState::abias, 0.002 * dt);

    return cov ;
}

BOOST_AUTO_TEST_CASE( USCKF )
{

    localization::MTK_FEATURE_TYPE( vec3 ) featuresVO(24, 1.34*Eigen::Vector3d::Identity() );
    localization::MTK_FEATURE_TYPE( vec3 ) featuresICP(24, 3.45*Eigen::Vector3d::Identity() );

    WVectorState vstate;
    WVectorState verror;

    std::cout<<"vstate::DOF is "<<vstate.DOF<<"\n";
    std::cout<<"vstate.statek is "<<vstate.statek<<"\n";
    std::cout<<"vstate.statek_l is "<<vstate.statek_l<<"\n";
    std::cout<<"vstate.statek_i is "<<vstate.statek_i<<"\n";
    std::cout<<"vstate: "<<vstate<<"\n";

    //std::vector< Eigen::Matrix <double, 3, 1> , Eigen::aligned_allocator < Eigen::Matrix <double, 3, 1> > > featuresVO;

    featuresVO = featuresICP;
    featuresVO.resize(5);
    std::cout<<"size of featuresk: "<<vstate.featuresk.size()<<"\n";
    vstate.featuresk.resize(featuresVO.size());
    vstate.featuresk = featuresVO;
    std::cout<<"size of featuresk: "<<vstate.featuresk.size()<<"\n";
    vstate.featuresk_l.resize(featuresICP.size());
    vstate.featuresk_l = featuresICP;
    std::cout<<"size of featuresk_l: "<<vstate.featuresk_l.size()<<"\n";
    std::cout<<"vstate: "<<vstate<<"\n";

    const double dt = 0.01; /** 100Hz */

    /** Initial covariance matrix **/
    localization::Usckf<WVectorState, WSingleState>::VectorStateCovariance P0;
    MTK::subblock (P0, &WVectorState::statek, &WVectorState::statek) = 0.0025 * localization::Usckf<WVectorState, WSingleState>::SingleStateCovariance::Identity();
    MTK::subblock (P0, &WVectorState::statek_l, &WVectorState::statek_l) = 0.0035 * localization::Usckf<WVectorState, WSingleState>::SingleStateCovariance::Identity();
    MTK::subblock (P0, &WVectorState::statek_i, &WVectorState::statek_i) = 0.0045 * localization::Usckf<WVectorState, WSingleState>::SingleStateCovariance::Identity();

    MTK::subblock (P0, &WVectorState::statek, &WVectorState::statek_l) = 0.0011 * localization::Usckf<WVectorState, WSingleState>::SingleStateCovariance::Identity();
    MTK::subblock (P0, &WVectorState::statek, &WVectorState::statek_i) = 0.0012 * localization::Usckf<WVectorState, WSingleState>::SingleStateCovariance::Identity();
    MTK::subblock (P0, &WVectorState::statek_l, &WVectorState::statek) = MTK::subblock (P0, &WVectorState::statek, &WVectorState::statek_l);
    MTK::subblock (P0, &WVectorState::statek_i, &WVectorState::statek) = MTK::subblock (P0, &WVectorState::statek, &WVectorState::statek_i);

    MTK::subblock (P0, &WVectorState::statek_l, &WVectorState::statek_i) = 0.0021 * localization::Usckf<WVectorState, WSingleState>::SingleStateCovariance::Identity();
    MTK::subblock (P0, &WVectorState::statek_i, &WVectorState::statek_l) = MTK::subblock (P0, &WVectorState::statek_l, &WVectorState::statek_i);

    //std::cout<<"[INITIALIZATION] P0:\n"<<P0<<std::endl;

    localization::Usckf<WVectorState, WSingleState> filter(static_cast<const WVectorState> (vstate),
                                                        static_cast<const WVectorState> (verror),
                                                        static_cast<const localization::Usckf<WVectorState,
                                                        WSingleState>::VectorStateCovariance> (P0));

    for (register int i=0; i<1; ++i)
    {
        Eigen::Vector3d acc, gyro;
        Eigen::Quaternion<double> orientq = filter.muState().statek_i.orient;
        acc << 100.00, 0.00, 0.00;
        gyro << (100.00*localization::D2R), 0.00, 0.00;
        std::cout<<"IN_LOOP ["<<i<<"]\n";
        localization::Usckf <WVectorState, WSingleState>::SingleStateCovariance myCov = processNoiseCov(dt);
        filter.predict(boost::bind(processModel, _1 , acc , gyro, orientq, dt), myCov);
        //std::cout<<"IN_LOOP Pk+i|k+l|k\n"<<filter.PkVectorState();
    }


    //std::cout<<"vstate.featuresk is "<<vstate.featuresk<<"\n";
    //std::cout<<"vstate.featuresk_l is "<<vstate.featuresk_l<<"\n";




}
