#define BOOST_TEST_MODULE template_for_test_test
#include <boost/test/included/unit_test.hpp>

/** Library **/
#include <localization/filters/Usckf.hpp> /** USCKF class with Manifolds */
#include <localization/filters/MtkWrap.hpp> /** USCKF wrapper for the state vector */
#include <localization/filters/State.hpp> /** Filters State */
#include <localization/Configuration.hpp> /** Constant values of the library */

/** Eigen **/
#include <Eigen/Core> /** Core */
#include <Eigen/StdVector> /** For STL container with Eigen types **/

/** Standard libs **/
#include <iostream>
#include <vector>

std::vector<int> hola(100, 2);
typedef std::vector< MTK::vect<3, double> > MTKintFeature;

MTKintFeature adios(100, Eigen::Vector3d::Zero() );

/** Wrap the AugmentedState **/
typedef localization::MtkWrap<localization::State> WSingleState;
typedef localization::MtkDynamicWrap<localization::AugmentedState> WAugmentedState;
typedef localization::Usckf<WAugmentedState, WSingleState > StateFilter;



WSingleState processModel (const WSingleState &state,  const Eigen::Vector3d &velocity, const Eigen::Vector3d &angular_velocity, double dt)
{
    WSingleState s2; /** Propagated state */

    /** Apply Rotation **/
    Eigen::Vector3d scaled_axis = angular_velocity * dt;
    localization::SO3 rot = localization::SO3::exp (scaled_axis);
    s2.orient = state.orient * rot ;
    s2.angvelo = angular_velocity;

    /** Apply Translation **/
    s2.velo = velocity;
    s2.pos = state.pos + state.velo * dt;

    return s2;
};

StateFilter::SingleStateCovariance processNoiseCov (double dt)
{
    StateFilter::SingleStateCovariance cov = StateFilter::SingleStateCovariance::Zero();
    MTK::setDiagonal (cov, &WSingleState::pos, 0.1 * dt);
    MTK::setDiagonal (cov, &WSingleState::orient, 0.1 * dt);
    MTK::setDiagonal (cov, &WSingleState::velo,  0.1 * dt);
    MTK::setDiagonal (cov, &WSingleState::angvelo,  0.1 * dt);

    return cov ;
}

BOOST_AUTO_TEST_CASE( STATES )
{

    localization::MeasurementType featuresVO, featuresICP;
    featuresVO.resize(4,1);
    featuresVO<<3.34, 3.34, 3.34, 3.34;
    featuresICP.resize(10,1);
    featuresICP<<1.34, 1.34, 1.34, 1.34,1.34, 1.34, 1.34, 1.34, 1.34, 1.34;

    WAugmentedState vstate;

    std::cout<<"[STATE] vstate::DOF is "<<vstate.DOF<<"\n";
    std::cout<<"[STATE] vstate.statek is "<<vstate.statek<<"\n";
    std::cout<<"[STATE] vstate.statek_l is "<<vstate.statek_l<<"\n";
    std::cout<<"[STATE] vstate.statek_i is "<<vstate.statek_i<<"\n";
    std::cout<<"[STATE] vstate: "<<vstate<<"\n";
    std::cout<<"[STATE] vectorized state:\n"<<vstate.getVectorizedState()<<"\n";
    std::cout<<"[STATE] featuresICP.size():\n"<<featuresICP.size() <<"\n";
    std::cout<<"[STATE] featuresICP:\n"<<featuresICP.matrix()<<"\n";

    //std::vector< Eigen::Matrix <double, 3, 1> , Eigen::aligned_allocator < Eigen::Matrix <double, 3, 1> > > featuresVO;

    std::cout<<"[STATE] size of featuresk: "<<vstate.featuresk.size()<<"\n";
    vstate.featuresk = featuresVO;
    std::cout<<"[STATE] size of featuresk: "<<vstate.featuresk.size()<<"\n";
    vstate.featuresk_l = featuresICP;
    std::cout<<"[STATE] size of featuresk_l: "<<vstate.featuresk_l.size()<<"\n";
    std::cout<<"[STATE] vstate: "<<vstate<<"\n";
    std::cout<<"[STATE] vectorized state:\n"<<vstate.getVectorizedState()<<"\n";

    localization::AugmentedState vstatebis(vstate.statek, vstate.statek_l, vstate.statek_i, featuresVO, featuresICP);
    std::cout<<"[STATE] vstatebis::DOF is "<<vstatebis.getDOF()<<"\n";
    std::cout<<"[STATE] vstatebis: "<<vstatebis<<"\n";

    localization::AugmentedState vstatebis2;
    std::cout<<"[STATE] vstatebis2::DOF is "<<vstatebis2.getDOF()<<"\n";
    vstatebis2 = vstatebis;
    std::cout<<"[STATE] vstatebis2::DOF is "<<vstatebis2.getDOF()<<"\n";
    std::cout<<"[STATE] vstatebis2: "<<vstatebis2<<"\n";

}

BOOST_AUTO_TEST_CASE( USCKF )
{
    localization::MeasurementType featuresVO, featuresICP;
    featuresVO.resize(4,1);
    featuresVO<<3.34, 3.34, 3.34, 3.34;
    featuresICP.resize(10,1);
    featuresICP<<1.34, 1.34, 1.34, 1.34,1.34, 1.34, 1.34, 1.34, 1.34, 1.34;

    WAugmentedState vstate;
    vstate.featuresk = featuresVO;
    vstate.featuresk_l = featuresICP;

    std::cout<<"[USCKF] vstate: "<<vstate<<"\n";

    const double dt = 0.01; /** 100Hz */

    /** Initial covariance matrix **/
    StateFilter::AugmentedStateCovariance P0;
    StateFilter::MultiStateCovariance P0_states;

    MTK::subblock (P0_states, &WAugmentedState::statek, &WAugmentedState::statek) = 0.0025 * StateFilter::SingleStateCovariance::Identity();
    MTK::subblock (P0_states, &WAugmentedState::statek_l, &WAugmentedState::statek_l) = 0.0035 * StateFilter::SingleStateCovariance::Identity();
    MTK::subblock (P0_states, &WAugmentedState::statek_i, &WAugmentedState::statek_i) = 0.0045 * StateFilter::SingleStateCovariance::Identity();

    MTK::subblock (P0_states, &WAugmentedState::statek, &WAugmentedState::statek_l) = 0.0011 * StateFilter::SingleStateCovariance::Identity();
    MTK::subblock (P0_states, &WAugmentedState::statek, &WAugmentedState::statek_i) = 0.0012 * StateFilter::SingleStateCovariance::Identity();
    MTK::subblock (P0_states, &WAugmentedState::statek_l, &WAugmentedState::statek) = MTK::subblock (P0_states, &WAugmentedState::statek, &WAugmentedState::statek_l);
    MTK::subblock (P0_states, &WAugmentedState::statek_i, &WAugmentedState::statek) = MTK::subblock (P0_states, &WAugmentedState::statek, &WAugmentedState::statek_i);

    MTK::subblock (P0_states, &WAugmentedState::statek_l, &WAugmentedState::statek_i) = 0.0021 * StateFilter::SingleStateCovariance::Identity();
    MTK::subblock (P0_states, &WAugmentedState::statek_i, &WAugmentedState::statek_l) = MTK::subblock (P0_states, &WAugmentedState::statek_l, &WAugmentedState::statek_i);

    P0.resize(vstate.getDOF(), vstate.getDOF());
    P0.block(0, 0, P0_states.rows(), P0_states.cols()) = P0_states;

    //std::cout<<"[USCKF] P0_states:\n"<<P0_states<<std::endl;
    std::cout<<"[USCKF] P0_states is of size "<<P0.rows() <<" x "<<P0.cols()<<std::endl;
    std::cout<<"[USCKF] P0_states is of size "<<P0_states.rows() <<" x "<<P0_states.cols()<<std::endl;

    StateFilter filter(static_cast<const WAugmentedState> (vstate), static_cast<const StateFilter::AugmentedStateCovariance> (P0));

    /***************************/
    /** PROPAGATION / PREDICT **/
    /***************************/
    for (register int i=0; i<2; ++i)
    {
        Eigen::Vector3d velo, angular_velo;
        velo << 100.00, 0.00, 0.00;
        angular_velo << (100.00*localization::D2R), (100.00*localization::D2R), (100.00*localization::D2R);
        std::cout<<"[USCKF] IN_LOOP ["<<i<<"]\n";
        StateFilter::SingleStateCovariance myCov = processNoiseCov(dt);
        filter.predict(boost::bind(processModel, _1 , velo , angular_velo, dt), myCov);
        std::cout<<"[USCKF] IN_LOOP Pk+i|k+l|k\n"<<filter.PkAugmentedState()<<"\n";
    }


    /***************************/
    /** CORRECTION / UPDATE   **/
    /***************************/

    vstate = filter.muState();
    Eigen::Matrix <double,localization::NUMAXIS,1> euler; /** In euler angles **/
    euler[2] = vstate.statek_i.orient.toRotationMatrix().eulerAngles(2,1,0)[0];//Yaw
    euler[1] = vstate.statek_i.orient.toRotationMatrix().eulerAngles(2,1,0)[1];//Pitch
    euler[0] = vstate.statek_i.orient.toRotationMatrix().eulerAngles(2,1,0)[2];//Roll
    std::cout<<"[USCKF] Result Roll: "<<euler[0]*localization::R2D<<" Pitch: "<<euler[1]*localization::R2D<<" Yaw: "<<euler[2]*localization::R2D<<"\n";
    std::cout<<"[USCKF] vectorized state: "<<vstate.getVectorizedState()<<"\n";

    Eigen::AngleAxisd angleAxis (vstate.statek_i.orient);
    std::cout<<"[USCKF] The angle of rotation is: "<<angleAxis.angle()<<"\n";
    std::cout<<"[USCKF] The angle of rotation is: "<<(angleAxis.axis() * angleAxis.angle()).norm()<<"\n";
    std::cout<<"[USCKF] The angle of rotation is (degrees): "<<angleAxis.angle()*localization::R2D<<"\n";
    std::cout<<"[USCKF] The axis of rotation is:\n"<<angleAxis.axis()<<"\n";

    euler[2] = 0.00 * localization::D2R;
    euler[1] = 0.00 * localization::D2R;
    euler[0] = 0.30 * localization::D2R;

    vstate.statek_i.orient.boxplus(euler);
    euler[2] = vstate.statek_i.orient.toRotationMatrix().eulerAngles(2,1,0)[0];//Yaw
    euler[1] = vstate.statek_i.orient.toRotationMatrix().eulerAngles(2,1,0)[1];//Pitch
    euler[0] = vstate.statek_i.orient.toRotationMatrix().eulerAngles(2,1,0)[2];//Roll
    std::cout<<"[USCKF] Result(Euler) Roll: "<<euler[0]*localization::R2D<<" Pitch: "<<euler[1]*localization::R2D<<" Yaw: "<<euler[2]*localization::R2D<<"\n";

    euler[2] = 0.05 * localization::D2R;
    euler[1] = 0.05 * localization::D2R;
    euler[0] = 0.05 * localization::D2R;

    vstate.statek_i.orient.boxplus(euler);
    euler[2] = vstate.statek_i.orient.toRotationMatrix().eulerAngles(2,1,0)[0];//Yaw
    euler[1] = vstate.statek_i.orient.toRotationMatrix().eulerAngles(2,1,0)[1];//Pitch
    euler[0] = vstate.statek_i.orient.toRotationMatrix().eulerAngles(2,1,0)[2];//Roll
    std::cout<<"[USCKF] Result(Euler2) Roll: "<<euler[0]*localization::R2D<<" Pitch: "<<euler[1]*localization::R2D<<" Yaw: "<<euler[2]*localization::R2D<<"\n";

    std::cout<<"[USCKF] vectorized state: "<<vstate.getVectorizedState()<<"\n";


    //std::cout<<"vstate.featuresk is "<<vstate.featuresk<<"\n";
    //std::cout<<"vstate.featuresk_l is "<<vstate.featuresk_l<<"\n";

}
