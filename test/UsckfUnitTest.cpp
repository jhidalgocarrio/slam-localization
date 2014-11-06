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
};

localization::MeasurementType measurementModelVO (const WAugmentedState &wastate)
{
    WSingleState delta_state, statek, statek_i; /** Propagated state */
    localization::MeasurementType z_hat;
    z_hat = wastate.featuresk;
    statek = wastate.statek;
    statek_i = wastate.statek_i;

    delta_state = statek - statek_i;
    Eigen::Affine3d delta_transform (delta_state.orient);
    delta_transform.translation() = delta_state.pos;

    for (register unsigned int i = 0; i < z_hat.size(); i+=3)
    {
        Eigen::Vector3d coord;
        coord<<wastate.featuresk[i], wastate.featuresk[i+1], wastate.featuresk[i+2];
        coord = delta_transform * coord;
        z_hat[i] = coord[0];
        z_hat[i+1] = coord[1];
        z_hat[i+2] = coord[2];
    }
//    std::cout<<"z_hat "<<z_hat<<"\n";

    return z_hat;
};

BOOST_AUTO_TEST_CASE( STATES )
{

    localization::MeasurementType featuresVO, featuresICP;
    featuresVO.resize(3,1);
    featuresVO<<3.34, 3.34, 3.34;
    featuresICP.resize(9,1);
    featuresICP<<1.34, 1.34, 1.34, 1.34,1.34, 1.34, 1.34, 1.34, 1.34;

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

    WAugmentedState vstatebis3;
    vstatebis3 = vstatebis;
    std::cout<<"[STATE] vstatebis3: "<<vstatebis3<<"\n";
}

BOOST_AUTO_TEST_CASE( OPERATIONS )
{
    localization::MeasurementType featuresVO, featuresICP;
    featuresVO.resize(3,1);
    featuresVO<<3.34, 3.34, 3.34;
    featuresICP.resize(9,1);
    featuresICP<<1.34, 1.34, 1.34, 1.34,1.34, 1.34, 1.34, 1.34, 1.34;

    WAugmentedState vstate;
    vstate.featuresk = featuresVO;
    vstate.featuresk_l = featuresICP;

    std::cout<<"[OPERATIONS] vstate: "<<vstate<<"\n";

    WAugmentedState vstatebis;
    vstatebis = vstate;
    vstatebis.statek_i.pos<< 20, 3.0, -4.00;

    Eigen::Matrix <double,localization::NUMAXIS,1> euler; /** In euler angles **/
    euler[2] = 1.00 * localization::D2R;
    euler[1] = 1.00 * localization::D2R;
    euler[0] = 1.00 * localization::D2R;

    vstatebis.statek_i.orient.boxplus(euler);
    vstatebis.statek = vstatebis.statek_i;
    vstatebis.statek_l = vstatebis.statek_i;
    std::cout<<"[OPERATIONS] vstatebis.getDOF(): "<<vstatebis.getDOF()<<"\n";
    std::cout<<"[OPERATIONS] vstatebis: "<<vstatebis<<"\n";

    /** Operation with states **/
    WAugmentedState sumstate, resstate;
    resstate = vstate - vstatebis;
    sumstate = vstate + vstatebis;
    std::cout<<"[OPERATIONS] vstate - vstatebis\n"<< vstate - vstatebis <<"\n";
    std::cout<<"[OPERATIONS] vstate + vstatebis\n"<< vstate + vstatebis <<"\n";

    std::cout<<"[OPERATIONS] resstate\n"<< resstate <<"\n";
    std::cout<<"[OPERATIONS] sumstate\n"<< sumstate <<"\n";
    std::cout<<"[OPERATIONS] resstate+sumstate\n"<< resstate+sumstate <<"\n";
    std::cout<<"[OPERATIONS] sumstate-resstate\n"<< sumstate-resstate <<"\n";
}

BOOST_AUTO_TEST_CASE( USCKF )
{
    WAugmentedState vstate;
    WSingleState state_single;

    std::cout<<"[USCKF] state_single: "<<state_single<<"\n";

    const double dt = 0.01; /** 100Hz */

    /** Initial covariance matrix **/
    StateFilter::SingleStateCovariance P0_single;
    P0_single = 0.0025 * StateFilter::SingleStateCovariance::Identity();

    std::cout<<"[USCKF] P0_single is of size "<<P0_single.rows() <<" x "<<P0_single.cols()<<std::endl;
    std::cout<<"[USCKF] P0_single:\n"<<P0_single<<std::endl;

    StateFilter filter(static_cast<const WSingleState> (state_single), static_cast<const StateFilter::SingleStateCovariance> (P0_single));

    std::cout<<"[USCKF] filter state: "<<filter.muState()<<"\n";
    std::cout<<"[USCKF] P0 is of size "<<filter.PkAugmentedState().rows() <<" x "<<filter.PkAugmentedState().cols()<<std::endl;
    std::cout<<"[USCKF] P0:\n"<<filter.PkAugmentedState()<<std::endl;

    /** First Measurements **/
    localization::MeasurementType featuresVO, featuresICP;
    featuresVO.resize(3,1);
    featuresVO<<3.34, 3.34, 3.34;
    featuresICP.resize(9,1);
    featuresICP<<1.34, 1.34, 1.34, 1.34,1.34, 1.34, 1.34, 1.34, 1.34;

    Eigen::Matrix<StateFilter::ScalarType, Eigen::Dynamic, Eigen::Dynamic> featuresVOCov, featuresICPCov;
    featuresVOCov.resize(featuresVO.size(), featuresVO.size());
    featuresVOCov.setIdentity(); featuresVOCov = 0.008 * featuresVOCov;
    featuresICPCov.resize(featuresICP.size(), featuresICP.size());
    featuresICPCov.setIdentity(); featuresICPCov = 0.008 * featuresICPCov;

    filter.setMeasurement<localization::MeasurementType, Eigen::Matrix<StateFilter::ScalarType, Eigen::Dynamic, Eigen::Dynamic> >(localization::STATEK, featuresVO, featuresVOCov);

    std::cout<<"[USCKF] filter state: "<<filter.muState()<<"\n";
    std::cout<<"[USCKF] P0 is of size "<<filter.PkAugmentedState().rows() <<" x "<<filter.PkAugmentedState().cols()<<std::endl;
    std::cout<<"[USCKF] P0:\n"<<filter.PkAugmentedState()<<std::endl;

    filter.setMeasurement<localization::MeasurementType, Eigen::Matrix<StateFilter::ScalarType, Eigen::Dynamic, Eigen::Dynamic> >(localization::STATEK_L, featuresICP, featuresICPCov);

    std::cout<<"[USCKF] filter state: "<<filter.muState()<<"\n";
    std::cout<<"[USCKF] P0 is of size "<<filter.PkAugmentedState().rows() <<" x "<<filter.PkAugmentedState().cols()<<std::endl;
    std::cout<<"[USCKF] P0:\n"<<filter.PkAugmentedState()<<std::endl;

    featuresVO<<3.35, 3.35, 3.35;
    featuresVOCov.setIdentity(); featuresVOCov = 0.05 * featuresVOCov;

    filter.setMeasurement<localization::MeasurementType, Eigen::Matrix<StateFilter::ScalarType, Eigen::Dynamic, Eigen::Dynamic> >(localization::STATEK, featuresVO, featuresVOCov);

    std::cout<<"[USCKF] filter state: "<<filter.muState()<<"\n";
    std::cout<<"[USCKF] P0 is of size "<<filter.PkAugmentedState().rows() <<" x "<<filter.PkAugmentedState().cols()<<std::endl;
    std::cout<<"[USCKF] P0:\n"<<filter.PkAugmentedState()<<std::endl;


    /***************************/
    /** PROPAGATION / PREDICT **/
    /***************************/
    std::cout<<"***************\n";
    std::cout<<"*** PREDICT ***\n";
    std::cout<<"***************\n";

    for (register int i=0; i<2; ++i)
    {
        Eigen::Vector3d velo, angular_velo;
        velo << 100.00, 0.00, 0.00;
        angular_velo << (100.00*localization::D2R), (100.00*localization::D2R), (100.00*localization::D2R);
        std::cout<<"[USCKF] IN_LOOP ["<<i<<"]\n";
        StateFilter::SingleStateCovariance myCov = processNoiseCov(dt);
        filter.predict(boost::bind(processModel, _1 , velo , angular_velo, dt), myCov);
       // std::cout<<"[USCKF] IN_LOOP Pk+i|k+l|k\n"<<filter.PkAugmentedState()<<"\n";
    }

    std::cout<<"[USCKF] vstate after propagation: "<<filter.muState()<<"\n";

    /***************************/
    /** CORRECTION / UPDATE   **/
    /***************************/
    std::cout<<"**************\n";
    std::cout<<"*** UPDATE ***\n";
    std::cout<<"**************\n";

    /** Measurement model **/
    localization::MeasurementType featuresVO_hat;
    featuresVO_hat = measurementModelVO(filter.muState());
    std::cout<<"[USCKF] featuresVO_hat: "<<featuresVO_hat<<"\n";

    /** Measurement **/
    Eigen::Matrix<StateFilter::ScalarType, Eigen::Dynamic, 1> measurementVO;
    std::vector< Eigen::Matrix<double, 3, 1> > vector_featuresVO;
    vector_featuresVO.resize(1);
    vector_featuresVO[0]<<2.33, 3.35, 3.35;
    measurementVO.resize(vector_featuresVO[0].rows()*vector_featuresVO.size(), 1);
    measurementVO.setZero();
    std::cout<<"[USCKF] MeasurementVO size is "<<measurementVO.size()<<"\n";
    register int k = 0;
    for (register size_t i=0; i<vector_featuresVO.size(); ++i)
    {
        measurementVO.block(0+k, 0, vector_featuresVO[i].rows(), vector_featuresVO[i].cols()) = vector_featuresVO[i];
        k=k+vector_featuresVO[i].rows();
    }

    std::cout<<"[USCKF] Measurement:\n"<<measurementVO<<"\n";
    Eigen::Matrix<StateFilter::ScalarType, Eigen::Dynamic, Eigen::Dynamic> measurementNoiseVO;
    measurementNoiseVO.resize(measurementVO.size(), measurementVO.size());
    measurementNoiseVO.setIdentity();
    measurementNoiseVO = 0.01 * measurementNoiseVO;
    filter.update(measurementVO, boost::bind(measurementModelVO, _1), measurementNoiseVO);

    /***********/
    /** END   **/
    /***********/

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
