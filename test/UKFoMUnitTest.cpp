#define BOOST_TEST_MODULE template_for_test_test
#include <boost/test/included/unit_test.hpp>

// MTK's pose and orientation definition:
#include <mtk/types/pose.hpp>
#include <mtk/types/SOn.hpp>
#include <mtk/build_manifold.hpp>

//UKF class with Manifolds
#include <ukfom/ukf.hpp>

//UKF wrapper for the state vector
#include <ukfom/mtkwrap.hpp>

/** Standard libs **/
#include <iostream>
#include <vector>

#ifndef D2R
#define D2R M_PI/180.00 /** Convert degree to radian **/
#endif
#ifndef R2D
#define R2D 180.00/M_PI /** Convert radian to degree **/
#endif

// We can't use types having a comma inside AutoConstruct macros :(
typedef MTK::vect<3, double> vec3;
typedef MTK::SO3<double> SO3;


MTK_BUILD_MANIFOLD ( mtk_state ,
(( vec3, pos ))
(( SO3, orient ))
(( vec3, vel ))
);


// Alternatively use predefined trafo-Type (already implements transformations):
typedef MTK::trafo<MTK::SO2<double> > Pose;

typedef ukfom::mtkwrap<mtk_state> state;

//a is acceleration
//w is angular velocity (gyros)
state process_model ( const state &s, const Eigen::Matrix<double, 3, 1> &a , const MTK::vect <3,double > &w, double dt)
{
    state s2 ;

    // apply rotation
    std::cout<<"gyros are: "<<w<<"\n";
    //MTK::vect <3, double> scaled_axis = w * dt ;
    Eigen::Vector3d scaled_axis = w * dt ;
    s2.orient.boxplus(scaled_axis); // this is the same than the two lines bellow
    //MTK::SO3 <double> rot = SO3::exp (scaled_axis);
    //s2.orient = s.orient * rot;


    // accelerate with gravity
    //Eigen::Matrix<double, 3, 1> gravity;
    MTK::vect<3, double> gravity;
    gravity <<0,0,9.81;
    s2.vel = s.vel + (s.orient * a + gravity) * dt;

    // translate
    s2.pos = s.pos + s.vel * dt;

    std::cout<<"[PROCESS_MODEL] pos:"<<s2.pos<<" orient:"<<s2.orient<<" vel:"<<s2.vel<<"\n";

    return s2;
}


ukfom::ukf <state>::cov process_noise_cov (double dt)
{
    ukfom::ukf<state>::cov cov = ukfom::ukf<state>::cov::Zero();
    MTK::setDiagonal (cov, &state::pos, 0);
    MTK::setDiagonal (cov, &state::orient,  0.0001 * dt);
    MTK::setDiagonal (cov, &state::vel, 0.0002 * dt);
    return cov ;
}

MTK::vect <3, double> gps_measurement_model ( const state &s )
{
    return s.pos ;
}

Eigen::Matrix<double, 3, 3> gps_measurement_noise_cov ()
{
    return 0.00000001 * Eigen::Matrix<double, 3, 3>::Identity();
}


BOOST_AUTO_TEST_CASE( UKFOM )
{

    const state ukf_state;

    std::cout<<"ukf_state::DOF is "<<ukf_state.DOF<<"\n";
    std::cout<<"ukf_state.pos is "<<ukf_state.pos<<"\n";
    std::cout<<"ukf_state.orient is "<<ukf_state.orient<<"\n";
    std::cout<<"ukf_state.vel is "<<ukf_state.vel<<"\n";

    const ukfom::ukf<state>::cov init_cov = 0.001 * ukfom::ukf<state>::cov::Identity();
    ukfom::ukf<state> filter(ukf_state, init_cov);
    Eigen::Vector3d acc, gyro, gps;
    gyro<<10.0*D2R, 0.0*D2R, 0.0*D2R;
    acc<<0.0, 0.0, 0.0;
    gps<<1.0, 0.0, 0.0;
    for (register int i=0; i<1; ++i)
    {
        std::cout<<"IN_LOOP["<<i<<"]\n";
        //boost::bind(process_model, _1 , acc, gyro, 0.01);
        filter.predict(boost::bind(process_model, _1 , acc , gyro, 0.01), process_noise_cov(0.01));
        filter.update (gps, boost::bind(gps_measurement_model, _1), gps_measurement_noise_cov);
    }

    state vstate = filter.mu();
    Eigen::Matrix <double, 3, 1> euler; /** In Euler angles **/
    euler[2] = vstate.orient.toRotationMatrix().eulerAngles(2,1,0)[0];//Yaw
    euler[1] = vstate.orient.toRotationMatrix().eulerAngles(2,1,0)[1];//Pitch
    euler[0] = vstate.orient.toRotationMatrix().eulerAngles(2,1,0)[2];//Roll
    std::cout<<"Result Roll: "<<euler[0]*R2D<<" Pitch: "<<euler[1]*R2D<<" Yaw: "<<euler[2]*R2D<<"\n";
    ::MTK::SO3<double> orientation = vstate.orient;
    std::cout<<"Orientation is: "<<vstate.orient<<"\n";
    std::cout<<"Orientation is: "<<orientation<<"\n";

}
