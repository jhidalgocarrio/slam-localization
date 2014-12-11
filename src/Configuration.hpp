/**\file configuration.hpp
 * Header function file and defines
 */

#ifndef _CONFIGURATION_HPP_
#define _CONFIGURATION_HPP_

namespace localization	
{

    /** General constants values **/
    static const int QUATERSIZE = 4; /** Number of parameters of a quaternion **/
    static const int NUMAXIS = 3; /** Number of axis sensed by the IMU **/

    /** WGS-84 ellipsoid constants (Nominal Gravity Model and Earth angular velocity) **/
    static const int Re = 6378137; /** Equatorial radius in meters **/
    static const int Rp = 6378137; /** Polar radius in meters **/
    static const double ECC = 0.0818191908426; /** First eccentricity **/
    static const double GRAVITY = 9.79766542; /** Mean value of gravity value in m/s^2 **/
    static const double GWGS0 = 9.7803267714; /** Gravity value at the equator in m/s^2 **/
    static const double GWGS1 = 0.00193185138639; /** Gravity formula constant **/
    static const double EARTHW = 7.292115e-05; /** Earth angular velocity in rad/s **/

    /** Magnetic declination **/
    enum DECLINATION_CONSTS {
         EAST = 1, /** EAST is 1 and means positive magnetic declination **/
         WEST = 2 /** WEST is 2 and means negative magnetic declination **/
    };

    /** Variables for the attitude estimation inside the algorithm **/
    static const int M1 = 1; /** Parameter for adaptive algorithm (to estimate Uk with is not directly observale) */
    static const int M2 = 5; /** Parameter for adaptive algorithm (to prevent falsering entering in no-external acc mode) */
    static const double GAMMA = 0.0005; /** Parameter for adaptive algorithm (only entering when Qstart is greater than RHR'+Ra) */
    static const int R2COUNT = 100; /** Parameter for adaptive algorithm */

    static const double D2R = M_PI/180.00; /** Convert degree to radian **/
    static const double R2D = 180.00/M_PI; /** Convert radian to degree **/

    static const double ZERO_UNCERTAINTY = 1.0e-10; /** Set as default zero uncertainty **/


    /** Integration of the delayed windows **/
    static const int INTEGRATION_XAXIS_WINDOW_SIZE = 1;//100; /** Windows size of the delay integration **/

    static const int INTEGRATION_YAXIS_WINDOW_SIZE = 1;//100; /** Windows size of the delay integration **/

    static const int INTEGRATION_ZAXIS_WINDOW_SIZE = 1;//10; /** Windows size of the delay integration **/

    static const int ANGVELO_WINDOW_SIZE = INTEGRATION_XAXIS_WINDOW_SIZE; /** Windows size of the delay integration **/
	
    static const unsigned int NUMBER_OF_WHEELS = 4; /** Rover number of wheels **/

    static const unsigned int NUMBER_OF_PASSIVE_JOINTS = 1; /** Rover chassis number of passive joints **/

    static const unsigned int NUMBER_OF_ACTIVE_JOINTS = 0; /** Rover chassis number of actuated joints **/

    static const unsigned int ENCODERS_VECTOR_SIZE = NUMBER_OF_ACTIVE_JOINTS+NUMBER_OF_PASSIVE_JOINTS+NUMBER_OF_WHEELS; /** Vector of rover sensed joints with encoders **/

    static const unsigned int SLIP_VECTOR_SIZE = NUMAXIS * NUMBER_OF_WHEELS; /** Size of slip vector of the mobil robot **/

    static const unsigned int NORDER_BESSEL_FILTER = 8; /** Order of the IIR Bessel filter **/

}

#endif //
