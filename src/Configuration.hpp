/**\file configuration.hpp
 * Header function file and defines
 */

#ifndef _CONFIGURATION_HPP_
#define _CONFIGURATION_HPP_

namespace localization	
{

    /** General defines **/
    #ifndef OK
    #define OK	0  /** Integer value in order to return when everything is all right. */
    #endif
    #ifndef ERROR_OUT
    #define ERROR_OUT -1  /** Integer value in order to return when an error occured. */
    #endif
    
    #ifndef QUATERSIZE
    #define QUATERSIZE 4 /** Number of parameters of a quaternion **/
    #endif
    
    
    /** WGS-84 ellipsoid constants (Nominal Gravity Model and Earth angular velocity) **/
    #ifndef Re
    #define Re	6378137 /** Equatorial radius in meters **/
    #endif
    #ifndef Rp
    #define Rp	6378137 /** Polar radius in meters **/
    #endif
    #ifndef ECC
    #define ECC  0.0818191908426 /** First eccentricity **/
    #endif
    #ifndef GRAVITY
    #define GRAVITY 9.79766542 /** Mean value of gravity value in m/s^2 **/
    #endif
    #ifndef GWGS0
    #define GWGS0 9.7803267714 /** Gravity value at the equator in m/s^2 **/
    #endif
    #ifndef GWGS1
    #define GWGS1 0.00193185138639 /** Gravity formula constant **/
    #endif
    #ifndef EARTHW
    #define EARTHW  7.292115e-05 /** Earth angular velocity in rad/s **/
    #endif
    
    #ifndef EAST
    #define EAST 1 /** EAST is 1 and means positive magnetic declination **/
    #endif
    
    #ifndef WEST
    #define WEST 2 /** WEST is 2 and means negative magnetic declination **/
    #endif
    
    /** Inertial Sensors constant parameters **/
    #ifndef NUMAXIS
    #define NUMAXIS 3 /** Number of axis sensed by the IMU **/
    #endif
    
    /** Variables for the attitude estimation inside the algorithm **/
    #define M1 3 /** Parameter for adaptive algorithm */
    #define M2 3 /** Parameter for adaptive algorithm (to prevent falsering entering in no-external acc mode) */
    #define GAMMA 0.01 /** Parameter for adaptive algorithm (only entering when Qstart is greater than RHR'+Ra)*/
    #define R2COUNT 100 /** Parameter for adaptive algorithm */

    #ifndef D2R
    #define D2R M_PI/180.00 /** Convert degree to radian **/
    #endif
    
    #ifndef R2D
    #define R2D 180.00/M_PI /** Convert radian to degree **/
    #endif
    
    #ifndef ZERO_UNCERTAINTY
    #define ZERO_UNCERTAINTY 1.0e-10 /** Set as default zero uncertainty **/
    #endif
    
    static const unsigned int NUMBER_OF_WHEELS = 4; /** Rover number of wheels **/
    
    static const unsigned int NUMBER_OF_PASSIVE_JOINTS = 1; /** Rover chassis number of passive joints **/
    
    static const unsigned int NUMBER_OF_ACTIVE_JOINTS = 0; /** Rover chassis number of actuated joints **/
    
    static const unsigned int ENCODERS_VECTOR_SIZE = NUMBER_OF_ACTIVE_JOINTS+NUMBER_OF_PASSIVE_JOINTS+NUMBER_OF_WHEELS; /** Vector of rover sensed joints with encoders **/
    
    static const unsigned int SLIP_VECTOR_SIZE = NUMAXIS * NUMBER_OF_WHEELS; /** Size of slip vector of the mobil robot **/
    
}

#endif // 