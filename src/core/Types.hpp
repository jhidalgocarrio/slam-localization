/**\file Datatypes.hpp
 * Header function file and defines
 */

#ifndef _DATATYPES_HPP_
#define _DATATYPES_HPP_

#include <base/time.h> /** For the timestamp **/
#include <base/eigen.h> /** For the matrices for orogen **/
#include <base/samples/rigid_body_state.h> /** Body State and Pose **/


namespace localization	
{
    /**************************************/
    /** Data struct for (internal) ports **/
    /**************************************/

    //Inertial sensor and extra information
    struct InertialState
    {
        /** Time stamp */
        base::Time time;

        /** Theoretical gravity value Usually computed from model */
        double theoretical_g;

        /** Experimental gravity value computed at initialization time (no-moving robot) */
        double estimated_g;

        /** On/Off initial accelerometers bias */
        base::Vector3d abias_onoff;

        /** On/Off initial gyro bias */
        base::Vector3d gbias_onoff;

    };

    /**************************************/
    /** Data struct for (internal) ports **/
    /**************************************/

    //Optimization problem information
    struct StateEstimation
    {
	base::Time time;
	base::VectorXd statek_i;
	base::VectorXd errork_i;
        base::Orientation orientation;
	base::MatrixXd Pki;
	base::MatrixXd K;
	base::MatrixXd Qk;
	base::MatrixXd Rk;
	base::VectorXd innovation;
	base::Matrix3d Hellinger;
	base::Matrix3d Threshold;
        double mahalanobis;
        base::Vector3d abias;
        base::Vector3d gbias;
        base::Vector3d accModel;
        base::Matrix3d accModelCov;
        base::Vector3d accInertial;
        base::Matrix3d accInertialCov;
        base::Vector3d accError;
        base::Matrix3d accErrorCov;
        base::Vector3d deltaVeloCommon;
        base::Matrix3d deltaVeloCommonCov;
    };

}
#endif
