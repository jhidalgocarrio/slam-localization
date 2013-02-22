/**\file Datatypes.hpp
 * Header function file and defines
 */

#ifndef _DATATYPES_HPP_
#define _DATATYPES_HPP_

#include <base/time.h> /** For the timestamp **/
#include <base/eigen.h> /** For the matrices for orogen **/


namespace localization	
{
    struct FilterInfo
    {
	
	base::Time time;
	base::VectorXd xki_k;
	base::MatrixXd Pki_k;
	base::MatrixXd K;
	base::MatrixXd Qk;
	base::MatrixXd Rk;
	base::MatrixXd Hk;
	base::VectorXd zki;
	base::VectorXd innovation;
	base::Vector3d eposition;
	base::Vector3d evelocity;
	base::Matrix3d evelocity_cov;
	base::Vector3d bghat;
	base::Vector3d bahat;
	base::Matrix3d Hellinger;
	base::Vector3d eacceleration;
	double Mahalanobis;
    };
    
    struct MeasurementGenerationInfo
    {
	base::Time time;
	double wlsNavigation;
	double wlsSlip;
    };
    
    struct SlipInfo
    {
	base::Time time;
	base::VectorXd slip_vector;
	base::MatrixXd Cov;
    };

}
#endif