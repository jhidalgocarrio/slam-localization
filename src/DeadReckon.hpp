/**\file DeadReckon.hpp
 * Header function file and defines
 */

#ifndef _DEADRECKON_HPP_
#define _DEADRECKON_HPP_

#include <Eigen/Core> /** Core methods of Eigen implementation **/
#include <Eigen/Dense> /** for the algebra and transformation matrices **/
#include <Eigen/Geometry> /** Eigen data type for Matrix, Quaternion, etc... */
#include <envire/core/Transform.hpp> /** Envire module which has transformation with uncertainty **/
#include "Configuration.hpp" /** For the localization framework constant and configuration values **/

namespace localization	
{
    
    /** Class to perform Dead Reckoning assuming constant acceleration **/
    class DeadReckon
    {
    
    public:
	
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW
	/** Dead reckoning variables **/
	Eigen::Matrix <double, NUMAXIS, 2> vState; /** robot state linear velocities time nT col(0) and (n-1)T col(1) (order is x,y,z) **/
	envire::TransformWithUncertainty actualPose; /** robot pose (rbsBC) in transformation with uncertainty form **/
	

    public:
	/** \Brief Default constructor
	*/
	DeadReckon();
	
	/** \Brief Constructor
	*/
	DeadReckon(base::samples::RigidBodyState initPose);
	
	/** \Brief Initialization method
	*/
	void setInitPose(base::samples::RigidBodyState initPose);
	
	/** \Brief Performs the time integration of delta pose updates
	 * 
	 * @return void
	 */
	base::samples::RigidBodyState updatePose(Eigen::Matrix< double, NUMAXIS , 1  > linvelo, Eigen::Quaternion <double> delta_q,
						     Eigen::Matrix< double, NUMAXIS , NUMAXIS  > covlinvelo, Eigen::Matrix< double, NUMAXIS , NUMAXIS  > covdelta_q,
						     Eigen::Matrix< double, NUMAXIS , 1  > linvelo_error, Eigen::Quaternion <double> delta_qe,
						     Eigen::Matrix< double, NUMAXIS , NUMAXIS  > covlinvelo_error, Eigen::Matrix< double, NUMAXIS , NUMAXIS  > covdelta_qe,
						     base::Time timeStamp, double delta_t);

    };

}

#endif