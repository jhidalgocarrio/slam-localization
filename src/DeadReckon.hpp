/**\file DeadReckon.hpp
 * Header function file and defines
 */

#ifndef _DEAD_RECKON_HPP_
#define _DEAD_RECKON_HPP_

#include <Eigen/Core> /** Core methods of Eigen implementation **/
#include <Eigen/Dense> /** for the algebra and transformation matrices **/
#include <Eigen/Geometry> /** Eigen data type for Matrix, Quaternion, etc... */
#include <Eigen/StdVector> /** For STL container with Eigen types **/
#include <envire/core/Transform.hpp> /** Envire module which has transformation with uncertainty **/
#include "Configuration.hpp" /** For the localization framework constant and configuration values **/

#include <rover_localization/Util.hpp> /**Helper class of the framework **/

//#define DEAD_RECKON_DEBUG_PRINTS 1

namespace localization	
{

    /** Class to perform Dead Reckoning assuming constant acceleration **/
    class DeadReckon
    {

    public:
	
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    public:
	
	/** \Brief Performs the time integration of delta pose updates
	 *
	 * @return delta pose of the change in pose
         *
	 */
	static base::samples::RigidBodyState updatePose(const double delta_t,
                            const std::vector< Eigen::Matrix <double, 2*NUMAXIS, 1> , Eigen::aligned_allocator < Eigen::Matrix <double, 2*NUMAXIS, 1> > > &cartesianVelocities,
                            const Eigen::Matrix <double, 2*NUMAXIS, 2*NUMAXIS> &cartesianVelCov,
                            const base::samples::RigidBodyState &prevPose,
                            base::samples::RigidBodyState &postPose,
                            bool useTranforWithUncertainty = false)
        {
            base::samples::RigidBodyState deltaPose; /** rbs form of the computation **/
            envire::TransformWithUncertainty tfDeltaPose; /** delta transformation between current pose and next pose **/
            envire::TransformWithUncertainty tfPrevPose; /** Previous pose**/
            envire::TransformWithUncertainty tfPostPose; /** Posterior pose **/

            deltaPose.invalidate();

            /** Calculate the delta position from velocity (dead reckoning) assuming constant acceleration **/
            deltaPose.position = (delta_t/2.0) * (cartesianVelocities[0].block<NUMAXIS,1> (0,0) + cartesianVelocities[1].block<NUMAXIS,1>(0,0));

            /** Calculate the delta quaternion **/
            std::vector< Eigen::Matrix < double, NUMAXIS, 1 > , Eigen::aligned_allocator
                < Eigen::Matrix <double, NUMAXIS, 1> > > angularVelocities (cartesianVelocities.size(), Eigen::Matrix < double, NUMAXIS, 1 >::Zero());
            Eigen::Matrix < double, NUMAXIS, NUMAXIS > angularVelCov;
            angularVelocities[0] = cartesianVelocities[0].block<NUMAXIS, 1> (NUMAXIS, 0);
            angularVelocities[1] = cartesianVelocities[1].block<NUMAXIS, 1> (NUMAXIS, 0);
            angularVelCov = cartesianVelCov.block<NUMAXIS, NUMAXIS> (NUMAXIS, NUMAXIS);

            deltaPose.orientation = updateAttitude(delta_t, angularVelocities, angularVelCov);

            /** Set the uncertainty matrices **/
            if (Util::isnotnan(cartesianVelCov))
            {
                //deltaPose.cov_position <<  1.38911051e-04,   5.85928263e-08,  3.73528732e-05,
                //                           5.85928263e-08,   1.39249391e-04,   2.03173653e-08,
                //                           3.73528732e-05,   2.03173653e-08,   1.39795352e-04;
                deltaPose.cov_position = cartesianVelCov.block<NUMAXIS, NUMAXIS>(0, 0) * delta_t * delta_t;
                deltaPose.cov_orientation = cartesianVelCov.block<NUMAXIS, NUMAXIS>(NUMAXIS, NUMAXIS) * delta_t * delta_t;
            }
            else
            {
                deltaPose.cov_position.setZero();
                deltaPose.cov_orientation.setZero();
            }

            #ifdef DEAD_RECKON_DEBUG_PRINTS
            std::cout<<"[DR] delta_t: "<<delta_t<<"\n";
            std::cout<<"[DR] cartesianVelocities[0]:\n"<<cartesianVelocities[0]<<"\n";
            std::cout<<"[DR] cartesianVelocities[1]:\n"<<cartesianVelocities[1]<<"\n";
            std::cout<<"[DR] deltaPose.position\n" <<deltaPose.position<<"\n";
            std::cout<<"[DR] deltaPose.cov_position\n" <<deltaPose.cov_position<<"\n";
            Eigen::Matrix <double,localization::NUMAXIS,1> euler; /** In euler angles **/
            euler[2] = deltaPose.orientation.toRotationMatrix().eulerAngles(2,1,0)[0];//Yaw
            euler[1] = deltaPose.orientation.toRotationMatrix().eulerAngles(2,1,0)[1];//Pitch
            euler[0] = deltaPose.orientation.toRotationMatrix().eulerAngles(2,1,0)[2];//Roll
            std::cout<<"[DR] delta Roll: "<<euler[0]*localization::R2D<<" Pitch: "<<euler[1]*localization::R2D<<" Yaw: "<<euler[2]*localization::R2D<<"\n";
            std::cout<<"[DR] deltaPose.cov_orientation\n" <<deltaPose.cov_orientation<<"\n";
            std::cout<<"[DR] prevPose.position\n" <<prevPose.position<<"\n";
            std::cout<<"[DR] prevPose.cov_position\n" <<prevPose.cov_position<<"\n";
            euler[2] = prevPose.orientation.toRotationMatrix().eulerAngles(2,1,0)[0];//Yaw
            euler[1] = prevPose.orientation.toRotationMatrix().eulerAngles(2,1,0)[1];//Pitch
            euler[0] = prevPose.orientation.toRotationMatrix().eulerAngles(2,1,0)[2];//Roll
            std::cout<<"[DR] prev Roll: "<<euler[0]*localization::R2D<<" Pitch: "<<euler[1]*localization::R2D<<" Yaw: "<<euler[2]*localization::R2D<<"\n";
            std::cout<<"[DR] prevPose.cov_orientation\n" <<prevPose.cov_orientation<<"\n";
            #endif

            if (useTranforWithUncertainty)
            {
                /** Prev Pose in TF form **/
                tfPrevPose = prevPose;

                /** Create the transformation from the delta position and the actual position **/
                tfDeltaPose = deltaPose;

                /** To perform the transformation **/
                tfPostPose = tfPrevPose * tfDeltaPose;

                /** Fill the rigid body state **/
                tfPostPose.copyToRigidBodyState(postPose);
            }
            else
            {
                /** Fill the instantaneous velocities in the rigid body state **/
                postPose.position += prevPose.orientation * deltaPose.position;
                postPose.cov_position += deltaPose.cov_position;
                postPose.orientation = prevPose.orientation * deltaPose.orientation;
                postPose.cov_orientation += deltaPose.cov_orientation;
            }

            /** Update the velocity information in the posterior **/
            postPose.velocity = cartesianVelocities[0].block<NUMAXIS, 1> (0,0); //!Velocity in body frame
            postPose.cov_velocity = cartesianVelCov.block<NUMAXIS, NUMAXIS> (0,0);
            postPose.angular_velocity = cartesianVelocities[0].block<NUMAXIS, 1> (NUMAXIS, 0);
            postPose.cov_angular_velocity = cartesianVelCov.block<NUMAXIS, NUMAXIS> (NUMAXIS, NUMAXIS);

            #ifdef DEAD_RECKON_DEBUG_PRINTS
            std::cout<<"[DR] Distance\n"<<postPose.orientation * cartesianVelocities[0].block<NUMAXIS, 1> (0,0) * delta_t << "\n";
            std::cout<<"[DR] postPose.position\n" <<postPose.position<<"\n";
            std::cout<<"[DR] postPose.cov_position\n" <<postPose.cov_position<<"\n";
            std::cout<<"[DR] postPose.velocity\n" <<postPose.velocity<<"\n";
            std::cout<<"[DR] postPose.cov_velocity\n" <<postPose.cov_velocity<<"\n";
            euler[2] = postPose.orientation.toRotationMatrix().eulerAngles(2,1,0)[0];//Yaw
            euler[1] = postPose.orientation.toRotationMatrix().eulerAngles(2,1,0)[1];//Pitch
            euler[0] = postPose.orientation.toRotationMatrix().eulerAngles(2,1,0)[2];//Roll
            std::cout<<"[DR] post Roll: "<<euler[0]*localization::R2D<<" Pitch: "<<euler[1]*localization::R2D<<" Yaw: "<<euler[2]*localization::R2D<<"\n";
            std::cout<<"[DR] postPose.cov_orientation\n" <<postPose.cov_orientation<<"\n";
            #endif

            return deltaPose;
        }

        /** \Brief Perform Quaternion integration assuming constant angular acceleration
        *
        * @return delta quaternion of the change in Attitude.
        *
        * */
        static Eigen::Quaternion<double> updateAttitude (const double dt,
                                const std::vector< Eigen::Matrix < double, NUMAXIS, 1> , Eigen::aligned_allocator < Eigen::Matrix <double, NUMAXIS, 1> > > &angularVelocities,
                                const Eigen::Matrix <double, NUMAXIS, NUMAXIS> &angularVelCov)
        {
            Eigen::Quaternion<double> deltaq; /** Instantaneous change in attitude **/
            Eigen::Matrix <double,QUATERSIZE,QUATERSIZE> omega4, oldomega4; /** Angular velocity matrix */
            Eigen::Matrix <double,QUATERSIZE,1> quat; /** Quaternion integration matrix */

            quat<< 1.00, 0.00, 0.00, 0.00; /**Identity quaternion */

            /** Discrete quaternion integration of the angular velocity **/
            omega4 << 0,-angularVelocities[0](0), -angularVelocities[0](1), -angularVelocities[0](2),
	        angularVelocities[0](0), 0, angularVelocities[0](2), -angularVelocities[0](1),
	        angularVelocities[0](1), -angularVelocities[0](2), 0, angularVelocities[0](0),
	        angularVelocities[0](2), angularVelocities[0](1), -angularVelocities[0](0), 0;

            oldomega4 << 0,-angularVelocities[1](0), -angularVelocities[1](1), -angularVelocities[1](2),
	        angularVelocities[1](0), 0, angularVelocities[1](2), -angularVelocities[1](1),
	        angularVelocities[1](1), -angularVelocities[1](2), 0, angularVelocities[1](0),
	        angularVelocities[1](2), angularVelocities[1](1), -angularVelocities[1](0), 0;


            /** Quaternion integration (third order linearization) **/
            quat = (Eigen::Matrix<double,QUATERSIZE,QUATERSIZE>::Identity() + (0.75 * omega4 *dt)- (0.25 * oldomega4 * dt) -
            ((1/6) * angularVelocities[0].squaredNorm() * pow(dt,2) *  Eigen::Matrix<double,QUATERSIZE,QUATERSIZE>::Identity()) -
            ((1/24) * omega4 * oldomega4 * pow(dt,2)) - ((1/48) * angularVelocities[0].squaredNorm() * omega4 * pow(dt,3))) * quat;

            /** Store in a quaternion form **/
            deltaq.w() = quat(0);
            deltaq.x() = quat(1);
            deltaq.y() = quat(2);
            deltaq.z() = quat(3);
            deltaq.normalize();

            #ifdef DEAD_RECKON_DEBUG_PRINTS
            std::cout <<"[DR] deltaq is "<<deltaq.w()<<" "<<deltaq.x()<<" "<<deltaq.y()<<" "<<deltaq.z()<<"\n";
            #endif


            return deltaq;
        }



    };

}

#endif
