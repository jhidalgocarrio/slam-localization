/**\file DeadReckon.cpp
 *
 * This class has the methods to compute a Dead Reckoning process
 * which inputs are average velocities in the sampling interval
 * 
 * @author Javier Hidalgo Carrio | DFKI RIC Bremen | javier.hidalgo_carrio@dfki.de
 * @date December 2012.
 * @version 1.0.
 */

#include "DeadReckon.hpp" 


using namespace localization;

/** \Brief Default Constructor
*/
DeadReckon::DeadReckon()
{
    base::samples::RigidBodyState rbsInit;
    
    rbsInit.invalidate();
    rbsInit.sourceFrame = "Body Frame";
    rbsInit.targetFrame = "Geographic_Frame (North-West-Up)";
    
    *this = DeadReckon(rbsInit);
}



/** \Brief Constructor
 */
DeadReckon::DeadReckon(base::samples::RigidBodyState initPose)
{
    vState.setZero();
    
    /** Set the initial pose in the uncertainty variable **/
    this->actualPose = initPose;
}


/** \Brief Performs the time integration of the delta pose updates.
*/
base::samples::RigidBodyState DeadReckon::updatePose(Eigen::Matrix< double, NUMAXIS , 1  > linvelo, Eigen::Quaternion <double> delta_q,
						     Eigen::Matrix< double, NUMAXIS , NUMAXIS  > covlinvelo, Eigen::Matrix< double, NUMAXIS , NUMAXIS  > covdelta_q,
						     Eigen::Matrix< double, NUMAXIS , 1  > linvelo_error, Eigen::Quaternion <double> delta_qe,
						     Eigen::Matrix< double, NUMAXIS , NUMAXIS  > covlinvelo_error, Eigen::Matrix< double, NUMAXIS , NUMAXIS  > covdelta_qe,
						     base::Time timeStamp, double delta_t)
{
    base::samples::RigidBodyState rbsDeltaPose, rbsBC;
    envire::TransformWithUncertainty deltaPose;
    
    /** Prepare the vState variable for the dead reckoning **/
    vState.block<NUMAXIS,1>(0,1) = vState.block<NUMAXIS,1>(0,0); // move the previous state to the col(1)

    /** Update the internal variables with the new measurements **/
    vState.block<NUMAXIS,1>(0,0) = linvelo; // x,y and z

    /** Calculate the delta position from velocity (dead reckoning) asuming constant acceleration **/
    rbsDeltaPose.position = ((delta_t/2.0) * (vState.block<3,1>(0,1) + vState.block<3,1>(0,0)));
    
    /** Create the transformation from the delta position and the actual position **/
    deltaPose = rbsDeltaPose;
    
    #ifdef DEBUG_PRINTS
    std::cout<<"[DR] actualPose(before)\n" <<actualPose;
    #endif
    
    /** To perform the transformation **/
    actualPose = actualPose * deltaPose;
    
    /** Fill the rigid body state **/
    rbsBC.time = timeStamp;
    actualPose.copyToRigidBodyState(rbsBC);
    rbsBC.velocity = vState.block<3,1>(0,0);
    
    #ifdef DEBUG_PRINTS
    std::cout<<"[DR] Delta_t"<<delta_t<< "\n";
    std::cout<<"[DR] Distance\n"<<rbsBC.orientation * vState.block<3,1>(0,0) * this->delta_t << "\n";
    #endif
    
    return rbsBC;
}
