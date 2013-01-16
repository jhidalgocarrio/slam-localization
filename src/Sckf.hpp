/**\file Sckf.hpp
 * Header function file and defines
 */

#ifndef _SCKF_HPP_
#define _SCKF_HPP_

#include <iostream> /** IO C++ Standard library */
#include <algorithm> /** Algorithm C++ Standard library */
#include <Eigen/LU> /** Lineal algebra of Eigen */
#include <Eigen/SVD> /** Singular Value Decomposition (SVD) of Eigen */
#include <Eigen/Geometry> /** Eigen data type for Matrix, Quaternion, etc... */
#include "Configuration.hpp" /** For the localization framework constant and configuration values **/
#include "Measurement.hpp" /** For the Measurement Generation module of the framework **/

namespace localization	
{
//     #ifndef EIGEN_NO_AUTOMATIC_RESIZING
//     #define EIGEN_NO_AUTOMATIC_RESIZING
//     #endif
    
    class Sckf
    {
	
    public:
	
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW
	
	/** CONSTANT VALUES TO THE CLASS**/
	
	/** Constants for the process model **/
	static const int A_STATE_VECTOR_SIZE = 9; /** Attitude state vector error (3 for attitude error, 3 for gyro bias and 3 for acc bias) **/
	static const int X_STATE_VECTOR_SIZE = (2*NUMAXIS + A_STATE_VECTOR_SIZE); /** Position, velocity and attitude state vector error **/

	/** Constants for the measurement model **/
	static const int X_OBSERVATION_VECTOR_SIZE = X_STATE_VECTOR_SIZE; /** Observation vector for the correction of vector state **/
	static const int X_MEASUREMENT_VECTOR_SIZE = SLIP_VECTOR_SIZE + (2*NUMAXIS); /** Measurement vector for computation of the correction measurement **/
	
	/** MEASUREMENT GENERATION **/
	
	/** Object of the measurement class **/
	Measurement filtermeasurement;
		
    private:
	
	/** FILTER VARIABLES (UPPER CASE MATRICES, LOWER CASE VECTORS) **/
	Eigen::Matrix <double,Eigen::Dynamic,1> xk_k; /** State vector Xk|k at the lastest extereoceptive measurement recorded */
	Eigen::Matrix <double,Eigen::Dynamic,1> xki_k; /** State vector Xk+i|k at the lastest proprioceptive measurement recorded (Robot's current state) */
	Eigen::Quaternion <double> q4;  /** Current robot attitude quaternion (integration) */
	Eigen::Quaternion <double> prev_q4;  /** Previous robot attitude quaternion (integration) (for computing the delta quaternion) */
	Eigen::Matrix <double,QUATERSIZE, QUATERSIZE> oldomega4; /** Quaternion integration matrix */
	
	/** System matrix **/
	Eigen::Matrix <double,Eigen::Dynamic,Eigen::Dynamic> A; /** Attitude System matrix */
	Eigen::Matrix <double,Eigen::Dynamic,Eigen::Dynamic> Fki; /** System matrix associated to vector Xk+i|k*/
	
	/** Error process matrix **/
	Eigen::Matrix <double,Eigen::Dynamic,Eigen::Dynamic> Qk; /** System Process noise convariance matrix of vector Xk+i|k and therefore Xk|k*/
	
	/** System covariance matrices **/
	Eigen::Matrix <double,Eigen::Dynamic,Eigen::Dynamic> Pk_k; /** Error covariance matrix of vector Xk|k*/
	Eigen::Matrix <double,Eigen::Dynamic,Eigen::Dynamic> Pki_k; /** Error covariance matrix of vector Xk+i|k*/
	
	/** Kalman gain matrix **/
	Eigen::Matrix <double,Eigen::Dynamic,Eigen::Dynamic> K; /** Kalman gain associated to xki_k*/
	
	/** Measurement observation matrices **/
	Eigen::Matrix <double,Eigen::Dynamic,Eigen::Dynamic> Hk; /** Measurement Observation matrix for xki_k */
	Eigen::Matrix <double,NUMAXIS,Eigen::Dynamic> H1a; /** Measurement 1 Observation matrix for attitude */
	Eigen::Matrix <double,NUMAXIS,Eigen::Dynamic> H2a; /** Measurement 2 Observation matrix for attitude */
	
	/** Measurement vector **/
	Eigen::Matrix <double,Eigen::Dynamic,1> zki; /** Measurement vector assodicate at time k+i **/
	
	/** Error measurement matrices **/
	Eigen::Matrix <double,NUMAXIS,NUMAXIS*M1> RHist; /** History of M1 measurement noise convariance matrix (for the adaptive algorithm) */
	Eigen::Matrix <double,NUMAXIS,NUMAXIS> Rg; /** Measurement and system noise convariance matrix for gyros (gyros are used in both: predict and update) */
	Eigen::Matrix <double,NUMAXIS,NUMAXIS> Ra; /** Measurement noise convariance matrix for acc */
	Eigen::Matrix <double,NUMAXIS,NUMAXIS> Rat; /** Measurement noise convariance matrix for the attitude (the estimation of the gravity vector) */
	Eigen::Matrix <double,NUMAXIS,NUMAXIS> Rm; /** Measurement noise convariance matrix for mag */
	Eigen::Matrix <double,Eigen::Dynamic,Eigen::Dynamic> Rk; /** Measurement noise convariance matrix for the measurement vector of proprioceptive sensors */
	
	Eigen::Matrix <double,Eigen::Dynamic,1> innovation; /** Internal variable of the filter innovation (zki - Hk*xki_k **/
	
	/** Tilde position **/
	Eigen::Matrix <double,NUMAXIS,1> eposition; /** Position error */
	
	/** Tilde velocity **/
// 	Eigen::Matrix <double,NUMAXIS,1> evelocity; /** Velocity error */
	DataModel evelocity; /** Velocity error */
	
	/** For the attitude computation **/
	Eigen::Matrix <double,NUMAXIS,1> gtilde; /** gravitation acceleration in world frame */
	Eigen::Matrix <double,NUMAXIS,1> mtilde; /** Magnetic dip angle in world frame */
	Eigen::Matrix <double,NUMAXIS,1> bghat; /** Estimated bias for gyroscope */
	Eigen::Matrix <double,NUMAXIS,1> bahat; /** Estimated bias for accelerometer */
	
	/** Adaptive external acceleration variables **/
	unsigned int r1count; /** Variable used in the adaptive algorithm, to compute the Uk matrix for SVD*/
	double r2count; /** Variable used in the adaptive algorithm, to compute the final Qstart cov. matrix*/
	
	/** Vel model at k-1 **/
	DataModel veloModelk_1, veloTruth;
	DataModel increVeloError;
	
    public:

	/**
	* Print a welcome to stdout
	* \return nothing
	*/
	void welcome();
	
	
	/**
	* @brief Gets the current state vector (x) of the filter
	* 
	* @author Javier Hidalgo Carrio.
	*
	* @return State Vector
	*
	*/
	Eigen::Matrix <double,Eigen::Dynamic,1> getStatex();
	
	
	/**
	* @brief Gets the current orientation in Quaternion
	* 
	* @author Javier Hidalgo Carrio.
	*
	* @return Quaternion with the current orientation.
	*
	*/
	Eigen::Quaternion <double> getAttitude();
	
	/**
	* @brief Gets the gravity value being used by the filter
	* 
	* @author Javier Hidalgo Carrio.
	*
	* @return Gravity value
	*
	*/
	 double getGravity();
	
	
	/**
	* @brief Gets the current orientation in Euler angles
	* 
	* @author Javier Hidalgo Carrio.
	*
	* @return Current orientation in Euler angles.
	*
	*/
	Eigen::Matrix <double, NUMAXIS, 1> getEuler();
	
	/**
	* @brief Gets Noise covariance matrix associated at vector-state x
	* 
	* @author Javier Hidalgo Carrio.
	*
	* @return Matrix P of the covariance of the state vector
	*
	*/
	Eigen::Matrix <double,Eigen::Dynamic,Eigen::Dynamic> getCovariancex();
	
	/**
	* @brief Gets Noise covariance matrix attitude part
	* 
	* @author Javier Hidalgo Carrio.
	*
	* @return Submatrix  of P 
	*
	*/
	Eigen::Matrix <double,A_STATE_VECTOR_SIZE,A_STATE_VECTOR_SIZE> getCovarianceAttitude();
	
	/**
	* @brief Return the filter Kalman Gain
	*
	* @author Javier Hidalgo Carrio.
	*
	* @return the filter kalman gain
	*
	*/
	Eigen::Matrix <double,Eigen::Dynamic, Eigen::Dynamic> getKalmanGain();
	
	/**
	* @brief Return the Kalman Gain of the attityde part 
	*
	* @author Javier Hidalgo Carrio.
	*
	* @return the attitude part of kalman gain
	*
	*/
	Eigen::Matrix <double,Eigen::Dynamic, Eigen::Dynamic> getAttitudeKalmanGain();
	
	/**
	* @brief Return the filter innovation (zki - Hk*xki_k) 
	*
	* @author Javier Hidalgo Carrio.
	*
	* @return the innovation vector
	*
	*/
	Eigen::Matrix <double,Eigen::Dynamic, 1> getInnovation();
	
	/**
	* @brief Gets the delta orientation in Quaternion
	* 
	* Delta quaternion is the ortation difference between
	* the previous attitude and the curret one.
	* 
	* @author Javier Hidalgo Carrio.
	*
	* @return Quaternion with the delta orientation (q(k-1) to q(k)) .
	*
	*/
	Eigen::Quaternion <double> deltaQuaternion();
	
	/**
	* @brief This function Initialize Attitude
	* 
	* Initial orientation value before start the IKF 
	*
	* @author Javier Hidalgo Carrio.
	*
	* @param[in] initq quaternion with the initial orientation
	*
	* @return OK is everything all right. ERROR on other cases.
	*
	*/
	int setAttitude (Eigen::Quaternion <double> &initq);
	
	/**
	* @brief Set gravity value
	* 
	* Set the gravity constant value to be used by the filter
	*
	* @author Javier Hidalgo Carrio.
	*
	* @param[in] g gravity value.
	*
	* @return nothing
	*
	*/
	void setGravity (double g);
	
	/**
	* @brief This function Initialize the State vector
	* 
	* The state vector is formed by defaulst for 16 + 9 element.
	* 16 elelements that governt the rover wheels interaction with the ground
	* (0-2) -> 3D slip vector of the contact current point
	* (3)   -> point contact angle
	* 
	* This 4 elements by 4 wheels of the rover (This is variable by the NUMBER_OF_WHEELS static
	* const variable).
	* 
	* 9 elements for the attitude estimation
	* (0-2) -> the vector patr of a error quaternion
	* (3-5) -> gyroscope bias estimation
	* (6-8) -> accelerometer bias estimation
	*
	* @param[in] x_0 a initial/desired state vector
	*
	* @return OK is everything all right. ERROR on other cases.
	*
	*/
	void setStatex (Eigen::Matrix <double,Eigen::Dynamic,1> &x_0);
	
	
	/**
	* @brief This function Initialize the Bias offset
	* 
	* Initial value of the gyros and acc bias offset
	*
	* @param[in] 
	*
	* @return OK is everything all right. ERROR on other cases.
	*
	*/
	void setBiasOffset (Eigen::Matrix <double,NUMAXIS,1> gbias, Eigen::Matrix <double,NUMAXIS,1> abias);
	
	/**
	* @brief This function set the initial Omega matrix
	* 
	* Initial Omega matrix with angular velocity for 
	* quaternion integration.
	*
	* @author Javier Hidalgo Carrio.
	*
	* @param[in] u pointer to vector with the angular velocity
	*
	* @return OK is everything all right. ERROR on other cases.
	*
	*/
	int setOmega (Eigen::Matrix <double,NUMAXIS,1>  &u);
	
	/**
	* @brief Set Yaw angle
	* 
	* Set the heading(yaw) angle from an external source
	*
	* @author Javier Hidalgo Carrio.
	*
	* @param[in] yaw yaw value in radians
	*
	* @return OK is everything all right. ERROR on other cases.
	*
	*/
	void setHeading (double yaw);
	
	/**
	* @brief This function Initilize the vectors and matrix of the Filter
	* 
	* This method receives the measurement noise matrix of the sensors
	* The theoretical gravity value and the Dip angle of the location.
	*
	* @author Javier Hidalgo Carrio.
	*
	* @param[in] P_0 Initial state covariance matrix
	* @param[in] Rg process noise matrix of Gyroscopes.
	* @param[in] Qbg process noise matrix of the gyroscopes bias (bias instability)
	* @param[in] Qba process noise matrix of the accelerometers bias (bias instability)
	* 
	*
	* @param[in] Ren measurement noise matrix of the encoders
	* @param[in] Rcontact measurement noise matrix of the contact angle information
	* @param[in] Ra measurement noise matrix of Accelerometers.
	* @param[in] Rm measurement noise matrix of Magnetometers.
	*
	* @param[in] g local gravitational value.
	* @param[in] alpha Dip angle
	*
	* @return void
	*
	*/
	void Init(Eigen::Matrix< double, Eigen::Dynamic, Eigen::Dynamic >& P_0,
		Eigen::Matrix< double, NUMAXIS , NUMAXIS  >& Rg,
		Eigen::Matrix< double, NUMAXIS , NUMAXIS  >& Qbg,
		Eigen::Matrix< double, NUMAXIS , NUMAXIS  >& Qba,
		Eigen::Matrix< double, NUMAXIS , NUMAXIS  >& Ra,
		Eigen::Matrix< double, NUMAXIS , NUMAXIS  >& Rat,
		Eigen::Matrix< double, NUMAXIS , NUMAXIS  >& Rm,
		double g, double alpha);
	
	/**
	* @brief Performs the prediction step of the filter.
	* 
	* It computes the discrete version of the matrix F to propagate forward
	* the state vector x. It computes the Q and P matrix as well as the 
	* quaternion integration from the input vector u and the delta time.
	*
	* @author Javier Hidalgo Carrio.
	*
	* @param[in] u vector with the angular velocity
	* @param[in] v vector with the linear acceleration
	* @param[in] dt delta time between samples
	*
	* @return void
	*
	*/
	void predict(Eigen::Matrix <double,NUMAXIS,1>  &u, Eigen::Matrix <double,NUMAXIS,1>  &v, double dt);
	
	/**
	* @brief Performs the measurement and correction steps of the filter.
	* 
	* The IKf: 
	*  1. Measurement of the system (velocity, encoders, etc...)
	*  2. Measurement step to correct Pitch and Roll from accelerometers.	
	*  3. Measurement step to correct Yaw angle from magnetometers (optional).
	* 
	* The first measurement step is dynamics. The noise covariamce matrix
	* of the update is dynamic depending on external accelerations felt on
	* the accelerometers. That means the variance noise increase or decrease
	* depending on the external acceleration. Thas is the main different between 
	* normal EKF.
	* 
	*
	* @author Javier Hidalgo Carrio.
	*
	* @param[in] He observation submatrix of the vector xk+1|k related to the encoders.
	* @param[in] Hme measurement submatrix of the measurement vector z related to the encoders.
	* @param[in] slip_error vector with the velocity error in slip.
	* @param[in] acc vector with accelerations
	* @param[in] gyro vector with angular velocities
	* @param[in] magn vector with magnetometers
	* @param[in] magn_on_off boolean value to connect or disconnect the magnetometers correction
	* @return void
	*
	*/
	void update(Eigen::Matrix <double,NUMAXIS,NUMAXIS> &Hme, Eigen::Matrix <double,NUMAXIS,NUMAXIS> &Rme,
		  Eigen::Matrix< double, NUMAXIS, 1  >& slip_error,
		  Eigen::Matrix< double, NUMAXIS , 1  >& acc,Eigen::Matrix< double, NUMAXIS , 1  >& mag, double dt, bool magn_on_off);
	
	
	/**
	* @brief It calls the measurement class to perform the proprioceptive measurement.
	* 
	* It computes the nav kinematics and slip kinematics
	* 
	*/
	void measurementGeneration (const Eigen::Matrix< double, Eigen::Dynamic, Eigen::Dynamic > &Anav, const Eigen::Matrix< double, Eigen::Dynamic, Eigen::Dynamic > &Bnav,
				    const Eigen::Matrix< double, Eigen::Dynamic, Eigen::Dynamic > &Aslip, const Eigen::Matrix< double, Eigen::Dynamic, Eigen::Dynamic > &Bslip,
				    const Eigen::Matrix< double, Eigen::Dynamic, 1  > &vjoints, Eigen::Matrix< double, SLIP_VECTOR_SIZE, 1> &slip_error,
				    Eigen::Matrix< double, SLIP_VECTOR_SIZE, SLIP_VECTOR_SIZE> &slip_errorCov, double dt);
	
	/**
	* @brief It calls the measurement class to perform the proprioceptive measurement.
	* 
	* It computes the nav kinematics and incremental velocity error
	* 
	*/
	void measurementGeneration(const Eigen::Matrix< double, Eigen::Dynamic, Eigen::Dynamic >& Anav, const Eigen::Matrix< double, Eigen::Dynamic, Eigen::Dynamic >& Bnav,
				 const Eigen::Matrix< double, Eigen::Dynamic, 1  > &vjoints, Eigen::Matrix< double, NUMAXIS, 1> &velo_error,
				 Eigen::Matrix< double, NUMAXIS, NUMAXIS> &vel_errorCov, double dt);
	
	/**
	* @brief It resets the state vector error
	* 
	* It reset the part of the error state that are independent for the next propagation
	* 
	*/
	void resetStateVector ();
	
	
	DataModel getVeloTruth ();
	
	DataModel getVeloError ();
	
	DataModel getIncreVeloError();
	
	/**
	* @brief Save the current filter status
	* 
	* Save the current filter status in FilterInfo form
	* 
	*/
	void toFilterInfo(localization::FilterInfo &finfo);
	 
    };
    
    
} // end namespace localization

#endif // 
