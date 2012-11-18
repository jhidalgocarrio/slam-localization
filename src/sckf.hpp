/**\file sckf.hpp
 * Header function file and defines
 */

#ifndef _SCKF_HPP_
#define _SCKF_HPP_

#include <iostream>
#include <Eigen/Geometry> /** Eigen data type for Matrix, Quaternion, etc... */


namespace localization	
{
//     #ifndef EIGEN_NO_AUTOMATIC_RESIZING
//     #define EIGEN_NO_AUTOMATIC_RESIZING
//     #endif
    
    /** General defines **/
    #ifndef OK
    #define OK	0  /** Integer value in order to return when everything is all right. */
    #endif
    #ifndef ERROR
    #define ERROR -1  /** Integer value in order to return when an error occured. */
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
    #define M1 5 /**< Parameter for adaptive algorithm */
    #define M2 3 /**< Parameter for adaptive algorithm */
    #define GAMMA 0.1 /**< Parameter for adaptive algorithm */
    #define R2COUNT 100 /**< Parameter for adaptive algorithm */

    #ifndef D2R
    #define D2R M_PI/180.00 /** Convert degree to radian **/
    #endif
    
    #ifndef R2D
    #define R2D 180.00/M_PI /** Convert radian to degree **/
    #endif
    
    class sckf
    {
	
    public:
	
	/** CONSTANT VALUES TO THE CLASS**/
	static const int NUMBER_OF_WHEELS = 4;
	static const int E_STATE_VECTOR_SIZE = (NUMAXIS + 1); /** Rover position state vector error (3 elements for slip + 1 for contact angle) **/
	static const int A_STATE_VECTOR_SIZE = 9; /** Attitude state vector error **/
	static const int X_STATE_VECTOR_SIZE = ((E_STATE_VECTOR_SIZE*NUMBER_OF_WHEELS) + A_STATE_VECTOR_SIZE); /** State vector error **/
	static const int Y_MEASUREMENT_VECTOR_SIZE = ((NUMBER_OF_WHEELS+1)+(2*NUMAXIS)); /** Measurement vector for computation of the position (5 is 4 motor encoders + 1 passive joint) **/
	static const int E_MEASUREMENT_VECTOR_SIZE = (sckf::NUMBER_OF_WHEELS*(2*NUMAXIS)); /** Measurement vector for the correction of the rover position state error **/
	static const int Z_MEASUREMENT_VECTOR_SIZE = (E_MEASUREMENT_VECTOR_SIZE + NUMAXIS); /** Whole rover measurement vector (6D for rover velocities + 4 of wheels encoders + 1 of passive joint) + IMU sensor **/

	
    private:
	
	/** FILTER VARIABLES (UPPER CASE MATRICES LOWER CASE VECTORS) **/
	Eigen::Matrix <double,Eigen::Dynamic,1> xk_k; /** State vector Xk|k at the lastest extereoceptive measurement recorded */
	Eigen::Matrix <double,Eigen::Dynamic,1> xki_k; /** State vector Xk+i|k at the lastest proprioceptive measurement recorded (Robot's current state) */
	Eigen::Quaternion <double> q4;  /** Current robot attitude quaternion (integration) */
	Eigen::Matrix <double,QUATERSIZE, QUATERSIZE> oldomega4; /** Quaternion integration matrix */
	
	/** System matrix **/
	Eigen::Matrix <double,Eigen::Dynamic,Eigen::Dynamic> A; /** Attitude System matrix */
	Eigen::Matrix <double,Eigen::Dynamic,Eigen::Dynamic> Fki; /** System matrix associated to vector Xk+i|k*/
	
	/** Error process matrix **/
	Eigen::Matrix <double,Eigen::Dynamic,Eigen::Dynamic> Qk; /** System Process noise convariance matrix of vector Xk+i|k and therefore Xk|k*/
	
	/** System covariance matrices **/
	Eigen::Matrix <double,Eigen::Dynamic,Eigen::Dynamic> Pk_k; /** Error covariance matrix of vector Xk|k*/
	Eigen::Matrix <double,Eigen::Dynamic,Eigen::Dynamic> Pki_k; /** Error covariance matrix of vector Xk+i|k*/
	
	/** Measurement observation matrices **/
	Eigen::Matrix <double,NUMAXIS,Eigen::Dynamic> H1a; /** Measurement 1 Observation matrix for attitude */
	Eigen::Matrix <double,NUMAXIS,Eigen::Dynamic> H2a; /** Measurement 2 Observation matrix for attitude */
	Eigen::Matrix <double,Eigen::Dynamic,Eigen::Dynamic> Hk; /** Measurement Observation matrix for the vector Xk+i|k and therefore Xk|k */
	
	/** Measurement vector **/
	Eigen::Matrix <double,Eigen::Dynamic,1> zki; /** Measurement vector assodicate at time k+i **/
	
	/** Error measurement matrices **/
	Eigen::Matrix <double,NUMAXIS,NUMAXIS*M1> RHist; /** History of M1 measurement noise convariance matrix (for the adaptive algorithm) */
	Eigen::Matrix <double,NUMAXIS,NUMAXIS> Rv; /** Measurement noise convariance matrix for linear velocities */
	Eigen::Matrix <double,NUMAXIS,NUMAXIS> Rg; /** Measurement and system noise convariance matrix for gyros (gyros are used in both: predict and update) */
	Eigen::Matrix <double,NUMAXIS,NUMAXIS> Ra; /** Measurement noise convariance matrix for acc */
	Eigen::Matrix <double,Eigen::Dynamic,Eigen::Dynamic> Ren; /** Measurement noise convariance matrix for joint and motor encoders */
	Eigen::Matrix <double,NUMAXIS,NUMAXIS> Rm; /** Measurement noise convariance matrix for mag */
	Eigen::Matrix <double,Eigen::Dynamic,Eigen::Dynamic> Rk; /** Noise convariance matrix for the measurement vector of the filter */
	Eigen::Matrix <double,Eigen::Dynamic,1> innovation; /** Internal variable of the filter innovation (zki - Hk*xki_k **/
	
	/** Kalman Gain matrix **/
	Eigen::Matrix <double,Eigen::Dynamic,Eigen::Dynamic> K; /** Kalman gain associted to the vector Xk+i|k */
	
	/** For the slip kinematics (each column is a wheel defined by a wheel_idx) **/
	Eigen::Matrix <double,NUMAXIS,NUMBER_OF_WHEELS> slipMatrix;
	
	/** For the contact angle (information of the contact angles, each row is a wheel) **/
	Eigen::Matrix <double,NUMBER_OF_WHEELS, 1> acontact;
	
	/** Linear velocities computed from acc information **/
	Eigen::Matrix <double, NUMAXIS, 1> linvelocity, angvelocity;
	
	/** For the attitude computation **/
	Eigen::Matrix <double,NUMAXIS,1> gtilde; /** gravitation acceleration in world frame */
	Eigen::Matrix <double,NUMAXIS,1> mtilde; /** Magnetic dip angle in world frame */
	Eigen::Matrix <double,NUMAXIS,1> bghat; /** Estimated bias for gyroscope */
	Eigen::Matrix <double,NUMAXIS,1> bahat; /** Estimated bias for accelerometer */
	
	unsigned int r1count; /** Variable used in the adaptive algorithm, to compute the Uk matrix for SVD*/
	double r2count; /** Variable used in the adaptive algorithm, to compute the final Qstart cov. matrix*/
	Eigen::Matrix <double,NUMAXIS,1> eccx, eccy, eccz; /** Accelerometers excentricity with respect to the body center of the robot **/
	
	
	
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
	* @brief Gets Linear velocities
	* 
	* @author Javier Hidalgo Carrio.
	*
	* @return Subvector of ye
	*
	*/
	Eigen::Matrix <double,NUMAXIS,1> getLinearVelocities();
	
	/**
	* @brief Gets Angular velocities
	* 
	* @author Javier Hidalgo Carrio.
	*
	* @return Subvector of ye
	*
	*/
	Eigen::Matrix <double,NUMAXIS,1> getAngularVelocities();
	
	/**
	* @brief Gets 3x1 slip vector
	* 
	* @param[in] wheel_idx wheel index between 0 and NUMBER_OF_WHEELS-1
	* 
	* @author Javier Hidalgo Carrio.
	*
	* @return the slip vector
	*
	*/
	Eigen::Matrix <double,NUMAXIS,1> getSlipVector(int wheel_idx);
	
	/**
	* @brief Gets contact angle vector
	*
	* @author Javier Hidalgo Carrio.
	*
	* @return the contact angles
	*
	*/
	Eigen::Matrix <double,Eigen::Dynamic,1> getContactAngles();
	
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
	* @param[in] *x_0 a initial/desired state vector
	*
	* @return OK is everything all right. ERROR on other cases.
	*
	*/
	void setStatex (Eigen::Matrix <double,Eigen::Dynamic,1> &x_0);
	
	/**
	* @brief This function set the initial Omega matrix
	* 
	* Initial Omega matrix with angular velocity for 
	* quaternion integration.
	*
	* @author Javier Hidalgo Carrio.
	*
	* @param[in] *u pointer to vector with the angular velocity
	*
	* @return OK is everything all right. ERROR on other cases.
	*
	*/
	int setOmega (Eigen::Matrix <double,NUMAXIS,1>  &u);
	
	/**
	* @brief This function set the Accelerometers excentricity
	* 
	* Here the eccentricity is defined as the vector of the distance
	* between the Body center and the accelerometer of the IMU.
	*
	* @author Javier Hidalgo Carrio.
	*
	* @param[in] eccx vector of the distance in meters for Accelerometers X axis
	* @param[in] eccy vector of the distance in meters for Accelerometers Y axis
	* @param[in] eccz vector of the distance in meters for Accelerometers Z axis
	*
	* @return OK is everything all right. ERROR on other cases.
	*
	*/
	void setEccentricity (Eigen::Matrix <double,NUMAXIS,1>  &eccx, Eigen::Matrix <double,NUMAXIS,1>  &eccy, Eigen::Matrix <double,NUMAXIS,1>  &eccz);
	
	/**
	* @brief This function Initilize the vectors and matrix of the Filter
	* 
	* This method receives the measurement noise matrix of the sensors
	* The theoretical gravity value and the Dip angle of the location.
	*
	* @author Javier Hidalgo Carrio.
	*
	* @param[in] P_0 Initial state covariance matrix
	* @param[in] Qec process noise matrix of slip vector and contact angles.
	* @param[in] Qbg process noise matrix of the gyroscopes bias
	* @param[in] Qba process noise matrix of the accelerometers bias 
	* 
	* @param[in] Rg process noise matrix of Gyroscopes.
	*
	* @param[in] Rz measurement noise matrix of measurement vector z.
	* @param[in] Ra measurement noise matrix of Accelerometers.
	* @param[in] Rm measurement noise matrix of Magnetometers.
	*
	* @param[in] ecc Inertial sensor eccentricity location (x, y and z).
	* @param[in] g local gravitational value.
	* @param[in] alpha Dip angle
	*
	* @return void
	*
	*/
	void Init(Eigen::Matrix< double, Eigen::Dynamic, Eigen::Dynamic >& P_0, Eigen::Matrix< double, Eigen::Dynamic, Eigen::Dynamic >& Qec,
		Eigen::Matrix< double, NUMAXIS , NUMAXIS  >& Qbg, Eigen::Matrix< double, NUMAXIS , NUMAXIS  >& Qba,
		Eigen::Matrix< double, NUMAXIS, NUMAXIS >& Rv, Eigen::Matrix< double, NUMAXIS , NUMAXIS  >& Rg,
		Eigen::Matrix <double,Eigen::Dynamic,Eigen::Dynamic> &Ren,
		Eigen::Matrix< double, NUMAXIS , NUMAXIS  >& Ra, Eigen::Matrix< double, NUMAXIS , NUMAXIS  >& Rm,
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
	* @param[in] dt delta time between samples
	*
	* @return void
	*
	*/
	void predict(Eigen::Matrix <double,NUMAXIS,1>  &u, double dt);
	
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
	* @param[in] Be measurement submatrix of the measurement vector z related to the encoders.
	* @param[in] encoders vector with the measurement values from the encoders.
	* @param[in] acc vector with accelerations
	* @param[in] gyro vector with angular velocities
	* @param[in] magn vector with magnetometers
	* @param[in] magn_on_off boolean value to connect or disconnect the magnetometers correction
	* @return void
	*
	*/
	void update(Eigen::Matrix <double,Eigen::Dynamic,Eigen::Dynamic> &He, Eigen::Matrix <double,Eigen::Dynamic,Eigen::Dynamic> &Be,
		    Eigen::Matrix <double,Eigen::Dynamic,1>  &encoders, Eigen::Matrix <double,NUMAXIS,1>  &acc, Eigen::Matrix <double,NUMAXIS,1>  &gyro,
		    Eigen::Matrix <double,NUMAXIS,1>  &mag, double dt, bool magn_on_off);
	    
    };
    
    /**
    * @brief This computes the theoretical gravity value according to the WGS-84 ellipsoid earth model.
    *
    * @author Javier Hidalgo Carrio.
    *
    * @param[in] latitude double the latitude value in radian
    * @param[in] altitude double with the altitude value in meters
    *
    * @return double. the theoretical value of the local gravity
    *
    */
    double GravityModel(double latitude, double altitude);
    
    /**
    * @brief Substract the Earth rotation from the gyroscopes readout
    *
    * This function computes the substraction of the rotation of the Earth (EARTHW)
    * from the gyroscope values. This function uses quaternion of transformation from
    * the body to the geographic frame and the latitude in radians.
    *
    * @author Javier Hidalgo Carrio.
    *
    * @param[in, out] *u pointer to angular velocity
    * @param[in] *qb_g quaternion from body frame to geographic frame
    * @param[in] latitude location latitude angle in radians
    *
    * @return void
    *
    */
    void SubstractEarthRotation(Eigen::Matrix <double, NUMAXIS, 1> *u, Eigen::Quaternion <double> *qb_g, double latitude);
    
    /**
    * @brief Correct the magnetic declination of the North
    *
    * Magnetic North and geographic North (Ertah rotation axis)
    * are different depending on geograohic location according
    * to a Declination Map. The function correct this bias.
    * See: http://www.magnetic-declination.com for futher information
    * about the declination angle of your location.
    *
    * @author Javier Hidalgo Carrio.
    *
    * @param[in, out] *quat pointer to quaternion with the orientation 
    * @param[in] double magnetic declination angle in radians
    * @param[in] mode. EAST or WEST depending on the magnetic declination
    *
    * @return OK is everything all right. ERROR on other cases.
    *
    */
    int CorrectMagneticDeclination (Eigen::Quaternion <double> *quat, double magnetic_declination,  int mode);

} // end namespace localization

#endif // 
