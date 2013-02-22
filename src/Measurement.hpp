/**\file measurement.hpp
 * Header function file and defines
 */

#ifndef _MEASUREMENT_HPP_
#define _MEASUREMENT_HPP_

#include <iostream> /** IO C++ Standard library */
#include <vector>
#include <algorithm> /** Algorithm C++ Standard library */
#include <Eigen/Geometry> /** Eigen data type for Matrix, Quaternion, etc... */
#include <Eigen/Core> /** Core methods of Eigen implementation **/
#include <Eigen/Dense> /** for the algebra and transformation matrices **/
#include <Eigen/Cholesky> /** For the Cholesky decomposition **/
#include <boost/circular_buffer.hpp> /** For circular_buffer **/
#include "DataModel.hpp" /** For the quantities models **/
#include "Configuration.hpp" /** For the localization framework constant and configuration values **/
#include "DataTypes.hpp" /** Orogen compatible export data types **/

namespace localization	
{
    class Measurement
    {
	
    public:
	
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW
	
    private:

	Eigen::Matrix <double,ENCODERS_VECTOR_SIZE,ENCODERS_VECTOR_SIZE> Rencoders; /** Measurement noise convariance matrix for joint velocities */
	
	/** Linear and angular velocities (from IMU) **/
	Eigen::Matrix <double, NUMAXIS, 1> linvelocity;
	
	/** Sensed encoders velocities **/
	Eigen::Matrix <double, ENCODERS_VECTOR_SIZE, 1>encodersvelocity;
	
	/** Offset quaternion to set the uneven weight distribution of the rover **/
	Eigen::Quaternion <double> q_weight_distribution;
	
	/** Velocity model from navigation kinematics **/
	DataModel velModel, prevVelModel;
	
	/** Increment in velocity model from navigation kinematics **/
	DataModel increVelModel;
	
	/** For the contact angle **/
	DataModel acontact;
	
	/** Slip vector **/
	DataModel slipModel, slipVector, slipError;
	
	/** For the slip kinematics (each column is a wheel defined by a wheel_idx) **/
	Eigen::Matrix <double,NUMAXIS,NUMBER_OF_WHEELS> slipMatrix;
	
	/** Translation distances for the accelerometers w.r.t to rover body **/
	Eigen::Matrix <double,NUMAXIS,1> eccx, eccy, eccz; /** Accelerometers excentricity with respect to the body center of the robot **/
	
	/** Circular Vector of accelerations for integral method **/
	Eigen::Matrix <double,NUMAXIS,1> cbAcc, cbAngvelo;
	
	/** Array of past rover velocity model **/
	Eigen::Matrix <double,NUMAXIS,1> cbVelModel;
	
	/** Least squared error of the navigation kinematics **/
	double leastSquaredNavigation;
	
	/** Least squared error of the slip kinematics **/
	double leastSquaredSlip;
	
	boost::circular_buffer<double> testAcc;
	
    public:
	
	/** Measurement contructor
         */
        Measurement();
	
	/** Measurement default descontructor
         */
        ~Measurement();
	
	/**
	* Print a welcome to stdout
	* \return nothing
	*/
	void welcome();
	
	/**
	* @brief Initialization method
	* 
	* Class initialization.
	* 
	*/
	void Init ( Eigen::Matrix< double, ENCODERS_VECTOR_SIZE , ENCODERS_VECTOR_SIZE  >& Ren,
		    Eigen::Matrix< double, NUMBER_OF_WHEELS , NUMBER_OF_WHEELS  >& Rcont,
		    Eigen::Quaternion <double> q_weight);
	
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
	* @brief Set the current linear velocities
	* 
	* @author Javier Hidalgo Carrio.
	*
	* @return void
	*
	*/
	void setLinearVelocities(Eigen::Matrix <double,NUMAXIS,1> &linvelo);
	
	/**
	* @brief Gets Linear velocities
	* 
	* @author Javier Hidalgo Carrio.
	*
	* @return the linear velocities
	*
	*/
	Eigen::Matrix <double,NUMAXIS,1> getLinearVelocities();
	
	/**
	* @brief Set the current angular velocity
	* 
	* @author Javier Hidalgo Carrio.
	*
	* @return void
	*
	*/
	void setAngularVelocities(Eigen::Matrix <double,NUMAXIS,1> &angvelo);
	
	/**
	* @brief Gets Angular velocities
	* 
	* @author Javier Hidalgo Carrio.
	*
	* @return angular velocity
	*
	*/
	Eigen::Matrix <double,NUMAXIS,1> getAngularVelocities();
	
	/**
	* @brief Set the linear Acceleration
	* 
	* @author Javier Hidalgo Carrio.
	*
	* @return Linear acceleration
	*
	*/
	void setLinearAcceleration(Eigen::Matrix <double,NUMAXIS,1> &linacc);
	
	/**
	* @brief Gets Linear Acceleration
	* 
	* @author Javier Hidalgo Carrio.
	*
	* @return Linear acceleration
	*
	*/
	Eigen::Matrix <double,NUMAXIS,1> getLinearAcceleration();
	
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
	Eigen::Matrix <double,NUMAXIS,1> getSlipVector(const unsigned int wheel_idx);
	
	/**
	* @brief Gets rover slip vector
	* 
	* 
	* @author Javier Hidalgo Carrio.
	*
	* @return the complete slip vector
	*
	*/
	Eigen::Matrix <double,SLIP_VECTOR_SIZE,1> getSlipVector();
	
	/**
	* @brief Gets rover slip vector Covariance
	* 
	* 
	* @author Javier Hidalgo Carrio.
	*
	* @return the complete slip vector cov matrix
	*
	*/
	Eigen::Matrix <double,SLIP_VECTOR_SIZE, SLIP_VECTOR_SIZE> getSlipVectorCov();
	
	/**
	* @brief Gets rover slip error vector
	* 
	* Return the slip error vector for the complete number of wheels
	* 
	* @author Javier Hidalgo Carrio.
	*
	* @return the rover slip error vector
	*
	*/
	Eigen::Matrix <double,SLIP_VECTOR_SIZE,1> getSlipErrorVector();
	
	/**
	* @brief Gets rover slip error covariance matrix
	* 
	* Return the slip error covariance matrix for the complete number of wheels
	* 
	* @author Javier Hidalgo Carrio.
	*
	* @return the rover slip error covariance matrix
	*
	*/
	Eigen::Matrix <double,SLIP_VECTOR_SIZE, SLIP_VECTOR_SIZE> getSlipErrorVectorCovariance();
	
	/**
	* @brief Set the current contact angles
	* 
	* Set the contact angle for the current wheel/foot point in contact
	* 
	* @param[in] contact_angles the vector of contact angles
	* 
	*/
	void setContactAnglesVelocity(Eigen::Matrix<double, NUMBER_OF_WHEELS, 1> &contact_angles,
				      Eigen::Matrix<double, NUMBER_OF_WHEELS, NUMBER_OF_WHEELS> &Cov);
	
	/**
	* @brief Gets contact angle vector
	*
	* @author Javier Hidalgo Carrio.
	*
	* @return the contact angles
	*
	*/
	Eigen::Matrix <double,Eigen::Dynamic,1> getContactAnglesVelocity();
	
	/**
	* @brief Set the current encoders velocities
	* 
	* Set the encoders velocities for the kinematics
	* 
	* @param[in] vjoints the vector of joints encoders velocities
	* 
	*/
	void setEncodersVelocity(const Eigen::Matrix<double, Eigen::Dynamic, 1> &vjoints);
	
	/**
	* @brief Get the velocity from the odometry model
	* 
	* Rover velocity from pure odometry model (navigation kinematics)
	* 
	* @return current rover velocity from odometry
	* 
	*/
	Eigen::Matrix<double, NUMAXIS, 1 > getCurrentVeloModel();
	
	/**
	* @brief Get the covariance of the velocity from the odometry model
	* 
	* Covariance matrix of the estimated quantity
	* 
	* @return the covariance matrix
	* 
	*/
	Eigen::Matrix<double, NUMAXIS, NUMAXIS > getCurrentVeloModelCovariance();
	
	/**
	* @brief Set the current velocity
	* 
	* It stores the current rover velocity from the odometry model
	* 
	* @param[in] velocity  current rover velocity
	* 
	*/
	void setCurrentVeloModel(Eigen::Matrix<double, NUMAXIS, 1> &velocity);
	
	
	/**
	* @brief Get the increment in velocity
	* 
	* Rover increment velocity from odometry model (navigation kinematics)
	* 
	* @return current rover incremet in velocity from odometry
	* 
	*/
	Eigen::Matrix<double, NUMAXIS, 1 > getIncrementalVeloModel();
	
	/**
	* @brief Get the covariance of the increment in velocity
	* 
	* Covariance matrix of the estimated quantity
	* 
	* @return the covariance matrix
	* 
	*/
	Eigen::Matrix<double, NUMAXIS, NUMAXIS > getIncrementalVeloModelCovariance();
	
	/**
	 * @brief Returns the covariance noise matrix
	 */
	Eigen::Matrix <double,ENCODERS_VECTOR_SIZE,ENCODERS_VECTOR_SIZE> getEncodersVelocityCovariance();
	
	/**
	 * @brief Returns the covariance noise matrix
	 */
	Eigen::Matrix <double,NUMBER_OF_WHEELS,NUMBER_OF_WHEELS>  getContactAnglesVelocityCovariance();
	
	/**
	 * @brief Returns the quaternion with the projection for the nadir vector
	 */
	Eigen::Quaternion <double>  getLevelWeightDistribution();
	
	
	/**
	* @brief Simpson's rule for numerical integration
	* 
	* It computes Simpson's numerical integration method
	*
	* @author Javier Hidalgo Carrio.
	*
	* @param[in] fa sample at time t-2
	* @param[in] fm sample at time t-1
	* @param[in] fb sample at time t
	* @param[in] delta_ab time between samples b and a
	*
	* @return OK is everything all right. ERROR on other cases.
	*
	*/
	inline double simpsonsIntegral (double fa, double fm, double fb, double delta_ab);
	
	/**
	* @brief Perform the accelerometers integration
	* 
	* Integration of accelerometers for the window defined
	* in INTEGRATION_WINDOWS_SIZE
	*/
	Eigen::Matrix<double, NUMAXIS,1> accIntegrationWindow(double dt);
	
	/**
	* @brief Return the buffer size for x,y and z
	* 
	* Integration size of the windows for each coordinates
	*/
	Eigen::Matrix< double, NUMAXIS , 1  > getIntegrationWindowSize();
	
	/**
	* @brief Least-Squares navigation kinematics solution
	* 
	* Computes the rover velocity, slip vector due to non-holonomic constraints
	* and contact angle velocities.
	* 
	* @param[in] A the matrix with the non-sensed values
	* @param[in] B the matrix with the sensed values
	* @param[in] R covariance matrix of the sensed values.
	* 
	* @return relative least-squares error
	* 
	*/
	double navigationKinematics (const Eigen::Matrix< double, Eigen::Dynamic, Eigen::Dynamic > &A, const Eigen::Matrix< double, Eigen::Dynamic, Eigen::Dynamic > &B,
				     const Eigen::Matrix< double, Eigen::Dynamic, Eigen::Dynamic > &R, const Eigen::Matrix< double, Eigen::Dynamic, Eigen::Dynamic > &W);
	
	/**
	* @brief Least-Squares slip kinematics solution
	* 
	* Computes the rover slip vector for the number of wheels
	* due to loose of traction with the ground
	* 
	* 
	* @param[in] A the matrix with the non-sensed values
	* @param[in] B the matrix with the sensed values
	* @param[in] R covariance matrix of the sensed values.
	* @param[in] dt delta integration step for the acc
	* 
	* @return relative least-squares error
	* 
	*/
	double slipKinematics (const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> &A, const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> &B,
				const Eigen::Matrix< double, Eigen::Dynamic, Eigen::Dynamic > &R);
	
	
	/**
	* @brief Save slip info to the Orogen-compatible DataType
	* 
	* 
	* @param[out] sinfo the slip information
	* 
	* @return void
	* 
	*/
	void toSlipInfo (localization::SlipInfo &sinfo);
	
	
	/**
	* @brief Save some useful measurement generation info.
	* 
	* 
	* @param[out] measurementinfo some debung info the measurement generation
	* 
	* @return void
	* 
	*/
	void toMeasurementGenerationInfo (localization::MeasurementGenerationInfo &measurementInfo);

	
	/**
	* @brief Computes the Bhattacharyya coefficient
	* 
	* @return BC
	* 
	*/
	Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> bhattacharyya (DataModel &data1, DataModel &data2);
	
	/**
	* @brief Computes the Mahalanobis distance
	* 
	* @return MH
	* 
	*/
	double mahalanobis (DataModel &data1);
	
	/**
	* @brief Convert data value in the range of MinValues..MaxValues to the range 350 - 650
	*/
	double getWaveLenghtFromValue (const double value, const double max, const double min);

	/**
	* @brief Convert data value in the range of MinValues..MaxValues to the range 350 - 650
	*/
	Eigen::Matrix<double, 3, 1> waveLenghtToColor (const double wavelength);
	
	/**
	* @brief Convert data value in the range of MinValues..MaxValues
	*/
	Eigen::Matrix<double, 3, 1> waveLenghtToColorv2 (const double wavelength);

	/**
	* @brief 
	*/
	Eigen::Matrix<double, 3, 1> valueToColor (const double value, const double max, const double min);
	
	
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
	static double GravityModel(double latitude, double altitude)
	{
	    double g; /** g magnitude at zero altitude **/

	    /** Nominal Gravity model **/
	    g = GWGS0*((1+GWGS1*pow(sin(latitude),2))/sqrt(1-pow(ECC,2)*pow(sin(latitude),2)));

	    /** Gravity affects by the altitude (aprox the value r = Re **/
	    g = g*pow(Re/(Re+altitude), 2);

	    std::cout<<"Theoretical gravity for this location (WGS-84 ellipsoid model): "<< g<<" [m/s^2]\n";

	    return g;
	};
	
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
	static void SubstractEarthRotation(Eigen::Matrix <double, NUMAXIS, 1> *u, Eigen::Quaternion <double> *qb_g, double latitude)
	{
	    Eigen::Matrix <double, NUMAXIS, 1> v (EARTHW*cos(latitude), 0, EARTHW*sin(latitude)); /** vector of earth rotation components expressed in the geografic frame according to the latitude **/

	    /** Compute the v vector expressed in the body frame **/
	    v = (*qb_g) * v;
	    
	    #ifdef DEBUG_PRINTS
	    std::cout<<"Earth Rotation:"<<v<<"\n";
	    #endif

	    /** Subtract the earth rotation to the vector of inputs (u = u-v**/
	    (*u)  = (*u) - v;
	    
	    return;
	};
	
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
	static int CorrectMagneticDeclination (Eigen::Quaternion <double> *quat, double magnetic_declination,  int mode)
	{
	    Eigen::Matrix <double, NUMAXIS, 1> euler;
	
	    euler[2] = quat->toRotationMatrix().eulerAngles(2,1,0)[0];//YAW
	    euler[1] = quat->toRotationMatrix().eulerAngles(2,1,0)[1];//PITCH
	    euler[0] = quat->toRotationMatrix().eulerAngles(2,1,0)[2];//ROLL
	    
	    if (mode == EAST)
	    {
		#ifdef DEBUG_PRINTS
		std::cout << "[EAST] magnetic declination\n";
		#endif
		euler[2] -= magnetic_declination; /** Magnetic declination is positive **/
	    }
	    else if (mode == WEST)
	    {
		#ifdef DEBUG_PRINTS
		std::cout << "[WEST] magnetic declination\n";
		#endif
		euler[2] += magnetic_declination; /** Magnetic declination is negative **/
	    }
	    else
	    {
		#ifdef DEBUG_PRINTS
		std::cerr << "[ERROR] In the correction of the magnetic declination\n";
		#endif
		return ERROR_OUT;
	    }
	    
	    *quat = Eigen::Quaternion <double> (Eigen::AngleAxisd(euler[0], Eigen::Vector3d::UnitX())*
				Eigen::AngleAxisd(euler[1], Eigen::Vector3d::UnitY()) *
				Eigen::AngleAxisd(euler[2], Eigen::Vector3d::UnitZ()));
	    
	    return OK_LOCALIZATION;
	};
	
	/**
	* @brief Finite difference
	*
	* Its applies the maximum precission possible
	* depending of the vector size.
	* data[0] -> the newest value
	* data[size-1] -> the oldest value
	*
	* @author Javier Hidalgo Carrio.
	*
	* @param[in] Eigen::Matrix <double, Eigen::Dynamic, 1> &data vector of values
	* @param[in] delta_t sampling interval
	*
	* @return double with the result of the difference. NaN if error occured
	*
	*/
	static double finiteDifference (Eigen::Matrix <double, Eigen::Dynamic, 1> &data, double delta_t)
	{
	    double result;
	    
	    std::cout << "[Difference] data.rows(): "<<data.rows()<<"\n";
	    std::cout << "[Difference] data\n"<<data<<"\n";
	    
	    if ((data.rows() < 4)&&(data.rows() > 1))
		result = (data[0] - data[1])/delta_t;
	    else if (data.rows() < 7)
		result = - ((-11.0/6.0)*data[0] + 3.0*data[1] - (1.5)*data[2] + (1.0/3.0)*data[3])/(delta_t);
	    else if (data.rows() >= 7)
		result = (2.45*data[0] - 6.0*data[1] + 7.5*data[2] - (20.0/3.0)*data[3]+ 3.75*data[4] - 1.2*data[5] + (1.0/6.0)*data[6])/delta_t;
	    else
		result = std::numeric_limits<double>::quiet_NaN();
	    
	    std::cout << "[Difference] result: "<<result<<"\n";
	    if (result < 0.00)
		std::cout << "[Difference] Negative velocity!!\n";

	    return result;
	};
	  
    }; //end of measurement class

}//end of namespace localization

#endif