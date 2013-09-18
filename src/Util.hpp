/**\file measurement.hpp
 * Header function file and defines
 */

#ifndef _MEASUREMENT_ITEM_HPP_
#define _MEASUREMENT_ITEM_HPP_

#include <iostream> /** IO C++ Standard library */
#include <cmath> /** Math C++ Standard library */
#include <queue> /** FIFO queues **/
#include <deque> /**  Dynamic-sized FIFO queues to have pop_back and other methods **/
#include <algorithm> /** Algorithm C++ Standard library */
#include <Eigen/Geometry> /** Eigen data type for Matrix, Quaternion, etc... */
#include <Eigen/Core> /** Core methods of Eigen implementation **/
#include <Eigen/Dense> /** for the algebra and transformation matrices **/
#include <Eigen/Cholesky> /** For the Cholesky decomposition **/

/** Localization library  headers **/
#include "DataModel.hpp" /** For the quantities models **/
#include "Configuration.hpp" /** For the localization framework constant and configuration values **/
#include "DataTypes.hpp" /** Orogen compatible export data types **/

//#define UTIL_DEBUG_PRINTS 1

namespace localization	
{
    class Util
    {
	
    public:
	
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW
	
    private:

	/** **/
	/** MOTION MODEL VARIABLES (UPPER CASE MATRICES, LOWER CASE VECTORS) **/
	unsigned int ibuffer_size; /** Size of the buffered velocities and inertial values **/
	Eigen::Matrix <double, Eigen::Dynamic, 1> besselBCoeff, besselACoeff; /** IIR Bessel filter Coefficients **/
	Eigen::Matrix <double, Eigen::Dynamic, Eigen::Dynamic> arrayVelMM, arrayVelMMIIR; /** Array of past rover velocities from the motion model (for Bessel filter) **/
	Eigen::Quaternion <double> q_weight_distribution; /** Offset quaternion to set the uneven weight distribution of the rover **/
	std::deque < DataModel<double, 3> > velMM; /** Velocity model from navigation kinematics (already filtered by Bessel IIR) **/
	std::deque < DataModel<double, 3> > increVelMM; /** Increment in velocity model from navigation kinematics (already filtered by Bessel IIR) **/
	DataModel<double,3> dcontactAngle; /** delta contact angle from the navigation kinematics **/
	
	Eigen::Matrix <double,ENCODERS_VECTOR_SIZE,ENCODERS_VECTOR_SIZE> Rencoders; /** Measurement noise convariance matrix for joint velocities */
	Eigen::Matrix <double,NUMAXIS,NUMAXIS> Rangvelo; /** Measurement noise convariance matrix for angular velocities */
	
	/** **/
	/** PROPRIOCEPTIVE SENSORS VARIABLES **/
	Eigen::Matrix <double, ENCODERS_VECTOR_SIZE, 1>encodersvelocity; /** Sensed encoders velocities **/
	std::deque < Eigen::Vector3d, Eigen::aligned_allocator < std::pair < const int, Eigen::Vector3d > > > acc, angvelo; /** Vector of accelerations and angular velocities **/
	
	/** **/
	/** MOTION MODEL EXTRA VARIABLES (INFORMATION DEBUG) **/
	double leastSquaredNavigation; /** Least squared error of the navigation kinematics **/
	double leastSquaredSlip; /** Least squared error of the slip kinematics **/
	
	/** **/
	/** KK VARIABLES **/
	/** For the slip kinematics (each column is a wheel defined by a wheel_idx) **/
	Eigen::Matrix <double,NUMAXIS,NUMBER_OF_WHEELS> slipMatrix;
	
	/** Translation distances for the accelerometers w.r.t to rover body **/
	Eigen::Matrix <double,NUMAXIS,1> eccx, eccy, eccz; /** Accelerometers excentricity with respect to the body center of the robot **/
		
	/** Slip vector **/
	DataModel<double,3> slipModel, slipVector, slipError;
	

    public:
	
	/** Util contructor
         */
        Util();
	
	/** Util default descontructor
         */
        ~Util();
	
	
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

            #ifdef UTIL_DEBUG_PRINTS
	    std::cout<<"[UTIL_CLASS] Theoretical gravity for this location (WGS-84 ellipsoid model): "<< g<<" [m/s^2]\n";
            #endif

	    return g;
	};
	
	/**
	* @brief Substract the Earth rotation from the gyroscopes readout
	*
	* This function computes the substraction of the rotation of the Earth (EARTHW)
	* from the gyroscope values. This function uses quaternion of transformation from
	* the geographici to body frame and the latitude in radians.
	*
	* @author Javier Hidalgo Carrio.
	*
	* @param[in, out] *u pointer to angular velocity
	* @param[in] *q quaternion from geographic to body frame vb = qbg * vg (qbg is transf a vector  from g to b)
	* @param[in] latitude location latitude angle in radians
	*
	* @return void
	*
	*/
	static void SubstractEarthRotation(Eigen::Matrix <double, NUMAXIS, 1, Eigen::DontAlign> *u, Eigen::Quaternion <double, Eigen::DontAlign> *q, double latitude)
	{
	    Eigen::Matrix <double, NUMAXIS, 1> v (EARTHW*cos(latitude), 0, EARTHW*sin(latitude)); /** vector of earth rotation components expressed in the geografic frame according to the latitude **/

	    /** Compute the v vector expressed in the body frame **/
	    v = (*q) * v;

	    #ifdef UTIL_DEBUG_PRINTS
	    std::cout<<"[UTIL_CLASS] Earth Rotation:"<<v<<"\n";
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
	static bool CorrectMagneticDeclination (Eigen::Quaternion <double> *quat, double magnetic_declination,  int mode)
	{
	    Eigen::Matrix <double, NUMAXIS, 1> euler;
	
	    euler[2] = quat->toRotationMatrix().eulerAngles(2,1,0)[0];//YAW
	    euler[1] = quat->toRotationMatrix().eulerAngles(2,1,0)[1];//PITCH
	    euler[0] = quat->toRotationMatrix().eulerAngles(2,1,0)[2];//ROLL

	    if (mode == EAST)
	    {
		#ifdef UTIL_DEBUG_PRINTS
		std::cout << "[EAST] magnetic declination\n";
		#endif
		euler[2] -= magnetic_declination; /** Magnetic declination is positive **/
	    }
	    else if (mode == WEST)
	    {
		#ifdef UTIL_DEBUG_PRINTS
		std::cout << "[WEST] magnetic declination\n";
		#endif
		euler[2] += magnetic_declination; /** Magnetic declination is negative **/
	    }
	    else
	    {
		#ifdef UTIL_DEBUG_PRINTS
		std::cerr << "[ERROR] In the correction of the magnetic declination\n";
		#endif
		return false;
	    }

	    *quat = Eigen::Quaternion <double> (Eigen::AngleAxisd(euler[0], Eigen::Vector3d::UnitX())*
				Eigen::AngleAxisd(euler[1], Eigen::Vector3d::UnitY()) *
				Eigen::AngleAxisd(euler[2], Eigen::Vector3d::UnitZ()));

	    return true;
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

	    //std::cout << "[Difference] data.rows(): "<<data.rows()<<"\n";
	    //std::cout << "[Difference] data\n"<<data<<"\n";

	    if ((data.rows() < 4)&&(data.rows() > 1))
		result = (data[0] - data[1])/delta_t;
	    else if (data.rows() < 7)
		result = ((11.0/6.0)*data[0] - 3.0*data[1] + (1.5)*data[2] - (1.0/3.0)*data[3])/(delta_t);
	    else if (data.rows() >= 7)
		result = (2.45*data[0] - 6.0*data[1] + 7.5*data[2] - (20.0/3.0)*data[3]+ 3.75*data[4] - 1.2*data[5] + (1.0/6.0)*data[6])/delta_t;
	    else
		result = std::numeric_limits<double>::quiet_NaN();

	    //std::cout << "[Difference] result: "<<result<<"\n";
	    //if (result < 0.00)
	    //	std::cout << "[Difference] Negative velocity!!\n";

	    return result;
	};
	
	/**
	* @brief IIR filter
	*
	* Infinite Impulse Response (IIR) Filter for a Eigen Vector
	*
	* @author Javier Hidalgo Carrio.
	*
	* @param[in] norder integer number with the filter order
	* @param[in] b array of feedforward filter coefficients
	* @param[in] a array of feedback filter coefficients
	* @param[in] x array of inputs
	* @param[in] y array of outputs
	*
	* @return double with the result y[n]
	*
	*/
	static Eigen::Matrix <double, Eigen::Dynamic, 1> iirFilter (const unsigned int norder,
				    const Eigen::Matrix <double, Eigen::Dynamic, 1> &b, const Eigen::Matrix <double, Eigen::Dynamic, 1> &a,
				    const Eigen::Matrix <double, Eigen::Dynamic, Eigen::Dynamic> &x, Eigen::Matrix <double, Eigen::Dynamic, Eigen::Dynamic> &y)
	{
	    int i=norder;

	    std::cout<<"Filter n.cols:"<<x.cols()<<" n.cols:"<<y.cols()<<"\n";

	    if ((x.cols() >= norder+1) && (y.cols() >= norder+1))
	    {
		std::cout<<"Performing filter\n";
		/** Perform the filter **/
		y.col(i) = 1.0/a[0] * (b[0]*x.col(i) + b[1]*x.col(i-1) + b[2]*x.col(i-2)
				    + b[3]*x.col(i-3) + b[4]*x.col(i-4) + b[5]*x.col(i-5)
				    + b[6]*x.col(i-6) + b[7]*x.col(i-7) + b[8]*x.col(i-8)
				    - a[1]*y.col(i-1) - a[2]*y.col(i-2)
				    - a[3]*y.col(i-3) - a[4]*y.col(i-4)
				    - a[5]*y.col(i-5) - a[6]*y.col(i-6)
				    - a[7]*y.col(i-7) - a[8]*y.col(i-8));
	    }

	    std::cout<<"[IIR]Result:\n"<<y.col(i)<<"\n";

	    return y.col(i);
	};

        /**
        * @brief Computes the Bhattacharyya coefficient
        */
        template<typename _Scalar, int _DIM>
        static Eigen::Matrix<_Scalar, _DIM, _DIM> bhattacharyya (const DataModel<_Scalar, _DIM> &data1, const DataModel<_Scalar, _DIM> &data2)
        {
            _Scalar number_expo = std::numeric_limits<_Scalar>::quiet_NaN();
            Eigen::Matrix<_Scalar, _DIM, _DIM> Inv_half;
            Eigen::Matrix<_Scalar, _DIM, _DIM> BC;
            Eigen::Matrix<_Scalar, _DIM, _DIM> Cov1llt;
            Eigen::Matrix<_Scalar, _DIM, _DIM> Cov2llt;
            DataModel<_Scalar, _DIM> aux;

            /** Substraction **/
            aux = data1 - data2;
            Inv_half = (aux.Cov * 0.5).inverse();

            Eigen::Matrix<_Scalar, 1, 1> expo = -0.125 * (aux.data.transpose() * Inv_half * aux.data);
            number_expo = exp(expo[0]);

            Cov1llt = data1.Cov.llt().matrixL();
            Cov2llt = data2.Cov.llt().matrixL();
            Cov1llt = Cov1llt.llt().matrixL();
            Cov2llt = Cov2llt.llt().matrixL();

            BC = (Cov1llt*Cov2llt) * Inv_half.llt().matrixL();

            #ifdef UTIL_DEBUG_PRINTS
            std::cout << "[BHATTACHARYYA] number_expo is:\n" << number_expo << std::endl;
            std::cout << "[BHATTACHARYYA] Cov1llt is:\n" << Cov1llt << std::endl;
            std::cout << "[BHATTACHARYYA] Cov2llt is:\n" << Cov2llt << std::endl;
            std::cout << "[BHATTACHARYYA] Inv_half is:\n" << Inv_half << std::endl;
            std::cout << "[BHATTACHARYYA] BC(without expo) is:\n" << BC << std::endl;
            #endif

            BC = BC * number_expo;

            #ifdef UTIL_DEBUG_PRINTS
            std::cout << "[BHATTACHARYYA] BC is:\n" << BC << std::endl;
            #endif

            return BC;
        };

        template <typename _Scalar, int _DIM>
        static _Scalar mahalanobis (const DataModel<_Scalar, _DIM> &data1, const DataModel<_Scalar, _DIM> &data2)
        {
            _Scalar distance;
            DataModel<_Scalar, _DIM> aux = data1 - data2;

            distance = aux.data.transpose() * aux.Cov.inverse() * aux.data;

            #ifdef UTIL_DEBUG_PRINTS
            std::cout << "[MAHALANOBIS] Distance is:" << distance << std::endl;
            #endif

            return distance;
        }

        template <typename _MatrixType>
        static _MatrixType guaranteeSPD (const _MatrixType &A)
        {
            _MatrixType spdA;
            Eigen::VectorXd s;
            s.resize(A.rows(), 1);

            /**
             * Single Value Decomposition
            */
            Eigen::JacobiSVD <Eigen::MatrixXd > svdOfA (A, Eigen::ComputeThinU | Eigen::ComputeThinV);

            s = svdOfA.singularValues(); //!eigenvalues

            #ifdef UTIL_DEBUG_PRINTS
            std::cout<<"[SPD-SVD] s: \n"<<s<<"\n";
            std::cout<<"[SPD-SVD] svdOfA.matrixU():\n"<<svdOfA.matrixU()<<"\n";
            std::cout<<"[SPD-SVD] svdOfA.matrixV():\n"<<svdOfA.matrixV()<<"\n";

            Eigen::EigenSolver<_MatrixType> eig(A);
            std::cout << "[SPD-SVD] BEFORE: eigen values: " << eig.eigenvalues().transpose() << std::endl;
            #endif

            for (register int i=0; i<s.size(); ++i)
            {
                #ifdef UTIL_DEBUG_PRINTS
                std::cout<<"[SPD-SVD] i["<<i<<"]\n";
                #endif

                if (s(i) < 0.00)
                    s(i) = 0.00;
            }
            spdA = svdOfA.matrixU() * s.matrix().asDiagonal() * svdOfA.matrixV();

            #ifdef UTIL_DEBUG_PRINTS
            Eigen::EigenSolver<_MatrixType> eigSPD(spdA);
            if (eig.eigenvalues() == eigSPD.eigenvalues())
                std::cout<<"[SPD-SVD] EQUAL!!\n";

            std::cout << "[SPD-SVD] AFTER: eigen values: " << eigSPD.eigenvalues().transpose() << std::endl;
            #endif

            return spdA;
        };

        /**
         * @brief Check if NaN values
         */
        template<typename _Derived>
        static inline bool isnotnan(const Eigen::MatrixBase<_Derived>& x)
        {
            return ((x.array() == x.array())).all();
        };

        template<typename _Derived>
        static inline bool isfinite(const Eigen::MatrixBase<_Derived>& x)
        {
            return isnotnan(x - x);
        };

    }; //end of util class

}//end of namespace localization

#endif
