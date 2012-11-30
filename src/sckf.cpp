/**\file sckf.cpp
 *
 * This class has the primitive methods for an Stochastic Cloning Indirect Kalman Filter implementation
 * the stet vector are formed by the errors. Therefore the name of indirect kalman filter
 * The measurement 
 * 
 * 
 * @author Javier Hidalgo Carrio | DFKI RIC Bremen | javier.hidalgo_carrio@dfki.de
 * @date November 2012.
 * @version 1.0.
 */

#include <iostream> /**< IO C++ Standard library */
#include <algorithm> /**< Algorithm C++ Standard library */
#include <Eigen/LU> /**< Lineal algebra of Eigen */
#include <Eigen/SVD> /**< Singular Value Decomposition (SVD) of Eigen */
#include "sckf.hpp"

// #define DEBUG_PRINTS 1

using namespace localization;
using namespace Eigen;

void sckf::welcome()
{
	std::cout << "You successfully compiled and executed SCFK. Welcome!" << std::endl;
}


/**
* @brief Gets the current vector x
*/
Eigen::Matrix< double, Eigen::Dynamic, 1  > sckf::getStatex()
{
    return xki_k;
}

/**
* @brief Gets the current orientation of the robot in Quaternion
*/
Eigen::Quaternion< double > sckf::getAttitude()
{
    return this->q4;
}

/**
* @brief Gets the gravity value
*/
double sckf::getGravity()
{
    return this->gtilde.norm();
}


/**
* @brief Gets the current orientation of the robot in Euler angles
*/
Eigen::Matrix< double, NUMAXIS , 1  > sckf::getEuler()
{
    Eigen::Matrix <double, NUMAXIS, 1> euler;
    
    //std::cout << Eigen::Matrix3d(q4) << std::endl; 
    Eigen::Vector3d e = Eigen::Matrix3d(q4).eulerAngles(2,1,0);
    euler(0) = e[2]; 
    euler(1) = e[1]; 
    euler(2) = e[0]; 
    
    std::cout << "Attitude (getEuler): "<< euler(0)<<" "<<euler(1)<<" "<<euler(2)<<"\n";
    std::cout << "Attitude in degrees (getEuler): "<< euler(0)*R2D<<" "<<euler(1)*R2D<<" "<<euler(2)*R2D<<"\n";
    
    return euler;

}

/**
* @brief Gets Noise covariance matrix
*/
Eigen::Matrix< double, Eigen::Dynamic, Eigen::Dynamic > sckf::getCovariancex()
{
    return Pki_k;
}

/**
* @brief Gets Noise covariance matrix for attitude estimation
*/
Eigen::Matrix< double, sckf::A_STATE_VECTOR_SIZE, sckf::A_STATE_VECTOR_SIZE > sckf::getCovarianceAttitude()
{
    return Pki_k.block<sckf::A_STATE_VECTOR_SIZE, sckf::A_STATE_VECTOR_SIZE> ((sckf::E_STATE_VECTOR_SIZE*sckf::NUMBER_OF_WHEELS), (sckf::E_STATE_VECTOR_SIZE*sckf::NUMBER_OF_WHEELS));
}

/**
* @brief Return linear velocities
*/
Eigen::Matrix< double, NUMAXIS , 1  > sckf::getLinearVelocities()
{
    return this->linvelocity;
}

/**
* @brief Return angular velocities
*/
Eigen::Matrix< double, NUMAXIS , 1  > sckf::getAngularVelocities()
{
    return this->angvelocity;
}

/**
* @brief Return linear acceleration
*/
Eigen::Matrix< double, 3 , 1  > sckf::getLinearAcceleration()
{
    return this->accSimps.col(0);
}


/**
* @brief Return slip vector for wheel_idx
*/
Eigen::Matrix< double, NUMAXIS , 1  > sckf::getSlipVector(int wheel_idx)
{
    if (wheel_idx < NUMBER_OF_WHEELS)
	return this->slipMatrix.col(wheel_idx);
    
    return Eigen::Matrix <double, NUMAXIS, 1>::Zero();
}

/**
* @brief Return the contact angles
*/

Eigen::Matrix< double, Eigen::Dynamic, 1  > sckf::getContactAngles()
{
    return this->acontact;
}

/**
* @brief Return the K matrix
*/
Eigen::Matrix< double, Eigen::Dynamic, Eigen::Dynamic > sckf::getKalmanGain()
{
    return this->K;
}

/**
* @brief Return the K associated to the attitude
*/
Eigen::Matrix< double, Eigen::Dynamic, Eigen::Dynamic > sckf::getAttitudeKalmanGain()
{
    return this->K.block<sckf::A_STATE_VECTOR_SIZE, NUMAXIS> ((sckf::E_STATE_VECTOR_SIZE*sckf::NUMBER_OF_WHEELS), sckf::E_MEASUREMENT_VECTOR_SIZE);
}

/**
* @brief Return the innovation vector
*/
Eigen::Matrix< double, Eigen::Dynamic, 1  > localization::sckf::getInnovation()
{
    return this->innovation;
}


/**
* @brief This function Initialize Attitude
*/
int sckf::setAttitude(Eigen::Quaternion< double >& initq)
{

    
    /** Initial orientation **/
    q4 = initq;
	
    return OK;
    
}

/**
* @brief This function sets the gravity value
*/
void sckf::setGravity(double g)
{
    this->gtilde << 0, 0, g;
}


/**
* @brief This function set the initial Omega matrix
*/
int sckf::setOmega(Eigen::Matrix< double, NUMAXIS , 1  >& u)
{
    if (&u != NULL)
    {
	/** Initialization for quaternion integration **/
	oldomega4 << 0,-u(0), -u(1), -u(2),
	    u(0), 0, u(2), -u(1),
	    u(1), -u(2), 0, u(0),
	    u(2), u(1), -u(0), 0;
	
	return OK;
    }

    return ERROR;

}

/**
* @brief Gets the current state vector of the filter
*/
void sckf::setStatex(Eigen::Matrix< double, Eigen::Dynamic, 1  > &x_0)
{
    x_0.resize(sckf::X_STATE_VECTOR_SIZE,1);
    this->xki_k = x_0;
}

/**
* @brief This function set the Accelerometers excentricity
*/
void sckf::setEccentricity(Eigen::Matrix <double,NUMAXIS,1>  &eccx, Eigen::Matrix <double,NUMAXIS,1>  &eccy, Eigen::Matrix <double,NUMAXIS,1>  &eccz)
{
    this->eccx = eccx;
    this->eccy = eccy;
    this->eccz = eccz;
    
    #ifdef DEBUG_PRINTS
    std::cout<< "eccx is of size "<<eccx.rows()<<"x"<<eccx.cols()<<"\n";
    std::cout<< "eccx:\n"<<eccx<<"\n";
    std::cout<< "eccy is of size "<<eccy.rows()<<"x"<<eccy.cols()<<"\n";
    std::cout<< "eccy:\n"<<eccy<<"\n";
    std::cout<< "eccz is of size "<<eccz.rows()<<"x"<<eccz.cols()<<"\n";
    std::cout<< "eccz:\n"<<eccz<<"\n";
    #endif
}

/**
* @brief This function set the Accelerometers excentricity
*/
double sckf::simpsonsIntegral(double fa, double fm, double fb, double delta_ab)
{
    return delta_ab/6.0 * (fa + (4.0*fm) + fb);
}

/**
* @brief Set the Heading angle
*/
void sckf::setHeading(double yaw)
{
    Eigen::Matrix< double, NUMAXIS , 1> euler;
    Eigen::Vector3d e = Eigen::Matrix3d(q4).eulerAngles(2,1,0);
    
    euler(0) = e[2]; 
    euler(1) = e[1]; 
    euler(2) = yaw; 
    
    q4 = Eigen::Quaternion <double> (Eigen::AngleAxisd(euler[2], Eigen::Vector3d::UnitZ())*
	    Eigen::AngleAxisd(euler[1], Eigen::Vector3d::UnitY()) *
	    Eigen::AngleAxisd(euler[0], Eigen::Vector3d::UnitX()));
    
}


 /**
* @brief This function Initilize the vectors and matrices
*/
void sckf::Init(Eigen::Matrix< double, Eigen::Dynamic, Eigen::Dynamic >& P_0, Eigen::Matrix< double, Eigen::Dynamic, Eigen::Dynamic >& Qec,
		Eigen::Matrix< double, NUMAXIS , NUMAXIS  >& Qbg, Eigen::Matrix< double, NUMAXIS , NUMAXIS  >& Qba,
		Eigen::Matrix< double, NUMAXIS, NUMAXIS >& Rv, Eigen::Matrix< double, NUMAXIS , NUMAXIS  >& Rg,
		Eigen::Matrix <double,Eigen::Dynamic,Eigen::Dynamic> &Ren,
		Eigen::Matrix< double, NUMAXIS , NUMAXIS  >& Ra, Eigen::Matrix< double, NUMAXIS , NUMAXIS  >& Rm,
		double g, double alpha)
{
    
    /** Set the matrix and vector dimension to the static values of the class **/
    xk_k.resize(sckf::X_STATE_VECTOR_SIZE,1);
    xki_k.resize(sckf::X_STATE_VECTOR_SIZE,1);
    
    A.resize(sckf::A_STATE_VECTOR_SIZE,sckf::A_STATE_VECTOR_SIZE);
    Fki.resize(sckf::X_STATE_VECTOR_SIZE,sckf::X_STATE_VECTOR_SIZE);
    
    Qk.resize(sckf::X_STATE_VECTOR_SIZE,sckf::X_STATE_VECTOR_SIZE);
    
    Pk_k.resize(sckf::X_STATE_VECTOR_SIZE,sckf::X_STATE_VECTOR_SIZE);
    Pki_k.resize(sckf::X_STATE_VECTOR_SIZE,sckf::X_STATE_VECTOR_SIZE);
    
    H1a.resize(NUMAXIS, sckf::A_STATE_VECTOR_SIZE);
    H2a.resize(NUMAXIS, sckf::A_STATE_VECTOR_SIZE);
    Hk.resize(sckf::Z_MEASUREMENT_VECTOR_SIZE, sckf::X_STATE_VECTOR_SIZE);
    
    zki.resize(Z_MEASUREMENT_VECTOR_SIZE, 1);
    
    Rk.resize(sckf::Z_MEASUREMENT_VECTOR_SIZE, sckf::Z_MEASUREMENT_VECTOR_SIZE);
    
    innovation.resize(sckf::Z_MEASUREMENT_VECTOR_SIZE, 1);
    
    K.resize(sckf::X_STATE_VECTOR_SIZE, sckf::Z_MEASUREMENT_VECTOR_SIZE);
    
    /** Resizing dynamic arguments to the correct dimension to avoid matrices errors **/
    P_0.resize(sckf::X_STATE_VECTOR_SIZE, sckf::X_STATE_VECTOR_SIZE);
    Qec.resize((sckf::E_STATE_VECTOR_SIZE*sckf::NUMBER_OF_WHEELS), (sckf::E_STATE_VECTOR_SIZE*sckf::NUMBER_OF_WHEELS));
    Ren.resize(1+sckf::NUMBER_OF_WHEELS, 1+sckf::NUMBER_OF_WHEELS);
    
    /** Gravitation acceleration **/
    gtilde << 0, 0, g;

    /** Dip angle (alpha is in rad) **/
    mtilde(0) = cos(alpha);
    mtilde(1) = 0;
    mtilde(2) = -sin(alpha);

    
    /** Kalman filter state, error covariance and process noise covariance **/
    xki_k = Matrix <double,sckf::X_STATE_VECTOR_SIZE,1>::Zero();
    xk_k = xki_k;
    
    Qk = Matrix <double,sckf::X_STATE_VECTOR_SIZE,sckf::X_STATE_VECTOR_SIZE>::Zero();
    Qk.block <(sckf::E_STATE_VECTOR_SIZE*sckf::NUMBER_OF_WHEELS), (sckf::E_STATE_VECTOR_SIZE*sckf::NUMBER_OF_WHEELS)> (0,0) = Qec;
    Qk.block <NUMAXIS, NUMAXIS> ((sckf::E_STATE_VECTOR_SIZE*sckf::NUMBER_OF_WHEELS),(sckf::E_STATE_VECTOR_SIZE*sckf::NUMBER_OF_WHEELS)) = 0.25 * Rg;
    Qk.block <NUMAXIS, NUMAXIS> ((sckf::E_STATE_VECTOR_SIZE*sckf::NUMBER_OF_WHEELS)+NUMAXIS,(sckf::E_STATE_VECTOR_SIZE*sckf::NUMBER_OF_WHEELS)+NUMAXIS) = Qbg;
    Qk.block <NUMAXIS, NUMAXIS> ((sckf::E_STATE_VECTOR_SIZE*sckf::NUMBER_OF_WHEELS)+(2*NUMAXIS),(sckf::E_STATE_VECTOR_SIZE*sckf::NUMBER_OF_WHEELS)+(2*NUMAXIS)) = Qba;
    

    /** Assign the initial values **/
    Pki_k = P_0;
    
    /** Asign the initial value for the measurement matrix **/
    Hk = Matrix <double,sckf::Z_MEASUREMENT_VECTOR_SIZE, sckf::X_STATE_VECTOR_SIZE>::Zero();
    H1a = Matrix <double,NUMAXIS,sckf::A_STATE_VECTOR_SIZE>::Zero();
    H2a = Matrix <double,NUMAXIS,sckf::A_STATE_VECTOR_SIZE>::Zero();
    H1a(0,6) = 1; H1a(1,7) = 1; H1a(2,8) = 1;
    
    
    /** System matrix A **/
    A = Matrix <double,sckf::A_STATE_VECTOR_SIZE, sckf::A_STATE_VECTOR_SIZE>::Zero();      
    A(0,3) = -0.5;A(1,4) = -0.5;A(2,5) = -0.5;
    
    /** Initial measurement noise **/
    Rk = Matrix <double,sckf::Z_MEASUREMENT_VECTOR_SIZE,sckf::Z_MEASUREMENT_VECTOR_SIZE>::Zero();
    RHist = Eigen::Matrix <double,NUMAXIS,NUMAXIS*M1>::Zero();
    
    /** Fill matrix Rv, Rg, Re, Ra and Rm **/
    this->Rv = Rv;
    this->Rg = Rg;
    this->Ren = Ren;
    this->Ra = Ra;
    this->Rm = Rm;
    
    /** Initial bias **/
    bghat = Matrix <double,NUMAXIS,1>::Zero();
    bahat = Matrix <double,NUMAXIS,1>::Zero();
    
    
    /** Default omega matrix **/
    oldomega4 << 0 , 0 , 0 , 0,
	0 , 0 , 0 , 0,
	0 , 0 , 0 , 0,
	0 , 0 , 0 , 0;
    
    /** Initial quaternion in Init is NaN**/
    q4.w() = std::numeric_limits<double>::quiet_NaN();
    q4.x() = std::numeric_limits<double>::quiet_NaN();
    q4.y() = std::numeric_limits<double>::quiet_NaN();
    q4.z() = std::numeric_limits<double>::quiet_NaN();
    
    
    /** Default initial bias **/
    bghat << 0.00, 0.00, 0.00;
    bahat << 0.00, 0.00, 0.00;
    
    /** Initial contact angle **/
    acontact << 0.00, 0.00, 0.00, 0.00;
    
    /** Initial slip matrix **/
    slipMatrix = Eigen::Matrix <double, NUMAXIS, NUMBER_OF_WHEELS>::Zero();
    
    /** Velocities form the filter (acce and bias substracted correctly) **/
    linvelocity = Eigen::Matrix <double, NUMAXIS, 1>::Zero();
    angvelocity = Eigen::Matrix <double, NUMAXIS, 1>::Zero();
    
    /** Variable in the adaptive algorithm **/
    r1count = 0;
    r2count = R2COUNT;
    
    /** Circular vector for the integration **/
    cbAccX.set_capacity(INTEGRATION_XAXIS_WINDOW_SIZE);
    cbAccY.set_capacity(INTEGRATION_YAXIS_WINDOW_SIZE);
    cbAccZ.set_capacity(INTEGRATION_ZAXIS_WINDOW_SIZE);
    
    cbAngveloX.set_capacity(ANGVELO_WINDOW_SIZE);
    cbAngveloY.set_capacity(ANGVELO_WINDOW_SIZE);
    cbAngveloZ.set_capacity(ANGVELO_WINDOW_SIZE);
    
    /** Print filter information **/
//     #ifdef DEBUG_PRINTS
    std::cout<< "xki_k is of size "<<xki_k.rows()<<"x"<<xki_k.cols()<<"\n";
    std::cout<< "xki_k:\n"<<xki_k<<"\n";
    std::cout<< "xk_k is of size "<<xk_k.rows()<<"x"<<xk_k.cols()<<"\n";
    std::cout<< "xk_k:\n"<<xk_k<<"\n";
    std::cout<< "A is of size "<<A.rows()<<"x"<<A.cols()<<"\n";
    std::cout<< "A:\n"<<A<<"\n";
    std::cout<< "Fki is of size "<<Fki.rows()<<"x"<<Fki.cols()<<"\n";
    std::cout<< "Fki:\n"<<Fki<<"\n";
    std::cout<< "Pk+i|k is of size "<<Pki_k.rows()<<"x"<<Pki_k.cols()<<"\n";
    std::cout<< "Pk+i|k:\n"<<Pki_k<<"\n";
    std::cout<< "Pk|k is of size "<<Pk_k.rows()<<"x"<<Pk_k.cols()<<"\n";
    std::cout<< "Pk|k:\n"<<Pk_k<<"\n";
    std::cout<< "Qk|k is of size "<<Qk.rows()<<"x"<<Qk.cols()<<"\n";
    std::cout<< "Qk|k:\n"<<Qk<<"\n";
    std::cout<< "Rk is of size "<<Rk.rows()<<"x"<<Rk.cols()<<"\n";
    std::cout<< "Rk:\n"<<Rk<<"\n";
    std::cout<< "H1a is of size "<<H1a.rows()<<"x"<<H1a.cols()<<"\n";
    std::cout<< "H1a:\n"<<H1a<<"\n";
    std::cout<< "H2a is of size "<<H2a.rows()<<"x"<<H2a.cols()<<"\n";
    std::cout<< "H2a:\n"<<H2a<<"\n";
    std::cout<< "Hk is of size "<<Hk.rows()<<"x"<<Hk.cols()<<"\n";
    std::cout<< "Hk:\n"<<Hk<<"\n";
    std::cout<< "zki is of size "<<zki.rows()<<"x"<<zki.cols()<<"\n";
    std::cout<< "zki:\n"<<zki<<"\n";
    std::cout<< "RHist is of size "<<RHist.rows()<<"x"<<RHist.cols()<<"\n";
    std::cout<< "RHist:\n"<<RHist<<"\n";
    std::cout<< "Rv is of size "<<Rv.rows()<<"x"<<Rv.cols()<<"\n";
    std::cout<< "Rv:\n"<<Rv<<"\n";
    std::cout<< "Rg is of size "<<Rg.rows()<<"x"<<Rg.cols()<<"\n";
    std::cout<< "Rg:\n"<<Rg<<"\n";
    std::cout<< "Ra is of size "<<Ra.rows()<<"x"<<Ra.cols()<<"\n";
    std::cout<< "Ra:\n"<<Ra<<"\n";
    std::cout<< "Ren is of size "<<Ren.rows()<<"x"<<Ren.cols()<<"\n";
    std::cout<< "Ren:\n"<<Ren<<"\n";
    std::cout<< "Rm is of size "<<Rm.rows()<<"x"<<Rm.cols()<<"\n";
    std::cout<< "Rm:\n"<<Rm<<"\n";
    std::cout<< "Rk is of size "<<Rk.rows()<<"x"<<Rk.cols()<<"\n";
    std::cout<< "Rk:\n"<<Rk<<"\n";
    std::cout<< "K is of size "<<K.rows()<<"x"<<K.cols()<<"\n";
    std::cout<< "K:\n"<<K<<"\n";
    std::cout<<"\n";
    std::cout<< "mtilde is of size "<<mtilde.rows()<<"x"<<mtilde.cols()<<"\n";
    std::cout<< "mtilde:\n"<<mtilde<<"\n";
    std::cout<< "gtilde is of size "<<gtilde.rows()<<"x"<<gtilde.cols()<<"\n";
    std::cout<< "gtilde:\n"<<gtilde<<"\n";
    std::cout<< "bghat is of size "<<bghat.rows()<<"x"<<bghat.cols()<<"\n";
    std::cout<< "bghat:\n"<<bghat<<"\n";
    std::cout<< "bahat is of size "<<bahat.rows()<<"x"<<bahat.cols()<<"\n";
    std::cout<< "bahat:\n"<<bahat<<"\n";
//     #endif

}

/**
* @brief Performs the prediction step of the filter.
*/
void sckf::predict(Eigen::Matrix< double, NUMAXIS , 1  >& u, double dt)
{
    
    Eigen::Matrix <double,NUMAXIS,NUMAXIS> vec2product; /** Vec 2 product  matrix */
    Eigen::Matrix <double,NUMAXIS,1> angvelo; /** Angular velocity */
    Eigen::Matrix <double,QUATERSIZE,QUATERSIZE> omega4; /** Quaternion integration matrix */
    Eigen::Matrix <double,QUATERSIZE,1> quat; /** Quaternion integration matrix */
    Eigen::Matrix <double,NUMAXIS, NUMAXIS> Fs; /** System matrix of a slip vector */
    Eigen::Matrix <double,NUMBER_OF_WHEELS, NUMBER_OF_WHEELS> Fc; /** System matrix of contact angles */
    Eigen::Matrix <double,sckf::X_STATE_VECTOR_SIZE,sckf::X_STATE_VECTOR_SIZE> dFki; /** Discrete System matrix */
    Eigen::Matrix <double,X_STATE_VECTOR_SIZE,X_STATE_VECTOR_SIZE> Qdk; /** Discrete Qk matrix */

    /** Compute the vector2product matrix with the angular velocity **/
    angvelo = u - bghat; /** Eliminate the Bias **/
    angvelocity = angvelo;
    
    #ifdef DEBUG_PRINTS
    std::cout<<"[Predict] angevelo:\n"<<angvelo<<"\n";
    #endif

    vec2product << 0, -angvelo(2), angvelo(1),
		angvelo(2), 0, -angvelo(0),
		-angvelo(1), angvelo(0), 0;
		
    /** Compute the dA Matrix of the attitude part **/
    A.block<NUMAXIS, NUMAXIS> (0,0) = -vec2product;
    
    /** Form a single wheel position error and contact angle **/
    Fs = Eigen::Matrix <double,NUMAXIS, NUMAXIS>::Zero();
    Fc = Eigen::Matrix <double,NUMBER_OF_WHEELS, NUMBER_OF_WHEELS>::Zero();
   
    /** Form the complete system model matrix **/
    for (int i = 0; i < NUMBER_OF_WHEELS; i++)
	Fki.block<NUMAXIS, NUMAXIS> (i*NUMAXIS,i*NUMAXIS) = Fs;
    
    Fki.block<sckf::NUMBER_OF_WHEELS,sckf::NUMBER_OF_WHEELS> (sckf::NUMBER_OF_WHEELS*NUMAXIS, sckf::NUMBER_OF_WHEELS*NUMAXIS) = Fc;

    Fki.block<sckf::A_STATE_VECTOR_SIZE, sckf::A_STATE_VECTOR_SIZE>(sckf::NUMBER_OF_WHEELS*sckf::E_STATE_VECTOR_SIZE, sckf::NUMBER_OF_WHEELS*sckf::E_STATE_VECTOR_SIZE) = A;
    
    dFki = Eigen::Matrix<double,sckf::X_STATE_VECTOR_SIZE,sckf::X_STATE_VECTOR_SIZE>::Identity() + Fki * dt + Fki * Fki * pow(dt,2)/2;
    
    #ifdef DEBUG_PRINTS
    std::cout<< "[Predict] xki|k is of size "<<xki_k.rows()<<"x"<<xki_k.cols()<<"\n";
    std::cout<< "[Predict] xki_k:\n"<<xki_k<<"\n";
    #endif
    
    /** Propagate the vector through the system **/
    xki_k = dFki * xki_k;
    
    #ifdef DEBUG_PRINTS
    std::cout<< "[After Predict] xki_k:\n"<<xki_k<<"\n";
    std::cout<< "[Predict] Fki is of size "<<Fki.rows()<<"x"<<Fki.cols()<<"\n";
    std::cout<< "[Predict] Fki:\n"<<Fki<<"\n";
    std::cout<< "[Predict] dFki is of size "<<dFki.rows()<<"x"<<dFki.cols()<<"\n";
    std::cout<< "[Predict] dFki:\n"<<dFki<<"\n";
    #endif
    
    /** Form the system noise matrix for the attitude **/
    Qdk = Qk*dt + 0.5*dt*dt*Fki*Qk + 0.5*dt*dt*Qk*Fki.transpose();
    Qdk = 0.5*(Qdk + Qdk.transpose());
    
    #ifdef DEBUG_PRINTS
    std::cout<< "[Predict] Qdk is of size "<<Qdk.rows()<<"x"<<Qdk.cols()<<"\n";
    std::cout<< "[Predict] Qdk:\n"<<Qdk<<"\n";
    #endif
    
    Pki_k = dFki*Pki_k*dFki.transpose() + Qdk;
    
    #ifdef DEBUG_PRINTS
    std::cout<< "[Predict] Pki_k:\n"<<Pki_k<<"\n";
    #endif
        
    omega4 << 0,-angvelo(0), -angvelo(1), -angvelo(2),
	    angvelo(0), 0, angvelo(2), -angvelo(1),
	    angvelo(1), -angvelo(2), 0, angvelo(0),
	    angvelo(2), angvelo(1), -angvelo(0), 0;
	    
    quat(0) = q4.w();
    quat(1) = q4.x();
    quat(2) = q4.y();
    quat(3) = q4.z();

    /** Quaternion integration **/
    quat = (Matrix<double,QUATERSIZE,QUATERSIZE>::Identity() +(0.75 * omega4 *dt)- (0.25 * oldomega4 * dt) -
    ((1/6) * angvelo.squaredNorm() * pow(dt,2) *  Matrix<double,QUATERSIZE,QUATERSIZE>::Identity()) -
    ((1/24) * omega4 * oldomega4 * pow(dt,2)) - ((1/48) * angvelo.squaredNorm() * omega4 * pow(dt,3))) * quat;

    /** Store in a quaternion form **/
    q4.w() = quat(0);
    q4.x() = quat(1);
    q4.y() = quat(2);
    q4.z() = quat(3);
    q4.normalize();

    oldomega4 = omega4;

    return;

}

/**
* @brief Conversion Quaternion to DCM (Direct Cosine Matrix) (Alternative to Eigen)
*/
void Quaternion2DCM(Eigen::Quaternion< double >* q, Eigen::Matrix< double, NUMAXIS, NUMAXIS  >*C)
{
    double q0, q1, q2, q3;

    if (C != NULL)
    {
    /** Take the parameters of the quaternion */
    q0 = q->w();
    q1 = q->x();
    q2 = q->y();
    q3 = q->z();
    
    /** Create the DCM matrix from the actual quaternion */
    (*C)(0,0) = 2 * q0 * q0 + 2 * q1 * q1 - 1;
    (*C)(0,1) = 2 * q1 * q2 + 2 * q0 * q3;
    (*C)(0,2) = 2 * q1 * q3 - 2 * q0 * q2;
    (*C)(1,0) = 2 * q1 * q2 - 2 * q0 * q3;
    (*C)(1,1) = 2 * q0 * q0 + 2 * q2 * q2 - 1;
    (*C)(1,2) = 2 * q2 * q3 + 2 * q0 * q1;
    (*C)(2,0) = 2 * q1 * q3 + 2 * q0 * q2;
    (*C)(2,1) = 2 * q2 * q3 - 2 * q0 * q1;
    (*C)(2,2) = 2 * q0 * q0 + 2 * q3 * q3 - 1;	
    }
    
    return;
}


Eigen::Matrix< double, 3 , 1  > sckf::accIntegrationWindow(double dt)
{
    Eigen::Matrix <double, NUMAXIS, 1> localvelocity;
    Eigen::Matrix <double,NUMAXIS,NUMAXIS> gyros2product; /** Vec 2 product  matrix for the gyroscopes (angular velocity) */
    
    localvelocity.setZero();
        
    for (unsigned int i=0; i<cbAccX.size(); i++)
    {
	gyros2product << 0, -cbAngveloZ[i], cbAngveloY[i],
		cbAngveloZ[i], 0, -cbAngveloX[i],
		-cbAngveloY[i], cbAngveloX[i], 0;
		
	localvelocity[0] += (cbAccX[i] * dt) - (gyros2product.row(0) * eccx);
    }
    
    for (unsigned int i=0; i<cbAccY.size(); i++)
    {
	gyros2product << 0, -cbAngveloZ[i], cbAngveloY[i],
		cbAngveloZ[i], 0, -cbAngveloX[i],
		-cbAngveloY[i], cbAngveloX[i], 0;
		
	localvelocity[1] += (cbAccY[i] * dt) - (gyros2product.row(1) * eccy);
    }
    
    for (unsigned int i=0; i<cbAccZ.size(); i++)
    {
	gyros2product << 0, -cbAngveloZ[i], cbAngveloY[i],
		cbAngveloZ[i], 0, -cbAngveloX[i],
		-cbAngveloY[i], cbAngveloX[i], 0;
		
	localvelocity[2] += (cbAccZ[i] * dt) - (gyros2product.row(2) * eccz);
    }
    
    return localvelocity;
}



void sckf::update(Eigen::Matrix <double,Eigen::Dynamic,Eigen::Dynamic> &He, Eigen::Matrix <double,Eigen::Dynamic,Eigen::Dynamic> &Be,
		  Eigen::Matrix< double, Eigen::Dynamic, 1  >& encoders, Eigen::Matrix <double,NUMAXIS, 1> &vel_model, Eigen::Matrix< double, NUMAXIS , 1  >& acc,
		  Eigen::Matrix< double, NUMAXIS , 1  >& gyro, Eigen::Matrix< double, NUMAXIS , 1  >& mag, double dt, bool magn_on_off)
{
    
    register int j;
    Eigen::Matrix <double,NUMAXIS,NUMAXIS> Cq; /** Rotational matrix */
    Eigen::Matrix <double,NUMAXIS,NUMAXIS> gtilde2product; /** Vec 2 product  matrix for the gravioty vector in body frame*/
    Eigen::Matrix <double,NUMAXIS,NUMAXIS> gyros2product; /** Vec 2 product  matrix for the gyroscopes (angular velocity) */
    Eigen::Matrix <double,NUMAXIS,1> angvelo; /** Angular velocity */
    Eigen::Matrix <double,NUMAXIS,NUMAXIS> fooR2; /**  Measurement noise matrix from accelerometers matrix Ra*/
    Eigen::Matrix <double,A_STATE_VECTOR_SIZE,1> xa_k; /** Attitude part of the state vector xk+i|k */
    Eigen::Matrix <double,A_STATE_VECTOR_SIZE,A_STATE_VECTOR_SIZE> P1a; /** Error convariance matrix for measurement 1 of the attitude */
    Eigen::Matrix <double,A_STATE_VECTOR_SIZE,A_STATE_VECTOR_SIZE> P2a; /** Error convariance matrix for measurement 2 of the attitude */
    Eigen::Matrix <double,A_STATE_VECTOR_SIZE,A_STATE_VECTOR_SIZE> auxM; /** Auxiliar matrix for computing Kalman gain in measurement */
    Eigen::Matrix <double,A_STATE_VECTOR_SIZE, NUMAXIS> K1a; /** Kalman Gain matrix for measurement 1 */
    Eigen::Matrix <double,A_STATE_VECTOR_SIZE, NUMAXIS> K2a; /** Kalman Gain matrix for measurement 2 */
    Eigen::Matrix <double,NUMAXIS,NUMAXIS> R1a; /** Measurement noise convariance matrix for measurement 1 */
    Eigen::Matrix <double,NUMAXIS,NUMAXIS> Uk; /** Uk measurement noise convariance matrix for the adaptive algorithm */
    Eigen::Matrix <double,NUMAXIS,NUMAXIS> Qstar; /** External acceleration covariance matrix */
    Eigen::Quaternion <double> qe;  /** Attitude error quaternion */
    Eigen::Matrix <double,NUMAXIS,1> gtilde_body; /** Gravitation in the body frame */
    Eigen::Matrix <double,NUMAXIS,1> mtilde_body; /** Magnetic field in the body frame */
    Eigen::Matrix <double,NUMAXIS,NUMAXIS> u; /** Unitary matrix U for the SVD decomposition */
    Eigen::Matrix <double,NUMAXIS,1> s; /** Unitary matrix V for the SVD decomposition */
    Eigen::Matrix <double,NUMAXIS,1> lambda; /** Lambda vector for the adaptive algorithm */
    Eigen::Matrix <double,NUMAXIS,1> mu; /** mu vector for the adaptive algorithm */
    Eigen::Matrix <double,NUMAXIS,1> z1a; /** Measurement vector 1 Acc */
    Eigen::Matrix <double,NUMAXIS,1> z2a; /** Measurement vector 2 Mag */
    Eigen::Matrix <double,Y_MEASUREMENT_VECTOR_SIZE,1> ye; /** Measurement vectors for the rover position error ze = Be*ye */
    Eigen::Matrix <double,NUMAXIS,1> auxvector; /** Auxiliar vector variable */
    
    
    /** Print filter information **/
    #ifdef DEBUG_PRINTS
    std::cout<<"[Update] Be is of size "<<Be.rows()<<"x"<<Be.cols()<<"\n";
    std::cout<<"[Update] Be:\n"<<Be<<"\n";
    std::cout<<"[Update] He is of size "<<He.rows()<<"x"<<He.cols()<<"\n";
    std::cout<<"[Update] He:\n"<<He<<"\n";
    #endif
    
    /** First measurement step (Pitch and Roll correction from Acc) **/
    
    /** Copy the attitude part of the state vector and covariance matrix **/
    xa_k = xki_k.block<sckf::A_STATE_VECTOR_SIZE, 1> ((sckf::E_STATE_VECTOR_SIZE*sckf::NUMBER_OF_WHEELS), 0);
    P1a = Pki_k.block<sckf::A_STATE_VECTOR_SIZE, sckf::A_STATE_VECTOR_SIZE> ((sckf::E_STATE_VECTOR_SIZE*sckf::NUMBER_OF_WHEELS), (sckf::E_STATE_VECTOR_SIZE*sckf::NUMBER_OF_WHEELS));
    
//     this->getEuler();
//     
//     /** Calculate the gravity vector in the body frame **/
//     gtilde_body = q4.toRotationMatrix()*gtilde;
//     gtilde2product << 0, -gtilde_body(2), gtilde_body(1),
// 		gtilde_body(2), 0, -gtilde_body(0),
// 		-gtilde_body(1), gtilde_body(0), 0;
// 		
//     #ifdef DEBUG_PRINTS
//     std::cout<<"[Update] gtilde_body of size "<<gtilde_body.rows()<<"x"<<gtilde_body.cols()<<"\n";
//     std::cout<<"[Update] g in body_frame:\n"<<gtilde_body<<"\n";
//     #endif
	
    /** Create the orientation matrix from the quaternion **/
    Quaternion2DCM (&q4, &Cq);
    
    /** Calculate the gravity vector in the body frame **/
    gtilde_body = Cq * gtilde;
    
    gtilde2product << 0, -gtilde_body(2), gtilde_body(1),
		    gtilde_body(2), 0, -gtilde_body(0),
		    -gtilde_body(1), gtilde_body(0), 0;

    #ifdef DEBUG_PRINTS
    std::cout<<"[Update] gtilde_body of size "<<gtilde_body.rows()<<"x"<<gtilde_body.cols()<<"\n";
    std::cout<<"[Update] g in body_frame:\n"<<gtilde_body<<"\n";
    #endif
      
    /** Eliminate the Bias from gyros**/
    angvelo = gyro - bghat; 
    angvelocity = angvelo;
    
    /** Push it in the blobal circular vector buffer **/
    cbAngveloX.push_back(angvelo[0]);
    cbAngveloY.push_back(angvelo[1]);
    cbAngveloZ.push_back(angvelo[2]);

    /** Compute the vector2product matrix with the angular velocity **/
    /** in order to remove the centripetal velocities **/
    gyros2product << 0, -angvelo(2), angvelo(1),
		angvelo(2), 0, -angvelo(0),
		-angvelo(1), angvelo(0), 0;
		
    /** Form the observation matrix Hk **/
    Hk.block<sckf::E_MEASUREMENT_VECTOR_SIZE, (sckf::E_STATE_VECTOR_SIZE*sckf::NUMBER_OF_WHEELS)> (0,0) = He;
    
    /** Form the matrix for the measurement 1 of the attitude (acc correction) **/
    H1a.block<NUMAXIS, NUMAXIS> (0,0) = 2*gtilde2product;
    
    /** Copy to the whole observation matrix **/
    Hk.block<NUMAXIS, sckf::A_STATE_VECTOR_SIZE> (sckf::E_MEASUREMENT_VECTOR_SIZE,(sckf::E_STATE_VECTOR_SIZE*sckf::NUMBER_OF_WHEELS)) = H1a;
    
    #ifdef DEBUG_PRINTS
    std::cout<<"[Update] Hk is of size "<<Hk.rows()<<"x"<<Hk.cols()<<"\n";
    std::cout<<"[Update] Hk:\n"<<Hk<<"\n";
    #endif
    
    /** Form the measurement vector z1a for the attitude (the real acceleration value) **/
    z1a = acc - bahat - gtilde_body;
    
    /** Push the acceleration value in the circular vector array **/
    cbAccX.push_back(z1a(0));
    cbAccY.push_back(z1a(1));
    cbAccZ.push_back(z1a(2));
    
    /** Store the z1a value in the global variable for the Simpson's rule integration **/
    accSimps.col(2) = accSimps.col(1); //t-1 -> t-2
    accSimps.col(1) = accSimps.col(0); //t -> t-1
    accSimps.col(0) = z1a; //actual sample -> t
    
    /** Compute the velocities and store them in the global variables **/
    linvelocity = accIntegrationWindow(dt);
    
//     linvelocity[0] = (z1a[0] * dt) - (gyros2product.row(0) * eccx);//simpsonsIntegral(accSimps(0,2), accSimps(0,1), accSimps(0,0), 2*dt) - (gyros2product.row(0) * eccx);
//     linvelocity[1] = (z1a[1] * dt) - (gyros2product.row(1) * eccy);//simpsonsIntegral(accSimps(1,2), accSimps(1,1), accSimps(1,0), 2*dt) - (gyros2product.row(1) * eccy);
//     linvelocity[2] = (z1a[2] * dt) - (gyros2product.row(2) * eccz);//simpsonsIntegral(accSimps(2,2), accSimps(2,1), accSimps(2,0), 2*dt) - (gyros2product.row(2) * eccz);
    
    
    /** Form the measurement vector ye of the rover position error (Be*ye) **/
    ye.block<NUMAXIS, 1> (0, 0) = vel_model + linvelocity;
    ye.block<NUMAXIS, 1> (NUMAXIS, 0) = angvelo;
    ye.block<1 + sckf::NUMBER_OF_WHEELS, 1> ((2*NUMAXIS), 0) = encoders;
    
    
    #ifdef DEBUG_PRINTS
    std::cout<<"[Update] linvelocity(simpson's rule:\n"<<linvelocity<<"\n";
    std::cout<<"[Update] linvelocity(before):\n"<<ye.block<NUMAXIS, 1> (0, 0)<<"\n";
    std::cout<<"[Update] acc is of size "<<acc.rows()<<"x"<<acc.cols()<<"\n";
    std::cout<<"[Update] acc:\n"<<acc<<"\n";
    std::cout<<"[Update] bahat:\n"<<bahat<<"\n";
    std::cout<<"[Update] dt:\n"<<dt<<"\n";
    std::cout<<"[Update] gyros2product.row(0) * eccx:\n"<<gyros2product.row(0) * eccx<<"\n";
    std::cout<<"[Update] ye is of size "<<ye.rows()<<"x"<<ye.cols()<<"\n";
    std::cout<<"[Update] ye:\n"<<ye<<"\n";
    std::cout<<"[Update] ze:\n"<<Be*ye<<"\n";
    #endif
    
    /** Form the complete zk vector **/
    zki.block<sckf::E_MEASUREMENT_VECTOR_SIZE, 1> (0,0)= Be*ye;
    zki.block<NUMAXIS, 1> (sckf::E_MEASUREMENT_VECTOR_SIZE, 0)= z1a;
    
    #ifdef DEBUG_PRINTS
    std::cout<<"[Update] zki is of size "<<zki.rows()<<"x"<<zki.cols()<<"\n";
    std::cout<<"[Update] zki:\n"<<zki<<"\n";
    #endif
   
    
    /** The adaptive algorithm for the attitude, the Uk matrix and SVD part **/
    R1a = (z1a - H1a*xa_k) * (z1a - H1a*xa_k).transpose();
    RHist.block <NUMAXIS, NUMAXIS> (0, (r1count*(M1-1))%M1) = R1a;
    
    /** r1count + 1 modulus the number of history M1 **/
    r1count++; 

    /** Starting the Uk is R **/
    Uk = R1a;
    for (j=0; j<M1; j++)
    {
	Uk = Uk + RHist.block <NUMAXIS, NUMAXIS> (0,NUMAXIS*j);
    }
    
    Uk = Uk/static_cast<double>(M1);
    
    fooR2 = H1a*P1a*H1a.transpose() + Ra;
    
    /**
    * Single Value Decomposition
    */
    JacobiSVD <MatrixXd > svdOfR1a(Uk, ComputeThinU);

    s = svdOfR1a.singularValues();
    u = svdOfR1a.matrixU();
    
    lambda << s(0), s(1), s(2);
    
    mu(0) = (u.transpose().row(0) * fooR2).dot(u.col(0));
    mu(1) = (u.transpose().row(1) * fooR2).dot(u.col(1));
    mu(2) = (u.transpose().row(2) * fooR2).dot(u.col(2));
    
    if ((lambda - mu).maxCoeff() > GAMMA)
    {
	#ifdef DEBUG_PRINTS
	std::cout<<"[Update] "<<(lambda - mu).maxCoeff()<<" Bigger than GAMMA\n";
	#endif
	
	r2count = 0;
	auxvector(0) = std::max(lambda(0)-mu(0),(double)0.00);
	auxvector(1) = std::max(lambda(1)-mu(1),(double)0.00);
	auxvector(2) = std::max(lambda(2)-mu(2),(double)0.00);
	
	Qstar = auxvector(0) * u.col(0) * u.col(0).transpose() + auxvector(1) * u.col(1) * u.col(1).transpose() + auxvector(2) * u.col(2) * u.col(2).transpose();
    }
    else
    {
	#ifdef DEBUG_PRINTS
	std::cout<<"[Update] r2count: "<<r2count<<"\n";
	#endif
	
	r2count ++;
	if (r2count < M2)
	    Qstar = auxvector(0) * u.col(0) * u.col(0).transpose() + auxvector(1) * u.col(1) * u.col(1).transpose() + auxvector(2) * u.col(2) * u.col(2).transpose();
	else
	    Qstar = Matrix<double, NUMAXIS, NUMAXIS>::Zero();
    }
    
    /** Form the Rk matrix **/
    for (int i = 0; i<NUMBER_OF_WHEELS; i++)
    {
	Rk.block<NUMAXIS, NUMAXIS>(i*(2*NUMAXIS),i*(2*NUMAXIS)) = Rv; /** For the linear velocity **/
	Rk.block<NUMAXIS, NUMAXIS>(NUMAXIS+(i*(2*NUMAXIS)), NUMAXIS+(i*(2*NUMAXIS))) = Rg; /** For the angular velocity **/
    }
    
    Rk.block<NUMAXIS, NUMAXIS>(sckf::E_MEASUREMENT_VECTOR_SIZE,sckf::E_MEASUREMENT_VECTOR_SIZE) = Ra + Qstar; /** For the attitude correction **/
    
    /** Compute the Kalman Gain Matrix **/
    K = Pki_k * Hk.transpose() * (Hk * Pki_k * Hk.transpose() + Rk).inverse();
    
    #ifdef DEBUG_PRINTS
    std::cout<<"[Update] Hk*xki_k:\n"<<(Hk*xki_k)<<"\n";
    std::cout<<"[Update] innovation:\n"<<(zki - Hk*xki_k)<<"\n";
    #endif
    
    innovation = (zki - Hk*xki_k);
    
    /** Update the state vector and the covariance matrix **/
    xki_k = xki_k + K * (zki - Hk*xki_k);
        
    Pki_k = (Matrix<double,X_STATE_VECTOR_SIZE,X_STATE_VECTOR_SIZE>::Identity()-K*Hk)*Pki_k*(Matrix<double,X_STATE_VECTOR_SIZE,X_STATE_VECTOR_SIZE>::Identity()-K*Hk).transpose() + K*Rk*K.transpose();
    Pki_k = 0.5 * (Pki_k + Pki_k.transpose());
    
    #ifdef DEBUG_PRINTS
    std::cout<< "[Update] Qstar is of size "<<Qstar.rows()<<"x"<<Qstar.cols()<<"\n";
    std::cout<< "[Update] Qstar:\n"<<Qstar<<"\n";
    std::cout<< "[Update] Rk is of size "<<Rk.rows()<<"x"<<Rk.cols()<<"\n";
    std::cout<< "[Update] Rk:\n"<<Rk<<"\n";
    std::cout<< "[Update] K is of size "<<K.rows()<<"x"<<K.cols()<<"\n";
    std::cout<< "[Update] K:\n"<<K<<"\n";
    std::cout<< "[Update] xki_k is of size "<<xki_k.rows()<<"x"<<xki_k.cols()<<"\n";
    std::cout<< "[Update] xki_k:\n"<<xki_k<<"\n";
    std::cout<< "[Update] Pki_k is of size "<<Pki_k.rows()<<"x"<<Pki_k.cols()<<"\n";
    std::cout<< "[Update] Pki_k:\n"<<Pki_k<<"\n";
    #endif
         
    /** Update the quaternion with the Indirect approach **/
    qe.w() = 1;
    qe.x() = xki_k((sckf::E_STATE_VECTOR_SIZE*sckf::NUMBER_OF_WHEELS));
    qe.y() = xki_k((sckf::E_STATE_VECTOR_SIZE*sckf::NUMBER_OF_WHEELS)+1);
    qe.z() = xki_k((sckf::E_STATE_VECTOR_SIZE*sckf::NUMBER_OF_WHEELS)+2);
    q4 = q4 * qe;
    
    #ifdef DEBUG_PRINTS
    Eigen::Matrix <double,NUMAXIS,1> euler; /** In euler angles **/
    euler[2] = qe.toRotationMatrix().eulerAngles(2,1,0)[0];//YAW
    euler[1] = qe.toRotationMatrix().eulerAngles(2,1,0)[1];//PITCH
    euler[0] = qe.toRotationMatrix().eulerAngles(2,1,0)[2];//ROLL
    std::cout<< "[Update] Error quaternion in euler\n";
    std::cout<< "[Update] Roll: "<<euler[0]*R2D<<" Pitch: "<<euler[1]*R2D<<" Yaw: "<<euler[2]*R2D<<"\n";
    #endif
    
    /** Normalize quaternion **/
    q4.normalize();

    
    /** Reset the quaternion part of the state vector **/
    xki_k.block<NUMAXIS,1>((sckf::E_STATE_VECTOR_SIZE*sckf::NUMBER_OF_WHEELS),0) = Matrix<double, NUMAXIS, 1>::Zero();
    
    /**---------------------------- **/
    /** Reset the rest of the state **/
    /**---------------------------- **/
    for (int i = 0; i<sckf::NUMBER_OF_WHEELS; i++)
    {
	slipMatrix.col(i) = xki_k.block<NUMAXIS, 1> (i*NUMAXIS,0);
	xki_k.block<NUMAXIS, 1> (i*NUMAXIS,0) = Eigen::Matrix <double, NUMAXIS, 1>::Zero();
    }
    
    acontact = xki_k.block<NUMBER_OF_WHEELS, 1> (NUMAXIS*NUMBER_OF_WHEELS,0);
    xki_k.block<sckf::NUMBER_OF_WHEELS, 1> (NUMAXIS*sckf::NUMBER_OF_WHEELS,0) = Eigen::Matrix <double, sckf::NUMBER_OF_WHEELS, 1>::Zero();
    
    bghat = bghat + xki_k.block<NUMAXIS, 1> ((sckf::E_STATE_VECTOR_SIZE*sckf::NUMBER_OF_WHEELS) + NUMAXIS,0);
    xki_k.block<NUMAXIS, 1> ((sckf::E_STATE_VECTOR_SIZE*sckf::NUMBER_OF_WHEELS) + NUMAXIS,0) = Matrix <double, NUMAXIS, 1>::Zero();
    
    bahat = bahat + xki_k.block<NUMAXIS, 1> ((sckf::E_STATE_VECTOR_SIZE*sckf::NUMBER_OF_WHEELS) + (2*NUMAXIS),0);
    xki_k.block<NUMAXIS, 1> ((sckf::E_STATE_VECTOR_SIZE*sckf::NUMBER_OF_WHEELS) + (2*NUMAXIS),0) = Matrix <double, NUMAXIS, 1>::Zero();
    
    #ifdef DEBUG_PRINTS
    std::cout<< "[Update After Reset] xki_k is of size "<<xki_k.rows()<<"x"<<xki_k.cols()<<"\n";
    std::cout<< "[Update After Reset] xki_k:\n"<<xki_k<<"\n";
    #endif
    
    return;

}


/**
* @brief This computes the theoretical gravity value according to the WGS-84 ellipsoid earth model.
*/
double localization::GravityModel(double latitude, double altitude)
{
    double g; /**< g magnitude at zero altitude **/

    /** Nominal Gravity model **/
    g = GWGS0*((1+GWGS1*pow(sin(latitude),2))/sqrt(1-pow(ECC,2)*pow(sin(latitude),2)));

    /** Gravity affects by the altitude (aprox the value r = Re **/
    g = g*pow(Re/(Re+altitude), 2);

    std::cout<<"Theoretical gravity for this location (WGS-84 ellipsoid model): "<< g<<" [m/s^2]\n";

    return g;

}

/**
* @brief Substract the Earth rotation from the gyroscopes readout
*/
void localization::SubstractEarthRotation(Eigen::Matrix <double, NUMAXIS, 1> *u, Eigen::Quaternion <double> *qb_g, double latitude)
{
    Eigen::Matrix <double, NUMAXIS, 1> v (EARTHW*cos(latitude), 0, EARTHW*sin(latitude)); /**< vector of earth rotation components expressed in the geografic frame according to the latitude **/

    /** Compute the v vector expressed in the body frame **/
    v = (*qb_g) * v;
    
    #ifdef DEBUG_PRINTS
    std::cout<<"Earth Rotation:"<<v<<"\n";
    #endif

    /** Subtract the earth rotation to the vector of inputs (u = u-v**/
    (*u)  = (*u) - v;
    
    return;
}

/**
* @brief Correct the magnetic declination of the North 
*/
int localization::CorrectMagneticDeclination(Eigen::Quaternion< double >* quat, double magnetic_declination, int mode)
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
	return ERROR;
    }
    
    *quat = Eigen::Quaternion <double> (Eigen::AngleAxisd(euler[0], Eigen::Vector3d::UnitX())*
			Eigen::AngleAxisd(euler[1], Eigen::Vector3d::UnitY()) *
			Eigen::AngleAxisd(euler[2], Eigen::Vector3d::UnitZ()));
    
    return OK;
}

