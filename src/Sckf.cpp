/**\file Sckf.cpp
 *
 * This class has the primitive methods for an Stochastic Cloning Indirect Kalman Filter implementation
 * the state vector are formed by the errors. Therefore the name of indirect kalman filter 
 * 
 * 
 * @author Javier Hidalgo Carrio | DFKI RIC Bremen | javier.hidalgo_carrio@dfki.de
 * @date November 2012.
 * @version 1.0.
 */

#include "Sckf.hpp"

#define DEBUG_PRINTS 1

using namespace localization;

void Sckf::welcome()
{
	std::cout << "You successfully compiled and executed SCFK. Welcome!" << std::endl;
}


/**
* @brief Gets the current vector x
*/
Eigen::Matrix< double, Eigen::Dynamic, 1  > Sckf::getStatex()
{
    return xki_k;
}

/**
* @brief Gets the current orientation of the robot in Quaternion
*/
Eigen::Quaternion< double > Sckf::getAttitude()
{
    return this->q4;
}

/**
* @brief Gets the gravity value
*/
double Sckf::getGravity()
{
    return this->gtilde.norm();
}


/**
* @brief Gets the current orientation of the robot in Euler angles
*/
Eigen::Matrix< double, NUMAXIS , 1  > Sckf::getEuler()
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
Eigen::Matrix< double, Eigen::Dynamic, Eigen::Dynamic > Sckf::getCovariancex()
{
    return Pki_k;
}

/**
* @brief Gets Noise covariance matrix for attitude estimation
*/
Eigen::Matrix< double, Sckf::A_STATE_VECTOR_SIZE, Sckf::A_STATE_VECTOR_SIZE > Sckf::getCovarianceAttitude()
{
    return Pki_k.block<Sckf::A_STATE_VECTOR_SIZE, Sckf::A_STATE_VECTOR_SIZE> (2*NUMAXIS, 2*NUMAXIS);
}

/**
* @brief Return the K matrix
*/
Eigen::Matrix< double, Eigen::Dynamic, Eigen::Dynamic > Sckf::getKalmanGain()
{
    return this->K;
}

/**
* @brief Get linear velocity covariance from IMU
*/
Eigen::Matrix< double, NUMAXIS, NUMAXIS > Sckf::getLinearCovarianceVelocities(double dt)
{
   return (this->Raup * dt);
}


/**
* @brief Return the K associated to the attitude
*/
Eigen::Matrix< double, Eigen::Dynamic, Eigen::Dynamic > Sckf::getAttitudeKalmanGain()
{
    return this->K.block<Sckf::A_STATE_VECTOR_SIZE, A_STATE_VECTOR_SIZE> (2*NUMAXIS, 2*NUMAXIS);
}

/**
* @brief Return the innovation vector
*/
Eigen::Matrix< double, Eigen::Dynamic, 1  > localization::Sckf::getInnovation()
{
    return this->innovation;
}

/**
* @brief Gets the delta orientation in Quaternion
*/
Eigen::Quaternion< double > Sckf::deltaQuaternion()
{
    return prev_q4.inverse() * q4;
}


/**
* @brief This function Initialize Attitude
*/
bool Sckf::setAttitude(Eigen::Quaternion< double >& initq)
{
    /** Initial orientation **/
    q4 = initq;
    prev_q4 = initq;
	
    return true;
    
}

/**
* @brief This function sets the gravity value
*/
void Sckf::setGravity(double g)
{
    this->gtilde << 0, 0, g;
    #ifdef DEBUG_PRINTS
    std::cout<<"Set gravity to: "<<gtilde[2]<<"\n";
    #endif
}

/**
* @brief This function set the initial Omega matrix
*/
bool Sckf::setOmega(Eigen::Matrix< double, NUMAXIS , 1  >& u)
{
    if (&u != NULL)
    {
	/** Initialization for quaternion integration **/
	oldomega4 << 0,-u(0), -u(1), -u(2),
	    u(0), 0, u(2), -u(1),
	    u(1), -u(2), 0, u(0),
	    u(2), u(1), -u(0), 0;
	
	return true;
    }

    return false;

}

/**
* @brief Gets the current state vector of the filter
*/
void Sckf::setStatex(Eigen::Matrix< double, Eigen::Dynamic, 1  > &x_0)
{
    x_0.resize(Sckf::X_STATE_VECTOR_SIZE,1);
    this->xki_k = x_0;
}

/**
* @brief This function Initialize the Bias offset
*/
void Sckf::setBiasOffset(Eigen::Matrix< double, NUMAXIS , 1  > gbias, Eigen::Matrix< double, NUMAXIS , 1  > abias)
{
    this->bghat = gbias;
    this->bahat = abias;
}


/**
* @brief Set the Heading angle
*/
void Sckf::setHeading(double yaw)
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


 /**
* @brief This function Initilize the vectors and matrices
*/
void Sckf::Init(Eigen::Matrix< double, Eigen::Dynamic, Eigen::Dynamic >& P_0,
		Eigen::Matrix< double, NUMAXIS , NUMAXIS  >& Rg,
		Eigen::Matrix< double, NUMAXIS , NUMAXIS  >& Qbg,
		Eigen::Matrix< double, NUMAXIS , NUMAXIS  >& Qba, 
		Eigen::Matrix< double, NUMAXIS , NUMAXIS  >& Rapr,
		Eigen::Matrix< double, NUMAXIS , NUMAXIS  >& Raup,
		Eigen::Matrix< double, NUMAXIS , NUMAXIS  >& Rat,
		Eigen::Matrix< double, NUMAXIS , NUMAXIS  >& Rm,
		double g, double alpha)
{
    
    /** Set the matrix and vector dimension to the static values of the class **/
    xk_k.resize(Sckf::X_STATE_VECTOR_SIZE,1);
    xki_k.resize(Sckf::X_STATE_VECTOR_SIZE,1);
    
    A.resize(Sckf::A_STATE_VECTOR_SIZE,Sckf::A_STATE_VECTOR_SIZE);
    Fki.resize(Sckf::X_STATE_VECTOR_SIZE,Sckf::X_STATE_VECTOR_SIZE);
    
    Qk.resize(Sckf::X_STATE_VECTOR_SIZE,Sckf::X_STATE_VECTOR_SIZE);
    
    Pk_k.resize(Sckf::X_STATE_VECTOR_SIZE,Sckf::X_STATE_VECTOR_SIZE);
    Pki_k.resize(Sckf::X_STATE_VECTOR_SIZE,Sckf::X_STATE_VECTOR_SIZE);
    
    K.resize(Sckf::X_STATE_VECTOR_SIZE, 2*NUMAXIS);
    
    H1a.resize(NUMAXIS, Sckf::A_STATE_VECTOR_SIZE);
    H2a.resize(NUMAXIS, Sckf::A_STATE_VECTOR_SIZE);
    Hk.resize(2*NUMAXIS, Sckf::X_STATE_VECTOR_SIZE);
    
    /** Resizing dynamic arguments to the correct dimension to avoid matrices errors **/
    P_0.resize(Sckf::X_STATE_VECTOR_SIZE, Sckf::X_STATE_VECTOR_SIZE);
    
    /** Gravitation acceleration **/
    gtilde << 0, 0, g;

    /** Dip angle (alpha is in rad) **/
    mtilde(0) = cos(alpha);
    mtilde(1) = 0;
    mtilde(2) = -sin(alpha);

    /** Kalman filter state, system matrix, error covariance and process noise covariance **/
    xki_k = Eigen::Matrix <double,Sckf::X_STATE_VECTOR_SIZE,1>::Zero();
    xk_k = xki_k;
    
    /** System matrix F **/
    Fki = Eigen::Matrix <double,Sckf::X_STATE_VECTOR_SIZE, Sckf::X_STATE_VECTOR_SIZE>::Zero();
    
    /** System matrix A **/
    A = Eigen::Matrix <double,Sckf::A_STATE_VECTOR_SIZE, Sckf::A_STATE_VECTOR_SIZE>::Zero();
    A(0,3) = -0.5;A(1,4) = -0.5;A(2,5) = -0.5;
    
    /** Process noise **/
    Qk = Eigen::Matrix <double,Sckf::X_STATE_VECTOR_SIZE,Sckf::X_STATE_VECTOR_SIZE>::Zero();
    Qk.block <NUMAXIS, NUMAXIS> (NUMAXIS,NUMAXIS) = Rapr;
    Qk.block <NUMAXIS, NUMAXIS> (2*NUMAXIS,2*NUMAXIS) = 0.25 * Rg;
    Qk.block <NUMAXIS, NUMAXIS> ((2*NUMAXIS)+NUMAXIS,(2*NUMAXIS)+NUMAXIS) = Qbg;
    Qk.block <NUMAXIS, NUMAXIS> ((4*NUMAXIS),(4*NUMAXIS)) = Qba;

    /** Assign the initial values **/
    Pki_k = P_0;
    Pk_k = Pki_k;
    
    /** Assign the initial value for the measurement matrix of the attitude **/
    H1a = Eigen::Matrix <double,NUMAXIS,Sckf::A_STATE_VECTOR_SIZE>::Zero();
    H2a = Eigen::Matrix <double,NUMAXIS,Sckf::A_STATE_VECTOR_SIZE>::Zero();
    H1a(0,6) = 1; H1a(1,7) = 1; H1a(2,8) = 1;

    /** Set the history of noise for the attitude **/
    RHist.resize(M1);
    
    /** Fill noise measurement matrices Rg, Ra, Rat, Rm , etc.. **/
    this->Rg = Rg;
    this->Rapr = Rapr;
    this->Raup = Raup;
    this->Rat = Rat;
    this->Rm = Rm;
    
    /** Resize the noise measurement matrix to the correct dimension **/
    Rk.resize (2*NUMAXIS, 2*NUMAXIS);
    
    /** Measurement vector **/
    zki.resize(2*NUMAXIS,1);
    zki = Eigen::Matrix<double, 2*NUMAXIS, 1>::Zero();
    
    /** Innovation **/
    innovation.resize(2*NUMAXIS, 1);
    innovation.setZero();
    
    /** Initial bias **/
    bghat = Eigen::Matrix <double,NUMAXIS,1>::Zero();
    bahat = Eigen::Matrix <double,NUMAXIS,1>::Zero();
        
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
    
    /** Initialize the previous q4 **/
    prev_q4 = q4;
    
    /** Default initial position errors **/
    eposition<<0.00, 0.00, 0.00;
    
    /** Default initial velocity errors **/
    evelocity<<0.00, 0.00, 0.00;
    evelocity_cov.setZero();
    filtered_evelocity.setZero();
//     evelocity = DataModel(NUMAXIS);
    
    /** Default initial bias **/
    bghat << 0.00, 0.00, 0.00;
    bahat << 0.00, 0.00, 0.00;
    
    /** Variable in the adaptive algorithm **/
    r1count = 0;
    r2count = R2COUNT;
    
//     veloModelk_1 = DataModel(NUMAXIS); //! Initialize vel from model at k-1
//     veloTruth = DataModel(NUMAXIS);
//     increVeloError = DataModel(NUMAXIS);
        
    /** Print filter information **/
    #ifdef DEBUG_PRINTS
    std::cout<< "[FILTER_INIT] xki_k is of size "<<xki_k.rows()<<"x"<<xki_k.cols()<<"\n";
    std::cout<< "[FILTER_INIT] xki_k:\n"<<xki_k<<"\n";
    std::cout<< "[FILTER_INIT] xk_k is of size "<<xk_k.rows()<<"x"<<xk_k.cols()<<"\n";
    std::cout<< "[FILTER_INIT] xk_k:\n"<<xk_k<<"\n";
    std::cout<< "[FILTER_INIT] A is of size "<<A.rows()<<"x"<<A.cols()<<"\n";
    std::cout<< "[FILTER_INIT] A:\n"<<A<<"\n";
    std::cout<< "[FILTER_INIT] Fki is of size "<<Fki.rows()<<"x"<<Fki.cols()<<"\n";
    std::cout<< "[FILTER_INIT] Fki:\n"<<Fki<<"\n";
    std::cout<< "[FILTER_INIT] Pk+i|k is of size "<<Pki_k.rows()<<"x"<<Pki_k.cols()<<"\n";
    std::cout<< "[FILTER_INIT] Pk+i|k:\n"<<Pki_k<<"\n";
    std::cout<< "[FILTER_INIT] Pk|k is of size "<<Pk_k.rows()<<"x"<<Pk_k.cols()<<"\n";
    std::cout<< "[FILTER_INIT] Pk|k:\n"<<Pk_k<<"\n";
    std::cout<< "[FILTER_INIT] Qk|k is of size "<<Qk.rows()<<"x"<<Qk.cols()<<"\n";
    std::cout<< "[FILTER_INIT] Qk|k:\n"<<Qk<<"\n";
    std::cout<< "[FILTER_INIT] H1a is of size "<<H1a.rows()<<"x"<<H1a.cols()<<"\n";
    std::cout<< "[FILTER_INIT] H1a:\n"<<H1a<<"\n";
    std::cout<< "[FILTER_INIT] H2a is of size "<<H2a.rows()<<"x"<<H2a.cols()<<"\n";
    std::cout<< "[FILTER_INIT] H2a:\n"<<H2a<<"\n";
    std::cout<< "[FILTER_INIT] zki is of size "<<zki.rows()<<"x"<<zki.cols()<<"\n";
    std::cout<< "[FILTER_INIT] zki:\n"<<zki<<"\n";
    std::cout<< "[FILTER_INIT] RHist is of size "<<RHist.size()<<"\n";
    std::cout<< "[FILTER_INIT] Rg is of size "<<Rg.rows()<<"x"<<Rg.cols()<<"\n";
    std::cout<< "[FILTER_INIT] Rg:\n"<<Rg<<"\n";
    std::cout<< "[FILTER_INIT] Rapr is of size "<<Rapr.rows()<<"x"<<Rapr.cols()<<"\n";
    std::cout<< "[FILTER_INIT] Rapr:\n"<<Rapr<<"\n";
    std::cout<< "[FILTER_INIT] Raup is of size "<<Raup.rows()<<"x"<<Raup.cols()<<"\n";
    std::cout<< "[FILTER_INIT] Raup:\n"<<Raup<<"\n";
    std::cout<< "[FILTER_INIT] Rat is of size "<<Rat.rows()<<"x"<<Rat.cols()<<"\n";
    std::cout<< "[FILTER_INIT] Rat:\n"<<Rat<<"\n";
    std::cout<< "[FILTER_INIT] Rm is of size "<<Rm.rows()<<"x"<<Rm.cols()<<"\n";
    std::cout<< "[FILTER_INIT] Rm:\n"<<Rm<<"\n";
    std::cout<< "[FILTER_INIT] mtilde is of size "<<mtilde.rows()<<"x"<<mtilde.cols()<<"\n";
    std::cout<< "[FILTER_INIT] mtilde:\n"<<mtilde<<"\n";
    std::cout<< "[FILTER_INIT] gtilde is of size "<<gtilde.rows()<<"x"<<gtilde.cols()<<"\n";
    std::cout<< "[FILTER_INIT] gtilde:\n"<<gtilde<<"\n";
    std::cout<< "[FILTER_INIT] bghat is of size "<<bghat.rows()<<"x"<<bghat.cols()<<"\n";
    std::cout<< "[FILTER_INIT] bghat:\n"<<bghat<<"\n";
    std::cout<< "[FILTER_INIT] bahat is of size "<<bahat.rows()<<"x"<<bahat.cols()<<"\n";
    std::cout<< "[FILTER_INIT] bahat:\n"<<bahat<<"\n";
    #endif

}

/**
* @brief Performs the prediction step of the filter.
*/
void Sckf::predict(Eigen::Matrix< double, NUMAXIS , 1  >& u, Eigen::Matrix< double, NUMAXIS , 1  >& v, double dt)
{
    Eigen::Matrix <double,NUMAXIS,NUMAXIS> Cq; /** Rotational matrix */
    Eigen::Matrix <double,NUMAXIS,NUMAXIS> velo2product; /** Vec 2 product  matrix */
    Eigen::Matrix <double,NUMAXIS,NUMAXIS> acc2product; /** Vec 2 product  matrix */
    Eigen::Matrix <double,NUMAXIS,1> angvelo; /** Angular velocity */
    Eigen::Matrix <double,NUMAXIS,1> linacc; /** Linear acceleration */
    Eigen::Matrix <double,NUMAXIS,1> gtilde_body; /** Gravitation in the body frame */
    Eigen::Matrix <double,QUATERSIZE,QUATERSIZE> omega4; /** Quaternion integration matrix */
    Eigen::Matrix <double,QUATERSIZE,1> quat; /** Quaternion integration matrix */
    Eigen::Matrix <double,Sckf::X_STATE_VECTOR_SIZE,Sckf::X_STATE_VECTOR_SIZE> dFki; /** Discrete System matrix */
    Eigen::Matrix <double,X_STATE_VECTOR_SIZE,X_STATE_VECTOR_SIZE> Qdk; /** Discrete Qk matrix */

    /** Compute the cross product matrix with the angular velocity **/
    angvelo = u - bghat; /** Eliminate the Bias **/

    /** In cross product form **/
    velo2product << 0, -angvelo(2), angvelo(1),
		angvelo(2), 0, -angvelo(0),
		-angvelo(1), angvelo(0), 0;
		
    /** Create the orientation matrix from the quaternion (q_body2world) **/
    Quaternion2DCM (&q4, &Cq);
    
    /** Calculate the gravity vector in the body frame **/
    gtilde_body = Cq * gtilde;
		
    /** Compute the cross product matrix with the linear acceleration **/
    linacc = v - bahat - gtilde_body; /** Eliminate the Bias and the local gravity vector **/
    
    #ifdef DEBUG_PRINTS
    std::cout<<"[Predict] angevelo:\n"<<angvelo<<"\n";
    std::cout<<"[Predict] gtilde_body of size "<<gtilde_body.rows()<<"x"<<gtilde_body.cols()<<"\n";
    std::cout<<"[Predict] g in body_frame:\n"<<gtilde_body<<"\n";
    std::cout<<"[Predict] linacc:\n"<<linacc<<"\n";
    #endif

    /** In cross product form **/
    acc2product << 0, -linacc(2), linacc(1),
		linacc(2), 0, -linacc(0),
		-linacc(1), linacc(0), 0;
		
    /** Compute the dA Matrix of the attitude part **/
    A.block<NUMAXIS, NUMAXIS> (0,0) = -velo2product;
    
    /** Form the complete system model matrix (position error) **/
    Fki.block<NUMAXIS, NUMAXIS> (0,NUMAXIS) = Cq.inverse();//Eigen::Matrix <double,NUMAXIS, NUMAXIS>::Identity();//!Cq.inverse() if velocity is in body_frame
    
    /** Velocity part **/
    Fki.block<NUMAXIS, NUMAXIS> (NUMAXIS, 2*NUMAXIS) = -/*Cq.inverse() */ acc2product;
    Fki.block<NUMAXIS, NUMAXIS> (NUMAXIS, 4*NUMAXIS) = -Eigen::Matrix <double,NUMAXIS, NUMAXIS>::Identity();//-Cq.inverse();
    
    /** Attitude part **/
    Fki.block<Sckf::A_STATE_VECTOR_SIZE, Sckf::A_STATE_VECTOR_SIZE>(2*NUMAXIS, 2*NUMAXIS) = A;

    /** Discretization of the linear system **/
    dFki = Eigen::Matrix<double,Sckf::X_STATE_VECTOR_SIZE,Sckf::X_STATE_VECTOR_SIZE>::Identity() + Fki * dt + Fki * Fki * pow(dt,2)/2.0;
    
    #ifdef DEBUG_PRINTS
    std::cout<< "[Predict] xki|k is of size "<<xki_k.rows()<<"x"<<xki_k.cols()<<"\n";
    std::cout<< "[Predict] xki_k:\n"<<xki_k<<"\n";
    #endif
    
    /** Propagate the vector through the system **/
    xki_k = dFki * xki_k;
    
    /** The process noise covariance matrix **/
    Qk.block <NUMAXIS, NUMAXIS> (0,0) = Cq.inverse() * Rapr * dt;
    Qk.block <NUMAXIS, NUMAXIS> (NUMAXIS,NUMAXIS) = /*Cq.inverse() */ this->Rapr;
    
    /** Form the system noise matrix (discretization) **/
    Qdk = Qk*dt + 0.5*dt*Fki*Qk + 0.5*dt*Qk*Fki.transpose();
    Qdk = 0.5*(Qdk + Qdk.transpose());
    
    /** Propagate the P covariance matrix **/
    Pki_k = dFki*Pki_k*dFki.transpose() + Qdk;
    
    /** Discrete quaternion integration of the angular velocity **/
    omega4 << 0,-angvelo(0), -angvelo(1), -angvelo(2),
	    angvelo(0), 0, angvelo(2), -angvelo(1),
	    angvelo(1), -angvelo(2), 0, angvelo(0),
	    angvelo(2), angvelo(1), -angvelo(0), 0;
	    
    /** Copy the quaternion in the previous one before performing the rotation given by the gyros **/
    prev_q4 = q4;
	    
    /** Copy in a vector form **/
    quat(0) = q4.w();
    quat(1) = q4.x();
    quat(2) = q4.y();
    quat(3) = q4.z();

    /** Quaternion integration **/
    quat = (Eigen::Matrix<double,QUATERSIZE,QUATERSIZE>::Identity() +(0.75 * omega4 *dt)- (0.25 * oldomega4 * dt) -
    ((1/6) * angvelo.squaredNorm() * pow(dt,2) *  Eigen::Matrix<double,QUATERSIZE,QUATERSIZE>::Identity()) -
    ((1/24) * omega4 * oldomega4 * pow(dt,2)) - ((1/48) * angvelo.squaredNorm() * omega4 * pow(dt,3))) * quat;

    /** Store in a quaternion form **/
    q4.w() = quat(0);
    q4.x() = quat(1);
    q4.y() = quat(2);
    q4.z() = quat(3);
    q4.normalize();

    /** Copy omega **/
    oldomega4 = omega4;
    
    /** Store the corrected values in the proprioceptive measurement object **/
    measurement.setInertialValues(linacc, angvelo);
    
    #ifdef DEBUG_PRINTS
    std::cout<< "[After Predict] xki_k:\n"<<xki_k<<"\n";
    std::cout<< "[Predict] Fki is of size "<<Fki.rows()<<"x"<<Fki.cols()<<"\n";
    std::cout<< "[Predict] Fki:\n"<<Fki<<"\n";
    std::cout<< "[Predict] dFki is of size "<<dFki.rows()<<"x"<<dFki.cols()<<"\n";
    std::cout<< "[Predict] dFki:\n"<<dFki<<"\n";
    std::cout<< "[Predict] Qdk is of size "<<Qdk.rows()<<"x"<<Qdk.cols()<<"\n";
    std::cout<< "[Predict] Qdk:\n"<<Qdk<<"\n";
    std::cout<< "[Predict] Pki_k:\n"<<Pki_k<<"\n";
    #endif
    
    #ifdef DEBUG_PRINTS
    std::cout<<"[Attitude] Inputs u(angular_velocity):\n"<<u<<"\nv(acc):\n"<<v<<"\n";
    std::cout<<"[Attitude] angevelo:\n"<<angvelo<<"\n";
    std::cout<<"[Attitude] gtilde_body of size "<<gtilde_body.rows()<<"x"<<gtilde_body.cols()<<"\n";
    std::cout<<"[Attitude] g in body_frame:\n"<<gtilde_body<<"\n";
    std::cout<<"[Attitude] g in body_frame(quat):\n"<<q4*gtilde<<"\n";
    std::cout<<"[Attitude] g in body_frame(quat.inverse):\n"<<q4.inverse()*gtilde<<"\n";
    std::cout<<"[Attitude] g in body_frame(Rot):\n"<<q4.toRotationMatrix()*gtilde<<"\n";
    std::cout<<"[Attitude] g in body_frame(Rot.inverse):\n"<<q4.toRotationMatrix().inverse()*gtilde<<"\n";
    std::cout<<"[Attitude] g in world_frame(Cq.inverse):\n"<<Cq.inverse() * gtilde_body<<"\n";
    std::cout<<"[Attitude] linacc:\n"<<linacc<<"\n";
    #endif
    
    return;

}


void Sckf::update(Eigen::Matrix <double,NUMAXIS,NUMAXIS> &Hme, Eigen::Matrix <double,NUMAXIS,NUMAXIS> &Rme,
		  Eigen::Matrix< double, NUMAXIS, 1  >& slip_error,
		  Eigen::Matrix< double, NUMAXIS , 1  >& acc,Eigen::Matrix< double, NUMAXIS , 1  >& mag, double dt, bool magn_on_off)
{
    Eigen::Matrix <double,NUMAXIS,NUMAXIS> Cq; /** Rotational matrix */
    Eigen::Matrix <double,NUMAXIS,NUMAXIS> gtilde2product; /** Vec 2 product  matrix for the gravity vector in body frame */
    Eigen::Matrix <double,NUMAXIS,NUMAXIS> fooR2; /**  Measurement noise matrix from accelerometers matrix Ra */
    Eigen::Matrix <double,A_STATE_VECTOR_SIZE,1> xa_k; /** Attitude part of the state vector xk+i|k */
    Eigen::Matrix <double,A_STATE_VECTOR_SIZE,A_STATE_VECTOR_SIZE> P1a; /** Error convariance matrix for measurement 1 of the attitude */
    Eigen::Matrix <double,A_STATE_VECTOR_SIZE,A_STATE_VECTOR_SIZE> P2a; /** Error convariance matrix for measurement 2 of the attitude */
    Eigen::Matrix <double,A_STATE_VECTOR_SIZE,A_STATE_VECTOR_SIZE> auxM; /** Auxiliar matrix for computing Kalman gain in measurement */
    Eigen::Matrix <double,NUMAXIS,NUMAXIS> Rve; /** Measurement noise convariance matrix for velocity error */
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
    Eigen::Matrix <double,NUMAXIS,1> auxvector; /** Auxiliar vector variable */
    Eigen::Matrix <double,NUMAXIS,1> z1a; /** Measurement vector 1 Acc */
    Eigen::Matrix <double,NUMAXIS,1> z2a; /** Measurement vector 2 Mag */
    
    
    /** Copy the attitude part of the state vector and covariance matrix **/
    xa_k = xki_k.block<Sckf::A_STATE_VECTOR_SIZE, 1> (2*NUMAXIS, 0);
    P1a = Pki_k.block<Sckf::A_STATE_VECTOR_SIZE, Sckf::A_STATE_VECTOR_SIZE> (2*NUMAXIS, 2*NUMAXIS);
    
    #ifdef DEBUG_PRINTS
    std::cout<<"[Update] xa_k is of size "<<xa_k.rows()<<"x"<<xa_k.cols()<<"\n";
    std::cout<<"[Update] xa_k:\n"<<xa_k<<"\n";
    std::cout<<"[Update] P1a is of size "<<P1a.rows()<<"x"<<P1a.cols()<<"\n";
    std::cout<<"[Update] P1a:\n"<<P1a<<"\n";
    #endif
    
    /** Create the orientation matrix from the quaternion **/
    Quaternion2DCM (&q4, &Cq);
    
    /** Calculate the gravity vector in the body frame **/
    gtilde_body = Cq * gtilde;
    
    gtilde2product << 0, -gtilde_body(2), gtilde_body(1),
		    gtilde_body(2), 0, -gtilde_body(0),
		    -gtilde_body(1), gtilde_body(0), 0;

    #ifdef DEBUG_PRINTS
    std::cout<<"[Update] gtilde_body is of size "<<gtilde_body.rows()<<"x"<<gtilde_body.cols()<<"\n";
    std::cout<<"[Update] g in body_frame:\n"<<gtilde_body<<"\n";
    #endif

    /** Form the matrix for the measurement 1 of the attitude (acc correction) **/
    H1a.block<NUMAXIS, NUMAXIS> (0,0) = 2*gtilde2product;
    
    /** Form the measurement vector z1a for the attitude (the real acceleration value) **/
    z1a = acc;// - bahat - gtilde_body;
    
    /** The adaptive algorithm for the attitude, the Uk matrix and SVD part **/
    R1a = (z1a - H1a*xa_k) * (z1a - H1a*xa_k).transpose();
    RHist[r1count] = R1a;
    
    /** r1count + 1 modulus the number of history M1 **/
    r1count = (r1count+1)%(M1); 

    Uk.setZero();
    /** Starting the Uk is R **/
    for (register int j=0; j<M1; j++)
    {
	Uk += RHist[j];
    }
    
    Uk = Uk/static_cast<double>(M1);
    
    fooR2 = H1a*P1a*H1a.transpose() + Raup;
    
    /**
    * Single Value Decomposition
    */
    Eigen::JacobiSVD <Eigen::MatrixXd > svdOfUk(Uk, Eigen::ComputeThinU);

    s = svdOfUk.singularValues(); //!eigenvalues
    u = svdOfUk.matrixU();//!eigenvectors
    
    lambda << s(0), s(1), s(2);
    
    mu(0) = u.col(0).transpose() * fooR2 * u.col(0);
    mu(1) = u.col(1).transpose() * fooR2 * u.col(1);
    mu(2) = u.col(2).transpose() * fooR2 * u.col(2);
    
    #ifdef DEBUG_PRINTS
    std::cout<<"[Update] (lambda - mu) is: "<<(lambda - mu)<<"\n";
    #endif
    
    if ((lambda - mu).maxCoeff() > GAMMA)
    {
	
	#ifdef DEBUG_PRINTS
	std::cout<<"[Update] Bigger than GAMMA("<<GAMMA<<")\n";
	#endif
    
	r2count = 0;
	auxvector(0) = std::max(lambda(0)-mu(0),static_cast<double>(0.00));
	auxvector(1) = std::max(lambda(1)-mu(1),static_cast<double>(0.00));
	auxvector(2) = std::max(lambda(2)-mu(2),static_cast<double>(0.00));
	
	Qstar = auxvector(0) * u.col(0) * u.col(0).transpose() + auxvector(1) * u.col(1) * u.col(1).transpose() + auxvector(2) * u.col(2) * u.col(2).transpose();
    }
    else
    {
	#ifdef DEBUG_PRINTS
	std::cout<<"[Update] Lower than GAMMA("<<GAMMA<<") r2count: "<<r2count<<"\n";
	#endif
	
	r2count ++;
	if (r2count < M2)
	{
	    auxvector(0) = std::max(lambda(0)-mu(0),static_cast<double>(0.00));
	    auxvector(1) = std::max(lambda(1)-mu(1),static_cast<double>(0.00));
	    auxvector(2) = std::max(lambda(2)-mu(2),static_cast<double>(0.00));
	    
	    Qstar = auxvector(0) * u.col(0) * u.col(0).transpose() + auxvector(1) * u.col(1) * u.col(1).transpose() + auxvector(2) * u.col(2) * u.col(2).transpose();
	}
	else
	    Qstar = Eigen::Matrix<double, NUMAXIS, NUMAXIS>::Zero();
    }
    
    /** Measurement vector **/
    zki.block<NUMAXIS, 1> (0,0) = Hme * slip_error;//xki_k.block<NUMAXIS,1> (3,0);
    zki.block<NUMAXIS, 1> (NUMAXIS,0) = z1a;
    
    /** Form the observation matrix Hk **/
    Hk = Eigen::Matrix<double, 2*NUMAXIS, X_STATE_VECTOR_SIZE>::Zero();
    Hk.block<NUMAXIS, NUMAXIS>(0,NUMAXIS) = Eigen::Matrix <double,NUMAXIS, NUMAXIS>::Identity();//Cq;
    Hk.block<NUMAXIS, A_STATE_VECTOR_SIZE>(NUMAXIS,2*NUMAXIS) = H1a;
    
    /** Define the velocity vector measurement noise from the slip_error vector measurement**/
    Rve = Hme * Rme * Hme.transpose();
    
    /** Define the attitude measurement noise **/
    R1a = Rat + Qstar; //!Qstart is the external acceleration covariance
    
    /** Composition of the measurement noise matrix **/
    Rk.setZero();
    Rk.block<NUMAXIS, NUMAXIS> (0,0) = Rve;
    Rk.block<NUMAXIS, NUMAXIS> (NUMAXIS,NUMAXIS) = R1a;
    
    /** Set the innovation global variable **/
    innovation = (zki - Hk * xki_k);
    
    #ifdef DEBUG_PRINTS
    std::cout<< "[Update] zki is of size "<<zki.rows()<<"x"<<zki.cols()<<"\n";
    std::cout<< "[Update] zki:\n"<<zki<<"\n";
    std::cout<< "[Update] observation:\n"<<Hk * xki_k<<"\n";
    std::cout<< "[Update] innovation is of size "<<innovation.rows()<<"x"<<innovation.cols()<<"\n";
    std::cout<< "[Update] innovation:\n"<<innovation<<"\n";
    #endif
    
    /** Compute the Kalman Gain Matrix **/
    Eigen::Matrix<double, 2*NUMAXIS, 2*NUMAXIS> S = Hk * Pki_k * Hk.transpose() + Rk;
//     K = Pki_k * Hk.transpose() * ((S.transpose() * S).inverse() * S.transpose()); //!Calculte K using the pseudoinverse of S
    K = Pki_k * Hk.transpose() * S.inverse(); //!Calculate K using the inverse of S

    /** Update the state vector and the covariance matrix **/
    xki_k = xki_k + K * (zki - Hk * xki_k);
    Pki_k = (Eigen::Matrix<double,Sckf::X_STATE_VECTOR_SIZE,Sckf::X_STATE_VECTOR_SIZE>::Identity()-K*Hk)*Pki_k*(Eigen::Matrix<double,Sckf::X_STATE_VECTOR_SIZE,Sckf::X_STATE_VECTOR_SIZE>::Identity()-K*Hk).transpose() + K*Rk*K.transpose();
    Pki_k = 0.5 * (Pki_k + Pki_k.transpose());
    
    #ifdef DEBUG_PRINTS
    std::cout<< "[Update] RHist is of size "<<RHist.size()<<"r1count is: "<<r1count<<"\n";
    std::cout<< "[Update] Uk is of size "<<Uk.rows()<<"x"<<Uk.cols()<<"\n";
    std::cout<< "[Update] Uk:\n"<<Uk<<"\n";
    std::cout<< "[Update] fooR2 is of size "<<fooR2.rows()<<"x"<<fooR2.cols()<<"\n";
    std::cout<< "[Update] fooR2:\n"<<fooR2<<"\n";
    std::cout<< "[Update] Qstar is of size "<<Qstar.rows()<<"x"<<Qstar.cols()<<"\n";
    std::cout<< "[Update] Qstar:\n"<<Qstar<<"\n";
    std::cout<< "[Update] R1a is of size "<<R1a.rows()<<"x"<<R1a.cols()<<"\n";
    std::cout<< "[Update] R1a:\n"<<R1a<<"\n";
    std::cout<< "[Update] Rk is of size "<<Rk.rows()<<"x"<<Rk.cols()<<"\n";
    std::cout<< "[Update] Rk:\n"<<Rk<<"\n";
    std::cout<< "[Update] Hk is of size "<<Hk.rows()<<"x"<<Hk.cols()<<"\n";
    std::cout<< "[Update] Hk:\n"<<Hk<<"\n";
    std::cout<< "[Update] xki_k is of size "<<xki_k.rows()<<"x"<<xki_k.cols()<<"\n";
    std::cout<< "[Update] xki_k:\n"<<xki_k<<"\n";
    std::cout<< "[Update] Pki_k is of size "<<Pki_k.rows()<<"x"<<Pki_k.cols()<<"\n";
    std::cout<< "[Update] Pki_k:\n"<<Pki_k<<"\n";
    std::cout<< "[Update] K is of size "<<K.rows()<<"x"<<K.cols()<<"\n";
    std::cout<< "[Update] K:\n"<<K<<"\n";
    
    #endif
    
    /** Update the quaternion with the Indirect approach **/
    qe.w() = 1;
    qe.x() = xki_k(6);
    qe.y() = xki_k(7);
    qe.z() = xki_k(8);
    
//     Eigen::Matrix <double,NUMAXIS,1> error_euler; /** The error in euler angles **/
//     error_euler[2] = 0.00; //Dont update the yaw
//     error_euler[1] = qe.toRotationMatrix().eulerAngles(2,1,0)[1];//PITCH
//     error_euler[0] = qe.toRotationMatrix().eulerAngles(2,1,0)[2];//ROLL
//     
//     qe = Eigen::Quaternion <double> (Eigen::AngleAxisd(error_euler[2], Eigen::Vector3d::UnitZ())*
// 			Eigen::AngleAxisd(error_euler[1], Eigen::Vector3d::UnitY()) *
// 			Eigen::AngleAxisd(error_euler[0], Eigen::Vector3d::UnitX()));
//     
//     xki_k(6) = qe.x();
//     xki_k(7) = qe.y();
//     xki_k(8) = qe.z();
    
    /** Correct the attitude using the error quaternion **/
    q4 = q4 * qe;
    
    /** Normalize quaternion **/
    q4.normalize();
    
    /** Update the position error **/
    eposition = eposition + xki_k.block<NUMAXIS, 1> (0,0);
    
    /** Update the velocity error **/
    evelocity = evelocity + xki_k.block<NUMAXIS, 1> (NUMAXIS,0);
    
    double Tc = 1.0/5.0;
    filtered_evelocity.col(0) = filtered_evelocity.col(1) + (dt/Tc) * (evelocity - filtered_evelocity.col(1));    
    
    /** Slip detector **/
    slipdetector = Eigen::Matrix<double,NUMAXIS,1>::Ones()*std::numeric_limits<double>::quiet_NaN();
    double cov = this->Pki_k(3,3);
    
    if (fabs(filtered_evelocity.col(0)[0] - filtered_evelocity.col(1)[0]) > fabs(2.0*sqrt(cov)))
    {
	#ifdef DEBUG_PRINTS
	std::cout<< "[Update] Slip detected in X axis: "<<filtered_evelocity.col(0)[0] - filtered_evelocity.col(1)[0]<<" cov:"<<fabs(sqrt(cov))<<"\n";
	#endif
	slipdetector[0] = filtered_evelocity.col(0)[0];
	
    }
	
    cov = this->Pki_k(4,4);
    
    if (fabs(filtered_evelocity.col(0)[1] - filtered_evelocity.col(1)[1]) > fabs(2.0*sqrt(cov)))
    {
	#ifdef DEBUG_PRINTS
	std::cout<< "[Update] Slip detected in Y axis: "<<filtered_evelocity.col(0)[1] - filtered_evelocity.col(1)[1]<<" cov:"<<fabs(sqrt(cov))<<"\n";
	#endif
	slipdetector[1] = filtered_evelocity.col(0)[1];
    }
    
    /** Update the bias with the bias error **/
    bghat = bghat + xki_k.block<NUMAXIS, 1> (3*NUMAXIS,0);
    bahat = bahat + xki_k.block<NUMAXIS, 1> (4*NUMAXIS,0);
    
    #ifdef DEBUG_PRINTS
    Eigen::Matrix <double,NUMAXIS,1> euler; /** In euler angles **/
    euler[2] = qe.toRotationMatrix().eulerAngles(2,1,0)[0];//YAW
    euler[1] = qe.toRotationMatrix().eulerAngles(2,1,0)[1];//PITCH
    euler[0] = qe.toRotationMatrix().eulerAngles(2,1,0)[2];//ROLL
    std::cout<< "[Update] Error quaternion in euler\n";
    std::cout<< "[Update] Roll: "<<euler[0]*R2D<<" Pitch: "<<euler[1]*R2D<<" Yaw: "<<euler[2]*R2D<<"\n";
    std::cout<< "[Update] bghat is of size "<<bghat.rows()<<"x"<<bghat.cols()<<"\n";
    std::cout<< "[Update] bghat:\n"<<bghat<<"\n";
    std::cout<< "[Update] bahat is of size "<<bahat.rows()<<"x"<<bahat.cols()<<"\n";
    std::cout<< "[Update] bahat:\n"<<bahat<<"\n";
    #endif

    return;

}

void Sckf::motionModel(const Eigen::Matrix< double, Eigen::Dynamic, Eigen::Dynamic >& Anav, const Eigen::Matrix< double, Eigen::Dynamic, Eigen::Dynamic >& Bnav, const Eigen::Matrix< double, Eigen::Dynamic, 1  >& vjoints, double dt)
{
    
    Eigen::Matrix <double, Eigen::Dynamic, Eigen::Dynamic> R; /** Measurement Noise matrix of the observation vector of the LS **/
    Eigen::Matrix <double, Eigen::Dynamic, Eigen::Dynamic> W; /** Wheel-weighting matrix of the LS **/
    Eigen::Matrix <double, NUMAXIS+ENCODERS_VECTOR_SIZE, NUMAXIS+ENCODERS_VECTOR_SIZE> Rm1; /** Measurement Noise matrix of the measurement vector of the LS **/
    Eigen::Matrix <double, 2*NUMAXIS, 2*NUMAXIS> Ri; /** Measurement Noise matrix of the measurement vector of the LS **/
    Eigen::Matrix <double, NUMAXIS, NUMAXIS> Rv; /** Measurement Noise matrix of the linear velocity from IMU **/
    Eigen::Matrix< double, NUMAXIS , 1  > linvelo; /** Linear velocities from acc integration **/
    Eigen::Matrix< double, NUMAXIS , 1  > nadir; /** The nadir vector (gravity /norm(gravity) **/
    Eigen::Matrix< double, NUMAXIS , 1  > wheel_weighting; /** The nadir vector applied to teh actual attitude * q_weight distribution **/
    Eigen::Quaternion <double> attitude; /** Pitch and roll **/
    
    R.resize(NUMBER_OF_WHEELS*(2*NUMAXIS),NUMBER_OF_WHEELS*(2*NUMAXIS));
    W.resize(NUMBER_OF_WHEELS*(2*NUMAXIS),NUMBER_OF_WHEELS*(2*NUMAXIS));
    R.setZero(); Rm1.setZero(); W.setIdentity();
    
    /** Get the attitude with Yaw to Zero **/
    Eigen::Vector3d e = Eigen::Matrix3d(q4).eulerAngles(2,1,0);
    attitude = Eigen::Quaternion <double> (Eigen::AngleAxisd(0.00 * D2R, Eigen::Vector3d::UnitZ())*
    Eigen::AngleAxisd(e[1] * D2R, Eigen::Vector3d::UnitY()) *
    Eigen::AngleAxisd(e[0] * D2R, Eigen::Vector3d::UnitX()));
    
    nadir << 0.0, 0.0, 1.0; //!nadir in navigation frame
    
    /** Form the noise matrix for the velocity model (navigationKinematics) **/
    Rm1.block<NUMAXIS, NUMAXIS>(0,0) = measurement.getRangvelo();
    Rm1.block<ENCODERS_VECTOR_SIZE, ENCODERS_VECTOR_SIZE>(NUMAXIS,NUMAXIS) = measurement.getEncodersVelocityCovariance();
    
    Ri = Bnav.block<2*NUMAXIS, ENCODERS_VECTOR_SIZE>(0,NUMAXIS) * measurement.getEncodersVelocityCovariance() * Bnav.block<2*NUMAXIS, ENCODERS_VECTOR_SIZE>(0,NUMAXIS).transpose();
    
    for (unsigned int i=0; i<NUMBER_OF_WHEELS; ++i)
    {    
	R(i*(2*NUMAXIS),i*(2*NUMAXIS)) = Ri(0,0);
	R(i*(2*NUMAXIS)+1,i*(2*NUMAXIS)+1) = Ri(0,0);
	R(i*(2*NUMAXIS)+2,i*(2*NUMAXIS)+2) = Ri(0,0);
	R.block<NUMAXIS, NUMAXIS>(NUMAXIS+i*(2*NUMAXIS), NUMAXIS+i*(2*NUMAXIS)) = Rg;
    }
    
    /** Compute the determinant of the first block matrix of R **/
    double idet = R.block<2*NUMAXIS, 2*NUMAXIS>(0,0).inverse().determinant();
    
    /** Compute the wheel-weighting factor **/
    wheel_weighting = (attitude * measurement.getLevelWeightDistribution()) * nadir;
    if (wheel_weighting[0] < 0.00)
	wheel_weighting[0] = 1.0 + wheel_weighting[0];
    
    /** Wheel-weighting matrix **/
    for (unsigned int i=0; i<NUMBER_OF_WHEELS; ++i)
    {
	if ((i==0)||(i==1)) //!Rear Wheels
	    W.block<(2*NUMAXIS),(2*NUMAXIS)> (i*(2*NUMAXIS),i*(2*NUMAXIS)) = (idet*wheel_weighting[0]) * Eigen::Matrix<double, (2*NUMAXIS), (2*NUMAXIS)>::Identity();
	else //!Front Wheels
	    W.block<(2*NUMAXIS),(2*NUMAXIS)> (i*(2*NUMAXIS),i*(2*NUMAXIS)) = (idet*(1.0 - wheel_weighting[0])) * Eigen::Matrix<double, (2*NUMAXIS), (2*NUMAXIS)>::Identity();
    }

    #ifdef DEBUG_PRINTS
    typedef Eigen::Matrix <double, NUMBER_OF_WHEELS*(2*NUMAXIS),NUMBER_OF_WHEELS*(2*NUMAXIS)> matrixRType;
    Eigen::FullPivLU<matrixRType> lu_decompR(R);
    std::cout << "The rank of R is " << lu_decompR.rank() << std::endl;
    
    std::cout<< "[MotionModel] idet:\n"<<idet<<"\n";
    std::cout<< "[MotionModel] wheel_weighting:\n"<<wheel_weighting<<"\n";
    std::cout<< "[MotionModel] Rm1 is of size "<<Rm1.rows()<<"x"<<Rm1.cols()<<"\n";
    std::cout<< "[MotionModel] Rm1:\n"<<Rm1<<"\n";
    std::cout<< "[MotionModel] R is of size "<<R.rows()<<"x"<<R.cols()<<"\n";
    std::cout<< "[MotionModel] R:\n"<<R<<"\n";
    std::cout<< "[MotionModel] R.inverse():\n"<<R.inverse()<<"\n";
    std::cout<< "[MotionModel] W is of size "<<W.rows()<<"x"<<W.cols()<<"\n";
    std::cout<< "[MotionModel] W:\n"<<W<<"\n";
    std::cout<< "[MotionModel] W+R.inverse():\n"<<W*R.inverse()<<"\n";
    std::cout<< "[MotionModel] (W+R.inverse()).inverse():\n"<<(W*R.inverse()).inverse()<<"\n";
    #endif
    
    /** Set encoders velocity **/
    measurement.setEncodersVelocity(vjoints);
    
    /** Call the Navigation kinematics to know the velocity from odometry, slip vector and contact angles velocities **/
    measurement.navigationKinematics(Anav, Bnav, R, W);
    
    return;
}

void Sckf::velocityError(Eigen::Matrix< double, NUMAXIS, 1  >& velo_error, Eigen::Matrix< double, NUMAXIS , NUMAXIS  >& vel_errorCov, double dt)
{
    DataModel delta_VeloModel, delta_veloIMU;
    DataModel increVeloError;
       
    /** Increment in velocity from the Motion Model **/
    delta_VeloModel.data = measurement.getIncrementalVeloModel();
    delta_VeloModel.Cov = measurement.getIncrementalVeloModelCovariance();
    
    /** Increment in velocity from the IMU **/
    delta_veloIMU.data = measurement.getLinearVelocities(dt);
    delta_veloIMU.Cov = (this->Raup * dt);
    
    
    /** Check to not have negative uncertainty values in the covariance matrices **/
    for (register int i=0; i<delta_VeloModel.size();i++)
    {
	for (register int j=0; j<delta_VeloModel.size();j++)
	{
	    if (delta_VeloModel.Cov(i,j) < 0.00)
		delta_VeloModel.Cov(i,j) = 0.00;
	    if (delta_veloIMU.Cov(i,j) < 0.00)
		delta_veloIMU.Cov(i,j) = 0.00;
	}
	
    }
    
    /** Compute the error in the incremental velocity **/
    increVeloError = delta_veloIMU - delta_VeloModel;
    
    /** Compute the Bhattacharyya distance **/
    Eigen::Matrix<double, NUMAXIS, NUMAXIS> BC;
    
    /** There is only error in velocity if there is enough statistical evidence **/
    BC = measurement.bhattacharyya(delta_veloIMU, increVeloError);
    
    Hellinger = Eigen::Matrix<double, NUMAXIS, NUMAXIS>::Identity() - BC;
    
    /** Store the value in the argument variables **/
    velo_error = Hellinger * increVeloError.data;
    vel_errorCov =  increVeloError.Cov;
    
    velo_error[2] = 0.00; //!No error in Z axis
        
    #ifdef DEBUG_PRINTS
    std::cout<< "[Measurement] incre_vel_model:\n"<<delta_VeloModel.data<<"\n";
    std::cout<< "[Measurement] incre_vel_model_cov:\n"<<delta_VeloModel.Cov<<"\n";
    std::cout<< "[Measurement] incre_vel_imu:\n"<<delta_veloIMU.data<<"\n";
    std::cout<< "[Measurement] incre_vel_imu_cov:\n"<<delta_veloIMU.Cov<<"\n";
    std::cout<< "[Measurement] Squared Hellinger distance:\n"<<Hellinger<<"\n";
    std::cout<< "[Measurement] vel_error:\n"<<velo_error<<"\n";
    std::cout<< "[Measurement] vel_error_cov(Rme):\n"<<vel_errorCov<<"\n";
    #endif
    
    
    return;
}


void Sckf::measurementGeneration(const Eigen::Matrix< double, Eigen::Dynamic, Eigen::Dynamic >& Anav, const Eigen::Matrix< double, Eigen::Dynamic, Eigen::Dynamic >& Bnav,
				 const Eigen::Matrix< double, Eigen::Dynamic, 1  > &vjoints, Eigen::Matrix< double, NUMAXIS, 1> &velo_error,
				 Eigen::Matrix< double, NUMAXIS, NUMAXIS> &vel_errorCov, double dt)
{
    Eigen::Matrix <double, Eigen::Dynamic, Eigen::Dynamic> R; /** Measurement Noise matrix of the observation vector of the LS **/
    Eigen::Matrix <double, Eigen::Dynamic, Eigen::Dynamic> W; /** Wheel-weighting matrix of the LS **/
    Eigen::Matrix <double, NUMAXIS+ENCODERS_VECTOR_SIZE, NUMAXIS+ENCODERS_VECTOR_SIZE> Rm1; /** Measurement Noise matrix of the measurement vector of the LS **/
    Eigen::Matrix <double, 2*NUMAXIS, 2*NUMAXIS> Ri; /** Measurement Noise matrix of the measurement vector of the LS **/
    Eigen::Matrix <double, NUMAXIS, NUMAXIS> Rv; /** Measurement Noise matrix of the linear velocity from IMU **/
    Eigen::Matrix< double, NUMAXIS , 1  > linvelo; /** Linear velocities from acc integration **/
    Eigen::Matrix< double, NUMAXIS , 1  > nadir; /** The nadir vector (gravity /norm(gravity) **/
    Eigen::Matrix< double, NUMAXIS , 1  > wheel_weighting; /** The nadir vector applied to teh actual attitude * q_weight distribution **/
    Eigen::Quaternion <double> attitude; /** Pitch and roll **/
    
    R.resize(NUMBER_OF_WHEELS*(2*NUMAXIS),NUMBER_OF_WHEELS*(2*NUMAXIS));
    W.resize(NUMBER_OF_WHEELS*(2*NUMAXIS),NUMBER_OF_WHEELS*(2*NUMAXIS));
    R.setZero(); Rm1.setZero(); W.setIdentity();
    
    /** Get the attitude with Yaw to Zero **/
    Eigen::Vector3d e = Eigen::Matrix3d(q4).eulerAngles(2,1,0);
    attitude = Eigen::Quaternion <double> (Eigen::AngleAxisd(0.00 * D2R, Eigen::Vector3d::UnitZ())*
    Eigen::AngleAxisd(e[1] * D2R, Eigen::Vector3d::UnitY()) *
    Eigen::AngleAxisd(e[0] * D2R, Eigen::Vector3d::UnitX()));
    
    nadir << 0.0, 0.0, 1.0; //!nadir in navigation frame
    
    #ifdef DEBUG_PRINTS
    std::cout<<"********** MEASUREMENT GENERATION *****************\n";
    #endif
	
    /** Form the noise matrix for the velocity model (navigationKinematics) **/
    Rm1.block<NUMAXIS, NUMAXIS>(0,0) = this->Rg;
    Rm1.block<ENCODERS_VECTOR_SIZE, ENCODERS_VECTOR_SIZE>(NUMAXIS,NUMAXIS) = measurement.getEncodersVelocityCovariance();
    
    Ri = Bnav.block<2*NUMAXIS, ENCODERS_VECTOR_SIZE>(0,NUMAXIS) * measurement.getEncodersVelocityCovariance() * Bnav.block<2*NUMAXIS, ENCODERS_VECTOR_SIZE>(0,NUMAXIS).transpose();
    
    for (unsigned int i=0; i<NUMBER_OF_WHEELS; ++i)
    {    
	R(i*(2*NUMAXIS),i*(2*NUMAXIS)) = Ri(0,0);
	R(i*(2*NUMAXIS)+1,i*(2*NUMAXIS)+1) = Ri(0,0);
	R(i*(2*NUMAXIS)+2,i*(2*NUMAXIS)+2) = Ri(0,0);
	R.block<NUMAXIS, NUMAXIS>(NUMAXIS+i*(2*NUMAXIS), NUMAXIS+i*(2*NUMAXIS)) = Rg;
    }
    
    /** Compute the determinant of the first block matrix of R **/
    double idet = R.block<2*NUMAXIS, 2*NUMAXIS>(0,0).inverse().determinant();
    
    /** Compute the wheel-weithing factor **/
    wheel_weighting = (attitude * measurement.getLevelWeightDistribution()) * nadir;
    if (wheel_weighting[0] < 0.00)
	wheel_weighting[0] = 1.0 + wheel_weighting[0];
    
    /** Wheel-weighting matrix **/
    for (unsigned int i=0; i<NUMBER_OF_WHEELS; ++i)
    {
	if ((i==0)||(i==1)) //!Rear Wheels
	    W.block<(2*NUMAXIS),(2*NUMAXIS)> (i*(2*NUMAXIS),i*(2*NUMAXIS)) = (idet*wheel_weighting[0]) * Eigen::Matrix<double, (2*NUMAXIS), (2*NUMAXIS)>::Identity();
	else //!Front Wheels
	    W.block<(2*NUMAXIS),(2*NUMAXIS)> (i*(2*NUMAXIS),i*(2*NUMAXIS)) = (idet*(1.0 - wheel_weighting[0])) * Eigen::Matrix<double, (2*NUMAXIS), (2*NUMAXIS)>::Identity();
    }

    #ifdef DEBUG_PRINTS
    typedef Eigen::Matrix <double, NUMBER_OF_WHEELS*(2*NUMAXIS),NUMBER_OF_WHEELS*(2*NUMAXIS)> matrixRType;
    Eigen::FullPivLU<matrixRType> lu_decompR(R);
    std::cout << "The rank of R is " << lu_decompR.rank() << std::endl;
    
    std::cout<< "[Measurement] idet:\n"<<idet<<"\n";
    std::cout<< "[Measurement] wheel_weighting:\n"<<wheel_weighting<<"\n";
    std::cout<< "[Measurement] Rm1 is of size "<<Rm1.rows()<<"x"<<Rm1.cols()<<"\n";
    std::cout<< "[Measurement] Rm1:\n"<<Rm1<<"\n";
    std::cout<< "[Measurement] R is of size "<<R.rows()<<"x"<<R.cols()<<"\n";
    std::cout<< "[Measurement] R:\n"<<R<<"\n";
    std::cout<< "[Measurement] R.inverse():\n"<<R.inverse()<<"\n";
    std::cout<< "[Measurement] W is of size "<<W.rows()<<"x"<<W.cols()<<"\n";
    std::cout<< "[Measurement] W:\n"<<W<<"\n";
    std::cout<< "[Measurement] W+R.inverse():\n"<<W*R.inverse()<<"\n";
    std::cout<< "[Measurement] (W+R.inverse()).inverse():\n"<<(W*R.inverse()).inverse()<<"\n";
    #endif
    
    /** Set encoders velocity **/
    measurement.setEncodersVelocity(vjoints);
    
    /** Call the Navigation kinematics to know the velocity from odometry, slip vector and contact angles velocities **/
    measurement.navigationKinematics(Anav, Bnav, R, W);
    
    /** Integrate accelerometers to have velocity **/
//     linvelo = measurement.accIntegrationWindow(dt);
    
    /** Set the IMU velocity to the linear ones **/
//     measurement.setLinearVelocities(linvelo);
    
    /** Compute the velocity error **/
    DataModel veloModel;
    DataModel delta_VeloModel, delta_veloIMU;
    DataModel increVeloError;
    
    veloModel.data = measurement.getCurrentVeloModel();
    veloModel.Cov = measurement.getCurrentVeloModelCovariance();
    delta_VeloModel.data = measurement.getIncrementalVeloModel();
    delta_VeloModel.Cov = measurement.getIncrementalVeloModelCovariance();
    delta_veloIMU.data  = measurement.getLinearVelocities(dt);
    delta_veloIMU.Cov = (this->Raup * dt);
    
    for (register int i=0; i<veloModel.size();i++)
    {
	for (register int j=0; j<veloModel.size();j++)
	{
	    if (veloModel.Cov(i,j) < 0.00)
		veloModel.Cov(i,j) = 0.00;
	    if (delta_VeloModel.Cov(i,j) < 0.00)
		delta_VeloModel.Cov(i,j) = 0.00;
	    if (delta_veloIMU.Cov(i,j) < 0.00)
		delta_veloIMU.Cov(i,j) = 0.00;
	}
	
    }
    
    
    /** Compute the error in the incremental velocity **/
    increVeloError = delta_veloIMU - delta_VeloModel;
    
    /** Mahalanobis **/
    Mahalanobis = measurement.mahalanobis(increVeloError);
    
    /** Compute the Bhattacharyya distance **/
    Eigen::Matrix<double, NUMAXIS, NUMAXIS> BC;
    
    /*BC = measurement.bhattacharyya(delta_veloIMU, delta_VeloModel);
    std::cout<< "[Measurement] Bhattacharyya coeff IMU-Model(function):\n"<<M<<"\n";
    
    HellingerI_M = Eigen::Matrix<double, NUMAXIS, NUMAXIS>::Identity() - BC;*/
    
    /** DANGER **/
    BC = measurement.bhattacharyya(delta_veloIMU, increVeloError);
    
    Hellinger = Eigen::Matrix<double, NUMAXIS, NUMAXIS>::Identity() - BC;
    /** DANGER **/
    
    
//     /** Estimate if there is statictical significant difference for the incremental velocity error **/
//     Eigen::Matrix<double, Eigen::Dynamic, 1> eigen_values;
//     Eigen::JacobiSVD <Eigen::MatrixXd > svdOfCov(increVeloError.Cov, Eigen::ComputeThinU); //!SVD
//     
//     /** Eigen values are the axis of the ellipsoid **/
//     eigen_values.resize(increVeloError.data.rows());
//     eigen_values = svdOfCov.singularValues();
//     
//     for (register int i=0; i < eigen_values.rows();++i)
//     {
// 	/** If the 3-sigma is bigger that the difference, discard the measurement **/
// 	if (1.0 * sqrt(eigen_values[i]) > fabs(increVeloError.data[i]))
// 	{
// 	    velo_error[i] = 0.00;
// 	    vel_errorCov.row(i).setZero();
// 	    vel_errorCov.col(i).setZero();
// 	}
// 	else
// 	{
// 	    velo_error[i] = increVeloError.data[i];
// 	    vel_errorCov.row(i) = increVeloError.Cov.row(i);
// 	    vel_errorCov.col(i) = increVeloError.Cov.col(i);
// 	}
// 	    
//     }
    
    
//     Eigen::Matrix<double, NUMAXIS, 1> velrange, modelrange;
//     velrange = delta_veloIMU.Cov.diagonal().array().cwiseSqrt();
//     velrange[0] = velrange[0] * 2.0;
//     modelrange= delta_VeloModel.Cov.diagonal().array().cwiseSqrt();
//     
//     
//     #ifdef DEBUG_PRINTS
//     std::cout<< "[Measurement] velrange:\n"<<velrange<<"\n";
//     #endif
//     
//     for (register int i=0; i < NUMAXIS;++i)
//     {
// 	if (delta_veloIMU.data[i] <= velrange[i])
// 	{
// 	    velo_error[i] = 0.00;
// 	    vel_errorCov.row(i).setZero();
// 	    vel_errorCov.col(i).setZero();
// 	 
// // 	    if (delta_VeloModel.data[i] <= modelrange[i])
// // 	    {
// // 		/** Experimental **/
// // 		#ifdef DEBUG_PRINTS
// // 		std::cout<< "[Update Reset] evelocity_goes_zero ["<<i<<"]\n";
// // 		#endif
// // 		this->evelocity[i] = 0.00;
// // 		this->evelocity_cov.col(i).setZero();
// // 		this->evelocity_cov.row(i).setZero();
// // 	    }
// 	}
//     }
    
    /** DANGER **/
    
    velo_error = Hellinger * increVeloError.data;
    vel_errorCov =  increVeloError.Cov;
    
//     evelocity = evelocity + velo_error;
    
       
//     for (register int i=0; i < Hellinger.rows();++i)
//     {
// 	if (Hellinger(i,i) < 0.68)
// 	{
// 	    this->evelocity[i] = 0.00;
// 	    this->evelocity_cov.col(i).setZero();
// 	    this->evelocity_cov.row(i).setZero();
// 
// 	
// 	}
//     }
    /** DANGER **/
    
    
    velo_error[2] = 0.00;
    
    
    #ifdef DEBUG_PRINTS
    std::cout<< "[Measurement] incre_vel_imu:\n"<<delta_veloIMU.data<<"\n";
    std::cout<< "[Measurement] incre_vel_imu_cov:\n"<<delta_veloIMU.Cov<<"\n";
    std::cout<< "[Measurement] incre_vel_model:\n"<<delta_VeloModel.data<<"\n";
    std::cout<< "[Measurement] incre_vel_model_cov:\n"<<delta_VeloModel.Cov<<"\n";
    std::cout<< "[Measurement] vel_model:\n"<<veloModel.data<<"\n";
    std::cout<< "[Measurement] vel_model_cov:\n"<<veloModel.Cov<<"\n";
//     std::cout<< "[Measurement] sqrt(eigen_values):\n"<<eigen_values.cwiseSqrt()<<"\n";
    std::cout<< "[Measurement] Squared Hellinger distance:\n"<<Hellinger<<"\n";
    std::cout<< "[Measurement] vel_error:\n"<<velo_error<<"\n";
    std::cout<< "[Measurement] vel_error_cov(Rme):\n"<<vel_errorCov<<"\n";
    #endif
    
    return;
}

void Sckf::computeSlipVector(const Eigen::Matrix< double, Eigen::Dynamic, Eigen::Dynamic >& Aslip,
			     const Eigen::Matrix< double, Eigen::Dynamic, Eigen::Dynamic >& Bslip,
			     const Eigen::Matrix< double, NUMAXIS, 1 >& linvelo,
			     Eigen::Matrix< double, SLIP_VECTOR_SIZE, 1  >& slip_vector,
			     Eigen::Matrix< double, SLIP_VECTOR_SIZE, SLIP_VECTOR_SIZE >& slip_vectorCov, double dt)
{
    Eigen::Matrix <double, Eigen::Dynamic, Eigen::Dynamic> R; /** Measurement Noise matrix for the LS **/
    Eigen::Matrix <double, 2*NUMAXIS + ENCODERS_VECTOR_SIZE + NUMBER_OF_WHEELS, 2*NUMAXIS + ENCODERS_VECTOR_SIZE + NUMBER_OF_WHEELS> Robs; /** Measurement Noise matrix of the measurement vector of the LS **/
    R.resize(NUMBER_OF_WHEELS*(2*NUMAXIS),NUMBER_OF_WHEELS*(2*NUMAXIS));
    
//     /** Set the linear velocity to the corrected one **/
//     measurement.setLinearVelocities(this->velocity);
    
    /** Set the angular velocity to the corrected one **/
    Eigen::Matrix <double,NUMAXIS,1> angevelo;
    
    angevelo = (Eigen::Matrix3d(this->deltaQuaternion()).eulerAngles(2,1,0))/dt; //Yaw, Pitch and Roll
    
    /** Change the order Yaw is in vector[0] **/
    double aux = angevelo[0];
    angevelo[0] = angevelo[2];
    angevelo[2] = aux;
    
//     measurement.setAngularVelocities(angevelo); //Roll Pitch and Yaw
    
    /** Form the observation noise matrix **/
    Robs.setZero();
    Robs.block<NUMAXIS, NUMAXIS>(0,0) = measurement.getCurrentVeloModelCovariance() + evelocity_cov; //!Linear Velo
    Robs.block<NUMAXIS, NUMAXIS>(NUMAXIS,NUMAXIS) = Pki_k.block<NUMAXIS, NUMAXIS> (2*NUMAXIS,2*NUMAXIS);//!Angular velo
    Robs.block<ENCODERS_VECTOR_SIZE, ENCODERS_VECTOR_SIZE>(2*NUMAXIS,2*NUMAXIS) = measurement.getEncodersVelocityCovariance();
    Robs.block<NUMBER_OF_WHEELS, NUMBER_OF_WHEELS>(2*NUMAXIS + ENCODERS_VECTOR_SIZE,2*NUMAXIS + ENCODERS_VECTOR_SIZE) = measurement.getContactAnglesVelocityCovariance();
    
    /** Form the Noise matrix **/
    Eigen::Matrix <double, 2*NUMAXIS + ENCODERS_VECTOR_SIZE + NUMBER_OF_WHEELS, 1> robs;
    Eigen::Matrix <double, NUMBER_OF_WHEELS*(2*NUMAXIS), 1> rdiagonal;
    
    robs = Robs.diagonal();
    R = Bslip * robs.asDiagonal() * Bslip.transpose();
    rdiagonal = R.diagonal();
    R = rdiagonal.asDiagonal();
    
    #ifdef DEBUG_PRINTS
    typedef Eigen::Matrix <double, NUMBER_OF_WHEELS*(2*NUMAXIS),NUMBER_OF_WHEELS*(2*NUMAXIS)> matrixRType;
    Eigen::FullPivLU<matrixRType> lu_decompR(R);
    std::cout << "The rank of R is " << lu_decompR.rank() << std::endl;
    
    std::cout<< "[Measurement] Robs is of size "<<Robs.rows()<<"x"<<Robs.cols()<<"\n";
    std::cout<< "[Measurement] Robs:\n"<<Robs<<"\n";
    std::cout<< "[Measurement] R is of size "<<R.rows()<<"x"<<R.cols()<<"\n";
    std::cout<< "[Measurement] R:\n"<<R<<"\n";
    std::cout<< "[Measurement] R.inverse():\n"<<R.inverse()<<"\n";
    #endif
    
    /** Computes the slip kinematics **/
    measurement.slipKinematics(Aslip, Bslip, R, linvelo);
    
    /** Get the estimated slip error vector and the covariance **/
    slip_vector = measurement.getSlipVector();
    slip_vectorCov = measurement.getSlipVectorCov();
    
    return;

}



void Sckf::resetStateVector()
{

    /**------------------------------ **/
    /** Reset some parts of the state **/
    /**------------------------------ **/
    
    /** Reset the position part of the state **/
    xki_k.block<NUMAXIS,1>(0,0) = Eigen::Matrix<double, NUMAXIS, 1>::Zero();
    
    /** Reset the velocity part of the state **/
    xki_k.block<NUMAXIS,1>(NUMAXIS,0) = Eigen::Matrix<double, NUMAXIS, 1>::Zero();
    
    /** Reset the quaternion part of the state vector (the error quaternion) (q4 is corrected in the update step) **/
    xki_k.block<NUMAXIS,1>(2*NUMAXIS,0) = Eigen::Matrix<double, NUMAXIS, 1>::Zero();
    
    /** Reset the gyroscopes bias (there is a bghat variable) **/
    xki_k.block<NUMAXIS, 1> ((3*NUMAXIS),0) = Eigen::Matrix <double, NUMAXIS, 1>::Zero();
    
    /** Reset the accelerometers bias (there is a bahat variable) **/
    xki_k.block<NUMAXIS, 1> ((4*NUMAXIS),0) = Eigen::Matrix <double, NUMAXIS, 1>::Zero();
    
    #ifdef DEBUG_PRINTS
    std::cout<< "[Update After Reset] xki_k is of size "<<xki_k.rows()<<"x"<<xki_k.cols()<<"\n";
    std::cout<< "[Update After Reset] xki_k:\n"<<xki_k<<"\n";
    std::cout<< "[Update] bghat is of size "<<bghat.rows()<<"x"<<bghat.cols()<<"\n";
    std::cout<< "[Update] bghat:\n"<<bghat<<"\n";
    std::cout<< "[Update] bahat is of size "<<bahat.rows()<<"x"<<bahat.cols()<<"\n";
    std::cout<< "[Update] bahat:\n"<<bahat<<"\n";
    #endif
    
    /** Update the evelocity covariance matrix **/
    evelocity_cov += Pki_k.block<NUMAXIS, NUMAXIS> (NUMAXIS,NUMAXIS);
    
    return;
}


Eigen::Matrix <double,NUMAXIS,1> Sckf::getVeloError()
{
//     return evelocity;
    return filtered_evelocity.col(0);
}


Eigen::Matrix <double,NUMAXIS,1> Sckf::getVelocity()
{
    return evelocity;
}


/**
* @brief Save the current filter status
*/
void Sckf::toFilterInfo(localization::FilterInfo &finfo, double dt)
{
    
    finfo.xki_k = this->xki_k;
    finfo.Pki_k = this->Pki_k;
    finfo.K = this->K;
    finfo.Qk = this->Qk;
    finfo.Rk = this->Rk;
    finfo.Hk = this->Hk;
    finfo.zki = this->zki;
    finfo.innovation = this->innovation;
    finfo.eposition = this->eposition;
    finfo.evelocity = this->evelocity;
    finfo.evelocity_cov = this->evelocity_cov;
    finfo.bghat = this->bghat;
    finfo.bahat = this->bahat;
    Eigen::VectorXd aux(2);
    aux = filtered_evelocity.row(0);
    finfo.eacceleration[0] = MeasurementItem::finiteDifference(aux, dt);
    aux = filtered_evelocity.row(1);
    finfo.eacceleration[1] = MeasurementItem::finiteDifference(aux, dt);
    aux = filtered_evelocity.row(2);
    finfo.eacceleration[2] = MeasurementItem::finiteDifference(aux, dt);
    finfo.Hellinger = this->Hellinger;
    finfo.Mahalanobis = this->Mahalanobis;
    finfo.slipdetector = this->slipdetector;
    
    /** For next step **/
    filtered_evelocity.col(1) = filtered_evelocity.col(0);
    
    return;
}

