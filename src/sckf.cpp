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
    
    //std::cout << "Attitude (getEuler): "<< euler(0)<<" "<<euler(1)<<" "<<euler(2)<<"\n";
    //std::cout << "Attitude in degrees (getEuler): "<< euler(0)*R2D<<" "<<euler(1)*R2D<<" "<<euler(2)*R2D<<"\n";
    
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
* @brief This function Initialize Attitude
*/
int sckf::setAttitude(Eigen::Quaternion< double >& initq)
{

    
    /** Initial orientation **/
    q4 = initq;
	
    return OK;
    
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
    x_0.resize(sckf::XSTATEVECTORSIZE,1);
    this->xki_k = x_0;
}

void sckf::setEccentricity(Eigen::Matrix <double,NUMAXIS,1>  &eccx, Eigen::Matrix <double,NUMAXIS,1>  &eccy, Eigen::Matrix <double,NUMAXIS,1>  &eccz)
{
    this->eccx = eccx;
    this->eccy = eccy;
    this->eccz = eccz;
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
    xk_k.resize(sckf::XSTATEVECTORSIZE,1);
    xki_k.resize(sckf::XSTATEVECTORSIZE,1);
    
    A.resize(sckf::ASTATEVECTORSIZE,sckf::ASTATEVECTORSIZE);
    Fki.resize(sckf::XSTATEVECTORSIZE,sckf::XSTATEVECTORSIZE);
    
    Qk.resize(sckf::XSTATEVECTORSIZE,sckf::XSTATEVECTORSIZE);
    
    Pk_k.resize(sckf::XSTATEVECTORSIZE,sckf::XSTATEVECTORSIZE);
    Pki_k.resize(sckf::XSTATEVECTORSIZE,sckf::XSTATEVECTORSIZE);
    
    H1a.resize(NUMAXIS, sckf::ASTATEVECTORSIZE);
    H2a.resize(NUMAXIS, sckf::ASTATEVECTORSIZE);
    Hk.resize(sckf::ZMEASUREMENTVECTORSIZE, sckf::XSTATEVECTORSIZE);
    
    zki.resize(ZMEASUREMENTVECTORSIZE, 1);
    
    Rk.resize(sckf::ZMEASUREMENTVECTORSIZE, sckf::ZMEASUREMENTVECTORSIZE);
    
    K.resize(sckf::XSTATEVECTORSIZE, sckf::ZMEASUREMENTVECTORSIZE);
    
    /** Resizing dynamic arguments to the correct dimension to avoid matrices errors **/
    P_0.resize(sckf::XSTATEVECTORSIZE, sckf::XSTATEVECTORSIZE);
    Qec.resize((sckf::ESTATEVECTORSIZE*sckf::NUMBEROFWHEELS), (sckf::ESTATEVECTORSIZE*sckf::NUMBEROFWHEELS));
    Ren.resize(sckf::EMEASUREMENTVECTORSIZE-(2*NUMAXIS), sckf::EMEASUREMENTVECTORSIZE-(2*NUMAXIS));
    
    /** Gravitation acceleration **/
    gtilde << 0, 0, g;

    /** Dip angle (alpha is in rad) **/
    mtilde(0) = cos(alpha);
    mtilde(1) = 0;
    mtilde(2) = -sin(alpha);

    
    /** Kalman filter state, error covariance and process noise covariance **/
    xki_k = Matrix <double,sckf::XSTATEVECTORSIZE,1>::Zero();
    xk_k = xki_k;
    
    Qk = Matrix <double,sckf::XSTATEVECTORSIZE,sckf::XSTATEVECTORSIZE>::Zero();
    Qk.block <(sckf::ESTATEVECTORSIZE*sckf::NUMBEROFWHEELS), (sckf::ESTATEVECTORSIZE*sckf::NUMBEROFWHEELS)> (0,0) = Qec;
    Qk.block <NUMAXIS, NUMAXIS> ((sckf::ESTATEVECTORSIZE*sckf::NUMBEROFWHEELS),(sckf::ESTATEVECTORSIZE*sckf::NUMBEROFWHEELS)) = 0.25 * Rg;
    Qk.block <NUMAXIS, NUMAXIS> ((sckf::ESTATEVECTORSIZE*sckf::NUMBEROFWHEELS)+NUMAXIS,(sckf::ESTATEVECTORSIZE*sckf::NUMBEROFWHEELS)+NUMAXIS) = Qbg;
    Qk.block <NUMAXIS, NUMAXIS> ((sckf::ESTATEVECTORSIZE*sckf::NUMBEROFWHEELS)+(2*NUMAXIS),(sckf::ESTATEVECTORSIZE*sckf::NUMBEROFWHEELS)+(2*NUMAXIS)) = Qba;

    /** Assign the initial values **/
    Pki_k = P_0;
    
    /** Asign the initial value for the measurement matrix **/
    Hk = Matrix <double,sckf::ZMEASUREMENTVECTORSIZE, sckf::XSTATEVECTORSIZE>::Zero();
    H1a = Matrix <double,NUMAXIS,sckf::ASTATEVECTORSIZE>::Zero();
    H2a = Matrix <double,NUMAXIS,sckf::ASTATEVECTORSIZE>::Zero();
    H1a(0,6) = 1; H1a(1,7) = 1; H1a(2,8) = 1;
    
    
    /** System matrix A **/
    A = Matrix <double,sckf::ASTATEVECTORSIZE, sckf::ASTATEVECTORSIZE>::Zero();      
    A(0,3) = -0.5;A(1,4) = -0.5;A(2,5) = -0.5;
    
    /** Initial measurement noise **/
    Rk = Matrix <double,sckf::ZMEASUREMENTVECTORSIZE,sckf::ZMEASUREMENTVECTORSIZE>::Zero();
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
    
	
    /** Variable in the adaptive algorithm **/
    r1count = 0;
    r2count = R2COUNT;
    
    /** Print filter information **/
    std::cout<< "Pk+i|k:\n"<<Pki_k<<"\n";
    std::cout<< "Qk|k:\n"<<Qk<<"\n";
    std::cout<< "Rk:\n"<<Rk<<"\n";
    std::cout<< "H1a:\n"<<H1a<<"\n";
    std::cout<< "H2a:\n"<<H2a<<"\n";
    std::cout<< "A:\n"<<A<<"\n";
    std::cout<< "mtilde:\n"<<mtilde<<"\n";
    std::cout<< "gtilde:\n"<<gtilde<<"\n";
    std::cout<< "Rv:\n"<<Rv<<"\n";
    std::cout<< "Rg:\n"<<Rg<<"\n";
    std::cout<< "Ra:\n"<<Ra<<"\n";
    std::cout<< "Ren:\n"<<Ren<<"\n";
    std::cout<< "Rm:\n"<<Rm<<"\n";

}

/**
* @brief Performs the prediction step of the filter.
*/
void sckf::predict(Eigen::Matrix< double, 3 , 1  >& u, double dt)
{
    
    Eigen::Matrix <double,NUMAXIS,NUMAXIS> vec2product; /** Vec 2 product  matrix */
    Eigen::Matrix <double,NUMAXIS,1> angvelo; /** Angular velocity */
    Eigen::Matrix <double,QUATERSIZE,QUATERSIZE> omega4; /** Quaternion integration matrix */
    Eigen::Matrix <double,QUATERSIZE,1> quat; /** Quaternion integration matrix */
    Eigen::Matrix <double,ESTATEVECTORSIZE, ESTATEVECTORSIZE> Fe; /** System matrix of a single wheel position error */
    Eigen::Matrix <double,ASTATEVECTORSIZE,ASTATEVECTORSIZE> dA; /** Discrete System matrix */
    Eigen::Matrix <double,XSTATEVECTORSIZE,XSTATEVECTORSIZE> Qdk; /** Discrete Qk matrix */

    /** Compute the vector2product matrix with the angular velocity **/
    angvelo = u - bghat; /** Eliminate the Bias **/

    vec2product << 0, -angvelo(2), angvelo(1),
		angvelo(2), 0, -angvelo(0),
		-angvelo(1), angvelo(0), 0;
		
    /** Compute the dA Matrix of the attitude part **/
    A.block<NUMAXIS, NUMAXIS> (0,0) = -vec2product;
    dA = Eigen::Matrix<double,ASTATEVECTORSIZE,ASTATEVECTORSIZE>::Identity() + A * dt + A * A * pow(dt,2)/2;
    
    /** Form a single wheel position error and contact angle **/
    Fe = Eigen::Matrix <double,ESTATEVECTORSIZE,ESTATEVECTORSIZE>::Zero();
    Fe(ESTATEVECTORSIZE-1, ESTATEVECTORSIZE-1)  = 1.0; /** The contact angle is modelled as a Wiener process (gaussian processes) **/
   
    /** Form the complete system model matrix **/
    for (int i = 0; i < NUMBEROFWHEELS; i++)
	Fki.block<ESTATEVECTORSIZE, ESTATEVECTORSIZE> (i*ESTATEVECTORSIZE,i*ESTATEVECTORSIZE) = Fe;

    Fki.block<ASTATEVECTORSIZE, ASTATEVECTORSIZE>(NUMBEROFWHEELS*ESTATEVECTORSIZE, NUMBEROFWHEELS*ESTATEVECTORSIZE) = dA;
    
    /** Propagate the vector through the system **/
    xki_k = Fki * xki_k;
    
    /** Form the system noise matrix for the attitude **/
    Qdk = Qk*dt + 0.5*Fki*Qk + 0.5*Fki*A.transpose();
    Qdk = 0.5*(Qdk + Qdk.transpose());
    Pki_k = Fki*Pki_k*Fki.transpose() + Qdk;
        
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

void sckf::update(Eigen::Matrix <double,Eigen::Dynamic,Eigen::Dynamic> &He, Eigen::Matrix <double,Eigen::Dynamic,Eigen::Dynamic> &Be,
		  Eigen::Matrix< double, Eigen::Dynamic, 1  >& encoders, Eigen::Matrix< double, 3 , 1  >& acc,
		  Eigen::Matrix< double, 3 , 1  >& gyro, Eigen::Matrix< double, 3 , 1  >& mag, double dt, bool magn_on_off)
{
    
    register int j;
    Eigen::Matrix <double,NUMAXIS,NUMAXIS> Cq; /** Rotational matrix */
    Eigen::Matrix <double,NUMAXIS,NUMAXIS> gtilde2product; /** Vec 2 product  matrix for the gravioty vector in body frame*/
    Eigen::Matrix <double,NUMAXIS,NUMAXIS> gyros2product; /** Vec 2 product  matrix for the gyroscopes (angular velocity) */
    Eigen::Matrix <double,NUMAXIS,1> angvelo; /** Angular velocity */
    Eigen::Matrix <double,NUMAXIS,NUMAXIS> fooR2; /**  Measurement noise matrix from accelerometers matrix Ra*/
    Eigen::Matrix <double,ASTATEVECTORSIZE,1> xa_k; /** Attitude part of the state vector xk+i|k */
    Eigen::Matrix <double,ASTATEVECTORSIZE,ASTATEVECTORSIZE> P1a; /** Error convariance matrix for measurement 1 of the attitude */
    Eigen::Matrix <double,ASTATEVECTORSIZE,ASTATEVECTORSIZE> P2a; /** Error convariance matrix for measurement 2 of the attitude */
    Eigen::Matrix <double,ASTATEVECTORSIZE,ASTATEVECTORSIZE> auxM; /** Auxiliar matrix for computing Kalman gain in measurement */
    Eigen::Matrix <double,ASTATEVECTORSIZE, NUMAXIS> K1a; /** Kalman Gain matrix for measurement 1 */
    Eigen::Matrix <double,ASTATEVECTORSIZE, NUMAXIS> K2a; /** Kalman Gain matrix for measurement 2 */
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
    Eigen::Matrix <double,EMEASUREMENTVECTORSIZE,1> ye; /** Measurement vectors for the rover position error ze = Be*ye */
    Eigen::Matrix <double,NUMAXIS,1> auxvector; /** Auxiliar vector variable */
    
    
    /** First measurement step (Pitch and Roll correction from Acc) **/
    
    /** Copy the attitude part of the state vector and covariance matrix **/
    xa_k = xki_k.block<ASTATEVECTORSIZE, 1> ((ESTATEVECTORSIZE*NUMBEROFWHEELS), 0);
    P1a = Pki_k.block<ASTATEVECTORSIZE, ASTATEVECTORSIZE> ((ESTATEVECTORSIZE*NUMBEROFWHEELS), (ESTATEVECTORSIZE*NUMBEROFWHEELS));
    
    /** Calculate the gravity vector in the body frame **/
    gtilde_body = q4 * gtilde;
    gtilde2product << 0, -gtilde_body(2), gtilde_body(1),
		gtilde_body(2), 0, -gtilde_body(0),
		-gtilde_body(1), gtilde_body(0), 0;
	
    /** Eliminate the Bias from gyros**/
    angvelo = gyro - bghat; 

    /** Compute the vector2product matrix with the angular velocity **/
    /** in order to In order to remove the centripetal velocities **/
    gyros2product << 0, -angvelo(2), angvelo(1),
		angvelo(2), 0, -angvelo(0),
		-angvelo(1), angvelo(0), 0;
		
    /** Form the observation matrix Hk **/
    Hk.block<EMEASUREMENTVECTORSIZE, XSTATEVECTORSIZE> (0,0) = He;
    
    /** Form the matrix for the measurement 1 of the attitude (acc correction) **/
    H1a.block<NUMAXIS, NUMAXIS> (0,0) = 2*gtilde2product;
    
    /** Copy to the whole observation matrix **/
    Hk.block<NUMAXIS, ASTATEVECTORSIZE> (EMEASUREMENTVECTORSIZE, XSTATEVECTORSIZE) = H1a;
    
    /** Form the measurement vector z1a for the attitude **/
    z1a = acc - bahat - gtilde_body;
    
    /** Form the measurement vector ye of the rover position error (Be*ye) **/
    ye(0,0) = (z1a[0] * dt) - (gyros2product.row(0) * eccx);
    ye(1,0) = (z1a[1] * dt) - (gyros2product.row(1) * eccy);
    ye(2,0) = (z1a[2] * dt) - (gyros2product.row(2) * eccz);
    ye.block<NUMAXIS, 1> (NUMAXIS, 0) = angvelo;
    ye.block<EMEASUREMENTVECTORSIZE-(2*NUMAXIS), 1> ((2*NUMAXIS), 0) = encoders;
    
    /** Form the complete zk vector **/
    zki.block<EMEASUREMENTVECTORSIZE, 1> (0,0)= Be*ye;
    zki.block<NUMAXIS, 1> (EMEASUREMENTVECTORSIZE, 0)= z1a;
   
    
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
    
    Uk = Uk / (M1);
    
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
	r2count = 0;
	auxvector(0) = std::max(lambda(0)-mu(0),(double)0.00);
	auxvector(1) = std::max(lambda(1)-mu(1),(double)0.00);
	auxvector(2) = std::max(lambda(2)-mu(2),(double)0.00);
	
	Qstar = auxvector(0) * u.col(0) * u.col(0).transpose() + auxvector(1) * u.col(1) * u.col(1).transpose() + auxvector(2) * u.col(2) * u.col(2).transpose();
    }
    else
    {
	r2count ++;
	if (r2count < M2)
	    Qstar = auxvector(0) * u.col(0) * u.col(0).transpose() + auxvector(1) * u.col(1) * u.col(1).transpose() + auxvector(2) * u.col(2) * u.col(2).transpose();
	else
	    Qstar = Matrix<double, NUMAXIS, NUMAXIS>::Zero();
    }
    
    /** Form the Rk matrix **/
    Rk.block<NUMAXIS, NUMAXIS>(0,0) = Rv; /** For the linear velocity **/
    Rk.block<NUMAXIS, NUMAXIS>(NUMAXIS,NUMAXIS) = Rg; /** For the angular velocity **/
    Rk.block<EMEASUREMENTVECTORSIZE-(2*NUMAXIS), EMEASUREMENTVECTORSIZE-(2*NUMAXIS)>((2*NUMAXIS),(2*NUMAXIS)) = Ren; /** For the encoders **/
    Rk.block<NUMAXIS, NUMAXIS>(EMEASUREMENTVECTORSIZE,EMEASUREMENTVECTORSIZE) = Ra + Qstar; /** For the attitude correction **/
    
    /** Compute the Kalman Gain Matrix **/
    K = Pki_k * Hk.transpose() * (Hk * Pki_k * Hk.transpose() + Rk).inverse();
    
    /** Update the state vector and the covariance matrix **/
    xki_k = xki_k + K * (zki - Hk*xki_k);
        
    Pki_k = (Matrix<double,XSTATEVECTORSIZE,XSTATEVECTORSIZE>::Identity()-K*Hk)*Pki_k*(Matrix<double,XSTATEVECTORSIZE,XSTATEVECTORSIZE>::Identity()-K*Hk).transpose() + K*Rk*K.transpose();
    Pki_k = 0.5 * (Pki_k + Pki_k.transpose());
         
    /** Update the quaternion with the Indirect approach **/
    qe.w() = 1;
    qe.x() = xki_k((ESTATEVECTORSIZE*NUMBEROFWHEELS));
    qe.y() = xki_k((ESTATEVECTORSIZE*NUMBEROFWHEELS)+1);
    qe.z() = xki_k((ESTATEVECTORSIZE*NUMBEROFWHEELS)+2);
    q4 = q4 * qe;
    
    /** Normalize quaternion **/
    q4.normalize();

    
    /** Reset the quaternion part of the state vector **/
    xki_k.block<NUMAXIS,1>((ESTATEVECTORSIZE*NUMBEROFWHEELS),0) = Matrix<double, NUMAXIS, 1>::Zero();
    
    /**---------------------------- **/
    /** Reset the rest of the state **/
    /**---------------------------- **/
    bghat = bghat + xki_k.block<NUMAXIS, 1> ((ESTATEVECTORSIZE*NUMBEROFWHEELS) + NUMAXIS,0);
    xki_k.block<NUMAXIS, 1> ((ESTATEVECTORSIZE*NUMBEROFWHEELS) + NUMAXIS,0) = Matrix <double, NUMAXIS, 1>::Zero();
    
    bahat = bahat + xki_k.block<NUMAXIS, 1> ((ESTATEVECTORSIZE*NUMBEROFWHEELS) + (2*NUMAXIS),0);
    xki_k.block<NUMAXIS, 1> ((ESTATEVECTORSIZE*NUMBEROFWHEELS) + (2*NUMAXIS),0) = Matrix <double, NUMAXIS, 1>::Zero();
    
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
    
//     std::cout<<"Earth Rotation:"<<v<<"\n";

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
	std::cout << "[EAST] magnetic declination\n";
	euler[2] -= magnetic_declination; /** Magnetic declination is positive **/
    }
    else if (mode == WEST)
    {
	std::cout << "[WEST] magnetic declination\n";
	euler[2] += magnetic_declination; /** Magnetic declination is negative **/
    }
    else
    {
	std::cerr << "[ERROR] In the correction of the magnetic declination\n";
	return ERROR;
    }
    
    *quat = Eigen::Quaternion <double> (Eigen::AngleAxisd(euler[0], Eigen::Vector3d::UnitX())*
			Eigen::AngleAxisd(euler[1], Eigen::Vector3d::UnitY()) *
			Eigen::AngleAxisd(euler[2], Eigen::Vector3d::UnitZ()));
    
    return OK;
}






