/**\file measurement.cpp
 *
 * This class has the primitive methods for the Measurement Generation of the localization framework.
 * The class perform proprioceptive measurements according to the kinematics jacobian
 * The Jacobian should externally being provided by the Navigation Kinematics implementation
 * which is particular of each mobile robot chassis.
 * 
 * Typical proprioceptive sensors are inertial sensors and encoders. The main output of this class
 * if the rover velocity, contact angle of the point in contact and slip vectors of each wheels.
 * 
 * @author Javier Hidalgo Carrio | DFKI RIC Bremen | javier.hidalgo_carrio@dfki.de
 * @date November 2012.
 * @version 1.0.
 */

#include "measurement.hpp"

#define DEBUG_PRINTS 1

using namespace localization;
using namespace Eigen;

void measurement::welcome()
{
	std::cout << "You successfully compiled and executed Measurement Generation. Welcome!" << std::endl;
}

/**
* @brief Set the current linear velocities
*/
void measurement::setLinearVelocities(Eigen::Matrix< double, 3 , 1  > linvelo)
{
    this->linvelocity = linvelo;
    
    return;
}

/**
* @brief Return linear velocities
*/
Eigen::Matrix< double, NUMAXIS , 1  > measurement::getLinearVelocities()
{
    return this->linvelocity;
}

/**
* @brief Set the current angular velocity
*/
void measurement::setAngularVelocities(Eigen::Matrix< double, NUMAXIS , 1  > angvelo)
{
    this->angvelocity = angvelo;
    
    return;
}

/**
* @brief Return angular velocities
*/
Eigen::Matrix< double, NUMAXIS , 1  > measurement::getAngularVelocities()
{
    return this->angvelocity;
}

/**
* @brief Set the linear Acceleration
*/
void measurement::setLinearAcceleration(Eigen::Matrix< double, NUMAXIS , 1  > linacc)
{
    this->linacceleration.col(0) = linacc;

    return;
}


/**
* @brief Return linear acceleration
*/
Eigen::Matrix< double, NUMAXIS , 1  > measurement::getLinearAcceleration()
{
    return this->linacceleration.col(0);
}


/**
* @brief Return slip vector for wheel_idx
*/
Eigen::Matrix< double, NUMAXIS , 1  > measurement::getSlipVector(int wheel_idx)
{
    if (wheel_idx < NUMBER_OF_WHEELS)
	return this->slipMatrix.col(wheel_idx);
    
    return Eigen::Matrix <double, NUMAXIS, 1>::Zero();
}

/**
* @brief Return the contact angles
*/

Eigen::Matrix< double, Eigen::Dynamic, 1  > measurement::getContactAnglesVelocity()
{
    return this->acontact;
}

/**
* @brief Set the current velocity model
*/
void measurement::setCurrentVeloModel(Eigen::Matrix< double, NUMAXIS , 1  > velocity)
{
    cbVelModelX.push_back(velocity(0));
    cbVelModelY.push_back(velocity(1));
    cbVelModelZ.push_back(velocity(2));
}

/**
* @brief Return the current velocity from odometry model
*/
Eigen::Matrix< double, NUMAXIS , 1  > measurement::getCurrentVeloModel()
{
    Eigen::Matrix< double, NUMAXIS , 1  > velocity;
    
    velocity[0] = cbVelModelX[cbVelModelX.size()-1];
    velocity[1] = cbVelModelY[cbVelModelY.size()-1];
    velocity[2] = cbVelModelZ[cbVelModelZ.size()-1];
    
    return velocity;
}

/**
* @brief This function set the Accelerometers excentricity
*/
void measurement::setEccentricity(Eigen::Matrix <double,NUMAXIS,1>  &eccx, Eigen::Matrix <double,NUMAXIS,1>  &eccy, Eigen::Matrix <double,NUMAXIS,1>  &eccz)
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
* @brief Initialization method
*/
void measurement::Init(Eigen::Matrix< double, NUMAXIS , NUMAXIS  >& Ra,
		       Eigen::Matrix< double, NUMAXIS , NUMAXIS  >& Rg,
		       Eigen::Matrix< double, ENCODERS_VECTOR_SIZE , ENCODERS_VECTOR_SIZE  >& Ren,
		       Eigen::Matrix< double, NUMBER_OF_WHEELS , NUMBER_OF_WHEELS  >& Rcont)
{
    /** Initial encoders velocities **/
    encodersvelocity.setZero();
    
    /** Initial contact angle **/
    acontact.setZero();
    
    /** Slip velocities structures **/
    slipModel.resize(SLIP_VECTOR_SIZE);
    slipInertial.resize(SLIP_VECTOR_SIZE);
    slipError.resize(SLIP_VECTOR_SIZE);
    
    /** Initial slip matrix **/
    slipMatrix = Eigen::Matrix <double, NUMAXIS, NUMBER_OF_WHEELS>::Zero();
    
    /** Velocities for the filter (acce and bias substracted correctly) **/
    linvelocity = Eigen::Matrix <double, NUMAXIS, 1>::Zero();
    angvelocity = Eigen::Matrix <double, NUMAXIS, 1>::Zero();
    linacceleration = Eigen::Matrix <double, NUMAXIS, 3>::Zero();
    
    /** Circular vector for the integration **/
    cbAccX.set_capacity(INTEGRATION_XAXIS_WINDOW_SIZE);
    cbAccY.set_capacity(INTEGRATION_YAXIS_WINDOW_SIZE);
    cbAccZ.set_capacity(INTEGRATION_ZAXIS_WINDOW_SIZE);
    
    cbAngveloX.set_capacity(ANGVELO_WINDOW_SIZE);
    cbAngveloY.set_capacity(ANGVELO_WINDOW_SIZE);
    cbAngveloZ.set_capacity(ANGVELO_WINDOW_SIZE);
    
    /** Resize the circular buffer for the velocities model **/
    cbVelModelX.set_capacity(INTEGRATION_XAXIS_WINDOW_SIZE);
    cbVelModelY.set_capacity(INTEGRATION_YAXIS_WINDOW_SIZE);
    cbVelModelZ.set_capacity(INTEGRATION_ZAXIS_WINDOW_SIZE);

    this->Ra = Ra;
    this->Rg = Rg;
    this->Rencoders = Ren;
    this->Rcontact = Rcont;
    
    /** Print filter information **/
    #ifdef DEBUG_PRINTS
    std::cout<< "Ra is of size "<<Ra.rows()<<"x"<<Ra.cols()<<"\n";
    std::cout<< "Ra:\n"<<Ra<<"\n";
    std::cout<< "Rg is of size "<<Rg.rows()<<"x"<<Rg.cols()<<"\n";
    std::cout<< "Rg:\n"<<Rg<<"\n";
    std::cout<< "Ren is of size "<<Rencoders.rows()<<"x"<<Rencoders.cols()<<"\n";
    std::cout<< "Ren:\n"<<Rencoders<<"\n";
    std::cout<< "Rcont is of size "<<Rcontact.rows()<<"x"<<Rcontact.cols()<<"\n";
    std::cout<< "Rcontact:\n"<<Rcontact<<"\n";
    std::cout<< "acontact is of size "<<acontact.rows()<<"x"<<acontact.cols()<<"\n";
    std::cout<< "acontact:\n"<<acontact<<"\n";
    #endif
}


/**
* @brief Set the current contact angles
*/
void measurement::setContactAnglesVelocity(Eigen::Matrix<double, localization::NUMBER_OF_WHEELS, 1> contact_angles)
{
    this->acontact = contact_angles;
}



/**
* @brief This perform simpson rule for integration
*/
double measurement::simpsonsIntegral(double fa, double fm, double fb, double delta_ab)
{
    return delta_ab/6.0 * (fa + (4.0*fm) + fb);
}

/**
* @brief To obtain the linear velocity 
*/
Eigen::Matrix< double, NUMAXIS , 1  > measurement::accIntegrationWindow(double dt)
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
    
    localvelocity[0] += cbVelModelX[0];
    
    for (unsigned int i=0; i<cbAccY.size(); i++)
    {
	gyros2product << 0, -cbAngveloZ[i], cbAngveloY[i],
		cbAngveloZ[i], 0, -cbAngveloX[i],
		-cbAngveloY[i], cbAngveloX[i], 0;
		
	localvelocity[1] += (cbAccY[i] * dt) - (gyros2product.row(1) * eccy);
    }
    
    localvelocity[1] += cbVelModelY[0];
    
    for (unsigned int i=0; i<cbAccZ.size(); i++)
    {
	gyros2product << 0, -cbAngveloZ[i], cbAngveloY[i],
		cbAngveloZ[i], 0, -cbAngveloX[i],
		-cbAngveloY[i], cbAngveloX[i], 0;
		
	localvelocity[2] += (cbAccZ[i] * dt) - (gyros2product.row(2) * eccz);
    }
    
    localvelocity[2] += cbVelModelZ[0];
    
    return localvelocity;
}

void measurement::navigationKinematics(const Eigen::Matrix< double, Eigen::Dynamic, Eigen::Dynamic > &A, const Eigen::Matrix< double, Eigen::Dynamic, Eigen::Dynamic > &B,
					const Eigen::Matrix< double, Eigen::Dynamic, Eigen::Dynamic > &R)
{
    double leastSquaresError;
    Eigen::Matrix< double, NUMAXIS , 1  > velocity;
    Eigen::Matrix <double, SLIP_VECTOR_SIZE, 1> slip;
    
    /** Navigation Vectors **/
    Eigen::Matrix< double, Eigen::Dynamic, Eigen::Dynamic > x;
    Eigen::Matrix< double, Eigen::Dynamic, Eigen::Dynamic > y;
    Eigen::Matrix< double, Eigen::Dynamic, Eigen::Dynamic > b;
    
    /** Resize vector to the correct size **/
    x.resize(NUMAXIS+(2*NUMBER_OF_WHEELS), 1);
    y.resize(NUMAXIS+ENCODERS_VECTOR_SIZE, 1);
    b.resize(NUMBER_OF_WHEELS*(2*NUMAXIS), 1);
    
    /** Form the sensed vector **/
    y.block<NUMAXIS, 1> (0,0) = angvelocity;
    y.block<ENCODERS_VECTOR_SIZE, 1> (NUMAXIS,0) = encodersvelocity;
    
    #ifdef DEBUG_PRINTS
    std::cout<<"[VM] R is of size "<<R.rows()<<"x"<<R.cols()<<"\n";
    std::cout<<"[VM] R:\n" << R << std::endl;
    std::cout<<"[VM] A is of size "<<A.rows()<<"x"<<A.cols()<<"\n";
    std::cout<<"[VM] A:\n" << A << std::endl;
    std::cout<<"[VM] B is of size "<<B.rows()<<"x"<<B.cols()<<"\n";
    std::cout<<"[VM] B:\n" << B << std::endl;
    std::cout<<"[VM] y is of size " <<y.rows()<<"x"<< y.cols()<<"\n";
    std::cout<<"[VM] The sensed vector y\n"<<y<<"\n";
    #endif
    
    b = B*y;
    
    /** DEBUG OUTPUT **/
    #ifdef DEBUG_PRINTS
    typedef Eigen::Matrix <double, NUMBER_OF_WHEELS*(2*NUMAXIS), 11> matrixAType;
    typedef Eigen::Matrix <double, NUMBER_OF_WHEELS*(2*NUMAXIS), 8> matrixBType;
    typedef Eigen::Matrix <double, NUMBER_OF_WHEELS*(2*NUMAXIS), 12> matrixConjType; // columns are columns of A + 1
    
    matrixConjType Conj;
    
    Eigen::FullPivLU<matrixAType> lu_decompA(A);
    std::cout << "The rank of A is " << lu_decompA.rank() << std::endl;
	    
    Eigen::FullPivLU<matrixBType> lu_decompB(B);
    std::cout << "The rank of B is " << lu_decompB.rank() << std::endl;
    
    Conj.block<NUMBER_OF_WHEELS*(2*NUMAXIS),11>(0,0) = A;
    Conj.block<NUMBER_OF_WHEELS*(2*NUMAXIS), 1>(0,11) = b;
    Eigen::FullPivLU<matrixConjType> lu_decompConj(Conj);
    std::cout << "The rank of A|B*y is " << lu_decompConj.rank() << std::endl;
    std::cout << "Pseudoinverse of A\n" << (A.transpose() * A).inverse() << std::endl;
    std::cout << "Inverse of R\n" << R.inverse() << std::endl;
    #endif
    /** **/
    
    x = (A.transpose() * R.inverse() * A).inverse() * A.transpose() * R.inverse() * b;
    
    if (b.norm() != 0.00)
	leastSquaresError = (A*x - b).norm() / b.norm();
    #ifdef DEBUG_PRINTS
    std::cout << "[VM] The relative error is:\n" << leastSquaresError << std::endl;
    std::cout << "[VM] The solution is:\n" << x << std::endl;
    #endif
    
    
    /** Save the results **/
    velocity = x.block<NUMAXIS,1>(0,0);
    setCurrentVeloModel(velocity);
    acontact = x.block<NUMBER_OF_WHEELS,1>(NUMAXIS+NUMBER_OF_WHEELS,0);
    setContactAnglesVelocity(acontact);
    
    
    return;
}

void measurement::calculateNavigationKinematics(const Eigen::Matrix< double, Eigen::Dynamic, Eigen::Dynamic >& A, const Eigen::Matrix< double, Eigen::Dynamic, Eigen::Dynamic >& B)
{
    Eigen::Matrix <double, Eigen::Dynamic, Eigen::Dynamic> R;
    
    R.resize(NUMBER_OF_WHEELS*(2*NUMAXIS), NUMBER_OF_WHEELS*(2*NUMAXIS));
    
    /** Form the covariance matrix **/
//     R.setZero();
//     R.block<NUMAXIS, NUMAXIS>(0,0) = this->Rg;
//     R.block<ENCODERS_VECTOR_SIZE, ENCODERS_VECTOR_SIZE>(NUMAXIS,NUMAXIS) = this->Rencoders;
    R.setIdentity();
    
    /** Peform the LS solution **/
    this->navigationKinematics(A, B, R);
    
    return;
}


void measurement::slipKinematics(Eigen::Matrix< double, Eigen::Dynamic, Eigen::Dynamic > A, Eigen::Matrix< double, Eigen::Dynamic, Eigen::Dynamic > B,
				    Eigen::Matrix< double, Eigen::Dynamic, Eigen::Dynamic > R, double dt)
{
    double leastSquaresError;
    Eigen::Matrix< double, NUMAXIS , 1  > linvelo; /** Linear velocities from acc integration **/
    
    /** Navigation Vectors **/
    Eigen::Matrix< double, Eigen::Dynamic, Eigen::Dynamic > x;
    Eigen::Matrix< double, Eigen::Dynamic, Eigen::Dynamic > y;
    Eigen::Matrix< double, Eigen::Dynamic, Eigen::Dynamic > b;
    
    /** Resize vector to the correct size **/
    x.resize(SLIP_VECTOR_SIZE, 1);
    y.resize(2*NUMAXIS + ENCODERS_VECTOR_SIZE + NUMBER_OF_WHEELS, 1);
    
    /** Vector b **/
    b.resize(NUMBER_OF_WHEELS*(2*NUMAXIS), 1);
    
    /** Integrate the lin acceleration **/
    linvelo = accIntegrationWindow(dt);
    setLinearVelocities(linvelo);
    
    /** Form the sensed vector **/
    y.block<NUMAXIS, 1> (0,0) = linvelocity;
    y.block<NUMAXIS, 1> (NUMAXIS,0) = angvelocity;
    y.block<ENCODERS_VECTOR_SIZE, 1> (2*NUMAXIS,0) = encodersvelocity;
    y.block<NUMBER_OF_WHEELS, 1> ((2*NUMAXIS) + ENCODERS_VECTOR_SIZE,0) = acontact;
    
    #ifdef DEBUG_PRINTS
    std::cout<<"[SK] A is of size "<<A.rows()<<"x"<<A.cols()<<"\n";
    std::cout<<"[SK] A:\n" << A << std::endl;
    std::cout<<"[SK] B is of size "<<B.rows()<<"x"<<B.cols()<<"\n";
    std::cout<<"[SK] B:\n" << B << std::endl;
    std::cout<<"[SK] y is of size " <<y.rows()<<"x"<< y.cols()<<"\n";
    std::cout<<"[SK] The sensed vector y\n"<<y<<"\n";
    #endif
    
    b = B*y;
    
    /** DEBUG OUTPUT **/
    #ifdef DEBUG_PRINTS
    typedef Eigen::Matrix <double, NUMBER_OF_WHEELS*(2*NUMAXIS), SLIP_VECTOR_SIZE> matrixAType;
    typedef Eigen::Matrix <double, NUMBER_OF_WHEELS*(2*NUMAXIS), 2*NUMAXIS + ENCODERS_VECTOR_SIZE + NUMBER_OF_WHEELS> matrixBType;
    typedef Eigen::Matrix <double, NUMBER_OF_WHEELS*(2*NUMAXIS), SLIP_VECTOR_SIZE+1> matrixConjType; // columns are columns of A + 1
    
    matrixConjType Conj;
    
    Eigen::FullPivLU<matrixAType> lu_decompA(A);
    std::cout << "The rank of A is " << lu_decompA.rank() << std::endl;
	    
    Eigen::FullPivLU<matrixBType> lu_decompB(B);
    std::cout << "The rank of B is " << lu_decompB.rank() << std::endl;
    
    Conj.block<NUMBER_OF_WHEELS*(2*NUMAXIS),SLIP_VECTOR_SIZE>(0,0) = A;
    Conj.block<NUMBER_OF_WHEELS*(2*NUMAXIS), 1>(0,SLIP_VECTOR_SIZE) = b;
    Eigen::FullPivLU<matrixConjType> lu_decompConj(Conj);
    std::cout << "The rank of A|B*y is " << lu_decompConj.rank() << std::endl;
    std::cout << "Pseudoinverse of A\n" << (A.transpose() * A).inverse() << std::endl;
    #endif
    /** **/
    
    x = (A.transpose() * R.inverse() * A).inverse() * A.transpose() * R.inverse() * b;
    
    leastSquaresError = (A*x - b).norm() / b.norm();
    #ifdef DEBUG_PRINTS
    std::cout << "[SK] The relative error is:\n" << leastSquaresError << std::endl;
    std::cout << "[SK] The solution is:\n" << x << std::endl;
    #endif
    
    
    /** Save the results **/
    
    
    return;

}

