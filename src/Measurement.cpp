/**\file Measurement.cpp
 *
 * This class has the primitive methods for the Measurement Generation of the localization framework.
 * The class perform proprioceptive Measurements according to the kinematics jacobian
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

#include "Measurement.hpp"

#define DEBUG_PRINTS 1

using namespace localization;
using namespace Eigen;

/**
* @brief Default constructor
*/
Measurement::Measurement()
{

    this->testAcc = boost::circular_buffer<double>(2);
}

/**
* @brief Default desconstructor
*/
Measurement::~Measurement()
{

}


void Measurement::welcome()
{
	std::cout << "You successfully compiled and executed Measurement Generation. Welcome!" << std::endl;
}

/**
* @brief Set the current linear velocities
*/
void Measurement::setLinearVelocities(Eigen::Matrix< double, NUMAXIS , 1  > &linvelo)
{
    this->linvelocity = linvelo;
    
    return;
}

/**
* @brief Return linear velocities
*/
Eigen::Matrix< double, NUMAXIS , 1  > Measurement::getLinearVelocities()
{
    return this->linvelocity;
}

/**
* @brief Set the current angular velocity
*/
void Measurement::setAngularVelocities(Eigen::Matrix< double, NUMAXIS , 1  > &angvelo)
{
    this->cbAngvelo = angvelo;
    
    return;
}

/**
* @brief Return angular velocities
*/
Eigen::Matrix< double, NUMAXIS , 1  > Measurement::getAngularVelocities()
{
    
    return this->cbAngvelo;
}

/**
* @brief Set the linear Acceleration
*/
void Measurement::setLinearAcceleration(Eigen::Matrix< double, NUMAXIS , 1  > &linacc)
{
    
    this->cbAcc = linacc;
    
    return;
}


/**
* @brief Return linear acceleration
*/
Eigen::Matrix< double, NUMAXIS , 1  > Measurement::getLinearAcceleration()
{
        
    return this->cbAcc;
}


/**
* @brief Return complete slip vector 
*/
Eigen::Matrix< double, SLIP_VECTOR_SIZE , 1  > Measurement::getSlipVector()
{
    return slipVector.data;
}

/**
* @brief Return complete slip vector Covariance 
*/
Eigen::Matrix< double, SLIP_VECTOR_SIZE , SLIP_VECTOR_SIZE  > Measurement::getSlipVectorCov()
{
    return slipVector.Cov;
}


/**
* @brief Return slip vector for wheel_idx
*/
Eigen::Matrix< double, NUMAXIS , 1  > Measurement::getSlipVector(const unsigned int wheel_idx)
{
    if (wheel_idx < NUMBER_OF_WHEELS)
	return this->slipMatrix.col(wheel_idx);
    
    return Eigen::Matrix <double, NUMAXIS, 1>::Zero();
}

/**
* @brief Gets rover slip error vector
*/
Eigen::Matrix< double, SLIP_VECTOR_SIZE, 1  > Measurement::getSlipErrorVector()
{
    return slipError.data;
}

/**
* @brief Gets rover slip error covariance
*/
Eigen::Matrix< double, SLIP_VECTOR_SIZE, SLIP_VECTOR_SIZE > Measurement::getSlipErrorVectorCovariance()
{
    return slipError.Cov;
}


/**
* @brief Set the current contact angles
*/
void Measurement::setContactAnglesVelocity(Eigen::Matrix<double, NUMBER_OF_WHEELS, 1> &contact_angles, Eigen::Matrix<double, NUMBER_OF_WHEELS, NUMBER_OF_WHEELS> &Cov)
{
    this->acontact.data = contact_angles;
    this->acontact.Cov = Cov;
}


/**
* @brief Return the contact angles
*/
Eigen::Matrix< double, Eigen::Dynamic, 1  > Measurement::getContactAnglesVelocity()
{
    return this->acontact.data;
}

/**
* @brief Set the current velocity model
*/
void Measurement::setCurrentVeloModel(Eigen::Matrix< double, NUMAXIS , 1  > &velocity)
{
    this->cbVelModel = velocity;
}

/**
* @brief Return the current velocity from odometry model
*/
Eigen::Matrix< double, NUMAXIS , 1  > Measurement::getCurrentVeloModel()
{
    
    return this->cbVelModel;
}

/**
* @brief Get the covariance of the velocity from the odometry model
*/
Eigen::Matrix< double, NUMAXIS , NUMAXIS> Measurement::getCurrentVeloModelCovariance()
{
    return velModel.Cov;
}

/**
* @brief Get the increment in velocity
*/
Eigen::Matrix< double, NUMAXIS, 1  > Measurement::getIncrementalVeloModel()
{
    
    return increVelModel.data;
}

/**
* @brief Get the covariance of the increment in velocity
*/
Eigen::Matrix< double, NUMAXIS, NUMAXIS> Measurement::getIncrementalVeloModelCovariance()
{
    return increVelModel.Cov;
}

/**
* @brief Set the current encoders velocities
*/
void Measurement::setEncodersVelocity(const Eigen::Matrix< double, Eigen::Dynamic, 1  > &vjoints)
{
    this->encodersvelocity = vjoints;
    
    #ifdef DEBUG_PRINTS
    std::cout<<"encodersvelocity\n"<<this->encodersvelocity<<"\n";
    #endif
    
    return;
}

/**
* @brief Returns the covariance noise matrix
*/
Eigen::Matrix< double, ENCODERS_VECTOR_SIZE, ENCODERS_VECTOR_SIZE > Measurement::getEncodersVelocityCovariance()
{
    return Rencoders;
}

/**
* @brief Returns the covariance noise matrix
*/
Eigen::Matrix< double, NUMBER_OF_WHEELS, NUMBER_OF_WHEELS > Measurement::getContactAnglesVelocityCovariance()
{
    return acontact.Cov;
}

/**
* @brief This function set the Accelerometers excentricity
*/
void Measurement::setEccentricity(Eigen::Matrix <double,NUMAXIS,1>  &eccx, Eigen::Matrix <double,NUMAXIS,1>  &eccy, Eigen::Matrix <double,NUMAXIS,1>  &eccz)
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
* @brief Returns the quaternion with the projection for the nadir vector
*/
Eigen::Quaternion< double > Measurement::getLevelWeightDistribution()
{
    return q_weight_distribution;
}


/**
* @brief Initialization method
*/
void Measurement::Init(Eigen::Matrix< double, ENCODERS_VECTOR_SIZE , ENCODERS_VECTOR_SIZE  >& Ren,
		       Eigen::Matrix< double, NUMBER_OF_WHEELS , NUMBER_OF_WHEELS  >& Rcont,
		       Eigen::Quaternion <double> q_weight)
{
    /** Initial encoders velocities **/
    encodersvelocity.setZero();
    
    /** Initial contact angle **/
    acontact = DataModel(NUMBER_OF_WHEELS);
    
    /** Velocity model from navigation kinematics **/
    velModel = DataModel(NUMAXIS);
    prevVelModel = DataModel(NUMAXIS);
    
    /** Incremental velocity model from navigation kinematics **/
    increVelModel = DataModel(NUMAXIS);
    
    /** Slip model **/
    slipModel = DataModel(SLIP_VECTOR_SIZE);
    slipVector = DataModel(SLIP_VECTOR_SIZE);
    slipError = DataModel(SLIP_VECTOR_SIZE);
    
    /** Initial slip matrix **/
    slipMatrix = Eigen::Matrix <double, NUMAXIS, NUMBER_OF_WHEELS>::Zero();
    
    /** Velocities for the filter (acceleration and bias substracted correctly) **/
    linvelocity = Eigen::Matrix <double, NUMAXIS, 1>::Zero();
    
    /** Set the eneven weight quaternion of the nadir vector to the rover platform **/
    q_weight_distribution = q_weight;
    
    /** Vector for the integration **/
    cbAcc.setZero();
    cbAngvelo.setZero();
    
    /** Resize the circular buffer for the velocities model **/
    cbVelModel.setZero();
    
    testAcc.set_capacity(INTEGRATION_XAXIS_WINDOW_SIZE);
    
    this->Rencoders = Ren;
    
    this->leastSquaredNavigation = std::numeric_limits<double>::quiet_NaN();
    this->leastSquaredSlip = std::numeric_limits<double>::quiet_NaN();
    
    /** Print filter information **/
    #ifdef DEBUG_PRINTS
    std::cout<< "Ren is of size "<<Rencoders.rows()<<"x"<<Rencoders.cols()<<"\n";
    std::cout<< "Ren:\n"<<Rencoders<<"\n";
    std::cout<< "acontact is of size "<<acontact.data.size()<<"\n";
    std::cout<< "acontact:\n"<<acontact<<"\n";
    std::cout<< "Rcontact is of size "<<acontact.Cov.rows()<<"x"<<acontact.Cov.cols()<<"\n";
    std::cout<< "Rcontact:\n"<<acontact.Cov<<"\n";
    #endif
}





/**
* @brief This perform simpson rule for integration
*/
double Measurement::simpsonsIntegral(double fa, double fm, double fb, double delta_ab)
{
    return delta_ab/6.0 * (fa + (4.0*fm) + fb);
}

/**
* @brief To obtain the linear velocity 
*/
Eigen::Matrix< double, NUMAXIS , 1  > Measurement::accIntegrationWindow(double dt)
{
    Eigen::Matrix <double, NUMAXIS, 1> localvelocity;
    Eigen::Matrix <double,NUMAXIS,NUMAXIS> gyros2product; /** Vec 2 product  matrix for the gyroscopes (angular velocity) */
    
    localvelocity.setZero();
    
    #ifdef DEBUG_PRINTS
    std::cout<< "IMU Velocity\n";
    #endif
    
    gyros2product << 0, -cbAngvelo[2], cbAngvelo[1],
	cbAngvelo[2], 0, -cbAngvelo[0],
	-cbAngvelo[1], cbAngvelo[0], 0;
		
    localvelocity[0] += (cbAcc[0] * dt);// - (gyros2product.row(0) * eccx);
    
//     localvelocity[0] += cbVelModel[0];
    
    #ifdef DEBUG_PRINTS
    std::cout<< "localvelocity[0] "<<localvelocity[0]<<"\n";
    #endif
    
    	
    localvelocity[1] += (cbAcc[1] * dt);// - (gyros2product.row(1) * eccy);
    
//     localvelocity[1] += cbVelModel[1];
    
    #ifdef DEBUG_PRINTS
    std::cout<< "localvelocity[1] "<<localvelocity[1]<<"\n";
    #endif
    
    localvelocity[2] += (cbAcc[2] * dt);// - (gyros2product.row(2) * eccz);
    
//     localvelocity[2] += cbVelModel[2];
    
    #ifdef DEBUG_PRINTS
    std::cout<< "localvelocity[2] "<<localvelocity[2]<<"\n";
    #endif
    
    return localvelocity;
}

/**
* @brief Return the buffer size for x,y and z
*/
Eigen::Matrix< double, NUMAXIS , 1  > Measurement::getIntegrationWindowSize()
{
    Eigen::Matrix< double, NUMAXIS , 1  > windowsSize;
    
//     windowsSize << cbAccX.size(), cbAccY.size(), cbAccZ.size();
    windowsSize.setIdentity();
    
    return windowsSize;
}

double Measurement::navigationKinematics(const Eigen::Matrix< double, Eigen::Dynamic, Eigen::Dynamic > &A, const Eigen::Matrix< double, Eigen::Dynamic, Eigen::Dynamic > &B,
					const Eigen::Matrix< double, Eigen::Dynamic, Eigen::Dynamic > &R, const Eigen::Matrix< double, Eigen::Dynamic, Eigen::Dynamic > &W)
{
    double leastSquaresError = std::numeric_limits<double>::quiet_NaN();
    Eigen::Matrix< double, NUMAXIS , 1  > velocity;
    Eigen::Matrix <double, SLIP_VECTOR_SIZE, 1> slip;
    Eigen::Matrix <double, NUMAXIS+(2*NUMBER_OF_WHEELS), NUMAXIS+(2*NUMBER_OF_WHEELS)> Iinverse;
    
    /** Navigation Vectors **/
    Eigen::Matrix< double, Eigen::Dynamic, Eigen::Dynamic > x;
    Eigen::Matrix< double, Eigen::Dynamic, Eigen::Dynamic > y;
    Eigen::Matrix< double, Eigen::Dynamic, Eigen::Dynamic > b;
    
    /** Resize vector to the correct size **/
    x.resize(NUMAXIS+(2*NUMBER_OF_WHEELS), 1);
    y.resize(NUMAXIS+ENCODERS_VECTOR_SIZE, 1);
    b.resize(NUMBER_OF_WHEELS*(2*NUMAXIS), 1);
    
    /** Form the sensed vector **/
    y.block<NUMAXIS, 1> (0,0) = this->getAngularVelocities();
    y.block<ENCODERS_VECTOR_SIZE, 1> (NUMAXIS,0) = encodersvelocity;
    
    #ifdef DEBUG_PRINTS    
    std::cout<<"[NK] A is of size "<<A.rows()<<"x"<<A.cols()<<"\n";
    std::cout<<"[NK] A:\n" << A << std::endl;
    std::cout<<"[NK] B is of size "<<B.rows()<<"x"<<B.cols()<<"\n";
    std::cout<<"[NK] B:\n" << B << std::endl;
    std::cout<<"[NK] y is of size " <<y.rows()<<"x"<< y.cols()<<"\n";
    std::cout<<"[NK] The sensed vector y\n"<<y<<"\n";
    std::cout<<"[NK] R is of size "<<R.rows()<<"x"<<R.cols()<<"\n";
    std::cout<<"[NK] R:\n" << R << std::endl;
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
    std::cout << "Pseudoinverse of A\n" << (A.transpose() * (W + R.inverse()) * A).inverse() << std::endl;
    #endif
    /** **/
   
    /** Least-Squares solution **/
    Iinverse = (A.transpose() * (W + R.inverse()) * A).inverse();
    x = Iinverse * A.transpose() * (W + R.inverse()) * b;
    
//     Eigen::MatrixXd M; /** dynamic memory matrix for the solution **/
//     M = A;
//     x = M.jacobiSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(b);
    
    Iinverse = (A.transpose() * (R.inverse()) * A).inverse();
    
    Eigen::Matrix<double, 1,1> squaredError = (((A*x - b).transpose() * R * (A*x - b)));
    Eigen::Matrix<double, 2*NUMAXIS, NUMBER_OF_WHEELS*(2*NUMAXIS)> idBlock;
    for (register unsigned int i=0; i<NUMBER_OF_WHEELS;i++)
    {
	idBlock.block<2*NUMAXIS, 2*NUMAXIS>(0,i*(2*NUMAXIS)).setIdentity();// = Eigen::Matrix<double, 2*NUMAXIS, 2*NUMAXIS>::Identity();
    }
    
    Eigen::Matrix<double, 2*NUMAXIS,1> vectorError = (idBlock * (R * (A*x - b)));
    
    if (b.norm() != 0.00)
	leastSquaresError = sqrt(squaredError[0]) / b.norm();
    #ifdef DEBUG_PRINTS
    std::cout << "[NK] Identity 6x24:\n" << idBlock << std::endl;
    std::cout << "[NK] The relative error is:\n" << leastSquaresError << std::endl;
    std::cout << "[NK] The absolute squared error is:\n" << squaredError << std::endl;
    std::cout << "[NK] The vector error is:\n" << vectorError << std::endl;
    std::cout << "[NK] The vector error/norm is:\n" << vectorError/b.norm() << std::endl;
    std::cout << "[NK] The sqrt(vectorerror)/norm is:\n" << vectorError.array().cwiseSqrt()/b.norm() << std::endl;
    std::cout << "[NK] The solution is:\n" << x << std::endl;
    std::cout << "[NK] The uncertainty is:\n" << Iinverse << std::endl;
    #endif
    
    /** Save the velocity model results **/
    velocity = x.block<NUMAXIS,1>(0,0);
    velModel.data = velocity;
    velModel.Cov = pow(leastSquaresError,2) * Eigen::Matrix<double, NUMAXIS, NUMAXIS>::Identity() +  Iinverse.block<NUMAXIS, NUMAXIS> (0,0);
    setCurrentVeloModel(velocity);
    
    #ifdef DEBUG_PRINTS
    std::cout << "[NK] velModel.data:\n" << velModel.data << std::endl;
    std::cout << "[NK] velModel.Cov:\n" << velModel.Cov << std::endl;
    #endif
    
    /** Compute the increment in velocity **/
    increVelModel = velModel - prevVelModel;
    
    /** Set the prev value to the actual one for the next iteration **/
    prevVelModel = velModel;
    
    /** Save the z-axis slip vector results **/
    slip.setZero();
    slip(2) = x(3); slip(5) = x(4); slip(8) = x(5); slip(11) = x(6);
    slipModel.data = slip;
    slipModel.Cov(2,2) = Iinverse(3,3);
    slipModel.Cov(5,5) = Iinverse(4,4);
    slipModel.Cov(8,8) = Iinverse(5,5);
    slipModel.Cov(11,11) = Iinverse(6,6);
    
    /** Store the acontact variable **/
    acontact.data = x.block<NUMBER_OF_WHEELS,1>(NUMAXIS+NUMBER_OF_WHEELS,0);
    acontact.Cov = Iinverse.bottomRightCorner<NUMBER_OF_WHEELS,NUMBER_OF_WHEELS>();
    
    this->leastSquaredNavigation = leastSquaresError;
    
    return leastSquaresError;
}

double Measurement::slipKinematics(const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> &A, const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> &B,
				const Eigen::Matrix< double, Eigen::Dynamic, Eigen::Dynamic > &R)
{
    double leastSquaresError = std::numeric_limits<double>::quiet_NaN();
    
    /** Navigation Vectors **/
    Eigen::Matrix< double, Eigen::Dynamic, Eigen::Dynamic > x;
    Eigen::Matrix< double, Eigen::Dynamic, Eigen::Dynamic > y;
    Eigen::Matrix< double, Eigen::Dynamic, Eigen::Dynamic > b;
    Eigen::Matrix <double, SLIP_VECTOR_SIZE, SLIP_VECTOR_SIZE> Iinverse;
    
    /** Resize vector to the correct size **/
    x.resize(SLIP_VECTOR_SIZE, 1);
    y.resize(2*NUMAXIS + ENCODERS_VECTOR_SIZE + NUMBER_OF_WHEELS, 1);
    b.resize(NUMBER_OF_WHEELS*(2*NUMAXIS), 1);
    
    
    /** Form the sensed vector **/
    y.block<NUMAXIS, 1> (0,0) = linvelocity;
    y.block<NUMAXIS, 1> (NUMAXIS,0) = this->getAngularVelocities();
    y.block<ENCODERS_VECTOR_SIZE, 1> (2*NUMAXIS,0) = encodersvelocity;
    y.block<NUMBER_OF_WHEELS, 1> ((2*NUMAXIS) + ENCODERS_VECTOR_SIZE,0) = acontact.data;
    

    #ifdef DEBUG_PRINTS
    std::cout<<"[SK] A is of size "<<A.rows()<<"x"<<A.cols()<<"\n";
    std::cout<<"[SK] A:\n" << A << std::endl;
    std::cout<<"[SK] B is of size "<<B.rows()<<"x"<<B.cols()<<"\n";
    std::cout<<"[SK] B:\n" << B << std::endl;
    std::cout<<"[SK] y is of size " <<y.rows()<<"x"<< y.cols()<<"\n";
    std::cout<<"[SK] The sensed vector y\n"<<y<<"\n";
    std::cout<<"[SK] R is of size "<<R.rows()<<"x"<<R.cols()<<"\n";
    std::cout<<"[SK] R:\n" << R << std::endl;
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
    std::cout << "Pseudoinverse of A\n" << (A.transpose() * R.inverse() * A).inverse() << std::endl;
    #endif
    /** **/
    
    Iinverse = (A.transpose() * R.inverse() * A).inverse();
    x =  Iinverse * A.transpose() * R.inverse() * b;
    
    Eigen::Matrix<double, 1,1> squaredError = (((A*x - b).transpose() * R * (A*x - b)));
    
    if (b.norm() != 0.00)
	leastSquaresError = sqrt(squaredError[0]) / b.norm();
    
    /** Save the results **/
    slipVector.data = x;
    slipVector.Cov = Iinverse;
    
    /** Compute the slip vector error **/
    slipError.data = slipVector.data - slipModel.data;
    slipError.Cov = slipVector.Cov;
    
    #ifdef DEBUG_PRINTS
    std::cout << "[SK] The relative error is:\n" << leastSquaresError << std::endl;
    std::cout << "[SK] The absolute squared error is:\n" << squaredError << std::endl;
    std::cout << "[SK] The solution is:\n" << x << std::endl;
    std::cout << "[SK] The uncertainty is:\n" << Iinverse << std::endl;
    std::cout << "[SK] SlipInertial is:\n" << slipVector << std::endl;
    std::cout << "[SK] SlipModel is:\n" << slipModel << std::endl;
    std::cout << "[SK] SlipError(+) is:\n" << slipError << std::endl;
    std::cout << "[SK] SlipError(-) is:\n" << slipVector.data - slipModel.data << std::endl;
    std::cout << "[SK] SlipError(+) is:\n" << slipVector.data + slipModel.data << std::endl;
    #endif
    
    this->leastSquaredSlip = leastSquaresError;
    
    return leastSquaresError;

}

/**
* @brief Save slip info to the Orogen-compatible DataType
*/
void Measurement::toSlipInfo(SlipInfo& sinfo)
{

    sinfo.slip_vector = slipVector.data;
    sinfo.Cov = slipVector.Cov;
    
    return;
}

void Measurement::toMeasurementGenerationInfo(MeasurementGenerationInfo& measurementInfo)
{
    measurementInfo.wlsNavigation = this->leastSquaredNavigation;
    measurementInfo.wlsSlip = this->leastSquaredSlip;

    return;
}



/**
* @brief Computes the Bhattacharyya coefficient
*/
Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> Measurement::bhattacharyya(DataModel& data1, DataModel& data2)
{
    double number_expo = std::numeric_limits<double>::quiet_NaN();
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> Inv_half(data1.size(), data1.size());
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> BC(data1.size(), data1.size());
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> Cov1llt(data1.size(), data1.size());
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> Cov2llt(data2.size(), data2.size());
    DataModel aux(data1.size());
    
    /** Substraction **/
    aux = data1 - data2;
    Inv_half = (aux.Cov * 0.5).inverse();
    
    Eigen::Matrix<double, 1, 1> expo = -0.125 * (aux.data.transpose() * Inv_half * aux.data);
    number_expo = exp(expo[0]);
    
    Cov1llt = data1.Cov.llt().matrixL();
    Cov2llt = data2.Cov.llt().matrixL();
    Cov1llt = Cov1llt.llt().matrixL();
    Cov2llt = Cov2llt.llt().matrixL();
        
    BC = (Cov1llt*Cov2llt) * Inv_half.llt().matrixL();
    
    #ifdef DEBUG_PRINTS
    std::cout << "[BHATTACHARYYA] number_expo is:\n" << number_expo << std::endl;
    std::cout << "[BHATTACHARYYA] Cov1llt is:\n" << Cov1llt << std::endl;
    std::cout << "[BHATTACHARYYA] Cov2llt is:\n" << Cov2llt << std::endl;
    std::cout << "[BHATTACHARYYA] Inv_half is:\n" << Inv_half << std::endl;
    std::cout << "[BHATTACHARYYA] BC(without expo) is:\n" << BC << std::endl;
    #endif
    
    BC = BC * number_expo;
    
    #ifdef DEBUG_PRINTS
    std::cout << "[BHATTACHARYYA] BC is:\n" << BC << std::endl;
    #endif
    
    return BC;
}

/**
* @brief Computes the Mahalanobis distance
*/
double Measurement::mahalanobis(DataModel& data1)
{
    Eigen::Matrix<double, 1, 1> result;
    
    result[0] = data1.data.transpose() * data1.Cov.inverse() * data1.data;
    return sqrt(result[0]);
}



/**
* @brief Convert data value in the range of MinValues..MaxValues to the range 350 - 650
*/
double Measurement::getWaveLenghtFromValue (const double value, const double max, const double min)
{
    const double minVisibleLength = 350.0;
    const double maxVisibleLength = 650.0;
    
    return (value - min)/(max-min)*(maxVisibleLength-minVisibleLength) + minVisibleLength;
    
}


int colorAdjust (const double colorvalue, const double factor)
{
    const double gamma        =   0.80;
    const double intensitymax = 255.0;
    
    if (colorvalue == 0.00)
	return 0;
    else
	return round(intensitymax * pow(colorvalue * factor, gamma));
    
}

/**
* @brief Convert data value in the range of MinValues..MaxValues to the range 350 - 650
*/
Eigen::Matrix<double, 3, 1> Measurement::waveLenghtToColor (const double wavelength)
{
    const double intensitymax = 255.0;
    
    Eigen::Matrix<double, 3, 1> color;
    double blue, factor, green, red;
    
    if ((wavelength > 379.00)&&(wavelength < 440.00))
    {
	    red = -(wavelength - 440.0) / (440.0 - 380.0);
	    green = 0.00;
	    blue = 1.00;
    }    
    else if ((wavelength > 439.00)&&(wavelength < 490.00))
    {
	    red = 0.00;
	    green = (wavelength - 440.0) / (490.0 - 440.0);
	    blue = 1.00;
    }    
    else if ((wavelength > 489.00)&&(wavelength < 510.00))
    {
	    red = 0.00;
	    green = 1.00;
	    blue = -(wavelength - 510.0) / (510 - 490);
    }    
    else if ((wavelength > 509.00)&&(wavelength < 580.00))
    {
	    red = (wavelength - 510.0) / (580 - 510);
	    green = 1.00;
	    blue = 0.00;
    }    
    else if ((wavelength > 579.00)&&(wavelength < 645.00))
    {
	    red = 1.00;
	    green = -(wavelength - 645.0) / (645 - 580);
	    blue = 0.00;
    }    
    else if ((wavelength > 644.00)&&(wavelength < 781.00))
    {
	    red = 1.00;
	    green = 0.00;
	    blue = 0.00;
    }
    else if (wavelength >781.00)
    {
	    red = 1.00;
	    green = 0.0;
	    blue = 0.0;
    }
    else
    {
	    red = 0.0;
	    green = 0.0;
	    blue = 1.00;
    }
    
    if ((wavelength > 379.00)&&(wavelength < 420.00))
    {
	    factor = 0.3 + 0.7*(wavelength - 380) / (420 - 380);
    }
    else if ((wavelength > 419.00)&&(wavelength < 701.00))
    {
	    factor = 1.0;
    }
    else if ((wavelength > 700.00)&&(wavelength < 781.00))
    {
	    factor = 0.3 + 0.7*(780.00 - wavelength) / (780 - 700);
    }
    else
    {
	    factor = 0.00;
    }
    
    color[0] = colorAdjust(red, factor) / intensitymax;
    color[1] = colorAdjust(green, factor) / intensitymax;
    color[2] = colorAdjust(blue, factor) / intensitymax;
    
    return color;
}


/**
* @brief Convert data value in the range of MinValues..MaxValues to the range 350 - 650
*/
Eigen::Matrix<double, 3, 1> Measurement::waveLenghtToColorv2 (const double wavelength)
{
    const double intensitymax = 255.0;
    
    Eigen::Matrix<double, 3, 1> color;
    double blue, factor, green, red;
    
    if ((wavelength > 379.00)&&(wavelength < 580.00))
    {
	    red = (wavelength - 510.0) / (580 - 510);
	    green = 1.00;
	    blue = 0.00;
    }    
    else if ((wavelength > 579.00)&&(wavelength < 645.00))
    {
	    red = 1.00;
	    green = -(wavelength - 645.0) / (645 - 580);
	    blue = 0.00;
    }    
    else if ((wavelength > 644.00)&&(wavelength < 781.00))
    {
	    red = 1.00;
	    green = 0.00;
	    blue = 0.00;
    }
    else if (wavelength >781.00)
    {
	    red = 1.00;
	    green = 0.0;
	    blue = 0.0;
    }
    else
    {
	    red = 0.0;
	    green = 1.0;
	    blue = 0.00;
    }
    
    if ((wavelength > 379.00)&&(wavelength < 420.00))
    {
	    factor = 0.3 + 0.7*(wavelength - 380) / (420 - 380);
    }
    else if ((wavelength > 419.00)&&(wavelength < 701.00))
    {
	    factor = 1.0;
    }
    else if ((wavelength > 700.00)&&(wavelength < 781.00))
    {
	    factor = 0.3 + 0.7*(780.00 - wavelength) / (780 - 700);
    }
    else
    {
	    factor = 0.00; //! Balck color in envire
    }
    
    color[0] = colorAdjust(red, factor) / intensitymax;
    color[1] = colorAdjust(green, factor) / intensitymax;
    color[2] = colorAdjust(blue, factor) / intensitymax;
    
    return color;
    
}
/**
* @brief 
*/
Eigen::Matrix<double, 3, 1> Measurement::valueToColor (const double value, const double max, const double min)
{
    double wavelength;//! In nanometers
    
    wavelength = this->getWaveLenghtFromValue(value, max, min);
    #ifdef DEBUG_PRINTS
    std::cout<< "[SLIP2COLOR] wavelength "<<wavelength<<"\n";
    #endif
    
    /** Restrict the wavelenght to be between Green and Red **/
    if (wavelength < 550.00)
	wavelength = 550.00;
    else if (wavelength > 780.00)
	wavelength = 780.00;
    
    return this->waveLenghtToColor(wavelength);
}

