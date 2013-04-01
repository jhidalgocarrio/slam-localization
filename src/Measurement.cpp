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
    /** In constructor set the buffer size to NaN **/
    this->ibuffer_size = std::numeric_limits<double>::quiet_NaN();

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


void Measurement::setInertialValues(const Eigen::Matrix< double, NUMAXIS , 1  >& acceleration, const Eigen::Matrix< double, NUMAXIS , 1  >& angular_velocity)
{

    /** Check if the buffer is full **/
    if (acc.size() >= this->ibuffer_size)
	acc.pop_back();
    
    if (angvelo.size() >= this->ibuffer_size)
	angvelo.pop_back();
    
    /** Store the new values **/
    acc.push_front(acceleration);
    angvelo.push_front(angular_velocity);
    
    #ifdef DEBUG_PRINTS
    std::cout << "[SETIN] angvelo is of size "<< angvelo.size()<<"\n";
    std::cout << "[SETIN] acc is of size "<< acc.size()<<"\n";
    #endif
    
    return;
}

/**
* @brief Gets Inertial Values (acc and angvelo)
*/
void Measurement::getInertialValues(Eigen::Matrix< double, NUMAXIS , 1  >& acc, Eigen::Matrix< double, NUMAXIS , 1  >& angvelo)
{
    acc = this->acc.front();
    angvelo = this->angvelo.front();
    
    return;
}

/**
* @brief Gets Avg Inertial Values (acc and angvelo)
*/
void Measurement::getStepInertialValues(Eigen::Matrix< double, NUMAXIS , 1  >& acc, Eigen::Matrix< double, NUMAXIS , 1  >& angvelo)
{
    Eigen::Matrix< double, NUMAXIS , 1  > a, g;
    int queue_size = std::min<int>(this->acc.size(), this->angvelo.size());
    a.setZero(); g.setZero();

    /** Get the mean of the stores values **/
    for (register int i=0; i<queue_size; i++)
    {
	a += this->acc[i];
	g += this->angvelo[i];
    }
    
    acc = (a/queue_size);
    angvelo = (g/queue_size);
    
    #ifdef DEBUG_PRINTS
    std::cout << "[GETAVGIN] Number of elements:\n" << queue_size<< std::endl;
    std::cout << "[GETAVGIN] The angevelo is:\n" << angvelo << std::endl;
    std::cout << "[GETAVGIN] The acc is:\n" << acc << std::endl;
    #endif
    
    return;
}

/**
* @brief Return increment in linear velocities
*/
Eigen::Matrix< double, NUMAXIS , 1  > Measurement::getLinearVelocities(double dt)
{
    Eigen::Matrix< double, NUMAXIS , 1  > a;
    int queue_size = this->acc.size();
    a.setZero();

    /** Get the mean of the stores values **/
    for (register int i=0; i<queue_size; i++)
	a += this->acc[i];
    
    return (a/queue_size)*dt;
}

/**
* @brief Return the current velocity from odometry model
*/
Eigen::Matrix< double, NUMAXIS , 1  > Measurement::getCurrentVeloModel()
{
    
    return this->velMM.front().data;
}

/**
* @brief Get the covariance of the velocity from the odometry model
*/
Eigen::Matrix< double, NUMAXIS , NUMAXIS> Measurement::getCurrentVeloModelCovariance()
{
    return velMM.front().Cov;
}

/**
* @brief Get the Step increment in velocity
*/
Eigen::Matrix< double, NUMAXIS, 1  > Measurement::getIncrementalVeloModel()
{
    Eigen::Matrix< double, NUMAXIS , 1  > iVelo;
    int queue_size = this->increVelMM.size();
    iVelo.setZero();

    /** Get the mean of the stores values **/
    for (register int i=0; i<queue_size; i++)
	iVelo += this->increVelMM[i].data;
    
    return iVelo;
}

/**
* @brief Get the covariance of the step increment in velocity
*/
Eigen::Matrix< double, NUMAXIS, NUMAXIS> Measurement::getIncrementalVeloModelCovariance()
{
    Eigen::Matrix< double, NUMAXIS , NUMAXIS  > iVeloCov;
    int queue_size = this->increVelMM.size();
    iVeloCov.setZero();

    /** Get the mean of the stores values **/
    for (register int i=0; i<queue_size; i++)
	iVeloCov += this->increVelMM[i].Cov;
    
    return iVeloCov;
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
* @brief Return the contact angles
*/
Eigen::Matrix< double, Eigen::Dynamic, 1  > Measurement::getContactAnglesVelocity()
{
    return this->dcontactAngle.data;
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
    return dcontactAngle.Cov;
}


/**
* @brief Returns the quaternion with the projection for the nadir vector
*/
Eigen::Quaternion< double > Measurement::getLevelWeightDistribution()
{
    return q_weight_distribution;
}

/**
* @brief Get angular velocity covariance matrix
*/
Eigen::Matrix< double, NUMAXIS , NUMAXIS  > Measurement::getRangvelo()
{
    return this->Rangvelo;
}

/**
* @brief Get joints encoders velocity covariance matrix
*/
Eigen::Matrix< double, ENCODERS_VECTOR_SIZE , ENCODERS_VECTOR_SIZE  > Measurement::getRencoders()
{
    return this->Rencoders;
}


/**
* @brief Initialization method
*/
void Measurement::Init(Eigen::Matrix< double, ENCODERS_VECTOR_SIZE , ENCODERS_VECTOR_SIZE  >& Ren,
		       Eigen::Matrix< double, NUMAXIS , NUMAXIS  >& Rangvelo,
		       Eigen::Quaternion <double> q_weight, int buffer_size)
{
        
    /** Initial contact angle **/
    dcontactAngle = DataModel(NUMBER_OF_WHEELS);
    
    /** Set the eneven weight quaternion of the nadir vector to the rover platform **/
    q_weight_distribution = q_weight;
       
    /** Initial encoders velocities **/
    encodersvelocity.setZero();
    
    /** Noise matrices for the proprioceptive measurements in the motion model **/
    this->Rencoders = Ren;
    this->Rangvelo = Rangvelo;
    
    /** Set the ibuffer_size variable **/
    this->ibuffer_size = buffer_size;
    
    /** Set the size of the Filter coeffcients **/
    besselBCoeff.resize(localization::NORDER_BESSEL_FILTER+1, 1);
    besselACoeff.resize(localization::NORDER_BESSEL_FILTER+1, 1);
    
    /** IIR filter coefficients **/ /** TO_DO: move to method a argument **/
    besselBCoeff[0] = 0.00467048;
    besselBCoeff[1] = 0.03736385;
    besselBCoeff[2] = 0.13077349;
    besselBCoeff[3] = 0.26154698;
    besselBCoeff[4] = 0.32693372;
    besselBCoeff[5] = 0.26154698;
    besselBCoeff[6] = 0.13077349;
    besselBCoeff[7] = 0.03736385;
    besselBCoeff[8] = 0.00467048;
    
    besselACoeff[0] = 1.00000000e+00;
    besselACoeff[1] = -3.87747570e-01;
    besselACoeff[2] = 7.13520818e-01;
    besselACoeff[3] = -2.49594003e-01;
    besselACoeff[4] = 1.47736180e-01;
    besselACoeff[5] = -3.59003821e-02;
    besselACoeff[6] = 8.56259334e-03;
    besselACoeff[7] = -9.97047726e-04;
    besselACoeff[8] = 6.27404353e-05;
         
    /** Initialize the LS-Error **/
    this->leastSquaredNavigation = std::numeric_limits<double>::quiet_NaN();
    this->leastSquaredSlip = std::numeric_limits<double>::quiet_NaN();
    
    /** Set the size for the Velocity model array buffer for bessel filter **/
    this->arrayVelMM.resize(NUMAXIS, localization::NORDER_BESSEL_FILTER+1);
    this->arrayVelMMIIR.resize(NUMAXIS, localization::NORDER_BESSEL_FILTER+1);
    
    /** Initialize the VelModel array for bessel filter **/
    this->arrayVelMM.setZero();
    this->arrayVelMMIIR.setZero();
    
    /** KK VARIABLES **/
    /** Slip model **/
    slipModel = DataModel(SLIP_VECTOR_SIZE);
    slipVector = DataModel(SLIP_VECTOR_SIZE);
    slipError = DataModel(SLIP_VECTOR_SIZE);
    
    /** Initial slip matrix **/
    slipMatrix = Eigen::Matrix <double, NUMAXIS, NUMBER_OF_WHEELS>::Zero();
    
    /** Print filter information **/
    #ifdef DEBUG_PRINTS
    std::cout<< "[MEASU_INIT] ibuffer_size "<<ibuffer_size<<"\n";
    std::cout<< "[MEASU_INIT] Ren is of size "<<Rencoders.rows()<<"x"<<Rencoders.cols()<<"\n";
    std::cout<< "[MEASU_INIT] Ren:\n"<<Rencoders<<"\n";
    std::cout<< "[MEASU_INIT] Rangvelo is of size "<<Rangvelo.rows()<<"x"<<Rangvelo.cols()<<"\n";
    std::cout<< "[MEASU_INIT] Rangvelo:\n"<<Rangvelo<<"\n";
    std::cout<< "[MEASU_INIT] dcontactAngle is of size "<<dcontactAngle.data.size()<<"\n";
    std::cout<< "[MEASU_INIT] dcontactAngle:\n"<<dcontactAngle<<"\n";
    std::cout<< "[MEASU_INIT] besselBCoeff is of size "<<besselBCoeff.rows()<<"x"<<besselBCoeff.cols()<<"\n";
    std::cout<< "[MEASU_INIT] besselBCoeff:\n"<<besselBCoeff<<"\n";
    std::cout<< "[MEASU_INIT] besselACoeff is of size "<<besselACoeff.rows()<<"x"<<besselACoeff.cols()<<"\n";
    std::cout<< "[MEASU_INIT] besselACoeff:\n"<<besselACoeff<<"\n";
    #endif
}


double Measurement::navigationKinematics(const Eigen::Matrix< double, Eigen::Dynamic, Eigen::Dynamic > &A, const Eigen::Matrix< double, Eigen::Dynamic, Eigen::Dynamic > &B,
					const Eigen::Matrix< double, Eigen::Dynamic, Eigen::Dynamic > &R, const Eigen::Matrix< double, Eigen::Dynamic, Eigen::Dynamic > &W)
{
    double leastSquaresError = std::numeric_limits<double>::quiet_NaN();
    DataModel velM = DataModel(NUMAXIS);
    DataModel ivelM = DataModel(NUMAXIS);
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
    y.block<NUMAXIS, 1> (0,0) = this->angvelo.front();
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
    
    /** Get the current velocity model from the LS solution **/
    velocity = x.block<NUMAXIS,1>(0,0);
    
    /** FILTERING THE VELOCITY **/
    /** Copy to the array for filtering **/
    arrayVelMM.block<NUMAXIS, NORDER_BESSEL_FILTER>(0,0) = arrayVelMM.block<NUMAXIS, NORDER_BESSEL_FILTER>(0,1);
    arrayVelMM.col(localization::NORDER_BESSEL_FILTER) = velocity;
    
    /** IIR Low-pass Bessel Filter **/
    velocity = Measurement::iirFilter(localization::NORDER_BESSEL_FILTER, besselBCoeff, besselACoeff, arrayVelMM, arrayVelMMIIR);
    /** END OF FILTERING THE VELOCITY **/
    
    /** Store in a DataModel variable **/
    velM.data = velocity;
    velM.Cov = pow(leastSquaresError,2) * Eigen::Matrix<double, NUMAXIS, NUMAXIS>::Identity() +  Iinverse.block<NUMAXIS, NUMAXIS> (0,0);
    
    /** INCREMENT IN VELOCITY **/
    /** Check the buffer size of the nav. kinematics velocities **/
    if (increVelMM.size() >= this->ibuffer_size)
	increVelMM.pop_back();
    
    /** Store the previous velocity in the local variable **/
    if (velMM.empty())
	ivelM = velM; //! If empty the increment is the actual velocity
    else
	ivelM = velM - velMM.front();
    
    /** Push the increment in velocities in the array of increments vel. Model **/
    increVelMM.push_front(ivelM);
    
    /** VELOCITY ARRAY **/
    /** Check the buffer size of the nav. kinematics velocities **/
    if (velMM.size() >= this->ibuffer_size)
	velMM.pop_back();
    
    /** Push the DataModel in the array of velocities **/
    velMM.push_front(velM);
    
    /** FILTER ARRANGEMENT FOR NEXT ITERATION **/
    /** Move back one column for the next bessel filter iteration **/
    arrayVelMMIIR.block<NUMAXIS, NORDER_BESSEL_FILTER>(0,0) = arrayVelMMIIR.block<NUMAXIS, NORDER_BESSEL_FILTER>(0,1);
   
    /** SLIP VECTOR (Z COMPONENT) **/
    /** Save the z-axis slip vector results **/
    slip.setZero();
    slip(2) = x(3); slip(5) = x(4); slip(8) = x(5); slip(11) = x(6);
    slipModel.data = slip;
    slipModel.Cov(2,2) = Iinverse(3,3);
    slipModel.Cov(5,5) = Iinverse(4,4);
    slipModel.Cov(8,8) = Iinverse(5,5);
    slipModel.Cov(11,11) = Iinverse(6,6);
    
    
    /** DELTA CHANGE OF THE CONTACT ANGLE **/
    /** Store in the dcontactAngle variable **/
    dcontactAngle.data = x.block<NUMBER_OF_WHEELS,1>(NUMAXIS+NUMBER_OF_WHEELS,0);
    dcontactAngle.Cov = Iinverse.bottomRightCorner<NUMBER_OF_WHEELS,NUMBER_OF_WHEELS>();
    
    /** Error in the least-Squares **/
    this->leastSquaredNavigation = leastSquaresError;
    
     #ifdef DEBUG_PRINTS
    std::cout << "[NK] velMM has size "<<velMM.size()<<"\n";
    std::cout << "[NK] velModel.data:\n" << velM.data << std::endl;
    std::cout << "[NK] velModel.Cov:\n" << velM.Cov << std::endl;
    std::cout << "[NK] increVelMM has size "<<increVelMM.size()<<"\n";
    std::cout << "[NK] increVelModel.data:\n" << increVelMM.front().data << std::endl;
    std::cout << "[NK] increVelModel.Cov:\n" << increVelMM.front().Cov << std::endl;
    std::cout << "[NK] ** IIR Filter Variables **"<< std::endl;
    std::cout << "[NK] arrayVelMM is of size "<< arrayVelMM.rows()<<"x"<<arrayVelMM.cols()<<"\n";
    std::cout << "[NK] The arrayVelMM is:\n" << arrayVelMM << std::endl;
    std::cout << "[NK] arrayVelMMIIR is of size "<< arrayVelMMIIR.rows()<<"x"<<arrayVelMMIIR.cols()<<"\n";
    std::cout << "[NK] The arrayVelMMIIR is:\n" << arrayVelMMIIR << std::endl;
    #endif
    
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
//     y.block<NUMAXIS, 1> (0,0) = linvelocity;
    y.block<NUMAXIS, 1> (NUMAXIS,0) = this->angvelo.front();
    y.block<ENCODERS_VECTOR_SIZE, 1> (2*NUMAXIS,0) = encodersvelocity;
    y.block<NUMBER_OF_WHEELS, 1> ((2*NUMAXIS) + ENCODERS_VECTOR_SIZE,0) = dcontactAngle.data;
    

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
* @brief This perform simpson rule for integration
*/
double Measurement::simpsonsIntegral(double fa, double fm, double fb, double delta_ab)
{
    return delta_ab/6.0 * (fa + (4.0*fm) + fb);
}


/**
* @brief Save slip info to the Orogen-compatible DataType
*/
void Measurement::toSlipInfo(localization::SlipInfo &sinfo)
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

