#ifndef _PROCESS_MODELS_HPP_
#define _PROCESS_MODELS_HPP_


#define PROCESS_MODEL_DEBUG_PRINTS 1

namespace localization
{
//    template<typename _SingleState>
//    _SingleState processModel (const _SingleState &serror, const Eigen::Matrix<double, 3, 1> &acc,
//                    const Eigen::Matrix<double, 3, 1> &angvelo,
//                    const Eigen::Quaternion<double> &orientq, double dt)
//    {
//        _SingleState s2; /** Propagated state */
//        Eigen::Matrix <double,3,3> dFki; /** Discrete form of subsytem models */
//
//        //std::cout<<"[PROCESS_MODEL] acc:\n"<<acc<<"\n";
//        //std::cout<<"[PROCESS_MODEL] angevelo:\n"<<angvelo<<"\n";
//
//        /** Propagate the position (discretization approx) **/
//        dFki = Eigen::Matrix<double, 3, 3>::Identity() * dt +
//            0.5 * Eigen::Matrix<double, 3, 3>::Identity() * Eigen::Matrix<double, 3, 3>::Identity() * pow(dt,2);
//        s2.pos = serror.pos + dFki * serror.vel;
//
//        /** Propagate the velocity (discretization approx) **/
//        Eigen::Matrix <double,3, 1> deltaVel;
//        deltaVel = /*orientq.inverse() */ (acc * dt); /** TO-DO: Increment in velocity in world frame */
//        s2.vel = serror.vel +  serror.orient * deltaVel - orientq.inverse() * (serror.abias * dt);
//
//        /** Propagate the error quaternion **/
//        Eigen::Matrix <double,3, 1> deltaAngle = (angvelo - 0.5 * serror.gbias) * dt;
//        MTK::SO3 <double> rot = MTK::SO3<double>::exp (deltaAngle);
//        s2.orient = (serror.orient * rot);
//
//        /** The bias (gyros and acc) update **/
//        s2.gbias = serror.gbias;
//        s2.abias = serror.abias;
//
//        return s2;
//    }

    /**@brief ProcessModel assuming corrected values
     */
    template <typename _SingleState>
    _SingleState processModel (const _SingleState &serror, const Eigen::Matrix<double, 3, 1> &deltaVel, const Eigen::Quaterniond &delta_orient,
                            const Eigen::Quaterniond &orient, double dt)
    {
        _SingleState s2; /** Propagated state */
        Eigen::Matrix <double,3,3> dFki; /** Discrete form of subsystem models */

        //std::cout<<"[PROCESS_MODEL] acc:\n"<<acc<<"\n";
        //std::cout<<"[PROCESS_MODEL] angevelo:\n"<<angvelo<<"\n";

        /** Propagate the error position (discretization approx) **/
        dFki = Eigen::Matrix<double, 3, 3>::Identity() * dt +
            0.5 * Eigen::Matrix<double, 3, 3>::Identity() * Eigen::Matrix<double, 3, 3>::Identity() * pow(dt,2);

        s2.pos = serror.pos + dFki * serror.vel; //!Position (error velocity comes in world frame)

        /** Propagate the error velocity (discretization approx) **/
        s2.vel = serror.vel; // serror orient is zero +  serror.orient * deltaVel; //! Error velocity in world frame

        /** Propagate the error quaternion **/
        //s2.orient = (serror.orient * delta_orient);

        /** The bias (gyros and acc) update **/
        s2.gbias = serror.gbias;
        s2.abias = serror.abias;

        return s2;
    }

    template <typename _SingleState, typename _ProcessModelMatrix>
    _ProcessModelMatrix processModel (const Eigen::Vector3d &angvelo, const Eigen::Vector3d &acc,
                                    const Eigen::Quaterniond &orient, double dt)
    {
        _ProcessModelMatrix F, dF; /** Process matrix and discrete process Matrix */
        Eigen::Matrix3d velo2product; /** Vec 2 product  matrix */
        Eigen::Matrix3d acc2product; /** Vec 2 product  matrix */

        /** In cross product form **/
        velo2product << 0, -angvelo(2), angvelo(1),
                     angvelo(2), 0, -angvelo(0),
                     -angvelo(1), angvelo(0), 0;

        /** In cross product form **/
        acc2product << 0, -acc(2), acc(1),
                    acc(2), 0, -acc(0),
                    -acc(1), acc(0), 0;

        /** Form the process model matrix **/
        F.setZero();dF.setZero();
        F.template block<3,3>(0,3) = orient.matrix(); //! Position part (velocity is comming in body, position in world)
        F.template block<3,3>(3,6) = -acc2product; //! Velocity part
        F.template block<3,3>(3,3*4) = -Eigen::Matrix3d::Identity(); //! Velocity part
        F.template block<3,3>(6,6) = -velo2product; //! Attitude part
        F.template block<3,3>(6,3*3) = -0.5*Eigen::Matrix3d::Identity(); //! Attitude part

        /** Discrete the matrix  (ZOH) **/
        dF = Eigen::Matrix<double, _SingleState::DOF, _SingleState::DOF>::Identity() + F * dt + F * F * pow(dt,2)/2.0;

        #ifdef PROCESS_MODEL_DEBUG_PRINTS
        std::cout<<"[PROCESS_MODEL] F is\n"<<F<<"\n";
        std::cout<<"[PROCESS_MODEL] dF is\n"<<dF<<"\n";
        #endif
        return dF;
    }


    template <typename _SingleState, typename _SingleStateCovariance>
    _SingleStateCovariance processNoiseCov (const Eigen::Matrix<double, _SingleState::DOF, _SingleState::DOF> &F,
                    const base::Vector3d &accrw,
                    const base::Vector3d &gyrorw,
                    const base::Vector3d &gbiasins,
                    const base::Vector3d &abiasins,
                    const Eigen::Quaterniond &orient, double &delta_t)
    {
        //TO-DO Eigen::Matrix<double, 3, 3> Cq;
        double sqrtdelta_t = sqrt(delta_t); /** Square root of delta time interval */

        /** Dimension is for one single error state **/
        _SingleStateCovariance Cov = _SingleStateCovariance::Zero();// Covariance matrix
        _SingleStateCovariance dCov = _SingleStateCovariance::Zero();// Discrete Cov matrix

        /** Noise for error in velocity **/
        Eigen::Matrix3d Qa;
        Qa.setZero();
        Qa(0,0) = pow(accrw[0]/sqrtdelta_t,2);
        Qa(1,1) = pow(accrw[1]/sqrtdelta_t,2);
        Qa(2,2) = pow(accrw[2]/sqrtdelta_t,2);
        ::MTK::subblock (Cov, &_SingleState::vel) = /*Cq.inverse */ Qa;

        /** Noise for error in position **/
        ::MTK::subblock (Cov, &_SingleState::pos) = orient.matrix() * Qa * delta_t;//Eigen::Matrix3d::Zero();

        /** Noise for error in orientation **/
        Eigen::Matrix3d Qg;
        Qg.setZero();
        Qg(0,0) = pow(gyrorw[0]/sqrtdelta_t,2);
        Qg(1,1) = pow(gyrorw[1]/sqrtdelta_t,2);
        Qg(2,2) = pow(gyrorw[2]/sqrtdelta_t,2);
        ::MTK::subblock (Cov, &_SingleState::orient) = 0.25 * Qg;

        /** Noise for error in gyros bias instability **/
        Eigen::Matrix3d Qgbias;
        Qgbias.setZero();
        Qgbias(0,0) = pow(gbiasins[0],2);
        Qgbias(1,1) = pow(gbiasins[1],2);
        Qgbias(2,2) = pow(gbiasins[2],2);
        ::MTK::subblock (Cov, &_SingleState::gbias) = Qgbias;

        /** Noise for error in acc bias instability **/
        Eigen::Matrix3d Qabias;
        Qabias.setZero();
        Qabias(0,0) = pow(abiasins[0],2);
        Qabias(1,1) = pow(abiasins[1],2);
        Qabias(2,2) = pow(abiasins[2],2);
        ::MTK::subblock (Cov, &_SingleState::abias) = Qabias;

        /** Form the system noise matrix (discretization) **/
        dCov = Cov*delta_t + 0.5*delta_t*F*Cov + 0.5*delta_t*Cov*F.transpose();
        dCov = 0.5*(dCov + dCov.transpose());

        #ifdef PROCESS_MODEL_DEBUG_PRINTS
        std::cout<<"[PROCESS_NOISE_COV] Cov is\n"<<Cov<<"\n";
        std::cout<<"[PROCESS_NOISE_COV] dCov is\n"<<dCov<<"\n";
        #endif

        return dCov ;
    }
}

#endif /** End of ProcessModels.hpp **/

