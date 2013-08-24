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


    template <typename _SingleState, typename _SingleStateCovariance>
    _SingleStateCovariance processNoiseCov (const base::Vector3d &accrw,
                    const base::Vector3d &gyrorw,
                    const base::Vector3d &gbiasins,
                    const base::Vector3d &abiasins,
                    const Eigen::Quaterniond &robotorient, double &delta_t)
    {
        //TO-DO Eigen::Matrix<double, 3, 3> Cq;
        double sqrtdelta_t = sqrt(delta_t); /** Square root of delta time interval */

        /** Dimension is for one single error state **/
        _SingleStateCovariance cov = _SingleStateCovariance::Zero();

        /** Noise for error in position **/
        ::MTK::subblock (cov, &_SingleState::pos) = Eigen::Matrix3d::Zero();

        /** Noise for error in velocity **/
        Eigen::Matrix3d Qa;
        Qa.setZero();
        Qa(0,0) = pow(accrw[0]/sqrtdelta_t,2);
        Qa(1,1) = pow(accrw[1]/sqrtdelta_t,2);
        Qa(2,2) = pow(accrw[2]/sqrtdelta_t,2);
        ::MTK::subblock (cov, &_SingleState::vel) = /*Cq.inverse */ Qa;

        /** Noise for error in orientation **/
        Eigen::Matrix3d Qg;
        Qg.setZero();
        Qg(0,0) = pow(gyrorw[0]/sqrtdelta_t,2);
        Qg(1,1) = pow(gyrorw[1]/sqrtdelta_t,2);
        Qg(2,2) = pow(gyrorw[2]/sqrtdelta_t,2);
        ::MTK::subblock (cov, &_SingleState::orient) = 0.25 * Qg;

        /** Noise for error in gyros bias instability **/
        Eigen::Matrix3d Qgbias;
        Qgbias.setZero();
        Qgbias(0,0) = pow(gbiasins[0],2);
        Qgbias(1,1) = pow(gbiasins[1],2);
        Qgbias(2,2) = pow(gbiasins[2],2);
        ::MTK::subblock (cov, &_SingleState::gbias) = Qgbias;

        /** Noise for error in acc bias instability **/
        Eigen::Matrix3d Qabias;
        Qabias.setZero();
        Qabias(0,0) = pow(abiasins[0],2);
        Qabias(1,1) = pow(abiasins[1],2);
        Qabias(2,2) = pow(abiasins[2],2);
        ::MTK::subblock (cov, &_SingleState::abias) = Qabias;

        return cov ;
    }
}

#endif /** End of ProcessModels.hpp **/

