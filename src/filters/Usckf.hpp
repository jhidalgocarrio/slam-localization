#ifndef _USCKF_HPP_
#define _USCKF_HPP_

/** Standard libraries **/
#include <vector> /** std::vector */
#include <algorithm> /** std::transform */
#include <numeric>

/** Boost **/
#include <boost/bind.hpp>

/** Eigen **/
#include <Eigen/Core>
#include <Eigen/LU>
#include <Eigen/QR>
#if EIGEN_VERSION_AT_LEAST(3,0,0)
    #include <Eigen/Eigenvalues>
#endif

#include <Eigen/Cholesky>
#include <Eigen/SVD> /** Singular Value Decomposition (SVD) of Eigen */

/** UKFOM library **/
#include <ukfom/lapack/cholesky.hpp>
#include <ukfom/traits/dof.hpp>
#include <ukfom/util.hpp>

// MTK's pose and orientation definition:
#include <mtk/startIdx.hpp>


#define USCKF_DEBUG_PRINTS 1

namespace localization
{

    /** Different cloning mechanism **/
    enum CloningMode
    {
        STATEK_I = 1,
        STATEK_L = 2,
        STATEK = 3
    };

    template <typename _AugmentedState, typename _SingleState>
    class Usckf
    {
        typedef Usckf self;

        public:
            enum
            {
                    DOF_AUGMENTED_STATE = _AugmentedState::DOF
            };
            enum
            {
                    DOF_SINGLE_STATE = _SingleState::DOF
            };


            typedef typename _AugmentedState::scalar_type ScalarType;
            typedef typename _AugmentedState::vectorized_type VectorizedAugmentedState;
            typedef typename _SingleState::vectorized_type VectorizedSingleState;
            typedef Eigen::Matrix<ScalarType, int(_AugmentedState::DOF), int(_AugmentedState::DOF)> AugmentedStateCovariance;
            typedef Eigen::Matrix<ScalarType, int(_SingleState::DOF), int(_SingleState::DOF)> SingleStateCovariance;
            typedef std::vector<_AugmentedState> AugmentedStateSigma;
            typedef std::vector<_SingleState> SingleStateSigma;


        private:

            _AugmentedState mu_state; /** Mean of the state vector **/
            AugmentedStateCovariance Pk; /** Covariance of the State vector **/

        public:
            /**@brief Constructor
             */
            Usckf(const _AugmentedState &state, const AugmentedStateCovariance &P0)
                : mu_state(state), Pk(P0)
            {
            }

            /**@brief Filter prediction step
             */
            template<typename _ProcessModel>
            void predict(_ProcessModel f, const SingleStateCovariance &Q)
            {
                predict(f, boost::bind(ukfom::id<SingleStateCovariance>, Q));
            }

            template<typename _ProcessModel, typename _ProcessNoiseCovariance>
            void predict(_ProcessModel f, _ProcessNoiseCovariance Q)
            {
                /** Get the current state vector to propagate **/
                _SingleState statek_i = mu_state.statek_i;
                SingleStateCovariance Pk_i = MTK::subblock (Pk, &_AugmentedState::statek_i);

                #ifdef USCKF_DEBUG_PRINTS
                std::cout<<"[USCKF_PREDICT] statek_i(k|k):\n"<<statek_i<<"\n";
                std::cout<<"[USCKF_PREDICT] Pk_i(k|k):\n"<<Pk_i<<"\n";
                #endif

                /** Propagation only uses the current state (DOF_SINGLE_STATE dimension ) **/
                SingleStateSigma X(2 * DOF_SINGLE_STATE + 1);

                /** Generates the sigma Points of a Single State **/
                generateSigmaPoints(statek_i, Pk_i, X);

                /** Create a copy before the transformation **/
                SingleStateSigma XCopy(2 * DOF_SINGLE_STATE + 1);
                XCopy = X;

                /*****************************/
                /** Process Model Transform **/
                /*****************************/

                /** Apply the non-linear transformation of the process model **/
                std::transform(X.begin(), X.end(), X.begin(), f);

                #ifdef USCKF_DEBUG_PRINTS
                printSigmaPoints<SingleStateSigma>(X);
                #endif

                /** Compute the mean **/
                mu_state.statek_i = meanSigmaPoints(X);

                /** Compute the cross-covariance matrix (cross because is between X and X before the transform) **/
                SingleStateCovariance Pxy;
                Pxy = crossCovSigmaPoints<_SingleState, _SingleState::DOF, SingleStateSigma,
                    _SingleState> (statek_i, mu_state.statek_i, XCopy, X);
                SingleStateCovariance Fk = Pxy.transpose()*Pk_i.inverse();

                #ifdef USCKF_DEBUG_PRINTS
                std::cout<<"[USCKF_PREDICT] Fk:\n"<< Fk << std::endl;
                #endif

                /*************************/
                /** Discretization of Q **/
                /*************************/

                SingleStateCovariance Qk = Q();

                /** Just for debugging purpose print **/
                /** Pyy =  Fk*Pk_i*Fk^T and  Pk_i = Pyy+Q **/
                /** After propagation/prediction, which means Pyy+Q = P(k+1|k) **/
                #ifdef USCKF_DEBUG_PRINTS
                std::cout<<"[USCKF_PREDICT] Fk*Pk_i*Fk^T + Qk:\n"<< Fk*Pk_i*Fk.transpose() + Qk<< std::endl;
                #endif

                /************************/
                /** Covariance Matrix  **/
                /************************/

                /** Compute the Process model Covariance **/
                Pk_i = covSigmaPoints<_SingleState::DOF, _SingleState>(mu_state.statek_i, X) + Qk;

                /** Store the subcovariance matrix for statek_i **/
                MTK::subblock (Pk, &_AugmentedState::statek_i) = Pk_i;

                /** Compute the Cross Cov for the Copy States of the AugmentedState **/
                SingleStateCovariance Pkkk;

                /************************/
                /**  Cross-Cov Matrix  **/
                /************************/

                /** Covariance between state(k) and state(k+i) **/
                Pkkk = MTK::subblock (Pk, &_AugmentedState::statek, &_AugmentedState::statek_i);
                Pkkk = Pkkk * Fk.transpose();
                MTK::subblock (Pk, &_AugmentedState::statek, &_AugmentedState::statek_i) = Pkkk;

                /** Covariance between state(k+l) and state(k+i) **/
                Pkkk = MTK::subblock (Pk, &_AugmentedState::statek_l, &_AugmentedState::statek_i);
                Pkkk = Pkkk * Fk.transpose();
                MTK::subblock (Pk, &_AugmentedState::statek_l, &_AugmentedState::statek_i) = Pkkk;

                /** Covariance between state(k+i) and state(k) **/
                Pkkk = MTK::subblock (Pk, &_AugmentedState::statek_i, &_AugmentedState::statek);
                Pkkk = Fk * Pkkk;
                MTK::subblock (Pk, &_AugmentedState::statek_i, &_AugmentedState::statek) = Pkkk;

                /** Covariance between state(k+i) and state(k+l) **/
                Pkkk = MTK::subblock (Pk, &_AugmentedState::statek_i, &_AugmentedState::statek_l);
                Pkkk = Fk * Pkkk;
                MTK::subblock (Pk, &_AugmentedState::statek_i, &_AugmentedState::statek_l) = Pkkk;

                /**********/
                /** TO-DO: cross-cov with the features (Dynamic size part of the filter) **/

                #ifdef  USCKF_DEBUG_PRINTS
                std::cout << "[USCKF_PREDICT] statek_i(k+1|k):" << std::endl << mu_state.statek_i << std::endl;
                std::cout << "[USCKF_PREDICT] Pk(k+1|k):"<< std::endl << Pk_i << std::endl;
                std::cout << "[USCKF_PREDICT] Process Noise Cov Q(k):"<< std::endl << Q() << std::endl;
                #endif
            }

            template<typename _Measurement, typename _MeasurementModel, typename _MeasurementNoiseCovariance>
            void update(const _Measurement &z, _MeasurementModel h, _MeasurementNoiseCovariance R)
            {
                    update(z, h, R, ukfom::accept_any_mahalanobis_distance<ScalarType>);
            }

            template<typename _Measurement, typename _MeasurementModel>
            void update(const _Measurement &z, _MeasurementModel h,
                        const Eigen::Matrix<ScalarType, ukfom::dof<_Measurement>::value, ukfom::dof<_Measurement>::value> &R)
            {
                    typedef Eigen::Matrix<ScalarType, ukfom::dof<_Measurement>::value, ukfom::dof<_Measurement>::value> measurement_cov;
                    update(z, h, boost::bind(ukfom::id<measurement_cov>, R), ukfom::accept_any_mahalanobis_distance<ScalarType>);
            }

            template<typename _Measurement, typename _MeasurementModel,
                    typename _MeasurementNoiseCovariance, typename _SignificanceTest>
            void update(const _Measurement &z, _MeasurementModel h,
                        _MeasurementNoiseCovariance R, _SignificanceTest mt)
            {
                    const static int measurement_rows = ukfom::dof<_Measurement>::value;
                    typedef _Measurement Measurement;
                    typedef Eigen::Matrix<ScalarType, measurement_rows, 1> VectorizedMeasurement;
                    typedef std::vector<Measurement> measurement_vector;
                    typedef Eigen::Matrix<ScalarType, measurement_rows, measurement_rows> MeasurementCov;
                    typedef Eigen::Matrix<ScalarType, _AugmentedState::DOF, measurement_rows> CrossCov;

                    AugmentedStateSigma X(2 * DOF_AUGMENTED_STATE + 1);
                    generateSigmaPoints(mu_state, Pk, X);

                    std::vector<Measurement> Z(X.size());
                    std::transform(X.begin(), X.end(), Z.begin(), h);

                    const Measurement meanZ = meanSigmaPoints(Z);
                    const MeasurementCov S = covSigmaPoints<measurement_rows>(meanZ, Z) + R();
                    const CrossCov covXZ = crossCovSigmaPoints<measurement_rows>(mu_state, meanZ, X, Z);

                    MeasurementCov S_inverse;
                    S_inverse = S.inverse();

                    const CrossCov K = covXZ * S_inverse;

                    const VectorizedMeasurement innovation = z - meanZ;

                    const ScalarType mahalanobis2 = (innovation.transpose() * S_inverse * innovation)(0);

                    if (mt(mahalanobis2))
                    {
                            Pk -= K * S * K.transpose();
                            //applyDelta(K * innovation);
                            mu_state = mu_state + K * innovation;
                    }

                    #ifdef USCKF_DEBUG_PRINTS
                    std::cout << "[USCKF_UPDATE] innovation:" << std::endl << innovation << std::endl;
                    std::cout << "[USCKF_UPDATE] mu_state':" << std::endl << mu_state << std::endl;
                    #endif
            }

            template<typename _Measurement, typename _MeasurementModel,
                     typename _MeasurementModelCovariance, typename _MeasurementCrossCovariance,
                     typename _MeasurementNoiseCovariance>
            void delayUpdate(const _Measurement &z, _MeasurementModel h,
                    _MeasurementModelCovariance Pzz,
                    _MeasurementCrossCovariance Pxz,
                    _MeasurementNoiseCovariance R)
            {
                    delayUpdate(z, h, Pzz, Pxz, R, accept_mahalanobis_distance<ScalarType>);//ukfom::accept_any_mahalanobis_distance<ScalarType>);
            }

            template<typename _Measurement, typename _MeasurementModel>
            void delayUpdate(const _Measurement &z, _MeasurementModel h,
                        const Eigen::Matrix<ScalarType, ukfom::dof<_Measurement>::value, ukfom::dof<_Measurement>::value> &Pzz,
                        const Eigen::Matrix<ScalarType, ukfom::dof<_AugmentedState>::value, ukfom::dof<_Measurement>::value> &Pxz,
                        const Eigen::Matrix<ScalarType, ukfom::dof<_Measurement>::value, ukfom::dof<_Measurement>::value> &R)
            {
                    typedef Eigen::Matrix<ScalarType, ukfom::dof<_AugmentedState>::value, ukfom::dof<_Measurement>::value> MeasurementCrossCov;
                    typedef Eigen::Matrix<ScalarType, ukfom::dof<_Measurement>::value, ukfom::dof<_Measurement>::value> MeasurementCov;
                    delayUpdate(z, h, boost::bind(ukfom::id<MeasurementCov>, Pzz),
                                 boost::bind(ukfom::id<MeasurementCrossCov>, Pxz),
                                 boost::bind(ukfom::id<MeasurementCov>, R),
                                 accept_mahalanobis_distance<ScalarType>);
                                 //ukfom::accept_any_mahalanobis_distance<ScalarType>);
            }

            template<typename _Measurement, typename _MeasurementModel,
                     typename _MeasurementModelCovariance, typename _MeasurementCrossCovariance,
                     typename _MeasurementNoiseCovariance, typename _SignificanceTest>
            void delayUpdate(const _Measurement &z, _MeasurementModel h,
                        _MeasurementModelCovariance Pzz,
                        _MeasurementCrossCovariance Pxz,
                        _MeasurementNoiseCovariance R,
                        _SignificanceTest mt)
            {
                const static int DOF_MEASUREMENT = ukfom::dof<_Measurement>::value;
                typedef _Measurement Measurement;
                typedef Eigen::Matrix<ScalarType, DOF_MEASUREMENT, 1> VectorizedMeasurement;
                typedef std::vector<Measurement> measurement_vector;
                typedef Eigen::Matrix<ScalarType, DOF_MEASUREMENT, DOF_MEASUREMENT> MeasurementCov;
                typedef Eigen::Matrix<ScalarType, _AugmentedState::DOF, DOF_MEASUREMENT> CrossCov;

                const Measurement meanZ = h(mu_state);
                const MeasurementCov S = Pzz() + R(); //Measurement Model Covariance + Measurement Noise Covariance
                const CrossCov covXZ = Pxz();

                MeasurementCov S_inverse;
                S_inverse = S.inverse();

                const CrossCov K = covXZ * S_inverse;

                const VectorizedMeasurement innovation = z - meanZ;

                const ScalarType mahalanobis2 = (innovation.transpose() * S_inverse * innovation)(0);
                #ifdef USCKF_DEBUG_PRINTS
                std::cout<<"[DELAY_UPDATE] mu_state(k+1|k):\n"<<mu_state<<"\n";
                std::cout<<"[DELAY_UPDATE] Pk(k+1|k):\n"<<Pk<<"\n";
                std::cout<<"[DELAY_UPDATE] K is of size:" <<K.rows()<<"x"<<K.cols()<<std::endl;
                std::cout<<"[DELAY_UPDATE] K:\n"<<K<<"\n";
                std::cout<<"[DELAY_UPDATE] mahalanobis2:" << std::endl << mahalanobis2 << std::endl;
                #endif

                if (mt(mahalanobis2, innovation.size()-1))
                {
                    std::cout<<"[DELAY_UPDATE]: OK - FUSION\n";

                    Pk -= K * S * K.transpose();
                    mu_state = mu_state + K * innovation;
                }
                else
                {
                    std::cout<<"[DELAY_UPDATE]: FAILED - DEFUSION\n";
                    Pk += K * (innovation * innovation.transpose()) * K.transpose();
                }

                #ifdef USCKF_DEBUG_PRINTS
                std::cout << "[DELAY_UPDATE] mu_state(K+1|K+1):" << std::endl << mu_state << std::endl;
                std::cout << "[DELAY_UPDATE] Pk(k+1|k+1):\n"<<Pk<<"\n";
                std::cout << "[DELAY_UPDATE] innovation:" << std::endl << innovation << std::endl;
                std::cout << "[DELAY_UPDATE] R is of size:" <<R().rows()<<"x"<<R().cols()<<std::endl;
                std::cout << "[DELAY_UPDATE] R:\n" << R() <<std::endl;
                std::cout << "[DELAY_UPDATE] Pzz is of size:" <<Pzz().rows()<<"x"<<Pzz().cols()<<std::endl;
                std::cout << "[DELAY_UPDATE] Pzz:\n" << Pzz() <<std::endl;
                #endif

            }


            template<typename _Measurement, typename _MeasurementModel, typename _MeasurementNoiseCovariance>
            void singleUpdate(const _Measurement &z, _MeasurementModel h, _MeasurementNoiseCovariance R)
            {
                    singleUpdate(z, h, R, ukfom::accept_any_mahalanobis_distance<ScalarType>);
            }

            template<typename _Measurement, typename _MeasurementModel>
            void singleUpdate(const _Measurement &z, _MeasurementModel h,
                        const Eigen::Matrix<ScalarType, ukfom::dof<_Measurement>::value, ukfom::dof<_Measurement>::value> &R, int order = STATEK_I)
            {
                    typedef Eigen::Matrix<ScalarType, ukfom::dof<_Measurement>::value, ukfom::dof<_Measurement>::value> MeasurementCov;
                    singleUpdate(z, h, boost::bind(ukfom::id<MeasurementCov>, R), ukfom::accept_any_mahalanobis_distance<ScalarType>, order);
            }


            template<typename _Measurement, typename _MeasurementModel,
                    typename _MeasurementNoiseCovariance, typename _SignificanceTest>
            void singleUpdate(const _Measurement &z, _MeasurementModel h,
                        _MeasurementNoiseCovariance R, _SignificanceTest mt, int order = STATEK_I)
            {
                    const static int DOF_MEASUREMENT = ukfom::dof<_Measurement>::value;
                    typedef _Measurement Measurement;
                    typedef Eigen::Matrix<ScalarType, DOF_MEASUREMENT, 1> VectorizedMeasurement;
                    typedef std::vector<Measurement> measurement_vector;
                    typedef Eigen::Matrix<ScalarType, DOF_MEASUREMENT, DOF_MEASUREMENT> MeasurementCov;
                    typedef Eigen::Matrix<ScalarType, _SingleState::DOF, DOF_MEASUREMENT> CrossCov;

                    /** Get the current state vector to propagate **/
                    _SingleState statek_i = muSingleState(order);
                    SingleStateCovariance Pk_i = PkSingleState(order);

                    #ifdef USCKF_DEBUG_PRINTS
                    std::cout<<"[SINGLE_UPDATE] statek_i(k+1|k):\n"<<statek_i<<"\n";
                    std::cout<<"[SINGLE_UPDATE] Pk_i(k+1|k):\n"<<Pk_i<<"\n";
                    #endif

                    SingleStateSigma X(2 * DOF_SINGLE_STATE + 1);
                    generateSigmaPoints(statek_i, Pk_i, X);

                    std::vector<Measurement> Z(X.size());
                    std::transform(X.begin(), X.end(), Z.begin(), h);

                    const Measurement meanZ = meanSigmaPoints(Z);
                    const MeasurementCov S = covSigmaPoints<DOF_MEASUREMENT, Measurement>(meanZ, Z) + R();
                    const CrossCov covXZ = crossCovSigmaPoints<_SingleState, DOF_MEASUREMENT,
                                            SingleStateSigma, Measurement>(statek_i, meanZ, X, Z);

                    MeasurementCov S_inverse;
                    S_inverse = S.inverse();

                    const CrossCov K = covXZ * S_inverse;

                    const VectorizedMeasurement innovation = z - meanZ;

                    const ScalarType mahalanobis2 = (innovation.transpose() * S_inverse * innovation)(0);

                    if (mt(mahalanobis2))
                    {
                            Pk_i -= K * S * K.transpose();//temporary covariance still relative to the old mean
                            applyDeltaSingleState(statek_i, Pk_i, K * innovation);
                    }

                    /** Store the sub-covariance matrix for statek_i **/
                    setPkSingleState(Pk_i, order);

                    /** Update the state in the general variable **/
                    setSingleState(statek_i, order);

                    #ifdef USCKF_DEBUG_PRINTS
                    std::cout <<"[SINGLE_UPDATE] innovation:" << std::endl << innovation << std::endl;
                    std::cout <<"[SINGLE_UPDATE] statek_i(k+1|k+1):" << std::endl << statek_i<< std::endl;
                    std::cout <<"[SINGLE_UPDATE] Pk_i(k+1|k+1):\n"<<Pk_i<<"\n";
                    std::cout <<"[SINGLE_UPDATE] R is of size:" <<R().rows()<<"x"<<R().cols()<<std::endl;
                    std::cout <<"[SINGLE_UPDATE] R:\n" << R() <<std::endl;
                    #endif
            }


            /** @brief Standard EKF Update of the state
             */
            template <typename _Measurement, typename _MeasurementModel, typename _MeasurementNoiseCovariance>
            void indirectUpdate(const _Measurement &z, _MeasurementModel H, _MeasurementNoiseCovariance R)
            {
                indirectUpdate(z, H, R, ukfom::accept_any_mahalanobis_distance<ScalarType>);

            }

            /** @brief Indirect Single Update of the state.
             */
            template <typename _Measurement, typename _MeasurementModel,
                    typename _MeasurementNoiseCovariance,  typename _SignificanceTest>
            void indirectUpdate(const _Measurement &z, _MeasurementModel H,
                    _MeasurementNoiseCovariance R, _SignificanceTest mt)
            {
                const static int DOF_MEASUREMENT = ukfom::dof<_Measurement>::value; /** Dimension of the measurement */

                /** Get the state in vector form **/
                VectorizedSingleState xk_i; xk_i.setZero();

                /** statek_i covariance matrix **/
                SingleStateCovariance Pk_i = MTK::subblock (Pk, &_AugmentedState::statek_i);

                #ifdef USCKF_DEBUG_PRINTS
                std::cout << "[INDIRECT_UPDATE] x_error(before):\n" << xk_i <<std::endl;
                std::cout << "[INDIRECT_UPDATE] Pk_i(before):\n" << Pk_i <<std::endl;
                #endif

                /** Compute the Kalman Gain Matrix **/
                Eigen::Matrix<ScalarType, DOF_MEASUREMENT, DOF_MEASUREMENT> S, S_inverse;
                Eigen::Matrix<ScalarType, DOF_SINGLE_STATE, DOF_MEASUREMENT> K;
                S = H * Pk_i* H.transpose() + R; //!Calculate the covariance of the innovation
                S_inverse = S.inverse();
                K = Pk_i * H.transpose() * S_inverse; //!Calculate K using the inverse of S

                /** Innovation **/
                const _Measurement innovation = (z - H * xk_i);
                const ScalarType mahalanobis2 = (innovation.transpose() *  S_inverse * innovation)(0);

                /** Update the state vector and the covariance matrix */
                if (mt(mahalanobis2))
                {

                    /** Update the state vector and the covariance matrix */
                    #ifdef USCKF_DEBUG_PRINTS
                    std::cout << "[INDIRECT_UPDATE] TRUE Update"<<std::endl;
                    std::cout << "[INDIRECT_UPDATE] H*xk_i:\n" << H*xk_i <<std::endl;
                    #endif

                    xk_i = xk_i + K * innovation;
                    Pk_i = (Eigen::Matrix<ScalarType, DOF_SINGLE_STATE, DOF_SINGLE_STATE>::Identity()
                            -K * H) * Pk_i *(Eigen::Matrix<ScalarType, DOF_SINGLE_STATE, DOF_SINGLE_STATE>::Identity()
                            -K * H).transpose() + K * R * K.transpose();
                    Pk_i = 0.5 * (Pk_i + Pk_i.transpose()); //! Guarantee symmetry
                }

                /** Store the sub-covariance matrix for statek_i **/
                MTK::subblock (Pk, &_AugmentedState::statek_i) = Pk_i;

                #ifdef USCKF_DEBUG_PRINTS
                std::cout << "[INDIRECT_UPDATE] xk_i(after):\n" << xk_i <<std::endl;
                std::cout << "[INDIRECT_UPDATE] Pk(after):\n" << Pk_i <<std::endl;
                std::cout << "[INDIRECT_UPDATE] K:\n" << K <<std::endl;
                std::cout << "[INDIRECT_UPDATE] S:\n" << S <<std::endl;
                std::cout << "[INDIRECT_UPDATE] z:\n" << z <<std::endl;
                std::cout << "[INDIRECT_UPDATE] innovation:\n" << innovation <<std::endl;
                std::cout << "[INDIRECT_UPDATE] R is of size:" <<R.rows()<<"x"<<R.cols()<<std::endl;
                std::cout << "[INDIRECT_UPDATE] R:\n" << R <<std::endl;
                #endif

                /**************************/
                /** Apply the Corrections */
                /**************************/

                Eigen::Quaterniond qe;

                /** Update the quaternion with the Indirect approach **/
                qe.w() = 1;
                qe.x() = xk_i(3);
                qe.y() = xk_i(4);
                qe.z() = xk_i(5);

                /** Apply correction **/
                mu_state.statek_i.pos += xk_i.template block<3, 1>(0,0);
                mu_state.statek_i.orient = (mu_state.statek_i.orient * qe);
                mu_state.statek_i.orient.normalize();
                mu_state.statek_i.velo += xk_i.template block<3, 1>(6, 0);
                mu_state.statek_i.angvelo += xk_i.template block<3, 1>(9, 0);

            }

            void cloning(int mode)
            {

                SingleStateCovariance Pk_i, Pk_l;

                switch (mode)
                {
                case STATEK_I:
                    /** Augmented state cloning **/
                    mu_state.statek_l = mu_state.statek_i;

                    /** Covariance state cloning **/
                    Pk_i = MTK::subblock (Pk, &_AugmentedState::statek_i);
                    MTK::subblock (Pk, &_AugmentedState::statek_l, &_AugmentedState::statek_l) = Pk_i;
                    MTK::subblock (Pk, &_AugmentedState::statek_l, &_AugmentedState::statek_i) = Pk_i;
                    MTK::subblock (Pk, &_AugmentedState::statek_i, &_AugmentedState::statek_l) = Pk_i;
                    break;

                case STATEK_L:
                    /** Augmented state cloning **/
                    mu_state.statek = mu_state.statek_l;

                    /** Covariance state cloning **/
                    Pk_l = MTK::subblock (Pk, &_AugmentedState::statek_l);
                    MTK::subblock (Pk, &_AugmentedState::statek, &_AugmentedState::statek) = Pk_l;
                    MTK::subblock (Pk, &_AugmentedState::statek, &_AugmentedState::statek_l) = Pk_l;
                    MTK::subblock (Pk, &_AugmentedState::statek_l, &_AugmentedState::statek) = Pk_l;
                    break;

                case STATEK_I + STATEK_L:
                    /** Augmented state cloning **/
                    mu_state.statek_l = mu_state.statek_i;
                    mu_state.statek = mu_state.statek_l;

                    /** Covariance state cloning **/
                    Pk_i = MTK::subblock (Pk, &_AugmentedState::statek_i);
                    MTK::subblock (Pk, &_AugmentedState::statek_l, &_AugmentedState::statek_l) = Pk_i;
                    MTK::subblock (Pk, &_AugmentedState::statek_l, &_AugmentedState::statek_i) = Pk_i;
                    MTK::subblock (Pk, &_AugmentedState::statek_i, &_AugmentedState::statek_l) = Pk_i;
                    MTK::subblock (Pk, &_AugmentedState::statek, &_AugmentedState::statek) = Pk_i;
                    MTK::subblock (Pk, &_AugmentedState::statek, &_AugmentedState::statek_l) = Pk_i;
                    MTK::subblock (Pk, &_AugmentedState::statek_l, &_AugmentedState::statek) = Pk_i;
                    break;

                default:
                    break;
                }
            }

            void setSingleState(const _SingleState & state, int order = STATEK_I)
            {
                switch (order)
                {
                case STATEK_I:
                    mu_state.statek_i = state;
                    break;
                default:
                    break;
                }

            }

            _SingleState muSingleState(int state = STATEK_I)
            {
                _SingleState statek_i;

                switch(state)
                {
                case STATEK_I:
                    statek_i = mu_state.statek_i;
                    break;
                case STATEK_L:
                    statek_i =  mu_state.statek_l;
                    break;
                case STATEK:
                    statek_i = mu_state.statek;
                    break;
                default:
                    statek_i = mu_state.statek_i;
                    break;
                }

                return statek_i;
            }

            void setPkSingleState(const SingleStateCovariance & Pk_i, int order = STATEK_I)
            {
                switch (order)
                {
                case STATEK_I:
                    MTK::subblock (Pk, &_AugmentedState::statek_i) = Pk_i;
                    break;
                default:
                    break;
                }

            }

            SingleStateCovariance PkSingleState(int state = STATEK_I)
            {
                SingleStateCovariance Pk_i;

                switch(state)
                {
                case STATEK_I:
                    Pk_i = MTK::subblock (Pk, &_AugmentedState::statek_i);
                    break;
                case STATEK_L:
                    Pk_i = MTK::subblock (Pk, &_AugmentedState::statek_l);
                    break;
                case STATEK:
                    Pk_i = MTK::subblock (Pk, &_AugmentedState::statek);
                    break;
                default:
                    Pk_i = MTK::subblock (Pk, &_AugmentedState::statek_i);
                    break;
                }

                return Pk_i;
            }

            const _AugmentedState& muState() const
            {
                return mu_state;
            }

            const AugmentedStateCovariance &PkAugmentedState() const
            {
                return Pk;
            }

    private:

            /**@brief Sigma Point Calculation for the complete Augmented State
            */
            void generateSigmaPoints(const _AugmentedState &mu, const AugmentedStateCovariance &sigma, AugmentedStateSigma &X) const
            {
                    generateSigmaPoints(mu, VectorizedAugmentedState::Zero(), sigma, X);
            }

            /**@brief Sigma Point Calculation for the complete Augmented State
             */
            void generateSigmaPoints(const _AugmentedState &mu, const VectorizedAugmentedState &delta,
                                    const AugmentedStateCovariance &sigma, AugmentedStateSigma &X) const
            {
                    assert(X.size() == 2 * DOF_AUGMENTED_STATE + 1);

                    ukfom::lapack::cholesky<DOF_AUGMENTED_STATE> L(sigma);

                    if (!L.isSPD())
                    {
                            std::cerr << std::endl << "sigma is not SPD:" << std::endl
                                              << sigma << std::endl
                                              << "---" << std::endl;
                            Eigen::EigenSolver<AugmentedStateCovariance> eig(sigma);
                            std::cerr << "eigen values: " << eig.eigenvalues().transpose() << std::endl;
                    }

                    assert(L.isSPD());


                    /*std::cout << ">> L" << std::endl
                                      << L.getL() << std::endl
                                      << "<< L" << std::endl;
                    */

                    X[0] = mu + delta;
                    for (std::size_t i = 1, j = 0; j < DOF_AUGMENTED_STATE; ++j)
                    {
                            //std::cout << "L.col(" << j << "): " << L.getL().col(j).transpose() << std::endl;
                            X[i++] = mu + (delta + L.getL().col(j));
                            X[i++] = mu + (delta - L.getL().col(j));
                    }
                    #ifdef USCKF_DEBUG_PRINTS
                    printSigmaPoints<AugmentedStateSigma>(X);
                    #endif
            }

            /**@brief Sigma Point Calculation for the Single State
            */
            void generateSigmaPoints(const _SingleState &mu, const SingleStateCovariance &sigma, SingleStateSigma &X) const
            {
                    generateSigmaPoints(mu, VectorizedSingleState::Zero(), sigma, X);
            }

            /**@brief Sigma Point Calculation for the Single State
             */
            void generateSigmaPoints(const _SingleState &mu, const VectorizedSingleState &delta,
                                    const SingleStateCovariance &sigma, SingleStateSigma &X) const
            {
                    assert(X.size() == 2 * DOF_SINGLE_STATE + 1);

                    Eigen::LLT< SingleStateCovariance > lltOfSigma(sigma); // compute the Cholesky decomposition of A
                    SingleStateCovariance L = lltOfSigma.matrixL(); // retrieve factor L  in the decomposition


                    /*std::cout << ">> L" << std::endl
                                      << L << std::endl
                                      << "<< L" << std::endl;
                     std::cout<<"L*L^T:\n"<< L * L.transpose()<<"\n";*/


                    X[0] = mu + delta;
                    for (std::size_t i = 1, j = 0; j < DOF_SINGLE_STATE; ++j)
                    {
                            //std::cout << "L.col(" << j << "): " << L.getL().col(j).transpose() << std::endl;
                            X[i++] = mu + (delta + L.col(j));
                            X[i++] = mu + (delta - L.col(j));
                    }

                    #ifdef USCKF_DEBUG_PRINTS
                    printSigmaPoints<SingleStateSigma>(X);
                    #endif
            }

            // manifold mean
            template<typename _Manifold>
            _Manifold meanSigmaPoints(const std::vector<_Manifold> &X) const
            {
                    _Manifold reference = X[0];
                    typename _Manifold::vectorized_type mean_delta;
                    const static std::size_t max_it = 10000;

                    std::size_t i = 0;
                    do {
                            mean_delta.setZero();
                            for (typename std::vector<_Manifold>::const_iterator Xi = X.begin(); Xi != X.end(); ++Xi)
                            {
                                    mean_delta += *Xi - reference;
                            }
                            mean_delta /= X.size();
                            reference += mean_delta;
                    } while (mean_delta.norm() > 1e-6
                                     && ++i < max_it);

                    if (i >= max_it)
                    {
                            std::cerr << "ERROR: meanSigmaPoints() did not converge. norm(mean_delta)=" << mean_delta.norm() << std::endl;
                            assert(false);
                    }

                    return reference;
            }

            // vector mean
            template<int _MeasurementRows>
            Eigen::Matrix<ScalarType, _MeasurementRows, 1>
            meanSigmaPoints(const std::vector<Eigen::Matrix<ScalarType, _MeasurementRows, 1> > &Z) const
            {
                    typedef Eigen::Matrix<ScalarType, _MeasurementRows, 1> Measurement;

                    return std::accumulate(Z.begin(), Z.end(), Measurement(Measurement::Zero())) / Z.size();
            }

#ifdef VECT_H_
            // MTK vector mean
            template<int _MeasurementRows>
            MTK::vect<_MeasurementRows, ScalarType>
            meanSigmaPoints(const std::vector<MTK::vect<_MeasurementRows, ScalarType> > &Z) const
            {
                    typedef MTK::vect<_MeasurementRows, ScalarType> Measurement;

                    return std::accumulate(Z.begin(), Z.end(), Measurement(Measurement::Zero())) / Z.size();
            }
#endif // VECT_H_

            template<int _CovSize, typename T>
            Eigen::Matrix<ScalarType, _CovSize, _CovSize>
            covSigmaPoints(const T &mean, const std::vector<T> &V) const
            {
                    typedef Eigen::Matrix<ScalarType, _CovSize, _CovSize> CovMat;
                    typedef Eigen::Matrix<ScalarType, _CovSize, 1> CovCol;

                    CovMat c(CovMat::Zero());

                    for (typename std::vector<T>::const_iterator Vi = V.begin(); Vi != V.end(); ++Vi)
                    {
                            CovCol d = *Vi - mean;
                            c += d * d.transpose();
                    }

                    return 0.5 * c;
            }

            template<typename _State, int _MeasurementRows, typename _SigmaPoints, typename _Measurement>
            Eigen::Matrix<ScalarType, _State::DOF, _MeasurementRows>
            crossCovSigmaPoints(const _State &meanX, const _Measurement &meanZ,
                                const _SigmaPoints &X, const std::vector<_Measurement> &Z) const
            {
                    assert(X.size() == Z.size());

                    typedef Eigen::Matrix<ScalarType, _State::DOF, _MeasurementRows> CrossCov;

                    CrossCov c(CrossCov::Zero());

                    {
                            typename _SigmaPoints::const_iterator Xi = X.begin();
                            typename std::vector<_Measurement>::const_iterator Zi = Z.begin();
                            for (;Zi != Z.end(); ++Xi, ++Zi)
                            {
                                    c += (*Xi - meanX) * (*Zi - meanZ).transpose();
                            }
                    }

                    return 0.5 * c;
            }

            void applyDeltaAugmentedState(const VectorizedAugmentedState &delta)
            {
                    AugmentedStateSigma X(2 * DOF_AUGMENTED_STATE + 1);
                    generateSigmaPoints(mu_state, delta, Pk, X);

                    mu_state = meanSigmaPoints(X);
                    Pk = covSigmaPoints<_AugmentedState::DOF>(mu_state, X);
            }

            void applyDeltaSingleState(_SingleState &statek_i, SingleStateCovariance &Pk_i, const VectorizedSingleState &delta)
            {
                    SingleStateSigma X(2 * DOF_SINGLE_STATE + 1);
                    generateSigmaPoints(statek_i, delta, Pk_i, X);

                    statek_i = meanSigmaPoints(X);
                    Pk_i = covSigmaPoints<_SingleState::DOF>(statek_i, X);
            }

            // for debugging only
            template <typename _SigmaType>
            void printSigmaPoints(const _SigmaType &X) const
            {
                    std::cout << "generated sigma points:" << std::endl;
                    for (typename _SigmaType::const_iterator Xi = X.begin(); Xi != X.end(); ++Xi)
                    {
                            std::cout << *Xi << std::endl << "***" << std::endl;
                    }
            }

    public:
            void checkSigmaPoints()
            {
                AugmentedStateSigma X(2 * DOF_AUGMENTED_STATE + 1);
                generateSigmaPoints(mu_state, Pk, X);

                _AugmentedState muX = meanSigmaPoints(X);

                AugmentedStateCovariance Pktest = covSigmaPoints<_AugmentedState::DOF>(muX, X);
                if((Pktest - Pk).cwise().abs().maxCoeff()>1e-6){
                        std::cerr << Pktest << "\n\n" << Pk;
                        assert(false);
                }

                if (mu_state != muX)
                {
                    //std::cout << "mu_:" << mu_ << std::endl;
                    //std::cout << "muX:" << muX << std::endl;
                    std::cout << "norm:" << ((mu_state - muX).norm() > 0. ? ">" : "=") << std::endl;
                }
                assert (mu_state == muX);
            }

    public:
            EIGEN_MAKE_ALIGNED_OPERATOR_NEW

            template <typename _ScalarType>
            static bool accept_mahalanobis_distance(const _ScalarType &mahalanobis2, const int dof)
            {
                #ifdef USCKF_DEBUG_PRINTS
                std::cout << "[MAHALANOBIS_DISTANCE] mahalanobis2: " << mahalanobis2 <<std::endl;
                std::cout << "[MAHALANOBIS_DISTANCE] dof: " << dof <<std::endl;
                #endif


                /** Only significance of alpha = 5% is computed **/
                switch (dof)
                {
                    case 1:
                        if (mahalanobis2 < 3.84)
                            return true;
                        else
                            return false;
                    case 2:
                        if (mahalanobis2 < 5.99)
                            return true;
                        else
                            return false;
                    case 3:
                        if (mahalanobis2 < 7.81)
                            return true;
                        else
                            return false;
                    case 4:
                        if (mahalanobis2 < 9.49)
                            return true;
                        else
                            return false;
                    case 5:
                        if (mahalanobis2 < 11.07)
                            return true;
                        else
                            return false;
                    case 6:
                        if (mahalanobis2 < 12.59)
                            return true;
                        else
                            return false;
                    case 7:
                        if (mahalanobis2 < 14.07)
                            return true;
                        else
                            return false;
                    case 8:
                        if (mahalanobis2 < 15.51)
                            return true;
                        else
                            return false;
                    case 9:
                        if (mahalanobis2 < 16.92)
                            return true;
                        else
                            return false;
                    default:
                        std::cerr<<"[MAHALANOBIS-ERROR] DoF("<<dof<<") not supported"<<std::endl;
                        return false;
                }
            };
    };
	
} // namespace localization

#endif // __USCKF_HPP_
