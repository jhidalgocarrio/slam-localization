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
#include <mtk/types/pose.hpp>
#include <mtk/types/SOn.hpp>
#include <mtk/build_manifold.hpp>


#define USCKF_DEBUG_PRINTS 1

namespace localization
{
	
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
            _AugmentedState mu_error; /** Mean of the error State vector **/
            AugmentedStateCovariance Pk_error; /** Covariance of the error State vector **/

        public:
            /**@brief Constructor
             */
            Usckf(const _AugmentedState &state, const _AugmentedState &error, const AugmentedStateCovariance &P0)
                : mu_state(state), mu_error(error),Pk_error(P0)
            {
            }

            /**@brief Set current State mu_state
            */
            void setStatek_i(const _SingleState & state)
            {
                mu_state.statek_i = state;
            }

            /**@brief Prediction step in the linearized form of an EKF
             */
            template<typename _ProcessModel>
            void ekfPredict (const _ProcessModel &F, const SingleStateCovariance &Q)
            {
                SingleStateCovariance Pk = MTK::subblock (Pk_error, &_AugmentedState::statek_i);

                /** Propagate the vector through the system **/
                mu_error.statek_i.set(F * mu_error.statek_i.getVectorizedState (_SingleState::ERROR_QUATERNION));

                /** Propagate the P covariance matrix **/
                Pk = F*Pk*F.transpose() + Q;

                /** Store the subcovariance matrix for statek_i **/
                MTK::subblock (Pk_error, &_AugmentedState::statek_i) = Pk;

                /************************/
                /**  Cross-Cov Matrix  **/
                /************************/

                /** Compute the Cross Cov for the Copy States of the AugmentedState **/
                SingleStateCovariance Pkkk;

                /** Covariance between state(k) and state(k+i) **/
                Pkkk = MTK::subblock (Pk_error, &_AugmentedState::statek, &_AugmentedState::statek_i);
                Pkkk = Pkkk * F.transpose();
                MTK::subblock (Pk_error, &_AugmentedState::statek, &_AugmentedState::statek_i) = Pkkk;

                /** Covariance between state(k+l) and state(k+i) **/
                Pkkk = MTK::subblock (Pk_error, &_AugmentedState::statek_l, &_AugmentedState::statek_i);
                Pkkk = Pkkk * F.transpose();
                MTK::subblock (Pk_error, &_AugmentedState::statek_l, &_AugmentedState::statek_i) = Pkkk;

                /** Covariance between state(k+i) and state(k) **/
                Pkkk = MTK::subblock (Pk_error, &_AugmentedState::statek_i, &_AugmentedState::statek);
                Pkkk = F * Pkkk;
                MTK::subblock (Pk_error, &_AugmentedState::statek_i, &_AugmentedState::statek) = Pkkk;

                /** Covariance between state(k+i) and state(k+l) **/
                Pkkk = MTK::subblock (Pk_error, &_AugmentedState::statek_i, &_AugmentedState::statek_l);
                Pkkk = F * Pkkk;
                MTK::subblock (Pk_error, &_AugmentedState::statek_i, &_AugmentedState::statek_l) = Pkkk;

                /**********/
                /** TO-DO: cross-cov with the features (Dynamic size part of the filter) **/

                #ifdef  USCKF_DEBUG_PRINTS
                std::cout << "[EKF_PREDICT] statek_i(k+1|k):" << std::endl << mu_error.statek_i << std::endl;
                std::cout << "[EKF_PREDICT] Pk(k+1|k):"<< std::endl << Pk << std::endl;
                std::cout << "[EKF_PREDICT] Process Noise Cov Q(k):"<< std::endl << Q << std::endl;
                #endif

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
                /** Get the current state vector error to propagate **/
                _SingleState statek_i = mu_error.statek_i;
                SingleStateCovariance Pk = MTK::subblock (Pk_error, &_AugmentedState::statek_i);

                #ifdef USCKF_DEBUG_PRINTS
                std::cout<<"[USCKF_PREDICT] statek_i(k|k):\n"<<statek_i<<"\n";
                std::cout<<"[USCKF_PREDICT] P(k|k):\n"<<Pk<<"\n";
                #endif

                /** Propagation only uses the current state (DOF_SINGLE_STATE dimension ) **/
                SingleStateSigma X(2 * DOF_SINGLE_STATE + 1);

                /** Generates the sigma Points of a Single State **/
                generateSigmaPoints(statek_i, Pk, X);

                /** Create a copy before the transformation **/
                SingleStateSigma XCopy(2 * DOF_SINGLE_STATE + 1);
                XCopy = X;

                /*****************************/
                /** Process Model Transform **/
                /*****************************/

                /** Apply the non-linear transformation of the process model **/
                std::transform(X.begin(), X.end(), X.begin(), f);

                printSigmaPoints<SingleStateSigma>(X);

                /** Compute the mean **/
                mu_error.statek_i = meanSigmaPoints(X);

                /** Compute the cross-covariance matrix (cross because is between X and X before the transform) **/
                SingleStateCovariance Pxy;
                Pxy = crossCovSigmaPoints<_SingleState, _SingleState::DOF, SingleStateSigma,
                    _SingleState> (statek_i, mu_error.statek_i, XCopy, X);
                SingleStateCovariance Fk = Pxy.transpose()*Pk.inverse();

                #ifdef USCKF_DEBUG_PRINTS
                std::cout<<"[USCKF_PREDICT] Fk:\n"<< Fk << std::endl;
                #endif

                /*************************/
                /** Discretization of Q **/
                /*************************/

                SingleStateCovariance Qk;
                /** TO-DO change it **/
                Qk = Q();
                //Qk = Q()*dt + 0.5*dt*dt*Fk*Q() + 0.5*dt*dt*Q()*Fk.transpose();
                //Qk = 0.5*(Qk + Qk.transpose());

                /** Just for debugging purpose print **/
                /** Pyy =  Fk*Pk*Fk^T and  Pk = Pyy+Q **/
                /** After propagation/prediction, which means Pyy+Q = P(k+1|k) **/
                #ifdef USCKF_DEBUG_PRINTS
                std::cout<<"[USCKF_PREDICT] Fk*Pk*Fk^T + Qk:\n"<< Fk*Pk*Fk.transpose() + Qk<< std::endl;
                #endif

                /************************/
                /** Covariance Matrix  **/
                /************************/

                /** Compute the Process model Covariance **/
                Pk = covSigmaPoints<_SingleState::DOF, _SingleState>(mu_error.statek_i, X) + Qk;

                /** Store the subcovariance matrix for statek_i **/
                MTK::subblock (Pk_error, &_AugmentedState::statek_i) = Pk;

                /** Compute the Cross Cov for the Copy States of the AugmentedState **/
                SingleStateCovariance Pkkk;

                /************************/
                /**  Cross-Cov Matrix  **/
                /************************/

                /** Covariance between state(k) and state(k+i) **/
                Pkkk = MTK::subblock (Pk_error, &_AugmentedState::statek, &_AugmentedState::statek_i);
                Pkkk = Pkkk * Fk.transpose();
                MTK::subblock (Pk_error, &_AugmentedState::statek, &_AugmentedState::statek_i) = Pkkk;

                /** Covariance between state(k+l) and state(k+i) **/
                Pkkk = MTK::subblock (Pk_error, &_AugmentedState::statek_l, &_AugmentedState::statek_i);
                Pkkk = Pkkk * Fk.transpose();
                MTK::subblock (Pk_error, &_AugmentedState::statek_l, &_AugmentedState::statek_i) = Pkkk;

                /** Covariance between state(k+i) and state(k) **/
                Pkkk = MTK::subblock (Pk_error, &_AugmentedState::statek_i, &_AugmentedState::statek);
                Pkkk = Fk * Pkkk;
                MTK::subblock (Pk_error, &_AugmentedState::statek_i, &_AugmentedState::statek) = Pkkk;

                /** Covariance between state(k+i) and state(k+l) **/
                Pkkk = MTK::subblock (Pk_error, &_AugmentedState::statek_i, &_AugmentedState::statek_l);
                Pkkk = Fk * Pkkk;
                MTK::subblock (Pk_error, &_AugmentedState::statek_i, &_AugmentedState::statek_l) = Pkkk;

                /**********/
                /** TO-DO: cross-cov with the features (Dynamic size part of the filter) **/

                #ifdef  USCKF_DEBUG_PRINTS
                std::cout << "[USCKF_PREDICT] statek_i(k+1|k):" << std::endl << mu_error.statek_i << std::endl;
                std::cout << "[USCKF_PREDICT] Pk(k+1|k):"<< std::endl << Pk << std::endl;
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
                    generateSigmaPoints(mu_error, Pk_error, X);

                    std::vector<Measurement> Z(X.size());
                    std::transform(X.begin(), X.end(), Z.begin(), h);

                    const Measurement meanZ = meanSigmaPoints(Z);
                    const MeasurementCov S = covSigmaPoints<measurement_rows>(meanZ, Z) + R();
                    const CrossCov covXZ = crossCovSigmaPoints<measurement_rows>(mu_error, meanZ, X, Z);

                    MeasurementCov S_inverse;
                    S_inverse = S.inverse();

                    const CrossCov K = covXZ * S_inverse;

                    const VectorizedMeasurement innovation = z - meanZ;

                    const ScalarType mahalanobis2 = (innovation.transpose() * S_inverse * innovation)(0);

                    if (mt(mahalanobis2))
                    {
                            Pk_error -= K * S * K.transpose();
                            //applyDelta(K * innovation);
                            mu_error = mu_error + K * innovation;
                    }

                    #ifdef USCKF_DEBUG_PRINTS
                    std::cout << "[USCKF_UPDATE] innovation:" << std::endl << innovation << std::endl;
                    std::cout << "[USCKF_UPDATE] mu_error':" << std::endl << mu_error << std::endl;
                    #endif
            }

            template<typename _Measurement, typename _MeasurementModel>
            void singleUpdate(const _Measurement &z, _MeasurementModel h,
                        const Eigen::Matrix<ScalarType, ukfom::dof<_Measurement>::value, ukfom::dof<_Measurement>::value> &R)
            {
                    typedef Eigen::Matrix<ScalarType, ukfom::dof<_Measurement>::value, ukfom::dof<_Measurement>::value> measurement_cov;
                    singleUpdate(z, h, boost::bind(ukfom::id<measurement_cov>, R));
            }

            /** @brief Single UKF Update of the state
             */
            template<typename _Measurement, typename _MeasurementModel, typename _MeasurementNoiseCovariance>
            void singleUpdate(const _Measurement &z, _MeasurementModel h, _MeasurementNoiseCovariance R)
            {
                const static int DOF_MEASUREMENT = ukfom::dof<_Measurement>::value; /** Dimension of the measurement */
                typedef Eigen::Matrix<ScalarType, DOF_MEASUREMENT, DOF_MEASUREMENT> MeasurementCov;
                typedef Eigen::Matrix<ScalarType, _SingleState::DOF, DOF_MEASUREMENT> CrossCov;

                /** Get the current state vector error to propagate **/
                _SingleState statek_i = mu_error.statek_i;

                /** statek_i cov matrix **/
                SingleStateCovariance Pk = MTK::subblock (Pk_error, &_AugmentedState::statek_i);

                #ifdef USCKF_DEBUG_PRINTS
                std::cout << "[USCKF_SINGLE_UPDATE] statek_i(k+1|k):\n" << statek_i <<std::endl;
                std::cout << "[USCKF_SINGLE_UPDATE] Pk(k+1|k):\n" << Pk <<std::endl;
                #endif

                SingleStateSigma X(2 * DOF_SINGLE_STATE + 1);
                generateSigmaPoints(statek_i, Pk, X);

                std::vector<_Measurement> Z(X.size());
                std::transform(X.begin(), X.end(), Z.begin(), h);

                /** Mean of the measurement model **/
                const _Measurement meanZ = meanSigmaPoints(Z);

                /** The innovation covariance matrix **/
                const MeasurementCov S = covSigmaPoints<DOF_MEASUREMENT, _Measurement>(meanZ, Z) + R();

                /** The cross-correlation matrix **/
                const CrossCov covXZ = crossCovSigmaPoints<_SingleState, DOF_MEASUREMENT, SingleStateSigma, _Measurement>(statek_i, meanZ, X, Z);

                MeasurementCov S_inverse;
                S_inverse = S.inverse();

                const CrossCov K = covXZ * S_inverse;

                const _Measurement innovation = z - meanZ;

                Pk -= K * S * K.transpose();
                Pk = 0.5 * (Pk + Pk.transpose()); //! Guarantee symmetry
                statek_i = statek_i + K * innovation;

                /** Store the subcovariance matrix for statek_i **/
                MTK::subblock (Pk_error, &_AugmentedState::statek_i) = Pk;

                #ifdef USCKF_DEBUG_PRINTS
                std::cout << "[USCKF_SINGLE_UPDATE] statek_i(k+1|k+1):\n" << statek_i <<std::endl;
                std::cout << "[USCKF_SINGLE_UPDATE] Pk(k+1|k+1):\n" << Pk <<std::endl;
                std::cout << "[USCKF_SINGLE_UPDATE] K:\n" << K <<std::endl;
                std::cout << "[USCKF_SINGLE_UPDATE] S:\n" << S <<std::endl;
                std::cout << "[USCKF_SINGLE_UPDATE] z:\n" << z <<std::endl;
                std::cout << "[USCKF_SINGLE_UPDATE] meanZ:\n" << meanZ <<std::endl;
                std::cout << "[USCKF_SINGLE_UPDATE] innovation:\n" << innovation <<std::endl;
                std::cout << "[USCKF_SINGLE_UPDATE] R is of size:" <<R().rows()<<"x"<<R().cols()<<std::endl;
                std::cout << "[USCKF_SINGLE_UPDATE] R:\n" << R() <<std::endl;
                #endif

                /**************************/
                /** Apply the Corrections */
                /**************************/

                /** Apply correction **/
                mu_state.statek_i.pos += statek_i.pos;
                mu_state.statek_i.vel += statek_i.vel;
                mu_state.statek_i.orient = (mu_state.statek_i.orient * statek_i.orient);
                mu_state.statek_i.orient.normalize();
                mu_state.statek_i.gbias += statek_i.gbias;
                mu_state.statek_i.abias += statek_i.abias;

            }


            /** @brief Standard EKF Update of the state (TO-DO: improve this part using Manifold
             * and Unscented tranform (UKF))
             */
            template<typename _Measurement, typename _MeasurementModel, typename _MeasurementNoiseCovariance>
            void ekfUpdate(VectorizedSingleState &xk_i, const _Measurement &z, _MeasurementModel H, _MeasurementNoiseCovariance R)
            {
                const static int DOF_MEASUREMENT = ukfom::dof<_Measurement>::value; /** Dimension of the measurement */

                /** statek_i cov matrix **/
                SingleStateCovariance Pk = MTK::subblock (Pk_error, &_AugmentedState::statek_i);

                #ifdef USCKF_DEBUG_PRINTS
                std::cout << "[USCKF_EKF_UPDATE] xk_i(before):\n" << xk_i <<std::endl;
                std::cout << "[USCKF_EKF_UPDATE] Pk(before):\n" << Pk <<std::endl;
                #endif

                /** Compute the Kalman Gain Matrix **/
                Eigen::Matrix<ScalarType, DOF_MEASUREMENT, DOF_MEASUREMENT> S;
                Eigen::Matrix<ScalarType, DOF_SINGLE_STATE, DOF_MEASUREMENT> K;
                S = H * Pk * H.transpose() + R; //!Calculate the cov of the innovation
                K = Pk * H.transpose() * S.inverse(); //!Calculate K using the inverse of S

                /** Update the state vector and the covariance matrix */
                xk_i = xk_i + K * (z - H * xk_i);
                Pk = (Eigen::Matrix<ScalarType, DOF_SINGLE_STATE, DOF_SINGLE_STATE>::Identity()
                        -K * H) * Pk *(Eigen::Matrix<ScalarType, DOF_SINGLE_STATE, DOF_SINGLE_STATE>::Identity()
                        -K * H).transpose() + K * R * K.transpose();
                Pk = 0.5 * (Pk + Pk.transpose()); //! Guarantee symmetry

                /** Store the subcovariance matrix for statek_i **/
                MTK::subblock (Pk_error, &_AugmentedState::statek_i) = Pk;

                #ifdef USCKF_DEBUG_PRINTS
                std::cout << "[USCKF_EKF_UPDATE] xk_i(after):\n" << xk_i <<std::endl;
                std::cout << "[USCKF_EKF_UPDATE] Pk(after):\n" << Pk <<std::endl;
                std::cout << "[USCKF_EKF_UPDATE] K:\n" << K <<std::endl;
                std::cout << "[USCKF_EKF_UPDATE] S:\n" << S <<std::endl;
                std::cout << "[USCKF_EKF_UPDATE] z:\n" << z <<std::endl;
                std::cout << "[USCKF_EKF_UPDATE] H*xk_i:\n" << H*xk_i <<std::endl;
                std::cout << "[USCKF_EKF_UPDATE] innovation:\n" << (z - H * xk_i) <<std::endl;
                std::cout << "[USCKF_EKF_UPDATE] R is of size:" <<R.rows()<<"x"<<R.cols()<<std::endl;
                std::cout << "[USCKF_EKF_UPDATE] R:\n" << R <<std::endl;
                #endif

                /**************************/
                /** Apply the Corrections */
                /**************************/

                Eigen::Quaterniond qe;

                /** Update the quaternion with the Indirect approach **/
                qe.w() = 1;
                qe.x() = xk_i(6);
                qe.y() = xk_i(7);
                qe.z() = xk_i(8);

                /** Apply correction **/
                mu_state.statek_i.pos += xk_i.template block<3, 1>(0,0);
                mu_state.statek_i.vel += xk_i.template block<3, 1>(3,0);
                mu_state.statek_i.orient = (mu_state.statek_i.orient * qe);
                mu_state.statek_i.orient.normalize();
                mu_state.statek_i.gbias += xk_i.template block<3, 1>(9, 0);
                mu_state.statek_i.abias += xk_i.template block<3, 1>(12, 0);

            }

            void muErrorSingleReset()
            {
                mu_error.statek_i.pos.setZero();
                mu_error.statek_i.vel.setZero();
                mu_error.statek_i.orient.setIdentity();
                mu_error.statek_i.gbias.setZero();
                mu_error.statek_i.abias.setZero();
            }

            void clonning()
            {
                /** Augmented state clonnning (error and magnitudes state) **/
                mu_state.statek_l = mu_state.statek_i;
                mu_state.statek = mu_state.statek_l;

                mu_error.statek_l = mu_error.statek_i;
                mu_error.statek = mu_error.statek_l;

                /** Covariance  state clonnning **/
                SingleStateCovariance Pk = MTK::subblock (Pk_error, &_AugmentedState::statek_i);
                MTK::subblock (Pk_error, &_AugmentedState::statek, &_AugmentedState::statek) = Pk;
                MTK::subblock (Pk_error, &_AugmentedState::statek_l, &_AugmentedState::statek_l) = Pk;
                MTK::subblock (Pk_error, &_AugmentedState::statek_i, &_AugmentedState::statek_i) = Pk;

                MTK::subblock (Pk_error, &_AugmentedState::statek, &_AugmentedState::statek_l) = Pk;
                MTK::subblock (Pk_error, &_AugmentedState::statek, &_AugmentedState::statek_i) = Pk;
                MTK::subblock (Pk_error, &_AugmentedState::statek_l, &_AugmentedState::statek) = Pk;
                MTK::subblock (Pk_error, &_AugmentedState::statek_i, &_AugmentedState::statek) = Pk;
                MTK::subblock (Pk_error, &_AugmentedState::statek_l, &_AugmentedState::statek_i) = Pk;
                MTK::subblock (Pk_error, &_AugmentedState::statek_i, &_AugmentedState::statek_l) = Pk;

            }

            const _AugmentedState& muState() const
            {
                return mu_state;
            }

            const _AugmentedState& muError() const
            {
                return mu_error;
            }

            const AugmentedStateCovariance &PkAugmentedState() const
            {
                return Pk_error;
            }


    private:

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
                    printSigmaPoints<AugmentedStateSigma>(X);
            }
            /**@brief Sigma Point Calculation for the complete Augmented State
             */
            void generateSigmaPoints(const _AugmentedState &mu, const AugmentedStateCovariance &sigma, AugmentedStateSigma &X) const
            {
                    generateSigmaPoints(mu, VectorizedAugmentedState::Zero(), sigma, X);
            }

            /**@brief Sigma Point Calculation for the State (SingleState)
             */
            void generateSigmaPoints(const _SingleState &mu, 
                    const SingleStateCovariance &sigma, SingleStateSigma &X) const
            {
                    assert(X.size() == 2 * DOF_SINGLE_STATE + 1);

                    Eigen::LLT< SingleStateCovariance > lltOfSigma(sigma); // compute the Cholesky decomposition of A
                    SingleStateCovariance L = lltOfSigma.matrixL(); // retrieve factor L  in the decomposition


                    /*std::cout << ">> L" << std::endl
                                      << L << std::endl
                                      << "<< L" << std::endl;
                     std::cout<<"L*L^T:\n"<< L * L.transpose()<<"\n";*/


                    X[0] = mu;
                    for (std::size_t i = 1, j = 0; j < DOF_SINGLE_STATE; ++j)
                    {
                            //std::cout << "L.col(" << j << "): " << L.getL().col(j).transpose() << std::endl;
                            X[i++] = mu + L.col(j);
                            X[i++] = mu + (-L.col(j));
                    }
                    printSigmaPoints<SingleStateSigma>(X);
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

            void applyDelta(const VectorizedAugmentedState &delta)
            {
                    SingleStateSigma X(2 * DOF_AUGMENTED_STATE + 1);
                    generateSigmaPoints(mu_error, delta, Pk_error, X);

                    mu_error = meanSigmaPoints(X);
                    Pk_error = covSigmaPoints<_AugmentedState::DOF>(mu_error, X);
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
                generateSigmaPoints(mu_error, Pk_error, X);

                _AugmentedState muX = meanSigmaPoints(X);

                AugmentedStateCovariance Pk_errortest = covSigmaPoints<_AugmentedState::DOF>(muX, X);
                if((Pk_errortest - Pk_error).cwise().abs().maxCoeff()>1e-6){
                        std::cerr << Pk_errortest << "\n\n" << Pk_error;
                        assert(false);
                }

                if (mu_error != muX)
                {
                    //std::cout << "mu_:" << mu_ << std::endl;
                    //std::cout << "muX:" << muX << std::endl;
                    std::cout << "norm:" << ((mu_error - muX).norm() > 0. ? ">" : "=") << std::endl;
                }
                assert (mu_error == muX);
            }

            template <typename _MatrixType>
            _MatrixType guaranteeSPD (const _MatrixType &A)
            {
                _MatrixType spdA;
                Eigen::VectorXd s;
                s.resize(A.rows(), 1);

                /**
                 * Single Value Decomposition
                */
                Eigen::JacobiSVD <Eigen::MatrixXd > svdOfA (A, Eigen::ComputeThinU | Eigen::ComputeThinV);

                s = svdOfA.singularValues(); //!eigenvalues
                std::cout<<"[SVD] s: \n"<<s<<"\n";
                std::cout<<"[SVD] svdOfA.matrixU():\n"<<svdOfA.matrixU()<<"\n";
                std::cout<<"[SVD] svdOfA.matrixV():\n"<<svdOfA.matrixV()<<"\n";

                Eigen::EigenSolver<_MatrixType> eig(A);
                std::cout << "[SVD] BEFORE: eigen values: " << eig.eigenvalues().transpose() << std::endl;

                for (register int i=0; i<s.size(); ++i)
                {
                    std::cout<<"[SVD] i["<<i<<"]\n";

                    if (s(i) < 0.00)
                        s(i) = 0.00;
                }
                spdA = svdOfA.matrixU() * s.matrix().asDiagonal() * svdOfA.matrixV();

                Eigen::EigenSolver<_MatrixType> eigSPD(spdA);
                if (eig.eigenvalues() == eigSPD.eigenvalues())
                    std::cout<<"[SVD] EQUAL!!\n";

                std::cout << "[SVD] AFTER: eigen values: " << eigSPD.eigenvalues().transpose() << std::endl;

                Eigen::EigenSolver<_MatrixType> eigSPD2(spdA);
                std::cout << "[SVD] AGAIN: eigen values: " << eigSPD2.eigenvalues().transpose() << std::endl;

                return spdA;
            }

    public:
            EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    };
	
} // namespace rover_localization

#endif // __USCKF_HPP_
