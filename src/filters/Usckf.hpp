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
	
    template <typename _VectorState, typename _SingleState>
    class Usckf
    {
        typedef Usckf self;

        public:
            enum
            {
                    DOF_STATE = _VectorState::DOF
            };
            enum
            {
                    DOF_SINGLE_STATE = _SingleState::DOF
            };


            typedef typename _VectorState::scalar_type ScalarType;
            typedef typename _VectorState::vectorized_type VectorizedVectorState;
            typedef typename _SingleState::vectorized_type VectorizedSingleState;
            typedef Eigen::Matrix<ScalarType, int(_VectorState::DOF), int(_VectorState::DOF)> VectorStateCovariance;
            typedef Eigen::Matrix<ScalarType, int(_SingleState::DOF), int(_SingleState::DOF)> SingleStateCovariance;
            typedef std::vector<_VectorState> VectorStateSigma;
            typedef std::vector<_SingleState> SingleStateSigma;


        private:

            _VectorState mu_state; /** Mean of the state vector **/
            _VectorState mu_error; /** Mean of the error State vector **/
            VectorStateCovariance Pk_error; /** Covariance of the error State vector **/

        public:
            /**@brief Constructor
             */
            Usckf(const _VectorState &state, const _VectorState &error, const VectorStateCovariance &P0)
                : mu_state(state), mu_error(error),Pk_error(P0)
            {
            }

            /**@brief Propagate the state of the navigation quantities
             */
            void propagate(const _SingleState & delta_state)
            {
                mu_state.statek_i.pos += delta_state.pos; //! Increment in position
                mu_state.statek_i.vel += delta_state.vel; //! Increment in velocity
                mu_state.statek_i.orient = mu_state.statek_i.orient * delta_state.orient; //! Delta in orientation
            }

            /**@brief Set current State
             */
            void setStatek_i(const _SingleState & state)
            {
                mu_state.statek_i = state;
            }

            /**@brief Filter prediction step
             */
            template<typename ProcessModel>
            void predict(ProcessModel g, const SingleStateCovariance &Q)
            {
                predict(g, boost::bind(ukfom::id<SingleStateCovariance>, Q));
            }

            template<typename ProcessModel, typename ProcessNoiseCovariance>
            void predict(ProcessModel g, ProcessNoiseCovariance Q)
            {
                /** Get the current state vector error to propagate **/
                _SingleState statek_i = mu_error.statek_i;
                SingleStateCovariance Pk = MTK::subblock (Pk_error, &_VectorState::statek_i);

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

                /** Apply the non-linear tranformation of the process model **/
                std::transform(X.begin(), X.end(), X.begin(), g);

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
                Qk = Q();//Q()*dt + 0.5*dt*dt*Fk*Q() + 0.5*dt*dt*Q()*Fk.transpose();
                //Qk = 0.5*(Qk + Qk.transpose());

                /** Just for debugging purpose print, this should be equal Pk **/
                /** Because to Pyy =  Fk*Pk*Fk^T and  Pk = Pyy+Q **/
                #ifdef USCKF_DEBUG_PRINTS
                std::cout<<"[USCKF_PREDICT] Fk*Pk*Fk^T + Qk:\n"<< Fk*Pk*Fk.transpose() + Qk<< std::endl;
                #endif

                /************************/
                /** Covariance Matrix  **/
                /************************/

                /** Compute the Process model Covariance **/
                Pk = covSigmaPoints<_SingleState::DOF, _SingleState>(mu_error.statek_i, X) + Qk;

                /** Store the subcovariance matrix for statek_i **/
                MTK::subblock (Pk_error, &_VectorState::statek_i) = Pk;

                /** Compute the Cross Cov for the Copy States of the VectorState **/
                SingleStateCovariance Pkkk;

                /************************/
                /**  Cross-Cov Matrix  **/
                /************************/

                /** Covariance between state(k) and state(k+i) **/
                Pkkk = MTK::subblock (Pk_error, &_VectorState::statek, &_VectorState::statek_i);
                Pkkk = Pkkk * Fk.transpose();
                MTK::subblock (Pk_error, &_VectorState::statek, &_VectorState::statek_i) = Pkkk;

                /** Covariance between state(k+l) and state(k+i) **/
                Pkkk = MTK::subblock (Pk_error, &_VectorState::statek_l, &_VectorState::statek_i);
                Pkkk = Pkkk * Fk.transpose();
                MTK::subblock (Pk_error, &_VectorState::statek_l, &_VectorState::statek_i) = Pkkk;

                /** Covariance between state(k+i) and state(k) **/
                Pkkk = MTK::subblock (Pk_error, &_VectorState::statek_i, &_VectorState::statek);
                Pkkk = Fk * Pkkk;
                MTK::subblock (Pk_error, &_VectorState::statek_i, &_VectorState::statek) = Pkkk;

                /** Covariance between state(k+i) and state(k+l) **/
                Pkkk = MTK::subblock (Pk_error, &_VectorState::statek_i, &_VectorState::statek_l);
                Pkkk = Fk * Pkkk;
                MTK::subblock (Pk_error, &_VectorState::statek_i, &_VectorState::statek_l) = Pkkk;

                /**********/
                /** TO-DO: cross-cov with the features (Dynamic size part of the filter) **/

                #ifdef  USCKF_DEBUG_PRINTS
                std::cout << "[USCKF_PREDICT] statek_i(k+1|k):" << std::endl << mu_error.statek_i << std::endl;
                std::cout << "[USCKF_PREDICT] Pk(k+1|k):"<< std::endl << Pk << std::endl;
                std::cout << "[USCKF_PREDICT] Process Noise Cov Q(k):"<< std::endl << Q() << std::endl;
                #endif
            }

            template<typename Measurement, typename MeasurementModel, typename MeasurementNoiseCovariance>
            void update(const Measurement &z,MeasurementModel h, MeasurementNoiseCovariance R)
            {
                    update(z, h, R, ukfom::accept_any_mahalanobis_distance<ScalarType>);
            }

            template<typename Measurement, typename MeasurementModel>
            void update(const Measurement &z, MeasurementModel h,
                        const Eigen::Matrix<ScalarType, ukfom::dof<Measurement>::value, ukfom::dof<Measurement>::value> &R)
            {
                    typedef Eigen::Matrix<ScalarType, ukfom::dof<Measurement>::value, ukfom::dof<Measurement>::value> measurement_cov;
                    update(z, h, boost::bind(ukfom::id<measurement_cov>, R), ukfom::accept_any_mahalanobis_distance<ScalarType>);
            }

            template<typename Measurement, typename MeasurementModel,
                    typename MeasurementNoiseCovariance, typename SignificanceTest>
            void update(const Measurement &z, MeasurementModel h,
                        MeasurementNoiseCovariance R, SignificanceTest mt)
            {
                    const static int measurement_rows = ukfom::dof<Measurement>::value;
                    typedef Measurement Measurement;
                    typedef Eigen::Matrix<ScalarType, measurement_rows, 1> VectorizedMeasurement;
                    typedef std::vector<Measurement> measurement_vector;
                    typedef Eigen::Matrix<ScalarType, measurement_rows, measurement_rows> MeasurementCov;
                    typedef Eigen::Matrix<ScalarType, _VectorState::DOF, measurement_rows> CrossCov;

                    SingleStateSigma X(2 * DOF_STATE + 1);
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
                    std::cout << "[USCKF_PREDICT] innovation:" << std::endl << innovation << std::endl;
                    std::cout << "[USCKF_PREDICT] mu_error':" << std::endl << mu_error << std::endl;
                    #endif
            }

            /** @brief Standard EKF Update of the state (TO-DO: improve this part using Manifold
             * and Unscented tranform (UKF))
             */
            template<typename Measurement, typename MeasurementModel, typename MeasurementNoiseCovariance>
            void ekfUpdate(VectorizedSingleState &xk_i, const Measurement &z, MeasurementModel H, MeasurementNoiseCovariance R)
            {
                const static int DOF_MEASUREMENT = ukfom::dof<Measurement>::value; /** Dimension of the measurement */

                /** statek_i cov matrix **/
                SingleStateCovariance Pk = MTK::subblock (Pk_error, &_VectorState::statek_i);

                #ifdef USCKF_DEBUG_PRINTS
                std::cout << "[USCKF_EKF_UPDATE] xk_i:\n" << xk_i <<std::endl;
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
                MTK::subblock (Pk_error, &_VectorState::statek_i) = Pk;

                #ifdef USCKF_DEBUG_PRINTS
                std::cout << "[USCKF_EKF_UPDATE] xk_i:\n" << xk_i <<std::endl;
                std::cout << "[USCKF_EKF_UPDATE] Pk:\n" << Pk <<std::endl;
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

                /** Store the subcovariance matrix for statek_i **/
                MTK::subblock (Pk_error, &_VectorState::statek_i) = Pk;
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
                /** Vector state clonnning (error and magnitudes state) **/
                mu_state.statek_l = mu_state.statek_i;
                mu_state.statek = mu_state.statek_l;

                mu_error.statek_l = mu_error.statek_i;
                mu_error.statek = mu_error.statek_l;

                /** Covariance  state clonnning **/
                SingleStateCovariance Pk = MTK::subblock (Pk_error, &_VectorState::statek_i);
                MTK::subblock (Pk_error, &_VectorState::statek, &_VectorState::statek) = Pk;
                MTK::subblock (Pk_error, &_VectorState::statek_l, &_VectorState::statek_l) = Pk;
                MTK::subblock (Pk_error, &_VectorState::statek_i, &_VectorState::statek_i) = Pk;

                MTK::subblock (Pk_error, &_VectorState::statek, &_VectorState::statek_l) = Pk;
                MTK::subblock (Pk_error, &_VectorState::statek, &_VectorState::statek_i) = Pk;
                MTK::subblock (Pk_error, &_VectorState::statek_l, &_VectorState::statek) = Pk;
                MTK::subblock (Pk_error, &_VectorState::statek_i, &_VectorState::statek) = Pk;
                MTK::subblock (Pk_error, &_VectorState::statek_l, &_VectorState::statek_i) = Pk;
                MTK::subblock (Pk_error, &_VectorState::statek_i, &_VectorState::statek_l) = Pk;

            }

            const _VectorState& muState() const
            {
                return mu_state;
            }

            const _VectorState& muError() const
            {
                return mu_error;
            }

            const VectorStateCovariance &PkVectorState() const
            {
                return Pk_error;
            }


    private:

            /**@brief Sigma Point Calculation for the complete Vector State
             */
            void generateSigmaPoints(const _VectorState &mu, const VectorizedVectorState &delta,
                                    const VectorStateCovariance &sigma, VectorStateSigma &X) const
            {
                    assert(X.size() == 2 * DOF_STATE + 1);

                    ukfom::lapack::cholesky<DOF_STATE> L(sigma);

                    if (!L.isSPD())
                    {
                            std::cerr << std::endl << "sigma is not SPD:" << std::endl
                                              << sigma << std::endl
                                              << "---" << std::endl;
                            Eigen::EigenSolver<VectorStateCovariance> eig(sigma);
                            std::cerr << "eigen values: " << eig.eigenvalues().transpose() << std::endl;
                    }

                    assert(L.isSPD());


                    /*std::cout << ">> L" << std::endl
                                      << L.getL() << std::endl
                                      << "<< L" << std::endl;
                    */

                    X[0] = mu + delta;
                    for (std::size_t i = 1, j = 0; j < DOF_STATE; ++j)
                    {
                            //std::cout << "L.col(" << j << "): " << L.getL().col(j).transpose() << std::endl;
                            X[i++] = mu + (delta + L.getL().col(j));
                            X[i++] = mu + (delta - L.getL().col(j));
                    }
                    printSigmaPoints<VectorStateSigma>(X);
            }
            /**@brief Sigma Point Calculation for the complete Vector State
             */
            void generateSigmaPoints(const _VectorState &mu, const VectorStateCovariance &sigma, VectorStateSigma &X) const
            {
                    generateSigmaPoints(mu, VectorizedVectorState::Zero(), sigma, X);
            }

            /**@brief Sigma Point Calculation for the State (SingleState)
             */
            void generateSigmaPoints(const _SingleState &mu, 
                    const SingleStateCovariance &sigma, SingleStateSigma &X) const
            {
                    assert(X.size() == 2 * DOF_SINGLE_STATE + 1);

                    ukfom::lapack::cholesky<DOF_SINGLE_STATE> L(sigma);

                    if (!L.isSPD())
                    {
                            std::cerr << std::endl << "sigma is not SPD:" << std::endl
                                              << sigma << std::endl
                                              << "---" << std::endl;
                            Eigen::EigenSolver<SingleStateCovariance> eig(sigma);
                            std::cerr << "eigen values: " << eig.eigenvalues().transpose() << std::endl;
                    }

                    assert(L.isSPD());


                    /*std::cout << ">> L" << std::endl
                                      << L.getL() << std::endl
                                      << "<< L" << std::endl;
                    */

                    X[0] = mu;
                    for (std::size_t i = 1, j = 0; j < DOF_SINGLE_STATE; ++j)
                    {
                            //std::cout << "L.col(" << j << "): " << L.getL().col(j).transpose() << std::endl;
                            X[i++] = mu + L.getL().col(j);
                            X[i++] = mu + (-L.getL().col(j));
                    }
                    printSigmaPoints<SingleStateSigma>(X);
            }

            // manifold mean
            template<typename manifold>
            manifold meanSigmaPoints(const std::vector<manifold> &X) const
            {
                    manifold reference = X[0];
                    typename manifold::vectorized_type mean_delta;
                    const static std::size_t max_it = 10000;

                    std::size_t i = 0;
                    do {
                            mean_delta.setZero();
                            for (typename std::vector<manifold>::const_iterator Xi = X.begin(); Xi != X.end(); ++Xi)
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
            template<int measurement_rows>
            Eigen::Matrix<ScalarType, measurement_rows, 1>
            meanSigmaPoints(const std::vector<Eigen::Matrix<ScalarType, measurement_rows, 1> > &Z) const
            {
                    typedef Eigen::Matrix<ScalarType, measurement_rows, 1> measurement;

                    return std::accumulate(Z.begin(), Z.end(), measurement(measurement::Zero())) / Z.size();
            }

#ifdef VECT_H_
            // MTK vector mean
            template<int measurement_rows>
            MTK::vect<measurement_rows, ScalarType>
            meanSigmaPoints(const std::vector<MTK::vect<measurement_rows, ScalarType> > &Z) const
            {
                    typedef MTK::vect<measurement_rows, ScalarType> measurement;

                    return std::accumulate(Z.begin(), Z.end(), measurement(measurement::Zero())) / Z.size();
            }
#endif // VECT_H_

            template<int cov_size, typename T>
            Eigen::Matrix<ScalarType, cov_size, cov_size>
            covSigmaPoints(const T &mean, const std::vector<T> &V) const
            {
                    typedef Eigen::Matrix<ScalarType, cov_size, cov_size> cov_mat;
                    typedef Eigen::Matrix<ScalarType, cov_size, 1> cov_col;

                    cov_mat c(cov_mat::Zero());

                    for (typename std::vector<T>::const_iterator Vi = V.begin(); Vi != V.end(); ++Vi)
                    {
                            cov_col d = *Vi - mean;
                            c += d * d.transpose();
                    }

                    return 0.5 * c;
            }

            template<typename _State, int measurement_rows, typename _SigmaPoints, typename Measurement>
            Eigen::Matrix<ScalarType, _State::DOF, measurement_rows>
            crossCovSigmaPoints(const _State &meanX, const Measurement &meanZ,
                                const _SigmaPoints &X, const std::vector<Measurement> &Z) const
            {
                    assert(X.size() == Z.size());

                    typedef Eigen::Matrix<ScalarType, _State::DOF, measurement_rows> cross_cov;

                    cross_cov c(cross_cov::Zero());

                    {
                            typename _SigmaPoints::const_iterator Xi = X.begin();
                            typename std::vector<Measurement>::const_iterator Zi = Z.begin();
                            for (;Zi != Z.end(); ++Xi, ++Zi)
                            {
                                    c += (*Xi - meanX) * (*Zi - meanZ).transpose();
                            }
                    }

                    return 0.5 * c;
            }

            void applyDelta(const VectorizedVectorState &delta)
            {
                    SingleStateSigma X(2 * DOF_STATE + 1);
                    generateSigmaPoints(mu_error, delta, Pk_error, X);

                    mu_error = meanSigmaPoints(X);
                    Pk_error = covSigmaPoints<_VectorState::DOF>(mu_error, X);
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
                VectorStateSigma X(2 * DOF_STATE + 1);
                generateSigmaPoints(mu_error, Pk_error, X);

                _VectorState muX = meanSigmaPoints(X);

                VectorStateCovariance Pk_errortest = covSigmaPoints<_VectorState::DOF>(muX, X);
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

    public:
            EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    };
	
} // namespace rover_localization

#endif // __USCKF_HPP_
