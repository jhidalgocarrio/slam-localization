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

/** MTK's pose and orientation definition **/
#include <mtk/startIdx.hpp>

/** Localization includes **/
#include <localization/tools/Util.hpp>

//#define USCKF_DEBUG_PRINTS 1

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

            /** Types related to Single State **/
            typedef typename _SingleState::vectorized_type VectorizedSingleState;
            typedef Eigen::Matrix<ScalarType, int(_SingleState::DOF), int(_SingleState::DOF)> SingleStateCovariance;
            typedef std::vector<_SingleState> SingleStateSigma;

            /** Types related to Augmented State **/
            typedef typename _AugmentedState::vectorized_type VectorizedAugmentedState;
            typedef Eigen::Matrix<ScalarType, Eigen::Dynamic, Eigen::Dynamic> AugmentedStateCovariance;
            typedef std::vector<_AugmentedState> AugmentedStateSigma;

            /** Covariance only for the multi-state part **/
            typedef Eigen::Matrix<ScalarType, int(_AugmentedState::DOF), int(_AugmentedState::DOF)> MultiStateCovariance;

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
                MultiStateCovariance Pk_states(Pk.block(0, 0, Pk_states.rows(), Pk_states.cols()));
                SingleStateCovariance Pk_i = MTK::subblock (Pk_states, &_AugmentedState::statek_i);

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
                SingleStateCovariance Fk = Pxy.transpose()*Pk_i.inverse(); // Fk = Pyx * (Pk_i)^-1

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
                MTK::subblock (Pk_states, &_AugmentedState::statek_i) = Pk_i;

                /** Compute the Cross Cov for the Copy States of the AugmentedState **/
                SingleStateCovariance Pkkk;

                /*******************************/
                /**  Cross-Covariance Matrix  **/
                /*******************************/

                /** Covariance between state(k) and state(k+i), Pk|k+i**/
                Pkkk = MTK::subblock (Pk_states, &_AugmentedState::statek, &_AugmentedState::statek_i);
                Pkkk = Pkkk * Fk.transpose();
                MTK::subblock (Pk_states, &_AugmentedState::statek, &_AugmentedState::statek_i) = Pkkk;

                /** Covariance between state(k+l) and state(k+i), Pk+l|k+i **/
                Pkkk = MTK::subblock (Pk_states, &_AugmentedState::statek_l, &_AugmentedState::statek_i);
                Pkkk = Pkkk * Fk.transpose();
                MTK::subblock (Pk_states, &_AugmentedState::statek_l, &_AugmentedState::statek_i) = Pkkk;

                /** Covariance between state(k+i) and state(k), Pk+i|k **/
                Pkkk = MTK::subblock (Pk_states, &_AugmentedState::statek_i, &_AugmentedState::statek);
                Pkkk = Fk * Pkkk;
                MTK::subblock (Pk_states, &_AugmentedState::statek_i, &_AugmentedState::statek) = Pkkk;

                /** Covariance between state(k+i) and state(k+l), Pk+i|Pk+l **/
                Pkkk = MTK::subblock (Pk_states, &_AugmentedState::statek_i, &_AugmentedState::statek_l);
                Pkkk = Fk * Pkkk;
                MTK::subblock (Pk_states, &_AugmentedState::statek_i, &_AugmentedState::statek_l) = Pkkk;

                /** Store the Pk_states into the complete Pk **/
                Pk.block(0, 0, Pk_states.rows(), Pk_states.cols()) = Pk_states;

                /*******************/
                /** Cross-Covariance with the features (Dynamic size part of the state vector) **/
                /*******************/

                Eigen::Matrix<ScalarType, _SingleState::DOF, Eigen::Dynamic> Pzk, Pzkl;
                Pzk.resize(_SingleState::DOF, mu_state.featuresk.size());
                Pzkl.resize(_SingleState::DOF, mu_state.featuresk_l.size());

                /** Covariance between state(k+i) and features(k), Px_(k+i),z_k **/
                Pzk = Pk.block(_AugmentedState::DOF-_SingleState::DOF, _AugmentedState::DOF, Pzk.rows(), Pzk.cols());
                Pzk = Fk * Pzk;
                Pk.block(_AugmentedState::DOF-_SingleState::DOF, _AugmentedState::DOF, Pzk.rows(), Pzk.cols()) = Pzk;

                /** Covariance between features(k) and state(k+i), Pz_(k),x_(k+i) **/
                Pk.block(_AugmentedState::DOF, _AugmentedState::DOF-_SingleState::DOF, Pzk.cols(), Pzk.rows()) = Pzk.transpose();

                /** Covariance between state(k+i) and features(k+l), Px_(k+i),z_(k+l) **/
                Pzkl = Pk.block(_AugmentedState::DOF-_SingleState::DOF, _AugmentedState::DOF+Pzk.cols(), Pzkl.rows(), Pzkl.cols());
                Pzkl = Fk * Pzkl;
                Pk.block(_AugmentedState::DOF-_SingleState::DOF, _AugmentedState::DOF+Pzk.cols(), Pzkl.rows(), Pzkl.cols()) = Pzkl;

                /** Covariance between features(k+l) and state(k+i), Pz_(k+l),x_(k+i)**/
                Pk.block(_AugmentedState::DOF+Pzk.transpose().rows(), _AugmentedState::DOF-_SingleState::DOF, Pzkl.cols(), Pzkl.rows()) = Pzkl.transpose();

                /** No cross-Covariance across measurements **/

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
                    typedef Eigen::Matrix<ScalarType, Eigen::Dynamic, measurement_rows> CrossCov;

                    AugmentedStateSigma X(2 * mu_state.getDOF() + 1);
                    VectorizedAugmentedState mu_delta(mu_state.getDOF(), 1);
                    mu_delta.setZero();
                    std::cout<<"[USCKF_UPDATE] X.size(): "<< X.size() <<"\n";
                    std::cout<<"[USCKF_UPDATE] mu_delta.size(): "<< mu_delta.size() <<"\n";
                    generateSigmaPoints(mu_state, mu_delta, Pk, X);

                    std::vector<Measurement> Z(X.size());
                    std::transform(X.begin(), X.end(), Z.begin(), h);

                    const Measurement meanZ = meanSigmaPoints(Z);
                    std::cout<<"meanZ: "<<meanZ<<"\n";
                    std::cout<<"measurement_rows: "<<measurement_rows<<"\n";

                    const MeasurementCov S = covSigmaPoints(meanZ, Z) + R();
                    const CrossCov covXZ = crossCovSigmaPoints(mu_state, meanZ, X, Z);

                    MeasurementCov S_inverse;
                    S_inverse = S.inverse();

                    const CrossCov K = covXZ * S_inverse;

                    const VectorizedMeasurement innovation = z - meanZ;

                    const ScalarType mahalanobis2 = (innovation.transpose() * S_inverse * innovation)(0);

                    if (mt(mahalanobis2))
                    {
                            Pk -= K * S * K.transpose();
                            //applyDelta(K * innovation);
                            std::cout<<"K "<<K.rows() <<" x "<<K.cols()<<"\n";
                            _AugmentedState innovation_state;
                            innovation_state.set(K * innovation, mu_state.featuresk.size(), mu_state.featuresk_l.size());
                            mu_state = mu_state + innovation_state;
                    }

                    //#ifdef USCKF_DEBUG_PRINTS
                    std::cout << "[USCKF_UPDATE] innovation:" << std::endl << innovation << std::endl;
                    std::cout << "[USCKF_UPDATE] mu_state':" << std::endl << mu_state << std::endl;
                    //#endif
            }

            template<typename _Measurement>
            void pushMeasurementIntoState(_Measurement &z_k_i)
            {

               /** State k+l is now state k and the associated covariances **/
               this->cloning(STATEK_L);

               /** State k+i is now state k+l and the associated covariances **/
               this->cloning(STATEK_I);

               /** Push a new set of features measurements **/
               mu_state.featuresk = mu_state.featuresk_l;
               mu_state.featuresk_l = z_k_i;

               /** Resize the process state covariance matrix to the new dimension **/
               Pk.resize(_AugmentedState::DOF + mu_state.featuresk.size() + mu_state.featuresk_l.size(),
                         _AugmentedState::DOF + mu_state.featuresk.size() + mu_state.featuresk_l.size());

               /** Get block for the state|measurement z(k+l) cross covariance matrices **/
               Eigen::Matrix<ScalarType, _AugmentedState::DOF, Eigen::Dynamic> Pzkl_block;
               Pzkl_block.resize(_AugmentedState::DOF, mu_state.featuresk_l.size());

               /** Move block for the state|measurement z(k+l) to the z(k) **/
               Pk.block(0, _AugmentedState::DOF, _AugmentedState::DOF, mu_state.featuresk_l.size()) = Pzkl_block;

               /** Move block for the measurement|state z(k+l) to the z(k) **/
               Pk.block(_AugmentedState::DOF, 0, mu_state.featuresk_l.size(), _AugmentedState::DOF) = Pzkl_block.transpose();

            }

            void cloning(int mode)
            {
                SingleStateCovariance Pk_i, Pk_l, Pk_l_i;

                switch (mode)
                {
                case STATEK_I:
                    /** Augmented state cloning, state k+l = state k+i**/
                    mu_state.statek_l = mu_state.statek_i;

                    /** Covariance state cloning, Pk+l = Pk+i, Pk+l|k+i = Pk+i, Pk+i|k+l = Pk+i **/
                    Pk_i = MTK::subblock (Pk, &_AugmentedState::statek_i);
                    Pk_l_i = MTK::subblock(Pk, &_AugmentedState::statek_l, &_AugmentedState::statek_i);
                    MTK::subblock (Pk, &_AugmentedState::statek_l, &_AugmentedState::statek_l) = Pk_i;
                    MTK::subblock (Pk, &_AugmentedState::statek_l, &_AugmentedState::statek_i) = Pk_i;
                    MTK::subblock (Pk, &_AugmentedState::statek_i, &_AugmentedState::statek_l) = Pk_i;

                    /** Covariance state cloning, Pk|k+i = Pk+l|k+i, Pk+i|k = Pk+i|k+l **/
                    MTK::subblock (Pk, &_AugmentedState::statek, &_AugmentedState::statek_i) = Pk_l_i;
                    MTK::subblock (Pk, &_AugmentedState::statek_i, &_AugmentedState::statek) = Pk_l_i.transpose();

                    break;

                case STATEK_L:
                    /** Augmented state cloning, state k = state k+l **/
                    mu_state.statek = mu_state.statek_l;

                    /** Covariance state cloning, Pk = Pk+l, Pk|k+l = Pk+l, Pk+l|k = Pk+l **/
                    Pk_l = MTK::subblock (Pk, &_AugmentedState::statek_l);
                    MTK::subblock (Pk, &_AugmentedState::statek, &_AugmentedState::statek) = Pk_l;
                    MTK::subblock (Pk, &_AugmentedState::statek, &_AugmentedState::statek_l) = Pk_l;
                    MTK::subblock (Pk, &_AugmentedState::statek_l, &_AugmentedState::statek) = Pk_l;
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

                case STATEK_L:
                    mu_state.statek_l = state;
                    break;

                case STATEK:
                    mu_state.statek_l = state;
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
                MultiStateCovariance Pk_states;
                Pk_states = Pk.block(0, 0, Pk_states.rows(), Pk_states.cols());

                switch(state)
                {
                case STATEK_I:
                    Pk_i = MTK::subblock (Pk_states, &_AugmentedState::statek_i);
                    break;
                case STATEK_L:
                    Pk_i = MTK::subblock (Pk_states, &_AugmentedState::statek_l);
                    break;
                case STATEK:
                    Pk_i = MTK::subblock (Pk_states, &_AugmentedState::statek);
                    break;
                default:
                    Pk_i = MTK::subblock (Pk_states, &_AugmentedState::statek_i);
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
            void generateSigmaPoints(const _AugmentedState &mu, const VectorizedAugmentedState &delta,
                                    const AugmentedStateCovariance &sigma, AugmentedStateSigma &X) const
            {
                    assert(X.size() == static_cast<unsigned int> (2 * mu_state.getDOF() + 1));

                    Eigen::LLT< AugmentedStateCovariance > lltOfSigma(sigma); // compute the Cholesky decomposition of A
                    AugmentedStateCovariance L = lltOfSigma.matrixL(); // retrieve factor L  in the decomposition

                    /*std::cout << "L is of size "<<L.rows()<<" x "<<L.cols()<<"\n";
                    std::cout << ">> L" << std::endl
                                      << L << std::endl
                                      << "<< L" << std::endl;
                     std::cout<<"L*L^T:\n"<< L * L.transpose()<<"\n";*/

                    _AugmentedState delta_state;
                    delta_state.set(delta, mu.featuresk.size(), mu.featuresk_l.size());

                    X[0] = mu + delta_state;
                    for (register unsigned int i = 1, j = 0; j < mu.getDOF(); ++j)
                    {
                            //std::cout << "L.col(" << j << "): " << L.col(j).transpose() << std::endl;
                            _AugmentedState l_state;
                            l_state.set(L.col(j), mu.featuresk.size(), mu.featuresk_l.size());
                            X[i++] = mu + (delta_state + l_state);
                            X[i++] = mu + (delta_state - l_state);
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

                    Measurement measurement_zero(Z[0].size(), 1);
                    measurement_zero.setZero();

                    return std::accumulate(Z.begin(), Z.end(), measurement_zero) / Z.size();
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

            Eigen::Matrix<ScalarType, Eigen::Dynamic, Eigen::Dynamic>
            covSigmaPoints(const Eigen::Matrix<ScalarType, Eigen::Dynamic, 1>  &mean,
                    const std::vector< Eigen::Matrix<ScalarType, Eigen::Dynamic, 1> > &V) const
            {
                    typedef Eigen::Matrix<ScalarType, Eigen::Dynamic, Eigen::Dynamic> CovMat;
                    typedef Eigen::Matrix<ScalarType, Eigen::Dynamic, 1> CovCol;

                    CovMat c(mean.size(), mean.size());
                    c.setZero();

                    for (typename std::vector< Eigen::Matrix<ScalarType, Eigen::Dynamic, 1> >::const_iterator Vi = V.begin(); Vi != V.end(); ++Vi)
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

            template<typename _State, typename _SigmaPoints>
            Eigen::Matrix<ScalarType, Eigen::Dynamic, Eigen::Dynamic>
            crossCovSigmaPoints(const _State &meanX, const Eigen::Matrix<ScalarType, Eigen::Dynamic, 1> &meanZ,
                                const _SigmaPoints &X, const std::vector< Eigen::Matrix<ScalarType, Eigen::Dynamic, 1> > &Z) const
            {
                    assert(X.size() == Z.size());

                    typedef Eigen::Matrix<ScalarType, Eigen::Dynamic, Eigen::Dynamic> CrossCov;

                    CrossCov c(mu_state.getDOF(), meanZ.size());
                    c.setZero();

                    {
                            typename _SigmaPoints::const_iterator Xi = X.begin();
                            typename std::vector< Eigen::Matrix<ScalarType, Eigen::Dynamic, 1> >::const_iterator Zi = Z.begin();
                            for (;Zi != Z.end(); ++Xi, ++Zi)
                            {
                                    _State tempXi (*Xi - meanX);
                                    c += tempXi.getVectorizedState() * (*Zi - meanZ).transpose();
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
