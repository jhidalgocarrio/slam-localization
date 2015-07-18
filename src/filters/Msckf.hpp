#ifndef _MSCKF_HPP_
#define _MSCKF_HPP_

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

#define MSCKF_DEBUG_PRINTS 1

namespace localization
{

    template <typename _MultiState, typename _SingleState>
    class Msckf
    {
        typedef Msckf self;

        public:
            enum
            {
                    DOF_MULTI_STATE = _MultiState::DOF
            };

            enum
            {
                    DOF_SINGLE_STATE = _MultiState::SingleState::DOF
            };

            enum
            {
                    SENSOR_DOF = _MultiState::SENSOR_DOF
            };



            typedef typename _MultiState::scalar_type ScalarType;

            /** Types related to State **/
            typedef typename _MultiState::vectorized_type VectorizedMultiState;
            typedef Eigen::Matrix<ScalarType, int(_MultiState::DOF), int(_MultiState::DOF)> MultiStateCovariance;
            typedef std::vector<_MultiState> MultiStateSigma;


            typedef typename _SingleState::vectorized_type VectorizedSingleState;
            typedef Eigen::Matrix<ScalarType, int(_SingleState::DOF), int(_SingleState::DOF)> SingleStateCovariance;
            typedef std::vector<_SingleState> SingleStateSigma;

        private:

            _MultiState mu_state; /** Mean of the state and sensors pose vector **/
            MultiStateCovariance Pk; /** Covariance of the State and sensor vector **/

        public:
            /**@brief Constructor
             */
            Msckf(const _MultiState &state, const MultiStateCovariance &P0)
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
                _SingleState statek_i = this->mu_state.statek;
                SingleStateCovariance Pk_i = MTK::subblock (this->Pk, &_MultiState::statek);

                #ifdef MSCKF_DEBUG_PRINTS
                std::cout<<"[MSCKF_PREDICT] statek_i(k|k):\n"<<statek_i<<"\n";
                std::cout<<"[MSCKF_PREDICT] Pk_i(k|k):\n"<<Pk_i<<"\n";
                #endif

                /** Propagation only uses the current state (DOF_SINGLE_STATE dimension ) **/
                SingleStateSigma X(2 * DOF_SINGLE_STATE + 1);

                /** Generates the sigma Points of a Single State **/
                this->generateSigmaPoints(statek_i, Pk_i, X);

                /** Create a copy before the transformation **/
                SingleStateSigma XCopy(2 * DOF_SINGLE_STATE + 1);
                XCopy = X;

                /*****************************/
                /** Process Model Transform **/
                /*****************************/

                /** Apply the non-linear transformation of the process model **/
                std::transform(X.begin(), X.end(), X.begin(), f);

                #ifdef MSCKF_DEBUG_PRINTS
                //this->printSigmaPoints<SingleStateSigma>(X);
                #endif

                /** Compute the mean **/
                this->mu_state.statek = this->meanSigmaPoints(X);

                /** Compute the cross-covariance matrix (cross because is between X and X before the transform) **/
                SingleStateCovariance Pxy;
                Pxy = this->crossCovSigmaPoints<_SingleState, DOF_SINGLE_STATE, SingleStateSigma, _SingleState>
                    (statek_i, mu_state.statek, XCopy, X);
                SingleStateCovariance Fk = Pxy.transpose()*Pk_i.inverse(); // Fk = Pyx * (Pk_i)^-1

                #ifdef MSCKF_DEBUG_PRINTS
                std::cout<<"[MSCKF_PREDICT] Fk:\n"<< Fk << std::endl;
                #endif

                /*************************/
                /** Discretization of Q **/
                /*************************/

                SingleStateCovariance Qk = Q();

                /** Just for debugging purpose print **/
                /** Pyy =  Fk*Pk_i*Fk^T and  Pk_i = Pyy+Q **/
                /** After propagation/prediction, which means Pyy+Q = P(k+1|k) **/
                #ifdef MSCKF_DEBUG_PRINTS
                std::cout<<"[MSCKF_PREDICT] Fk*Pk_i*Fk^T + Qk:\n"<< Fk*Pk_i*Fk.transpose() + Qk<< std::endl;
                #endif

                /************************/
                /** Covariance Matrix  **/
                /************************/

                /** Compute the Process model Covariance **/
                Pk_i = this->covSigmaPoints<DOF_SINGLE_STATE, _SingleState>(mu_state.statek, X) + Qk;

                /** Store the subcovariance matrix for statek **/
                MTK::subblock (this->Pk, &_MultiState::statek) = Pk_i;

                /*******************************/
                /**  Cross-Covariance Matrix  **/
                /*******************************/

                Eigen::Matrix<ScalarType, _SingleState::DOF, _MultiState::SENSOR_DOF> Pkk;

                /** Covariance between state and sensor poses **/
                Pkk = Pk.template block<_SingleState::DOF, _MultiState::SENSOR_DOF>(0, _SingleState::DOF);
                Pkk = Fk * Pkk;
                Pk.template block<_SingleState::DOF, _MultiState::SENSOR_DOF>(0, _SingleState::DOF) = Pkk;

                /** Covariance between sensor poses and state **/
                Pkk.transpose() = Pk.template block<_MultiState::SENSOR_DOF,_SingleState::DOF>(_SingleState::DOF, 0);
                Pkk.transpose() = Pkk.transpose() * Fk.transpose();
                Pk.template block<_MultiState::SENSOR_DOF, _SingleState::DOF>(_SingleState::DOF, 0) = Pkk.transpose();

                #ifdef  MSCKF_DEBUG_PRINTS
                std::cout << "[MSCKF_PREDICT] statek_i(k+1|k):" << std::endl << mu_state.statek << std::endl;
                std::cout << "[MSCKF_PREDICT] Pk(k+1|k):"<< std::endl << Pk_i << std::endl;
                std::cout << "[MSCKF_PREDICT] Process Noise Cov Q(k):"<< std::endl << Q() << std::endl;
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
                    typedef Eigen::Matrix<ScalarType, measurement_rows, measurement_rows> MeasurementCov;
                    typedef Eigen::Matrix<ScalarType, Eigen::Dynamic, measurement_rows> CrossCov;

                    MultiStateSigma X(2 * mu_state.getDOF() + 1);
                    this->generateSigmaPoints(mu_state, Pk, X);

                    VectorizedMultiState mu_delta(mu_state.getDOF(), 1);
                    mu_delta.setZero();

                    std::vector<Measurement> Z(X.size());
                    std::transform(X.begin(), X.end(), Z.begin(), h);

                    const Measurement meanZ = meanSigmaPoints(Z);

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
                            _MultiState innovation_state;
                            innovation_state.set(K * innovation, mu_state.featuresk.size(), mu_state.featuresk_l.size());
                            mu_state = mu_state + innovation_state;
                    }

                    #ifdef MSCKF_DEBUG_PRINTS
                    std::cout << "[MSCKF_UPDATE] innovation:" << std::endl << innovation << std::endl;
                    std::cout << "[MSCKF_UPDATE] mu_state':" << std::endl << mu_state << std::endl;
                    #endif
            }


            void muSingleState(const _SingleState & state)
            {
                mu_state.statek = state;
            }

            _SingleState muSingleState()
            {
                _SingleState statek_i = mu_state.statek;

                return statek_i;
            }

            void setPkState(const SingleStateCovariance & Pk_i)
            {
                MTK::subblock (Pk, &_MultiState::statek) = Pk_i;
            }

            SingleStateCovariance getPkState()
            {
                SingleStateCovariance Pk_i;
                Pk_i = Pk.block(0, 0, Pk_i.rows(), Pk_i.cols());

                return Pk_i;
            }

            const _MultiState& muState() const
            {
                return mu_state;
            }

            const MultiStateCovariance &getPk() const
            {
                return Pk;
            }

    private:
            /**@brief Sigma Point Calculation for the complete Multi State
            */
            void generateSigmaPoints(const _MultiState &mu, const MultiStateCovariance &sigma, MultiStateSigma &X) const
            {
                    generateSigmaPoints(mu, VectorizedMultiState::Zero(), sigma, X);
            }

            /**@brief Sigma Point Calculation for the complete Multi State
             */
            void generateSigmaPoints(const _MultiState &mu, const VectorizedMultiState &delta,
                                    const MultiStateCovariance &sigma, MultiStateSigma &X) const
            {
                    assert(X.size() == static_cast<unsigned int> (2 * mu_state.getDOF() + 1));

                    Eigen::LLT< MultiStateCovariance > lltOfSigma(sigma); // compute the Cholesky decomposition of A
                    MultiStateCovariance L = lltOfSigma.matrixL(); // retrieve factor L  in the decomposition

                    /*std::cout << "L is of size "<<L.rows()<<" x "<<L.cols()<<"\n";
                    std::cout << ">> L" << std::endl
                                      << L << std::endl
                                      << "<< L" << std::endl;
                     std::cout<<"L*L^T:\n"<< L * L.transpose()<<"\n";*/

                    X[0] = mu + delta;
                    for (register unsigned int i = 1, j = 0; j < mu.getDOF(); ++j)
                    {
                            //std::cout << "L.col(" << j << "): " << L.col(j).transpose() << std::endl;
                            X[i++] = mu + (delta + L.col(j));
                            X[i++] = mu + (delta - L.col(j));
                    }
                    #ifdef MSCKF_DEBUG_PRINTS
                    this->printSigmaPoints<MultiStateSigma>(X);
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

                    #ifdef MSCKF_DEBUG_PRINTS
                    this->printSigmaPoints<SingleStateSigma>(X);
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

            void applyDeltaMultiState(const VectorizedMultiState &delta)
            {
                    MultiStateSigma X(2 * DOF_MULTI_STATE + 1);
                    generateSigmaPoints(mu_state, delta, Pk, X);

                    mu_state = meanSigmaPoints(X);
                    Pk = covSigmaPoints<_MultiState::DOF>(mu_state, X);
            }

            void applyDeltaState(_SingleState &statek_i, SingleStateCovariance &Pk_i, const VectorizedSingleState &delta)
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
                MultiStateSigma X(2 * DOF_MULTI_STATE + 1);
                generateSigmaPoints(mu_state, Pk, X);

                _MultiState muX = meanSigmaPoints(X);

                MultiStateCovariance Pktest = covSigmaPoints<_MultiState::DOF>(muX, X);
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
                #ifdef MSCKF_DEBUG_PRINTS
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

#endif // __MSCKF_HPP_
