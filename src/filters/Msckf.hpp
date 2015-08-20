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

/** Base types **/
#include <base/Float.hpp>
#include <base/Eigen.hpp>

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
                    DOF_SINGLE_STATE = _MultiState::SingleState::DOF
            };

            enum
            {
                    SENSOR_DOF = _MultiState::SENSOR_DOF
            };



            typedef typename _MultiState::scalar_type ScalarType;

            /** Types related to Single State **/
            typedef typename _SingleState::vectorized_type VectorizedSingleState;
            typedef Eigen::Matrix<ScalarType, int(_SingleState::DOF), int(_SingleState::DOF)> SingleStateCovariance;
            typedef std::vector<_SingleState> SingleStateSigma;

            /** Types related to Multi State **/
            typedef typename _MultiState::vectorized_type VectorizedMultiState;
            typedef Eigen::Matrix<ScalarType, Eigen::Dynamic, Eigen::Dynamic> MultiStateCovariance;
            typedef std::vector<_MultiState> MultiStateSigma;


        private:

            _MultiState mu_state; /** Mean of the state and sensors pose vector **/
            MultiStateCovariance Pk; /** Covariance of the State and sensor vector **/

        public:
            /**@brief Constructor
             */
            Msckf(const _MultiState &state, const MultiStateCovariance &P0)
                : mu_state(state)
            {
                this->Pk.resize(P0.rows(), P0.cols());
                this->Pk = P0;
            }

            /**@brief Filter prediction step
             */
            template<typename _ProcessModel>
            void predict(_ProcessModel f, const SingleStateCovariance &Q)
            {
                Eigen::Matrix<ScalarType, _SingleState::DOF, 4> Nk;
                Nk = base::NaN<double>() * Eigen::Matrix<ScalarType, _SingleState::DOF, 4>::Identity();
                predict(f, boost::bind(ukfom::id<SingleStateCovariance>, Q), Nk);
            }

            template<typename _ProcessModel, typename _ProcessNoiseCovariance, typename _NullSpaceMatrix>
            void predict(_ProcessModel f, _ProcessNoiseCovariance Q, _NullSpaceMatrix Nk)
            {

                /** Get the current state vector to propagate **/
                _SingleState statek_i = this->mu_state.statek;
                SingleStateCovariance Pk_i = this->Pk.block(0, 0, _SingleState::DOF, _SingleState::DOF);

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
                this->Pk.block(0, 0, _SingleState::DOF, _SingleState::DOF) = Pk_i;

                /*******************************/
                /**  Cross-Covariance Matrix  **/
                /*******************************/

                //Eigen::Matrix<ScalarType, _SingleState::DOF, Eigen::Dynamic> Pkk;
                //Pkk.resize(_SingleState::DOF, _MultiState::SENSOR_DOF * mu_state.sensorsk.size());

                ///** Covariance between state and sensor poses **/
                //Pkk = Pk.block(0, _SingleState::DOF, _SingleState::DOF, _MultiState::SENSOR_DOF*mu_state.sensorsk.size());
                //Pkk = Fk * Pkk;
                //Pk.block(0, _SingleState::DOF, _SingleState::DOF, _MultiState::SENSOR_DOF*mu_state.sensorsk.size()) = Pkk;

                ///** Covariance between sensor poses and state **/
                //Pkk.transpose() = Pk.block(_SingleState::DOF, 0, _MultiState::SENSOR_DOF*mu_state.sensorsk.size(), _SingleState::DOF);
                //Pkk.transpose() = Pkk.transpose() * Fk.transpose();
                //Pk.block(_SingleState::DOF, 0, _MultiState::SENSOR_DOF*mu_state.sensorsk.size(), _SingleState::DOF) = Pkk.transpose();

                #ifdef  MSCKF_DEBUG_PRINTS
                std::cout << "[MSCKF_PREDICT] statek_i(k+1|k):" << std::endl << mu_state.statek << std::endl;
                std::cout << "[MSCKF_PREDICT] Pk(k+1|k):"<< std::endl << Pk_i << std::endl;
                std::cout << "[MSCKF_PREDICT] Process Noise Cov Q(k):"<< std::endl << Q() << std::endl;
                #endif
            }

            template<typename _Measurement, typename _MeasurementModel, typename _MeasurementNoiseCovariance>
            unsigned int update(const _Measurement &z, _MeasurementModel h, _MeasurementNoiseCovariance &R)
            {
                    return update(z, h, R, accept_mahalanobis_distance<ScalarType>);
            }

            template<typename _Measurement, typename _MeasurementModel>
            unsigned int update(const _Measurement &z, _MeasurementModel h,
                        const Eigen::Matrix<ScalarType, Eigen::Dynamic, Eigen::Dynamic> &R)
            {
                    typedef Eigen::Matrix<ScalarType, ukfom::dof<_Measurement>::value, ukfom::dof<_Measurement>::value> measurement_cov;
                    return update(z, h, R, accept_mahalanobis_distance<ScalarType>);
            }

            template<typename _Measurement, typename _MeasurementModel,
                    typename _MeasurementNoiseCovariance, typename _SignificanceTest>
            unsigned int update(const _Measurement &z, _MeasurementModel &h,
                        _MeasurementNoiseCovariance &R, _SignificanceTest mt)
            {
                    typedef Eigen::Matrix<ScalarType, Eigen::Dynamic, 1> VectorXd;
                    typedef Eigen::Matrix<ScalarType, Eigen::Dynamic, Eigen::Dynamic> MatrixXd;

                    MultiStateSigma X(2 * mu_state.getDOF() + 1);
                    generateSigmaPoints(mu_state, Pk, X);

                    std::vector<VectorXd> Z(X.size());
                    std::transform(X.begin(), X.end(), Z.begin(), h);

                    const VectorXd mean_z = meanSigmaPoints(Z);

                    VectorXd innovation = z - mean_z;

                    MatrixXd covXZ = this->crossCovSigmaPoints(mu_state, mean_z, X, Z);

                    MatrixXd H = covXZ.transpose() * Pk.inverse(); //H = Pzx * (Pk)^-1

                    //std::cout<<"[MSCKF_UKF_UPDATE] H size "<<H.rows()<<" x "<<H.cols()<<"\n";
                    const unsigned int number_outliers = removeOutliers (innovation, H, Pk, R, mt, 2);

                    #ifdef MSCKF_DEBUG_PRINTS
                    std::cout<<"[MSCKF_UKF_UPDATE] H size "<<H.rows()<<" x "<<H.cols()<<"\n";
                    std::cout<<"[MSCKF_UKF_UPDATE] H \n"<<H<<"\n";
                    std::cout<<"[MSCKF_UKF_UPDATE] Pk size "<<Pk.rows()<<" x "<<Pk.cols()<<"\n";
                    std::cout<<"[MSCKF_UKF_UPDATE] R size "<<R.rows()<<" x "<<R.cols()<<"\n";
                    #endif

                    if (innovation.rows() > 0)
                    {
                        #ifdef MSCKF_DEBUG_PRINTS
                        std::cout << "[MSCKF_UKF_UPDATE] innovation size "<<innovation.rows()<<" x "<<innovation.cols()<<"\n";
                        std::cout << "[MSCKF_UKF_UPDATE] innovation\n"<<innovation<<"\n";
                        #endif

                        reduceDimension (innovation, H, R);
                        const MatrixXd S = H * Pk * H.transpose() + R;
                        const MatrixXd K = Pk * H.transpose() * S.inverse();
                        #ifdef MSCKF_DEBUG_PRINTS
                        std::cout << "[MSCKF_UKF_UPDATE] innovation\n"<<innovation<<"\n";
                        #endif

                        Pk -= K * S * K.transpose();
                        this->applyDelta(K * innovation);

                        #ifdef MSCKF_DEBUG_PRINTS
                        std::cout<<"[MSCKF_UKF_UPDATE] K "<<K.rows() <<" x "<<K.cols()<<"\n";
                        std::cout<<"[MSCKF_UKF_UPDATE] Pk "<<Pk.rows() <<" x "<<Pk.cols()<<"\n";
                        #endif
                    }

                    #ifdef MSCKF_DEBUG_PRINTS
                    std::cout << "[MSCKF_UKF_UPDATE] mu_state':" << std::endl << mu_state << std::endl;
                    std::cout << "[MSCKF_UKF_UPDATE] Pk':" << std::endl << Pk << std::endl;
                    #endif

                    return number_outliers;
            }

            /**@brief update
             *
             * EKF update
             *
             */
            template<typename _Measurement, typename _MeasurementModel>
            unsigned int update(const _Measurement &z, _MeasurementModel h,
                    Eigen::Matrix<ScalarType, Eigen::Dynamic, Eigen::Dynamic> &H,
                    Eigen::Matrix<ScalarType, Eigen::Dynamic, Eigen::Dynamic> &R)
            {
                    return update(z, h, H, R, accept_mahalanobis_distance<ScalarType>);
            }

            /**@brief update
             *
             * EKF update
             *
             */
            template<typename _Measurement, typename _MeasurementModel,
                    typename _MeasurementNoiseCovariance, typename _SignificanceTest>
            unsigned int update(const _Measurement &z, _MeasurementModel &h,
                        Eigen::Matrix<ScalarType, Eigen::Dynamic, Eigen::Dynamic> &H,
                        _MeasurementNoiseCovariance &R, _SignificanceTest mt)
            {
                    typedef Eigen::Matrix<ScalarType, Eigen::Dynamic, 1> VectorXd;
                    typedef Eigen::Matrix<ScalarType, Eigen::Dynamic, Eigen::Dynamic> MatrixXd;

                    const VectorXd mean_z = h(this->mu_state, H);

                    VectorXd innovation = z - mean_z;

                    const unsigned int number_outliers = removeOutliers (innovation, H, Pk, R, mt, 2);
                    #ifdef MSCKF_DEBUG_PRINTS
                    std::cout<<"[MSCKF_EKF_UPDATE] H size "<<H.rows()<<" x "<<H.cols()<<"\n";
                    std::cout<<"[MSCKF_EKF_UPDATE] H \n"<<H<<"\n";
                    std::cout<<"[MSCKF_EKF_UPDATE] Pk size "<<Pk.rows()<<" x "<<Pk.cols()<<"\n";
                    std::cout<<"[MSCKF_EKF_UPDATE] R size "<<R.rows()<<" x "<<R.cols()<<"\n";
                    #endif

                    if (innovation.rows() > 0)
                    {
                        #ifdef MSCKF_DEBUG_PRINTS
                        std::cout << "[MSCKF_EKF_UPDATE] innovation size "<<innovation.rows()<<" x "<<innovation.cols()<<"\n";
                        std::cout << "[MSCKF_EKF_UPDATE] innovation\n"<<innovation<<"\n";
                        #endif

                        reduceDimension (innovation, H, R);
                        const MatrixXd S = H * Pk * H.transpose() + R;
                        const MatrixXd K = Pk * H.transpose() * S.inverse();
                        #ifdef MSCKF_DEBUG_PRINTS
                        std::cout << "[MSCKF_EKF_UPDATE] innovation\n"<<innovation<<"\n";
                        #endif

                        Pk -= K * S * K.transpose();
                        this->mu_state = this->mu_state + K * innovation;

                        #ifdef MSCKF_DEBUG_PRINTS
                        std::cout<<"[MSCKF_EKF_UPDATE] K "<<K.rows() <<" x "<<K.cols()<<"\n";
                        std::cout<<"[MSCKF_EKF_UPDATE] Pk "<<Pk.rows() <<" x "<<Pk.cols()<<"\n";
                        #endif

                        base::guaranteeSPD(Pk);
                    }

                    #ifdef MSCKF_DEBUG_PRINTS
                    std::cout << "[MSCKF_EKF_UPDATE] mu_state':" << std::endl << mu_state << std::endl;
                    std::cout << "[MSCKF_EKF_UPDATE] Pk':" << std::endl << Pk << std::endl;
                    #endif

                    return number_outliers;
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

            void setPkSingleState(const SingleStateCovariance & Pk_i)
            {
                this->Pk.block(0, 0, _SingleState::DOF, _SingleState::DOF) = Pk_i;
            }

            SingleStateCovariance getPkSingleState()
            {
                SingleStateCovariance Pk_i;
                Pk_i = Pk.block(0, 0, Pk_i.rows(), Pk_i.cols());

                return Pk_i;
            }

            const _MultiState& muState() const
            {
                return mu_state;
            }

            _MultiState& muState()
            {
                return mu_state;
            }

            const MultiStateCovariance &getPk() const
            {
                return Pk;
            }

            void setPk (const MultiStateCovariance &Pk_i)
            {
                Pk.resize(Pk_i.rows(), Pk_i.cols());
                this->Pk = Pk_i;
            }

    private:
            /**@brief Sigma Point Calculation for the complete Multi State
            */
            void generateSigmaPoints(const _MultiState &mu, const MultiStateCovariance &sigma, MultiStateSigma &X) const
            {
                    generateSigmaPoints(mu, VectorizedMultiState::Zero(mu.getDOF(), 1), sigma, X);
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

            // manifold mean for single state
            _SingleState meanSigmaPoints(const std::vector<_SingleState> &X) const
            {
                    _SingleState reference = X[0];
                    typename _SingleState::vectorized_type mean_delta;
                    const static std::size_t max_it = 10000;

                    std::size_t i = 0;
                    do {
                            mean_delta.setZero();
                            for (typename std::vector<_SingleState>::const_iterator Xi = X.begin(); Xi != X.end(); ++Xi)
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

            // manifold mean for multi state
            _MultiState meanSigmaPoints(const std::vector<_MultiState> &X) const
            {
                    _MultiState reference = X[0];
                    typename _MultiState::vectorized_type mean_delta;
                    mean_delta.resize(reference.getDOF(), 1);
                    const static std::size_t max_it = 10000;

                    std::size_t i = 0;
                    do {
                            mean_delta.setZero();
                            for (typename std::vector<_MultiState>::const_iterator Xi = X.begin(); Xi != X.end(); ++Xi)
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

            /*@brief covariance of sigma point when using static size manifold vector
             */
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

            /*@brief covariance of sigma points when using the _MultiState
             */
            Eigen::Matrix<ScalarType, Eigen::Dynamic, Eigen::Dynamic>
            covSigmaPoints(const _MultiState &mean, const std::vector<_MultiState> &V) const
            {
                    typedef Eigen::Matrix<ScalarType, Eigen::Dynamic, Eigen::Dynamic> CovMat;
                    typedef Eigen::Matrix<ScalarType, Eigen::Dynamic, 1> CovCol;

                    CovMat c(CovMat::Zero(mean.getDOF(), mean.getDOF()));

                    for (typename std::vector< _MultiState >::const_iterator Vi = V.begin(); Vi != V.end(); ++Vi)
                    {
                            CovCol d = *Vi - mean;
                            c += d * d.transpose();
                    }

                    return 0.5 * c;
            }

            /*@brief covariance of sigma points for the dynamic size measurement vector
             */
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
            crossCovSigmaPoints(const _State &mean_x, const _Measurement &mean_z,
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
                                    c += (*Xi - mean_x) * (*Zi - mean_z).transpose();
                            }
                    }

                    return 0.5 * c;
            }

            template<typename _State, typename _SigmaPoints>
            Eigen::Matrix<ScalarType, Eigen::Dynamic, Eigen::Dynamic>
            crossCovSigmaPoints(const _State &mean_x, const Eigen::Matrix<ScalarType, Eigen::Dynamic, 1> &mean_z,
                                const _SigmaPoints &X, const std::vector< Eigen::Matrix<ScalarType, Eigen::Dynamic, 1> > &Z) const
            {
                    assert(X.size() == Z.size());

                    typedef Eigen::Matrix<ScalarType, Eigen::Dynamic, Eigen::Dynamic> CrossCov;

                    CrossCov c(mu_state.getDOF(), mean_z.size());
                    c.setZero();

                    {
                            typename _SigmaPoints::const_iterator Xi = X.begin();
                            typename std::vector< Eigen::Matrix<ScalarType, Eigen::Dynamic, 1> >::const_iterator Zi = Z.begin();
                            for (;Zi != Z.end(); ++Xi, ++Zi)
                            {
                                    c += (*Xi - mean_x) * (*Zi - mean_z).transpose();
                            }
                    }

                    return 0.5 * c;
            }

            void applyDelta(const VectorizedMultiState &delta)
            {
                    MultiStateSigma X(2 * mu_state.getDOF() + 1);
                    generateSigmaPoints(mu_state, delta, Pk, X);

                    mu_state = meanSigmaPoints(X);
                    Pk = covSigmaPoints(mu_state, X);
            }

            void applyDelta(_SingleState &statek_i, SingleStateCovariance &Pk_i, const VectorizedSingleState &delta)
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

            static void removeRow(Eigen::Matrix<ScalarType, Eigen::Dynamic, 1> &vector, unsigned int rowToRemove)
            {
                unsigned int numRows = vector.rows()-1;

                if( rowToRemove < numRows )
                    vector.block(rowToRemove,0,numRows-rowToRemove, 1) =
                        vector.block(rowToRemove+1,0,numRows-rowToRemove, 1);

                vector.conservativeResize(numRows);
            }

            static void removeRow(Eigen::Matrix<ScalarType, Eigen::Dynamic, Eigen::Dynamic> &matrix, unsigned int rowToRemove)
            {
                unsigned int numRows = matrix.rows()-1;
                unsigned int numCols = matrix.cols();

                if( rowToRemove < numRows )
                    matrix.block(rowToRemove,0,numRows-rowToRemove,numCols) =
                        matrix.block(rowToRemove+1,0,numRows-rowToRemove,numCols);

                matrix.conservativeResize(numRows,numCols);
            }

            static void removeColumn(Eigen::Matrix<ScalarType, Eigen::Dynamic, Eigen::Dynamic>& matrix, unsigned int colToRemove)
            {
                unsigned int numRows = matrix.rows();
                unsigned int numCols = matrix.cols()-1;

                if( colToRemove < numCols )
                    matrix.block(0,colToRemove,numRows,numCols-colToRemove) =
                        matrix.block(0,colToRemove+1,numRows,numCols-colToRemove);

                matrix.conservativeResize(numRows,numCols);
            }


            template <typename _SignificanceTest>
            unsigned int removeOutliers(Eigen::Matrix<ScalarType, Eigen::Dynamic, 1> &innovation,
                    Eigen::Matrix<ScalarType, Eigen::Dynamic, Eigen::Dynamic> &h_matrix,
                    Eigen::Matrix<ScalarType, Eigen::Dynamic, Eigen::Dynamic> &p_matrix,
                    Eigen::Matrix<ScalarType, Eigen::Dynamic, Eigen::Dynamic> &r_matrix,
                    _SignificanceTest mt,
                    const unsigned int dof)
            {
                /** Information matrix **/
                Eigen::Matrix<ScalarType, Eigen::Dynamic, Eigen::Dynamic> information  =
                    ((h_matrix * p_matrix * h_matrix.transpose()) + r_matrix).inverse();

                /** Process the innovation by block of dof size **/
                unsigned int number_outliers = 0;
                register unsigned int i=0;
                while (i<innovation.size()/dof)
                {
                    const ScalarType mahalanobis2 = (innovation.block(dof*i, 0, dof, 1).transpose() * information.block(dof*i, dof*i, dof, dof) * innovation.block(dof*i, 0, dof, 1))[0];
                    //std::cout<<"feature["<<i<<"]:\n"<<innovation.block(dof*i, 0, dof, 1) <<"\n";
                    if (!mt(mahalanobis2, dof))
                    {
                        //std::cout<<"OUTLIER!!\n";
                        removeRow(innovation, dof*i); removeRow(innovation, (dof*i)+1);
                        removeRow(r_matrix, dof*i); removeColumn(r_matrix, dof*i);
                        removeRow(r_matrix, (dof*i)+1); removeColumn(r_matrix, (dof*i)+1);
                        removeRow(h_matrix, dof*i);removeRow(h_matrix, (dof*i)+1);
                        number_outliers++;
                    }
                    else
                    {
                        i++;
                    }
                }

                return number_outliers;
            }

            void reduceDimension(Eigen::Matrix<ScalarType, Eigen::Dynamic, 1> &innovation,
                    Eigen::Matrix<ScalarType, Eigen::Dynamic, Eigen::Dynamic> &h_matrix,
                    Eigen::Matrix<ScalarType, Eigen::Dynamic, Eigen::Dynamic> &r_matrix)
            {
                typedef Eigen::Matrix<ScalarType, Eigen::Dynamic, Eigen::Dynamic> MatrixXd;

                Eigen::HouseholderQR<MatrixXd> qr(h_matrix);
                MatrixXd Q = qr.householderQ();
                MatrixXd R = qr.matrixQR().template triangularView<Eigen::Upper>();

                MatrixXd thinQ(MatrixXd::Identity(h_matrix.rows(), h_matrix.cols()));
                thinQ = qr.householderQ() * thinQ;

                /** Reduced H matrix **/
                h_matrix = R.block(0, 0, this->mu_state.getDOF(), this->mu_state.getDOF());

                /** Reduced innovation **/
                innovation = thinQ.transpose() * innovation;

                /** Reduced noise matrix **/
                r_matrix = thinQ.transpose() * r_matrix * thinQ;

                return;
            }

    public:
            void checkSigmaPoints()
            {
                MultiStateSigma X(2 * mu_state.getDOF() + 1);
                generateSigmaPoints(mu_state, Pk, X);

                _MultiState muX = meanSigmaPoints(X);

                MultiStateCovariance Pktest = covSigmaPoints(muX, X);
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
