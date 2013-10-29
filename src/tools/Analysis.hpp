/**\file Analysis.hpp
 * Header function file and defines
 */

#ifndef _SENSITIVITY_ANALYSIS_HPP_
#define _SENSITIVITY_ANALYSIS_HPP_

#include <iostream> /** IO C++ Standard library */
#include <cmath> /** Math C++ Standard library */
#include <algorithm> /** Algorithm C++ Standard library */
#include <Eigen/Geometry> /** Eigen data type for Matrix, Quaternion, etc... */
#include <Eigen/Core> /** Core methods of Eigen implementation **/
#include <Eigen/Dense> /** for the algebra and transformation matrices **/
#include <Eigen/Cholesky> /** For the Cholesky decomposition **/


//#define SENSITIVITY_ANALYSIS_PRINTS 1

namespace localization	
{
    template <int _StateSize, int _DoF>
    class Analysis
    {
        private:
            std::vector< Eigen::Matrix <double, _StateSize, 1> , Eigen::aligned_allocator < Eigen::Matrix <double, _StateSize, 1> > >vstates;
            std::vector< Eigen::Matrix <double, _StateSize, _StateSize> , Eigen::aligned_allocator < Eigen::Matrix <double, _StateSize, _StateSize> > > vCovariances;
            std::vector< Eigen::Matrix <double, _DoF, 1> , Eigen::aligned_allocator < Eigen::Matrix <double, _DoF, 1> > > vparameters;

        public:

            /**@brief Constructor
             */
            Analysis()
            {
                vstates = std::vector< Eigen::Matrix <double, _StateSize, 1> , Eigen::aligned_allocator < Eigen::Matrix <double, _StateSize, 1> > > (2);
                vCovariances = std::vector< Eigen::Matrix <double, _StateSize, _StateSize> , Eigen::aligned_allocator < Eigen::Matrix <double, _StateSize, _StateSize> > > (2);
                vparameters = std::vector< Eigen::Matrix <double, _DoF, 1> , Eigen::aligned_allocator < Eigen::Matrix <double, _DoF, 1> > > (2);

                for (register unsigned int i=0; i<vstates.size(); ++i)
                {
                    vstates[i].setZero();
                    vCovariances[i].setZero();
                    vparameters[i].setZero();
                }

            }

            Eigen::Matrix <double, _DoF, 1> solve (const  Eigen::Matrix <double, _StateSize, 1> &state, const Eigen::Matrix <double, _StateSize, _StateSize> &Cov,
                                                   const Eigen::Matrix <double, _DoF, 1> &parameter,
                                                   Eigen::Matrix <double, _DoF, 1> &TCov)
            {
                Eigen::Matrix <double, _DoF, 1> Tstate;

                /** Push the new values that arrived **/

                /** New state **/
                vstates[1] = vstates[0];
                vstates[0] = state;

                /** New covariance **/
                vCovariances[1] = vCovariances[0];
                vCovariances[0] = Cov;

                /** New vector of parameters **/
                vparameters[1] = vparameters[0];
                vparameters[0] = parameter;

                #ifdef SENSITIVITY_ANALYSIS_PRINTS
                std::cout<<"[SENSITIVITY_ANALYSIS]vstates\n"<<vstates[0]<<"\n"<<vstates[1]<<"\n";
                std::cout<<"[SENSITIVITY_ANALYSIS]vCovariances\n"<<vCovariances[0]<<"\n"<<vCovariances[1]<<"\n";
                std::cout<<"[SENSITIVITY_ANALYSIS]vparameters\n"<<vparameters[0]<<"\n"<<vparameters[1]<<"\n";
                #endif

                /** Get the coefficient **/
                for (register int i=0; i<_DoF; ++i)
                {
                    double aux = vparameters[0][i] - vparameters[1][i];

                    if (aux!= 0.00)
                    {
                        Tstate[i] = ((vstates[0] - vstates[1])/aux).transpose() /* Cov.inverse()*/ * ((vstates[0] - vstates[1])/aux);
                        TCov[i] = ((vCovariances[0] - vCovariances[1])/aux).determinant() ;/*/ Cov.determinant();*/
                    }
                    else
                    {
                        Tstate[i] = 0.00;
                        TCov[i] = 0.00;
                    }
                    #ifdef SENSITIVITY_ANALYSIS_PRINTS
                    std::cout<<"[SENSITIVITY_ANALYSIS]division_state\n"<<(vstates[0] - vstates[1])/aux<<"\n";
                    std::cout<<"[SENSITIVITY_ANALYSIS]division_cov.determinant\n"<<((vCovariances[0] - vCovariances[1])/aux).determinant()<<"\n";
                    std::cout<<"[SENSITIVITY_ANALYSIS]Cov.inverse\n"<<Cov.inverse()<<"\n";
                    std::cout<<"[SENSITIVITY_ANALYSIS]Cov.determinant\n"<<Cov.determinant()<<"\n";
                    std::cout<<"[SENSITIVITY_ANALYSIS]aux\n"<<aux<<"\n";
                    std::cout<<"[SENSITIVITY_ANALYSIS]Tstate:"<<Tstate[i]<<"\n";
                    std::cout<<"[SENSITIVITY_ANALYSIS]TCov:"<<TCov[i]<<"\n";
                    #endif

                }

                return Tstate;
            }

    };
}
#endif

