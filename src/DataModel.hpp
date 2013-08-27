/**\file Datamodel.hpp
 * Header function file and defines
 */

#ifndef _DATAMODEL_HPP_
#define _DATAMODEL_HPP_

#include <iostream>

#include <Eigen/Core> /** Core methods of Eigen implementation **/
#include <Eigen/Dense> /** for the algebra and transformation matrices **/
#include "Configuration.hpp" /** For the localization framework constant and configuration values **/

namespace localization	
{

    /** Class for representing a 3D slip, linear or contact angle velocity vector and
     * its uncertainty in the estimation by Weighted Least-Squares**/
    template < typename _Scalar, int _DIM >
    class DataModel
    {
    
    public:
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW
	Eigen::Matrix<_Scalar, _DIM, 1> data; //! instantaneous slip velocity vector
	Eigen::Matrix<_Scalar, _DIM, _DIM> Cov; //! Covariance matrix

    public:
	
        /*@brief constructor
         */
	DataModel()
        {
            data.setZero();
            Cov = localization::ZERO_UNCERTAINTY * Eigen::Matrix<_Scalar, _DIM, _DIM>::Identity();
        }

	DataModel( const Eigen::Matrix<_Scalar, _DIM, 1> &initData, Eigen::Matrix<_Scalar, _DIM, _DIM> &initCov)
            : data(initData), Cov(initCov)
        {
        }

	inline int size ()
        {
            return (data.size());
        }

	void fusion(const DataModel<_Scalar, _DIM> &data2)
        {
            DataModel &data1 (*this);

            Eigen::Matrix <_Scalar, _DIM, _DIM> P;

            P = (data1.Cov.inverse() + data2.Cov.inverse()).inverse();
            data1.data = P*(data1.Cov.inverse() * data1.data + data2.Cov.inverse() * data2.data);

            data1.Cov = P;

            return;
        }

	void safeFusion(const DataModel<_Scalar, _DIM> &data2)
        {
            DataModel &data1 (*this);

            /** Variables for the safe fusion **/
            Eigen::Matrix<_Scalar, _DIM, _DIM> I1, I2;
            Eigen::Matrix<_Scalar, _DIM, _DIM> U1, sqrtD1, isqrtD1;
            Eigen::Matrix<_Scalar, _DIM, _DIM> U2, D2;
            Eigen::Matrix<_Scalar, _DIM, _DIM> D3;

            Eigen::Matrix<_Scalar, _DIM, 1> result;
            Eigen::Matrix<_Scalar, _DIM, 1> datatrans1, datatrans2;

            D3.setZero();

            I1 = data1.Cov.inverse();
            I2 = data2.Cov.inverse();

            /**
            * Single Value Decomposition
            */
            Eigen::JacobiSVD < Eigen::MatrixXd > svdOfI1(I1, Eigen::ComputeThinU);

            sqrtD1 = svdOfI1.singularValues().array().cwiseSqrt().matrix().asDiagonal();
            U1 = svdOfI1.matrixU();

        //     std::cout<<"sqrt(D1):\n"<<sqrtD1<<"\n";
        //     std::cout<<"U1:\n"<<U1<<"\n";

            isqrtD1 = sqrtD1.inverse();

        //     std::cout<<"sqrt(D1)^-1:\n"<<isqrtD1<<"\n";

            I2 = isqrtD1 * U1.transpose() * I2 * U1 * isqrtD1;

            Eigen::JacobiSVD <Eigen::MatrixXd > svdOfI2(I2, Eigen::ComputeThinU);

            D2 = svdOfI2.singularValues().array().matrix().asDiagonal();
            U2 = svdOfI2.matrixU();

        //     std::cout<<"D2:\n"<<D2<<"\n";

            Eigen::Matrix<_Scalar, 3, 3> T;

            T = U2.transpose() * sqrtD1 * U1;

            /** Safe transformation **/
            datatrans1 = T * data1.data;
            datatrans2 = T * data2.data;

            for (int i = 0; i<data1.data.size(); i++)
            {
                if (D2(i,i) < 1.0)
                {
                    result[i] = datatrans1[i];
                    D3(i,i) = 1.0;
                }
                else
                {
                    result[i] = datatrans2[i];
                    D3(i,i) =  D2(i,i);
                }

            }

            data1.data = T.inverse() * result;
            data1.Cov = T.inverse() * D3.inverse() * T.inverse().transpose();

        }

	DataModel operator+(const DataModel<_Scalar, _DIM> &data2) const
        {
            const DataModel &data1 (*this);
            DataModel result;

            result.data = data1.data + data2.data;
            result.Cov = data1.Cov + data2.Cov;

            return result;
        }

	DataModel operator-(const DataModel<_Scalar, _DIM> &data2) const
        {
            const DataModel &data1 (*this);
            DataModel result;

            result.data = data1.data - data2.data;
            result.Cov = data1.Cov + data2.Cov;

            return result;
        }

	DataModel& operator=(const DataModel<_Scalar, _DIM> &dmodel)
        {
            data = dmodel.data;
            Cov = dmodel.Cov;

            return *this;
        }

        /** Default std::cout function
        */
        friend std::ostream & operator<<(std::ostream &out, const  DataModel<_Scalar, _DIM> &m1)
        {
            out <<"\n" << m1.data << "\n";
            out <<"\n" << m1.Cov << "\n";

            return out;
        }

    };
}

#endif
