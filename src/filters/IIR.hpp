#ifndef _IIR_HPP_
#define _IIR_HPP_

#include <Eigen/Core> /** Core methods of Eigen implementation **/
#include <Eigen/Cholesky> /** For the Cholesky decomposition **/

//#define IIR_DEBUG_PRINTS 1

namespace localization
{
    /**@brief Implements an Infinite Impulse Response filter with order
     * _Order specified as template parameters and _DataDimension is the
     * Dimension of the data to filter (1-D, 2-D, 3-D, etc..) */
    template < unsigned int _Order, unsigned int _DataDimension > class IIR
    {

    public:
	
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW
	
    protected:

        /** Coefficients of a polinomia of order _Order */
        /** IIR filter Coefficients */
        Eigen::Matrix <double, _Order+1, 1> bCoeff, aCoeff;

        /** Data values for original and filtered measurements **/
        Eigen::Matrix <double, _DataDimension, _Order+1> originalData, filteredData;

        /** Variables for the uncertainty propagation through an IIR filter (is your original data has uncertainty) **/
        Eigen::Matrix <double, (_Order+1)*_DataDimension, (_Order+1)*_DataDimension> originalDataCov, filteredDataCov;
        Eigen::Matrix <double, _DataDimension, (_Order+1)*_DataDimension> bMatrix, aMatrix;
        Eigen::Matrix <double, (_Order+1)*_DataDimension, (_Order+1)*_DataDimension> crossCov;


    public:

        /**
         * @brief Constructor
         *
         * @param[in] b the feedforward filter coefficients
         * @param[in] a the feedback filter coefficients
         *
         */
        IIR(const Eigen::Matrix <double, _Order+1, 1> &b, const Eigen::Matrix <double, _Order+1, 1> &a)
            :bCoeff(b), aCoeff(a)
        {
            originalData.setZero();
            filteredData.setZero();
            originalDataCov.setZero();
            filteredDataCov.setZero();
            crossCov.setZero();

            /** Form the matrix A and B of the coefficients (Only used for the uncertainty propagation) **/
            aMatrix.setZero(); bMatrix.setZero();
            Eigen::Matrix<double, _DataDimension, _Order+1> aSubMatrix, bSubMatrix;

            for (register int i=0; i<static_cast<int>(_DataDimension); ++i)
            {
                aSubMatrix.row(i) = 1.0/a[0] * a.transpose(); //! 1 x _Order+1 row
                bSubMatrix.row(i) = 1.0/a[0] * b.transpose();
            }

            for (register int i=0; i<static_cast<int>(_DataDimension); ++i)
            {
                aMatrix.template block<_DataDimension, _Order+1> (0, i*(_Order+1)) = aSubMatrix;
                bMatrix.template block<_DataDimension, _Order+1> (0, i*(_Order+1)) = bSubMatrix;
            }

            /** The element for the a[0] is a zero matrix, because this is the value it is calculated by the IIR filter **/
            aMatrix.template block<_DataDimension, _DataDimension>(0,0) = Eigen::Matrix<double, _DataDimension, _DataDimension>::Zero();

        }

         /** Peform the filter calling the protected method of the class
         *
         * @param[in] data vector of dimension _DataDimension of the new sample.
         *
         * @return filtered data
         *
         */
        Eigen::Matrix<double, _DataDimension, 1> perform (const Eigen::Matrix <double, _DataDimension, 1> &data)
        {
            Eigen::Matrix <double, _DataDimension, _DataDimension> dataCov; dataCov.setZero();
            return this-> perform (data, dataCov, false);
        }

        /** Peform the filter calling the protected method of the class
         *
         * @param[in] data vector of dimension _DataDimension of the new sample.
         *
         * @return filtered data
         *
         */
        Eigen::Matrix<double, _DataDimension, 1> perform (const Eigen::Matrix <double, _DataDimension, 1> &data,
                                                        Eigen::Matrix <double, _DataDimension, _DataDimension> &dataCov,
                                                        const bool covariance = true)
        {
            Eigen::Matrix<double, _DataDimension, 1> result;

            /** Copy the new data samples to the array for filtering */
            originalData.template block<_DataDimension, _Order>(0,0) = originalData.template block<_DataDimension, _Order>(0,1);
            originalData.col(_Order) = data;

            #ifdef IIR_DEBUG_PRINTS
            std::cout<<"[IIR-FILTER] originalData.cols "<<originalData.cols()<<" filteredData.cols "<<filteredData.cols()<<"\n";
            #endif

            /** Perform the IIR filter with the designed coefficients */
            result = this->iirFilter(bCoeff, aCoeff, originalData, filteredData);

            /** Compute the uncertainty propagation if set to true **/
            if (covariance == true)
            {

                /** Check if NaN values in the matrix **/
                if (!base::isnotnan(dataCov))
                {
                    #ifdef IIR_DEBUG_PRINTS
                    std::cout<<"[IIR] dataCov has NaN values\n";
                    #endif

                    dataCov.setZero();
                }

                #ifdef IIR_DEBUG_PRINTS
                std::cout<<"[IIR] dataCov :\n"<<dataCov<<"\n";
                #endif

                /** Copy the cov matrix for the data samples **/
                originalDataCov.template block<_Order*_DataDimension, _Order*_DataDimension>(0,0) =
                    originalDataCov.template block<_Order*_DataDimension, _Order*_DataDimension> (_DataDimension, _DataDimension);
                originalDataCov.template block<_DataDimension, _DataDimension>(_Order*_DataDimension, _Order*_DataDimension) = dataCov;

                crossCov.template block<_Order*_DataDimension, _Order*_DataDimension>(0,0) =
                    crossCov.template block<_Order*_DataDimension, _Order*_DataDimension> (_DataDimension, _DataDimension);
                crossCov.template block<_DataDimension, _DataDimension>(_Order*_DataDimension, _Order*_DataDimension) = dataCov  - (result * result.transpose());

                dataCov = this->iirFilterCov(bMatrix, aMatrix, originalDataCov, filteredDataCov);

                #ifdef IIR_DEBUG_PRINTS
                std::cout<<"[IIR] ResultCov:\n"<<dataCov<<"\n";
                #endif

                /** Move back one block for the next filter iteration **/
                filteredDataCov.template block<_Order*_DataDimension, _Order*_DataDimension> (0,0) =
                    filteredDataCov.template block<_Order*_DataDimension, _Order*_DataDimension> (_DataDimension, _DataDimension);
            }

            #ifdef IIR_DEBUG_PRINTS
            std::cout<<"[IIR] Result:\n"<<result<<"\n";
            #endif

            /** Move back one column for the next filter iteration */
            filteredData.template block<_DataDimension, _Order>(0,0) = filteredData.template block<_DataDimension, _Order>(0,1);

            return result;
        }

    protected:

        /**
	* @brief IIR filter
	*
	* Infinite Impulse Response (IIR) Filter for an Eigen Vector
	*
	* @author Javier Hidalgo Carrio.
	*
	* @param[in] norder integer number with the filter order
	* @param[in] b array of feedforward filter coefficients
	* @param[in] a array of feedback filter coefficients
	* @param[in] x array of inputs
	* @param[in,out] y array of outputs
	*
	* @return double with the result y[n]
	*
	*/
	Eigen::Matrix <double, _DataDimension, 1> iirFilter (
				    const Eigen::Matrix <double, _Order+1, 1> &b,
                                    const Eigen::Matrix <double, _Order+1, 1> &a,
				    const Eigen::Matrix <double, _DataDimension, _Order+1> &x,
                                    Eigen::Matrix <double, _DataDimension, _Order+1> &y)
	{
            register int j;
	    register int i = static_cast<int>(_Order);

            Eigen::Matrix <double, _DataDimension, 1> tmpA, tmpB;
            tmpB.setZero(); tmpA.setZero();

            /** Perform the filter **/
            for(j=1; j<static_cast<int>(_Order+1); ++j)
            {
                tmpB += b[j-1]*x.col(i);
                tmpA += a[j]*y.col(i-1);

                i--;
            }
            tmpB += b[j-1]*x.col(i);

            y.col(_Order) = 1.0/a[0] * (tmpB - tmpA);

	    return y.col(_Order);
	};


        /**
	* @brief IIR filter Cov
	*
	* Infinite Impulse Response (IIR) Filter for the Cov Matrix of the data
	*
	* @author Javier Hidalgo Carrio.
	*
	* @param[in] norder integer number with the filter order
	* @param[in] b array of feedforward filter coefficients
	* @param[in] a array of feedback filter coefficients
	* @param[in] x array of inputs
	* @param[in,out] y array of outputs
	*
	* @return double with the result y[n]
	*
	*/
	Eigen::Matrix <double, _DataDimension, _DataDimension> iirFilterCov (
				    const Eigen::Matrix <double, _DataDimension, (_Order+1)*_DataDimension> &bMatrix,
                                    const Eigen::Matrix <double, _DataDimension, (_Order+1)*_DataDimension> &aMatrix,
				    const Eigen::Matrix <double, (_Order+1)*_DataDimension, (_Order+1)*_DataDimension> &xCov,
                                    Eigen::Matrix <double, (_Order+1)*_DataDimension, (_Order+1)*_DataDimension> &yCov)
	{
            //Eigen::Matrix <double, (_Order+1)*_DataDimension, (_Order+1)*_DataDimension> xSigma = xCov.llt().matrixL();
            //Eigen::Matrix <double, (_Order+1)*_DataDimension, (_Order+1)*_DataDimension> ySigma = yCov.llt().matrixL();

            #ifdef IIR_DEBUG_PRINTS
            std::cout<<"[IIR-FILTER] bMatrix is\n"<<bMatrix<<"\naMatrix is\n"<<aMatrix<<"\n";
            std::cout<<"[IIR-FILTER] xCov is\n"<<xCov<<"\nyCov is\n"<<yCov<<"\n";
            std::cout<<"[IIR-FILTER] xyCov is\n"<<crossCov<<"\n";
            #endif

            yCov.template block<_DataDimension, _DataDimension> (_Order*_DataDimension, _Order*_DataDimension) =
                bMatrix * xCov * bMatrix.transpose() + aMatrix * yCov * aMatrix.transpose() - bMatrix *  crossCov *  aMatrix.transpose() - aMatrix * crossCov * bMatrix.transpose();

	    return yCov.template block<_DataDimension, _DataDimension> (_Order*_DataDimension, _Order*_DataDimension);
	};

    };
}

#endif

