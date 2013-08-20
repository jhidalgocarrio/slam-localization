#ifndef _IIR_HPP_
#define _IIR_HPP_

#include <Eigen/Core> /** Core methods of Eigen implementation **/

//#define IIR_DEBUG_PRINTS 1

namespace localization
{
    /**@brief Implements an Infinite Impulse Response filter with order
     * _Order especified as template parameters and _DataDimension is the
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
            Eigen::Matrix<double, _DataDimension, 1> result;

            /** Copy the new data samples to the array for filtering */
            originalData.template block<_DataDimension, _Order>(0,0) = originalData.template block<_DataDimension, _Order>(0,1);
            originalData.col(_Order) = data;

            #ifdef IIR_DEBUG_PRINTS
            std::cout<<"[IIR-FILTER] originalData.cols "<<originalData.cols()<<" filteredData.cols "<<filteredData.cols()<<"\n";
            #endif

            /** Perform the IIR filter with the designed coefficients */
            result = this->iirFilter(bCoeff, aCoeff, originalData, filteredData);

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

	    if ((x.cols() >= _Order+1) && (y.cols() >= _Order+1))
	    {
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

	    }


	    return y.col(i);
	};

    };
}

#endif

