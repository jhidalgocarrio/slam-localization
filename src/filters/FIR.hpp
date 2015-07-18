#ifndef _FIR_HPP_
#define _FIR_HPP_

#include <Eigen/Core> /** Core methods of Eigen implementation **/
#include <Eigen/Cholesky> /** For the Cholesky decomposition **/

//#define FIR_DEBUG_PRINTS 1

namespace localization
{
    /**@brief Implements an Finite Impulse Response filter with order
     * _Order specified as template parameters and _DataDimension is the
     * Dimension of the data to filter (1-D, 2-D, 3-D, etc..) */
    template < unsigned int _Order, unsigned int _DataDimension > class FIR
    {

    public:
	
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW
	
    protected:

        /** Coefficients of a polinomia of order _Order */
        /** FIR filter Coefficients */
        Eigen::Matrix <double, _Order+1, 1> bCoeff;

        /** Data values for original and filtered measurements **/
        Eigen::Matrix <double, _DataDimension, _Order+1> originalData, filteredData;

        /** Variables for the uncertainty propagation through an FIR filter (is your original data has uncertainty) **/
        Eigen::Matrix <double, (_Order+1)*_DataDimension, (_Order+1)*_DataDimension> originalDataCov, filteredDataCov;
        Eigen::Matrix <double, _DataDimension, (_Order+1)*_DataDimension> bMatrix;
        Eigen::Matrix <double, (_Order+1)*_DataDimension, (_Order+1)*_DataDimension> crossCov;


    public:

        /**
         * @brief Constructor
         *
         * @param[in] b the feedforward filter coefficients
         *
         */
        FIR(const Eigen::Matrix <double, _Order+1, 1> &b)
            :bCoeff(b)
        {
            originalData.setZero();
            filteredData.setZero();
            originalDataCov.setZero();
            filteredDataCov.setZero();
            crossCov.setZero();

            /** Form the matrix B of the coefficients (Only used for the uncertainty propagation) **/
            bMatrix.setZero();
            Eigen::Matrix<double, _DataDimension, _Order+1> bSubMatrix;

            for (register int i=0; i<static_cast<int>(_DataDimension); ++i)
            {
                bSubMatrix.row(i) = 1.0 * b.transpose();
            }

            for (register int i=0; i<static_cast<int>(_DataDimension); ++i)
            {
                bMatrix.template block<_DataDimension, _Order+1> (0, i*(_Order+1)) = bSubMatrix;
            }

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

            #ifdef FIR_DEBUG_PRINTS
            std::cout<<"[FIR-FILTER] originalData.cols "<<originalData.cols()<<" filteredData.cols "<<filteredData.cols()<<"\n";
            #endif

            /** Perform the FIR filter with the designed coefficients */
            result = this->firFilter(bCoeff, originalData, filteredData);

            /** Compute the uncertainty propagation if set to true **/
            if (covariance == true)
            {

                /** Check if NaN values in the matrix **/
                if (!base::isnotnan(dataCov))
                {
                    #ifdef FIR_DEBUG_PRINTS
                    std::cout<<"[FIR] dataCov has NaN values\n";
                    #endif

                    dataCov.setZero();
                }

                #ifdef FIR_DEBUG_PRINTS
                std::cout<<"[FIR] dataCov :\n"<<dataCov<<"\n";
                #endif

              /** TO-DO **/

            }

            #ifdef FIR_DEBUG_PRINTS
            std::cout<<"[FIR] Result:\n"<<result<<"\n";
            #endif

            /** Move back one column for the next filter iteration */
            filteredData.template block<_DataDimension, _Order>(0,0) = filteredData.template block<_DataDimension, _Order>(0,1);

            return result;
        }

    protected:

        /**
	* @brief FIR filter
	*
	* Finite Impulse Response (FIR) Filter for an Eigen Vector
	*
	* @author Javier Hidalgo Carrio.
	*
	* @param[in] b array of feedforward filter coefficients
	* @param[in] x array of inputs
	* @param[in,out] y array of outputs
	*
	* @return double with the result y[n]
	*
	*/
	Eigen::Matrix <double, _DataDimension, 1> firFilter (
				    const Eigen::Matrix <double, _Order+1, 1> &b,
				    const Eigen::Matrix <double, _DataDimension, _Order+1> &x,
                                    Eigen::Matrix <double, _DataDimension, _Order+1> &y)
	{
            register int j;
	    register int i = static_cast<int>(_Order);

            Eigen::Matrix <double, _DataDimension, 1> tmpB;
            tmpB.setZero();

            /** Perform the filter **/
            for(j=1; j<static_cast<int>(_Order+1); ++j)
            {
                tmpB += b[j-1]*x.col(i);
                i--;
            }
            tmpB += b[j-1]*x.col(i);

            y.col(_Order) = (tmpB);

	    return y.col(_Order);
	};


        /**
	* @brief FIR filter Cov
	*
	* Finite Impulse Response (FIR) Filter for the Cov Matrix of the data
	*
	* @author Javier Hidalgo Carrio.
	*
	* @param[in] b array of feedforward filter coefficients
	* @param[in] x array of inputs
	* @param[in,out] y array of outputs
	*
	* @return double with the result y[n]
	*
	*/
	Eigen::Matrix <double, _DataDimension, _DataDimension> firFilterCov (
				    const Eigen::Matrix <double, _DataDimension, (_Order+1)*_DataDimension> &bMatrix,
				    const Eigen::Matrix <double, (_Order+1)*_DataDimension, (_Order+1)*_DataDimension> &xCov,
                                    Eigen::Matrix <double, (_Order+1)*_DataDimension, (_Order+1)*_DataDimension> &yCov)
	{

        };

    };
}

#endif

