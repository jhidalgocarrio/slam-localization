/**\file sckf.hpp
 * Header function file and defines
 */

#ifndef _SCKF_HPP_
#define _SCKF_HPP_

#include <iostream>
#include <Eigen/Geometry> /**< Eigen data type for Matrix, Quaternion, etc... */


namespace localization	
{
	class sckf
	{
		public: 
			/**
			* Print a welcome to stdout
			* \return nothing
			*/
			void welcome();
	};

} // end namespace localization

#endif // 
