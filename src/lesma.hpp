/**\file lesma.hpp
 * Header function file and defines
 */

#ifndef _LESMA_HPP_
#define _LESMA_HPP_

#include <iostream>
#include <Eigen/Geometry> /**< Eigen data type for Matrix, Quaternion, etc... */


namespace localization	
{
	class lesma
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
