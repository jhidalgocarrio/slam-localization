#define BOOST_TEST_MODULE template_for_test_test
#include <boost/test/included/unit_test.hpp>

/** Library **/
#include <rover_localization/Configuration.hpp> /** Constant values of the library */
#include <rover_localization/Util.hpp> /** Util function for localization **/
#include <rover_localization/DataModel.hpp> /** Data Model using Gaussian pdf **/

/** Eigen **/
#include <Eigen/Core> /** Core */

/** Standard libs **/
#include <iostream>
#include <vector>

bool increment(std::vector<int> & vector, int k)
{
  for (unsigned int i = 0; i < vector.size(); ++i)
  {
    int j = vector[i] + 1;
    if (j <= k) {
      vector[i] = j;
      return true; 
    }
    else vector[i] = 0;
    // and carry a 1 by looping back again
  }
  return false;
}

BOOST_AUTO_TEST_CASE( DATAMODEL )
{
   	localization::DataModel data1, data2;
	data1 = localization::DataModel(3);
	data2 = localization::DataModel(3);

	data1.data << 0.0124889, 0.00171945, -0.0138983;
	data2.data << 0.0168381, 0.000632167, -0.0235605;
	
// 	data1.data << 1.0, 1.0, 1.0;
// 	data2.data << 1.0, 1.0, 1.0;
	
//   	data2.Cov = 0.5 * data1.Cov;
	
	std::cout<<"data1 is: "<<data1<<"\n";
	std::cout<<"Inverse of data1 Cov:\n"<<data1.Cov.inverse()<<"\n";
	std::cout<<"data2 is: "<<data2<<"\n";
	std::cout<<"Inverse of data2 Cov:\n"<<data2.Cov.inverse()<<"\n";
	
	localization::DataModel data3;
	data3 = data1 + data2;
	
	std::cout<<"*** OPERATOR + ***\n";
	std::cout<<"data1 + data2 is in data3: "<<data3<<"\n";
	std::cout<<"Inverse of data3 Cov:\n"<<data3.Cov.inverse()<<"\n";
	
	data3 = data1 - data2;
	
	std::cout<<"*** OPERATOR - ***\n";
	std::cout<<"data1 - data2 is in data3: "<<data3<<"\n";
	std::cout<<"Inverse of data3 Cov:\n"<<data3.Cov.inverse()<<"\n";
	
	data3 = data2 - data1;
	
	std::cout<<"data2 - data1 is in data3: "<<data3<<"\n";
	std::cout<<"Inverse of data3 Cov:\n"<<data3.Cov.inverse()<<"\n";
	
	std::cout<<"*** FUSION *** \n";
	data3 = data1;
	data3.fusion(data2);
	
	std::cout<<"data1.fusion(data2) is in data3: "<<data3<<"\n";
	std::cout<<"Inverse of data3 Cov:\n"<<data3.Cov.inverse()<<"\n";
	
	std::cout<<"*** SAFE FUSION *** \n";
	data3 = data1;
	data3.safeFusion(data2);
	
	std::cout<<"data1.safeFusion(data2) is in data3: "<<data3<<"\n";
	std::cout<<"Inverse of data3 Cov:\n"<<data3.Cov.inverse()<<"\n";
	
	Eigen::Matrix<double, localization::NUMAXIS, 1> v1, v2;
	double offset = 18.26;
	Eigen::Matrix<double, localization::NUMAXIS, 1> euler;
	Eigen::Quaternion <double> attitude; /** Initial attitude in case no port in orientation is connected **/
	Eigen::Quaternion <double> q_weight; /** Identity if equaly weight distributed of rover **/
	
	v1 << 0.0, 0.0, 1.0;
	
 	euler<<0.00, 1.29348, .00;
	
	/** Set the initial attitude with the Yaw provided from the initial pose **/
	attitude = Eigen::Quaternion <double> (Eigen::AngleAxisd(euler[2] * localization::D2R, Eigen::Vector3d::UnitZ())*
	Eigen::AngleAxisd(euler[1] * localization::D2R, Eigen::Vector3d::UnitY()) *
	Eigen::AngleAxisd(euler[0] * localization::D2R, Eigen::Vector3d::UnitX()));
	
	q_weight = Eigen::Quaternion <double> (Eigen::AngleAxisd(0.00 * localization::D2R, Eigen::Vector3d::UnitZ())*
	Eigen::AngleAxisd(offset * localization::D2R, Eigen::Vector3d::UnitY()) *
	Eigen::AngleAxisd(0.0 * localization::D2R, Eigen::Vector3d::UnitX()));
	
	std::cout<<"*** Quaternions  *** \n";
	std::cout<<"attitude: "<<"w:"<<attitude.w()<<"x:"<<attitude.x()<<"  y:"<<attitude.y()<<"  z:"<<attitude.z()<<"\n";
	std::cout<<"q_weight: "<<"w:"<<q_weight.w()<<"x:"<<q_weight.x()<<"  y:"<<q_weight.y()<<"  z:"<<q_weight.z()<<"\n";
	
	v2 = (attitude * q_weight) * v1;
	
	std::cout<<"v1\n"<<v1<<"\n";
	std::cout<<"v2\n"<<v2<<"\n";
	

// 	Eigen::Matrix <double,localization::NUMAXIS,localization::NUMAXIS> angle2product;
// 	
// 	v1 << 0.0, 0.0, 1.0;
// 	
// 	angle2product << 0, -0.00, cos(offset),
// 		    0.00, 0, 0.00,
// 		    sin(-offset), 0.00, 0;
// 		    
// 	v2 = angle2product * v1;
// 	
// 	std::cout<<"v1\n"<<v1<<"\n";
// 	std::cout<<"v2\n"<<v2<<"\n";
//


}


BOOST_AUTO_TEST_CASE( PERMUTATIONS )
{
    	double n[] = {0,0,0,0};
	int k = 1;
	
	do{
	    for(int i = 0; i < k; i++)
	    {
		std::cout << n[i];
	    }
	    std::cout << '\n';
	} while (std::next_permutation(n, n + k));
	
// 	int myints[] = {0,0,0,0};
// 	std::vector<int> fifth (myints, myints + sizeof(myints) / sizeof(int) );
	std::vector<int> fifth (4, 0);
	
	
	std::cout<<"fifth size: "<<fifth.size()<<"\n";
	for (int l=0; l<16; ++l)
	{
	    
	    for( std::vector<int>::const_iterator i = fifth.begin(); i != fifth.end(); ++i)
		std::cout << *i << ' ';
	    
	    bool value = increment(fifth,k);
	    std::cout<<"\n bool: "<<value<<"\n";
	    
	}
	
	fifth = std::vector<int>(4, 0);
	
	std::cout<<"Incrementing in 5\n";
	for (int l=0; l<1; ++l)
	    increment(fifth,1);
	    
	for( std::vector<int>::const_iterator i = fifth.begin(); i != fifth.end(); ++i)
	    std::cout << *i << ' ';
	    
	std::cout<<"\n";

}
