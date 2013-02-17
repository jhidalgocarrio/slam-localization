/**\file DataModel.cpp
 *
 * This class has the primitive methods for the data Model 
 *  of the measurement generation
 * 
 * 
 * @author Javier Hidalgo Carrio | DFKI RIC Bremen | javier.hidalgo_carrio@dfki.de
 * @date December 2012.
 * @version 1.0.
 */

#include "DataModel.hpp" 

using namespace localization;


DataModel::DataModel()
{
    data.resize(1,1);
    data.setZero();
    Cov.resize(1,1);
    Cov.setIdentity();
    
    
    Cov = ZERO_UNCERTAINTY * Cov;
}

DataModel::DataModel(const unsigned int dim)
{
    data.resize(dim,1);
    data.setZero();
    Cov.resize(dim,dim);
    Cov.setIdentity();    
    
    Cov = ZERO_UNCERTAINTY * Cov;
    
}

DataModel::DataModel(Eigen::Matrix< double, Eigen::Dynamic, 1  >& data, Eigen::Matrix< double, Eigen::Dynamic, Eigen::Dynamic > &Cov)
{
    this->data.resize(data.size(),1);
    this->data = data;
    
    this->Cov.resize(Cov.rows(), Cov.cols());
    this->Cov = Cov;

}


int DataModel::size()
{
    return (this->data.rows());

}


void DataModel::fusion(const DataModel& data2)
{
    DataModel &data1 (*this);
    
    if ((data1.Cov.size() != 0 && data2.Cov.size() != 0) && (data1.Cov.size() == data2.Cov.size()))
    {
	Eigen::Matrix <double, Eigen::Dynamic, Eigen::Dynamic> P;
	
	P.resize(data1.Cov.rows(), data1.Cov.cols());
	
	P = (data1.Cov.inverse() + data2.Cov.inverse()).inverse();
	data1.data = P*(data1.Cov.inverse() * data1.data + data2.Cov.inverse() * data2.data);
	
	data1.Cov = P;
    }
    
    return;
}

void DataModel::safeFusion(const localization::DataModel& data2)
{
    DataModel &data1 (*this);
    
    /** Variables for the safe fusion **/
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> I1, I2;
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> U1, sqrtD1, isqrtD1;
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> U2, D2;
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> D3;
    
    Eigen::Matrix<double, Eigen::Dynamic, 1> result;
    Eigen::Matrix<double, Eigen::Dynamic, 1> datatrans1, datatrans2;
    
    /** Resize vectors and matrices **/
    I1.resize (data1.Cov.rows(), data1.Cov.cols());
    U1.resize (data1.Cov.rows(), data1.Cov.cols());
    sqrtD1.resize (data1.Cov.rows(), data1.Cov.cols());
    isqrtD1.resize (data1.Cov.rows(), data1.Cov.cols());
    
    I2.resize (data2.Cov.rows(), data2.Cov.cols());
    D2.resize (data2.Cov.rows(), data2.Cov.cols());
    
    D3.resize (data2.Cov.rows(), data2.Cov.cols());
    
    result.resize(data1.data.size(), 1);
    datatrans1.resize(data1.data.size(), 1);
    datatrans2.resize(data2.data.size(), 1);
    
    D3.setZero();
    
    I1 = data1.Cov.inverse();
    I2 = data2.Cov.inverse();
	    
    /**
    * Single Value Decomposition
    */
    Eigen::JacobiSVD <Eigen::MatrixXd > svdOfI1(I1, Eigen::ComputeThinU);
    
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
    
    Eigen::Matrix<double, NUMAXIS, NUMAXIS> T;
	    
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


DataModel DataModel::operator+(const DataModel& data2) const
{
    const DataModel &data1 (*this);
    DataModel result;
    
    result.data = data1.data + data2.data;
    result.Cov = data1.Cov + data2.Cov;
    
    return result;
}

DataModel DataModel::operator-(const DataModel& data2) const
{
    const DataModel &data1 (*this);
    DataModel result;
    
    result.data = data1.data - data2.data;
    result.Cov = data1.Cov + data2.Cov;
    
    return result;
}

DataModel& DataModel::operator=(const DataModel& dmodel)
{
    data = dmodel.data;
    Cov = dmodel.Cov;
    
    return *this;
}


std::ostream& localization::operator<<(std::ostream &out, const  DataModel& data1)
{
    out <<"\n" << data1.data << "\n";
    out <<"\n" << data1.Cov << "\n";
    
    return out;
}





