#if defined KOKKOS_DEVICE_CUDA
#include <cuda_runtime.h>
#include "cusparse.h"
#endif

#include <MatVecMult.hpp>

template<class Scalar, class DeviceType>
struct CRSMatVec;

#if defined KOKKOS_DEVICE_CUDA
//Specialize MatVecMult functor of CUDA device
template<class Scalar>
struct CRSMatVec<Scalar , Kokkos::DeviceCuda>
{
	typedef Kokkos::MultiVectorView<Scalar, Kokkos::DeviceCuda> 		scalar_vector;
	typedef Kokkos::MultiVectorView<int, Kokkos::DeviceCuda> 			int_vector;	 
	
	void MatVec(Scalar , Scalar);
	cusparseStatus_t status;
	cusparseHandle_t handle;
	cusparseMatDescr_t descra;
	
	Scalar * value, *x , *b ;
	int * row ,* col;
	int m;
	
	//Constructor that initalizes cusparse library and 
	//gets pointers on the device to use 
	CRSMatVec(scalar_vector & arg_value , int_vector & arg_row , int_vector arg_col,
					 scalar_vector & arg_x , scalar_vector & arg_b)
	{
		status = cusparseCreate(&handle);
		if(status != CUSPARSE_STATUS_SUCCESS)
		{
			std::cout<<"ERROR - Library Initialization failed"<<std::endl;
		}

		status = cusparseCreateMatDescr(&descra);
		if(status != CUSPARSE_STATUS_SUCCESS)
		{
			std::cout<<"ERROR - Matrix descriptor failed"<<std::endl;
		}

		cusparseSetMatType(descra , CUSPARSE_MATRIX_TYPE_GENERAL);
		cusparseSetMatIndexBase(descra , CUSPARSE_INDEX_BASE_ZERO);

		value = (Scalar*)arg_value.ptr_on_device();
		row = (int*)arg_row.ptr_on_device();
		col = (int*)arg_col.ptr_on_device();
		x = (Scalar*)arg_x.ptr_on_device();
		b = (Scalar*)arg_b.ptr_on_device();
		m = arg_row.length() - 1;
	}
};

//Compute equivalent MatVec operation using 
//cuda sparse library function
template<>
void CRSMatVec<float, Kokkos::DeviceCuda >::MatVec(float alpha , float beta)
{

	cusparseScsrmv(handle , CUSPARSE_OPERATION_NON_TRANSPOSE , m, m , alpha , 
				descra , value , row , col , x , beta , b);
}

template<>
void CRSMatVec<double, Kokkos::DeviceCuda>::MatVec(double alpha , double beta)
{
	cusparseDcsrmv(handle , CUSPARSE_OPERATION_NON_TRANSPOSE , m, m , alpha , 
				descra , value , row , col , x , beta , b);
}

#endif



template<class Scalar, class DeviceType>
struct CRSMatVec
{
	typedef DeviceType										device_type;
	typedef Kokkos::MultiVectorView<Scalar, device_type>	scalar_vector;
	typedef Kokkos::MultiVectorView<int, device_type> 		int_vector;	 
	
	void MatVec(Scalar alpha , Scalar beta);
	int rows;

	scalar_vector value ;
	int_vector row , col;
	scalar_vector x , b;
	CRSMatVec(scalar_vector & arg_value , int_vector & arg_row , int_vector & arg_col,
					 scalar_vector & arg_x , scalar_vector & arg_b)
			: value(arg_value) , row(arg_row) , col(arg_col) , x(arg_x) , b(arg_b)			
	{
		value = arg_value;
		row = arg_row;
		col = arg_col;
		x = arg_x;
		b = arg_b;
		rows = arg_row.length() - 1;
	}
};


template<class Scalar , class DeviceType>
void CRSMatVec<Scalar , DeviceType>::MatVec(Scalar alpha , Scalar beta)
{
	Kokkos::parallel_for(rows , MatVecMult<Scalar , DeviceType>(value , row , col , x , b , beta) );
}

