template<class Scalar , class DeviceType >
struct MatVecMult;

template<class Scalar >
struct MatVecMult<Scalar , KOKKOS_MACRO_DEVICE >
{
	
	typedef KOKKOS_MACRO_DEVICE 									device_type;
	typedef device_type::size_type									size_type;
	typedef typename Kokkos::MultiVectorView<Scalar, device_type> 	scalar_vector;
	typedef typename Kokkos::MultiVectorView<int, device_type> 		int_vector;	

	scalar_vector A_value;
	int_vector A_row;
	int_vector A_col;

	scalar_vector x;
	scalar_vector b;
	Scalar beta;
  
  	MatVecMult(	scalar_vector & arg_A_value, int_vector & arg_A_row , int_vector & arg_A_col , 
  			scalar_vector & arg_x , scalar_vector & arg_b, Scalar arg_beta) 
  	: A_value(arg_A_value) , A_row(arg_A_row) , A_col(arg_A_col) , x(arg_x) , b(arg_b) , beta(arg_beta)
  	{  	}
  	
  	
	// b = beta*b + A * x
  	KOKKOS_MACRO_DEVICE_FUNCTION
  	void operator()(int irow ) const 
  	{

		Scalar sum = beta*b(irow);

		const size_type end = A_row(irow+1);

		for(size_type i = A_row(irow) ; i < end; i++)
		{
			sum += A_value(i) * x(A_col(i));
		}
		
		b(irow) = sum;

  	}
  	

}; //MatVec


