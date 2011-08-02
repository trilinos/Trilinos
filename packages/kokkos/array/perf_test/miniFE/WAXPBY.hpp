template<class Scalar , class DeviceType >
struct WAXPBY;

template<class Scalar >
struct WAXPBY<Scalar , KOKKOS_MACRO_DEVICE >
{
	
	typedef KOKKOS_MACRO_DEVICE 							device_type;
	typedef device_type::size_type							size_type;
	typedef Kokkos::MultiVectorView<Scalar, device_type> 	scalar_vector;	
	typedef Kokkos::ValueView<Scalar , device_type>			value;


	value alpha;
	scalar_vector x;
	
	value beta;
	scalar_vector y;
	
	scalar_vector w;
	
  	WAXPBY(value arg_alpha ,  scalar_vector & arg_x , value arg_beta , scalar_vector & arg_y, 
  			 scalar_vector & arg_w )
  		:  alpha(arg_alpha) , x(arg_x), beta(arg_beta) , y(arg_y) , w(arg_w) { }
  	
  	
  	
  	KOKKOS_MACRO_DEVICE_FUNCTION
  	void operator()(int inode) const 
  	{
		w(inode) = (*alpha)*x(inode) + (*beta)*y(inode);
  	}

}; //WAXPBY

