template<class Scalar , class DeviceType >
struct Dot;

template<class Scalar >
struct Dot<Scalar , KOKKOS_MACRO_DEVICE >
{
	
	typedef KOKKOS_MACRO_DEVICE 								device_type;
	typedef device_type::size_type								size_type;
	typedef Kokkos::MultiVectorView<Scalar, device_type> 		scalar_vector;	
	typedef Scalar 												value_type;


	scalar_vector x;
	scalar_vector y;
	
  
  	Dot(scalar_vector & arg_x , scalar_vector & arg_y) :  x(arg_x) , y(arg_y) { }
  	
  	
  	
  	KOKKOS_MACRO_DEVICE_FUNCTION
  	void operator()(int iwork , value_type & update) const 
  	{
		update += x(iwork) * y(iwork);
  	}
  	
  	KOKKOS_MACRO_DEVICE_FUNCTION
  	static void join(volatile value_type & update , const volatile value_type & source)
  	{
  		update += source;  	
  	}
  	
  	KOKKOS_MACRO_DEVICE_FUNCTION
  	static void init( value_type & update )
  	{ 
  		update = 0 ; 
  	}

}; //Dot


