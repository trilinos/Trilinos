namespace Test {

template< typename Scalar , class DeviceType>
struct simpleFill;

template<typename Scalar>
struct simpleFill<Scalar, KOKKOS_MACRO_DEVICE>{

	typedef KOKKOS_MACRO_DEVICE     							device_type ;
  	typedef typename Kokkos::MDArrayView<Scalar,device_type> 	array_type ;

	array_type data;

  	simpleFill(	array_type & arg_data):		data(arg_data)	{ ;} 


	KOKKOS_MACRO_DEVICE_FUNCTION
  	void operator()( int ielem )const {

		for(unsigned int i = 0; i < data.dimension(1); i++){

			data(ielem, i) = 1.0;

		}

	}

}; // struct

} // namespace

