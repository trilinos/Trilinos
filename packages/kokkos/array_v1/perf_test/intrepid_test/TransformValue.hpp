
namespace Test {


template <class Scalar , class DeviceType >
struct TransformValue; 

template<class Scalar >
struct TransformValue<Scalar , KOKKOS_MACRO_DEVICE >
{
	typedef KOKKOS_MACRO_DEVICE 		device_type;
	typedef device_type::size_type 		size_type;
	
	typedef typename Kokkos::MDArrayView<Scalar,device_type> array_type;
	typedef typename Kokkos::MDArrayView<Scalar,Kokkos::DeviceHost> host_array;
	
  private:
	
	array_type output;
	array_type input;

	int numFields;
	int numPoints;
	
  public:
  
	TransformValue(	array_type 			& arg_output , 
				const array_type	& arg_input ) : output(arg_output) , input(arg_input) 
	{
		numFields = output.dimension(1);
		numPoints = output.dimension(2);
	}
  	
  	//Assume fields is rank 3, input is rank 2
  	KOKKOS_MACRO_DEVICE_FUNCTION
  	void operator()(size_type ielem) const 
  	{
		for(int bf = 0; bf < numFields; bf++) {
          for(int pt = 0; pt < numPoints; pt++) {
            output(ielem, bf, pt) = input(bf, pt);
          } // P-loop
        } // F-loop
  	}

};




} // namespace Test
