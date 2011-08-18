/*
* Finalize functor for dot product. 
* Returns both dot product result 
* as well quotient of divison operation.	
*/
template<class Scalar , class DeviceType >
struct Divide;

template<class Scalar >
struct Divide<Scalar , KOKKOS_MACRO_DEVICE>
{
	typedef KOKKOS_MACRO_DEVICE 									device_type;
	typedef device_type::size_type									size_type;
	typedef Kokkos::ValueView<Scalar , device_type>					value;



	value den;
	value val;
	value num;

  
  	Divide(value & arg_den, value & arg_val, value & arg_num) : 
  		 den(arg_den) , val(arg_val) , num(arg_num)  { }
  	
  	
  	
  	KOKKOS_MACRO_DEVICE_FUNCTION
  	void operator()(Scalar & result ) const 
  	{
  		*num = result;
		*val = (result)/(*den);
  	}
  		
}; //Divide
