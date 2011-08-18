
namespace Test {

//Assume operatorIntegral, output is rank 3 , fields is rank 4 (contractFieldVector)
template <class Scalar , class DeviceType >
struct Integrate; 

template<class Scalar >
struct Integrate<Scalar , KOKKOS_MACRO_DEVICE >
{
	typedef KOKKOS_MACRO_DEVICE 		device_type;
	typedef device_type::size_type 		size_type;
	
	typedef typename Kokkos::MDArrayView<Scalar,device_type> array_type;
	typedef typename Kokkos::MDArrayView<Scalar,Kokkos::DeviceHost> host_array;
	
  private:
	
	array_type output;
	array_type left;
	array_type right;  

	int numLeft;
	int numRight;
	int numPoints;
	int dim;
	
  public:
  
	Integrate(	array_type 			& arg_output , 
				const array_type	& arg_left ,
				const array_type	& arg_right  ) : output(arg_output) , left(arg_left) , right(arg_right)
	{
		numLeft = left.dimension(1);
		numRight = right.dimension(1);
		numPoints = left.dimension(2);
		dim = left.dimension(3);
		if(output.rank() == 2) numLeft = 1;
	}
  	
  	//Assume compEngine is COMP_CPP, sumInto = false
  	KOKKOS_MACRO_DEVICE_FUNCTION
  	void operator()(size_type ielem) const 
  	{
		for (int lbf = 0; lbf < numLeft; lbf++) {
        	for (int rbf = 0; rbf < numRight; rbf++) {
            	Scalar tmpVal(0);
              		for (int qp = 0; qp < numPoints; qp++) {
                		for (int iVec = 0; iVec < dim; iVec++) {
                  			tmpVal += left(ielem, lbf, qp, iVec)*right(ielem, rbf, qp, iVec);
		                } //D-loop
        		      } // P-loop
           	   output(ielem, lbf, rbf) = tmpVal;
	        } // R-loop
		} // L-loop
  	}

};




} // namespace Test
