
namespace Test {


template <class Scalar , class DeviceType >
struct Multiply; 

template<class Scalar >
struct Multiply<Scalar , KOKKOS_MACRO_DEVICE >
{
	typedef KOKKOS_MACRO_DEVICE 		device_type;
	typedef device_type::size_type 		size_type;
	
	typedef typename Kokkos::MDArrayView<Scalar,device_type> array_type;
	typedef typename Kokkos::MDArrayView<Scalar,Kokkos::DeviceHost> host_array;
	
  private:
	
	array_type output;
	array_type input;
	array_type fields;  
	int data_rank;
	int numDataPts;
	int in_rank;
	int out_rank;
	int numFields;
	int numPoints;
	int dim;
	
  public:
  
	Multiply(	array_type 			& arg_output , 
				const array_type	& arg_input ,
				const array_type	& arg_fields  ) : output(arg_output) , input(arg_input) , fields(arg_fields)
	{
		data_rank = input.rank();
		numDataPts = input.dimension(1);
		in_rank = fields.rank();
		out_rank = output.rank(); 
		numFields = output.dimension(1);
		numPoints = output.dimension(2);
		dim = output.dimension(3);		
	}
  	
  	//Assume no reciprocal
  	KOKKOS_MACRO_DEVICE_FUNCTION
  	void operator()(size_type ielem) const 
  	{
  		switch(in_rank) {
  		
  			case 4: {
				for(int bf = 0; bf < numFields; bf++) {
					for(int pt = 0; pt < numPoints; pt++) {
				  		for( int iVec = 0; iVec < dim; iVec++) {
							output(ielem, bf, pt, iVec) = fields(ielem, bf, pt, iVec)*input(ielem, pt);
				  		} // D1-loop
					} // P-loop
			  	} // F-loop
			}
			
			case 3: {
			 for(int bf = 0; bf < numFields; bf++) {
                for(int pt = 0; pt < numPoints; pt++) {
                  output(ielem, bf, pt) = fields(ielem, bf, pt)*input(ielem, pt);
                } // P-loop
              } // F-loop
            }
    	}
  	}

};




} // namespace Test
