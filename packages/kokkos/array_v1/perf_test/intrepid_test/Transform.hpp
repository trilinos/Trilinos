
namespace Test {


template <class Scalar , class DeviceType >
struct Transform; 

template<class Scalar >
struct Transform<Scalar , KOKKOS_MACRO_DEVICE >
{
	typedef KOKKOS_MACRO_DEVICE 		device_type;
	typedef device_type::size_type 		size_type;
	
	typedef typename Kokkos::MDArrayView<Scalar,device_type> array_type;
	typedef typename Kokkos::MDArrayView<Scalar,Kokkos::DeviceHost> host_array;
	
  private:
	
	array_type output;
	array_type input;
	array_type fields;  
	array_type basisGrads;
	int data_rank;
	int numDataPts;
	int in_rank;
	int numCells;
	int numFields;
	int numPoints;
	int dim;
	
  public:
  
	Transform(	array_type 			& arg_output , 
				const array_type	& arg_input ,
				const array_type	& arg_fields  ) : output(arg_output) , input(arg_input) , fields(arg_fields)
	{
		data_rank = input.rank();
		numDataPts = input.dimension(1);
		in_rank = fields.rank();
		numCells = output.dimension(0);
		numFields = output.dimension(1);
		numPoints = output.dimension(2);
		dim = output.dimension(3);		
	}
  	
  	//Assume fields is rank 3, input is rank 4 and transpose is on
  	KOKKOS_MACRO_DEVICE_FUNCTION
  	void operator()(size_type ielem) const 
  	{
		for(int field = 0; field < numFields; field++){
          for(int point = 0; point < numPoints; point++){
            for(int row = 0; row < dim; row++){
              output(ielem, field, point, row) = 0.0;
              for(int col = 0; col < dim; col++){
                output(ielem, field, point, row) += \
                  input(ielem, point, col, row)*fields(field, point, col);
              }// col
            } //row
          }// point
        }// field
  	}

};




} // namespace Test
