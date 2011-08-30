
namespace Test {


template <class Scalar , class DeviceType , class CellTraits >
struct Jacobian; 

//Specialized for the hexahedron and case 2 multiple jacobian for a single set of reference points
template<class Scalar >
struct Jacobian<Scalar , KOKKOS_MACRO_DEVICE , shards::Hexahedron<8> >
{
	typedef KOKKOS_MACRO_DEVICE 		device_type;
	typedef device_type::size_type 		size_type;
	
	typedef typename Kokkos::MDArrayView<Scalar,device_type> array_type;
	typedef typename Kokkos::MDArrayView<Scalar,Kokkos::DeviceHost> host_array;
	
  private:
  
    array_type jacobian ;
    array_type cellcoords;
	int spaceDim ;
	int numCells ;
	int numPoints ;
	int basisCardinality;
	array_type basisGrads_device;
  public:
  
  	Jacobian( 	array_type &	arg_jacobian ,
  				const host_array &	arg_points_host ,
  				const array_type &	arg_cellcoords )
  	: jacobian(arg_jacobian) , cellcoords(arg_cellcoords) 
  	{ 
  		
  		spaceDim = shards::Hexahedron<8>::dimension;
  		numCells = cellcoords.dimension(0);
  		numPoints = arg_points_host.dimension(0);
		
		
		//Specialized for hexahedron<8>
  		Intrepid::Basis_HGRAD_HEX_C1_FEM<Scalar, host_array > HGRAD_Basis;
  		
  		basisCardinality = HGRAD_Basis.getCardinality();
  		
  		//Create local temporary host MDArray to get basisGrad on host
  		host_array basisGrads = Kokkos::create_mdarray<host_array>(basisCardinality, numPoints , spaceDim);
  		
  		//Data shared among all calls
  		basisGrads_device = Kokkos::create_mdarray<array_type>(basisCardinality, numPoints , spaceDim);
  		
        HGRAD_Basis.getValues(basisGrads, arg_points_host, Intrepid::OPERATOR_GRAD);
		
		//Copy basisGrads onto device       
        Kokkos::deep_copy(basisGrads_device , basisGrads);

  	}
  	
  	KOKKOS_MACRO_DEVICE_FUNCTION
  	void operator()(size_type ielem) const 
  	{
		for(int pointOrd = 0; pointOrd < numPoints; pointOrd++) {
			for(int row = 0; row < spaceDim; row++){
         		for(int col = 0; col < spaceDim; col++){
                    
                    // The entry is computed by contracting the basis index. Number of basis functions and vertices must be the same.
         	    	for(int bfOrd = 0; bfOrd < basisCardinality; bfOrd++){
                  		jacobian(ielem, pointOrd, row, col) += cellcoords(ielem, bfOrd, row)*basisGrads_device(bfOrd, pointOrd, col);
                  	} // bfOrd
               	} // col
       		} // row
		} // pointOrd
  	}

};




} // namespace Test
