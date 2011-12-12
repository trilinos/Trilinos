namespace Test {

template< typename Scalar , class DeviceType>
struct computeCellMeasure;

template<typename Scalar>
struct computeCellMeasure<Scalar, KOKKOS_MACRO_DEVICE>{

	typedef KOKKOS_MACRO_DEVICE     							device_type ;
  	typedef typename Kokkos::MDArrayView<Scalar,device_type> 	array_type ;

	array_type inDet;
	array_type inWeights;
	array_type outVals;

  	computeCellMeasure(	array_type & arg_outVals,
						array_type & arg_inDet, 
						array_type & arg_inWeights):
		inDet		(arg_inDet), 
		inWeights	(arg_inWeights), 
		outVals		(arg_outVals)
	{}


	KOKKOS_MACRO_DEVICE_FUNCTION
  	void operator()( int ielem )const {

		for(unsigned int point = 0; point < outVals.dimension(1); point++){

			outVals(ielem, point) = inDet(ielem, point) * inWeights(point);
			
		}// for, cubature

		if(inDet(ielem, 0) < 0.0){			

			for(unsigned int point = 0; point < outVals.dimension(1); point++){

				outVals(ielem, point) *= -1;
		
			}// for, point

		}// if

	}

}; // struct

} // namespace

