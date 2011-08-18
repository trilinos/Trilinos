#include <impl/Kokkos_Preprocessing_macros.hpp>

namespace Test {

template <typename Scalar , class DeviceType >
struct HexFill;

template< typename Scalar >
struct HexFill< Scalar , KOKKOS_MACRO_DEVICE >
{
	typedef KOKKOS_MACRO_DEVICE		 device_type ;
	typedef device_type::size_type	 size_type ;
	
	typedef typename Kokkos::MDArrayView<Scalar,device_type> array_type;
	
	array_type coords;
	
	HexFill( const array_type & arg_coords) : coords(arg_coords) {}
	
	KOKKOS_MACRO_DEVICE_FUNCTION
	void operator()( size_type ielem) const 
	{
		coords(ielem,0,0) = 0.;
		coords(ielem,0,1) = 0.;
		coords(ielem,0,2) = 0.;
	
		coords(ielem,1,0) = 1.2;
		coords(ielem,1,1) = -0.1;
		coords(ielem,1,2) = 0.2;
		
		coords(ielem,2,0) = 1.1;
		coords(ielem,2,1) = 1.1;
		coords(ielem,2,2) = -0.1;
		
		coords(ielem,3,0) = 0.1;
		coords(ielem,3,1) = 0.9;
		coords(ielem,3,2) = 0.1;
		
		coords(ielem,4,0) = 0.;
		coords(ielem,4,1) = -0.1;
		coords(ielem,4,2) = 1.1;
		
		coords(ielem,5,0) = 1.1;
		coords(ielem,5,1) = 0.1;
		coords(ielem,5,2) = 1.2;
		
		coords(ielem,6,0) = 0.9;
		coords(ielem,6,1) = 0.9;
		coords(ielem,6,2) = 0.9;
		
		coords(ielem,7,0) = -0.1;
		coords(ielem,7,1) = 1.;
		coords(ielem,7,2) = 1.;
	}
};


template< class DeviceType>
void poisson_run(int beg , int end);

template<>
void poisson_run< KOKKOS_MACRO_DEVICE>(int beg , int end)
{
	std::ofstream outfile;
	
	std::string label_poisson;
	label_poisson.append("\"Poisson< double , ");
	label_poisson.append(KOKKOS_MACRO_TO_STRING(KOKKOS_MACRO_DEVICE) );
	label_poisson.append(" >\"");

	double seconds = 0.0;
	double temp = 0.0;
	double min_seconds = 0.0;
	double max_seconds = 0.0;
	double avg_seconds = 0.0;

	for(int i = beg ; i < end ; ++i) {
	int numCells = (1<<i)+7;

	for(int j = 0 ; j < 5 ; j++) {
	seconds = 0.0 ; temp = 0.0;
	typedef double Scalar;
	typedef KOKKOS_MACRO_DEVICE device_type;
	typedef device_type::size_type size_type;
	typedef Kokkos::MDArrayView<Scalar , device_type > array_type;
	typedef Kokkos::MDArrayView<Scalar , Kokkos::DeviceHost > host_array;
	
	
/*********************************************************/
/*               Define cell topology                    */
/*********************************************************/

	typedef shards::Hexahedron<8> CellTraits;

    // Get dimensions
    int numNodesPerCell = CellTraits::node_count;
    int spaceDim = CellTraits::dimension;
	
/*********************************************************/
/*    Create hexahedral cell and copy into batch         */
/*********************************************************/

	array_type cellcoords = Kokkos::create_mdarray<array_type>(numCells , numNodesPerCell , spaceDim );


	
	Kokkos::parallel_for(numCells, HexFill<Scalar , device_type>(cellcoords) , temp);
	
/*********************************************************/
/*               Define cell cubature                    */
/*********************************************************/
	host_array cubPointsHost, cubWeightsHost;
	array_type cubPointsDevice, cubWeightsDevice;
	int cubDim , numCubPoints;
	{
	// Define cubature 
    int cubDegree = 2;
    std::vector< Teuchos::RCP< Intrepid::Cubature<Scalar , host_array> > > lineCubs(3);
    lineCubs[0]  = Teuchos::rcp(new Intrepid::CubatureDirectLineGauss<Scalar , host_array>(cubDegree));
    lineCubs[1]  = Teuchos::rcp(new Intrepid::CubatureDirectLineGauss<Scalar , host_array>(cubDegree));
    lineCubs[2]  = Teuchos::rcp(new Intrepid::CubatureDirectLineGauss<Scalar , host_array>(cubDegree));
    Intrepid::CubatureTensor<Scalar , host_array >cellCubature(lineCubs);
	
    cubDim       = cellCubature.getDimension(); 
    numCubPoints = cellCubature.getNumPoints();
    

	cubPointsHost = Kokkos::create_mdarray<host_array>(numCubPoints,cubDim);
	cubWeightsHost = Kokkos::create_mdarray<host_array>(numCubPoints);
	
	//Device allocation
	cubWeightsDevice = Kokkos::create_mdarray<array_type>(numCubPoints);
	
	cellCubature.getCubature(cubPointsHost , cubWeightsHost);
	}
	//Deep Copy
	Kokkos::deep_copy(cubWeightsDevice , cubWeightsHost);
 
 
/*********************************************************/
/*                 Define cell basis                     */
/*********************************************************/

	Intrepid::Basis_HGRAD_HEX_C1_FEM<Scalar, host_array > hexHGradBasis;
  		
  	int numFields = hexHGradBasis.getCardinality();
  	
	host_array basisGrads = Kokkos::create_mdarray<host_array>(numFields , numCubPoints , spaceDim);
	host_array basisValues = Kokkos::create_mdarray<host_array>(numFields, numCubPoints);
	
	// Evaluate basis values and gradients at cubature points
    hexHGradBasis.getValues(basisValues, cubPointsHost, Intrepid::OPERATOR_VALUE);
    hexHGradBasis.getValues(basisGrads, cubPointsHost, Intrepid::OPERATOR_GRAD);
    
    array_type basisGrads_device = Kokkos::create_mdarray<array_type>(numFields , numCubPoints , spaceDim);
	array_type basisValues_device = Kokkos::create_mdarray<array_type>(numFields, numCubPoints);
	
	//Copy to device
	Kokkos::deep_copy(basisGrads_device , basisGrads);
	Kokkos::deep_copy(basisValues_device , basisValues);
	

/*********************************************************/
/*            Calculate cell Jacobians                   */
/*********************************************************/
    
    array_type worksetJacobian = Kokkos::create_mdarray<array_type>(numCells , numCubPoints , 
																	spaceDim , spaceDim );
	

	array_type worksetJacobianInv = Kokkos::create_mdarray<array_type>(numCells, numCubPoints, spaceDim, spaceDim);

	array_type worksetJacobianDet = Kokkos::create_mdarray<array_type>(numCells, numCubPoints);
	
	//  Set Jacobian, inverse and determinant
	Kokkos::parallel_for(numCells , Test::Jacobian<Scalar , device_type , CellTraits >(worksetJacobian , cubPointsHost , cellcoords) , temp );

	seconds = seconds + temp;

	
	Kokkos::parallel_for(numCells, Test::Invert<Scalar,device_type, 3>(worksetJacobian, worksetJacobianInv, 8), temp);
	
	seconds = seconds + temp;


	Kokkos::parallel_for(numCells , Test::Determinant<Scalar , device_type>(worksetJacobianDet , worksetJacobian), temp);
	
	seconds = seconds + temp;


/*********************************************************/
/*             Compute Stiffness Matrix                  */
/*********************************************************/


		// Containers for stiffness matrix
	array_type worksetStiffMatrix = Kokkos::create_mdarray<array_type>(numCells, numFields , numFields);
	array_type worksetCubWeights = Kokkos::create_mdarray<array_type>(numCells, numCubPoints);
	array_type worksetBasisGrads = Kokkos::create_mdarray<array_type>(numCells, numFields , numCubPoints , spaceDim);
	array_type worksetBasisGradsWeighted = Kokkos::create_mdarray<array_type>(numCells, numCubPoints , spaceDim);
		
	Kokkos::parallel_for(numCells , Test::Transform<Scalar , device_type>(worksetBasisGrads , worksetJacobianInv , basisGrads_device) , temp);	
	
		seconds = seconds + temp;
	
	Kokkos::parallel_for(numCells , Test::computeCellMeasure<Scalar , device_type>(worksetCubWeights , worksetJacobianDet , cubWeightsDevice) , temp);	
	
		seconds = seconds + temp;
	
	
	Kokkos::parallel_for(numCells , Test::Multiply<Scalar , device_type>(worksetBasisGradsWeighted , worksetCubWeights , worksetBasisGrads) , temp);
	
		seconds = seconds + temp;
		
		
	Kokkos::parallel_for(numCells , Test::Integrate<Scalar , device_type>(worksetStiffMatrix , worksetCubWeights , worksetBasisGrads) , temp);
	
		seconds = seconds + temp;
	

/**********************************************************************************/
/*                                   Compute RHS                                  */
/**********************************************************************************/	

	array_type worksetRHS = Kokkos::create_mdarray<array_type>(numCells , numFields);
	array_type worksetBasisValues  = Kokkos::create_mdarray<array_type>(numCells , numFields, numCubPoints);
	array_type worksetBasisValuesWeighted  = Kokkos::create_mdarray<array_type>(numCells , numFields, numCubPoints);

	array_type worksetSourceTerm  = Kokkos::create_mdarray<array_type>(numCells , numCubPoints);
	
	
	Kokkos::parallel_for(numCells , Test::simpleFill<Scalar , device_type>(worksetSourceTerm) , temp);
	
		seconds = seconds + temp;
	
	Kokkos::parallel_for(numCells , Test::TransformValue<Scalar , device_type>(worksetBasisValues ,basisValues_device ) , temp);
	
		seconds = seconds + temp;
	
	Kokkos::parallel_for(numCells , Test::Multiply<Scalar , device_type>(worksetBasisValuesWeighted,worksetCubWeights,worksetBasisValues) , temp);
	
		seconds = seconds + temp;
	
	Kokkos::parallel_for(numCells , Test::Integrate<Scalar , device_type>(worksetRHS ,worksetSourceTerm ,  worksetBasisValuesWeighted) , temp);
	
		seconds = seconds + temp;
	
		
	if( 0 == j) {
		min_seconds = seconds;
		max_seconds = seconds;
	}
	else {
		if( seconds < min_seconds ) min_seconds = seconds;
		if( seconds > max_seconds ) max_seconds = seconds;
	}

}
	avg_seconds +=seconds;

	avg_seconds /= 5.0;
	std::cout 	<< label_poisson 
				<< " , " <<numCells
				<< " , " <<min_seconds
				<< " , " <<(min_seconds/numCells)
				<< std::endl;

}
}

} //Test
