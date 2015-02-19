/*
//@HEADER
// ************************************************************************
// 
//   Kokkos: Manycore Performance-Portable Multidimensional Arrays
//              Copyright (2012) Sandia Corporation
// 
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
// 
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact  H. Carter Edwards (hcedwar@sandia.gov) 
// 
// ************************************************************************
//@HEADER
*/

#include <impl/Kokkos_Preprocessing_macros.hpp>

namespace Test {

template <typename Scalar , class DeviceType >
struct HexFill;

template< typename Scalar >
struct HexFill< Scalar , KOKKOS_MACRO_DEVICE >
{
	typedef KOKKOS_MACRO_DEVICE		 execution_space ;
	typedef execution_space::size_type	 size_type ;
	
	typedef typename Kokkos::MDArrayView<Scalar,execution_space> array_type;
	
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
	typedef KOKKOS_MACRO_DEVICE execution_space;
	typedef execution_space::size_type size_type;
	typedef Kokkos::MDArrayView<Scalar , execution_space > array_type;
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


	
	Kokkos::parallel_for(numCells, HexFill<Scalar , execution_space>(cellcoords) , temp);
	
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
	Kokkos::parallel_for(numCells , Test::Jacobian<Scalar , execution_space , CellTraits >(worksetJacobian , cubPointsHost , cellcoords) , temp );

	seconds = seconds + temp;

	
	Kokkos::parallel_for(numCells, Test::Invert<Scalar,execution_space, 3>(worksetJacobian, worksetJacobianInv, 8), temp);
	
	seconds = seconds + temp;


	Kokkos::parallel_for(numCells , Test::Determinant<Scalar , execution_space>(worksetJacobianDet , worksetJacobian), temp);
	
	seconds = seconds + temp;


/*********************************************************/
/*             Compute Stiffness Matrix                  */
/*********************************************************/


		// Containers for stiffness matrix
	array_type worksetStiffMatrix = Kokkos::create_mdarray<array_type>(numCells, numFields , numFields);
	array_type worksetCubWeights = Kokkos::create_mdarray<array_type>(numCells, numCubPoints);
	array_type worksetBasisGrads = Kokkos::create_mdarray<array_type>(numCells, numFields , numCubPoints , spaceDim);
	array_type worksetBasisGradsWeighted = Kokkos::create_mdarray<array_type>(numCells, numCubPoints , spaceDim);
		
	Kokkos::parallel_for(numCells , Test::Transform<Scalar , execution_space>(worksetBasisGrads , worksetJacobianInv , basisGrads_device) , temp);	
	
		seconds = seconds + temp;
	
	Kokkos::parallel_for(numCells , Test::computeCellMeasure<Scalar , execution_space>(worksetCubWeights , worksetJacobianDet , cubWeightsDevice) , temp);	
	
		seconds = seconds + temp;
	
	
	Kokkos::parallel_for(numCells , Test::Multiply<Scalar , execution_space>(worksetBasisGradsWeighted , worksetCubWeights , worksetBasisGrads) , temp);
	
		seconds = seconds + temp;
		
		
	Kokkos::parallel_for(numCells , Test::Integrate<Scalar , execution_space>(worksetStiffMatrix , worksetCubWeights , worksetBasisGrads) , temp);
	
		seconds = seconds + temp;
	

/**********************************************************************************/
/*                                   Compute RHS                                  */
/**********************************************************************************/	

	array_type worksetRHS = Kokkos::create_mdarray<array_type>(numCells , numFields);
	array_type worksetBasisValues  = Kokkos::create_mdarray<array_type>(numCells , numFields, numCubPoints);
	array_type worksetBasisValuesWeighted  = Kokkos::create_mdarray<array_type>(numCells , numFields, numCubPoints);

	array_type worksetSourceTerm  = Kokkos::create_mdarray<array_type>(numCells , numCubPoints);
	
	
	Kokkos::parallel_for(numCells , Test::simpleFill<Scalar , execution_space>(worksetSourceTerm) , temp);
	
		seconds = seconds + temp;
	
	Kokkos::parallel_for(numCells , Test::TransformValue<Scalar , execution_space>(worksetBasisValues ,basisValues_device ) , temp);
	
		seconds = seconds + temp;
	
	Kokkos::parallel_for(numCells , Test::Multiply<Scalar , execution_space>(worksetBasisValuesWeighted,worksetCubWeights,worksetBasisValues) , temp);
	
		seconds = seconds + temp;
	
	Kokkos::parallel_for(numCells , Test::Integrate<Scalar , execution_space>(worksetRHS ,worksetSourceTerm ,  worksetBasisValuesWeighted) , temp);
	
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
