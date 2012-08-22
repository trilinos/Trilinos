/*
//@HEADER
// ************************************************************************
// 
//          Kokkos: Node API and Parallel Node Kernels
//              Copyright (2008) Sandia Corporation
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ************************************************************************
//@HEADER
*/

#include <impl/KokkosArray_Preprocessing_macros.hpp>

namespace Test {

template <typename Scalar , class DeviceType >
struct HexFill;

template< typename Scalar >
struct HexFill< Scalar , KOKKOSARRAY_MACRO_DEVICE >
{
	typedef KOKKOSARRAY_MACRO_DEVICE		 device_type ;
	typedef device_type::size_type	 size_type ;
	
	typedef typename KokkosArray::MDArrayView<Scalar,device_type> array_type;
	
	array_type coords;
	
	HexFill( const array_type & arg_coords) : coords(arg_coords) {}
	
	KOKKOSARRAY_MACRO_DEVICE_FUNCTION
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
void poisson_run< KOKKOSARRAY_MACRO_DEVICE>(int beg , int end)
{
	std::ofstream outfile;
	
	std::string label_poisson;
	label_poisson.append("\"Poisson< double , ");
	label_poisson.append(KOKKOSARRAY_MACRO_TO_STRING(KOKKOSARRAY_MACRO_DEVICE) );
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
	typedef KOKKOSARRAY_MACRO_DEVICE device_type;
	typedef device_type::size_type size_type;
	typedef KokkosArray::MDArrayView<Scalar , device_type > array_type;
	typedef KokkosArray::MDArrayView<Scalar , KokkosArray::DeviceHost > host_array;
	
	
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

	array_type cellcoords = KokkosArray::create_mdarray<array_type>(numCells , numNodesPerCell , spaceDim );


	
	KokkosArray::parallel_for(numCells, HexFill<Scalar , device_type>(cellcoords) , temp);
	
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
    

	cubPointsHost = KokkosArray::create_mdarray<host_array>(numCubPoints,cubDim);
	cubWeightsHost = KokkosArray::create_mdarray<host_array>(numCubPoints);
	
	//Device allocation
	cubWeightsDevice = KokkosArray::create_mdarray<array_type>(numCubPoints);
	
	cellCubature.getCubature(cubPointsHost , cubWeightsHost);
	}
	//Deep Copy
	KokkosArray::deep_copy(cubWeightsDevice , cubWeightsHost);
 
 
/*********************************************************/
/*                 Define cell basis                     */
/*********************************************************/

	Intrepid::Basis_HGRAD_HEX_C1_FEM<Scalar, host_array > hexHGradBasis;
  		
  	int numFields = hexHGradBasis.getCardinality();
  	
	host_array basisGrads = KokkosArray::create_mdarray<host_array>(numFields , numCubPoints , spaceDim);
	host_array basisValues = KokkosArray::create_mdarray<host_array>(numFields, numCubPoints);
	
	// Evaluate basis values and gradients at cubature points
    hexHGradBasis.getValues(basisValues, cubPointsHost, Intrepid::OPERATOR_VALUE);
    hexHGradBasis.getValues(basisGrads, cubPointsHost, Intrepid::OPERATOR_GRAD);
    
    array_type basisGrads_device = KokkosArray::create_mdarray<array_type>(numFields , numCubPoints , spaceDim);
	array_type basisValues_device = KokkosArray::create_mdarray<array_type>(numFields, numCubPoints);
	
	//Copy to device
	KokkosArray::deep_copy(basisGrads_device , basisGrads);
	KokkosArray::deep_copy(basisValues_device , basisValues);
	

/*********************************************************/
/*            Calculate cell Jacobians                   */
/*********************************************************/
    
    array_type worksetJacobian = KokkosArray::create_mdarray<array_type>(numCells , numCubPoints , 
																	spaceDim , spaceDim );
	

	array_type worksetJacobianInv = KokkosArray::create_mdarray<array_type>(numCells, numCubPoints, spaceDim, spaceDim);

	array_type worksetJacobianDet = KokkosArray::create_mdarray<array_type>(numCells, numCubPoints);
	
	//  Set Jacobian, inverse and determinant
	KokkosArray::parallel_for(numCells , Test::Jacobian<Scalar , device_type , CellTraits >(worksetJacobian , cubPointsHost , cellcoords) , temp );

	seconds = seconds + temp;

	
	KokkosArray::parallel_for(numCells, Test::Invert<Scalar,device_type, 3>(worksetJacobian, worksetJacobianInv, 8), temp);
	
	seconds = seconds + temp;


	KokkosArray::parallel_for(numCells , Test::Determinant<Scalar , device_type>(worksetJacobianDet , worksetJacobian), temp);
	
	seconds = seconds + temp;


/*********************************************************/
/*             Compute Stiffness Matrix                  */
/*********************************************************/


		// Containers for stiffness matrix
	array_type worksetStiffMatrix = KokkosArray::create_mdarray<array_type>(numCells, numFields , numFields);
	array_type worksetCubWeights = KokkosArray::create_mdarray<array_type>(numCells, numCubPoints);
	array_type worksetBasisGrads = KokkosArray::create_mdarray<array_type>(numCells, numFields , numCubPoints , spaceDim);
	array_type worksetBasisGradsWeighted = KokkosArray::create_mdarray<array_type>(numCells, numCubPoints , spaceDim);
		
	KokkosArray::parallel_for(numCells , Test::Transform<Scalar , device_type>(worksetBasisGrads , worksetJacobianInv , basisGrads_device) , temp);	
	
		seconds = seconds + temp;
	
	KokkosArray::parallel_for(numCells , Test::computeCellMeasure<Scalar , device_type>(worksetCubWeights , worksetJacobianDet , cubWeightsDevice) , temp);	
	
		seconds = seconds + temp;
	
	
	KokkosArray::parallel_for(numCells , Test::Multiply<Scalar , device_type>(worksetBasisGradsWeighted , worksetCubWeights , worksetBasisGrads) , temp);
	
		seconds = seconds + temp;
		
		
	KokkosArray::parallel_for(numCells , Test::Integrate<Scalar , device_type>(worksetStiffMatrix , worksetCubWeights , worksetBasisGrads) , temp);
	
		seconds = seconds + temp;
	

/**********************************************************************************/
/*                                   Compute RHS                                  */
/**********************************************************************************/	

	array_type worksetRHS = KokkosArray::create_mdarray<array_type>(numCells , numFields);
	array_type worksetBasisValues  = KokkosArray::create_mdarray<array_type>(numCells , numFields, numCubPoints);
	array_type worksetBasisValuesWeighted  = KokkosArray::create_mdarray<array_type>(numCells , numFields, numCubPoints);

	array_type worksetSourceTerm  = KokkosArray::create_mdarray<array_type>(numCells , numCubPoints);
	
	
	KokkosArray::parallel_for(numCells , Test::simpleFill<Scalar , device_type>(worksetSourceTerm) , temp);
	
		seconds = seconds + temp;
	
	KokkosArray::parallel_for(numCells , Test::TransformValue<Scalar , device_type>(worksetBasisValues ,basisValues_device ) , temp);
	
		seconds = seconds + temp;
	
	KokkosArray::parallel_for(numCells , Test::Multiply<Scalar , device_type>(worksetBasisValuesWeighted,worksetCubWeights,worksetBasisValues) , temp);
	
		seconds = seconds + temp;
	
	KokkosArray::parallel_for(numCells , Test::Integrate<Scalar , device_type>(worksetRHS ,worksetSourceTerm ,  worksetBasisValuesWeighted) , temp);
	
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
