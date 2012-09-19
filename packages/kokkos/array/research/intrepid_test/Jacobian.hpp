/*
//@HEADER
// ************************************************************************
// 
//   KokkosArray: Manycore Performance-Portable Multidimensional Arrays
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


namespace Test {


template <class Scalar , class DeviceType , class CellTraits >
struct Jacobian; 

//Specialized for the hexahedron and case 2 multiple jacobian for a single set of reference points
template<class Scalar >
struct Jacobian<Scalar , KOKKOSARRAY_MACRO_DEVICE , shards::Hexahedron<8> >
{
	typedef KOKKOSARRAY_MACRO_DEVICE 		device_type;
	typedef device_type::size_type 		size_type;
	
	typedef typename KokkosArray::MDArrayView<Scalar,device_type> array_type;
	typedef typename KokkosArray::MDArrayView<Scalar,KokkosArray::DeviceHost> host_array;
	
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
  		host_array basisGrads = KokkosArray::create_mdarray<host_array>(basisCardinality, numPoints , spaceDim);
  		
  		//Data shared among all calls
  		basisGrads_device = KokkosArray::create_mdarray<array_type>(basisCardinality, numPoints , spaceDim);
  		
        HGRAD_Basis.getValues(basisGrads, arg_points_host, Intrepid::OPERATOR_GRAD);
		
		//Copy basisGrads onto device       
        KokkosArray::deep_copy(basisGrads_device , basisGrads);

  	}
  	
  	KOKKOSARRAY_MACRO_DEVICE_FUNCTION
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
