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

//Assume operatorIntegral, output is rank 3 , fields is rank 4 (contractFieldVector)
template <class Scalar , class DeviceType >
struct Integrate; 

template<class Scalar >
struct Integrate<Scalar , KOKKOSARRAY_MACRO_DEVICE >
{
	typedef KOKKOSARRAY_MACRO_DEVICE 		device_type;
	typedef device_type::size_type 		size_type;
	
	typedef typename KokkosArray::MDArrayView<Scalar,device_type> array_type;
	typedef typename KokkosArray::MDArrayView<Scalar,KokkosArray::DeviceHost> host_array;
	
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
  	KOKKOSARRAY_MACRO_DEVICE_FUNCTION
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
