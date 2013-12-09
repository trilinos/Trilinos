// @HEADER
// ***********************************************************************
//
//            Domi: Multidimensional Datastructures Package
//                 Copyright (2013) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
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
// Questions? Contact William F. Spotz (wfspotz@sandia.gov)
//
// ***********************************************************************
// @HEADER


#include "Domi_getValidParameters.hpp"
#include "Teuchos_Array.hpp"
#include "Domi_Utils.hpp"

namespace Domi

{

  Teuchos::RCP< const Teuchos::ParameterList > getValidParameters() 
  {
    static Teuchos::RCP< const Teuchos::ParameterList > result;

    if ( result.is_null() )
      {

	Teuchos::ParameterList* plist = new Teuchos::ParameterList;

	/* Parameter applies to MDComm */
	Teuchos::Array< int > axisCommSizes(1);
	axisCommSizes[0] = -1;
	plist->set("axisCommSizes", axisCommSizes, 
		   "An array of ints whose length is the number of dimensions "
		   "of the problem and whose entries specify the size of the "
		   "MDComm along each axis. A negative value tells Domi to "
		   "fill in a logical value based on the total number of "
		   "processors");

	/* Parameter applies to MDComm */
	Teuchos::Array< int > periodic(1);
	periodic[0] = 0;
	plist->set("periodic", periodic,
		   "An array of int flags specifying whether each axis is a "
		   "periodic axis. If the periodic array is shorter than the "
		   "length of axisCommSizes array, then the unspecified "
		   "entries are given a default value of zero (not "
		   "periodic).");

	/* Parameter applies to MDMap */
	Teuchos::Array< long long > dimensions(1);
	dimensions[0] = 0;
	plist->set("dimensions", dimensions,
		   "An array of ordinals specifying the global dimensions of "
		   "the MDMap. The length of this array should be the same as "
		   "the length of the axisCommSizes array. Note that this is "
		   "an array of long long and it will need to be copied to an "
		   "array of type GlobalOrd.");

	Teuchos::Array< int > pad;
	/* Parameter applies to MDMap */
	plist->set("boundaryPad", pad,
		  "An array of ints specifying the size of the boundary "
		   "padding along each axis. All unspecified entries are "
		   "assumed to be zero.");

	/* Parameter applies to MDMap */
	plist->set("commPad", pad,
		  "An array of ints specifying the size of the commmunication "
		   "padding along each axis. All unspecified entries are "
		   "assumed to be zero.");


	/* Parameter applies to MDMap */
	std::string layout = "Default";
	plist->set("layout", layout,
		   "A string indicating how the data is laid out in memory. "
		   "Default is currently set to Fortran order.");

	result.reset(plist);

      }
    

    return result;
 
  }

}
