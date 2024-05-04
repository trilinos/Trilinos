// @HEADER
// ************************************************************************
//
//                           Intrepid Package
//                 Copyright (2007) Sandia Corporation
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
// Questions? Contact Pavel Bochev  (pbboche@sandia.gov)
//                    Denis Ridzal  (dridzal@sandia.gov), or
//                    Kara Peterson (kjpeter@sandia.gov)
//
// ************************************************************************
// @HEADER

/** \file   Intrepid_ProductTopologyDef.hpp
    \brief  Implementation file for products of lines
    \author Created by R. Kirby
*/

namespace Intrepid
{
  void ProductTopology::lineProduct2d( const int dim0 ,
				       const int entity0 ,
				       const int dim1 ,
				       const int entity1 ,
				       int &resultdim ,
				       int &resultentity )
  {
    // two vertices
    if (dim0 == 0 && dim1 == 0) {
      resultdim = 0;
      if (entity0 == 0 && entity1 == 0) {
	resultentity = 0;
      }
      else if (entity0 == 0 && entity1 == 1) {
	resultentity = 3;
      }
      else if (entity0 == 1 && entity1 == 0) {
	resultentity = 1;
      }
      else if (entity0 == 1 && entity1 == 1) {
	resultentity = 2;
      } 
      else {
	TEUCHOS_TEST_FOR_EXCEPTION( true , std::invalid_argument ,
			    "Intrepid::ProductTopology::lineProduct2D: illegal inputs" );
      }
    }
    else if (dim0 == 0 && dim1 == 1) {
      resultdim = 1;
      if (entity0 == 0 && entity1 == 0) {
	resultentity = 3;
      }
      else if (entity0 == 1 && entity1 == 0) {
	resultentity = 1;
      }
      else {
	TEUCHOS_TEST_FOR_EXCEPTION( true , std::invalid_argument ,
			    "Intrepid::ProductTopology::lineProduct2D: illegal inputs" );
      }
    }
    else if (dim0 == 1 && dim1 == 0) {
      resultdim = 1;
      if (entity0 == 0 && entity1 == 0) {
	resultentity = 0;
      }
      else if (entity0 == 0 && entity1 == 1) {
	resultentity = 2;
      }
      else {
	TEUCHOS_TEST_FOR_EXCEPTION( true , std::invalid_argument ,
			    "Intrepid::ProductTopology::lineProduct2D: illegal inputs" );
      }
    }
    else if (dim0 == 1 && dim1 == 1) {
      resultdim = 2;
      if (entity0 == 0 && entity1 == 0) {
	resultentity = 0;
      }
      else {
	TEUCHOS_TEST_FOR_EXCEPTION( true , std::invalid_argument ,
			    "Intrepid::ProductTopology::lineProduct2D: illegal inputs" );
      }
    }
    else {
      TEUCHOS_TEST_FOR_EXCEPTION( true , std::invalid_argument ,
			  "Intrepid::ProductTopology::lineProduct2D: illegal inputs" );
    }

  }

  void ProductTopology::lineProduct3d( const int dim0 ,
				       const int entity0 ,
				       const int dim1 ,
				       const int entity1 ,
				       const int dim2 ,
				       const int entity2 ,
				       int &resultdim ,
				       int &resultentity ) 
  {
    // on vertex
    if (dim0 == 0 && dim1 == 0 && dim2 == 0) {
      resultdim = 0;
      if (entity0 == 0 && entity1 == 0 && entity2 == 0 ) {
	resultentity = 0;
      }
      else if (entity0 == 0 && entity1 == 0 && entity2 == 1 ) {
	resultentity = 4;
      }
      else if (entity0 == 0 && entity1 == 1 && entity2 == 0 ) {
	resultentity = 3;
      }
      else if (entity0 == 0 && entity1 == 1 && entity2 == 1 ) {
	resultentity = 7;
      }
      else if (entity0 == 1 && entity1 == 0 && entity2 == 0) {
	resultentity = 1;
      }
      else if (entity0 == 1 && entity1 == 0 && entity2 == 1) {
	resultentity = 5;
      }
      else if (entity0 == 1 && entity1 == 1 && entity2 == 0) {
	resultentity = 2;
      }
      else if (entity0 == 1 && entity1 == 1 && entity2 == 1) {
	resultentity = 6;
      }
      else {
	TEUCHOS_TEST_FOR_EXCEPTION( true , std::invalid_argument ,
			    "Intrepid::ProductTopology::lineProduct3D: illegal inputs" );
      }
    }
    // LINES
    // z coord is on line, other two on vertex, this makes an ascending vertical edge
    else if (dim0 == 0 && dim1 == 0 && dim2 == 1) {
      resultdim = 1;
      if (entity0 == 0 && entity1 == 0 && entity2 == 0) {
	resultentity = 8;
      }
      else if (entity0 == 0 && entity1 == 1 && entity2 == 0) {
	resultentity = 11;
      }
      else if (entity0 == 1 && entity1 == 0 && entity2 == 0) {
	resultentity = 9;
      }
      else if (entity0 == 1 && entity1 == 1 && entity2 == 0) {
	resultentity = 10;
      }
      else {
	TEUCHOS_TEST_FOR_EXCEPTION( true , std::invalid_argument ,
			    "Intrepid::ProductTopology::lineProduct3D: illegal inputs" );
      }
    }
    // only y coord is on line, other two on vertex, this makes line along y axis
    else if (dim0 == 0 && dim1 == 1 && dim2 == 0) {
      resultdim = 1;
      if (entity0 == 0 && entity1 == 0 && entity2 == 0) {
	resultentity = 3;
      }
      else if (entity0 == 0 && entity1 == 0 && entity2 == 1) {
	resultentity = 7;
      }
      else if (entity0 == 1 && entity1 == 0 && entity2 == 0) {
	resultentity = 1;
      }
      else if (entity0 == 1 && entity1 == 0 && entity2 == 1) {
 	resultentity = 5;
      }
      else {
	TEUCHOS_TEST_FOR_EXCEPTION( true , std::invalid_argument ,
			    "Intrepid::ProductTopology::lineProduct3D: illegal inputs" );
      }
    }
    // x dof is on line, others on vertex.  
    else if (dim0 == 1 && dim1 == 0 && dim2 == 0) {
      resultdim = 1;
      if (entity0 == 0 && entity1 == 0 && entity2 == 0) {
	resultentity = 0;
      }
      else if (entity0 == 0 && entity1 == 0 && entity2 == 1) {
	resultentity = 4;
      }
      else if (entity0 == 0 && entity1 == 1 && entity2 == 0) {
	resultentity = 2;
      }
      else if (entity0 == 0 && entity1 == 1 && entity2 == 1) {
	resultentity = 6;
      }
      else {
	TEUCHOS_TEST_FOR_EXCEPTION( true , std::invalid_argument ,
			    "Intrepid::ProductTopology::lineProduct3D: illegal inputs" );
      }
    }
    // FACES, these require two of the line dimensions to be 1
    else if (dim0 == 0 && dim1 == 1 && dim2 == 1) {
      resultdim = 2;
      if (entity0 == 0 && entity1 == 0 && entity2 == 0) { 
	resultentity = 3;
      }
      else if (entity0 == 1 && entity1 == 0 && entity2 == 0) { 
	resultentity = 1;
      }
      else {
	TEUCHOS_TEST_FOR_EXCEPTION( true , std::invalid_argument ,
			    "Intrepid::ProductTopology::lineProduct3D: illegal inputs" );
      }
    }
    else if (dim0 == 1 && dim1 == 0 && dim2 == 1) {
      resultdim = 2;
      if (entity0 == 0 && entity1 == 0 && entity2 == 0) { 
	resultentity = 0;
      }
      else if (entity0 == 0 && entity1 == 1 && entity2 == 0) { 
	resultentity = 2;
      }
      else {
	TEUCHOS_TEST_FOR_EXCEPTION( true , std::invalid_argument ,
			    "Intrepid::ProductTopology::lineProduct3D: illegal inputs" );
      }
    }
    else if (dim0 == 1 && dim1 == 1 && dim2 == 0) {
      resultdim = 2;
      if (entity0 == 0 && entity1 == 0 && entity2 == 0) { 
	resultentity = 4;
      }
      else if (entity0 == 0 && entity1 == 0 && entity2 == 1) { 
	resultentity = 5;
      }
      else {
	TEUCHOS_TEST_FOR_EXCEPTION( true , std::invalid_argument ,
			    "Intrepid::ProductTopology::lineProduct3D: illegal inputs" );
      }
    }
    // CELL ITSELF
    else if (dim0 == 1 && dim1 == 1 && dim2 == 1) {
      resultdim = 3;
      if (entity0 == 0 && entity1 == 0 && entity2 == 0) {
	resultentity = 0;
      }
      else {
	TEUCHOS_TEST_FOR_EXCEPTION( true , std::invalid_argument ,
			    "Intrepid::ProductTopology::lineProduct3D: illegal inputs" );
      }
    }
  }
}

#if defined(Intrepid_SHOW_DEPRECATED_WARNINGS)
#ifdef __GNUC__
#warning "The Intrepid package is deprecated"
#endif
#endif

