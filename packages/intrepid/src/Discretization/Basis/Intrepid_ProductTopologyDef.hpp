// @HEADER
// ************************************************************************
//
//                           Intrepid Package
//                 Copyright (2007) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Pavel Bochev (pbboche@sandia.gov) or
//                    Denis Ridzal (dridzal@sandia.gov) or
//                    Robert Kirby (robert.c.kirby@ttu.edu)
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
	TEST_FOR_EXCEPTION( true , std::invalid_argument ,
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
	TEST_FOR_EXCEPTION( true , std::invalid_argument ,
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
	TEST_FOR_EXCEPTION( true , std::invalid_argument ,
			    "Intrepid::ProductTopology::lineProduct2D: illegal inputs" );
      }
    }
    else if (dim0 == 1 && dim1 == 1) {
      resultdim = 2;
      if (entity0 == 0 && entity1 == 0) {
	resultentity = 0;
      }
      else {
	TEST_FOR_EXCEPTION( true , std::invalid_argument ,
			    "Intrepid::ProductTopology::lineProduct2D: illegal inputs" );
      }
    }
    else {
      TEST_FOR_EXCEPTION( true , std::invalid_argument ,
			  "Intrepid::ProductTopology::lineProduct2D: illegal inputs" );
    }

  }
}
