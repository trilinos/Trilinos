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
//                    Robert Kirby (robert.c.kirby@ttu.edu) or
//
// ************************************************************************
// @HEADER

/** \file   Intrepid_ProductTopology.hpp
    \brief  Header file for utility class to take Cartesion products,
            of cell topologies.  Currently just for products of lines.
    \author Created by R. Kirby
*/

#ifndef INTREPID_PRODUCTTOPOLOGY_HPP
#define INTREPID_PRODUCTTOPOLOGY_HPP

#include "Shards_CellTopology.hpp"
#include "Teuchos_TestForException.hpp"
#include <stdexcept>
#include <map>

namespace Intrepid {
  
  /** \class Intrepid::ProductTopology
      \brief Utility class that provides methods for calculating
             distributions of points on different cells
	     
  */
  class ProductTopology {
  public:
    /** \brief 
	\param dim0         [in]  - dimension of first term in cartesian product
	\param entity0      [in]  - entity id of first term in cartesian product
	\param dim1         [in]  - dimension of second term in cartesian product
	\param entity1      [in]  - entity id of second term in cartesian product
	\param resultdim    [out] - dimension of resulting entity
	\param resultentity [out] - id of resulting entity
    */    
    static void lineProduct2d( const int dim0 ,
			       const int entity0 ,
			       const int dim1 ,
			       const int entity1 ,
			       int &resultdim ,
			       int &resultentity )
    {
      // two vertices
      if (dim0 == 0 && dim1 == 0) 
	{
	  resultdim = 0;
	  if (entity0 == 0 && entity1 == 0) 
	    {
	      resultentity = 0;
	    }
	  else if (entity0 == 0 && entity1 == 1) 
	    {
	      resultentity = 3;
	    }
	  else if (entity0 == 1 && entity1 == 0) 
	    {
	      resultentity = 1;
	    }
	  else if (entity0 == 1 && entity1 == 1) 
	    {
	      resultentity = 2;
	    } 
	  else 
	    {
	      TEST_FOR_EXCEPTION( true , std::invalid_argument ,
				  "Intrepid::ProductTopology::lineProduct2D: illegal inputs" );
	    }
	}
      else if (dim0 == 0 && dim1 == 1) 
	{
	  resultdim = 1;
	  if (entity0 == 0 && entity1 == 0) 
	    {
	      resultentity = 3;
	    }
	  else if (entity0 == 1 && entity1 == 0) 
	    {
	      resultentity = 1;
	    }
	  else 
	    {
	      TEST_FOR_EXCEPTION( true , std::invalid_argument ,
				  "Intrepid::ProductTopology::lineProduct2D: illegal inputs" );
	    }
	}
      else if (dim0 == 1 && dim1 == 0) 
	{
	  resultdim = 1;
	  if (entity0 == 0 && entity1 == 0) 
	    {
	      resultentity = 0;
	    }
	  else if (entity0 == 0 && entity1 == 1) 
	    {
	      resultentity = 2;
	    }
	  else 
	    {
	      TEST_FOR_EXCEPTION( true , std::invalid_argument ,
				  "Intrepid::ProductTopology::lineProduct2D: illegal inputs" );
	    }
	}
      else if (dim0 == 1 && dim1 == 1) 
	{
	  resultdim = 2;
	  if (entity0 == 0 && entity1 == 0) 
	    {
	      resultentity = 0;
	    }
	  else 
	    {
	      TEST_FOR_EXCEPTION( true , std::invalid_argument ,
				  "Intrepid::ProductTopology::lineProduct2D: illegal inputs" );
	    }
	}
      else 
	{
	  TEST_FOR_EXCEPTION( true , std::invalid_argument ,
			      "Intrepid::ProductTopology::lineProduct2D: illegal inputs" );
	}
      
    }

    /** \brief 
	\param dim0         [in]  - dimension of first term in cartesian product
	\param entity0      [in]  - entity id of first term in cartesian product
	\param dim1         [in]  - dimension of second term in cartesian product
	\param entity1      [in]  - entity id of second term in cartesian product
	\param dim2         [in]  - dimension of third term in cartesian product
	\param entity2      [in]  - entity id of third term in cartesian product
	\param resultdim    [out] - dimension of resulting entity
	\param resultentity [out] - id of resulting entity
    */    
    static void lineProduct3d( const int dim0 ,
			       const int entity0 ,
			       const int dim1 ,
			       const int entity1 ,
			       const int dim2 ,
			       const int entity2 ,
			       int &resultdim ,
			       int &resultentity )
    {
      // on vertex
      if (dim0 == 0 && dim1 == 0 && dim2 == 0) 
	{
	  resultdim = 0;
	  if (entity0 == 0 && entity1 == 0 && entity2 == 0 ) 
	    {
	      resultentity = 0;
	    }
	  else if (entity0 == 0 && entity1 == 0 && entity2 == 1 ) 
	    {
	      resultentity = 4;
	    }
	  else if (entity0 == 0 && entity1 == 1 && entity2 == 0 ) 
	    {
	      resultentity = 3;
	    }
	  else if (entity0 == 0 && entity1 == 1 && entity2 == 1 ) 
	    {
	      resultentity = 7;
	    }
	  else if (entity0 == 1 && entity1 == 0 && entity2 == 0) 
	    {
	      resultentity = 1;
	    }
	  else if (entity0 == 1 && entity1 == 0 && entity2 == 1) 
	    {
	      resultentity = 5;
	    }
	  else if (entity0 == 1 && entity1 == 1 && entity2 == 0) 
	    {
	      resultentity = 2;
	    }
	  else if (entity0 == 1 && entity1 == 1 && entity2 == 1) 
	    {
	      resultentity = 6;
	    }
	  else 
	    {
	      TEST_FOR_EXCEPTION( true , std::invalid_argument ,
				  "Intrepid::ProductTopology::lineProduct3D: illegal inputs" );
	    }
	}
    // LINES
    // z coord is on line, other two on vertex, this makes an ascending vertical edge
    else if (dim0 == 0 && dim1 == 0 && dim2 == 1) 
      {
	resultdim = 1;
	if (entity0 == 0 && entity1 == 0 && entity2 == 0) 
	  {
	    resultentity = 8;
	  }
	else if (entity0 == 0 && entity1 == 1 && entity2 == 0) 
	  {
	    resultentity = 11;
	  }
	else if (entity0 == 1 && entity1 == 0 && entity2 == 0) 
	  {
	    resultentity = 9;
	  }
	else if (entity0 == 1 && entity1 == 1 && entity2 == 0) 
	  {
	    resultentity = 10;
	  }
	else 
	  {
	    TEST_FOR_EXCEPTION( true , std::invalid_argument ,
				"Intrepid::ProductTopology::lineProduct3D: illegal inputs" );
	  }
      }
      // only y coord is on line, other two on vertex, this makes line along y axis
    else if (dim0 == 0 && dim1 == 1 && dim2 == 0) 
      {
	resultdim = 1;
	if (entity0 == 0 && entity1 == 0 && entity2 == 0) 
	  {
	    resultentity = 3;
	  }
	else if (entity0 == 0 && entity1 == 0 && entity2 == 1) 
	  {
	    resultentity = 7;
	  }
	else if (entity0 == 1 && entity1 == 0 && entity2 == 0) 
	  {
	    resultentity = 1;
	  }
	else if (entity0 == 1 && entity1 == 0 && entity2 == 1) 
	  {
	    resultentity = 5;
	  }
	else 
	  {
	    TEST_FOR_EXCEPTION( true , std::invalid_argument ,
				"Intrepid::ProductTopology::lineProduct3D: illegal inputs" );
	  }
      }
    // x dof is on line, others on vertex.  
    else if (dim0 == 1 && dim1 == 0 && dim2 == 0) 
      {
	resultdim = 1;
	if (entity0 == 0 && entity1 == 0 && entity2 == 0) 
	  {
	    resultentity = 0;
	  }
	else if (entity0 == 0 && entity1 == 0 && entity2 == 1) 
	  {
	    resultentity = 4;
	  }
	else if (entity0 == 0 && entity1 == 1 && entity2 == 0) 
	  {
	    resultentity = 2;
	  }
	else if (entity0 == 0 && entity1 == 1 && entity2 == 1) 
	  {
	    resultentity = 6;
	  }
	else 
	  {
	    TEST_FOR_EXCEPTION( true , std::invalid_argument ,
				"Intrepid::ProductTopology::lineProduct3D: illegal inputs" );
	  }
      }
    // FACES, these require two of the line dimensions to be 1
    else if (dim0 == 0 && dim1 == 1 && dim2 == 1) 
      {
	resultdim = 2;
	if (entity0 == 0 && entity1 == 0 && entity2 == 0) 
	  { 
	    resultentity = 3;
	  }
	else if (entity0 == 1 && entity1 == 0 && entity2 == 0) 
	  { 
	    resultentity = 1;
	  }
	else 
	  {
	    TEST_FOR_EXCEPTION( true , std::invalid_argument ,
				"Intrepid::ProductTopology::lineProduct3D: illegal inputs" );
	  }
      }
    else if (dim0 == 1 && dim1 == 0 && dim2 == 1) 
      {
	resultdim = 2;
	if (entity0 == 0 && entity1 == 0 && entity2 == 0) 
	  { 
	    resultentity = 0;
	  }
	else if (entity0 == 0 && entity1 == 1 && entity2 == 0) 
	  { 
	    resultentity = 2;
	  }
      else 
	{
	  TEST_FOR_EXCEPTION( true , std::invalid_argument ,
			      "Intrepid::ProductTopology::lineProduct3D: illegal inputs" );
	}
      }
    else if (dim0 == 1 && dim1 == 1 && dim2 == 0) 
      {
	resultdim = 2;
	if (entity0 == 0 && entity1 == 0 && entity2 == 0) 
	  { 
	    resultentity = 4;
	  }
	else if (entity0 == 0 && entity1 == 0 && entity2 == 1) 
	  { 
	    resultentity = 5;
	  }
      else 
	{
	  TEST_FOR_EXCEPTION( true , std::invalid_argument ,
			      "Intrepid::ProductTopology::lineProduct3D: illegal inputs" );
	}
      }
     // CELL ITSELF
    else if (dim0 == 1 && dim1 == 1 && dim2 == 1) 
      {
	resultdim = 3;
	if (entity0 == 0 && entity1 == 0 && entity2 == 0) 
	  {
	    resultentity = 0;
	  }
	else 
	  {
	    TEST_FOR_EXCEPTION( true , std::invalid_argument ,
				"Intrepid::ProductTopology::lineProduct3D: illegal inputs" );
	  }
      }
    }

  };
}

#endif
