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
// Questions? Contact Robert Kirby (robert.c.kirby@ttu.edu) 
// ************************************************************************
// @HEADER

/** \file FIAT_Line.hpp
    \brief  The 1d reference element base class used by FIAT
    \author Created by R. Kirby.
*/

#ifndef FIAT_LINE_HPP
#define FIAT_LINE_HPP

#include "FIAT_ReferenceElement.hpp"

namespace FIAT
{
  /** \class Line
   */
  template<class Scalar>
  class Line: public ReferenceElement<Scalar>
  {
  public:
    Line( RefCountPtr<SerialDenseMatrix<int,Scalar> > verts ,
	  RefCountPtr<Array<Array<Array<int> > > > topology ) :
      ReferenceElement<Scalar>( 1 , FIAT_LINE , verts , topology ,
				_computeNormals( verts ) ,
				_computeScaledNormals( verts ) ,
				_computeTangents( verts ) ,
				_computeScaledTangents( verts ) ,
				_computeFaceTangents( ) ,
				_computeMeasure( verts ) )
    {
    }
    
    /** \brief	Constructs the coordinates of a lattice of points on the reference element.
	n is the order of the lattice, so that if n == 1, the lattice includes the vertices.
	If n == 2, the lattice includes the vertices and edge midpoints, and so on.
	If the interior flag is other than zero, then that many points from the boundary in
	are skipped.  For example, if n = 3 and interior = 1 on a triangle, only the triangle's
	midpoint is included.
	The return array stores the points id in the first index and
	the spatial coordinate in the second.  
	An exception is raised on an illegal facet number
     */
    virtual RefCountPtr<SerialDenseMatrix<int,Scalar> > 
    makeLattice( const int order ,
		 const int interior = 0 ) const;
    
    /** \brief Constructs the lattice of points on the interior of the given facet
	of the given topological dimension.  
	The resize method of pts may be invoked.
	An exception is raised on illegal facet number.
     */
    virtual RefCountPtr<SerialDenseMatrix<int,Scalar> >
    makePoints( const int topdim , const int facet , const int n ) const;


  private:
    RefCountPtr<SerialDenseMatrix<int,Scalar> > 
    _computeNormals( RefCountPtr<SerialDenseMatrix<int,Scalar> > verts ) ;
			
    RefCountPtr<SerialDenseMatrix<int,Scalar> > 
    _computeScaledNormals( RefCountPtr<SerialDenseMatrix<int,Scalar> > verts );

    RefCountPtr<SerialDenseMatrix<int,Scalar> >
    _computeTangents( RefCountPtr<SerialDenseMatrix<int,Scalar> > verts );

    RefCountPtr<SerialDenseMatrix<int,Scalar> >
    _computeScaledTangents( RefCountPtr<SerialDenseMatrix<int,Scalar> > verts );

    RefCountPtr<SerialDenseMatrix<int,Scalar> >
    _computeFaceTangents( );

    RefCountPtr<Array<Array<Scalar> > >
    _computeMeasure( RefCountPtr<SerialDenseMatrix<int,Scalar> > verts );
  };


}

#include "FIAT_LineDef.hpp"

#endif
