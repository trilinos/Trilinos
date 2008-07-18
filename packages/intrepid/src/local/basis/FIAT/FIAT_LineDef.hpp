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

/** \file FIAT_LineDef.hpp
    \brief  The 1d reference element base class used by FIAT,
        plus some particular lines
    \author Created by R. Kirby.
*/

#ifndef FIAT_LINEDEF_HPP
#define FIAT_LINEDEF_HPP

namespace FIAT
{
  template<class Scalar>
  RefCountPtr<SerialDenseMatrix<int,Scalar> > 
  Line<Scalar>::_computeNormals( RefCountPtr<SerialDenseMatrix<int,Scalar> >verts )
  {
    RefCountPtr<SerialDenseMatrix<int,Scalar> > normals =
      rcp( new SerialDenseMatrix<int,Scalar>(2,1) );
    
    if ( (*verts)(0,0) < (*verts)(1,0) ) {
      (*normals)(0,0) = -1.0;
      (*normals)(1,0) = 1.0;
    }
    else {
      (*normals)(0,0) = 1.0;
      (*normals)(1,0) = -1.0;
    }
    
    return normals;
  }

  template<class Scalar>
  RefCountPtr<SerialDenseMatrix<int,Scalar> > 
  Line<Scalar>::_computeScaledNormals( RefCountPtr<SerialDenseMatrix<int,Scalar> > verts )
  {
    RefCountPtr<SerialDenseMatrix<int,Scalar> > normals =
      rcp( new SerialDenseMatrix<int,Scalar>(2,1) );
    
    (*normals)(0,0) = (*verts)(0,0) - (*verts)(1,0);
    (*normals)(1,0) = (*verts)(1,0) - (*verts)(0,0);

    return normals;
  }

  template<class Scalar>
  RefCountPtr<SerialDenseMatrix<int,Scalar> > 
  Line<Scalar>::_computeTangents( RefCountPtr<SerialDenseMatrix<int,Scalar> > verts )
  {
    RefCountPtr<SerialDenseMatrix<int,Scalar> > tangents =
      rcp( new SerialDenseMatrix<int,Scalar>(1,1) );

    (*tangents)(0,0) = (*verts)(0,0) < (*verts)(1,0) ? 1.0 : -1.0;
    
    return tangents;
  }

  template<class Scalar>
  RefCountPtr<SerialDenseMatrix<int,Scalar> > 
  Line<Scalar>::_computeScaledTangents( RefCountPtr<SerialDenseMatrix<int,Scalar> > verts )
  {
    RefCountPtr<SerialDenseMatrix<int,Scalar> > tangents =
      rcp( new SerialDenseMatrix<int,Scalar>(1,1) );

    (*tangents)(0,0) = (*verts)(1,0) - (*verts)(0,0);
    
    return tangents;
  }

  template<class Scalar>
  RefCountPtr<SerialDenseMatrix<int,Scalar> >
  Line<Scalar>::_computeFaceTangents( )
  {
    return rcp( new SerialDenseMatrix<int,Scalar>(0,0) );
  }

  template<class Scalar>
  RefCountPtr<Array<Array<Scalar> > >
  Line<Scalar>::_computeMeasure( RefCountPtr<SerialDenseMatrix<int,Scalar> > verts )
  {
    RefCountPtr<Array<Array<Scalar> > > measure =
      rcp( new Array<Array<Scalar> >(2) );

    (*measure)[0].resize( 2 );
    (*measure)[1].resize( 1 );

    (*measure)[0][0] = 1.0;
    (*measure)[0][1] = 1.0;
    Scalar diff = (*verts)(1,0) - (*verts)(0,0);
    (*measure)[1][0] = diff > 0 ? 

    return measure;
  }

  template<class Scalar>
  RefCountPtr<SerialDenseMatrix<int,Scalar> > 
  Line<Scalar>::makeLattice( const int order ,
			     const int interior ) const
  {
    int numPoints = std::max<int>( order + 1 - 2 * interior , 0 );
    
    RefCountPtr<SerialDenseMatrix<int,Scalar> > pts = 
      rcp( new SerialDenseMatrix<int,Scalar>(numPoints,1) );

    RefCountPtr<SerialDenseMatrix<int,Scalar> > vs_ptr = this->getVertices();
    SerialDenseMatrix<int,Scalar> &vs = *vs_ptr;

    if (numPoints > 0) {
      Scalar h = (vs(1,0) - vs(0,0)) / (Scalar) order;
      
      for (int i=order;i<numPoints-order;i++) {
	(*pts)(i-order,0) = vs(0,0) + i * h;
      }
    }

    return pts;
   
  }

  template<class Scalar>
  RefCountPtr<SerialDenseMatrix<int,Scalar> > 
  Line<Scalar>::makePoints( const int topdim , const int facet , const int n )  const
  {
    const int sd = this->getSpatialDimension();

    TEST_FOR_EXCEPTION( topdim < 0 || topdim > sd ,
			std::invalid_argument ,
			">>> ERROR (FIAT::Line::makePoints): invalid dimension" );
    TEST_FOR_EXCEPTION( facet < 0 || facet >= this->getNumFacets( topdim ) ,
			std::invalid_argument ,
			">>> ERROR (FIAT::Line::makePoints): invalid facet number" );
    switch( topdim ) {
    case 0:
      return this->getVerticesOfSubcomplex( 0 , facet );
      break;
    case 1:
      return this->makeLattice( n , 1 );
      break;
    default:
      TEST_FOR_EXCEPTION( true ,
			  std::invalid_argument ,
			  ">>> ERROR (FIAT::Line::makePoints): invalid argument" );
    }
  }

}
#endif
