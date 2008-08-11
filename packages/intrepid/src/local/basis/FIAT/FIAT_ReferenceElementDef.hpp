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

/** \file FIAT_ReferenceElementDef.hpp
    \brief  Implementation of methods for ReferenceElement class
    \author Created by R. Kirby.
*/

namespace FIAT
{
  template<class Scalar>
  RefCountPtr<SerialDenseMatrix<int,Scalar> > 
  ReferenceElement<Scalar>::getVerticesOfSubcomplex( const int topdim ,
						    const int facet ) const
  {
    TEST_FOR_EXCEPTION( topdim < 0 || topdim > spatialDim_ ,
			std::invalid_argument ,
			">>> ERROR (FIAT::ReferenceElement::getVerticesOfSubcomplex): invalid dimension" );
    TEST_FOR_EXCEPTION( facet < 0 || facet >= getNumFacets( topdim ) ,
			std::invalid_argument ,
			">>> ERROR (FIAT::ReferenceElement::getVerticesOfSubcomplex): invalid facet number" );
    int numVerts = (*topology_)[topdim][facet].size();
    SerialDenseMatrix<int,Scalar> &verts = *verts_;
    RefCountPtr<SerialDenseMatrix<int,Scalar> > verts_sub_ptr
      = rcp( new SerialDenseMatrix<int,Scalar>( numVerts , spatialDim_ ) );
    SerialDenseMatrix<int,Scalar> &verts_sub = *verts_sub_ptr;
    for (int i=0;i<numVerts;i++) {
      int vert_cur = (*topology_)[topdim][facet][i];
      for (int j=0;j<spatialDim_;j++) {
	verts_sub(i,j) = verts( vert_cur , j );
      }
    }
    return verts_sub_ptr;
  }

  template<class Scalar>
  RefCountPtr<SerialDenseMatrix<int,Scalar> > 
  ReferenceElement<Scalar>::getFaceTangents( ) const
  {
    TEST_FOR_EXCEPTION( spatialDim_ != 3 , 
			std::invalid_argument ,
			">>> ERROR (FIAT::ReferenceElement::getFaceTangents): only works in 3d" );
    return faceTangents_;
  }

  template<class Scalar>
  int ReferenceElement<Scalar>::getNumFacets( const int topdim ) const
  {
    TEST_FOR_EXCEPTION( topdim < 0 || topdim > spatialDim_ ,
			std::invalid_argument , 
			">>> ERROR (ReferenceElement::getNumFacets) : invalid dimension" );
    return (*topology_)[topdim].size();
  }

  template<class Scalar>
  Scalar ReferenceElement<Scalar>::getMeasure( const int topdim , const int facet ) const
  {
    TEST_FOR_EXCEPTION( topdim < 0 || topdim > spatialDim_ ,
			std::invalid_argument ,
			">>> ERROR( FIAT::ReferenceElement::getMeasure ): invalid dimension" );
    TEST_FOR_EXCEPTION( topdim < 0 || topdim >= getNumFacets( topdim ) ,
			std::invalid_argument ,
			">>> ERROR( FIAT::ReferenceElement::getMeasure ): invalid facet" );
    return (*measure_)[ topdim ][ facet ];
    
  }


} 
