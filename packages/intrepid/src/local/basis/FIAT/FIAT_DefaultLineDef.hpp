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

/** \file FIAT_DefaultLineDef.hpp
    \brief  The [-1,1] reference element used by FIAT,
    \author Created by R. Kirby.
*/

#ifndef FIAT_DEFAULTLINEDEF_HPP
#define FIAT_DEFAULTLINEDEF_HPP

namespace FIAT
{
  template<class Scalar>
  RefCountPtr<SerialDenseMatrix<int,Scalar> > DefaultLine<Scalar>::_getVertices()
  {
    RefCountPtr<SerialDenseMatrix<int,Scalar> > vs = rcp( new SerialDenseMatrix<int,Scalar>(2,1) );
    (*vs)(0,0) = -1.0;
    (*vs)(1,0) = 1.0;
    return vs;
  }

  template<class Scalar>
  RefCountPtr<Array<Array<Array<int> > > > DefaultLine<Scalar>::_getTopology()
  {
    RefCountPtr<Array<Array<Array<int> > > > top = rcp( new Array<Array<Array<int> > >(2) );
    Array<Array<Array<int> > > &t = *top;

    // index 0 is vertices, of which there are two.
    t[0].resize( 2 );
    // each vertex has exactly one vertex
    t[0][0].resize( 1 );
    t[0][1].resize( 1 );
    t[0][0][0] = 0;
    t[0][1][0] = 1;

    // index 1 is for edges, of which there is 1
    t[1].resize( 1 );
    // that edge contains two vertices, which are vertices 0 and 1.
    t[1][0].resize( 2 );
    t[1][0][0] = 0;
    t[1][0][1] = 1;

    return top;
  }
}

#endif
