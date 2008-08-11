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

/** \file FIAT_DefaultLine.hpp
    \brief  The [-1,1] reference element used by FIAT,
    \author Created by R. Kirby.
*/

#ifndef FIAT_DEFAULTLINE_HPP
#define FIAT_DEFAULTLINE_HPP

#include "FIAT_Line.hpp"
#include "Teuchos_Array.hpp"
#include "Teuchos_RefCountPtr.hpp"


using Teuchos::Array;
using Teuchos::rcp;

namespace FIAT
{
  /** \class FIAT::DefaultLine
      \brief particular line with [-1,1] coordinates
   */
  template<class Scalar>
  class DefaultLine : public Line<Scalar>
  {
  public:
    /** \brief Constructor */
    DefaultLine( ) : Line<Scalar>( DefaultLine<Scalar>::_getVertices() , DefaultLine<Scalar>::_getTopology() ) {}
    virtual ~DefaultLine( ) {}

  private:
    /** \brief returns the vertices of the array, invoked by constructor */
    static RefCountPtr<SerialDenseMatrix<int,Scalar> > _getVertices();

    /** \brief provides the array structure containing topology, invoked by the constructor */
    static RefCountPtr<Array<Array<Array<int> > > > _getTopology();
  };
}

#include "FIAT_DefaultLineDef.hpp"

#endif
