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

/** \file FIAT_AffineMap.hpp
    \brief Defines affine maps between two sets of vertices
    \author Created by R. Kirby.
*/

#ifndef FIAT_AFFINEMAP_HPP
#define FIAT_AFFINEMAP_HPP

#include "Teuchos_SerialDenseMatrix.hpp"
#include "Teuchos_SerialDenseVector.hpp"
#include "Teuchos_RefCountPtr.hpp"
#include "Teuchos_SerialDenseSolver.hpp"

using Teuchos::SerialDenseMatrix;
using Teuchos::SerialDenseVector;
using Teuchos::RefCountPtr;
using Teuchos::SerialDenseSolver;

namespace FIAT 
{
  /** \class FIAT::AffineMapping
      \brief constructs affine mappings between two sets of vertices.  This is represented
      by a matrix A and a vector b such that 
      y = A * x + b 
      baseVerts is a m x n matrix, where m is the number of coordinates to define the simplex
         and n is the intrinsic spatial dimension of these coordinates.
         For example, if the base simplex is a triangle embedded in R^3, then m == 3 and n == 3
	 If the base simplex is a triangle in R^2, then m == 3 and n == 2
      targetVerts is an m' x n' matrix, where m' is the number of coordinates to define
         the simplex and n' is the spatial dimension of these coordinates.
      Note that m' == m is required for the transformation to be well-defined, but
         that n' need not equal n (for example, a mapping from the face of a tetrahedron in R^3
         to a reference triangle in R^2).
  */
  template<class Scalar>
  class AffineMap
  {
  public:
    /** \brief Constructor */
    AffineMap( RefCountPtr<SerialDenseMatrix<int,Scalar> > baseVerts ,
	       RefCountPtr<SerialDenseMatrix<int,Scalar> > targetVerts );

    /** \brief Returns a RefCountPtr to the matrix A of the transformation (the Jacobian)
     */
    RefCountPtr<SerialDenseMatrix<int,Scalar> >
    getMatrix( ) const 
    {
      return A_;
    }
 
    /** \brief Returns a RefCountPtr to the vector b of the transformation (the translation)
     */
    RefCountPtr<SerialDenseMatrix<int,Scalar> >
    getVector( ) const
    {
      return b_;
    }

  private:
    RefCountPtr<SerialDenseMatrix<int,Scalar> > A_;
    RefCountPtr<SerialDenseMatrix<int,Scalar> > b_;

  };
}

#include "FIAT_AffineMapDef.hpp"

#endif
