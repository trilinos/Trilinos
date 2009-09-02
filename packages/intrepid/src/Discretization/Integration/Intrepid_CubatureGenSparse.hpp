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
//                    Denis Ridzal (dridzal@sandia.gov).
//
// ************************************************************************
// @HEADER

/** \file   Intrepid_CubatureGenSparse.hpp
    \brief  Header file for the Intrepid::CubatureGenSparse class.
    \author Created by P. Bochev, D. Ridzal, and M. Keegan.
*/

#ifndef INTREPID_CUBATURE_GEN_SPARSE_HPP
#define INTREPID_CUBATURE_GEN_SPARSE_HPP

#include "Intrepid_ConfigDefs.hpp"
#include "Intrepid_Cubature.hpp"
#include "Intrepid_CubatureDirectLineGauss.hpp"
#include "Intrepid_CubatureSparseHelper.hpp"
#include "Teuchos_TestForException.hpp"

/** \def INTREPID_CUBATURE_SPARSE2D_GAUSS_MAX
  \brief The maximum degree of the polynomial that can be integrated exactly by
         a generalized sparse rule of the Gauss(-Legendre) type in 2D.
*/
#define INTREPID_CUBATURE_GENSPARSE_GAUSS_MAX 17


namespace Intrepid{


template<class Scalar, int dimension_, class ArrayType=FieldContainer<Scalar> >
class CubatureGenSparse : public Intrepid::Cubature<Scalar,ArrayType> {
  private:

  int numPoints_;

  const int degree_;

  SGNodes<Scalar, dimension_> grid;

  
  public:

  ~CubatureGenSparse() {}


  CubatureGenSparse(const int degree);

  /** \brief Returns cubature points and weights
             (return arrays must be pre-sized/pre-allocated).

      \param cubPoints       [out]     - Array containing the cubature points.
      \param cubWeights      [out]     - Array of corresponding cubature weights.
  */
  virtual void getCubature(ArrayType & cubPoints,
                           ArrayType & cubWeights) const;

  /** \brief Returns the number of cubature points.
  */
  virtual int getNumPoints() const;

  /** \brief Returns dimension of the integration domain.
  */
  virtual int getDimension() const;

  /** \brief Returns algebraic accuracy (e.g. max. degree of polynomial
             that is integrated exactly).
  */
  virtual void getAccuracy(std::vector<int> & accuracy) const;

}; // end class CubatureGenSparse

// helper functions
template<class Scalar>
inline Scalar Sum(Scalar* list, int first, int last)
{
  Scalar sum = 0;
  for(int i = first; i <= last; i++)
    sum += list[i];
  return sum;
}


} // end namespace Intrepid


// include templated definitions
#include <Intrepid_CubatureGenSparseDef.hpp>

#endif
