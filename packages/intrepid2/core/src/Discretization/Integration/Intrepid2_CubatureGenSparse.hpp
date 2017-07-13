// @HEADER
// ************************************************************************
//
//                           Intrepid2 Package
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
// Questions? Contact Kyungjoo Kim  (kyukim@sandia.gov), or
//                    Mauro Perego  (mperego@sandia.gov)
//
// ************************************************************************
// @HEADER

/** \file   Intrepid_CubatureGenSparse.hpp
    \brief  Header file for the Intrepid2::CubatureGenSparse class.
    \author Created by P. Bochev, D. Ridzal, and M. Keegan.
*/

#ifndef INTREPID2_CUBATURE_GEN_SPARSE_HPP
#define INTREPID2_CUBATURE_GEN_SPARSE_HPP

#include "Intrepid2_ConfigDefs.hpp"
#include "Intrepid2_Cubature.hpp"
#include "Intrepid2_CubatureDirectLineGauss.hpp"
#include "Intrepid2_CubatureSparseHelper.hpp"
#include "Teuchos_Assert.hpp"

/** \def INTREPID2_CUBATURE_SPARSE2D_GAUSS_MAX
  \brief The maximum degree of the polynomial that can be integrated exactly by
         a generalized sparse rule of the Gauss(-Legendre) type in 2D.
*/
#define INTREPID2_CUBATURE_GENSPARSE_GAUSS_MAX 17


namespace Intrepid2{


template<class Scalar, int dimension_, class ArrayPoint=FieldContainer<Scalar>, class ArrayWeight = ArrayPoint>
class CubatureGenSparse : public Intrepid2::Cubature<Scalar,ArrayPoint,ArrayWeight> {
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
  virtual void getCubature(ArrayPoint  & cubPoints,
                           ArrayWeight & cubWeights) const;

  /** \brief Returns cubature points and weights.
              Method for physical space cubature, throws an exception.

       \param cubPoints             [out]        - Array containing the cubature points.
       \param cubWeights            [out]        - Array of corresponding cubature weights.
       \param cellCoords             [in]        - Array of cell coordinates
  */
  virtual void getCubature(ArrayPoint& cubPoints,
                           ArrayWeight& cubWeights,
                           ArrayPoint& cellCoords) const;

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


} // end namespace Intrepid2


// include templated definitions
#include <Intrepid2_CubatureGenSparseDef.hpp>

#endif
