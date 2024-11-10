// @HEADER
// ************************************************************************
//
//                           Intrepid Package
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
// Questions? Contact Pavel Bochev  (pbboche@sandia.gov)
//                    Denis Ridzal  (dridzal@sandia.gov), or
//                    Kara Peterson (kjpeter@sandia.gov)
//
// ************************************************************************
// @HEADER

/** \file   Intrepid_CubaturePolylib.hpp
    \brief  Header file for the Intrepid::CubaturePolylib class.
    \author Created by P. Bochev and D. Ridzal.
*/

#ifndef INTREPID_CUBATURE_POLYLIB_HPP
#define INTREPID_CUBATURE_POLYLIB_HPP

#include "Intrepid_ConfigDefs.hpp"
#include "Intrepid_Cubature.hpp"
#include "Intrepid_Polylib.hpp"
#include "Teuchos_Assert.hpp"


namespace Intrepid {

/** \class Intrepid::CubaturePolylib
    \brief Utilizes cubature (integration) rules contained in the library Polylib
           (Spencer Sherwin, Aeronautics, Imperial College London) within Intrepid.

           They are based on zeros of Jacobi polynomials,
           e.g. Legendre (alpha=beta=0, default), Chebyshev (alpha=beta=-0.5), etc.
           They are given on the interval [-1,1] and are optimal with respect to
           the following requirements, yielding 4 subclasses:
           \li Gauss             - no restrictions (all integration points are contained in the open interval (-1,1))
           \li Gauss-Radau-Left  - left-most integration point is fixed at -1
           \li Gauss-Radau-Right - right-most integration point is fixed at +1
           \li Gauss-Lobatto     - left-most and right-most integration points are fixed at -1 and +1, respectively
*/
template<class Scalar, class ArrayPoint = FieldContainer<Scalar>, class ArrayWeight = ArrayPoint>
class CubaturePolylib : public Intrepid::Cubature<Scalar,ArrayPoint,ArrayWeight> {
  private:

  /** \brief The degree of polynomials that are integrated
             exactly by this cubature rule.
  */
  int degree_;

  /** \brief Dimension of integration domain.
  */
  int dimension_;

  /** \brief Type of integration points.
  */
  EIntrepidPLPoly poly_type_;

  /** \brief Jacobi parameter alpha.
  */
  Scalar alpha_;

  /** \brief Jacobi parameter beta.
  */
  Scalar beta_;

  /** \brief Cubature name.
  */
  static const char *cubature_name_;


  public:

  ~CubaturePolylib() {}

  /** \brief Constructor.

      \param degree           [in]     - The degree of polynomials that are integrated
                                         exactly by this cubature rule. Default: 0.
  */
  CubaturePolylib(int degree = 0, EIntrepidPLPoly pt_type = PL_GAUSS, Scalar alpha = 0.0, Scalar beta = 0.0);


  /** \brief Returns cubature points and weights
             (return arrays must be pre-sized/pre-allocated).

      \param cubPoints       [out]     - Array containing the cubature points.
      \param cubWeights      [out]     - Array of corresponding cubature weights.
  */
  void getCubature(ArrayPoint  & cubPoints,
                   ArrayWeight & cubWeights) const;

  /** \brief Returns cubature points and weights.
              Method for physical space cubature, throws an exception.

       \param cubPoints             [out]        - Array containing the cubature points.
       \param cubWeights            [out]        - Array of corresponding cubature weights.
       \param cellCoords             [in]        - Array of cell coordinates
  */
  void getCubature(ArrayPoint& cubPoints,
                   ArrayWeight& cubWeights,
                   ArrayPoint& cellCoords) const;

  /** \brief Returns the number of cubature points.
  */
  int getNumPoints() const;

  /** \brief Returns dimension of integration domain.
  */
  virtual int getDimension() const;

  /** \brief Returns max. degree of polynomials that are integrated exactly.
             The return vector has size 1.
  */
  void getAccuracy(std::vector<int> & accuracy) const;

  /** \brief Returns cubature name.
  */
  const char* getName() const;

}; // end class CubaturePolylib

} // end namespace Intrepid

// include templated definitions
#include <Intrepid_CubaturePolylibDef.hpp>

#endif

#if defined(Intrepid_SHOW_DEPRECATED_WARNINGS)
#ifdef __GNUC__
#warning "The Intrepid package is deprecated"
#endif
#endif

