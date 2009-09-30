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

/** \file   Intrepid_CubaturePolylib.hpp
    \brief  Header file for the Intrepid::CubaturePolylib class.
    \author Created by P. Bochev and D. Ridzal.
*/

#ifndef INTREPID_CUBATURE_POLYLIB_HPP
#define INTREPID_CUBATURE_POLYLIB_HPP

#include "Intrepid_ConfigDefs.hpp"
#include "Intrepid_Cubature.hpp"
#include "Intrepid_Polylib.hpp"
#include "Teuchos_TestForException.hpp"


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
