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

/** \file   Intrepid_CubatureTensor.hpp
\brief  Header file for the Intrepid::CubatureTensor class.
\author Created by P. Bochev and D. Ridzal.
*/

#ifndef INTREPID_CUBATURE_TENSOR_HPP
#define INTREPID_CUBATURE_TENSOR_HPP

#include "Intrepid_ConfigDefs.hpp"
#include "Intrepid_Cubature.hpp"
#include "Intrepid_CubatureDirect.hpp"
#include "Intrepid_MultiCell.hpp"
#include "Teuchos_TestForException.hpp"


namespace Intrepid {

/** \class Intrepid::CubatureTensor
  \brief Defines tensor-product cubature (integration) rules in Intrepid.

  Tensor-product rules are defined for quads, hexes, and triprisms.
  They are assembled as tensor products of direct cubature rules
  (see class CubatureDirect).

  All rules are defined on a reference cell and can be mapped to physical space
  cells by the methods available in the MultiCell class.
*/
template<class Scalar>
class CubatureTensor : public Intrepid::Cubature<Scalar> {
  private:
  
  public:

  ~CubatureTensor() {}

  /** \brief Returns number of cubature points, cubature points and weights.

    \param numCubPoints    [out]     - Number of cubature points.
    \param cubPoints       [out]     - Vector containing the cubature points.
    \param cubWeights      [out]     - Vector of corresponding cubature weights.
    \param cellType         [in]     - Type of cell on which the cubature rule is defined.
    \param degree           [in]     - In general, represents the degree of polynomials that are integrated
                                       exactly by this cubature rule. For certain derived classes,
                                       <var>degree</var> is a hash code, whose meaning is purely contextual,
                                       see classes CubatureTensorVar and CubatureTensorSparse.
  */
  virtual void getCubature(int &                            numCubPoints,
                           Teuchos::Array< Point<Scalar> >& cubPoints,
                           Teuchos::Array<Scalar>&          cubWeights,
                           const ECell                      cellType,
                           const int                        degree) const;

  /** \brief Returns the number of cubature points.

    \param cellType         [in]     - Type of cell on which the cubature rule is defined.
    \param degree           [in]     - In general, represents the degree of polynomials that are integrated
                                       exactly by this cubature rule. For certain derived classes,
                                       <var>degree</var> is a hash code, whose meaning is purely contextual,
                                       see classes CubatureTensorVar and CubatureTensorSparse.
  */
  virtual int getNumPoints(const ECell cellType,
                           const int   degree) const;

}; // end class CubatureTensor 


} // end namespace Intrepid


// include templated definitions
#include <Intrepid_CubatureTensorDef.hpp>

#endif
