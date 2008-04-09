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

/** \file   Intrepid_CubatureDirect.hpp
    \brief  Header file for the Intrepid::CubatureDirect class.
    \author Created by P. Bochev and D. Ridzal.
*/

#ifndef INTREPID_CUBATURE_DIRECT_HPP
#define INTREPID_CUBATURE_DIRECT_HPP

#include "Intrepid_ConfigDefs.hpp"
#include "Intrepid_Cubature.hpp"
#include "Intrepid_MultiCell.hpp"
#include "Teuchos_TestForException.hpp"


namespace Intrepid {

/** \class Intrepid::CubatureDirect
    \brief Defines direct cubature (integration) rules in Intrepid.

    Cubature template (rule) consists of cubature points and cubature weights.
    Intrepid provides a small collection of frequently used cubature rule templates
    for FEM reconstructions on simplices (edge, tri, tet) and pyramid cells.

    For quad, hex, and triprism rules, see tensor-product rules
    defined in the class CubatureTensor, and its derived classes.

    Cubature rules for simplices and the pyramid are stored in the
    <var>cubature_data_</var> array.

    All templates are defined on a reference cell and can be mapped to physical space
    cells by the methods available in the MultiCell class.
*/
template<class Scalar>
class CubatureDirect : public Intrepid::Cubature<Scalar> {
  private:

  /** \brief Type of cell on which the cubature rule is defined.
  */
  ECell cellType_;

  /** \brief The degree of polynomials that are integrated
             exactly by this cubature rule.
  */
  int degree_;
  
  /** \brief Complete set of data defining frequently used direct cubature rules.
  */
  static const CubatureTemplate cubature_data_[CUBATURE_MAX];

  /** \brief Names of templates for frequently used direct cubature rules.
  */
  static const char *cubature_names_[];


  public:

  ~CubatureDirect() {}

  /** \brief Constructor.

      \param cellType         [in]     - Type of cell on which the cubature rule is defined.
      \param degree           [in]     - The degree of polynomials that are integrated
                                         exactly by this cubature rule.
  */
  CubatureDirect(const ECell cellType, const int degree);

  /** \brief Returns number of cubature points, cubature points, and weights
             (return arrays will be sized and memory will be allocated).

      \param numCubPoints    [out]     - Number of cubature points.
      \param cubPoints       [out]     - Array containing the cubature points.
      \param cubWeights      [out]     - Array of corresponding cubature weights.
  */
  void getCubature(int &                            numCubPoints,
                   Teuchos::Array< Point<Scalar> >& cubPoints,
                   Teuchos::Array<Scalar>&          cubWeights) const;

  /** \brief Returns cubature points and weights
             (return arrays must be pre-sized/pre-allocated).

      \param cubPoints       [out]     - Array containing the cubature points.
      \param cubWeights      [out]     - Array of corresponding cubature weights.
  */
  void getCubature(Teuchos::Array< Point<Scalar> >& cubPoints,
                   Teuchos::Array<Scalar>&          cubWeights) const;

  /** \brief Returns the number of cubature points.
  */
  int getNumPoints() const;

  /** \brief Returns integration cell type.
  */
  ECell getCellType() const;

  /** \brief Returns max. degree of polynomials that are integrated exactly.
  */
  int getAccuracy() const;

  /** \brief Returns global index of a direct cubature rule.
  */
  int getIndex(const ECell cellType,
               const int   degree) const;

  /** \brief Returns cubature name.
  */
  const char* getName() const;

  /** \brief Exposes static cubature data.
  */
  static const CubatureTemplate (& exposeData())[CUBATURE_MAX];

}; // end class CubatureDirect 

template<class Scalar>
inline const CubatureTemplate (& CubatureDirect<Scalar>::exposeData())[CUBATURE_MAX] {
  return cubature_data_;
}

} // end namespace Intrepid


// include templated definitions
#include <Intrepid_CubatureDirectDef.hpp>

#endif
