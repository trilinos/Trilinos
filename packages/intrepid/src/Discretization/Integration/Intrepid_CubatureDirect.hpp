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
template<class Scalar, class ArrayType = FieldContainer<Scalar> >
class CubatureDirect : public Intrepid::Cubature<Scalar,ArrayType> {
  private:

  protected:

  /** \brief The degree of polynomials that are integrated
             exactly by this cubature rule.
  */
  int degree_;

  /** \brief Dimension of integration domain.
  */
  int dimension_;


  public:

  virtual ~CubatureDirect() {}

  CubatureDirect() {}  


  /** \brief Extracts points and weights from cubData.

      \param cubPoints       [out]     - Array containing the cubature points.
      \param cubWeights      [out]     - Array of corresponding cubature weights.
      \param cubData         [in]      - Pointer to raw cubature data.
  */
  virtual void getCubatureData(ArrayType &              cubPoints,
                               ArrayType &              cubWeights,
                               const CubatureTemplate * cubData) const;

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

  /** \brief Returns dimension of integration domain.
  */
  virtual int getDimension() const;

  /** \brief Returns max. degree of polynomials that are integrated exactly.
             The return vector has size 1.
  */
  virtual void getAccuracy(std::vector<int> & accuracy) const;

  /** \brief Returns cubature name.
  */
  virtual const char* getName() const = 0;

  /** \brief Exposes cubature data.
  */
  virtual const CubatureTemplate * exposeCubatureData() const = 0;

  /** \brief Returns maximum cubature accuracy.
  */
  virtual int getMaxAccuracy() const = 0;

}; // end class CubatureDirect 

} // end namespace Intrepid


// include templated definitions
#include <Intrepid_CubatureDirectDef.hpp>

#endif
