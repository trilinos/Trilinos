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

/** \file   Intrepid_Cubature.hpp
    \brief  Header file for the Intrepid::Cubature class.
    \author Created by P. Bochev and D. Ridzal.
*/

#ifndef INTREPID_CUBATURE_HPP
#define INTREPID_CUBATURE_HPP

#include "Intrepid_ConfigDefs.hpp"
#include "Intrepid_RealSpace.hpp"
#include "Intrepid_Types.hpp"
#include "Teuchos_Array.hpp"


namespace Intrepid {

/** \class Intrepid::Cubature
  \brief Defines the base class for cubature (integration) rules in Intrepid.

  Cubature template (rule) consists of cubature points and cubature weights.
  Intrepid provides a small collection of frequently used cubature rule templates
  for FEM reconstructions on simplices (edge, tri, tet) and the pyramid cell,
  defined in the derived class CubatureDirect.

  For quad, hex, and triprism cells cubature templates are tensor products of CubatureDirect
  templates. The tensor-product cubatures are defined in the derived class CubatureTensor.

  All templates are defined on a reference cell and can be mapped to physical space
  cells by the methods available in the MultiCell class.
*/
template<class Scalar>
class Cubature {
  private:

  public:

  Cubature() {}

  virtual ~Cubature() {}


  /** \brief Returns number of cubature points, cubature points, and weights
             (return arrays will be sized and memory will be allocated).

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
                           const int                        degree) const = 0;


  /** \brief Returns cubature points and weights
             (return arrays must be pre-sized/pre-allocated).

    \param cubPoints       [out]     - Vector containing the cubature points.
    \param cubWeights      [out]     - Vector of corresponding cubature weights.
    \param cellType         [in]     - Type of cell on which the cubature rule is defined.
    \param degree           [in]     - In general, represents the degree of polynomials that are integrated
                                       exactly by this cubature rule. For certain derived classes,
                                       <var>degree</var> is a hash code, whose meaning is purely contextual,
                                       see classes CubatureTensorVar and CubatureTensorSparse.
  */
  virtual void getCubature(Teuchos::Array< Point<Scalar> >& cubPoints,
                           Teuchos::Array<Scalar>&          cubWeights,
                           const ECell                      cellType,
                           const int                        degree) const = 0;


  /** \brief Returns the number of cubature points.

    \param cellType         [in]     - Type of cell on which the cubature rule is defined.
    \param degree           [in]     - In general, represents the degree of polynomials that are integrated
                                       exactly by this cubature rule. For certain derived classes,
                                       <var>degree</var> is a hash code, whose meaning is purely contextual,
                                       see classes CubatureTensorVar and CubatureTensorSparse.
   */
  virtual int getNumPoints(const ECell cell_type,
                           const int   degree) const = 0;

}; // class Cubature 

}// end namespace Intrepid

#endif
