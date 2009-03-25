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
#include "Intrepid_FieldContainer.hpp"
#include "Intrepid_Types.hpp"

namespace Intrepid {

/** \class Intrepid::Cubature
    \brief Defines the base class for cubature (integration) rules in Intrepid.

    Cubature template (rule) consists of cubature points and cubature weights.
    Intrepid provides a small collection of frequently used cubature rule templates
    for FEM reconstructions on simplices (edge, tri, tet) and the pyramid cell,
    defined in the derived classes of CubatureDirect.

    For quad, hex, and triprism cells cubature templates are tensor products of CubatureDirect
    templates. The tensor-product cubatures are defined in the derived class CubatureTensor.
*/
template<class Scalar, class ArrayType = FieldContainer<Scalar> >
class Cubature {
  private:

  public:

  Cubature() {}

  virtual ~Cubature() {}


  /** \brief Returns cubature points and weights
             (return arrays must be pre-sized/pre-allocated).

      \param cubPoints       [out]     - Array containing the cubature points.
      \param cubWeights      [out]     - Array of corresponding cubature weights.
  */
  virtual void getCubature(ArrayType & cubPoints,
                           ArrayType & cubWeights) const = 0;


  /** \brief Returns the number of cubature points.
  */
  virtual int getNumPoints() const = 0;


  /** \brief Returns dimension of the integration domain.
  */
  virtual int getDimension() const = 0;


  /** \brief Returns algebraic accuracy (e.g. max. degree of polynomial
             that is integrated exactly). For tensor-product or sparse
             rules, algebraic accuracy for each coordinate direction
             is returned.

             Since the size of the return argument need not be known
             ahead of time, other return options are possible, depending
             on the type of the cubature rule. 
  */
  virtual void getAccuracy(std::vector<int> & accuracy) const = 0;

}; // class Cubature 

}// end namespace Intrepid

#endif
