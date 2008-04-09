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

/** \file   Intrepid_F0_TRI_C1_FEM_DEFAULT.hpp
    \brief  Header file for the Intrepid::F0_TRI_C1_FEM_DEFAULT class.
    \author Created by P. Bochev and D. Ridzal.
*/

#ifndef INTREPID_F0_TRI_C1_FEM_DEFAULT_HPP
#define INTREPID_F0_TRI_C1_FEM_DEFAULT_HPP

#include "Intrepid_Basis.hpp"
#include "Intrepid_RealSpace.hpp"
#include "Teuchos_Array.hpp"

namespace Intrepid {
  
/** \class Intrepid::Basis_F0_TRI_C1_FEM_DEFAULT
    \brief Implementation of default FEM basis functions of degree 1 for 0-forms on 
           TRI cells. Reconstruction space type is COMPLETE, i.e., linear polynomials.
*/
template<class Scalar> 
class Basis_F0_TRI_C1_FEM_DEFAULT: public Basis<Scalar> {
  private:
  
  static Teuchos::Array<Teuchos::Array<Teuchos::Array<int> > > tagToEnum_;
  static Teuchos::Array<LocalDofTag> enumToTag_;
  static bool isSet_;

  public:

  /** \brief Constructor.
  */
  Basis_F0_TRI_C1_FEM_DEFAULT();
    
  void getValues(VarContainer<Scalar>&                  outputValues,
                 const Teuchos::Array< Point<Scalar> >& inputPoints,
                 const EOperator                        operatorType) const;
    
  void getValues(VarContainer<Scalar>&                  outputValues,
                 const Teuchos::Array< Point<Scalar> >& inputPoints,
                 const MultiCell<Scalar>&               mCell) const;

  int getLocalDofEnumeration(const LocalDofTag dofTag) const;

  LocalDofTag getLocalDofTag(int id) const;

  const Teuchos::Array<LocalDofTag> & getAllLocalDofTags() const;

  ECell getCellType() const;

  EBasis getBasisType() const;

  ECoordinates getCoordinateSystem() const;

  int getDegree() const;

};

}// namespace Intrepid

#include "Intrepid_F0_TRI_C1_FEM_DEFAULTDef.hpp"

#endif
