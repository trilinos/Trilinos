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

/** \file   Intrepid_Basis.hpp
\brief  Header file for the abstract base class Intrepid::Basis.
\author Created by P. Bochev and D. Ridzal.
*/

#ifndef INTREPID_BASIS_HPP
#define INTREPID_BASIS_HPP

#include "Intrepid_ConfigDefs.hpp"
#include "Intrepid_VarContainer.hpp"
#include "Intrepid_MultiCell.hpp"
#include "Teuchos_RCP.hpp"

namespace Intrepid {
  
/** \class Intrepid::Basis
    \brief An abstract base class that defines interface for concrete basis implementations for
           FEM and FVD reconstructions.
*/
template<class Scalar>
class Basis {
  private:

  public:
    
  /** \brief Destructor
  */
  virtual ~Basis() {}
    
  /** \brief Returns multi-indexed value representing values of an operator applied to FEM basis 
             functions, evaluated at a set of reference cell points. FEM reconstruction relies on COMPLETE   
             or INCOMPLETE polynomial spaces defined on reference cells, i.e., smooth local spaces. As a  
             result, admissible operators include VALUE and partial derivatives of all orders up to the  
             maximal order admitted in Intrepid (GRAD, DIV, CURL, D1,...,D10).

             This method returns multi-indexed quantity <var>outputValues</var> with
             multi-index {p,d0,...,dk}, where the indices are as follows:

             \arg p  - index of the evaluation point. Range 0 <= p < np, np = number of evaluation points
             \arg di - order of the ith partial derivative. Range 0 <= di < 2 (the space dimension)
             \arg k  - order of the total derivative. Range 0 <= k < 10 (upper range: VALUE = 0, D10 = 10)

             Rank of the multi-indexed quantity is 1 + k and is determined from the <var>operatorType</var>.
             Example:

             \li OPERATOR_VALUE -> {p}             scalar basis function evaluated at points
             \li OPERATOR_GRAD  -> {p,d0}          1st derivatives of scalar basis function evaluated at points
             \li OPERATOR_D2    -> {p,d0,d1}       2nd derivatives of scalar basis functions evaluated at points

             If total derivative order exceeds polynomial degree, output container is set to zero rank and zero value.
 
      \param outputValues   [out]         - VarContainer with the computed values (see implementation
                                            for index ordering)
      \param inputPoints     [in]         - evaluation points on the reference cell  
      \param operatorType    [in]         - the operator being applied to the basis function
  */    
  virtual void getValues(VarContainer<Scalar>&                  outputValues,
                         const Teuchos::Array< Point<Scalar> >& inputPoints,
                         const EOperator                        operatorType) = 0;
    
    
  /** \brief Returns multi-indexed value representing values of an FVD basis function, evaluated 
             at a set of physical cell points. FVD reconstruction relies on "basis" functions defined 
             directly on the physical cell using BROKEN polynomial spaces, i.e., discontinuous local 
             spaces. As a result, admissible operators are restricted to VALUE, which is the default for 
             this method. Because this method works on actual physical cells, it needs a MultiCell argument
             with their vertex coordinates.
    
      \param outputValues   [out]         - VarContainer with the computed values 
      \param inputPoints     [in]         - evaluation points on the physical cell  
      \param mCell           [in]         - the MultiCell containing the physical cells
  */    
  virtual void getValues(VarContainer<Scalar>&                  outputValues,
                         const Teuchos::Array< Point<Scalar> >& inputPoints,
                         const MultiCell<Scalar>&               mCell) = 0;

  
  /** \brief Returns the local enumeration (Id) of a degree of freedom with a given LocalDofTag.

      \param dofTag  [in]  - LocalDofTag encoding the localization of the degree of freedom.

      \return
              - local enumeration (Id)
  */
  virtual int getLocalDofEnumeration(const LocalDofTag dofTag) = 0;


  /** \brief Returns the local degree-of-freedom tag (LocalDofTag) for a given enumeration (Id).

      \param id  [in]  - Degree-of-freedom enumeration (id).

      \return
              - local degree-of-freedom tag
  */
  virtual LocalDofTag getLocalDofTag(int id) = 0;


  /** \brief Returns a Teuchos::Array containing local degree-of-freedom tags (LocalDofTag).

      \param dofTags  [out]  - Teuchos::Array of degree-of-freedom tags.
  */
  virtual void getAllLocalDofTags(Teuchos::Array<LocalDofTag>& dofTags) = 0;

};
  
}// namespace Intrepid

#endif
