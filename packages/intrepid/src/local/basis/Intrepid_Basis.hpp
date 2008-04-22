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
#include "Intrepid_Cell.hpp"

namespace Intrepid {
  
/** \class Intrepid::Basis
    \brief An abstract base class that defines interface for concrete basis implementations for
  FEM and FVD reconstructions. A FEM reconstruction approximates fields by COMPLETE or INCOMPLETE
  polynomial spaces; an FVD reconstruction typically uses a BROKEN polynomial space. A particular
  basis of a given reconstruction space is completely defined by a unisolvent set of linear 
  functionals called <strong>degrees-of-freedom</strong>, or DoF, by choosing the basis set to be
  the dual of the DoF set. Because the two sets are isomorphic, enumeration of the DoF set establishes
  enumeration of the basis set which assigns to every basis function and DoF functional a local
  <strong> DoF Id </strong> number. Every DoF Id is further assigned a unique <strong>DoF tag</strong>
  which associates that DoF ID with a particular subcell and establishes the order of the DoF with
  respect to that subcell; see LocalDofTag. Basis functions for FEM reconstructions are defined on 
  standard <strong>reference</strong> cells and the values returned need to be further transformed 
  to physical cells. Basis functions for FVD reconstructions are defined directly on the 
  <strong>physical</strong> cells and values returned don't have to be transformed. 
*/
template<class Scalar>
class Basis {
  private:

  public:
    
  /** \brief Destructor
  */
  virtual ~Basis() {}
    
  /** \brief Returns VarContainer with multi-indexed quantity representing the values of 
  an operator applied to a set of FEM basis functions, evaluated at a set of points in a <strong>reference</strong>
  cell. The rank of the return VarContainer argument depends on the number of points, number of basis, 
  functions, their rank, and the rank of the operator applied to them; see VarContainer::getEnumeration 
  for summary of the admissible index combinations. In particular, the admissible range of the 
  <var>operatorType</var> argument depends on the rank of the basis functions and the space dimension.
  Because FEM reconstruction relies on COMPLETE or INCOMPLETE polynomials, i.e., smooth local spaces,
  the admissible range of <var>operatorType</var> always includes VALUE, D0, D1,...,D10. If derivative 
  order exceeds polynomial degree, output container is set to zero rank and zero value.
 
      \param outputValues   [out]         - VarContainer with the computed values (see implementation
                                            for index ordering)
      \param inputPoints     [in]         - evaluation points on the reference cell  
      \param operatorType    [in]         - the operator being applied to the basis function
  */    
  virtual void getValues(VarContainer<Scalar>&                  outputValues,
                         const Teuchos::Array< Point<Scalar> >& inputPoints,
                         const EOperator                        operatorType) const = 0;
    
    
  /** \brief Returns VarContainer with multi-indexed quantity representing the values of an operator 
  applied to a set of FVD basis functions, evaluated at a set of points in a <strong>physical</strong> cell.  
  Because FVD reconstruction relies on "basis" functions defined directly on the physical cell using 
  BROKEN polynomial spaces, i.e., discontinuous local spaces, the admissible range of the 
  <var>operatorType<var> argument is restricted to VALUE. 
    
      \param outputValues   [out]         - VarContainer with the computed values 
      \param inputPoints     [in]         - evaluation points on the physical cell  
      \param cell            [in]         - the physical cell
  */    
  virtual void getValues(VarContainer<Scalar>&                  outputValues,
                         const Teuchos::Array< Point<Scalar> >& inputPoints,
                         const Cell<Scalar>&                     cell) const = 0;

  
  /** \brief Returns the number of basis functions in a basis. 
    
      \return
              - number of basis funcions in a concrete basis
  */
  virtual int getNumLocalDof() const = 0;
  
  /** \brief Returns the local enumeration (DoF Id) of a degree-of-freedom with a given LocalDofTag.

      \param dofTag  [in]  - The LocalDoF tag assigned to the degree-of-freedom 

      \return               
              - local enumeration (DoF Id) associated with the degree-of-freedom LocalDoFtag
  */
  virtual int getLocalDofEnumeration(const LocalDofTag dofTag) = 0;


  /** \brief Returns the local degree-of-freedom tag (LocalDofTag) for a given enumeration (DoF Id).

      \param id  [in]  - The local enumeration (DoF id) assigned to the degree-of-freedom.

      \return
              - local degree-of-freedom tag associated with the local DoF enumeration (DoF Id)
  */
  virtual LocalDofTag getLocalDofTag(const int id) = 0;


  /** \brief Returns a Teuchos::Array containing all local degree-of-freedom tags (LocalDofTag).

      \param dofTags  [out]  - Teuchos::Array of degree-of-freedom tags.
  */
  virtual const Teuchos::Array<LocalDofTag> & getAllLocalDofTags() = 0;


  /** \brief Returns cell type on which the basis is defined.

      \return
              - cell type
  */
  virtual ECell getCellType() const = 0;


  /** \brief Returns basis type.

      \return
              - basis type
  */
  virtual EBasis getBasisType() const = 0;


  /** \brief Returns type of coordinate system (Cartesian, polar, R-Z, etc.).

      \return
              - coordinate system
  */
  virtual ECoordinates getCoordinateSystem() const = 0;


  /** \brief Returns the degree of a polynomial basis.

      \return
              - if the basis is a polynomial basis, returns degree of max. complete polynomial
                that can be represented by the basis, otherwise -1
  */
  virtual int getDegree() const = 0;

};
  
}// namespace Intrepid

#endif
