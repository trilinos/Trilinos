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
#include "Intrepid_FieldContainer.hpp"
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
  <strong>DoF Id</strong> number. Every DoF Id is further assigned a unique <strong>DoF tag</strong>
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
    
  /** \brief Fills <var>outputValues</var> with values of <var>primOp</var>, acting on a set of FEM 
  basis functions, at a set of <strong>reference cell</strong> points <var>inputPoints</var>. Rank of the
  output FieldContainer depends on the ranks of the basis fields and the specified operator and
  ranges from 2 to 5. The first two dimensions are always the number of points in <var>inputPoints</var> 
  and the number of basis functions in the basis set. See FieldContainer::resize(int,int,EField,EOperator,int)
  for a complete list of the admissible combinations of primitive operators and fields and the 
  resulting container configurations.
  
  \remarks
  Because FEM reconstruction relies on COMPLETE or INCOMPLETE polynomials which are smooth functions,
  the admissible <var>primOp</var> always include VALUE, D0, D1,...,D10. When derivative order exceeds
  the polynomial degree, <var>outputValues</var> is filled with the appropriate number of zeros. 
 
      \param outputValues   [out]         - FieldContainer with the computed values (see implementation
                                            for index ordering)
      \param inputPoints     [in]         - evaluation points on the reference cell  
      \param primOp          [in]         - primitive operator acting on the basis function
  */    
  virtual void getValues(FieldContainer<Scalar>&                outputValues,
                         const Teuchos::Array< Point<Scalar> >& inputPoints,
                         const EOperator                        primOp) const = 0;
  
    
  /** \brief Fills <var>outputValues</var> with values of a set of FVD basis functions, at a set of 
  <strong>physical cell</strong> points <var>inputPoints</var>. Rank of the output FieldContainer 
  equals 2 + rank of the basis field and ranges from 2 to 4. The first two dimensions are always 
  the number of points in <var>inputPoints</var> and the number of basis functions in the basis set.  
  
  \remarks
  In this method the range of the admissible primitive operators is restricted to VALUE because FVD
  FVD reconstruction relies on "basis" functions defined directly on the physical cell using 
  BROKEN polynomial spaces, i.e., discontinuous functions.
    
      \param outputValues   [out]         - FieldContainer with the computed values 
      \param inputPoints     [in]         - evaluation points on the physical cell  
      \param cell            [in]         - the physical cell
  */    
  virtual void getValues(FieldContainer<Scalar>&                outputValues,
                         const Teuchos::Array< Point<Scalar> >& inputPoints,
                         const Cell<Scalar>&                    cell) const = 0;

  
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
