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

/** \file   Intrepid_VarContainer.hpp
\brief  Header file for utility class to provide variable containers. VarContainer is used to 
store values of vector and scalar fields and their derivatives at a set of points.
\author Created by P. Bochev and D. Ridzal.
*/

#ifndef INTREPID_VARCONTAINER_HPP
#define INTREPID_VARCONTAINER_HPP

#include "Intrepid_ConfigDefs.hpp"
#include "Intrepid_Types.hpp"
#include "Intrepid_LexContainer.hpp"
#include "Teuchos_Array.hpp"
#include "Teuchos_oblackholestream.hpp"
#include "Teuchos_TestForException.hpp"

namespace Intrepid {
  
  /** \class Intrepid::VarContainer
  \brief Implementation of a templated variable container to store and access values of an operator 
  applied to a <strong>set of fields</strong> (scalar, vector or tensor), evaluated at a 
  <strong>set of points</strong>  on a <\strong>single</strong> cell. 
  
  Provides more specialized storage option than a LexContainer by allowing access to values by point 
  index, function index, function coordinate, and enumeration of derivatives of order k. 
  Values can also be accessed by their enumeration value using an overloaded [] operator. 
  Uses LexContainer for actual storage.     
  */
  template<class Scalar>
  class VarContainer : public LexContainer<Scalar> {
  private:
    
    /** \brief Type of the field whose values will be stored
    */
    EField fieldType_;
    
    /** \brief Type of the operator that is applied to the field
    */
    EOperator operatorType_;
    
    /** \brief The space dimension in which this container is defined
    */
    int spaceDim_;
    
  public:
      
    /** \brief Default destructor.
    */
    ~VarContainer() {};
    
    
    /** \brief Default constructor: LexContainer with zero capacity, dummy field, operator & dim.
    */
    VarContainer() : LexContainer<Scalar>::LexContainer() {
      fieldType_    = FIELD_FORM_0;
      operatorType_ = OPERATOR_VALUE;
      spaceDim_     = 1;
    };
    
    
    /** \brief Copy constructor.
    */
    VarContainer(const VarContainer& right);
    
    
    /** \brief Define VarContainer to store values based on number of points, number of fields
               field type, applied operator type and space dimension. 
      
      \param numPoints       [in]        - number of evaluation points
      \param numFields       [in]        - number of fields that will be evaluated
      \param fieldType       [in]        - type of the field that will be evaluated
      \param operatorType    [in]        - type of the operator that will be applied to the field
      \param spaceDim        [in]        - dimension of the ambient space
    */
    VarContainer(const int         numPoints,
                 const int         numFields,
                 const EField      fieldType,
                 const EOperator   operatorType,
                 const int         spaceDim);
    
    
    /** \brief Returns number of evaluation points. Same as getIndexBound(0).
    */
    int getNumPoints() const;
    
    
    /** \brief Returns number of fields that are being evaluated. Same as getIndexBound(1).
    */
    int getNumFields() const;
    
    
    /** \brief Returns type of the fields that are being evaluated.
    */
    EField getFieldType() const;
    
    
    /** \brief Returns type of the operator that is applied to the fields.
    */
    EOperator getOperatorType() const;
    
    
    /** \brief Returns dimension of the ambient space where fields are defined.
    */
    int getSpaceDim() const;
    
        
    /** \brief Returns the enumeration of a value based on its indices. The number, type  
      and bounds of the indices depend on the number of points, number of fields, the rank of the
      field and the rank of the operator applied to the field. They are summarized below:
      \verbatim
      |--------------------|-------------------|-------------------|------------------------|
      |operator/field rank |       rank 0      | rank 1 2D/3D      | rank 2 2D/3D           |
      |--------------------|-------------------|-------------------|------------------------|
      |       VALUE        | [P][F]            | [P][F][D]         | [P][F][D][D]           |
      |--------------------|-------------------|-------------------|------------------------|
      |     GRAD, D1       | [P][F][D]         | [P][F][D][D]      | [P][F][D][D][D]        |
      |--------------------|-------------------|-------------------|------------------------|
      |        CURL        | [P][F][D](undef3D)| [P][F]/[P][F][D]  | [P][F][D]/[P][F][D][D] |
      |--------------------|-------------------|-------------------|------------------------|
      |        DIV         | [P][F][D](only 1D)| [P][F]            | [P][F][D]              |
      |--------------------|-------------------|-------------------|------------------------|
      |    D1,D2,..,D10    | [P][F][K]         | [P][F][D][K]      | [P][F][D][D][K]        |
      |--------------------|-------------------|-------------------|------------------------|

      Legend:
        P -> point index            range 0 <= P < numPoints
        F -> field index            range 0 <= F < numFields
        D -> field component index  range 0 <= D < spaceDim
        K -> enumeration of Dk      range 0 <= K < DkCardinality
      \endverbatim
      
      \note If operator is one of the k-th derivative operators Dk, the enumeration index K of
      a particular derivative in the set of all kth order derivatives can be obtained from the
      multiplicities of dx, dy and dz using global function getDkEnumeration.
     
      \note The rank of the VarContainer equals the number of indices needed to tag its values 
      and ranges from 2 to 5. The input array must contain exactly the same number of 
      indices, and each index value must be within its admissible bounds as defined above.      
      
      \param multiIndex      [in]     - a set of indices whose enumeration we want
    */
    int getEnumeration(const Teuchos::Array<int>& multiIndex) const;
    
    
    /** \brief Resets VarContainer to store values of an operator applied to a set of numFields
      fields, evaluated at numPoints points.
      
      \param numPoints       [in]        - number of evaluation points
      \param numFields       [in]        - number of fields that will be evaluated
      \param fieldType       [in]        - type of the field whose basis will be evaluated
      \param operatorType    [in]        - type of the operator that will be applied to the basis
      \param spaceDim        [in]        - dimension of the ambient space
    */
    void reset(const int       numPoints,
               const int       numFields,
               const EField    fieldType,
               const EOperator operatorType,
               const int       spaceDim);
    
    
    /** \brief Assignment operator <var>*this = right</var>.
    */
    VarContainer& operator  = (const VarContainer& right);
    
    
    /** \brief Prints basic facts about the container: size, rank, index bounds, number of
      points, number of fields, operator type and space dimension
    */
    void print(std::ostream& os) const;

    
    
  }; // end class VarContainer
  
  //===========================================================================//
  //                                                                           //
  //                Function declarations related to VarContainer              //
  //                                                                           //
  //===========================================================================//
  
  /** \relates VarContainer
  Outputs a formated stream with VarContainer data. For debugging purposes.
  */
  template<class Scalar>
  std::ostream& operator <<  (std::ostream& os, const VarContainer<Scalar>& container);
  
} // end namespace Intrepid

// include templated definitions
#include <Intrepid_VarContainerDef.hpp>

#endif
