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

/** \file   Intrepid_F0_TRI_C2_FEM_DEFAULT.hpp
    \brief  Header file for the Intrepid::F0_TRI_C2_FEM_DEFAULT class.
    \author Created by P. Bochev and D. Ridzal.
*/

#ifndef INTREPID_F0_TRI_C2_FEM_DEFAULT_HPP
#define INTREPID_F0_TRI_C2_FEM_DEFAULT_HPP

#include "Intrepid_Basis.hpp"
#include "Intrepid_RealSpace.hpp"
#include "Intrepid_Utils.hpp"
#include "Teuchos_Array.hpp"

namespace Intrepid {
  
/** \class Intrepid::Basis_F0_TRI_C2_FEM_DEFAULT
  \brief Implementation of default FEM basis functions of degree 2 for 0-forms on TRI cells. 
  Reconstruction space type is COMPLETE, i.e., 2D linear polynomials. Definition of the DoF set 
  for this basis, its enumeration and the associated local DoF tags are as follows,
                                                                                    
  \verbatim
  =================================================================================================
  |        |                degree-of-freedom-tag              |                                  |
  | DoF Id |---------------------------------------------------|          DoF definition          |
  |        |  subc dim  |  subc id   | subc DoFId |subc num DoF|                                  |
  |========|============|============|============|============|==================================|
  |    0   |     0      |     0      |     0      |     1      |         L_0(u) = u(v_0)          |
  |--------|------------|------------|------------|------------|----------------------------------|
  |    1   |     0      |     1      |     0      |     1      |         L_1(u) = u(v_1)          |
  |--------|------------|------------|------------|------------|----------------------------------|
  |    2   |     0      |     2      |     0      |     1      |         L_2(u) = u(v_2)          |
  |--------|------------|------------|------------|------------|----------------------------------|
  |    3   |     1      |     0      |     0      |     1      |         L_3(u) = u(v_01)         |
  |--------|------------|------------|------------|------------|----------------------------------|
  |    4   |     1      |     1      |     0      |     1      |         L_4(u) = u(v_12)         |
  |--------|------------|------------|------------|------------|----------------------------------|
  |    5   |     1      |     2      |     0      |     1      |         L_5(u) = u(v_20)         |
  |========|============|============|============|============|==================================|
  \endverbatim
  
  \remarks
  \li v_i  is the ith vertex of the TRI cell;
  \li v_ij is the edge midpoint of edge with endpoints {i,j}
  \li DefaultBasisFactory will select this class if the following parameters are specified:

  \verbatim
  |=======================|===================================|
  |  EField               |  FIELD_FORM_0                     |
  |-----------------------|-----------------------------------|
  |  ECell                |  CELL_TRI                         |
  |-----------------------|-----------------------------------|
  |  EReconstructionSpace |  RECONSTRUCTION_SPACE_COMPLETE    |
  |-----------------------|-----------------------------------|
  |  degree               |  2                                |
  |-----------------------|-----------------------------------|
  |  EBasis               |  BASIS_FEM_DEFAULT                |
  |-----------------------|-----------------------------------|
  |  ECoordinates         |  COORDINATES_CARTESIAN            |
  |=======================|===================================|
  \endverbatim
  
*/
template<class Scalar> 
class Basis_F0_TRI_C2_FEM_DEFAULT: public Basis<Scalar> {
  private:
  
  /** \brief Dimension of the space spanned by the basis = number of degrees of freedom.
  */
  int numDof_;
  
  /**\brief Lookup table for the DoF's local enumeration (DoF Id) by its local DoF tag
  */
  Teuchos::Array<Teuchos::Array<Teuchos::Array<int> > > tagToEnum_;
  
  /**\brief Lookup table for the DoF's local DoF tag by its local enumeration (DoF Id)
  */
  Teuchos::Array<LocalDofTag> enumToTag_;
  
  /**\brief "true" if both lookup arrays have been set by initialize()
  */
  bool isSet_;

  public:

  /** \brief Constructor.
  */
  Basis_F0_TRI_C2_FEM_DEFAULT() : numDof_(6), isSet_(false) {};
  
  /** \brief Initializes arrays needed for the lookup of the local enumeration (DoF Id) of a 
    degree-of-freedom by its local DoF tag and the reverse lookup of the DoF tag by DoF Id.
   */
  void initialize();
  
  
  /** \brief Returns a FieldContainer with the values of <var>opertorType</var> applied to the basis
    functions. For admissible <var>operatorType</var> arguments and the format of the output container 
    see LocalForm0::getOperator(FieldContainer<Scalar>&, const Teuchos::Array<Point<Scalar> >&, const EOperator)
    
    \param outputValues   [out]         - FieldContainer of rank 2 or 3 with the computed values
    \param inputPoints     [in]         - evaluation points on the reference cell  
    \param operatorType    [in]         - the operator being applied to the basis function    
    
    \remarks 
    \li Enumeration of Dk (derivatives of total order k) follows the lexicographical order of 
    the partial derivatives; see getDkEnumeration() for details.
    
    \li For linear basis functions all 2nd order and higher derivatives are identically zero. 
    Nevertheless, the output container for D2,...,D10 is still shaped using DkCardinality as an upper
    bound for the last index, i.e., the output container is filled with as many zeroes as there are
    partial derivatives of a particular order; see getDkCardinality. 
  */
  void getValues(FieldContainer<Scalar>&                outputValues,
                 const Teuchos::Array< Point<Scalar> >& inputPoints,
                 const EOperator                        operatorType) const;
  
  
  /** \brief This method is intended for FVD reconstructions and should not be used here. Its 
    invocation will throw an exception. 
  */  
  void getValues(FieldContainer<Scalar>&                outputValues,
                 const Teuchos::Array< Point<Scalar> >& inputPoints,
                 const Cell<Scalar>&                    cell) const;

  
  int getNumLocalDof() const;
  
  
  int getLocalDofEnumeration(const LocalDofTag dofTag);

  
  LocalDofTag getLocalDofTag(int id);

  
  const Teuchos::Array<LocalDofTag> & getAllLocalDofTags();

  
  ECell getCellType() const;

  
  EBasis getBasisType() const;

  
  ECoordinates getCoordinateSystem() const;

  
  int getDegree() const;

};

}// namespace Intrepid

#include "Intrepid_F0_TRI_C2_FEM_DEFAULTDef.hpp"

#endif
