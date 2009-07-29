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

/** \file   Intrepid_HCURL_TET_In_FEM.hpp
    \brief  Header file for the Intrepid::HCURL_TET_In_FEM class.
    \author Created by R. Kirby 
*/

#ifndef INTREPID_HCURL_TET_In_FEM_HPP
#define INTREPID_HCURL_TET_In_FEM_HPP

#include "Intrepid_Types.hpp"
#include "Intrepid_Basis.hpp"
#include "Intrepid_PointTools.hpp"
#include "Intrepid_HGRAD_TET_Cn_FEM_ORTH.hpp"
#include "Intrepid_CubatureDirectTetDefault.hpp"
#include "Intrepid_CubatureDirectLineGauss.hpp"
#include "Teuchos_SerialDenseMatrix.hpp"
#include "Teuchos_SerialDenseSolver.hpp"


namespace Intrepid {
  
/** \class  Intrepid::Basis_HCURL_TET_In_FEM
    \brief  Implementation of the default H(curl)-compatible Nedelec (first kind) 
            basis of arbitrary degree  on Tetrahedron cell.  The lowest order space
	    is indexted with 0 rather than one.
            Implements nodal basis of degree n (n>=1) on the reference Tetrahedron cell. The basis has
            cardinality n*(n+1)*(n+3)/2 and spans an INCOMPLETE polynomial space. Basis functions are dual 
            to a unisolvent set of degrees-of-freedom (DoF) defined and enumerated as follows:
  
  \verbatim
  \endverbatim
  
    \remarks
    \li     DefaultBasisFactory will select this class if the following parameters are specified:
  
  \verbatim
  |=======================|===================================|
  |  CellTopology         |  Tetrahedron                      |
  |-----------------------|-----------------------------------|
  |  EFunctionSpace       |  FUNCTION_SPACE_HCURL             |
  |-----------------------|-----------------------------------|
  |  EDiscreteSpace       |  DISCRETE_SPACE_INCOMPLETE        |
  |-----------------------|-----------------------------------|
  |  degree               |  n                                |
  |-----------------------|-----------------------------------|
  |  EBasis               |  BASIS_FEM_FIAT                   |
  |-----------------------|-----------------------------------|
  |  ECoordinates         |  COORDINATES_CARTESIAN            |
  |=======================|===================================|
  \endverbatim
  */
  
template<class Scalar, class ArrayScalar> 
class Basis_HCURL_TET_In_FEM: public Basis<Scalar, ArrayScalar> {
private:
  
  /** \brief  Initializes <var>tagToOrdinal_</var> and <var>ordinalToTag_</var> lookup arrays.
   */
  virtual void initializeTags();

  Basis_HGRAD_TET_Cn_FEM_ORTH<Scalar,FieldContainer<Scalar> > Phis_;
  FieldContainer<Scalar> coeffs_;

public:
  
  /** \brief  Constructor.
   */
  Basis_HCURL_TET_In_FEM( const int n , const EPointType pointType );    
  
  /** \brief  Evaluation of a FEM basis on a <strong>reference Tetrahedron</strong> cell. 
  
              Returns values of <var>operatorType</var> acting on FEM basis functions for a set of
              points in the <strong>reference Tetrahedron</strong> cell. For rank and dimensions of
              I/O array arguments see Section \ref basis_md_array_sec .

      \param  outputValues      [out] - variable rank array with the basis values
      \param  inputPoints       [in]  - rank-2 array (P,D) with the evaluation points
      \param  operatorType      [in]  - the operator acting on the basis functions    
   */
  void getValues(ArrayScalar &          outputValues,
                 const ArrayScalar &    inputPoints,
                 const EOperator        operatorType) const;
  
  
  /**  \brief  FVD basis evaluation: invocation of this method throws an exception.
   */
  void getValues(ArrayScalar &          outputValues,
                 const ArrayScalar &    inputPoints,
                 const ArrayScalar &    cellVertices,
                 const EOperator        operatorType = OPERATOR_VALUE) const;
};

}// namespace Intrepid

#include "Intrepid_HCURL_TET_In_FEMDef.hpp"

#endif
