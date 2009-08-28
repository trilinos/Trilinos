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
//                    Kara Peterson (kjpeter@sandia.gov).
//
// ************************************************************************
// @HEADER

/** \file   Intrepid_HCURL_HEX_I1_FEM.hpp
    \brief  Header file for the Intrepid::HCURL_HEX_I1_FEM class.
    \author Created by P. Bochev, D. Ridzal and K. Peterson.
 */

#ifndef INTREPID_HCURL_HEX_I1_FEM_HPP
#define INTREPID_HCURL_HEX_I1_FEM_HPP
#include "Intrepid_Basis.hpp"

namespace Intrepid {
  
/** \class  Intrepid::Basis_HCURL_HEX_I1_FEM
    \brief  Implementation of the default H(curl)-compatible FEM basis of degree 1 on Hexahedron cell 
  
            Implements Nedelec basis of degree 1 on the reference Hexahedron cell. The basis has
            cardinality 12 and spans an INCOMPLETE tri-linear polynomial space. Basis functions are dual 
            to a unisolvent set of degrees-of-freedom (DoF) defined and enumerated as follows:
  
  \verbatim
  ===================================================================================================
  |         |           degree-of-freedom-tag table                    |                            |
  |   DoF   |----------------------------------------------------------|       DoF definition       |
  | ordinal |  subc dim    | subc ordinal | subc DoF ord |subc num DoF |                            |
  |=========|==============|==============|==============|=============|============================|
  |    0    |       1      |       0      |       0      |      1      |   L_0(u) = (u.t)(0,-1,-1)  |
  |---------|--------------|--------------|--------------|-------------|----------------------------|
  |    1    |       1      |       1      |       0      |      1      |   L_1(u) = (u.t)(1,0,-1)   |
  |---------|--------------|--------------|--------------|-------------|----------------------------|
  |    2    |       1      |       2      |       0      |      1      |   L_2(u) = (u.t)(0,1,-1)   |
  |---------|--------------|--------------|--------------|-------------|----------------------------|
  |    3    |       1      |       3      |       0      |      1      |   L_3(u) = (u.t)(-1,0,-1)  |
  |---------|--------------|--------------|--------------|-------------|----------------------------|
  |    4    |       1      |       4      |       0      |      1      |   L_4(u) = (u.t)(0,-1,1)   |
  |---------|--------------|--------------|--------------|-------------|----------------------------|
  |    5    |       1      |       5      |       0      |      1      |   L_5(u) = (u.t)(1,0,1)    |
  |---------|--------------|--------------|--------------|-------------|----------------------------|
  |    6    |       1      |       6      |       0      |      1      |   L_6(u) = (u.t)(0,1,1)    |
  |---------|--------------|--------------|--------------|-------------|----------------------------|
  |    7    |       1      |       7      |       0      |      1      |   L_7(u) = (u.t)(-1,0,1)   |
  |---------|--------------|--------------|--------------|-------------|----------------------------|
  |    8    |       1      |       8      |       0      |      1      |   L_8(u) = (u.t)(-1,-1,0)  |
  |---------|--------------|--------------|--------------|-------------|----------------------------|
  |    9    |       1      |       9      |       0      |      1      |   L_9(u) = (u.t)(1,-1,0)   |
  |---------|--------------|--------------|--------------|-------------|----------------------------|
  |   10    |       1      |      10      |       0      |      1      |   L_10(u) = (u.t)(1,1,0)   |
  |---------|--------------|--------------|--------------|-------------|----------------------------|
  |   11    |       1      |      11      |       0      |      1      |   L_11(u) = (u.t)(-1,1,0)  |
  |=========|==============|==============|==============|=============|============================|
  |   MAX   |  maxScDim=1  |  maxScOrd=11 |  maxDfOrd=0  |      -      |                            |
  |=========|==============|==============|==============|=============|============================|
  \endverbatim
  
    \remarks
    \li       The edge tangent \c t in the DoF definition is normalized by the \s edge length.
              As a result, the DoF functional is the value of the tangent component of a vector field 
              at the edge midpoint times the edge length. The resulting basis is equivalent to
              a basis defined by using the edge circulation as a DoF functional.
  
    \li       DefaultBasisFactory will select this class if the following parameters are specified:
  
  \verbatim
  |=======================|===================================|
  |  CellTopology         |  Hexahedron                       |
  |-----------------------|-----------------------------------|
  |  EFunctionSpace       |  FUNCTION_SPACE_HCURL             |
  |-----------------------|-----------------------------------|
  |  EDiscreteSpace       |  DISCRETE_SPACE_INCOMPLETE        |
  |-----------------------|-----------------------------------|
  |  degree               |  1                                |
  |-----------------------|-----------------------------------|
  |  EBasis               |  BASIS_FEM_DEFAULT                |
  |-----------------------|-----------------------------------|
  |  ECoordinates         |  COORDINATES_CARTESIAN            |
  |=======================|===================================|
\endverbatim
 */
  
template<class Scalar, class ArrayScalar> 
class Basis_HCURL_HEX_I1_FEM : public Basis<Scalar, ArrayScalar> {
private:
  
  /** \brief  Initializes <var>tagToOrdinal_</var> and <var>ordinalToTag_</var> lookup arrays.
   */
  void initializeTags();
  
public:

  /** \brief  Constructor.
    */
  Basis_HCURL_HEX_I1_FEM();
  
    
  /** \brief  Evaluation of a FEM basis on a <strong>reference Hexahedron</strong> cell. 
    
              Returns values of <var>operatorType</var> acting on FEM basis functions for a set of
              points in the <strong>reference Hexahedron</strong> cell. For rank and dimensions of
              I/O array arguments see Section \ref basis_md_array_sec.
  
      \param  outputValues      [out] - rank-3 array with the computed basis values
      \param  inputPoints       [in]  - rank-2 array with dimensions (P,D) containing reference points  
      \param  operatorType      [in]  - operator applied to basis functions    
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

#include "Intrepid_HCURL_HEX_I1_FEMDef.hpp"

#endif
