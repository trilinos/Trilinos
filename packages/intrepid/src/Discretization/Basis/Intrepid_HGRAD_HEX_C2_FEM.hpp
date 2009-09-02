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

/** \file   Intrepid_HGRAD_HEX_C2_FEM.hpp
    \brief  Header file for the Intrepid::HGRAD_HEX_C2_FEM class.
    \author Created by P. Bochev and D. Ridzal.
 */

#ifndef INTREPID_HGRAD_HEX_C2_FEM_HPP
#define INTREPID_HGRAD_HEX_C2_FEM_HPP
#include "Intrepid_Basis.hpp"

namespace Intrepid {
  
/** \class  Intrepid::Basis_HGRAD_HEX_C2_FEM
    \brief  Implementation of the default H(grad)-compatible FEM basis of degree 2 on Hexahedron cell 
  
            Implements Lagrangian basis of degree 2 on the reference Hexahedron cell. The basis has
            cardinality 27 and spans a COMPLETE tri-quadratic polynomial space. Basis functions are dual 
            to a unisolvent set of degrees-of-freedom (DoF) defined and enumerated as follows:
  
  \verbatim
  =================================================================================================
  |         |           degree-of-freedom-tag table                    |                           |
  |   DoF   |----------------------------------------------------------|      DoF definition       |
  | ordinal |  subc dim    | subc ordinal | subc DoF ord |subc num DoF |                           |
  |=========|==============|==============|==============|=============|===========================|
  |    0    |       0      |       0      |       0      |      1      |   L_0(u) = u(-1,-1,-1)    |
  |---------|--------------|--------------|--------------|-------------|---------------------------|
  |    1    |       0      |       1      |       0      |      1      |   L_1(u) = u( 1,-1,-1)    |
  |---------|--------------|--------------|--------------|-------------|---------------------------|
  |    2    |       0      |       2      |       0      |      1      |   L_2(u) = u( 1, 1,-1)    |
  |---------|--------------|--------------|--------------|-------------|---------------------------|
  |    3    |       0      |       3      |       0      |      1      |   L_3(u) = u(-1, 1,-1)    |
  |---------|--------------|--------------|--------------|-------------|---------------------------|
  |    4    |       0      |       4      |       0      |      1      |   L_4(u) = u(-1,-1, 1)    |
  |---------|--------------|--------------|--------------|-------------|---------------------------|
  |    5    |       0      |       5      |       0      |      1      |   L_5(u) = u( 1,-1, 1)    |
  |---------|--------------|--------------|--------------|-------------|---------------------------|
  |    6    |       0      |       6      |       0      |      1      |   L_6(u) = u( 1, 1, 1)    |
  |---------|--------------|--------------|--------------|-------------|---------------------------|
  |    7    |       0      |       7      |       0      |      1      |   L_7(u) = u(-1, 1, 1)    |
  |---------|--------------|--------------|--------------|-------------|---------------------------|
  |---------|--------------|--------------|--------------|-------------|---------------------------|
  |    8    |       1      |       0      |       0      |      1      |   L_8(u) = u( 0,-1,-1)    |
  |---------|--------------|--------------|--------------|-------------|---------------------------|
  |    9    |       1      |       1      |       0      |      1      |   L_9(u) = u( 1, 0,-1)    |
  |---------|--------------|--------------|--------------|-------------|---------------------------|
  |   10    |       1      |       2      |       0      |      1      |   L_10(u) = u( 0, 1,-1)   |
  |---------|--------------|--------------|--------------|-------------|---------------------------|
  |   11    |       1      |       3      |       0      |      1      |   L_11(u) = u(-1, 0,-1)   |
  |---------|--------------|--------------|--------------|-------------|---------------------------|
  |   12    |       1      |       8      |       0      |      1      |   L_12(u) = u(-1,-1, 0)   |
  |---------|--------------|--------------|--------------|-------------|---------------------------|
  |   13    |       1      |       9      |       0      |      1      |   L_13(u) = u( 1,-1, 0)   |
  |---------|--------------|--------------|--------------|-------------|---------------------------|
  |   14    |       1      |      10      |       0      |      1      |   L_14(u) = u( 1, 1, 0)   |
  |---------|--------------|--------------|--------------|-------------|---------------------------|
  |   15    |       1      |      11      |       0      |      1      |   L_15(u) = u(-1, 1, 0)   |
  |---------|--------------|--------------|--------------|-------------|---------------------------|
  |   16    |       1      |       4      |       0      |      1      |   L_16(u) = u( 0,-1, 1)   |
  |---------|--------------|--------------|--------------|-------------|---------------------------|
  |   17    |       1      |       5      |       0      |      1      |   L_17(u) = u( 1, 0, 1)   |
  |---------|--------------|--------------|--------------|-------------|---------------------------|
  |   18    |       1      |       6      |       0      |      1      |   L_18(u) = u( 0, 1, 1)   |
  |---------|--------------|--------------|--------------|-------------|---------------------------|
  |   19    |       1      |       7      |       0      |      1      |   L_19(u) = u(-1, 0, 1)   |
  |---------|--------------|--------------|--------------|-------------|---------------------------|
  |---------|--------------|--------------|--------------|-------------|---------------------------|
  |   20    |       3      |       0      |       0      |      1      |   L_20(u) = u( 0, 0, 0)   |
  |---------|--------------|--------------|--------------|-------------|---------------------------|
  |---------|--------------|--------------|--------------|-------------|---------------------------|
  |   21    |       2      |       4      |       0      |      1      |   L_21(u) = u( 0, 0,-1)   |
  |---------|--------------|--------------|--------------|-------------|---------------------------|
  |   22    |       2      |       5      |       0      |      1      |   L_22(u) = u( 0, 0, 1)   |
  |---------|--------------|--------------|--------------|-------------|---------------------------|
  |   23    |       2      |       3      |       0      |      1      |   L_23(u) = u(-1, 0, 0)   |
  |---------|--------------|--------------|--------------|-------------|---------------------------|
  |   24    |       2      |       1      |       0      |      1      |   L_24(u) = u( 1, 0, 0)   |
  |---------|--------------|--------------|--------------|-------------|---------------------------|
  |   25    |       2      |       0      |       0      |      1      |   L_25(u) = u( 0,-1, 0)   |
  |---------|--------------|--------------|--------------|-------------|---------------------------|
  |   26    |       2      |       2      |       0      |      1      |   L_26(u) = u( 0, 1, 0)   |
  |=========|==============|==============|==============|=============|===========================|
  |   MAX   |  maxScDim=2  |  maxScOrd=12 |  maxDfOrd=0  |      -      |                           |
  |=========|==============|==============|==============|=============|===========================|
  \endverbatim
  
    \remark   Ordering of DoFs follows the node order in Hexahedron<27> topology. Note that node
              order in this topology does not follow the natural oder of k-subcells where the nodes
              are located, except for nodes 0 to 7 which coincide with the vertices of the base
              Hexahedrn <8> topology. As a result, L_0 to L_7 are associated with nodes 0 to 7, but
              L_8 to L_19 are not associated with edges 0 to 12 in that order.
  
    \remarks
    \li       DefaultBasisFactory will select this class if the following parameters are specified:
  
  \verbatim
  |=======================|===================================|
  |  CellTopology         |  Hexahedron                       |
  |-----------------------|-----------------------------------|
  |  EFunctionSpace       |  FUNCTION_SPACE_HGRAD             |
  |-----------------------|-----------------------------------|
  |  EDiscreteSpace       |  DISCRETE_SPACE_COMPLETE          |
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
class Basis_HGRAD_HEX_C2_FEM : public Basis<Scalar, ArrayScalar> {
private:
  
  /** \brief  Initializes <var>tagToOrdinal_</var> and <var>ordinalToTag_</var> lookup arrays.
   */
  void initializeTags();
  
public:

  /** \brief  Constructor.
    */
  Basis_HGRAD_HEX_C2_FEM();
  
    
  /** \brief  Evaluation of a FEM basis on a <strong>reference Hexahedron</strong> cell. 
    
              Returns values of <var>operatorType</var> acting on FEM basis functions for a set of
              points in the <strong>reference Hexahedron</strong> cell. For rank and dimensions of
              I/O array arguments see Section \ref basis_md_array_sec.
  
      \param  outputValues      [out] - rank-2 or 3 array with the computed basis values
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

#include "Intrepid_HGRAD_HEX_C2_FEMDef.hpp"

#endif
