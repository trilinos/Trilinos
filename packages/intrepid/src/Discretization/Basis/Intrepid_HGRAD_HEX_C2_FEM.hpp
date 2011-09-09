// @HEADER
// ************************************************************************
//
//                           Intrepid Package
//                 Copyright (2007) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Pavel Bochev  (pbboche@sandia.gov)
//                    Denis Ridzal  (dridzal@sandia.gov), or
//                    Kara Peterson (kjpeter@sandia.gov)
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
  
 */
  
template<class Scalar, class ArrayScalar> 
class Basis_HGRAD_HEX_C2_FEM : public Basis<Scalar, ArrayScalar>, public DofCoordsInterface<ArrayScalar> {
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

  /** \brief  Returns spatial locations (coordinates) of degrees of freedom on a
              <strong>reference Quadrilateral</strong>.

      \param  DofCoords      [out] - array with the coordinates of degrees of freedom,
                                     dimensioned (F,D)
  */
  void getDofCoords(ArrayScalar & DofCoords) const;
};
}// namespace Intrepid

#include "Intrepid_HGRAD_HEX_C2_FEMDef.hpp"

#endif
