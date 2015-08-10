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

/** \file   Intrepid_HDIV_TRI_I1_FEM.hpp
    \brief  Header file for the Intrepid::HDIV_TRI_I1_FEM class.
    \author Created by P. Bochev and D. Ridzal and K. Peterson.
 */

#ifndef INTREPID_HDIV_TRI_I1_FEM_HPP
#define INTREPID_HDIV_TRI_I1_FEM_HPP
#include "Intrepid_Basis.hpp"

namespace Intrepid {
  
/** \class  Intrepid::Basis_HDIV_TRI_I1_FEM
    \brief  Implementation of the default H(div)-compatible FEM basis of degree 1 on a Triangle cell.
  
            Implements Raviart-Thomas basis of degree 1 on the reference Triangle cell. The basis has
            cardinality 3 and spans an INCOMPLETE linear polynomial space. Basis functions are dual 
            to a unisolvent set of degrees-of-freedom (DoF) defined and enumerated as follows:
  
  \verbatim
  =========================================================================================================
  |         |           degree-of-freedom-tag table                    |                                  |
  |   DoF   |----------------------------------------------------------|       DoF definition             |
  | ordinal |  subc dim    | subc ordinal | subc DoF ord |subc num DoF |                                  |
  |=========|==============|==============|==============|=============|==================================|
  |    0    |       1      |       0      |       0      |      1      | L_0(u) = (u.n)(1/2,0)            |
  |---------|--------------|--------------|--------------|-------------|----------------------------------|
  |    1    |       1      |       1      |       0      |      1      | L_1(u) = (u.n)(1/2,1/2)          |
  |---------|--------------|--------------|--------------|-------------|----------------------------------|
  |    2    |       1      |       2      |       0      |      1      | L_2(u) = (u.n)(0,1/2)            |
  |=========|==============|==============|==============|=============|==================================|
  |   MAX   |  maxScDim=1  |  maxScOrd=2  |  maxDfOrd=0  |      -      |                                  |
  |=========|==============|==============|==============|=============|==================================|
  \endverbatim
  
    \remarks 
    \li     In the DOF functional \f${\bf n}=(t_2,-t_1)\f$ where \f${\bf t}=(t_1,t_2)\f$ 
            is the side (edge) tangent, i.e., the choice of normal direction is such that 
            the pair \f$({\bf n},{\bf t})\f$ is positively oriented. 
  
    \li     Direction of side tangents is determined by the vertex order of the sides in the
            cell topology and runs from side vertex 0 to side vertex 1, whereas their length is set 
            equal to the side length. For example, side 1 of all Triangle reference cells has vertex 
            order {1,2}, i.e., its tangent runs from vertex 1 of the reference Triangle to vertex 2 
            of that cell. On the reference Triangle the coordinates of these vertices are (1,0) and 
            (0,1), respectively. Therefore, the tangent to side 1 is (0,1)-(1,0) = (-1,1) and the 
            normal to that side is (1,1). Because its length already equals side length, no further 
            rescaling of the side tangent is needed.
        
    \li     The length of the side normal equals the length of the side. As a result, the 
            DoF functional is the value of the normal component of a vector field 
            at the side center times the side length. The resulting basis is equivalent to
            a basis defined by using the side flux as a DoF functional. Note that sides 0 and 2 of 
            reference Triangle<> cells have length 1 and side 1 has length Sqrt(2).
  
 */
  
template<class Scalar, class ArrayScalar> 
class Basis_HDIV_TRI_I1_FEM : public Basis<Scalar, ArrayScalar> {
private:
  
  /** \brief  Initializes <var>tagToOrdinal_</var> and <var>ordinalToTag_</var> lookup arrays.
   */
  void initializeTags();
  
public:

  /** \brief  Constructor.
    */
  Basis_HDIV_TRI_I1_FEM();
  
    
  /** \brief  Evaluation of a FEM basis on a <strong>reference Triangle</strong> cell. 
    
              Returns values of <var>operatorType</var> acting on FEM basis functions for a set of
              points in the <strong>reference Triangle</strong> cell. For rank and dimensions of
              I/O array arguments see Section \ref basis_md_array_sec.
  
      \param  outputValues      [out] - rank-3 or 4 array with the computed basis values
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

#include "Intrepid_HDIV_TRI_I1_FEMDef.hpp"

#endif
