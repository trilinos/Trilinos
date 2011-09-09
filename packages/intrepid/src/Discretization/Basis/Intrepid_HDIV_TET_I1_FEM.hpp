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

/** \file   Intrepid_HDIV_TET_I1_FEM.hpp
    \brief  Header file for the Intrepid::HDIV_TET_I1_FEM class.
    \author Created by P. Bochev and D. Ridzal and K. Petrson.
 */

#ifndef INTREPID_HDIV_TET_I1_FEM_HPP
#define INTREPID_HDIV_TET_I1_FEM_HPP
#include "Intrepid_Basis.hpp"

namespace Intrepid {
  
/** \class  Intrepid::Basis_HDIV_TET_I1_FEM
    \brief  Implementation of the default H(div)-compatible FEM basis of degree 1 on Tetrahedron cell 
  
            Implements Raviart-Thomas basis of degree 1 on the reference Tetrahedron cell. The basis has
            cardinality 4 and spans an INCOMPLETE linear polynomial space. Basis functions are dual 
            to a unisolvent set of degrees-of-freedom (DoF) defined and enumerated as follows:
  
  \verbatim
  =========================================================================================================
  |         |           degree-of-freedom-tag table                    |                                  |
  |   DoF   |----------------------------------------------------------|       DoF definition             |
  | ordinal |  subc dim    | subc ordinal | subc DoF ord |subc num DoF |                                  |
  |=========|==============|==============|==============|=============|==================================|
  |    0    |       2      |       0      |       0      |      1      | L_0(u) = (u.n)(1/3,0,1/3)        |
  |---------|--------------|--------------|--------------|-------------|----------------------------------|
  |    1    |       2      |       1      |       0      |      1      | L_1(u) = (u.n)(1/3,1/3,1/3)      |
  |---------|--------------|--------------|--------------|-------------|----------------------------------|
  |    2    |       2      |       2      |       0      |      1      | L_2(u) = (u.n)(0,1/3,1/3)        |
  |---------|--------------|--------------|--------------|-------------|----------------------------------|
  |    3    |       2      |       3      |       0      |      1      | L_3(u) = (u.n)(1/3,1/3,0)        |
  |=========|==============|==============|==============|=============|==================================|
  |   MAX   |  maxScDim=2  |  maxScOrd=3  |  maxDfOrd=0  |      -      |                                  |
  |=========|==============|==============|==============|=============|==================================|
  \endverbatim
  
    \remarks
    \li     In the DoF functional \f${\bf n}\f$ is a face normal. Direction of face normals 
            is determined by the right-hand rule applied to faces oriented by their vertex order
            in the cell topology, from face vertex 0 to last face vertex, whereas their length is
            set equal to face area (see http://mathworld.wolfram.com/Right-HandRule.html for definition 
            of right-hand rule). For example, face 1 of all Tetrahedron cells has vertex order {1,2,3} 
            and its right-hand rule normal can be computed, e.g., by the vector product of edge 
            tangents to edges {1,2} and {2,3}. On the reference Tetrahedron the coordinates of 
            face 1 vertices are (1,0,0), (0,1,0), and (0,0,1), the edge tangents are (-1,1,0) and 
            (0,-1,1) and the face normal direction is (-1,1,0) X (0,-1,1) = (1,1,1). Length of this 
            raw face normal is twice the face area of face 1 and so the final face normal to face 1 is
            obtained by scaling the raw normal by 1/2: (1/2,1/2,1/2).
  
    \li     The length of the face normal equals the face area. As a result, the DoF functional 
            is the value of the normal component of a vector field at the face center times the 
            face area. The resulting basis is equivalent to a basis defined by using the face 
            flux as a DoF functional. Note that faces 0, 2, and 3 of reference Tetrahedron<> 
            cells have area 1/2 and face 1 has area Sqrt(3)/2.
  
 */
  
template<class Scalar, class ArrayScalar> 
class Basis_HDIV_TET_I1_FEM : public Basis<Scalar, ArrayScalar> {
private:
  
  /** \brief  Initializes <var>tagToOrdinal_</var> and <var>ordinalToTag_</var> lookup arrays.
   */
  void initializeTags();
  
public:

  /** \brief  Constructor.
    */
  Basis_HDIV_TET_I1_FEM();
  
    
  /** \brief  Evaluation of a FEM basis on a <strong>reference Tetrahedron</strong> cell. 
    
              Returns values of <var>operatorType</var> acting on FEM basis functions for a set of
              points in the <strong>reference Tetrahedron</strong> cell. For rank and dimensions of
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

#include "Intrepid_HDIV_TET_I1_FEMDef.hpp"

#endif
