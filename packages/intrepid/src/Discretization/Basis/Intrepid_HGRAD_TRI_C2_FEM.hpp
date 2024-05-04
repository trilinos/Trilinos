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

/** \file   Intrepid_HGRAD_TRI_C2_FEM.hpp
    \brief  Header file for the Intrepid::HGRAD_TRI_C2_FEM class.
    \author Created by P. Bochev and D. Ridzal.
*/

#ifndef INTREPID_HGRAD_TRI_C2_FEM_HPP
#define INTREPID_HGRAD_TRI_C2_FEM_HPP

#include "Intrepid_Basis.hpp"

namespace Intrepid {
  
/** \class  Intrepid::Basis_HGRAD_TRI_C2_FEM
    \brief  Implementation of the default H(grad)-compatible FEM basis of degree 2 on Triangle cell 
  
            Implements Lagrangian basis of degree 2 on the reference Triangle cell. The basis has
            cardinality 6 and spans a COMPLETE quadratic polynomial space. Basis functions are dual 
            to a unisolvent set of degrees-of-freedom (DoF) defined and enumerated as follows:
  
  \verbatim
  =================================================================================================
  |         |           degree-of-freedom-tag table                    |                           |
  |   DoF   |----------------------------------------------------------|      DoF definition       |
  | ordinal |  subc dim    | subc ordinal | subc DoF ord |subc num DoF |                           |
  |=========|==============|==============|==============|=============|===========================|
  |    0    |       0      |       0      |       0      |      1      |   L_0(u) = u(0,0)         |
  |---------|--------------|--------------|--------------|-------------|---------------------------|
  |    1    |       0      |       1      |       0      |      1      |   L_1(u) = u(1,0)         |
  |---------|--------------|--------------|--------------|-------------|---------------------------|
  |    2    |       0      |       2      |       0      |      1      |   L_2(u) = u(0,1)         |
  |---------|--------------|--------------|--------------|-------------|---------------------------|
  |    3    |       1      |       0      |       0      |      1      |   L_3(u) = u(1/2,0)       |
  |---------|--------------|--------------|--------------|-------------|---------------------------|
  |    4    |       1      |       1      |       0      |      1      |   L_4(u) = u(1/2,1/2)     |
  |---------|--------------|--------------|--------------|-------------|---------------------------|
  |    5    |       1      |       2      |       0      |      1      |   L_5(u) = u(0,1/2)       |
  |=========|==============|==============|==============|=============|===========================|
  |   MAX   |  maxScDim=1  |  maxScOrd=2  |  maxDfOrd=0  |     -       |                           |
  |=========|==============|==============|==============|=============|===========================|
  \endverbatim
  */
  
template<class Scalar, class ArrayScalar> 
class Basis_HGRAD_TRI_C2_FEM: public Basis<Scalar, ArrayScalar>, public DofCoordsInterface<ArrayScalar> {
private:
  
  /** \brief  Initializes <var>tagToOrdinal_</var> and <var>ordinalToTag_</var> lookup arrays.
   */
  void initializeTags();
  
public:
  
  /** \brief  Constructor.
   */
  Basis_HGRAD_TRI_C2_FEM();  
  
  
  /** \brief  Evaluation of a FEM basis on a <strong>reference Triangle</strong> cell. 
  
              Returns values of <var>operatorType</var> acting on FEM basis functions for a set of
              points in the <strong>reference Triangle</strong> cell. For rank and dimensions of
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

  /** \brief  Returns spatial locations (coordinates) of degrees of freedom on a
              <strong>reference Quadrilateral</strong>.

      \param  DofCoords      [out] - array with the coordinates of degrees of freedom,
                                     dimensioned (F,D)
  */
  void getDofCoords(ArrayScalar & DofCoords) const;

};

}// namespace Intrepid

#include "Intrepid_HGRAD_TRI_C2_FEMDef.hpp"

#endif

#if defined(Intrepid_SHOW_DEPRECATED_WARNINGS)
#ifdef __GNUC__
#warning "The Intrepid package is deprecated"
#endif
#endif

