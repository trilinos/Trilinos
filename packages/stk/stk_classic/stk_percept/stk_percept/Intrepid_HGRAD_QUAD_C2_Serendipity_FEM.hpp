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

/** \file   Intrepid_G_QUAD_C2_Serendipity_FEM.hpp
    \brief  Header file for the Intrepid::G_QUAD_C2_Serendipity_FEM class.
    \author Created by P. Bochev and D. Ridzal.
 */

#ifndef INTREPID_HGRAD_QUAD_C2_Serendipity_FEM_HPP
#define INTREPID_HGRAD_QUAD_C2_Serendipity_FEM_HPP
#include "Intrepid_Basis.hpp"

namespace Intrepid {
  
/** \class  Intrepid::Basis_HGRAD_QUAD_C2_Serendipity_FEM
    \brief  Implementation of the default H(grad)-compatible FEM serendipity basis of incomplete degree 2 on Quadrilateral cell
  
            Implements Lagrangian serendipity basis of degree 2 on the reference Quadrilateral cell. The basis has
            cardinality 8 and spans an INCOMPLETE bi-quadratic polynomial space. Basis functions are dual 
            to a unisolvent set of degrees-of-freedom (DoF) defined and enumerated as follows:
  
  \verbatim
  =================================================================================================
  |         |           degree-of-freedom-tag table                    |                           |
  |   DoF   |----------------------------------------------------------|      DoF definition       |
  | ordinal |  subc dim    | subc ordinal | subc DoF ord |subc num DoF |                           |
  |=========|==============|==============|==============|=============|===========================|
  |    0    |       0      |       0      |       0      |      1      |   L_0(u) = u(-1,-1)       |
  |---------|--------------|--------------|--------------|-------------|---------------------------|
  |    1    |       0      |       1      |       0      |      1      |   L_1(u) = u( 1,-1)       |
  |---------|--------------|--------------|--------------|-------------|---------------------------|
  |    2    |       0      |       2      |       0      |      1      |   L_2(u) = u( 1, 1)       |
  |---------|--------------|--------------|--------------|-------------|---------------------------|
  |    3    |       0      |       3      |       0      |      1      |   L_3(u) = u(-1, 1)       |
  |---------|--------------|--------------|--------------|-------------|---------------------------|
  |    4    |       1      |       0      |       0      |      1      |   L_4(u) = u( 0,-1)       |
  |---------|--------------|--------------|--------------|-------------|---------------------------|
  |    5    |       1      |       1      |       0      |      1      |   L_5(u) = u( 1, 0)       |
  |---------|--------------|--------------|--------------|-------------|---------------------------|
  |    6    |       1      |       2      |       0      |      1      |   L_6(u) = u( 0, 1)       |
  |---------|--------------|--------------|--------------|-------------|---------------------------|
  |    7    |       1      |       3      |       0      |      1      |   L_7(u) = u(-1, 0)       |
  |=========|==============|==============|==============|=============|===========================|
  |   MAX   |  maxScDim=2  |  maxScOrd=3  |  maxDfOrd=0  |     -       |                           |
  |=========|==============|==============|==============|=============|===========================|
  \endverbatim
*/
template<class Scalar, class ArrayScalar> 
class Basis_HGRAD_QUAD_C2_Serendipity_FEM : public Basis<Scalar, ArrayScalar>, public DofCoordsInterface<ArrayScalar>  {
private:

  /** \brief Initializes <var>tagToOrdinal_</var> and <var>ordinalToTag_</var> lookup arrays.
   */
  void initializeTags();
  
public:

  /** \brief Constructor.
  */
  Basis_HGRAD_QUAD_C2_Serendipity_FEM();
  
    
  /** \brief  FEM basis evaluation on a <strong>reference Quadrilateral</strong> cell. 
    
              Returns values of <var>operatorType</var> acting on FEM basis functions for a set of
              points in the <strong>reference Quadrilateral</strong> cell. For rank and dimensions of 
              I/O array arguments see Section \ref basis_md_array_sec .
  
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

#include "Intrepid_HGRAD_QUAD_C2_Serendipity_FEMDef.hpp"

#endif
