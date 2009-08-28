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

/** \file   Intrepid_HDIV_WEDGE_I1_FEM.hpp
    \brief  Header file for the Intrepid::HDIV_WEDGE_I1_FEM class.
    \author Created by P. Bochev and D. Ridzal.
 */

#ifndef INTREPID_HDIV_WEDGE_I1_FEM_HPP
#define INTREPID_HDIV_WEDGE_I1_FEM_HPP
#include "Intrepid_Basis.hpp"

namespace Intrepid {
  
/** \class  Intrepid::Basis_HDIV_WEDGE_I1_FEM
    \brief  Implementation of the default H(div)-compatible FEM basis of degree 1 on Wedge cell 
  
            Implements Raviart-Thomas basis of degree 1 on the reference Wedge cell. The basis has
            cardinality 5 and spans an INCOMPLETE bi-linear polynomial space. Basis functions are dual 
            to a unisolvent set of degrees-of-freedom (DoF) defined and enumerated as follows:
  
  \verbatim
  ===================================================================================================
  |         |           degree-of-freedom-tag table                    |                            |
  |   DoF   |----------------------------------------------------------|       DoF definition       |
  | ordinal |  subc dim    | subc ordinal | subc DoF ord |subc num DoF |                            |
  |=========|==============|==============|==============|=============|============================|
  |    0    |       2      |       0      |       0      |      1      |  L_0(u) = (u.n)(1/2,0,0)   |
  |---------|--------------|--------------|--------------|-------------|----------------------------|
  |    1    |       2      |       1      |       0      |      1      |  L_1(u) = (u.n)(1/2,1/2,0) |
  |---------|--------------|--------------|--------------|-------------|----------------------------|
  |    2    |       2      |       2      |       0      |      1      |  L_2(u) = (u.n)(0,1/2,0)   |
  |---------|--------------|--------------|--------------|-------------|----------------------------|
  |    3    |       2      |       3      |       0      |      1      |  L_3(u) = (u.n)(1/3,1/3,-1)|
  |---------|--------------|--------------|--------------|-------------|----------------------------|
  |    4    |       2      |       4      |       0      |      1      |  L_4(u) = (u.n)(1/3,1/3,1) |
  |=========|==============|==============|==============|=============|============================|
  |   MAX   |  maxScDim=2  |  maxScOrd=4  |  maxDfOrd=0  |      -      |                            |
  |=========|==============|==============|==============|=============|============================|
  \endverbatim
  
    \remarks
    \li       The face outer normal \c n in the DoF definition is normalized by the \s face area.
              As a result, the DoF functional is the value of the normal component of a vector field 
              at the face center times the face area. The so defined basis is equivalent to a basis 
              defined by using the face flux as a DoF functional. Note that faces 0 and 2 of reference 
              Wedge<> cells have area 2, face 1 has area 2*Sqrt(2), and faces 3 and 4 have area 1/2.
  
    \li       DefaultBasisFactory will select this class if the following parameters are specified:
  
  \verbatim
  |=======================|===================================|
  |  CellTopology         |  Wedge                            |
  |-----------------------|-----------------------------------|
  |  EFunctionSpace       |  FUNCTION_SPACE_HDIV              |
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
class Basis_HDIV_WEDGE_I1_FEM : public Basis<Scalar, ArrayScalar> {
private:
  
  /** \brief  Initializes <var>tagToOrdinal_</var> and <var>ordinalToTag_</var> lookup arrays.
   */
  void initializeTags();
  
public:

  /** \brief  Constructor.
    */
  Basis_HDIV_WEDGE_I1_FEM();
  
    
  /** \brief  Evaluation of a FEM basis on a <strong>reference Wedge</strong> cell. 
    
              Returns values of <var>operatorType</var> acting on FEM basis functions for a set of
              points in the <strong>reference Wedge</strong> cell. For rank and dimensions of
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

#include "Intrepid_HDIV_WEDGE_I1_FEMDef.hpp"

#endif
