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
//                    Denis Ridzal (dridzal@sandia.gov)
//                    Kara Peterson (kjpeter@sandia.gov).
//
// ************************************************************************
// @HEADER

/** \file   Intrepid_HDIV_QUAD_I1_FEM.hpp
    \brief  Header file for the Intrepid::HDIV_QUAD_I1_FEM class.
    \author Created by P. Bochev and D. Ridzal and K. Petrson.
 */

#ifndef INTREPID_HDIV_QUAD_I1_FEM_HPP
#define INTREPID_HDIV_QUAD_I1_FEM_HPP
#include "Intrepid_Basis.hpp"

namespace Intrepid {
  
/** \class  Intrepid::Basis_HDIV_QUAD_I1_FEM
    \brief  Implementation of the default H(div)-compatible FEM basis of degree 1 on Quadrilateral cell 
  
            Implements Raviart-Thomas basis of degree 1 on the reference Quadrilateral cell. The basis has
            cardinality 4 and spans a INCOMPLETE bi-linear polynomial space. Basis functions are dual 
            to a unisolvent set of degrees-of-freedom (DoF) defined and enumerated as follows:
  
  \verbatim
  ===================================================================================================
  |         |           degree-of-freedom-tag table                    |                            |
  |   DoF   |----------------------------------------------------------|       DoF definition       |
  | ordinal |  subc dim    | subc ordinal | subc DoF ord |subc num DoF |                            |
  |=========|==============|==============|==============|=============|============================|
  |    0    |       1      |       0      |       0      |      1      |   L_0(u) = (u.n)(0,-1)     |
  |---------|--------------|--------------|--------------|-------------|----------------------------|
  |    1    |       1      |       1      |       0      |      1      |   L_1(u) = (u.n)(1,0)      |
  |---------|--------------|--------------|--------------|-------------|----------------------------|
  |    2    |       1      |       2      |       0      |      1      |   L_2(u) = (u.n)(0,1)      |
  |---------|--------------|--------------|--------------|-------------|----------------------------|
  |    3    |       1      |       3      |       0      |      1      |   L_3(u) = (u.n)(-1,0)     |
  |=========|==============|==============|==============|=============|============================|
  |   MAX   |  maxScDim=2  |  maxScOrd=5  |  maxDfOrd=0  |      -      |                            |
  |=========|==============|==============|==============|=============|============================|
  \endverbatim
  
    \remarks
    \li     In the DOF functional \f${\bf n}=(t_2,-t_1)\f$ where \f${\bf t}=(t_1,t_2)\f$ 
            is the side (edge) tangent, i.e., the choice of normal direction is such that 
            the pair \f$({\bf n},{\bf t})\f$ is positively oriented. 
  
    \li     Direction of side tangents is determined by the vertex order of the sides in the
            cell topology and runs from side vertex 0 to side vertex 1, whereas their length is set 
            equal to the side length. For example, side 1 of all Quadrilateral reference cells has 
            vertex order {1,2}, i.e., its tangent runs from vertex 1 of the reference Quadrilateral 
            to vertex 2 of that cell. On the reference Quadrilateral the coordinates of these vertices 
            are (1,-1) and (1,1), respectively. Therefore, the tangent to side 1 is (1,1)-(1,-1) = (0,2) 
            and the normal to that side is (2,0). Because its length already equals side length, no 
            further rescaling of the side tangent is needed.
  
    \li     The length of the side normal equals the length of the side. As a result, the 
            DoF functional is the value of the normal component of a vector field 
            at the side center times the side length. The resulting basis is equivalent to
            a basis defined by using the side flux as a DoF functional. Note that all sides of 
            the reference Quadrilateral<> cells have length 2.
  
 */
  
template<class Scalar, class ArrayScalar> 
class Basis_HDIV_QUAD_I1_FEM : public Basis<Scalar, ArrayScalar> {
private:
  
  /** \brief  Initializes <var>tagToOrdinal_</var> and <var>ordinalToTag_</var> lookup arrays.
   */
  void initializeTags();
  
public:

  /** \brief  Constructor.
    */
  Basis_HDIV_QUAD_I1_FEM();
  
    
  /** \brief  Evaluation of a FEM basis on a <strong>reference Quadrilateral</strong> cell. 
    
              Returns values of <var>operatorType</var> acting on FEM basis functions for a set of
              points in the <strong>reference Quadrilateral</strong> cell. For rank and dimensions of
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

#include "Intrepid_HDIV_QUAD_I1_FEMDef.hpp"

#endif
