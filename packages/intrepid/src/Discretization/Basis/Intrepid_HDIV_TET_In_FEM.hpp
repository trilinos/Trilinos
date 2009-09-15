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

/** \file   Intrepid_HDIV_TET_In_FEM.hpp
    \brief  Header file for the Intrepid::HDIV_TET_In_FEM class.
    \author Created by R. Kirby and P. Bochev and D. Ridzal.
*/

#ifndef INTREPID_HDIV_TET_In_FEM_HPP
#define INTREPID_HDIV_TET_In_FEM_HPP

#include "Intrepid_Types.hpp"
#include "Intrepid_Basis.hpp"
#include "Intrepid_PointTools.hpp"
#include "Intrepid_HGRAD_TET_Cn_FEM_ORTH.hpp"
#include "Intrepid_CubatureDirectTetDefault.hpp"
#include "Teuchos_SerialDenseMatrix.hpp"
#include "Teuchos_SerialDenseSolver.hpp"


namespace Intrepid {
  
/** \class  Intrepid::Basis_HDIV_TET_In_FEM
    \brief  Implementation of the default H(div)-compatible
            Raviart-Thomas basis of arbitrary degree  on Tetrahedron
            cell. The lowest order instance starts with n.  Implements
            the nodal basis of degree n the reference Tetrahedron cell. The basis has
            cardinality n(n+1)(n+3)/2 and spans an INCOMPLETE
            polynomial space of degree n. Basis functions are dual 
            to a unisolvent set of degrees-of-freedom (DoF) defined and enumerated as follows:
  
            \li The normal component on a lattice of order n+1 and
            offset 1 on each face (see PointTools). This gives one point per edge in
            the lowest-order case.  These are the first
            4 * (n*(n+1)/2) degrees of freedom.

            \li If n > 1, the x and y z components at a lattice of
            order n+2 and offset on the interior of the tetrahedron.  These are the rest
            of the degrees of freedom.  

            If the pointType argument to the constructor specifies equispaced points, then the face
            and interior points will be equispaced.  If
            the pointType argument specifies warp-blend points, the
            interior of a warp-blend lattice will be used on each face
            and also for the cell interior.


  */
  
template<class Scalar, class ArrayScalar> 
class Basis_HDIV_TET_In_FEM: public Basis<Scalar, ArrayScalar> {
private:
  
  /** \brief  Initializes <var>tagToOrdinal_</var> and <var>ordinalToTag_</var> lookup arrays.
   */
  virtual void initializeTags();

  /** \brief Orthogonal basis out of which the nodal basis is
      constructed */
  Basis_HGRAD_TET_Cn_FEM_ORTH<Scalar,FieldContainer<Scalar> > Phis_;

  /** \brief expansion coefficients of the nodal basis in terms of the
      orthgonal one */
  FieldContainer<Scalar> coeffs_;


public:
  
  /** \brief  Constructor.
   */
  Basis_HDIV_TET_In_FEM( const int n , const EPointType pointType );    
  
  /** \brief  Evaluation of a FEM basis on a <strong>reference Tetrahedron</strong> cell. 
  
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
};

}// namespace Intrepid

#include "Intrepid_HDIV_TET_In_FEMDef.hpp"

#endif
