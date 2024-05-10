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

/** \file   Intrepid_HGRAD_TRI_Cn_FEM.hpp
    \brief  Header file for the Intrepid::HGRAD_TRI_Cn_FEM class.
    \author Created by R. Kirby and P. Bochev and D. Ridzal.
*/

#ifndef INTREPID_HGRAD_TRI_Cn_FEM_HPP
#define INTREPID_HGRAD_TRI_Cn_FEM_HPP

#include "Intrepid_Types.hpp"
#include "Intrepid_Basis.hpp"
#include "Intrepid_PointTools.hpp"
#include "Intrepid_HGRAD_TRI_Cn_FEM_ORTH.hpp"
#include "Teuchos_SerialDenseMatrix.hpp"
#include "Teuchos_SerialDenseSolver.hpp"


namespace Intrepid {
  
/** \class  Intrepid::Basis_HGRAD_TRI_Cn_FEM
    \brief  Implementation of the default H(grad)-compatible Lagrange basis of arbitrary degree  on Triangle cell 
  
            Implements Lagrangian basis of degree n on the reference Triangle cell. The basis has
            cardinality (n+1)(n+2)/2 and spans a COMPLETE polynomial space of degree n.
            Basis functions are dual 
            to a unisolvent set of degrees-of-freedom (DoF) defined on a lattice of order n (see PointTools).
            In particular, the degrees of freedom are point evaluation at
            \li the vertices
            \li (n-1) points on each edge of the triangle
            \li max((n-1)(n-2)/2,0) points on the inside of the triangle.
  
            The distribution of these points is specified by the pointType argument to the class constructor.
            Currently, either equispaced lattice points or Warburton's warp-blend points are available.

            The dof are enumerated according to the ordering on the lattice (see PointTools).  In particular,
            dof number 0 is at the bottom left vertex (0,0).  The dof increase
            along the lattice with points along the lines of constant
            x adjacent in the enumeration. 
  */
  
template<class Scalar, class ArrayScalar> 
class Basis_HGRAD_TRI_Cn_FEM: public Basis<Scalar, ArrayScalar> {
private:
  
  /** \brief  Initializes <var>tagToOrdinal_</var> and <var>ordinalToTag_</var> lookup arrays.
   */
  virtual void initializeTags();

  /** \brief  The orthogonal basis on triangles, out of which the nodal basis is constructed
   */
  Basis_HGRAD_TRI_Cn_FEM_ORTH<Scalar,FieldContainer<Scalar> > Phis;
  /** \brief  The Vandermonde matrix with V_{ij} = phi_i(x_j), where x_j is the j_th point in the lattice
   */
  FieldContainer<Scalar> V;
  /** \brief  The inverse of V.  The columns of Vinv express the Lagrange basis in terms of the
      orthogonal basis
   */
  FieldContainer<Scalar> Vinv;
  /** \brief stores the points at which degrees of freedom are located. 
   */
  FieldContainer<Scalar> latticePts;

public:
  
  /** \brief  Constructor.
   */
  Basis_HGRAD_TRI_Cn_FEM( const int n , const EPointType pointType );    
  
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
};

}// namespace Intrepid

#include "Intrepid_HGRAD_TRI_Cn_FEMDef.hpp"

#endif

#if defined(Intrepid_SHOW_DEPRECATED_WARNINGS)
#ifdef __GNUC__
#warning "The Intrepid package is deprecated"
#endif
#endif

