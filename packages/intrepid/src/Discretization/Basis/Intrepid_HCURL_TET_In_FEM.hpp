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

/** \file   Intrepid_HCURL_TET_In_FEM.hpp
    \brief  Header file for the Intrepid::HCURL_TET_In_FEM class.
    \author Created by R. Kirby 
*/

#ifndef INTREPID_HCURL_TET_In_FEM_HPP
#define INTREPID_HCURL_TET_In_FEM_HPP

#include "Intrepid_Types.hpp"
#include "Intrepid_Basis.hpp"
#include "Intrepid_PointTools.hpp"
#include "Intrepid_HGRAD_TET_Cn_FEM_ORTH.hpp"
#include "Intrepid_CubatureDirectTetDefault.hpp"
#include "Intrepid_CubatureDirectLineGauss.hpp"
#include "Teuchos_SerialDenseMatrix.hpp"
#include "Teuchos_SerialDenseSolver.hpp"


namespace Intrepid {
  
/** \class  Intrepid::Basis_HCURL_TET_In_FEM
    \brief  Implementation of the default H(curl)-compatible Nedelec (first kind) 
            basis of arbitrary degree  on Tetrahedron cell.  The lowest order space
            is indexted with 1 rather than 0.
            Implements nodal basis of degree n (n>=1) on the reference Tetrahedron cell. The basis has
            cardinality n*(n+2)*(n+3)/2 and spans an INCOMPLETE
            polynomial space of degree n. Basis functions are dual 
            to a unisolvent set of degrees-of-freedom (DoF) defined by

            \li The tangential component of the vector field at n
            points per edge.  These are the first 6 * n dof.
            \li The components in the directions of the tangents to
            each face (see CellTools) on a lattice of order n+1 with
            offset 1 (see PointTools) on each face.  These are the
            next 4 * n*(n-1) dof.
            \li The x,y,z components on a lattice of order n+1 with
            offset 1 on the interior of the tetrahedron, the remaining
            degrees of freedom.
            
            If the pointType argument to the constructor specifies equispaced points, then the edge points
            will be equispaced on each edge and the interior points equispaced also.  If
            the pointType argument specifies warp-blend points, then Gauss-Lobatto points of order n
            are chosen on each edge and the interior of warp-blend
            lattice is chosen on each face and the interior.
  
  \verbatim
  \endverbatim
  
  */
  
template<class Scalar, class ArrayScalar> 
class Basis_HCURL_TET_In_FEM: public Basis<Scalar, ArrayScalar> {
private:
  
  /** \brief  Initializes <var>tagToOrdinal_</var> and <var>ordinalToTag_</var> lookup arrays.
   */
  virtual void initializeTags();


  /** \brief Orthogonal basis of ofder n, in terms of which the H(curl) basis functions are expressed */
  Basis_HGRAD_TET_Cn_FEM_ORTH<Scalar,FieldContainer<Scalar> > Phis_;
  /** \brief Array holding the expansion coefficients of the nodal basis in terms of Phis_ */
  FieldContainer<Scalar> coeffs_;


public:
  
  /** \brief  Constructor.
   */
  Basis_HCURL_TET_In_FEM( const int n , const EPointType pointType );    
  
  /** \brief  Evaluation of a FEM basis on a <strong>reference Tetrahedron</strong> cell. 
  
              Returns values of <var>operatorType</var> acting on FEM basis functions for a set of
              points in the <strong>reference Tetrahedron</strong> cell. For rank and dimensions of
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

#include "Intrepid_HCURL_TET_In_FEMDef.hpp"

#endif

#if defined(Intrepid_SHOW_DEPRECATED_WARNINGS)
#ifdef __GNUC__
#warning "The Intrepid package is deprecated"
#endif
#endif

