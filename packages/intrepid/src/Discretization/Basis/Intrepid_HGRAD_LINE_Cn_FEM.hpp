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

/** \file   Intrepid_HGRAD_LINE_Cn_FEM_JACOBI.hpp
    \brief  Header file for the Intrepid::HGRAD_LINE_Cn_FEM class.
    \author Created by R. Kirby and P. Bochev and D. Ridzal.
*/

#ifndef INTREPID_HGRAD_LINE_Cn_FEM_HPP
#define INTREPID_HGRAD_LINE_Cn_FEM_HPP

#include "Intrepid_Basis.hpp"
#include "Intrepid_PointTools.hpp"
#include "Intrepid_Polylib.hpp"
#include "Intrepid_HGRAD_LINE_Cn_FEM_JACOBI.hpp"
#include "Teuchos_SerialDenseMatrix.hpp"
#include "Teuchos_SerialDenseSolver.hpp"


namespace Intrepid {
  
/** \class  Intrepid::Basis_HGRAD_LINE_Cn_FEM
    \brief  Implementation of the locally H(grad)-compatible FEM basis of variable order
            on the [-1,1] reference line cell, using Lagrange polynomials. 
  
            Implements Lagrange basis of variable order \f$n\f$ on
            the reference [-1,1] line cell.  The distribution of the points
	    may be equispaced points with our without the endpoints, the Gauss-Legendre
	    or Gauss-Lobatto points.  These points are {x_i}_{i=0}^n. with x_i < x_{i+1}

            The basis has cardinality \f$n+1\f$ and spans a COMPLETE linear polynomial space.
            Basis functions are dual to a unisolvent set of degrees of freedom (DoF)
	    n_i( psi ) = psi(x_i).  The DoF are ordered by i.  The DoF at points 
	    -1 and 1 (if included in {x_i} are attached to the vertices, and the rest of
	    the DoF are attached to the edge itself.
*/
template<class Scalar, class ArrayScalar> 
class Basis_HGRAD_LINE_Cn_FEM: 
    public Basis<Scalar, ArrayScalar>, 
    public DofCoordsInterface<ArrayScalar> {
private:
  /** \brief Holds the points defining the Lagrange basis */
  FieldContainer<Scalar> latticePts_;

  /** \brief orthogonal basis */
  Basis_HGRAD_LINE_Cn_FEM_JACOBI<Scalar, FieldContainer<Scalar> > Phis_;
  
  /** \brief Generalized Vandermonde matrix V_{ij} = phis_i(x_j) */
  FieldContainer<Scalar> V_;

  /** \brief inverse of Generalized Vandermonde matrix, whose columns store the expansion
      coefficients of the nodal basis in terms of phis_ */
  FieldContainer<Scalar> Vinv_;

  /** \brief  Initializes <var>tagToOrdinal_</var> and <var>ordinalToTag_</var> lookup arrays.
  */
  void initializeTags();
  
public:
  /** \brief Destructor
   */
  ~Basis_HGRAD_LINE_Cn_FEM( ) { }

  /** \brief  Constructor.
  */
  Basis_HGRAD_LINE_Cn_FEM(int order , const ArrayScalar &pts );  

  /** \brief Constructor.
      \param  int order:        [in] polynomial degree of the basis
      \param  int pointType:    [in] type of points, either POINTTYPE_EQUISPACED or POINTTYPE_SPECTRAL */
  Basis_HGRAD_LINE_Cn_FEM(int order , const EPointType &pointType );  
  
  
  /** \brief  Evaluation of a FEM basis on a <strong>reference Line</strong> cell. 
  
              Returns values of <var>operatorType</var> acting on FEM basis functions for a set of
              points in the <strong>reference Line</strong> cell. For rank and dimensions of
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

  /** \brief implements the dofcoords interface */
  virtual void getDofCoords( ArrayScalar & DofCoords ) const;

};

}// namespace Intrepid

#include "Intrepid_HGRAD_LINE_Cn_FEMDef.hpp"

#endif
