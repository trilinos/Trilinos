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
