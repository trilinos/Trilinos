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

/** \file   Intrepid_HGRAD_TET_Cn_FEM.hpp
    \brief  Header file for the Intrepid::HGRAD_TET_Cn_FEM class.
    \author Created by R. Kirby and P. Bochev and D. Ridzal.
*/

#ifndef INTREPID_HGRAD_TET_Cn_FEM_HPP
#define INTREPID_HGRAD_TET_Cn_FEM_HPP

#include "Intrepid_Types.hpp"
#include "Intrepid_Basis.hpp"
#include "Intrepid_PointTools.hpp"
#include "Intrepid_HGRAD_TET_Cn_FEM_ORTH.hpp"
#include "Teuchos_SerialDenseMatrix.hpp"
#include "Teuchos_SerialDenseSolver.hpp"


namespace Intrepid {
  
/** \class  Intrepid::Basis_HGRAD_TET_Cn_FEM
    \brief  Implementation of the default H(grad)-compatible Lagrange basis of arbitrary degree  on Tetrahedron cell 
  
            Implements Lagrangian basis of degree n on the reference Tetrahedron cell. The basis has
            cardinality (n+1)(n+2)(n+3)/6 and spans a COMPLETE polynomial space of degree n. 
	    Nodal basis functions are dual to a unisolvent set of
            degrees-of-freedom (DoF) defined at a lattice of order n
            (see \ref PointTools).  In particular, the degrees of freedom
            are point evaluation at
	    \li The vertices
	    \li (n-1) points on each edge of the tetrahedron
	    \li max((n-1)(n-2)/2,0) points on each face of the
            tetrahedron
	    \li max((n-1)(n-2)(n-3)/6,0) points in the interior
	    of the tetrahedron.

	    The distribution of these points is specified by the pointType argument to the class constructor.
	    Currently, either equispaced lattice points or Warburton's warp-blend points are available.

	    The dof are enumerated according to the ordering on the lattice (see PointTools).  In particular,
	    dof number 0 is at the vertex (0,0,0).  The dof increase
	    along the lattice with points along the lines of constant
	    x adjacent in the enumeration. 


\remarks
  DefaultBasisFactory will select this class if the following parameters are specified:
  
  \verbatim
  |=======================|===================================|
  |  CellTopology         |  Tetrahedron                      |
  |-----------------------|-----------------------------------|
  |  EFunctionSpace       |  FUNCTION_SPACE_HGRAD             |
  |-----------------------|-----------------------------------|
  |  EDiscreteSpace       |  DISCRETE_SPACE_COMPLETE          |
  |-----------------------|-----------------------------------|
  |  degree               |  n                                |
  |-----------------------|-----------------------------------|
  |  EBasis               |  BASIS_FEM_FIAT                   |
  |-----------------------|-----------------------------------|
  |  ECoordinates         |  COORDINATES_CARTESIAN            |
  |=======================|===================================|
  \endverbatim
  */
  
template<class Scalar, class ArrayScalar> 
class Basis_HGRAD_TET_Cn_FEM: public Basis<Scalar, ArrayScalar> {
private:
  
  /** \brief  Initializes <var>tagToOrdinal_</var> and <var>ordinalToTag_</var> lookup arrays.
   */
  virtual void initializeTags();

  /** \brief  The orthogonal basis on triangles, out of which the nodal basis is constructed
   */
  Basis_HGRAD_TET_Cn_FEM_ORTH<Scalar,FieldContainer<Scalar> > Phis;
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
  Basis_HGRAD_TET_Cn_FEM( const int n , const EPointType pointType );    
  
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

#include "Intrepid_HGRAD_TET_Cn_FEMDef.hpp"

#endif
