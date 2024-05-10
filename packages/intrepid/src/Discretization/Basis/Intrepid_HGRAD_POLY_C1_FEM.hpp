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
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301
// USA
// Questions? Contact Pavel Bochev (pbboche@sandia.gov) or
//                    Denis Ridzal (dridzal@sandia.gov).
//
// ************************************************************************
// @HEADER

/** \file   Intrepid_HGRAD_POLY_C1_FEM.hpp
    \brief  Header file for the Intrepid::HGRAD_POLY_C1_FEM class.
    \author Created by P. Bochev and J. Lai.
*/

#ifndef INTREPID_HGRAD_POLY_C1_FEM_HPP
#define INTREPID_HGRAD_POLY_C1_FEM_HPP

#include "Intrepid_Basis.hpp"
#include "Shards_CellTopology.hpp"

namespace Intrepid{
  /** \class Intrepid::HGRAD_POLY_C1_FEM
      \brief Implementation of the default H(grad) compatible FEM basis of degree 1 on a polygon cell

             Implements Wachspress rational basis on polygons.  The basis has cardinality equal to the
	     number of vertices.  
  */
  template<class Scalar, class ArrayScalar> 
  class Basis_HGRAD_POLY_C1_FEM : public Basis<Scalar, ArrayScalar> {
  public:
    /** \brief  Constructor.
	
	\param  cellTopology          [in]     - the topology of the polygon
    */
    Basis_HGRAD_POLY_C1_FEM(const shards::CellTopology& cellTopology);
    
    /** \brief  FEM reference basis evaluation: invocation of this method throws an exception.
     */
    void getValues(ArrayScalar& outputValues,
		   const ArrayScalar& inputPoints,
		   const EOperator operatorType) const;

    /** \brief  Evaluation of a FEM basis on a <strong>physical polygon</strong> cell.

	\param  outputValues         [out]  - variable rank array with the basis values
	\param  inputPoints          [in]   - rank-2 array (P,D) with the evaluation points
	\param  cellVertices         [in]   - rank-2 array (V,D) with the vertices of the polygon
	\param  operatorType         [in]   - the operator acting on the basis functions
    */
    void getValues(ArrayScalar& outputValues,
		   const ArrayScalar& inputPoints,
		   const ArrayScalar& cellVertices,
		   const EOperator operatorType = OPERATOR_VALUE) const;
    
  private:
    
    /** \brief Initializes <var> tagToOrdinal_</var> and <var>ordinalToTag_</var> lookup arrays.
     */
    void initializeTags();

    
    /** \brief  Helper function to compute area of triangle formed by 3 points
     */
    template<class Scalar1, class ArrayScalar1>
    Scalar1 computeArea(const ArrayScalar1& p1,
		       const ArrayScalar1& p2,
		       const ArrayScalar1& p3) const;

    /** \brief  Evaluation of the Wachspress weight functions.
     */
    template<class Scalar1, class ArrayScalar1>
    void evaluateWeightFunctions(ArrayScalar1& outputValues,
				 const ArrayScalar1& inputValues,
				 const ArrayScalar1& cellVertices) const;
    
    
    
    /** \brief  Evaluation of Wachspress shape functions.
     */
    template<class Scalar1, class ArrayScalar1>
    void shapeFunctions(ArrayScalar1& outputValues,
			const ArrayScalar1& inputValues,
			const ArrayScalar1& cellVertices)const;
  }; // end class Basis_HGRAD_POLY_C1_FEM

} // namespace Intrepid

#include "Intrepid_HGRAD_POLY_C1_FEMDef.hpp"
#endif

#if defined(Intrepid_SHOW_DEPRECATED_WARNINGS)
#ifdef __GNUC__
#warning "The Intrepid package is deprecated"
#endif
#endif

