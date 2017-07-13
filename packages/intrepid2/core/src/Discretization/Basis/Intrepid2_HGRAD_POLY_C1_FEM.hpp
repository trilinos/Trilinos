// @HEADER
// ************************************************************************
//
//                           Intrepid2 Package
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
// Questions? Contact Kyungjoo Kim  (kyukim@sandia.gov), or
//                    Mauro Perego  (mperego@sandia.gov)
//
// ************************************************************************
// @HEADER

/** \file   Intrepid_HGRAD_POLY_C1_FEM.hpp
    \brief  Header file for the Intrepid2::HGRAD_POLY_C1_FEM class.
    \author Created by P. Bochev and J. Lai.
*/

#ifndef INTREPID2_HGRAD_POLY_C1_FEM_HPP
#define INTREPID2_HGRAD_POLY_C1_FEM_HPP

#include "Intrepid2_Basis.hpp"
#include "Shards_CellTopology.hpp"

namespace Intrepid2{
  /** \class Intrepid2::HGRAD_POLY_C1_FEM
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

} // namespace Intrepid2

#include "Intrepid2_HGRAD_POLY_C1_FEMDef.hpp"
#endif
