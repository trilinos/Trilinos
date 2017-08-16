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

/** \file   Intrepid_TensorProductSpaceTools.hpp
    \brief  Header file for the Intrepid2::TensorProductSpaceTools class.
    \author Created by R. Kirby
*/

#ifndef INTREPID2_TENSORPRODUCTSPACETOOLS_HPP
#define INTREPID2_TENSORPRODUCTSPACETOOLS_HPP

#include "Intrepid2_Basis.hpp"
#include "Intrepid2_FieldContainer.hpp"
#include "Teuchos_Array.hpp"
#include "Teuchos_RCP.hpp"

using Teuchos::Array;
using Teuchos::RCP;
using Intrepid2::FieldContainer;

namespace Intrepid2 {

/** \class Intrepid2::TensorProductSpaceTools
    \brief Defines expert-level interfaces for the evaluation, differentiation
           and integration of finite element-functions defined by tensor
           products of one-dimensional spaces.  These are useful in
           implementing spectral element methods.
*/
class TensorProductSpaceTools
{
public:
  /** \brief Computes point values of a set of polynomials expressed in a
             tensor product basis at output points.  The array
             <b>coeffs</b> is assumed to have dimensions (C,F1,F2),
	     where F1 runs over the number of different polynomials per cell and
             F2 runs over the coefficients run over a tensor product
             basis (lowest space dimension runs fastest).  The
	     Teuchos::Array of (pointers to) Arrays bases
	     have the one-dimensional bases tabulated at the
             one-dimensional points.  The output array is (C,F1,P).
	     
       \param vals  [out] - output point values of the discrete function
       \param coeffs [in] - coefficients of the input function
       \param bases  [in] - one-dimensional bases tabulated at points
   */
  template<class Scalar, class ArrayTypeOut, class ArrayTypeCoeffs,
	   class ArrayTypeBasis>
  static void evaluate(  ArrayTypeOut &vals ,
			 const ArrayTypeCoeffs &coeffs ,
			 const Array<RCP<ArrayTypeBasis> > &bases );

  // /** \brief Computes point values of a set of array-valued
  //            polynomials expressed in a  
  //            tensor product basis at output points.  The array
  //            <b>coeffs</b> is assumed to have dimensions (C,F1,F2),
  // 	     where F1 runs over the number of different polynomials per cell and
  //            F2 runs over the coefficients run over a tensor product
  //            basis (lowest space dimension runs fastest).  The
  // 	     Teuchos::Array of (pointers to) Arrays bases
  // 	     have the one-dimensional bases tabulated at the
  //            one-dimensional points.  The output array is (C,F1,P,D).
  // 	     This method assumes that the nodes for the basis coincide with
  // 	     the evaluation points, which leads to a big simplification.
	     
  //      \param vals  [out] - output point values of the discrete function
  //      \param coeffs [in] - coefficients of the input function
  //      \param bases  [in] - one-dimensional bases tabulated at points
  //  */
  // template<class Scalar, class ArrayTypeOut, class ArrayTypeCoeffs,
  // 	   class ArrayTypeBasis>
  // static void evaluateCollocated(  ArrayTypeOut &vals ,
  // 				   const ArrayTypeCoeffs &coeffs ,
  // 				   const Array<Array<RCP<ArrayTypeBasis> > > &bases );

  /** \brief Computes point values of a set of polynomials expressed in a
             tensor product basis at output points.  The array
             <b>coeffs</b> is assumed to have dimensions (C,F1,F2),
  	     where F1 runs over the number of different polynomials per cell and
             F2 runs over the coefficients run over a tensor product
             basis (lowest space dimension runs fastest).  The
  	     Teuchos::Array of (pointers to) Arrays bases
  	     have the one-dimensional bases tabulated at the
             one-dimensional points.  The output array is (C,F1,P).
  	     This method assumes that the nodes for the basis coincide with
  	     the evaluation points, which leads to a big simplification.
	     
       \param vals  [out] - output point values of the discrete function
       \param coeffs [in] - coefficients of the input function
       \param bases  [in] - one-dimensional bases tabulated at points
   */
  template<class Scalar, class ArrayTypeOut, class ArrayTypeCoeffs,
  	   class ArrayTypeBasis>
  static void evaluateCollocated(  ArrayTypeOut &vals ,
  				   const ArrayTypeCoeffs &coeffs ,
  				   const Array<RCP<ArrayTypeBasis> > &bases );

  /** \brief Given a polynomial expressed in a tensor product basis,
      evaluates the gradient at a tensor product of points.
      The array <b>coeffs</b> is assumed to have dimensions (C,F1,F2),
      where F1 runs over the number of different polynomials per cell and
      F2 runs over the coefficients run over a tensor product
      basis (lowest space dimension runs fastest).  The Teuchos::Array of
      (pointers to) Arrays bases and Dbases have the one-dimensional
      bases and their derivatives, respectively, tabulated at the
      one-dimensional points.  The output array is (C,F1,P).

      \param vals  [out] - output point values of the discrete function
      \param coeffs [in] - coefficients of the input function
      \param bases  [in] - one-dimensional bases tabulated at points
      \param Dbases  [in] - one-dimensional bases differentiated at points
   */
  template<class Scalar, class ArrayTypeOut, class ArrayTypeCoeffs,
	   class ArrayTypeBasis>
  static void evaluateGradient(  ArrayTypeOut &vals ,
				 const ArrayTypeCoeffs &coeffs ,
				 const Array<RCP<ArrayTypeBasis> > &bases ,
				 const Array<RCP<ArrayTypeBasis> > &Dbases );

  /** \brief Given a polynomial expressed in a tensor product basis,
      evaluates the gradient at a tensor product of points.
      The array <b>coeffs</b> is assumed to have dimensions (C,F1,F2),
      where F1 runs over the number of different polynomials per cell and
      F2 runs over the coefficients run over a tensor product
      basis (lowest space dimension runs fastest).  The Teuchos::Array of
      (pointers to) Arrays bases and Dbases have the one-dimensional
      bases and their derivatives, respectively, tabulated at the
      one-dimensional points.  The output array is (C,F1,P).
      This method assumes that the basis nodes and the output points coincide,
      which leads to a considerable simplification.

      \param vals  [out] - output point values of the discrete function
      \param coeffs [in] - coefficients of the input function
      \param bases  [in] - one-dimensional bases tabulated at points
      \param Dbases  [in] - one-dimensional bases differentiated at points
   */
  template<class Scalar, class ArrayTypeOut, class ArrayTypeCoeffs,
	   class ArrayTypeBasis>
  static void evaluateGradientCollocated( ArrayTypeOut &vals ,
					  const ArrayTypeCoeffs &coeffs ,
					  const Array<RCP<ArrayTypeBasis> > &bases ,
					  const Array<RCP<ArrayTypeBasis> > &Dbases );

  /** \brief Computes the moments of a set of data integrated against
      a basis tabulated at points.

      \param vals  [out]    - (C,F1,F2) output moments of the data against
                              the basis functions
      \param data [in]      - (C,F1,P) data tabulated at the tensor product
                              of points
      \param basisvals [in] - one-dimensional bases tabulated at points
      \param wts [in]       - array of one-dimensional quadrature weights
   */
  template<class Scalar, class ArrayTypeOut, class ArrayTypeData,
	   class ArrayTypeBasis, class ArrayTypeWeights>
  static void moments( ArrayTypeOut &vals ,
		       const ArrayTypeData &data ,
		       const Array<RCP<ArrayTypeBasis> > &basisVals ,
		       const Array<RCP<ArrayTypeWeights> > &wts );

  /** \brief Computes the moments of a set of data integrated against
      a basis tabulated at points, assuming that the basis nodes
      and integration points coincide.

      \param vals  [out]    - (C,F1,F2) output moments of the data against
                              the basis functions
      \param data [in]      - (C,F1,P) data tabulated at the tensor product
                              of points
      \param basisvals [in] - one-dimensional bases tabulated at points
      \param wts [in]       - array of one-dimensional quadrature weights
   */
  template<class Scalar, class ArrayTypeOut, class ArrayTypeData,
	   class ArrayTypeBasis, class ArrayTypeWeights>
  static void momentsCollocated( ArrayTypeOut &vals ,
				 const ArrayTypeData &data ,
				 const Array<RCP<ArrayTypeBasis> > &basisVals ,
		       const Array<RCP<ArrayTypeWeights> > &wts );

  /** \brief Computes the moments of a collection of F1 data integrated against
      a list of functions tabulated at points. F1 runs over the input data,
      F2 runs over the members of the basis.

      \param vals  [out]    - (C,F1,F2) output moments of the data against
                              the basis functions
      \param data [in]      - (C,F1,P,D) data tabulated at the tensor product
                              of points
      \param basisvals [in] - one-dimensional bases tabulated at points
      \param basisDvals [in] - one-dimensional bases differentated at points
      \param wts [in]     - one-dimensional quadrature weights
   */
  template<class Scalar, class ArrayTypeOut, class ArrayTypeData,
	   class ArrayTypeBasis, class ArrayTypeWeights>
  static void momentsGrad( ArrayTypeOut &vals ,
			   const ArrayTypeData &data ,
			   const Array<RCP<ArrayTypeBasis> > &basisVals ,
			   const Array<RCP<ArrayTypeBasis> > &basisDVals ,
			   const Array<RCP<ArrayTypeWeights> > &wts );

  /** \brief Computes the moments of a collection of F1 data integrated against
      a list of functions tabulated at points. F1 runs over the input data,
      F2 runs over the members of the basis.  This assumes the basis nodes
      and integration points coincide.

      \param vals  [out]    - (C,F1,F2) output moments of the data against
                              the basis functions
      \param data [in]      - (C,F1,P,D) data tabulated at the tensor product
                              of points
      \param basisvals [in] - one-dimensional bases tabulated at points
      \param basisDvals [in] - one-dimensional bases differentated at points
      \param wts [in]     - one-dimensional quadrature weights
   */
  template<class Scalar, class ArrayTypeOut, class ArrayTypeData,
	   class ArrayTypeBasis, class ArrayTypeWeights>
  static void momentsGradCollocated( ArrayTypeOut &vals ,
				     const ArrayTypeData &data ,
				     const Array<RCP<ArrayTypeBasis> > &basisVals ,
				     const Array<RCP<ArrayTypeBasis> > &basisDVals ,
				     const Array<RCP<ArrayTypeWeights> > &wts );


private:
  template<class Scalar, class ArrayTypeOut, class ArrayTypeCoeffs,
	   class ArrayTypeBasis>
  static void evaluate2D( ArrayTypeOut &vals ,
			  const ArrayTypeCoeffs &coeffs ,
			  const Array<RCP<ArrayTypeBasis> > &basisVals );

  template<class Scalar, class ArrayTypeOut, class ArrayTypeCoeffs,
	   class ArrayTypeBasis>
  static void evaluate3D( ArrayTypeOut &vals ,
			  const ArrayTypeCoeffs &coeffs ,
			  const Array<RCP<ArrayTypeBasis> > &basisDVals );

  template<class Scalar, class ArrayTypeOut, class ArrayTypeCoeffs,
	   class ArrayTypeBasis>
  static void evaluateCollocated2D( ArrayTypeOut &vals ,
				    const ArrayTypeCoeffs &coeffs ,
				    const Array<RCP<ArrayTypeBasis> > &basisVals );

  template<class Scalar, class ArrayTypeOut, class ArrayTypeCoeffs,
	   class ArrayTypeBasis>
  static void evaluateCollocated3D( ArrayTypeOut &vals ,
				    const ArrayTypeCoeffs &coeffs ,
				    const Array<RCP<ArrayTypeBasis> > &basisDVals );

  template<class Scalar, class ArrayTypeOut, class ArrayTypeCoeffs,
	   class ArrayTypeBasis>
  static void evaluateGradient2D( ArrayTypeOut &vals ,
				  const ArrayTypeCoeffs &coeffs ,
				  const Array<RCP<ArrayTypeBasis> > &basisVals ,
				  const Array<RCP<ArrayTypeBasis> > &basisDVals );

  template<class Scalar, class ArrayTypeOut, class ArrayTypeCoeffs,
	   class ArrayTypeBasis>
  static void evaluateGradient3D( ArrayTypeOut &vals ,
				  const ArrayTypeCoeffs &coeffs ,
				  const Array<RCP<ArrayTypeBasis> > &basisVals ,
				  const Array<RCP<ArrayTypeBasis> > &basisDVals );


  template<class Scalar, class ArrayTypeOut, class ArrayTypeCoeffs,
	   class ArrayTypeBasis>
  static void evaluateGradientCollocated2D( ArrayTypeOut &vals ,
				  const ArrayTypeCoeffs &coeffs ,
				  const Array<RCP<ArrayTypeBasis> > &basisVals ,
				  const Array<RCP<ArrayTypeBasis> > &basisDVals );

  template<class Scalar, class ArrayTypeOut, class ArrayTypeCoeffs,
	   class ArrayTypeBasis>
  static void evaluateGradientCollocated3D( ArrayTypeOut &vals ,
				  const ArrayTypeCoeffs &coeffs ,
				  const Array<RCP<ArrayTypeBasis> > &basisVals ,
				  const Array<RCP<ArrayTypeBasis> > &basisDVals );

  template<class Scalar, class ArrayTypeOut, class ArrayTypeData,
	   class ArrayTypeBasis, class ArrayTypeWeights>
  static void moments2D( ArrayTypeOut &vals ,
			 const ArrayTypeData &data ,
			 const Array<RCP<ArrayTypeBasis> > &basisVals ,
			 const Array<RCP<ArrayTypeWeights> > &wts );

  template<class Scalar, class ArrayTypeOut, class ArrayTypeData,
	   class ArrayTypeBasis, class ArrayTypeWeights>
  static void moments3D( ArrayTypeOut &vals ,
			 const ArrayTypeData &data ,
			 const Array<RCP<ArrayTypeBasis> > &basisVals ,
			 const Array<RCP<ArrayTypeWeights> > &wts );

  template<class Scalar, class ArrayTypeOut, class ArrayTypeData,
	   class ArrayTypeBasis, class ArrayTypeWeights>
  static void momentsCollocated2D( ArrayTypeOut &vals ,
			 const ArrayTypeData &data ,
			 const Array<RCP<ArrayTypeBasis> > &basisVals ,
			 const Array<RCP<ArrayTypeWeights> > &wts );

  template<class Scalar, class ArrayTypeOut, class ArrayTypeData,
	   class ArrayTypeBasis, class ArrayTypeWeights>
  static void momentsCollocated3D( ArrayTypeOut &vals ,
			 const ArrayTypeData &data ,
			 const Array<RCP<ArrayTypeBasis> > &basisVals ,
			 const Array<RCP<ArrayTypeWeights> > &wts );

  template<class Scalar, class ArrayTypeOut, class ArrayTypeData,
	   class ArrayTypeBasis, class ArrayTypeWeights>
  static void momentsGradCollocated2D( ArrayTypeOut &vals ,
			     const ArrayTypeData &data ,
			     const Array<RCP<ArrayTypeBasis> > &basisVals ,
			     const Array<RCP<ArrayTypeBasis> > &basisDVals ,
			     const Array<RCP<ArrayTypeWeights> > &wts );

  template<class Scalar, class ArrayTypeOut, class ArrayTypeData,
	   class ArrayTypeBasis, class ArrayTypeWeights>
  static void momentsGradCollocated3D( ArrayTypeOut &vals ,
			     const ArrayTypeData &data ,
			     const Array<RCP<ArrayTypeBasis> > &basisVals ,
			     const Array<RCP<ArrayTypeBasis> > &basisDVals ,
			     const Array<RCP<ArrayTypeWeights> > &wts );

 template<class Scalar, class ArrayTypeOut, class ArrayTypeData,
	   class ArrayTypeBasis, class ArrayTypeWeights>
  static void momentsGrad2D( ArrayTypeOut &vals ,
			     const ArrayTypeData &data ,
			     const Array<RCP<ArrayTypeBasis> > &basisVals ,
			     const Array<RCP<ArrayTypeBasis> > &basisDVals ,
			     const Array<RCP<ArrayTypeWeights> > &wts );

  template<class Scalar, class ArrayTypeOut, class ArrayTypeData,
	   class ArrayTypeBasis, class ArrayTypeWeights>
  static void momentsGrad3D( ArrayTypeOut &vals ,
			     const ArrayTypeData &data ,
			     const Array<RCP<ArrayTypeBasis> > &basisVals ,
			     const Array<RCP<ArrayTypeBasis> > &basisDVals ,
			     const Array<RCP<ArrayTypeWeights> > &wts );
};

} //end namespace Intrepid2

#include "Intrepid2_TensorProductSpaceToolsDef.hpp"
#endif 
