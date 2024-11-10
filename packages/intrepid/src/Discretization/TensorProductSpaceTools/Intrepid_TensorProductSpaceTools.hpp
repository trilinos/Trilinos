// @HEADER
// ************************************************************************
//
//                           Intrepid Package
//                 Copytest (2007) Sandia Corporation
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

/** \file   Intrepid_TensorProductSpaceTools.hpp
    \brief  Header file for the Intrepid::TensorProductSpaceTools class.
    \author Created by R. Kirby
*/

#ifndef INTREPID_TENSORPRODUCTSPACETOOLS_HPP
#define INTREPID_TENSORPRODUCTSPACETOOLS_HPP

#include "Intrepid_Basis.hpp"
#include "Intrepid_FieldContainer.hpp"
#include "Teuchos_Array.hpp"
#include "Teuchos_RCP.hpp"

using Teuchos::Array;
using Teuchos::RCP;
using Intrepid::FieldContainer;

namespace Intrepid {

/** \class Intrepid::TensorProductSpaceTools
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

} //end namespace Intrepid

#include "Intrepid_TensorProductSpaceToolsDef.hpp"
#endif 

#if defined(Intrepid_SHOW_DEPRECATED_WARNINGS)
#ifdef __GNUC__
#warning "The Intrepid package is deprecated"
#endif
#endif

