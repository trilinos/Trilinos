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

/** \file   Intrepid_TensorProductSpaceToolsDef.hpp
    \brief  Definition file for the Intrepid::TensorProductSpaceTools class.
    \author Created by R. Kirby
*/


namespace Intrepid {
  template<class Scalar, class ArrayTypeOut, class ArrayTypeCoeffs,
	   class ArrayTypeBasis>
  void TensorProductSpaceTools::evaluate( ArrayTypeOut &vals ,
					  const ArrayTypeCoeffs &coeffs ,
					  const Array<RCP<ArrayTypeBasis> > &bases )
  {
    const unsigned space_dim = bases.size();

#ifdef HAVE_INTREPID_DEBUG
    TEUCHOS_TEST_FOR_EXCEPTION( vals.rank() != 3 , std::invalid_argument ,
				">>> ERROR (TensorProductSpaceTools::evaluate): vals must be rank 3 array." );

    TEUCHOS_TEST_FOR_EXCEPTION( coeffs.rank() != 3 , std::invalid_argument ,
				">>> ERROR (TensorProductSpaceTools::evaluate): coeffs must be rank 3 array." );

    for (unsigned d=0;d<space_dim;d++)
      {
	TEUCHOS_TEST_FOR_EXCEPTION( bases[d]->rank() != 2 , std::invalid_argument ,
				    ">>> ERROR (TensorProductSpaceTools::evaluate): each tabulated basis must be rank 2 array." );
      }

    TEUCHOS_TEST_FOR_EXCEPTION( vals.dimension(0) != coeffs.dimension(0) , std::invalid_argument,
				">>> ERROR (TensorProductSpaceTools::evaluate): Number of cells for vals and coeffs must match." );

    TEUCHOS_TEST_FOR_EXCEPTION( vals.dimension(1) != coeffs.dimension(1) , std::invalid_argument,
				">>> ERROR (TensorProductSpaceTools::evaluate): Number of fields for vals and coeffs must match." );
    
    int product_of_basis_dimensions = 1;
    int product_of_basis_points = 1;
    for (unsigned d=0;d<space_dim;d++)
      {
	product_of_basis_dimensions *= bases[d]->dimension(0);
	product_of_basis_points *= bases[d]->dimension(1);
      }
    
    TEUCHOS_TEST_FOR_EXCEPTION( vals.dimension(2) != product_of_basis_points , std::invalid_argument,
				">>> ERROR (TensorProductSpaceTools::evaluate): Incompatible number of points in vals and bases ." );
    
    TEUCHOS_TEST_FOR_EXCEPTION( coeffs.dimension(2) != product_of_basis_dimensions , std::invalid_argument,
				">>> ERROR (TensorProductSpaceTools::evaluate): Incompatible number of basis functions in coeffs and bases ." );
#endif    
    

    switch (space_dim)
      {
      case 2:
	evaluate2D<Scalar,ArrayTypeOut,ArrayTypeCoeffs,ArrayTypeBasis>( vals , coeffs , bases );
	break;
      case 3:
	evaluate3D<Scalar,ArrayTypeOut,ArrayTypeCoeffs,ArrayTypeBasis>( vals , coeffs , bases );
	break;
      }

  }

  template<class Scalar, class ArrayTypeOut, class ArrayTypeCoeffs,
	   class ArrayTypeBasis>
  void TensorProductSpaceTools::evaluateCollocated( ArrayTypeOut &vals ,
						    const ArrayTypeCoeffs &coeffs ,
						    const Array<RCP<ArrayTypeBasis> > &bases )
  {
    const unsigned space_dim = bases.size();

#ifdef HAVE_INTREPID_DEBUG
    TEUCHOS_TEST_FOR_EXCEPTION( vals.rank() != 3 , std::invalid_argument ,
				">>> ERROR (TensorProductSpaceTools::evaluateCollocated): vals must be rank 3 array." );

    TEUCHOS_TEST_FOR_EXCEPTION( coeffs.rank() != 3 , std::invalid_argument ,
				">>> ERROR (TensorProductSpaceTools::evaluateCollocated): coeffs must be rank 3 array." );

    for (unsigned d=0;d<space_dim;d++)
      {
	TEUCHOS_TEST_FOR_EXCEPTION( bases[d]->rank() != 2 , std::invalid_argument ,
				    ">>> ERROR (TensorProductSpaceTools::evaluateCollocated): each tabulated basis must be rank 2 array." );
      }

    TEUCHOS_TEST_FOR_EXCEPTION( vals.dimension(0) != coeffs.dimension(0) , std::invalid_argument,
				">>> ERROR (TensorProductSpaceTools::evaluateCollocated): Number of cells for vals and coeffs must match." );

    TEUCHOS_TEST_FOR_EXCEPTION( vals.dimension(1) != coeffs.dimension(1) , std::invalid_argument,
				">>> ERROR (TensorProductSpaceTools::evaluateCollocated): Number of fields for vals and coeffs must match." );
    
    int product_of_basis_dimensions = 1;
    int product_of_basis_points = 1;
    for (unsigned d=0;d<space_dim;d++)
      {
	product_of_basis_dimensions *= bases[d]->dimension(0);
	product_of_basis_points *= bases[d]->dimension(1);
      }
    
    TEUCHOS_TEST_FOR_EXCEPTION( vals.dimension(2) != product_of_basis_points , std::invalid_argument,
				">>> ERROR (TensorProductSpaceTools::evaluateCollocated): Incompatible number of points in vals and bases ." );
    
    TEUCHOS_TEST_FOR_EXCEPTION( coeffs.dimension(2) != product_of_basis_dimensions , std::invalid_argument,
				">>> ERROR (TensorProductSpaceTools::evaluateCollocated): Incompatible number of basis functions in coeffs and bases ." );
#endif    
    

    switch (space_dim)
      {
      case 2:
	evaluateCollocated2D<Scalar,ArrayTypeOut,ArrayTypeCoeffs,ArrayTypeBasis>( vals , coeffs , bases );
	break;
      case 3:
	evaluateCollocated3D<Scalar,ArrayTypeOut,ArrayTypeCoeffs,ArrayTypeBasis>( vals , coeffs , bases );
	break;
      }

  }

//   template<class Scalar, class ArrayTypeOut, class ArrayTypeCoeffs,
// 	   class ArrayTypeBasis>
//   void TensorProductSpaceTools::evaluateCollocated( ArrayTypeOut &vals ,
// 						    const ArrayTypeCoeffs &coeffs ,
// 						    const Array<Array<RCP<ArrayTypeBasis> > > &bases )
//   {
//     const unsigned space_dim = bases.size();

// #ifdef HAVE_INTREPID_DEBUG
//     TEUCHOS_TEST_FOR_EXCEPTION( vals.rank() != 4 , std::invalid_argument ,
// 				">>> ERROR (TensorProductSpaceTools::evaluateCollocated): vals must be rank 3 array." );

//     TEUCHOS_TEST_FOR_EXCEPTION( coeffs.rank() != 3 , std::invalid_argument ,
// 				">>> ERROR (TensorProductSpaceTools::evaluateCollocated): coeffs must be rank 3 array." );

//     for (unsigned d0=0;d0<bases.size();d0++)
//       {
// 	for (unsigned d1=1;d1<space_dim;d1++)
// 	  {
// 	    TEUCHOS_TEST_FOR_EXCEPTION( bases[d0][d1]->rank() != 2 , std::invalid_argument ,
// 					">>> ERROR (TensorProductSpaceTools::evaluateCollocated): each tabulated basis must be rank 2 array." );
// 	  }
//       }

//     TEUCHOS_TEST_FOR_EXCEPTION( vals.dimension(0) != coeffs.dimension(0) , std::invalid_argument,
// 				">>> ERROR (TensorProductSpaceTools::evaluateCollocated): Number of cells for vals and coeffs must match." );

//     TEUCHOS_TEST_FOR_EXCEPTION( vals.dimension(1) != coeffs.dimension(1) , std::invalid_argument,
// 				">>> ERROR (TensorProductSpaceTools::evaluateCollocated): Number of fields for vals and coeffs must match." );
    
//     int product_of_basis_dimensions = 1;
//     int product_of_basis_points = 1;
//     for (unsigned d=0;d<space_dim;d++)
//       {
// 	product_of_basis_dimensions *= bases[0][d]->dimension(0);
// 	product_of_basis_points *= bases[0][d]->dimension(1);
//       }
    
//     TEUCHOS_TEST_FOR_EXCEPTION( vals.dimension(2) != product_of_basis_points , std::invalid_argument,
// 				">>> ERROR (TensorProductSpaceTools::evaluateCollocated): Incompatible number of points in vals and bases ." );
    
//     TEUCHOS_TEST_FOR_EXCEPTION( coeffs.dimension(2) != product_of_basis_dimensions , std::invalid_argument,
// 				">>> ERROR (TensorProductSpaceTools::evaluateCollocated): Incompatible number of basis functions in coeffs and bases ." );
// #endif    

//     TEUCHOS_TEST_FOR_EXCEPTION( vals.rank(3) != bases.size() , std::invalid_argument ,
// 				">>> ERROR (TensorProductSpaceTools::evaluateCollocated): wrong dimensions for vals");
    

//     switch (space_dim)
//       {
//       case 2:
// 	evaluateCollocated2D<Scalar,ArrayTypeOut,ArrayTypeCoeffs,ArrayTypeBasis>( vals , coeffs , bases );
// 	break;
//       case 3:
// 	evaluateCollocated3D<Scalar,ArrayTypeOut,ArrayTypeCoeffs,ArrayTypeBasis>( vals , coeffs , bases );
// 	break;
//       }

//   }

  template<class Scalar, class ArrayTypeOut, class ArrayTypeCoeffs,
	   class ArrayTypeBasis>
  void TensorProductSpaceTools::evaluateGradient( ArrayTypeOut &vals ,
						  const ArrayTypeCoeffs &coeffs ,
						  const Array<RCP<ArrayTypeBasis> > &bases ,
						  const Array<RCP<ArrayTypeBasis> > &Dbases )
  {
    const unsigned space_dim = bases.size();

#ifdef HAVE_INTREPID_DEBUG
    TEUCHOS_TEST_FOR_EXCEPTION( vals.rank() != 4 , std::invalid_argument ,
				">>> ERROR (TensorProductSpaceTools::evaluate): vals must be rank 3 array." );

    TEUCHOS_TEST_FOR_EXCEPTION( coeffs.rank() != 3 , std::invalid_argument ,
				">>> ERROR (TensorProductSpaceTools::evaluate): coeffs must be rank 3 array." );

    TEUCHOS_TEST_FOR_EXCEPTION( Dbases.size() != (int)space_dim , std::invalid_argument ,
				">>> ERROR (TensorProductSpaceTools::evaluate): bases and Dbases must be same size." );
   
    for (unsigned d=0;d<space_dim;d++)
      {
	TEUCHOS_TEST_FOR_EXCEPTION( bases[d]->rank() != 2 , std::invalid_argument ,
				    ">>> ERROR (TensorProductSpaceTools::evaluate): each tabulated basis must be rank 2 array." );

	TEUCHOS_TEST_FOR_EXCEPTION( Dbases[d]->rank() != 3 , std::invalid_argument ,
				    ">>> ERROR (TensorProductSpaceTools::evaluate): each tabulated Dbasis must be rank 3 array." );
      }

    TEUCHOS_TEST_FOR_EXCEPTION( vals.dimension(0) != coeffs.dimension(0) , std::invalid_argument,
				">>> ERROR (TensorProductSpaceTools::evaluate): Number of cells for vals and coeffs must match." );

    TEUCHOS_TEST_FOR_EXCEPTION( vals.dimension(1) != coeffs.dimension(1) , std::invalid_argument,
				">>> ERROR (TensorProductSpaceTools::evaluate): Number of fields for vals and coeffs must match." );
    
    int product_of_basis_dimensions = 1;
    int product_of_basis_points = 1;
    for (unsigned d=0;d<space_dim;d++)
      {
	product_of_basis_dimensions *= bases[d]->dimension(0);
	product_of_basis_points *= bases[d]->dimension(1);
      }
    
    TEUCHOS_TEST_FOR_EXCEPTION( vals.dimension(2) != product_of_basis_points , std::invalid_argument,
				">>> ERROR (TensorProductSpaceTools::evaluate): Incompatible number of points in vals and bases ." );
    
    TEUCHOS_TEST_FOR_EXCEPTION( coeffs.dimension(2) != product_of_basis_dimensions , std::invalid_argument,
				">>> ERROR (TensorProductSpaceTools::evaluate): Incompatible number of basis functions in coeffs and bases ." );

    for (unsigned d=0;d<space_dim;d++)
      {
	for (unsigned i=0;i<2;i++)
	  {
	    TEUCHOS_TEST_FOR_EXCEPTION( bases[d]->dimension(i) != Dbases[d]->dimension(i) , std::invalid_argument ,
					">>>ERROR (TensorProductSpaceTools::evaluate): bases and Dbases have incompatible shape." );
	  }
      }
#endif    

    switch (space_dim)
      {
      case 2:
	TensorProductSpaceTools::evaluateGradient2D<Scalar,ArrayTypeOut,ArrayTypeCoeffs,ArrayTypeBasis>( vals , coeffs , bases , Dbases );
	break;
      case 3:
	TensorProductSpaceTools::evaluateGradient3D<Scalar,ArrayTypeOut,ArrayTypeCoeffs,ArrayTypeBasis>( vals , coeffs , bases , Dbases );
	break;
      }

  }

  template<class Scalar, class ArrayTypeOut, class ArrayTypeCoeffs,
	   class ArrayTypeBasis>
  void TensorProductSpaceTools::evaluateGradientCollocated( ArrayTypeOut &vals ,
							    const ArrayTypeCoeffs &coeffs ,
							    const Array<RCP<ArrayTypeBasis> > &bases ,
							    const Array<RCP<ArrayTypeBasis> > &Dbases )
  {
    const unsigned space_dim = bases.size();

#ifdef HAVE_INTREPID_DEBUG
    TEUCHOS_TEST_FOR_EXCEPTION( vals.rank() != 4 , std::invalid_argument ,
				">>> ERROR (TensorProductSpaceTools::evaluate): vals must be rank 3 array." );

    TEUCHOS_TEST_FOR_EXCEPTION( coeffs.rank() != 3 , std::invalid_argument ,
				">>> ERROR (TensorProductSpaceTools::evaluate): coeffs must be rank 3 array." );

    TEUCHOS_TEST_FOR_EXCEPTION( Dbases.size() != (int)space_dim , std::invalid_argument ,
				">>> ERROR (TensorProductSpaceTools::evaluate): bases and Dbases must be same size." );
   
    for (unsigned d=0;d<space_dim;d++)
      {
	TEUCHOS_TEST_FOR_EXCEPTION( bases[d]->rank() != 2 , std::invalid_argument ,
				    ">>> ERROR (TensorProductSpaceTools::evaluate): each tabulated basis must be rank 2 array." );

	TEUCHOS_TEST_FOR_EXCEPTION( Dbases[d]->rank() != 3 , std::invalid_argument ,
				    ">>> ERROR (TensorProductSpaceTools::evaluate): each tabulated Dbasis must be rank 3 array." );
      }

    TEUCHOS_TEST_FOR_EXCEPTION( vals.dimension(0) != coeffs.dimension(0) , std::invalid_argument,
				">>> ERROR (TensorProductSpaceTools::evaluate): Number of cells for vals and coeffs must match." );

    TEUCHOS_TEST_FOR_EXCEPTION( vals.dimension(1) != coeffs.dimension(1) , std::invalid_argument,
				">>> ERROR (TensorProductSpaceTools::evaluate): Number of fields for vals and coeffs must match." );
    
    int product_of_basis_dimensions = 1;
    int product_of_basis_points = 1;
    for (unsigned d=0;d<space_dim;d++)
      {
	product_of_basis_dimensions *= bases[d]->dimension(0);
	product_of_basis_points *= bases[d]->dimension(1);
      }
    
    TEUCHOS_TEST_FOR_EXCEPTION( vals.dimension(2) != product_of_basis_points , std::invalid_argument,
				">>> ERROR (TensorProductSpaceTools::evaluate): Incompatible number of points in vals and bases ." );
    
    TEUCHOS_TEST_FOR_EXCEPTION( coeffs.dimension(2) != product_of_basis_dimensions , std::invalid_argument,
				">>> ERROR (TensorProductSpaceTools::evaluate): Incompatible number of basis functions in coeffs and bases ." );

    for (unsigned d=0;d<space_dim;d++)
      {
	for (unsigned i=0;i<2;i++)
	  {
	    TEUCHOS_TEST_FOR_EXCEPTION( bases[d]->dimension(i) != Dbases[d]->dimension(i) , std::invalid_argument ,
					">>>ERROR (TensorProductSpaceTools::evaluate): bases and Dbases have incompatible shape." );
	  }
      }
#endif    

    switch (space_dim)
      {
      case 2:
	evaluateGradientCollocated2D<Scalar,ArrayTypeOut,ArrayTypeCoeffs,ArrayTypeBasis>( vals , coeffs , bases , Dbases );
	break;
      case 3:
	evaluateGradientCollocated3D<Scalar,ArrayTypeOut,ArrayTypeCoeffs,ArrayTypeBasis>( vals , coeffs , bases , Dbases );
	break;
      }

  }

  template<class Scalar, class ArrayTypeOut, class ArrayTypeData,
	   class ArrayTypeBasis, class ArrayTypeWeights>
  void TensorProductSpaceTools::moments( ArrayTypeOut &vals ,
					 const ArrayTypeData &data ,
					 const Array<RCP<ArrayTypeBasis> > &basisVals ,
					 const Array<RCP<ArrayTypeWeights> > &wts )
  {
    const unsigned space_dim = basisVals.size();

#ifdef HAVE_INTREPID_DEBUG
    TEUCHOS_TEST_FOR_EXCEPTION( vals.rank() != 3 , std::invalid_argument ,
				">>> ERROR (TensorProductSpaceTools::evaluate): vals must be rank 3 array." );

    TEUCHOS_TEST_FOR_EXCEPTION( data.rank() != 3 , std::invalid_argument ,
				">>> ERROR (TensorProductSpaceTools::evaluate): coeffs must be rank 3 array." );

    TEUCHOS_TEST_FOR_EXCEPTION( basisVals.size() != (int)space_dim , std::invalid_argument ,
				">>> ERROR (TensorProductSpaceTools::evaluate): bases ill-sized." );

    TEUCHOS_TEST_FOR_EXCEPTION( wts.size() != (int)space_dim , std::invalid_argument ,
				">>> ERROR (TensorProductSpaceTools::evaluate):  quadrature weights ill-sized." );
   
    for (unsigned d=0;d<space_dim;d++)
      {
	TEUCHOS_TEST_FOR_EXCEPTION( basisVals[d]->rank() != 2 , std::invalid_argument ,
				    ">>> ERROR (TensorProductSpaceTools::evaluate): each tabulated basis must be rank 2 array." );

	TEUCHOS_TEST_FOR_EXCEPTION( wts[d]->rank() != 1 , std::invalid_argument ,
				    ">>> ERROR (TensorProductSpaceTools::evaluate): each tabulated Dbasis must be rank 2 array." );
      }

    TEUCHOS_TEST_FOR_EXCEPTION( vals.dimension(0) != data.dimension(0) , std::invalid_argument,
				">>> ERROR (TensorProductSpaceTools::evaluate): Number of cells for vals and coeffs must match." );

    TEUCHOS_TEST_FOR_EXCEPTION( vals.dimension(1) != data.dimension(1) , std::invalid_argument,
				">>> ERROR (TensorProductSpaceTools::evaluate): Number of fields for vals and coeffs must match." );
    
    int product_of_basis_dimensions = 1;
    int product_of_basis_points = 1;
    for (unsigned d=0;d<space_dim;d++)
      {
	product_of_basis_dimensions *= basisVals[d]->dimension(0);
	product_of_basis_points *= basisVals[d]->dimension(1);
      }
    
    TEUCHOS_TEST_FOR_EXCEPTION( vals.dimension(2) != product_of_basis_dimensions , std::invalid_argument,
				">>> ERROR (TensorProductSpaceTools::evaluate): Incompatible number of points in vals and bases ." );
    
    TEUCHOS_TEST_FOR_EXCEPTION( data.dimension(2) != product_of_basis_points , std::invalid_argument,
				">>> ERROR (TensorProductSpaceTools::evaluate): Incompatible number of basis functions in data and bases ." );

#endif    

    switch (space_dim)
      {
      case 2:
	moments2D<Scalar, ArrayTypeOut, ArrayTypeData, ArrayTypeBasis, ArrayTypeWeights>( vals , data , basisVals , wts );
	break;
      case 3:
	moments3D<Scalar, ArrayTypeOut, ArrayTypeData, ArrayTypeBasis, ArrayTypeWeights>( vals , data , basisVals , wts );
	break;
      }

  }

  template<class Scalar, class ArrayTypeOut, class ArrayTypeData,
	   class ArrayTypeBasis, class ArrayTypeWeights>
  void TensorProductSpaceTools::momentsCollocated( ArrayTypeOut &vals ,
						   const ArrayTypeData &data ,
						   const Array<RCP<ArrayTypeBasis> > &basisVals ,
						   const Array<RCP<ArrayTypeWeights> > &wts )
  {
    const unsigned space_dim = basisVals.size();

#ifdef HAVE_INTREPID_DEBUG
    TEUCHOS_TEST_FOR_EXCEPTION( vals.rank() != 3 , std::invalid_argument ,
				">>> ERROR (TensorProductSpaceTools::evaluate): vals must be rank 3 array." );

    TEUCHOS_TEST_FOR_EXCEPTION( data.rank() != 3 , std::invalid_argument ,
				">>> ERROR (TensorProductSpaceTools::evaluate): coeffs must be rank 3 array." );

    TEUCHOS_TEST_FOR_EXCEPTION( basisVals.size() != (int)space_dim , std::invalid_argument ,
				">>> ERROR (TensorProductSpaceTools::evaluate): bases ill-sized." );

    TEUCHOS_TEST_FOR_EXCEPTION( wts.size() != (int)space_dim , std::invalid_argument ,
				">>> ERROR (TensorProductSpaceTools::evaluate):  quadrature weights ill-sized." );
   
    for (unsigned d=0;d<space_dim;d++)
      {
	TEUCHOS_TEST_FOR_EXCEPTION( basisVals[d]->rank() != 2 , std::invalid_argument ,
				    ">>> ERROR (TensorProductSpaceTools::evaluate): each tabulated basis must be rank 2 array." );

	TEUCHOS_TEST_FOR_EXCEPTION( wts[d]->rank() != 1 , std::invalid_argument ,
				    ">>> ERROR (TensorProductSpaceTools::evaluate): each tabulated Dbasis must be rank 2 array." );
      }

    TEUCHOS_TEST_FOR_EXCEPTION( vals.dimension(0) != data.dimension(0) , std::invalid_argument,
				">>> ERROR (TensorProductSpaceTools::evaluate): Number of cells for vals and coeffs must match." );

    TEUCHOS_TEST_FOR_EXCEPTION( vals.dimension(1) != data.dimension(1) , std::invalid_argument,
				">>> ERROR (TensorProductSpaceTools::evaluate): Number of fields for vals and coeffs must match." );
    
    int product_of_basis_dimensions = 1;
    int product_of_basis_points = 1;
    for (unsigned d=0;d<space_dim;d++)
      {
	product_of_basis_dimensions *= basisVals[d]->dimension(0);
	product_of_basis_points *= basisVals[d]->dimension(1);
      }
    
    TEUCHOS_TEST_FOR_EXCEPTION( vals.dimension(2) != product_of_basis_dimensions , std::invalid_argument,
				">>> ERROR (TensorProductSpaceTools::evaluate): Incompatible number of points in vals and bases ." );
    
    TEUCHOS_TEST_FOR_EXCEPTION( data.dimension(2) != product_of_basis_points , std::invalid_argument,
				">>> ERROR (TensorProductSpaceTools::evaluate): Incompatible number of basis functions in data and bases ." );

#endif    

    switch (space_dim)
      {
      case 2:
	moments2D<Scalar, ArrayTypeOut, ArrayTypeData, ArrayTypeBasis, ArrayTypeWeights>( vals , data , basisVals , wts );
	break;
      case 3:
	moments3D<Scalar, ArrayTypeOut, ArrayTypeData, ArrayTypeBasis, ArrayTypeWeights>( vals , data , basisVals , wts );
	break;
      }

  }

  template<class Scalar, class ArrayTypeOut, class ArrayTypeData,
	   class ArrayTypeBasis, class ArrayTypeWeights>
  void TensorProductSpaceTools::momentsGrad( ArrayTypeOut &vals ,
					     const ArrayTypeData &data ,
					     const Array<RCP<ArrayTypeBasis> > &basisVals ,
					     const Array<RCP<ArrayTypeBasis> > &basisDVals ,
					     const Array<RCP<ArrayTypeWeights> > &wts )
  {
    const unsigned space_dim = basisVals.size();

#ifdef HAVE_INTREPID_DEBUG
    TEUCHOS_TEST_FOR_EXCEPTION( vals.rank() != 3 , std::invalid_argument ,
				">>> ERROR (TensorProductSpaceTools::evaluate): vals must be rank 3 array." );

    TEUCHOS_TEST_FOR_EXCEPTION( data.rank() != 4 , std::invalid_argument ,
				">>> ERROR (TensorProductSpaceTools::evaluate): coeffs must be rank 4 array." );

    TEUCHOS_TEST_FOR_EXCEPTION( basisVals.size() != (int)space_dim , std::invalid_argument ,
				">>> ERROR (TensorProductSpaceTools::evaluate): bases ill-sized." );

    TEUCHOS_TEST_FOR_EXCEPTION( basisDVals.size() != (int)space_dim , std::invalid_argument ,
				">>> ERROR (TensorProductSpaceTools::evaluate): bases ill-sized." );

    TEUCHOS_TEST_FOR_EXCEPTION( wts.size() != (int)space_dim , std::invalid_argument ,
				">>> ERROR (TensorProductSpaceTools::evaluate):  quadrature weights ill-sized." );
   
    for (unsigned d=0;d<space_dim;d++)
      {
	TEUCHOS_TEST_FOR_EXCEPTION( basisVals[d]->rank() != 2 , std::invalid_argument ,
				    ">>> ERROR (TensorProductSpaceTools::evaluate): each tabulated basis must be rank 2 array." );

	TEUCHOS_TEST_FOR_EXCEPTION( basisDVals[d]->rank() != 3 , std::invalid_argument ,
				    ">>> ERROR (TensorProductSpaceTools::evaluate): each tabulated derivative basis must be rank 3 array." );

	TEUCHOS_TEST_FOR_EXCEPTION( wts[d]->rank() != 1 , std::invalid_argument ,
				    ">>> ERROR (TensorProductSpaceTools::evaluate): each tabulated Dbasis must be rank 2 array." );
      }

    TEUCHOS_TEST_FOR_EXCEPTION( vals.dimension(0) != data.dimension(0) , std::invalid_argument,
				">>> ERROR (TensorProductSpaceTools::evaluate): Number of cells for vals and coeffs must match." );
    
    int product_of_basis_dimensions = 1;
    int product_of_basis_points = 1;
    for (unsigned d=0;d<space_dim;d++)
      {
	product_of_basis_dimensions *= basisVals[d]->dimension(0);
	product_of_basis_points *= basisVals[d]->dimension(1);
      }
    
    TEUCHOS_TEST_FOR_EXCEPTION( vals.dimension(2) != product_of_basis_dimensions , std::invalid_argument,
				">>> ERROR (TensorProductSpaceTools::evaluate): Incompatible number of points in vals and bases ." );
    
    TEUCHOS_TEST_FOR_EXCEPTION( data.dimension(2) != product_of_basis_points , std::invalid_argument,
				">>> ERROR (TensorProductSpaceTools::evaluate): Incompatible number of basis functions in data and bases ." );

#endif    

    switch (space_dim)
      {
      case 2:
	momentsGrad2D<Scalar, ArrayTypeOut, ArrayTypeData, ArrayTypeBasis, ArrayTypeWeights>( vals , data , basisVals , basisDVals , wts );
	break;
      case 3:
	momentsGrad3D<Scalar, ArrayTypeOut, ArrayTypeData, ArrayTypeBasis, ArrayTypeWeights>( vals , data , basisVals , basisDVals , wts );
	break;
      }

  }

  template<class Scalar, class ArrayTypeOut, class ArrayTypeData,
	   class ArrayTypeBasis, class ArrayTypeWeights>
  void TensorProductSpaceTools::momentsGradCollocated( ArrayTypeOut &vals ,
					     const ArrayTypeData &data ,
					     const Array<RCP<ArrayTypeBasis> > &basisVals ,
					     const Array<RCP<ArrayTypeBasis> > &basisDVals ,
					     const Array<RCP<ArrayTypeWeights> > &wts )
  {
    const unsigned space_dim = basisVals.size();

#ifdef HAVE_INTREPID_DEBUG
    TEUCHOS_TEST_FOR_EXCEPTION( vals.rank() != 3 , std::invalid_argument ,
				">>> ERROR (TensorProductSpaceTools::evaluate): vals must be rank 3 array." );

    TEUCHOS_TEST_FOR_EXCEPTION( data.rank() != 4 , std::invalid_argument ,
				">>> ERROR (TensorProductSpaceTools::evaluate): coeffs must be rank 4 array." );

    TEUCHOS_TEST_FOR_EXCEPTION( basisVals.size() != (int)space_dim , std::invalid_argument ,
				">>> ERROR (TensorProductSpaceTools::evaluate): bases ill-sized." );

    TEUCHOS_TEST_FOR_EXCEPTION( basisDVals.size() != (int)space_dim , std::invalid_argument ,
				">>> ERROR (TensorProductSpaceTools::evaluate): bases ill-sized." );

    TEUCHOS_TEST_FOR_EXCEPTION( wts.size() != (int)space_dim , std::invalid_argument ,
				">>> ERROR (TensorProductSpaceTools::evaluate):  quadrature weights ill-sized." );
   
    for (unsigned d=0;d<space_dim;d++)
      {
	TEUCHOS_TEST_FOR_EXCEPTION( basisVals[d]->rank() != 2 , std::invalid_argument ,
				    ">>> ERROR (TensorProductSpaceTools::evaluate): each tabulated basis must be rank 2 array." );

	TEUCHOS_TEST_FOR_EXCEPTION( basisDVals[d]->rank() != 3 , std::invalid_argument ,
				    ">>> ERROR (TensorProductSpaceTools::evaluate): each tabulated derivative basis must be rank 3 array." );

	TEUCHOS_TEST_FOR_EXCEPTION( wts[d]->rank() != 1 , std::invalid_argument ,
				    ">>> ERROR (TensorProductSpaceTools::evaluate): each tabulated Dbasis must be rank 2 array." );
      }

    TEUCHOS_TEST_FOR_EXCEPTION( vals.dimension(0) != data.dimension(0) , std::invalid_argument,
				">>> ERROR (TensorProductSpaceTools::evaluate): Number of cells for vals and coeffs must match." );
    
    int product_of_basis_dimensions = 1;
    int product_of_basis_points = 1;
    for (unsigned d=0;d<space_dim;d++)
      {
	product_of_basis_dimensions *= basisVals[d]->dimension(0);
	product_of_basis_points *= basisVals[d]->dimension(1);
      }
    
    TEUCHOS_TEST_FOR_EXCEPTION( vals.dimension(2) != product_of_basis_dimensions , std::invalid_argument,
				">>> ERROR (TensorProductSpaceTools::evaluate): Incompatible number of points in vals and bases ." );
    
    TEUCHOS_TEST_FOR_EXCEPTION( data.dimension(2) != product_of_basis_points , std::invalid_argument,
				">>> ERROR (TensorProductSpaceTools::evaluate): Incompatible number of basis functions in data and bases ." );

#endif    

    switch (space_dim)
      {
      case 2:
	momentsGradCollocated2D<Scalar, ArrayTypeOut, ArrayTypeData, ArrayTypeBasis, ArrayTypeWeights>( vals , data , basisVals , basisDVals , wts );
	break;
      case 3:
	momentsGradCollocated3D<Scalar, ArrayTypeOut, ArrayTypeData, ArrayTypeBasis, ArrayTypeWeights>( vals , data , basisVals , basisDVals , wts );
	break;
      }

  }


  template<class Scalar, class ArrayTypeOut, class ArrayTypeCoeffs,
	   class ArrayTypeBasis>
  void TensorProductSpaceTools::evaluate2D( ArrayTypeOut &vals ,
					    const ArrayTypeCoeffs &coeffs ,
					    const Array<RCP<ArrayTypeBasis> > &bases )
  {
    const int numBfx = bases[0]->dimension(0);
    const int numBfy = bases[1]->dimension(0);
    
    const int numPtsx = bases[0]->dimension(1);
    const int numPtsy = bases[1]->dimension(1);
    
    const int numCells = vals.dimension(0);
    const int numFields = vals.dimension(1);

    const ArrayTypeBasis &Phix = *bases[0];
    const ArrayTypeBasis &Phiy = *bases[1];
    
    FieldContainer<double> Xi(numCells,numBfy,numPtsx);
    
    // sum factorization step

    for (int f=0;f<numFields;f++)
      {
	for (int cell=0;cell<numCells;cell++)
	  {
	    for (int j=0;j<numBfy;j++)
	      {
		for (int i=0;i<numBfx;i++)
		  {
		    const int I = j * numBfx + i;
		    for (int k=0;k<numPtsx;k++)
		      {
			Xi(cell,j,k) += coeffs(cell,f,I) * Phix(i,k);
		      }
		  }
	      }
	  }
	
	for (int cell=0;cell<numCells;cell++)
	  {
	    for (int kpty=0;kpty<numPtsy;kpty++)
	      {
		for (int kptx=0;kptx<numPtsx;kptx++)
		  {
		    vals(cell,f,kptx+numPtsx*kpty) = 0.0;
		  }
	      }
	    
	    // evaluation step
	    for (int kpty=0;kpty<numPtsy;kpty++)
	      {
		for (int kptx=0;kptx<numPtsx;kptx++)
		  {
		    const int I=kptx+numPtsx*kpty;
		    
		    for (int j=0;j<numBfy;j++)
		      {
			vals(cell,f,I) += Phiy(j,kpty) * Xi(cell,j,kptx);
		      }
		  }
	      }
	  }
      }
    
    return;
  }

  template<class Scalar, class ArrayTypeOut, class ArrayTypeCoeffs,
	   class ArrayTypeBasis>
  void TensorProductSpaceTools::evaluateCollocated2D( ArrayTypeOut &vals ,
						      const ArrayTypeCoeffs &coeffs ,
						      const Array<RCP<ArrayTypeBasis> > &bases )
  {
    // just copy coeffs to vals!

    const int numBfx = bases[0]->dimension(0);
    const int numBfy = bases[1]->dimension(0);
    
    const int numCells = vals.dimension(0);
    const int numFields = vals.dimension(1);

    for (int cell=0;cell<numCells;cell++)
      {
	for (int f=0;f<numFields;f++)
	  {
	    for (int j=0;j<numBfy;j++)
	      {
		for (int i=0;i<numBfx;i++)
		  {
		    const int I = j * numBfx + i;
		    vals(cell,f,I) = coeffs(cell,f,I);
		  }
	      }
	  }
      }

    return;
  }

  // template<class Scalar, class ArrayTypeOut, class ArrayTypeCoeffs,
  // 	   class ArrayTypeBasis>
  // void TensorProductSpaceTools::evaluateCollocated2D( ArrayTypeOut &vals ,
  // 						      const ArrayTypeCoeffs &coeffs ,
  // 						      const Array<Array<RCP<ArrayTypeBasis> > > &bases )
  // {
  //   // just copy coeffs to vals!
  //   const int numCells = vals.dimension(0);
  //   const int numFields = vals.dimension(1);
    
  //   const int numBfx = bases[comp][0]->dimension(0);
  //   const int numBfy = bases[comp][1]->dimension(0);
    
    
  //   for (int cell=0;cell<numCells;cell++)
  //     {
  // 	for (int f=0;f<numFields;f++)
  // 	  {
  // 	    for (int j=0;j<numBfy;j++)
  // 	      {
  // 		for (int i=0;i<numBfx;i++)
  // 		  {			const int I = j * numBfx + i;
  // 		    vals(cell,f,I,comp) = coeffs(cell,f,I);
  // 		  }
  // 	      }
  // 	  }
  //     }

  //   return;
  // }

  template<class Scalar, class ArrayTypeOut, class ArrayTypeCoeffs,
	   class ArrayTypeBasis>
  void TensorProductSpaceTools::evaluate3D( ArrayTypeOut &vals ,
					    const ArrayTypeCoeffs &coeffs ,
					    const Array<RCP<ArrayTypeBasis> > &bases )  

  {
    // checks to do:
    // vals point dimension is product of sizes of point arrays
    // vals
    
    const int numBfx = bases[0]->dimension(0);
    const int numBfy = bases[1]->dimension(0);
    const int numBfz = bases[2]->dimension(0);
    
    const int numPtsx = bases[0]->dimension(1);
    const int numPtsy = bases[1]->dimension(1);
    const int numPtsz = bases[2]->dimension(1);
    
    const int numCells = vals.dimension(0);
    const int numFields = vals.dimension(1);
    
    const ArrayTypeBasis &Phix = *bases[0];
    const ArrayTypeBasis &Phiy = *bases[1];
    const FieldContainer<double> &Phiz = *bases[2];

    FieldContainer<double> Xi(numCells, 
			      numBfz, numBfy , numPtsx);
    
    FieldContainer<double> Theta(numCells, 
				 numBfz , numPtsy, numPtsx );
    
    
    for (int f=0;f<numFields;f++)
      {

	// Xi section
	for (int c=0;c<numCells;c++)
	  {
	    for (int k=0;k<numBfz;k++)
	      {
		for (int j=0;j<numBfy;j++)
		  {
		    for (int l=0;l<numPtsx;l++)
		      {
			for (int i=0;i<numBfx;i++)
			  {
			    const int coeffIndex = k*numBfy*numBfx + j * numBfx + i;
			    Xi(c,k,j,l) += Phix(i,l) * coeffs(c,f,coeffIndex);
			  }
		      }
		  }
	      }
	  }
    
	// Theta section
	for (int c=0;c<numCells;c++)
	  {
	    for (int k=0;k<numBfz;k++)
	      {
		for (int m=0;m<numPtsy;m++)
		  {
		    for (int l=0;l<numPtsx;l++)
		      {
			for (int j=0;j<numBfy;j++)
			  {
			    Theta(c,k,m,l) += Phiy(j,m) * Xi(c,k,j,l);
			  }
		      }
		  }
	      }
	  }
    
	// final section
	for (int c=0;c<numCells;c++)
	  {
	    for (int n=0;n<numPtsz;n++)
	      {
		for (int m=0;m<numPtsy;m++)
		  {
		    for (int l=0;l<numPtsx;l++)
		      {
			vals(c,f,n*numPtsx*numPtsy+m*numPtsx+l) = 0.0;
			for (int k=0;k<numBfz;k++)
			  {
			    vals(c,f,n*numPtsx*numPtsy+m*numPtsx+l) += Theta(c,k,m,l) * Phiz(k,n);
			  }
		      }
		  }
	      }
	  }
      }

    return;
    
  }

  template<class Scalar, class ArrayTypeOut, class ArrayTypeCoeffs,
	   class ArrayTypeBasis>
  void TensorProductSpaceTools::evaluateCollocated3D( ArrayTypeOut &vals ,
						      const ArrayTypeCoeffs &coeffs ,
						      const Array<RCP<ArrayTypeBasis> > &bases )  

  {
    // copy coeffs to vals
    
    const int numBfx = bases[0]->dimension(0);
    const int numBfy = bases[1]->dimension(0);
    const int numBfz = bases[2]->dimension(0);
    
    const int numCells = vals.dimension(0);
    const int numFields = vals.dimension(1);
    
    for (int cell=0;cell<numCells;cell++)
      {
	for (int field=0;field<numFields;field++)
	  {	
	    for (int k=0;k<numBfz;k++)
	      {
		for (int j=0;j<numBfy;j++)
		  {
		    for (int i=0;i<numBfx;i++)
		      {
			const int I = k*numBfy*numBfx + j * numBfx + i;
			vals(cell,field,I) = coeffs(cell,field,I);
		      }
		  }
	      }
	  }
      }

    return;
    
  }


  template<class Scalar, class ArrayTypeOut, class ArrayTypeCoeffs,
	   class ArrayTypeBasis>
  void TensorProductSpaceTools::evaluateGradient2D( ArrayTypeOut &vals ,
						    const ArrayTypeCoeffs &coeffs ,
						    const Array<RCP<ArrayTypeBasis> > &bases ,
						    const Array<RCP<ArrayTypeBasis> > &Dbases )
  {
    
    const int numBfx = bases[0]->dimension(0);
    const int numBfy = bases[1]->dimension(0);
    const int numPtsx = bases[0]->dimension(1);
    const int numPtsy = bases[1]->dimension(1);
    const int numCells = vals.dimension(0);
    const int numFields = vals.dimension(1);
    const ArrayTypeBasis &Phix = *bases[0];
    const ArrayTypeBasis &Phiy = *bases[1];
    const ArrayTypeBasis &DPhix = *Dbases[0];
    const ArrayTypeBasis &DPhiy = *Dbases[1];
    
    FieldContainer<double> Xi(numBfy,numPtsx);

    for (int f=0;f<numFields;f++)
      {
    
	for (int cell=0;cell<numCells;cell++)
	  {
	    // x derivative
	
	    // sum factorization step
	    for (int j=0;j<numBfy;j++)
	      {
		for (int k=0;k<numPtsx;k++)
		  {
		    Xi(j,k) = 0.0;
		  }
	      }
	
	    for (int j=0;j<numBfy;j++)
	      {
		for (int i=0;i<numBfx;i++)
		  {
		    const int I = j * numBfx + i;
		    for (int k=0;k<numPtsx;k++)
		      {
			Xi(j,k) += coeffs(cell,f,I) * DPhix(i,k,0);
		      }
		  }
	      }
	
	    for (int kpty=0;kpty<numPtsy;kpty++)
	      {
		for (int kptx=0;kptx<numPtsx;kptx++)
		  {
		    vals(cell,f,kptx+numPtsx*kpty,0) = 0.0;
		  }
	      }
	
	    // evaluation step
	    for (int kpty=0;kpty<numPtsy;kpty++)
	      {
		for (int kptx=0;kptx<numPtsx;kptx++)
		  {
		    const int I=kptx+numPtsx*kpty;
		
		    for (int j=0;j<numBfy;j++)
		      {
			vals(cell,f,I,0) += Phiy(j,kpty) * Xi(j,kptx);
		      }
		  }
	      }
	
	    // y derivative
	
	    // sum factorization step
	    for (int j=0;j<numBfy;j++)
	      {
		for (int k=0;k<numPtsx;k++)
		  {
		    Xi(j,k) = 0.0;
		  }
	      }
	
	    for (int j=0;j<numBfy;j++)
	      {
		for (int i=0;i<numBfx;i++)
		  {
		    const int I = j * numBfx + i;
		    for (int k=0;k<numPtsx;k++)
		      {
			Xi(j,k) += coeffs(cell,f,I) * Phix(i,k);
		      }
		  }
	      }
	
	    for (int kpty=0;kpty<numPtsy;kpty++)
	      {
		for (int kptx=0;kptx<numPtsx;kptx++)
		  {
		    vals(cell,f,kptx+numPtsx*kpty,1) = 0.0;
		  }
	      }
	
	    // evaluation step
	    for (int kpty=0;kpty<numPtsy;kpty++)
	      {
		for (int kptx=0;kptx<numPtsx;kptx++)
		  {
		    const int I=kptx+numPtsx*kpty;
		
		    for (int j=0;j<numBfy;j++)
		      {
			vals(cell,f,I,1) += DPhiy(j,kpty,0) * Xi(j,kptx);
		      }
		  }
	      }
	  }
      }    
    return;
  }

  template<class Scalar, class ArrayTypeOut, class ArrayTypeCoeffs,
	   class ArrayTypeBasis>
  void TensorProductSpaceTools::evaluateGradientCollocated2D( ArrayTypeOut &vals ,
							      const ArrayTypeCoeffs &coeffs ,
							      const Array<RCP<ArrayTypeBasis> > &bases ,
							      const Array<RCP<ArrayTypeBasis> > &Dbases )
  {
    
    const int numBfx = bases[0]->dimension(0);
    const int numBfy = bases[1]->dimension(0);
    const int numPtsx = bases[0]->dimension(1);
    const int numPtsy = bases[1]->dimension(1);
    const int numCells = vals.dimension(0);
    const int numFields = vals.dimension(1);
    const ArrayTypeBasis &DPhix = *Dbases[0];
    const ArrayTypeBasis &DPhiy = *Dbases[1];
    
    for (int cell=0;cell<numCells;cell++)
      {
	for (int field=0;field<numFields;field++)
	  {
	    for (int j=0;j<numPtsy;j++)
	      {
		for (int i=0;i<numPtsx;i++)
		  {
		    const int I = j * numPtsx + i;
		    vals(cell,field,I,0) = 0.0;
		    vals(cell,field,I,1) = 0.0;
		  }
	      }
	  }
      }
    
    // x derivative
    for (int cell=0;cell<numCells;cell++)
      {
	for (int field=0;field<numFields;field++)
	  {
	    for (int j=0;j<numPtsy;j++)
	      {
		for (int i=0;i<numPtsx;i++)
		  {
		    const int I = j * numPtsx + i;
		    for (int ell=0;ell<numBfx;ell++)
		      {
			const int Itmp = j*numPtsx + ell;
			vals(cell,field,I,0) +=  coeffs(cell,field,Itmp) * DPhix( ell , i );
		      }
		  }
	      }
	  }
      }

    // y derivative
    for (int cell=0;cell<numCells;cell++)
      {
	for (int field=0;field<numFields;field++)
	  {
	    for (int j=0;j<numPtsy;j++)
	      {
		for (int i=0;i<numPtsx;i++)
		  {
		    const int I = j * numPtsx + i;
		    for (int m=0;m<numBfy;m++)
		      {
			const int Itmp = m*numPtsx + i;
			vals(cell,field,I,1) +=  coeffs(cell,field,Itmp) * DPhiy( m , j );
		      }
		  }
	      }
	  }
      }

  }

  template<class Scalar, class ArrayTypeOut, class ArrayTypeCoeffs,
	   class ArrayTypeBasis>
  void TensorProductSpaceTools::evaluateGradient3D( ArrayTypeOut &vals ,
						    const ArrayTypeCoeffs &coeffs ,
						    const Array<RCP<ArrayTypeBasis> > &bases ,
						    const Array<RCP<ArrayTypeBasis> > &Dbases )
  
  {
    // checks to do:
    // vals point dimension is product of sizes of point arrays
    // vals

    const int numBfx = bases[0]->dimension(0);
    const int numBfy = bases[1]->dimension(0);
    const int numBfz = bases[2]->dimension(0);
    const int numPtsx = bases[0]->dimension(1);
    const int numPtsy = bases[1]->dimension(1);
    const int numPtsz = bases[2]->dimension(1);
    const int numCells = vals.dimension(0);
    const int numFields = vals.dimension(1);
    const ArrayTypeBasis &Phix = *bases[0];
    const ArrayTypeBasis &Phiy = *bases[1];
    const ArrayTypeBasis &Phiz = *bases[2];
    const ArrayTypeBasis &DPhix = *Dbases[0];
    const ArrayTypeBasis &DPhiy = *Dbases[1];
    const ArrayTypeBasis &DPhiz = *Dbases[2];

    FieldContainer<double> Xi(numCells, 
			      numBfz, numBfy , numPtsx , 3) ;

    FieldContainer<double> Theta(numCells, 
				 numBfz , numPtsy, numPtsx , 3);
  
    for (int f=0;f<numFields;f++)
      {

	// Xi section
	for (int c=0;c<numCells;c++)
	  {
	    for (int k=0;k<numBfz;k++)
	      {
		for (int j=0;j<numBfy;j++)
		  {
		    for (int l=0;l<numPtsx;l++)
		      {
			for (int i=0;i<numBfx;i++)
			  {
			    const int coeffIndex = k*numBfz*numBfx + j * numBfx + i;
			    Xi(c,k,j,l,0) += DPhix(i,l,0) * coeffs(c,f,coeffIndex);
			    Xi(c,k,j,l,1) += Phix(i,l) * coeffs(c,f,coeffIndex);
			    Xi(c,k,j,l,2) += Phix(i,l) * coeffs(c,f,coeffIndex);
			  }
		      }
		  }
	      }
	  }

	// Theta section
	for (int c=0;c<numCells;c++)
	  {
	    for (int k=0;k<numBfz;k++)
	      {
		for (int m=0;m<numPtsy;m++)
		  {
		    for (int l=0;l<numPtsx;l++)
		      {
			for (int j=0;j<numBfy;j++)
			  {
			    Theta(c,k,m,l,0) += Phiy(j,m) * Xi(c,k,j,l,0);
			    Theta(c,k,m,l,1) += DPhiy(j,m,0) * Xi(c,k,j,l,1);
			    Theta(c,k,m,l,2) += Phiy(j,m) * Xi(c,k,j,l,2);
			  }
		      }
		  }
	      }
	  }

	// final section
	for (int c=0;c<numCells;c++)
	  {
	    for (int n=0;n<numPtsz;n++)
	      {
		for (int m=0;m<numPtsy;m++)
		  {
		    for (int l=0;l<numPtsx;l++)
		      {
			vals(c,f,n*numPtsx*numPtsy+m*numPtsx+l,0) = 0.0;
			vals(c,f,n*numPtsx*numPtsy+m*numPtsx+l,1) = 0.0;
			vals(c,f,n*numPtsx*numPtsy+m*numPtsx+l,2) = 0.0;
			for (int k=0;k<numBfz;k++)
			  {
			    vals(c,f,n*numPtsx*numPtsy+m*numPtsx+l,0) += Theta(c,k,m,l,0) * Phiz(k,n);
			    vals(c,f,n*numPtsx*numPtsy+m*numPtsx+l,1) += Theta(c,k,m,l,1) * Phiz(k,n);
			    vals(c,f,n*numPtsx*numPtsy+m*numPtsx+l,2) += Theta(c,k,m,l,2) * DPhiz(k,n,0);

			  }
		      }
		  }
	      }
	  }
      }
    return;
  }

  template<class Scalar, class ArrayTypeOut, class ArrayTypeCoeffs,
	   class ArrayTypeBasis>
  void TensorProductSpaceTools::evaluateGradientCollocated3D( ArrayTypeOut &vals ,
						    const ArrayTypeCoeffs &coeffs ,
						    const Array<RCP<ArrayTypeBasis> > &bases ,
						    const Array<RCP<ArrayTypeBasis> > &Dbases )
  
  {
    const int numBfx = bases[0]->dimension(0);
    const int numBfy = bases[1]->dimension(0);
    const int numBfz = bases[2]->dimension(0);
    const int numPtsx = bases[0]->dimension(1);
    const int numPtsy = bases[1]->dimension(1);
    const int numPtsz = bases[2]->dimension(1);
    const int numCells = vals.dimension(0);
    const int numFields = vals.dimension(1);
    const ArrayTypeBasis &Phix = *bases[0];
    const ArrayTypeBasis &Phiy = *bases[1];
    const ArrayTypeBasis &Phiz = *bases[2];
    const ArrayTypeBasis &DPhix = *Dbases[0];
    const ArrayTypeBasis &DPhiy = *Dbases[1];
    const ArrayTypeBasis &DPhiz = *Dbases[2];

    for (int cell=0;cell<numCells;cell++)
      {
	for (int field=0;field<numFields;field++)
	  {
	    for (int k=0;k<numPtsz;k++)
	      {
		for (int j=0;j<numPtsy;j++)
		  {
		    for (int i=0;i<numPtsx;i++)
		      {
			const int I = k * numPtsx * numPtsy + j * numPtsx + i;
			vals(cell,field,I,0) = 0.0;
			vals(cell,field,I,1) = 0.0;
			vals(cell,field,I,2) = 0.0;
		      }
		  }
	      }
	  }
      }
    
    // x derivative
    for (int cell=0;cell<numCells;cell++)
      {
	for (int field=0;field<numFields;field++)
	  {
	    for (int k=0;k<numPtsz;k++)
	      {
		for (int j=0;j<numPtsy;j++)
		  {
		    for (int i=0;i<numPtsx;i++)
		      {
			const int I = k * numPtsx * numPtsy + j * numPtsx + i;
			for (int ell=0;ell<numBfx;ell++)
			  {
			    const int Itmp = k * numPtsx * numPtsy + j*numPtsx + ell;
			    vals(cell,field,I,0) +=  coeffs(cell,field,Itmp) * DPhix( ell , i );
			  }
		      }
		  }
	      }
	  }
      }

    // y derivative
    for (int cell=0;cell<numCells;cell++)
      {
	for (int field=0;field<numFields;field++)
	  {
	    for (int k=0;k<numPtsz;k++)
	      {
		for (int j=0;j<numPtsy;j++)
		  {
		    for (int i=0;i<numPtsx;i++)
		      {
			const int I = k * numPtsx * numPtsy + j * numPtsx + i;
			for (int m=0;m<numBfy;m++)
			  {
			    const int Itmp = k * numPtsx * numPtsy + m * numPtsx + i;
			    vals(cell,field,I,1) +=  coeffs(cell,field,Itmp) * DPhiy( m , j );
			  }
		      }
		  }
	      }
	  }
      }


    // z derivative
    for (int cell=0;cell<numCells;cell++)
      {
	for (int field=0;field<numFields;field++)
	  {
	    for (int k=0;k<numPtsz;k++)
	      {
		for (int j=0;j<numPtsy;j++)
		  {
		    for (int i=0;i<numPtsx;i++)
		      {
			const int I = k * numPtsx * numPtsy + j * numPtsx + i;
			for (int n=0;n<numBfz;n++)
			  {
			    const int Itmp = n * numPtsx * numPtsy + j * numPtsx + i;
			    vals(cell,field,I,2) +=  coeffs(cell,field,Itmp) * DPhiz( n , k );
			  }
		      }
		  }
	      }
	  }
      }


  
    return;
  }

  template<class Scalar, class ArrayTypeOut, class ArrayTypeData,
	   class ArrayTypeBasis, class ArrayTypeWeights>
  void TensorProductSpaceTools::moments2D( ArrayTypeOut &vals ,
					   const ArrayTypeData &data ,
					   const Array<RCP<ArrayTypeBasis> > &basisVals ,
					   const Array<RCP<ArrayTypeWeights> > &wts )
  {
    const int numBfx = basisVals[0]->dimension(0);
    const int numBfy = basisVals[1]->dimension(0);
    const int numPtsx = basisVals[0]->dimension(1);
    const int numPtsy = basisVals[1]->dimension(1);
    const int numCells = vals.dimension(0);
    const int numFields = vals.dimension(1);
    const ArrayTypeBasis &Phix = *basisVals[0];
    const ArrayTypeBasis &Phiy = *basisVals[1];
    
    FieldContainer<double> Xi(numBfx,numPtsy);
    
    const ArrayTypeWeights &wtsx = *wts[0];
    const ArrayTypeWeights &wtsy = *wts[1];

    for (int f=0;f<numFields;f++)
      {
	for (int cell=0;cell<numCells;cell++)
	  {
	    // sum factorization step
	    for (int i=0;i<numBfx;i++)
	      {
		for (int k=0;k<numPtsy;k++)
		  {
		    Xi(i,k) = 0.0;
		  }
	      }
	
	    for (int i=0;i<numBfx;i++)
	      {
		for (int l=0;l<numPtsy;l++)
		  {
		    for (int k=0;k<numPtsx;k++)
		      {
			Xi(i,l) += wtsx(k) * data(cell,f,l*numPtsx+k) * Phix(i,k);
		      }
		  }
	      }
	
	    for (int j=0;j<numBfy;j++)
	      {
		for (int i=0;i<numBfx;i++)
		  {
		    vals(cell,f,j*numBfx+i) = 0.0;
		  }
	      }

	    // evaluate moments with sum factorization
	    for (int j=0;j<numBfy;j++)
	      {
		for (int i=0;i<numBfx;i++)
		  {
		    for (int l=0;l<numPtsy;l++)
		      {
			vals(cell,f,j*numBfx+i) += wtsy(l) * Phiy(j,l) * Xi(i,l);
		      }
		  }
	      }
	  }
      }
    return;

  }

  template<class Scalar, class ArrayTypeOut, class ArrayTypeData,
	   class ArrayTypeBasis, class ArrayTypeWeights>
  void TensorProductSpaceTools::momentsCollocated2D( ArrayTypeOut &vals ,
					   const ArrayTypeData &data ,
					   const Array<RCP<ArrayTypeBasis> > &basisVals ,
					   const Array<RCP<ArrayTypeWeights> > &wts )
  {
    const int numBfx = basisVals[0]->dimension(0);
    const int numBfy = basisVals[1]->dimension(0);
    const int numPtsx = basisVals[0]->dimension(1);
    const int numPtsy = basisVals[1]->dimension(1);
    const int numCells = vals.dimension(0);
    const int numFields = vals.dimension(1);
    
    FieldContainer<double> Xi(numBfx,numPtsy);
    
    const ArrayTypeWeights &wtsx = *wts[0];
    const ArrayTypeWeights &wtsy = *wts[1];

    for (int cell=0;cell<numCells;cell++)
      {
	for (int field=0;field<numFields;field++)
	  {
	    for (int i=0;i<numBfx;i++)
	      {
		for (int j=0;j<numBfy;j++)
		  {
		    const int I = numBfy * i + j;
		    vals(cell,field,I) = wtsx(i) * wtsy(j) * data(cell,field,I); 
		  }
	      }
	  }
      }

    return;

  }

  template<class Scalar, class ArrayTypeOut, class ArrayTypeData,
	   class ArrayTypeBasis, class ArrayTypeWeights>
  void TensorProductSpaceTools::moments3D( ArrayTypeOut &vals ,
					   const ArrayTypeData &data ,
					   const Array<RCP<ArrayTypeBasis> > &basisVals ,
					   const Array<RCP<ArrayTypeWeights> > &wts )
  {
    const int numBfx = basisVals[0]->dimension(0);
    const int numBfy = basisVals[1]->dimension(0);
    const int numBfz = basisVals[2]->dimension(0);

    const int numPtsx = basisVals[0]->dimension(1);
    const int numPtsy = basisVals[1]->dimension(1);
    const int numPtsz = basisVals[2]->dimension(1);

    const int numCells = vals.dimension(0);
    const int numFields = vals.dimension(1);

    const ArrayTypeBasis &Phix = *basisVals[0];
    const ArrayTypeBasis &Phiy = *basisVals[1];
    const ArrayTypeBasis &Phiz = *basisVals[2];

    const ArrayTypeWeights &Wtx = *wts[0];
    const ArrayTypeWeights &Wty = *wts[1];
    const ArrayTypeWeights &Wtz = *wts[2];

    FieldContainer<double> Xi(numCells,numBfx,numPtsz,numPtsy);
    FieldContainer<double> Theta(numCells,numBfy,numBfx,numPtsz);

    for (int f=0;f<numFields;f++)
      {
	// Xi phase
	for (int c=0;c<numCells;c++)
	  {
	    for (int i=0;i<numBfx;i++)
	      {
		for (int n=0;n<numPtsz;n++)
		  {
		    for (int m=0;m<numPtsy;m++)
		      {
			for (int l=0;l<numPtsx;l++)
			  {
			    Xi(c,i,n,m) += Wtx(l) * Phix(i,l) * data(c,f,n*numPtsy*numPtsx+m*numPtsx+l);
			  }
		      }
		  }
	      }
	  }

	// Theta phase
	for (int c=0;c<numCells;c++)
	  {
	    for (int j=0;j<numBfy;j++)
	      {
		for (int i=0;i<numBfx;i++)
		  {
		    for (int n=0;n<numPtsz;n++)
		      {
			for (int m=0;m<numPtsy;m++)
			  {
			    Theta(c,j,i,n) += Wty(m) * Phiy(j,m) * Xi(c,i,n,m);
			  }
		      }
		  }
	      }
	  }

	// Final phase
	for (int c=0;c<numCells;c++)
	  {
	    for (int k=0;k<numBfz;k++)
	      {
		for (int j=0;j<numBfy;j++)
		  {
		    for (int i=0;i<numBfx;i++)
		      {
			const int momIdx = k*numBfx*numBfy+j*numBfx+i;
			for (int n=0;n<numPtsz;n++)
			  {
			    vals(c,f,momIdx) += Wtz(n) * Phiz(k,n) * Theta(c,j,i,n);
			  }
		      }
		  }
	      }
	  }

      }
    return;

  }

  template<class Scalar, class ArrayTypeOut, class ArrayTypeData,
	   class ArrayTypeBasis, class ArrayTypeWeights>
  void TensorProductSpaceTools::momentsCollocated3D( ArrayTypeOut &vals ,
					   const ArrayTypeData &data ,
					   const Array<RCP<ArrayTypeBasis> > &basisVals ,
					   const Array<RCP<ArrayTypeWeights> > &wts )
  {
    const int numBfx = basisVals[0]->dimension(0);
    const int numBfy = basisVals[1]->dimension(0);
    const int numBfz = basisVals[2]->dimension(0);

    const int numPtsx = basisVals[0]->dimension(1);
    const int numPtsy = basisVals[1]->dimension(1);
    const int numPtsz = basisVals[2]->dimension(1);

    const int numCells = vals.dimension(0);
    const int numFields = vals.dimension(1);
    
    const ArrayTypeWeights &Wtx = *wts[0];
    const ArrayTypeWeights &Wty = *wts[1];
    const ArrayTypeWeights &Wtz = *wts[2];

    for (int cell=0;cell<numCells;cell++)
      {
	for (int field=0;field<numFields;field++)
	  {
	    for (int k=0;k<numBfz;k++)
	      {
		for (int j=0;j<numBfy;j++)
		  {
		    for (int i=0;i<numBfx;i++)
		      {
			const int I = k * numBfy * numBfx + j * numBfx + i;
			vals(cell,field,I) = Wtx(i) * Wty(j) * Wtz(k) * data(cell,field,I);
		      }
		  }
	      }
	  }
      }

    return;

  }

  // data is now (C,P,D)
  // want to compute the moments against the gradients of the basis 
  // functions.

  template<class Scalar, class ArrayTypeOut, class ArrayTypeData,
	   class ArrayTypeBasis, class ArrayTypeWeights>
  void TensorProductSpaceTools::momentsGrad2D( ArrayTypeOut &vals ,
					       const ArrayTypeData &data ,
					       const Array<RCP<ArrayTypeBasis> > &basisVals ,
					       const Array<RCP<ArrayTypeBasis> > &Dbases ,
					       const Array<RCP<ArrayTypeWeights> > &wts )
  {
    
    const int numBfx = basisVals[0]->dimension(0);
    const int numBfy = basisVals[1]->dimension(0);
    const int numPtsx = basisVals[0]->dimension(1);
    const int numPtsy = basisVals[1]->dimension(1);
    const int numCells = vals.dimension(0);
    const int numFields = vals.dimension(1);
    const ArrayTypeBasis &Phix = *basisVals[0];
    const ArrayTypeBasis &Phiy = *basisVals[1];
    const ArrayTypeBasis &DPhix = *Dbases[0];
    const ArrayTypeBasis &DPhiy = *Dbases[1];
    const ArrayTypeWeights &wtsx = *wts[0];
    const ArrayTypeWeights &wtsy = *wts[1];
    
    FieldContainer<double> Xi(numBfx,numPtsy,2);

    for (int f=0;f<numFields;f++)
      {
	// Xi phase
	for (int c=0;c<numCells;c++)
	  {
	    for (int i=0;i<numBfx;i++)
	      {
		for (int m=0;m<numPtsy;m++)
		  {
		    for (int l=0;l<numPtsx;l++)
		      {
			Xi(i,m,0) += wtsx(l) * data(c,f,m*numPtsy+l,0) * DPhix(i,l);
			Xi(i,m,1) += wtsx(l) * data(c,f,m*numPtsy+l,1) * Phix(i,l);
		      }
		  }
	      }
	  }

	// main phase
	for (int c=0;c<numCells;c++)
	  {
	    for (int j=0;j<numBfy;j++)
	      {
		for (int i=0;i<numBfx;i++)
		  {
		    const int bfIdx = j*numBfx+i;
		    vals(c,f,bfIdx) = 0.0;
		    for (int m=0;m<numPtsy;m++)
		      {
			vals(c,f,bfIdx) += wtsy(m) * Xi(i,m,0) * Phiy(j,m);
			vals(c,f,bfIdx) += wtsy(m) * Xi(i,m,1) * DPhiy(j,m);
		      }
		  }
	      }
	  }
      }
    return;

  }

  template<class Scalar, class ArrayTypeOut, class ArrayTypeData,
	   class ArrayTypeBasis, class ArrayTypeWeights>
  void TensorProductSpaceTools::momentsGradCollocated2D( ArrayTypeOut &vals ,
					       const ArrayTypeData &data ,
					       const Array<RCP<ArrayTypeBasis> > &basisVals ,
					       const Array<RCP<ArrayTypeBasis> > &Dbases ,
					       const Array<RCP<ArrayTypeWeights> > &wts )
  {
    
    const int numBfx = basisVals[0]->dimension(0);
    const int numBfy = basisVals[1]->dimension(0);
    const int numPtsx = basisVals[0]->dimension(1);
    const int numPtsy = basisVals[1]->dimension(1);
    const int numCells = vals.dimension(0);
    const int numFields = vals.dimension(1);
    const ArrayTypeBasis &Phix = *basisVals[0];
    const ArrayTypeBasis &Phiy = *basisVals[1];
    const ArrayTypeBasis &DPhix = *Dbases[0];
    const ArrayTypeBasis &DPhiy = *Dbases[1];
    const ArrayTypeWeights &wtsx = *wts[0];
    const ArrayTypeWeights &wtsy = *wts[1];

    for (int cell=0;cell<numCells;cell++)
      {
	for (int field=0;field<numFields;field++)
	  {
	    for (int j=0;j<numBfy;j++)
	      {
		for (int i=0;i<numBfx;i++)
		  {
		    const int I=j*numBfx+i;
		    vals(cell,field,I) = 0.0;
		  }
	      }
	  }
      }

    for (int cell=0;cell<numCells;cell++)
      {
	for (int field=0;field<numFields;field++)
	  {
	    for (int j=0;j<numBfy;j++)
	      {
		for (int i=0;i<numBfx;i++)
		  {
		    const int I=j*numBfx+i;
		    for (int m=0;m<numPtsx;m++)
		      {
			const int Itmp = j * numBfy + m;
			vals(cell,field,Itmp) += wtsx(m) * wtsy(j) * DPhix(i,m);
		      }
		  }
	      }
	  }
      }

    for (int cell=0;cell<numCells;cell++)
      {
	for (int field=0;field<numFields;field++)
	  {
	    for (int j=0;j<numBfy;j++)
	      {
		for (int i=0;i<numBfx;i++)
		  {
		    const int I=j*numBfx+i;
		    for (int n=0;n<numPtsy;n++)
		      {
			const int Itmp = n * numBfy + i;
			vals(cell,field,Itmp) += wtsx(i) * wtsy(n) * DPhiy(j,n);
		      }
		  }
	      }
	  }
      }
    
  }

  template<class Scalar, class ArrayTypeOut, class ArrayTypeData,
	   class ArrayTypeBasis, class ArrayTypeWeights>
  void TensorProductSpaceTools::momentsGrad3D( ArrayTypeOut &vals ,
					       const ArrayTypeData &data ,
					       const Array<RCP<ArrayTypeBasis> > &basisVals ,
					       const Array<RCP<ArrayTypeBasis> > &basisDVals ,
					       const Array<RCP<ArrayTypeWeights> > &wts )
  {

    const int numBfx = basisVals[0]->dimension(0);
    const int numBfy = basisVals[1]->dimension(0);
    const int numBfz = basisVals[2]->dimension(0);
    const int numPtsx = basisVals[0]->dimension(1);
    const int numPtsy = basisVals[1]->dimension(1);
    const int numPtsz = basisVals[2]->dimension(1);
    const int numCells = vals.dimension(0);
    const int numFields = vals.dimension(1);
    const ArrayTypeBasis &Phix = *basisVals[0];
    const ArrayTypeBasis &Phiy = *basisVals[1];
    const ArrayTypeBasis &Phiz = *basisVals[2];
    const ArrayTypeBasis &DPhix = *basisDVals[0];
    const ArrayTypeBasis &DPhiy = *basisDVals[1];
    const ArrayTypeBasis &DPhiz = *basisDVals[2];
    const ArrayTypeWeights &wtsx = *wts[0];
    const ArrayTypeWeights &wtsy = *wts[1];
    const ArrayTypeWeights &wtsz = *wts[2];

    FieldContainer<double> Xi(numCells,numBfx,numPtsz,numPtsy,3);
    FieldContainer<double> Theta(numCells,numBfy,numBfx,numPtsz,3);
  
    // Xi phase
    for (int f=0;f<numFields;f++)
      {
	for (int c=0;c<numCells;c++)
	  {
	    for (int i=0;i<numBfx;i++)
	      {
		for (int n=0;n<numPtsz;n++)
		  {
		    for (int m=0;m<numPtsy;m++)
		      {
			for (int l=0;l<numPtsx;l++)
			  {
			    const int dataIdx = n * numPtsy * numPtsx + m * numPtsx + l;
			    Xi(c,i,n,m,0) += wtsx(l) * DPhix(i,l,0) * data(c,f,dataIdx,0);
			    Xi(c,i,n,m,1) += wtsx(l) * Phix(i,l) * data(c,f,dataIdx,1);
			    Xi(c,i,n,m,2) += wtsx(l) * Phix(i,l) * data(c,f,dataIdx,2);
			  }
		      }
		  }
	      }
	  }

	// Theta phase
	for (int c=0;c<numCells;c++)
	  {
	    for (int j=0;j<numBfy;j++)
	      {
		for (int i=0;i<numBfx;i++)
		  {
		    for (int n=0;n<numPtsz;n++)
		      {
			for (int m=0;m<numPtsy;m++)
			  {
			    Theta(c,j,i,n,0) += wtsy(j) * Phiy(j,m) * Xi(c,i,n,m,0);
			    Theta(c,j,i,n,1) += wtsy(j) * DPhiy(j,m,0) * Xi(c,i,n,m,1);
			    Theta(c,j,i,n,2) += wtsy(j) * Phiy(j,m) * Xi(c,i,n,m,2);
			  }
		      }
		  }
	      }
	  }

	// last phase
	for (int c=0;c<numCells;c++)
	  {
	    for (int k=0;k<numBfz;k++)
	      {
		for (int j=0;j<numBfy;j++)
		  {
		    for (int i=0;i<numBfx;i++)
		      {
			const int momIdx = k * numBfx * numBfy + j * numBfx + i;
			for (int n=0;n<numPtsz;n++)
			  {
			    vals(c,f,momIdx) += wtsz(n) * Theta(c,j,i,n,0) * Phiz(k,n);
			    vals(c,f,momIdx) += wtsz(n) * Theta(c,j,i,n,1) * Phiz(k,n);
			    vals(c,f,momIdx) += wtsz(n) * Theta(c,j,i,n,2) * DPhiz(k,n,0);
			  }
		      }
		  }
	      }
	  }
      }

}

  template<class Scalar, class ArrayTypeOut, class ArrayTypeData,
	   class ArrayTypeBasis, class ArrayTypeWeights>
  void TensorProductSpaceTools::momentsGradCollocated3D( ArrayTypeOut &vals ,
					       const ArrayTypeData &data ,
					       const Array<RCP<ArrayTypeBasis> > &basisVals ,
					       const Array<RCP<ArrayTypeBasis> > &basisDVals ,
					       const Array<RCP<ArrayTypeWeights> > &wts )
  {

    const int numBfx = basisVals[0]->dimension(0);
    const int numBfy = basisVals[1]->dimension(0);
    const int numBfz = basisVals[2]->dimension(0);
    const int numPtsx = basisVals[0]->dimension(1);
    const int numPtsy = basisVals[1]->dimension(1);
    const int numPtsz = basisVals[2]->dimension(1);
    const int numCells = vals.dimension(0);
    const int numFields = vals.dimension(1);
    const ArrayTypeBasis &Phix = *basisVals[0];
    const ArrayTypeBasis &Phiy = *basisVals[1];
    const ArrayTypeBasis &Phiz = *basisVals[2];
    const ArrayTypeBasis &DPhix = *basisDVals[0];
    const ArrayTypeBasis &DPhiy = *basisDVals[1];
    const ArrayTypeBasis &DPhiz = *basisDVals[2];
    const ArrayTypeWeights &wtsx = *wts[0];
    const ArrayTypeWeights &wtsy = *wts[1];
    const ArrayTypeWeights &wtsz = *wts[2];

    for (int cell=0;cell<numCells;cell++)
      {
	for (int field=0;field<numFields;field++)
	  {
	    // x component of data versus x derivative of bases
	    for (int k=0;k<numBfz;k++)
	      {
		for (int j=0;j<numBfy;j++)
		  {
		    for (int i=0;i<numBfx;i++)
		      {
			const int I = numBfy * numBfx * k + numBfy * j + i;
			for (int ell=0;ell<numPtsx;ell++)
			  {
			    const int Itmp = numBfy * numBfx * k + numBfy * j + ell;
			    vals(cell,field,I) += wtsx(ell) * wtsy(j) * wtsz(k) * DPhix( i , ell );
			  }
		      }
		  }
	      }
	    // y component of data versus y derivative of bases
	    for (int k=0;k<numBfz;k++)
	      {
		for (int j=0;j<numBfy;j++)
		  {
		    for (int i=0;i<numBfx;i++)
		      {
			const int I = numBfy * numBfx * k + numBfy * j + i;
			for (int m=0;m<numPtsy;m++)
			  {
			    const int Itmp = numBfy * numBfx * k + numBfy * m + i;
			    vals(cell,field,I) += wtsx(i) * wtsy(m) * wtsz(k) * DPhiy( j , m );
			  }
		      }
		  }
	      }
	    // z component of data versus z derivative of bases
	    for (int k=0;k<numBfz;k++)
	      {
		for (int j=0;j<numBfy;j++)
		  {
		    for (int i=0;i<numBfx;i++)
		      {
			const int I = numBfy * numBfx * k + numBfy * j + i;
			for (int n=0;n<numPtsz;n++)
			  {
			    const int Itmp = numBfy * numBfx * n + numBfy * j + i;
			    vals(cell,field,I) += wtsx(i) * wtsy(j) * wtsz(n) * DPhiz( k , n );
			  }
		      }
		  }
	      }
	  }
      }
  
}



} // end namespace Intrepid
