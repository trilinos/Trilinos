// @HEADER
// *****************************************************************************
//                           Intrepid2 Package
//
// Copyright 2007 NTESS and the Intrepid2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER


/** \file test_01.hpp
    \brief  Unit tests for the Intrepid2::HGRAD_TET_C2_FEM class.
    \author Created by P. Bochev, D. Ridzal, K. Peterson and Kyungjoo Kim.
*/

#include "Intrepid2_config.h"

#ifdef HAVE_INTREPID2_DEBUG
#define INTREPID2_TEST_FOR_DEBUG_ABORT_OVERRIDE_TO_CONTINUE
#endif

#include "Intrepid2_HGRAD_TET_C2_FEM.hpp"

#include "Teuchos_oblackholestream.hpp"
#include "Teuchos_RCP.hpp"

#include "packages/intrepid2/unit-test/Discretization/Basis/Macros.hpp"
#include "packages/intrepid2/unit-test/Discretization/Basis/Setup.hpp"

namespace Intrepid2 {

namespace Test {

  using HostSpaceType = Kokkos::DefaultHostExecutionSpace;

  template<typename ValueType, typename DeviceType>
  int HGRAD_TET_C2_FEM_Test01(const bool verbose) {

    //! Create an execution space instance.
    const auto space = Kokkos::Experimental::partition_space(typename DeviceType::execution_space {}, 1)[0];

    //! Setup test output stream.
    Teuchos::RCP<std::ostream> outStream = setup_output_stream<DeviceType>(
      verbose, "Basis_HGRAD_TET_C2_FEM", {
        "1) Conversion of Dof tags into Dof ordinals and back",
        "2) Basis values for VALUE, GRAD, and Dk operators"
    });

    Teuchos::oblackholestream oldFormatState;
    oldFormatState.copyfmt(std::cout);

    typedef Kokkos::DynRankView<ValueType,DeviceType> DynRankView;
    typedef Kokkos::DynRankView<ValueType,HostSpaceType>   DynRankViewHost;

    const ValueType tol = tolerence();
    int errorFlag = 0;
      
    // for virtual function, value and point types are declared in the class
    typedef ValueType outputValueType;
    typedef ValueType pointValueType;
    Basis_HGRAD_TET_C2_FEM<DeviceType,outputValueType,pointValueType> tetBasis;

  *outStream
    << "\n"
    << "===============================================================================\n"
    << "| TEST 1: Basis creation, exception testing                                   |\n"
    << "===============================================================================\n";
  
  try{
    ordinal_type nthrow = 0, ncatch = 0;
#ifdef HAVE_INTREPID2_DEBUG 

    DynRankView ConstructWithLabel(tetNodes, 10, 3);

    const ordinal_type numFields = tetBasis.getCardinality();
    const ordinal_type numPoints = tetNodes.extent(0);
    const ordinal_type spaceDim  = tetBasis.getBaseCellTopology().getDimension();

    DynRankView ConstructWithLabel(vals, numFields, numPoints);
    DynRankView ConstructWithLabel(vals_vec, numFields, numPoints, 4);
    {
    // exception #1: CURL cannot be applied to scalar functions
    // resize vals to rank-3 container with dimensions (num. points, num. basis functions, arbitrary)
      INTREPID2_TEST_ERROR_EXPECTED( tetBasis.getValues(vals_vec, tetNodes, OPERATOR_CURL) );
    }
    {
    // exception #2: DIV cannot be applied to scalar functions
    // resize vals to rank-2 container with dimensions (num. points, num. basis functions)
      INTREPID2_TEST_ERROR_EXPECTED( tetBasis.getValues(vals, tetNodes, OPERATOR_DIV) );
    }
        
    // Exceptions 3-7: all bf tags/bf Ids below are wrong and should cause getDofOrdinal() and 
    // getDofTag() to access invalid array elements thereby causing bounds check exception
    {
    // exception #3
      INTREPID2_TEST_ERROR_EXPECTED( tetBasis.getDofOrdinal(3,0,0) );
    // exception #4
      INTREPID2_TEST_ERROR_EXPECTED( tetBasis.getDofOrdinal(1,1,1) );
    // exception #5
      INTREPID2_TEST_ERROR_EXPECTED( tetBasis.getDofOrdinal(0,4,0) );
    // exception #6
      INTREPID2_TEST_ERROR_EXPECTED( tetBasis.getDofTag(numFields) );
    // exception #7
      INTREPID2_TEST_ERROR_EXPECTED( tetBasis.getDofTag(-1) );
    }
    
    // Exceptions 8-18 test exception handling with incorrectly dimensioned input/output arrays
    {
    // exception #8: input points array must be of rank-2
      DynRankView ConstructWithLabel(badPoints1, 4, 5, 3);
      INTREPID2_TEST_ERROR_EXPECTED( tetBasis.getValues(vals, badPoints1, OPERATOR_VALUE) );
    }
    { 
    // exception #9 dimension 1 in the input point array must equal space dimension of the cell
      DynRankView ConstructWithLabel(badPoints2, 4, spaceDim - 1);
      INTREPID2_TEST_ERROR_EXPECTED( tetBasis.getValues(vals, badPoints2, OPERATOR_VALUE) );
    }
    { 
    // exception #10 output values must be of rank-2 for OPERATOR_VALUE
      DynRankView ConstructWithLabel(badVals1, 4, 3, 1);
      INTREPID2_TEST_ERROR_EXPECTED( tetBasis.getValues(badVals1, tetNodes, OPERATOR_VALUE) );
    }
    { 
    // exception #11 output values must be of rank-3 for OPERATOR_GRAD
      DynRankView ConstructWithLabel(badVals2, 4, 3);
      INTREPID2_TEST_ERROR_EXPECTED( tetBasis.getValues(badVals2, tetNodes, OPERATOR_GRAD) );
    // exception #12 output values must be of rank-3 for OPERATOR_D1
      INTREPID2_TEST_ERROR_EXPECTED( tetBasis.getValues(badVals2, tetNodes, OPERATOR_D1) );
    // exception #13 output values must be of rank-3 for OPERATOR_D2
      INTREPID2_TEST_ERROR_EXPECTED( tetBasis.getValues(badVals2, tetNodes, OPERATOR_D2) );
    }
    {
    // exception #14 incorrect 0th dimension of output array (must equal number of basis functions)
      DynRankView ConstructWithLabel(badVals3, numFields + 1, numPoints);
      INTREPID2_TEST_ERROR_EXPECTED( tetBasis.getValues(badVals3, tetNodes, OPERATOR_VALUE) );
    }
    { 
    // exception #15 incorrect 1st dimension of output array (must equal number of points)
      DynRankView ConstructWithLabel(badVals4, numFields, numPoints + 1);
      INTREPID2_TEST_ERROR_EXPECTED( tetBasis.getValues(badVals4, tetNodes, OPERATOR_VALUE) );
    }
    { 
    // exception #16: incorrect 2nd dimension of output array (must equal the space dimension)
      DynRankView ConstructWithLabel(badVals5, numFields, numPoints, spaceDim + 1);
      INTREPID2_TEST_ERROR_EXPECTED( tetBasis.getValues(badVals5, tetNodes, OPERATOR_GRAD) );
    }
    { 
    // exception #17: incorrect 2nd dimension of output array (must equal D2 cardinality in 2D)
      DynRankView ConstructWithLabel(badVals6, numFields, numPoints, 40);
      INTREPID2_TEST_ERROR_EXPECTED( tetBasis.getValues(badVals6, tetNodes, OPERATOR_D1) );
    // exception #18: incorrect 2nd dimension of output array (must equal D3 cardinality in 2D)
      INTREPID2_TEST_ERROR_EXPECTED( tetBasis.getValues(badVals6, tetNodes, OPERATOR_D2) );
    }
#endif
    // Check if number of thrown exceptions matches the one we expect 
    if (nthrow != ncatch) {
      errorFlag++;
      *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
      *outStream << "# of catch ("<< ncatch << ") is different from # of throw (" << nthrow << ")\n";
    }
  } catch (std::logic_error &err) {
    *outStream << "UNEXPECTED ERROR !!! ----------------------------------------------------------\n";
    *outStream << err.what() << '\n';
    *outStream << "-------------------------------------------------------------------------------" << "\n\n";
    errorFlag = -1000;
  };
  
  *outStream
    << "\n"
    << "===============================================================================\n"
    << "| TEST 2: correctness of tag to enum and enum to tag lookups                  |\n"
    << "===============================================================================\n";
  
  try{
    const ordinal_type numFields = tetBasis.getCardinality();
    const auto allTags = tetBasis.getAllDofTags();
    
    // Loop over all tags, lookup the associated dof enumeration and then lookup the tag again
    const ordinal_type dofTagSize = allTags.extent(0);
    for (ordinal_type i = 0; i < dofTagSize; ++i) {
      const auto bfOrd  = tetBasis.getDofOrdinal(allTags(i,0), allTags(i,1), allTags(i,2));
      
      const auto myTag = tetBasis.getDofTag(bfOrd);
       if( !( (myTag(0) == allTags(i,0)) &&
              (myTag(1) == allTags(i,1)) &&
              (myTag(2) == allTags(i,2)) &&
              (myTag(3) == allTags(i,3)) ) ) {
        errorFlag++;
        *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
        *outStream << " getDofOrdinal( {" 
          << allTags(i,0) << ", " 
          << allTags(i,1) << ", " 
          << allTags(i,2) << ", " 
          << allTags(i,3) << "}) = " << bfOrd <<" but \n";   
        *outStream << " getDofTag(" << bfOrd << ") = { "
          << myTag(0) << ", " 
          << myTag(1) << ", " 
          << myTag(2) << ", " 
          << myTag(3) << "}\n";        
      }
    }
    
    // Now do the same but loop over basis functions
    for( ordinal_type bfOrd = 0; bfOrd < numFields; bfOrd++) {
      const auto myTag  = tetBasis.getDofTag(bfOrd);
      const auto myBfOrd = tetBasis.getDofOrdinal(myTag(0), myTag(1), myTag(2));
      if( bfOrd != myBfOrd) {
        errorFlag++;
        *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
        *outStream << " getDofTag(" << bfOrd << ") = { "
          << myTag(0) << ", " 
          << myTag(1) << ", " 
          << myTag(2) << ", " 
          << myTag(3) << "} but getDofOrdinal({" 
          << myTag(0) << ", " 
          << myTag(1) << ", " 
          << myTag(2) << ", " 
          << myTag(3) << "} ) = " << myBfOrd << "\n";
      }
    }
  }
  catch (std::logic_error &err){
    *outStream << err.what() << "\n\n";
    errorFlag = -1000;
  };
  
  *outStream
    << "\n"
    << "===============================================================================\n"
    << "| TEST 3: correctness of basis function values                                |\n"
    << "===============================================================================\n";
  
  outStream -> precision(20);
  
  // VALUE: in (F,P) format
  const ValueType basisValues[] = {
    1.00000, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.00000, 0, 0, 0, 0, 0, 0, 0, \
    0, 0, 0, 1.00000, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.00000, 0, 0, 0, 0, \
    0, 0, 0, 0, 0, 0, 1.00000, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.00000, 0, \
    0, 0, 0, 0, 0, 0, 0, 0, 0, 1.00000, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, \
    1.00000, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.00000, 0, 0, 0, 0, 0, 0, 0, \
    0, 0, 0, 1.00000 };
  
  // GRAD and D1: in (F,P,D) format
  const ValueType basisGrads[] = {
    -3.00000, -3.00000, -3.00000, 1.00000, 1.00000, 1.00000, 1.00000, \
    1.00000, 1.00000, 1.00000, 1.00000, 1.00000, -1.00000, -1.00000, \
    -1.00000, 1.00000, 1.00000, 1.00000, -1.00000, -1.00000, -1.00000, \
    -1.00000, -1.00000, -1.00000, 1.00000, 1.00000, 1.00000, 1.00000, \
    1.00000, 1.00000, -1.00000, 0, 0, 3.00000, 0, 0, -1.00000, 0, 0, \
    -1.00000, 0, 0, 1.00000, 0, 0, 1.00000, 0, 0, -1.00000, 0, 0, \
    -1.00000, 0, 0, 1.00000, 0, 0, -1.00000, 0, 0, 0, -1.00000, 0, 0, \
    -1.00000, 0, 0, 3.00000, 0, 0, -1.00000, 0, 0, -1.00000, 0, 0, \
    1.00000, 0, 0, 1.00000, 0, 0, -1.00000, 0, 0, -1.00000, 0, 0, \
    1.00000, 0, 0, 0, -1.00000, 0, 0, -1.00000, 0, 0, -1.00000, 0, 0, \
    3.00000, 0, 0, -1.00000, 0, 0, -1.00000, 0, 0, -1.00000, 0, 0, \
    1.00000, 0, 0, 1.00000, 0, 0, 1.00000, 4.00000, 0, 0, -4.00000, \
    -4.00000, -4.00000, 0, 0, 0, 0, 0, 0, 0, -2.00000, -2.00000, \
    -2.00000, -2.00000, -2.00000, 2.00000, 0, 0, 2.00000, 0, 0, -2.00000, \
    -2.00000, -2.00000, 0, 0, 0, 0, 0, 0, 0, 4.00000, 0, 4.00000, 0, 0, \
    0, 0, 0, 0, 2.00000, 0, 2.00000, 2.00000, 0, 2.00000, 0, 0, 0, 0, 0, \
    0, 2.00000, 0, 2.00000, 0, 0, 0, 4.00000, 0, 0, 0, 0, -4.00000, \
    -4.00000, -4.00000, 0, 0, 0, 0, 2.00000, 0, -2.00000, -2.00000, \
    -2.00000, -2.00000, 0, -2.00000, 0, 2.00000, 0, 0, 0, 0, -2.00000, \
    -2.00000, -2.00000, 0, 0, 4.00000, 0, 0, 0, 0, 0, 0, -4.00000, \
    -4.00000, -4.00000, 0, 0, 2.00000, 0, 0, 0, 0, 0, 2.00000, -2.00000, \
    -2.00000, 0, -2.00000, -2.00000, -2.00000, -2.00000, -2.00000, \
    -2.00000, 0, 0, 0, 0, 0, 4.00000, 0, 0, 0, 4.00000, 0, 0, 0, 0, \
    2.00000, 0, 0, 2.00000, 0, 0, 0, 2.00000, 0, 0, 2.00000, 0, 2.00000, \
    2.00000, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 4.00000, 0, 4.00000, 0, 0, 0, \
    0, 0, 0, 2.00000, 0, 0, 2.00000, 0, 2.00000, 0, 0, 2.00000, 0, 0, \
    2.00000, 2.00000};
  
  // D2 values in (F,P, Dk) format
  const ValueType basisD2[]={
    4.00000, 4.00000, 4.00000, 4.00000, 4.00000, 4.00000, 4.00000, \
    4.00000, 4.00000, 4.00000, 4.00000, 4.00000, 4.00000, 4.00000, \
    4.00000, 4.00000, 4.00000, 4.00000, 4.00000, 4.00000, 4.00000, \
    4.00000, 4.00000, 4.00000, 4.00000, 4.00000, 4.00000, 4.00000, \
    4.00000, 4.00000, 4.00000, 4.00000, 4.00000, 4.00000, 4.00000, \
    4.00000, 4.00000, 4.00000, 4.00000, 4.00000, 4.00000, 4.00000, \
    4.00000, 4.00000, 4.00000, 4.00000, 4.00000, 4.00000, 4.00000, \
    4.00000, 4.00000, 4.00000, 4.00000, 4.00000, 4.00000, 4.00000, \
    4.00000, 4.00000, 4.00000, 4.00000, 4.00000, 0, 0, 0, 0, 0, 4.00000, \
    0, 0, 0, 0, 0, 4.00000, 0, 0, 0, 0, 0, 4.00000, 0, 0, 0, 0, 0, \
    4.00000, 0, 0, 0, 0, 0, 4.00000, 0, 0, 0, 0, 0, 4.00000, 0, 0, 0, 0, \
    0, 4.00000, 0, 0, 0, 0, 0, 4.00000, 0, 0, 0, 0, 0, 4.00000, 0, 0, 0, \
    0, 0, 0, 0, 0, 4.00000, 0, 0, 0, 0, 0, 4.00000, 0, 0, 0, 0, 0, \
    4.00000, 0, 0, 0, 0, 0, 4.00000, 0, 0, 0, 0, 0, 4.00000, 0, 0, 0, 0, \
    0, 4.00000, 0, 0, 0, 0, 0, 4.00000, 0, 0, 0, 0, 0, 4.00000, 0, 0, 0, \
    0, 0, 4.00000, 0, 0, 0, 0, 0, 4.00000, 0, 0, 0, 0, 0, 0, 0, 4.00000, \
    0, 0, 0, 0, 0, 4.00000, 0, 0, 0, 0, 0, 4.00000, 0, 0, 0, 0, 0, \
    4.00000, 0, 0, 0, 0, 0, 4.00000, 0, 0, 0, 0, 0, 4.00000, 0, 0, 0, 0, \
    0, 4.00000, 0, 0, 0, 0, 0, 4.00000, 0, 0, 0, 0, 0, 4.00000, 0, 0, 0, \
    0, 0, 4.00000, -8.00000, -4.00000, -4.00000, 0, 0, 0, -8.00000, \
    -4.00000, -4.00000, 0, 0, 0, -8.00000, -4.00000, -4.00000, 0, 0, 0, \
    -8.00000, -4.00000, -4.00000, 0, 0, 0, -8.00000, -4.00000, -4.00000, \
    0, 0, 0, -8.00000, -4.00000, -4.00000, 0, 0, 0, -8.00000, -4.00000, \
    -4.00000, 0, 0, 0, -8.00000, -4.00000, -4.00000, 0, 0, 0, -8.00000, \
    -4.00000, -4.00000, 0, 0, 0, -8.00000, -4.00000, -4.00000, 0, 0, 0, \
    0, 4.00000, 0, 0, 0, 0, 0, 4.00000, 0, 0, 0, 0, 0, 4.00000, 0, 0, 0, \
    0, 0, 4.00000, 0, 0, 0, 0, 0, 4.00000, 0, 0, 0, 0, 0, 4.00000, 0, 0, \
    0, 0, 0, 4.00000, 0, 0, 0, 0, 0, 4.00000, 0, 0, 0, 0, 0, 4.00000, 0, \
    0, 0, 0, 0, 4.00000, 0, 0, 0, 0, 0, -4.00000, 0, -8.00000, -4.00000, \
    0, 0, -4.00000, 0, -8.00000, -4.00000, 0, 0, -4.00000, 0, -8.00000, \
    -4.00000, 0, 0, -4.00000, 0, -8.00000, -4.00000, 0, 0, -4.00000, 0, \
    -8.00000, -4.00000, 0, 0, -4.00000, 0, -8.00000, -4.00000, 0, 0, \
    -4.00000, 0, -8.00000, -4.00000, 0, 0, -4.00000, 0, -8.00000, \
    -4.00000, 0, 0, -4.00000, 0, -8.00000, -4.00000, 0, 0, -4.00000, 0, \
    -8.00000, -4.00000, 0, 0, 0, -4.00000, 0, -4.00000, -8.00000, 0, 0, \
    -4.00000, 0, -4.00000, -8.00000, 0, 0, -4.00000, 0, -4.00000, \
    -8.00000, 0, 0, -4.00000, 0, -4.00000, -8.00000, 0, 0, -4.00000, 0, \
    -4.00000, -8.00000, 0, 0, -4.00000, 0, -4.00000, -8.00000, 0, 0, \
    -4.00000, 0, -4.00000, -8.00000, 0, 0, -4.00000, 0, -4.00000, \
    -8.00000, 0, 0, -4.00000, 0, -4.00000, -8.00000, 0, 0, -4.00000, 0, \
    -4.00000, -8.00000, 0, 0, 4.00000, 0, 0, 0, 0, 0, 4.00000, 0, 0, 0, \
    0, 0, 4.00000, 0, 0, 0, 0, 0, 4.00000, 0, 0, 0, 0, 0, 4.00000, 0, 0, \
    0, 0, 0, 4.00000, 0, 0, 0, 0, 0, 4.00000, 0, 0, 0, 0, 0, 4.00000, 0, \
    0, 0, 0, 0, 4.00000, 0, 0, 0, 0, 0, 4.00000, 0, 0, 0, 0, 0, 0, 0, \
    4.00000, 0, 0, 0, 0, 0, 4.00000, 0, 0, 0, 0, 0, 4.00000, 0, 0, 0, 0, \
    0, 4.00000, 0, 0, 0, 0, 0, 4.00000, 0, 0, 0, 0, 0, 4.00000, 0, 0, 0, \
    0, 0, 4.00000, 0, 0, 0, 0, 0, 4.00000, 0, 0, 0, 0, 0, 4.00000, 0, 0, \
    0, 0, 0, 4.00000, 0
  };
  
  try{
    DynRankViewHost ConstructWithLabel(tetNodesHost, 10, 3);
    tetNodesHost(0,0) =  0.0;  tetNodesHost(0,1) =  0.0;  tetNodesHost(0,2) =  0.0;  
    tetNodesHost(1,0) =  1.0;  tetNodesHost(1,1) =  0.0;  tetNodesHost(1,2) =  0.0;  
    tetNodesHost(2,0) =  0.0;  tetNodesHost(2,1) =  1.0;  tetNodesHost(2,2) =  0.0;
    tetNodesHost(3,0) =  0.0;  tetNodesHost(3,1) =  0.0;  tetNodesHost(3,2) =  1.0;  
    
    tetNodesHost(4,0) =  0.5;  tetNodesHost(4,1) =  0.0;  tetNodesHost(4,2) =  0.0;
    tetNodesHost(5,0) =  0.5;  tetNodesHost(5,1) =  0.5;  tetNodesHost(5,2) =  0.0;  
    tetNodesHost(6,0) =  0.0;  tetNodesHost(6,1) =  0.5;  tetNodesHost(6,2) =  0.0;  
    tetNodesHost(7,0) =  0.0;  tetNodesHost(7,1) =  0.0;  tetNodesHost(7,2) =  0.5;  
    tetNodesHost(8,0) =  0.5;  tetNodesHost(8,1) =  0.0;  tetNodesHost(8,2) =  0.5;  
    tetNodesHost(9,0) =  0.0;  tetNodesHost(9,1) =  0.5;  tetNodesHost(9,2) =  0.5;  

    auto tetNodes = Kokkos::create_mirror_view(typename DeviceType::memory_space(), tetNodesHost);
    Kokkos::deep_copy(tetNodes, tetNodesHost);
        
    // Dimensions for the output arrays:
    const ordinal_type numFields = tetBasis.getCardinality();
    const ordinal_type numPoints = tetNodes.extent(0);
    const ordinal_type spaceDim  = tetBasis.getBaseCellTopology().getDimension();
    const ordinal_type D2cardinality = getDkCardinality(OPERATOR_D2, spaceDim);
    
    {
    // Check VALUE of basis functions: resize vals to rank-2 container:
    DynRankView ConstructWithLabel(vals, numFields, numPoints);
    tetBasis.getValues(space, vals, tetNodes, OPERATOR_VALUE);
    auto vals_host = Kokkos::create_mirror_view(typename HostSpaceType::memory_space(), vals);
    Kokkos::deep_copy(vals_host, vals);
    for (ordinal_type i = 0; i < numFields; ++i) {
      for (ordinal_type j = 0; j < numPoints; ++j) {
          const ordinal_type l =  i + j * numFields;
           if (std::abs(vals_host(i,j) - basisValues[l]) > tol) {
             errorFlag++;
             *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";

             // Output the multi-index of the value where the error is:
             *outStream << " At multi-index { ";
             *outStream << i << " ";*outStream << j << " ";
             *outStream << "}  computed value: " << vals_host(i,j)
               << " but reference value: " << basisValues[l] << "\n";
         }
      }
    }
    }

    {
    // Check GRAD of basis function: resize vals to rank-3 container
    DynRankView ConstructWithLabel(vals, numFields, numPoints, spaceDim);
    tetBasis.getValues(space, vals, tetNodes, OPERATOR_GRAD);
    auto vals_host = Kokkos::create_mirror_view(typename HostSpaceType::memory_space(), vals);
    Kokkos::deep_copy(vals_host, vals);
    for (ordinal_type i = 0; i < numFields; ++i) {
      for (ordinal_type j = 0; j < numPoints; ++j) {
        for (ordinal_type k = 0; k < spaceDim; ++k) {
 
          // basisGrads is (F,P,D), compute offset:
          const ordinal_type l = k + j * spaceDim + i * spaceDim * numPoints;
           if (std::abs(vals_host(i,j,k) - basisGrads[l]) > tol) {
             errorFlag++;
             *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";

             // Output the multi-index of the value where the error is:
             *outStream << " At multi-index { ";
             *outStream << i << " ";*outStream << j << " ";*outStream << k << " ";
             *outStream << "}  computed grad component: " << vals_host(i,j,k)
               << " but reference grad component: " << basisGrads[l] << "\n";
            }
         }
      }
    }

    // Check D1 of basis function (do not resize vals because it has the correct size: D1 = GRAD)
    tetBasis.getValues(space, vals, tetNodes, OPERATOR_D1);
    Kokkos::deep_copy(vals_host, vals);
    for (ordinal_type i = 0; i < numFields; ++i) {
      for (ordinal_type j = 0; j < numPoints; ++j) {
        for (ordinal_type k = 0; k < spaceDim; ++k) {
          
          // basisGrads is (F,P,D), compute offset:
          const ordinal_type l = k + j * spaceDim + i * spaceDim * numPoints;
           if (std::abs(vals_host(i,j,k) - basisGrads[l]) > tol) {
             errorFlag++;
             *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";

             // Output the multi-index of the value where the error is:
             *outStream << " At multi-index { ";
             *outStream << i << " ";*outStream << j << " ";*outStream << k << " ";
             *outStream << "}  computed D1 component: " << vals_host(i,j,k)
               << " but reference D1 component: " << basisGrads[l] << "\n";
            }
         }
      }
    }
    }

    {
    // Check D2 of basis function
    DynRankView ConstructWithLabel(vals, numFields, numPoints, D2cardinality);
    tetBasis.getValues(space, vals, tetNodes, OPERATOR_D2);
    auto vals_host = Kokkos::create_mirror_view(typename HostSpaceType::memory_space(), vals);
    Kokkos::deep_copy(vals_host, vals);
    for (ordinal_type i = 0; i < numFields; ++i) {
      for (ordinal_type j = 0; j < numPoints; ++j) {
        for (ordinal_type k = 0; k < D2cardinality; ++k) {
          
          // basisD2 is (F,P,Dk), compute offset:
          const ordinal_type l = k + j * D2cardinality + i * D2cardinality * numPoints;
          if (std::abs(vals_host(i,j,k) - basisD2[l]) > tol) {
            errorFlag++;
            *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
            
            // Output the multi-index of the value where the error is:
            *outStream << " At multi-index { ";
            *outStream << i << " ";*outStream << j << " ";*outStream << k << " ";
            *outStream << "}  computed D2 component: " << vals_host(i,j,k)
              << " but reference D2 component: " << basisD2[l] << "\n";
          }
        }
      }
    }
    }
    
    {
    // Check all higher derivatives - must be zero. 

      const EOperator ops[] = { OPERATOR_D3,
                                OPERATOR_D4,
                                OPERATOR_D5,
                                OPERATOR_D6,
                                OPERATOR_D7,
                                OPERATOR_D8,
                                OPERATOR_D9,
                                OPERATOR_D10,
                                OPERATOR_MAX };
      for (auto h=0;ops[h]!=OPERATOR_MAX;++h) {
        const auto op = ops[h];
      // The last dimension is the number of kth derivatives and needs to be resized for every Dk
        const ordinal_type DkCardin  = getDkCardinality(op, spaceDim);
        DynRankView vals("vals", numFields, numPoints, DkCardin);

        tetBasis.getValues(space, vals, tetNodes, op);
        auto vals_host = Kokkos::create_mirror_view(typename HostSpaceType::memory_space(), vals);
        Kokkos::deep_copy(vals_host, vals);

        for (ordinal_type i1 = 0; i1 < numFields; ++i1)
          for (ordinal_type i2 = 0; i2 < numPoints; ++i2)
            for (ordinal_type i3 = 0; i3 < DkCardin; ++i3) {
              if (std::abs(vals_host(i1,i2,i3)) > tol) {
                errorFlag++;
                *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
            
                // Get the multi-index of the value where the error is and the operator order
                int ord = Intrepid2::getOperatorOrder(op);
                *outStream << " At multi-index { "<<i1<<" "<<i2 <<" "<<i3;
                *outStream << "}  computed D"<< ord <<" component: " << vals_host(i1,i2,i3) 
                           << " but reference D" << ord << " component:  0 \n";
              }
            }
      }    
    }    
  } catch (std::logic_error &err) {
    *outStream << err.what() << "\n\n";
    errorFlag = -1000;
  };
    
  *outStream
  << "\n"
  << "===============================================================================\n"
  << "| TEST 4: Function Space is Correct                                           |\n"
  << "===============================================================================\n";
  
  try {
    const EFunctionSpace fs = tetBasis.getFunctionSpace();
    
    if (fs != FUNCTION_SPACE_HGRAD)
    {
      *outStream << std::setw(70) << "------------- TEST FAILURE! -------------" << "\n";
      
      // Output the multi-index of the value where the error is:
      *outStream << " Expected a function space of FUNCTION_SPACE_HGRAD (enum value " << FUNCTION_SPACE_HGRAD << "),";
      *outStream << " but got " << fs << "\n";
      if (fs == FUNCTION_SPACE_MAX)
      {
        *outStream << "Note that this matches the default value defined by superclass, FUNCTION_SPACE_MAX.  Likely the subclass has failed to set the superclass functionSpace_ field.\n";
      }
      errorFlag++;
    }
  } catch (std::logic_error &err){
    *outStream << "UNEXPECTED ERROR !!! ----------------------------------------------------------\n";
    *outStream << err.what() << '\n';
    *outStream << "-------------------------------------------------------------------------------" << "\n\n";
    errorFlag = -1000;
  }

  if (errorFlag != 0)
    std::cout << "End Result: TEST FAILED\n";
  else
    std::cout << "End Result: TEST PASSED\n";
  
  // reset format state of std::cout
  std::cout.copyfmt(oldFormatState);
  return errorFlag;
}

} //end namespace
} //end namespace
