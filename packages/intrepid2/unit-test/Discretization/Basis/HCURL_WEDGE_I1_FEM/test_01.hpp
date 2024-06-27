// @HEADER
// *****************************************************************************
//                           Intrepid2 Package
//
// Copyright 2007 NTESS and the Intrepid2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/** \file   test_01.hpp
    \brief  Unit tests for the Intrepid2::C_WEDGE_I1_FEM class.
    \author Created by P. Bochev, D. Ridzal, K. Peterson, and Kyungjoo Kim.
*/

#include "Intrepid2_config.h"

#ifdef HAVE_INTREPID2_DEBUG
#define INTREPID2_TEST_FOR_DEBUG_ABORT_OVERRIDE_TO_CONTINUE
#endif

#include "Intrepid2_Basis.hpp"
#include "Intrepid2_HCURL_WEDGE_I1_FEM.hpp"

#include "Teuchos_oblackholestream.hpp"
#include "Teuchos_RCP.hpp"

#include "packages/intrepid2/unit-test/Discretization/Basis/Macros.hpp"

namespace Intrepid2 {
  
namespace Test {

  template<typename ValueType, typename DeviceType>
  int HCURL_WEDGE_I1_FEM_Test01(const bool verbose) {
      
    Teuchos::RCP<std::ostream> outStream;
    Teuchos::oblackholestream bhs; // outputs nothing

    if (verbose)
      outStream = Teuchos::rcp(&std::cout, false);
    else
      outStream = Teuchos::rcp(&bhs, false);
      
    Teuchos::oblackholestream oldFormatState;
    oldFormatState.copyfmt(std::cout);
    using DeviceSpaceType = typename DeviceType::execution_space;       
    typedef typename
      Kokkos::DefaultHostExecutionSpace HostSpaceType ;
      
    *outStream << "DeviceSpace::  "; DeviceSpaceType().print_configuration(*outStream, false);
    *outStream << "HostSpace::    ";   HostSpaceType().print_configuration(*outStream, false);

  *outStream 
    << "\n"
    << "===============================================================================\n" 
    << "|                                                                             |\n" 
    << "|                 Unit Test (Basis_HCURL_WEDGE_I1_FEM)                        |\n" 
    << "|                                                                             |\n" 
    << "|     1) Conversion of Dof tags into Dof ordinals and back                    |\n" 
    << "|     2) Basis values for VALUE and CURL operators                            |\n" 
    << "|                                                                             |\n" 
    << "|  Questions? Contact  Pavel Bochev  (pbboche@sandia.gov),                    |\n" 
    << "|                      Denis Ridzal  (dridzal@sandia.gov),                    |\n" 
    << "|                      Kara Peterson (kjpeter@sandia.gov).                    |\n" 
    << "|                      Kyungjoo Kim  (kyukim@sandia.gov).                     |\n"
    << "|                                                                             |\n" 
    << "|  Intrepid's website: http://trilinos.sandia.gov/packages/intrepid           |\n" 
    << "|  Trilinos website:   http://trilinos.sandia.gov                             |\n" 
    << "|                                                                             |\n" 
    << "===============================================================================\n";

      typedef Kokkos::DynRankView<ValueType,DeviceType> DynRankView;
      typedef Kokkos::DynRankView<ValueType,HostSpaceType> DynRankViewHost;

      const ValueType tol = tolerence();

      int errorFlag = 0;

      // for virtual function, value and point types are declared in the class
      typedef ValueType outputValueType;
      typedef ValueType pointValueType;
      Basis_HCURL_WEDGE_I1_FEM<DeviceType,outputValueType,pointValueType> wedgeBasis;

  *outStream 
    << "\n"
    << "===============================================================================\n"
    << "| TEST 1: Basis creation, exception testing                                   |\n"
    << "===============================================================================\n";

  try{
    ordinal_type nthrow = 0, ncatch = 0;
#ifdef HAVE_INTREPID2_DEBUG
    // Define array containing the 4 vertices of the reference WEDGE and its center.
    DynRankView ConstructWithLabel(wedgeNodes, 12, 3);

    // Generic array for the output values; needs to be properly resized depending on the operator type
    const ordinal_type numFields = wedgeBasis.getCardinality();
    const ordinal_type numPoints = wedgeNodes.extent(0);
    const ordinal_type spaceDim  = wedgeBasis.getBaseCellTopology().getDimension();

    DynRankView vals ("vals", numFields, numPoints);
    DynRankView vals_vec ("vals", numFields, numPoints, spaceDim);

    {
    // exception #1: GRAD cannot be applied to HCURL functions 
    // resize vals to rank-3 container with dimensions (num. basis functions, num. points, arbitrary)
      INTREPID2_TEST_ERROR_EXPECTED( wedgeBasis.getValues(vals_vec, wedgeNodes, OPERATOR_GRAD) );

    // exception #2: DIV cannot be applied to HCURL functions
    // resize vals to rank-2 container with dimensions (num. basis functions, num. points)
      INTREPID2_TEST_ERROR_EXPECTED( wedgeBasis.getValues(vals, wedgeNodes, OPERATOR_DIV) );
    }   
    // Exceptions 3-7: all bf tags/bf Ids below are wrong and should cause getDofOrdinal() and 
    // getDofTag() to access invalid array elements thereby causing bounds check exception
    {
    // exception #3
      INTREPID2_TEST_ERROR_EXPECTED( wedgeBasis.getDofOrdinal(3,0,0) );
    // exception #4
      INTREPID2_TEST_ERROR_EXPECTED( wedgeBasis.getDofOrdinal(1,1,1) );
    // exception #5
      INTREPID2_TEST_ERROR_EXPECTED( wedgeBasis.getDofOrdinal(0,4,1) );
    // exception #6
      INTREPID2_TEST_ERROR_EXPECTED( wedgeBasis.getDofTag(numFields) );
    // exception #7
      INTREPID2_TEST_ERROR_EXPECTED( wedgeBasis.getDofTag(-1) );
    }

    // Exceptions 8- test exception handling with incorrectly dimensioned input/output arrays
    {
    // exception #8: input points array must be of rank-2
      DynRankView ConstructWithLabel(badPoints1, 4, 5, 3);
      INTREPID2_TEST_ERROR_EXPECTED( wedgeBasis.getValues(vals, badPoints1, OPERATOR_VALUE) );
    } 
    {
    // exception #9 dimension 1 in the input point array must equal space dimension of the cell
      DynRankView ConstructWithLabel(badPoints2, 4, 2);
      INTREPID2_TEST_ERROR_EXPECTED( wedgeBasis.getValues(vals, badPoints2, OPERATOR_VALUE) );
    }
    {
    // exception #10 output values must be of rank-3 for OPERATOR_VALUE
      DynRankView ConstructWithLabel(badVals1, 4, 3);
      INTREPID2_TEST_ERROR_EXPECTED( wedgeBasis.getValues(badVals1, wedgeNodes, OPERATOR_VALUE) );
    // exception #11 output values must be of rank-3 for OPERATOR_CURL
      INTREPID2_TEST_ERROR_EXPECTED( wedgeBasis.getValues(badVals1, wedgeNodes, OPERATOR_CURL) );
    }
    { 
    // exception #12 incorrect 0th dimension of output array (must equal number of basis functions)
      DynRankView ConstructWithLabel(badVals2, numFields + 1, numPoints, 3);
      INTREPID2_TEST_ERROR_EXPECTED( wedgeBasis.getValues(badVals2, wedgeNodes, OPERATOR_VALUE) );
    }
    {
    // exception #13 incorrect 1st dimension of output array (must equal number of points)
      DynRankView ConstructWithLabel(badVals3, numFields, numPoints + 1, 3);
      INTREPID2_TEST_ERROR_EXPECTED( wedgeBasis.getValues(badVals3, wedgeNodes, OPERATOR_VALUE) );
    }
    {
    // exception #14: incorrect 2nd dimension of output array (must equal the space dimension)
      DynRankView ConstructWithLabel(badVals4, numFields, numPoints, 4);
      INTREPID2_TEST_ERROR_EXPECTED( wedgeBasis.getValues(badVals4, wedgeNodes, OPERATOR_VALUE) );
    
    // exception #15: incorrect 2nd dimension of output array (must equal the space dimension)
      INTREPID2_TEST_ERROR_EXPECTED( wedgeBasis.getValues(badVals4, wedgeNodes, OPERATOR_CURL) );
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
    const ordinal_type numFields = wedgeBasis.getCardinality();
    const auto allTags = wedgeBasis.getAllDofTags();

    // Loop over all tags, lookup the associated dof enumeration and then lookup the tag again
    const ordinal_type dofTagSize = allTags.extent(0);
    for (ordinal_type i = 0; i < dofTagSize; ++i) {
      auto bfOrd  = wedgeBasis.getDofOrdinal(allTags(i,0), allTags(i,1), allTags(i,2));
      
      const auto myTag = wedgeBasis.getDofTag(bfOrd);
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
    for( ordinal_type bfOrd = 0; bfOrd < numFields; ++bfOrd) {
      const auto  myTag  = wedgeBasis.getDofTag(bfOrd);
      const auto myBfOrd = wedgeBasis.getDofOrdinal(myTag(0), myTag(1), myTag(2));
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
  } catch (std::logic_error &err){
    *outStream << err.what() << "\n\n";
    errorFlag = -1000;
  };
  
  *outStream
    << "\n"
    << "===============================================================================\n"
    << "| TEST 3: correctness of basis function values                                |\n"
    << "===============================================================================\n";
  
  outStream -> precision(20);
  
  // VALUE: Each row pair gives the 9x3 correct basis set values at an evaluation point: (P,F,D) layout
  const ValueType basisValues[] = {
    2.0, 0, 0, 0, 0, 0, 0, -2.0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, \
    0, 1.0, 0, 0, 0, 0, 0, 0, 2.0, 2.0, 0, 0, 2.0, 0, 0, \
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.0, 0, 0, 0, 0, \
    0, 0, -2.0, 0, 0, -2.0, -2.0, 0, 0, 0, 0, 0, 0, 0, 0, 0, \
    0, 0, 0, 0, 0, 0, 0, 0, 0, 1.0, 0, 0, 0, 0, 0, 0, 0, 0, 0, \
    2.0, 0, 0, 0, 0, 0, 0, -2.0, 0, 0, 0, 1.0, 0, 0, 0, 0, \
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2.0, 2.0, 0, 0, 2.0, 0, \
    0, 0, 0, 0, 0, 0, 0, 0, 1.0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, \
    0, 0, 0, -2.0, 0, 0, -2.0, -2.0, 0, 0, 0, 0, 0, 0, 0, 0, \
    0, 1.0, 1.0, 0.5, 0, -1.0, 0.5, 0, \
    -1.0, -1.5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.25, \
    0, 0, 0.25, 0, 0, 0.5, 0.75, 0.5, 0, -0.25, \
    0.5, 0, -0.25, -0.5, 0, 0.75, 0.5, 0, \
    -0.25, 0.5, 0, -0.25, -0.5, 0, 0, 0, 0.25, 0, \
    0, 0.5, 0, 0, 0.25, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.5, \
    0.5, 0, -0.5, 0.5, 0, -0.5, -1.5, 0, 0, 0, \
    0.5, 0, 0, 0.25, 0, 0, 0.25, 0.25, 0.0625, 0, 0, \
    0.0625, 0, 0, -0.1875, 0, 1.75, 0.4375, 0, 0, 0.4375, \
    0, 0, -1.3125, 0, 0, 0, 0.75, 0, 0, 0.25, 0, 0, 0, \
    0.625, 0, 0, -0.625, 0, 0, -0.625, -1.25, 0, 0.375, \
    0, 0, -0.375, 0, 0, -0.375, -0.75, 0, 0, 0, 0.5, 0, 0, \
    0, 0, 0, 0.5, 0.5, 0.5, 0, -0.5, 0.5, 0, \
    -0.5, -0.5, 0, 0.5, 0.5, 0, -0.5, 0.5, \
    0, -0.5, -0.5, 0, 0, 0, 0, 0, 0, 0.5, 0, 0, 0.5
  };
  
  // CURL: each row pair gives the 9x3 correct values of the curls of the 9 basis functions: (P,F,D) layout
  const ValueType basisCurls[] = {  
    0, -1.0, 4.0, 0, 0, 4.0, -1.0, 0, 4.0, 0, \
    1.0, 0, 0, 0, 0, 1.0, 0, 0, -1.0, 1.0, 0, 0, \
    -1.0, 0, 1.0, 0, 0, 1.0, -1.0, 4.0, 1.0, \
    0, 4.0, 0, 0, 4.0, -1.0, 1.0, 0, -1.0, 0, 0, \
    0, 0, 0, -1.0, 1.0, 0, 0, -1.0, 0, 1.0, 0, 0, 0, \
    0, 4.0, 0, 1.0, 4.0, -1.0, 1.0, 4.0, 0, 0, \
    0, 0, -1.0, 0, 1.0, -1.0, 0, -1.0, 1.0, 0, \
    0, -1.0, 0, 1.0, 0, 0, 0, -1.0, 0, 0, 0, 0, -1.0, \
    0, 0, 0, 1.0, 4.0, 0, 0, 4.0, 1.0, 0, 4.0, \
    -1.0, 1.0, 0, 0, -1.0, 0, 1.0, 0, 0, 1.0, \
    -1.0, 0, 1.0, 0, 0, 0, 0, 0, -1.0, 1.0, 4.0, \
    -1.0, 0, 4.0, 0, 0, 4.0, -1.0, 1.0, 0, 0, \
    -1.0, 0, 1.0, 0, 0, 0, 0, 0, 0, 1.0, 0, -1.0, \
    1.0, 0, 0, 0, 4.0, 0, -1.0, 4.0, 1.0, \
    -1.0, 4.0, -1.0, 1.0, 0, 0, -1.0, 0, \
    1.0, 0, 0, 0.25, -0.5, 4.0, 0.25, 0.5, \
    4.0, -0.75, 0.5, 4.0, -0.25, 0.5, 0, \
    -0.25, -0.5, 0, 0.75, -0.5, 0, -1.0, 1.0, \
    0, 0, -1.0, 0, 1.0, 0, 0, 0.5, -0.75, 2.0, \
    0.5, 0.25, 2.0, -0.5, 0.25, 2.0, -0.5, \
    0.75, 2.0, -0.5, -0.25, 2.0, 0.5, \
    -0.25, 2.0, -1.0, 1.0, 0, 0, -1.0, 0, \
    1.0, 0, 0, 0.25, -0.75, 0, 0.25, 0.25, 0, \
    -0.75, 0.25, 0, -0.25, 0.75, 4.0, -0.25, \
    -0.25, 4.0, 0.75, -0.25, 4.0, -1.0, \
    1.0, 0, 0, -1.0, 0, 1.0, 0, 0, 0.25, -1.0, \
    0.5, 0.25, 0, 0.5, -0.75, 0, 0.5, -0.25, \
    1.0, 3.5, -0.25, 0, 3.5, 0.75, 0, 3.5, \
    -1.0, 1.0, 0, 0, -1.0, 0, 1.0, 0, 0, 0, \
    -0.5, 2.5, 0, 0.5, 2.5, -1.0, 0.5, \
    2.5, 0, 0.5, 1.5, 0, -0.5, 1.5, 1.0, \
    -0.5, 1.5, -1.0, 1.0, 0, 0, -1.0, 0, \
    1.0, 0, 0, 0.5, -0.5, 2.0, 0.5, 0.5, \
    2.0, -0.5, 0.5, 2.0, -0.5, 0.5, 2.0, \
    -0.5, -0.5, 2.0, 0.5, -0.5, 2.0, \
    -1.0, 1.0, 0, 0, -1.0, 0, 1.0, 0, 0
  };
  
  try{
    
    DynRankViewHost ConstructWithLabel(wedgeNodesHost, 12, 3);
    wedgeNodesHost(0,0) =  0.0;  wedgeNodesHost(0,1) =  0.0;  wedgeNodesHost(0,2) = -1.0;
    wedgeNodesHost(1,0) =  1.0;  wedgeNodesHost(1,1) =  0.0;  wedgeNodesHost(1,2) = -1.0;
    wedgeNodesHost(2,0) =  0.0;  wedgeNodesHost(2,1) =  1.0;  wedgeNodesHost(2,2) = -1.0;
    wedgeNodesHost(3,0) =  0.0;  wedgeNodesHost(3,1) =  0.0;  wedgeNodesHost(3,2) =  1.0;
    wedgeNodesHost(4,0) =  1.0;  wedgeNodesHost(4,1) =  0.0;  wedgeNodesHost(4,2) =  1.0;
    wedgeNodesHost(5,0) =  0.0;  wedgeNodesHost(5,1) =  1.0;  wedgeNodesHost(5,2) =  1.0;
  
    wedgeNodesHost(6,0) =  0.25; wedgeNodesHost(6,1) =  0.5;  wedgeNodesHost(6,2) = -1.0;
    wedgeNodesHost(7,0) =  0.5;  wedgeNodesHost(7,1) =  0.25; wedgeNodesHost(7,2) =  0.0;
    wedgeNodesHost(8,0) =  0.25; wedgeNodesHost(8,1) =  0.25; wedgeNodesHost(8,2) =  1.0;
    wedgeNodesHost(9,0) =  0.25; wedgeNodesHost(9,1) =  0.0;  wedgeNodesHost(9,2) =  0.75;
    wedgeNodesHost(10,0)=  0.0;  wedgeNodesHost(10,1)=  0.5;  wedgeNodesHost(10,2)= -0.25;
    wedgeNodesHost(11,0)=  0.5;  wedgeNodesHost(11,1)=  0.5;  wedgeNodesHost(11,2)=  0.0;

    const auto wedgeNodes = Kokkos::create_mirror_view(typename DeviceType::memory_space(), wedgeNodesHost);
    Kokkos::deep_copy(wedgeNodes, wedgeNodesHost);
        
    // Dimensions for the output arrays:
    const ordinal_type numFields = wedgeBasis.getCardinality();
    const ordinal_type numPoints = wedgeNodes.extent(0);
    const ordinal_type spaceDim  = wedgeBasis.getBaseCellTopology().getDimension();
    
    {
    // Check VALUE of basis functions: resize vals to rank-3 container:
    DynRankView ConstructWithLabel(vals, numFields, numPoints, spaceDim);
    wedgeBasis.getValues(vals, wedgeNodes, OPERATOR_VALUE);
    auto vals_host = Kokkos::create_mirror_view(typename HostSpaceType::memory_space(), vals);
    Kokkos::deep_copy(vals_host, vals);
    for (ordinal_type i = 0; i < numFields; ++i) {
      for (ordinal_type j = 0; j < numPoints; ++j) {
        for (ordinal_type k = 0; k < spaceDim; ++k) {
          
          // compute offset for (P,F,D) data layout: indices are P->j, F->i, D->k
           const ordinal_type l = k + i * spaceDim + j * spaceDim * numFields;
           if (std::abs(vals_host(i,j,k) - basisValues[l]) > tol) {
             errorFlag++;
             *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";

             // Output the multi-index of the value where the error is:
             *outStream << " At multi-index { ";
             *outStream << i << " ";*outStream << j << " ";*outStream << k << " ";
             *outStream << "}  computed value: " << vals_host(i,j,k)
               << " but reference value: " << basisValues[l] << "\n";
            }
         }
      }
    }
    }

    {
    // Check CURL of basis function: resize vals to rank-3 container
    DynRankView ConstructWithLabel(vals, numFields, numPoints, spaceDim);
    wedgeBasis.getValues(vals, wedgeNodes, OPERATOR_CURL);
    auto vals_host = Kokkos::create_mirror_view(typename HostSpaceType::memory_space(), vals);
    Kokkos::deep_copy(vals_host, vals);
    for (ordinal_type i = 0; i < numFields; ++i) {
      for (ordinal_type j = 0; j < numPoints; ++j) {
        for (ordinal_type k = 0; k < spaceDim; ++k) {
          
          // compute offset for (P,F,D) data layout: indices are P->j, F->i, D->k
           const ordinal_type l = k + i * spaceDim + j * spaceDim * numFields;
           if (std::abs(vals_host(i,j,k) - basisCurls[l]) > tol) {
             errorFlag++;
             *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";

             // Output the multi-index of the value where the error is:
             *outStream << " At multi-index { ";
             *outStream << i << " ";*outStream << j << " ";*outStream << k << " ";
             *outStream << "}  computed curl component: " << vals_host(i,j,k)
               << " but reference curl component: " << basisCurls[l] << "\n";
            }
         }
      }
    }
    }
    
   }    
  
  // Catch unexpected errors
  catch (std::logic_error &err) {
    *outStream << err.what() << "\n\n";
    errorFlag = -1000;
  };

  *outStream
  << "\n"
  << "===============================================================================\n"
  << "| TEST 4: Function Space is Correct                                           |\n"
  << "===============================================================================\n";
  
  try {
    const EFunctionSpace fs = wedgeBasis.getFunctionSpace();
    
    if (fs != FUNCTION_SPACE_HCURL)
    {
      *outStream << std::setw(70) << "------------- TEST FAILURE! -------------" << "\n";
      
      // Output the multi-index of the value where the error is:
      *outStream << " Expected a function space of FUNCTION_SPACE_HCURL (enum value " << FUNCTION_SPACE_HCURL << "),";
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
