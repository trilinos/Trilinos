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

/** \file   test_01.hpp
    \brief  Unit tests for the Intrepid2::HCURL_TRI_I1_FEM class.
    \author Created by P. Bochev, D. Ridzal, K. Peterson, and Kyungjoo Kim.
*/


#include "Intrepid2_config.h"

#ifdef HAVE_INTREPID2_DEBUG
#define INTREPID2_TEST_FOR_DEBUG_ABORT_OVERRIDE_TO_CONTINUE
#endif

#include "Intrepid2_HCURL_TRI_I1_FEM.hpp"

#include "Teuchos_oblackholestream.hpp"
#include "Teuchos_RCP.hpp"

namespace Intrepid2 {
  
  namespace Test {
    
#define INTREPID2_TEST_ERROR_EXPECTED( S )                              \
    try {                                                               \
      ++nthrow;                                                         \
      S ;                                                               \
    }                                                                   \
    catch (std::logic_error err) {                                      \
      ++ncatch;                                                         \
      *outStream << "Expected Error ----------------------------------------------------------------\n"; \
      *outStream << err.what() << '\n';                                 \
      *outStream << "-------------------------------------------------------------------------------" << "\n\n"; \
    }
    
    template<typename ValueType, typename DeviceSpaceType>
    int HCURL_TRI_I1_FEM_Test01(const bool verbose) {
      
      Teuchos::RCP<std::ostream> outStream;
      Teuchos::oblackholestream bhs; // outputs nothing

      if (verbose)
        outStream = Teuchos::rcp(&std::cout, false);
      else
        outStream = Teuchos::rcp(&bhs, false);
      
      Teuchos::oblackholestream oldFormatState;
      oldFormatState.copyfmt(std::cout);
      
      typedef typename
        Kokkos::Impl::is_space<DeviceSpaceType>::host_mirror_space::execution_space HostSpaceType ;
      
      *outStream << "DeviceSpace::  "; DeviceSpaceType::print_configuration(*outStream, false);
      *outStream << "HostSpace::    ";   HostSpaceType::print_configuration(*outStream, false);
  
  *outStream 
    << "===============================================================================\n" 
    << "|                                                                             |\n" 
    << "|                 Unit Test (Basis_HCURL_TRI_I1_FEM)                          |\n" 
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

      typedef Kokkos::DynRankView<ValueType,DeviceSpaceType> DynRankView;
#define ConstructWithLabel(obj, ...) obj(#obj, __VA_ARGS__)
      const ValueType tol = Parameters::Tolerence;

      int errorFlag = 0;

      // for virtual function, value and point types are declared in the class
      typedef ValueType outputValueType;
      typedef ValueType pointValueType;
      Basis_HCURL_TRI_I1_FEM<DeviceSpaceType,outputValueType,pointValueType> triBasis;

    *outStream
    << "\n"
    << "===============================================================================\n"
    << "| TEST 1: Basis creation, exception testing                                   |\n"
    << "===============================================================================\n";
 /* 

  FieldContainer<double> triNodes(7, 2);
  triNodes(0,0) =  0.0;  triNodes(0,1) =  0.0;  
  triNodes(1,0) =  1.0;  triNodes(1,1) =  0.0;  
  triNodes(2,0) =  0.0;  triNodes(2,1) =  1.0;  
  // edge midpoints
  triNodes(3,0) =  0.5;  triNodes(3,1) =  0.0;  
  triNodes(4,0) =  0.5;  triNodes(4,1) =  0.5;  
  triNodes(5,0) =  0.0;  triNodes(5,1) =  0.5;  
  // Inside Triangle
  triNodes(6,0) =  0.25; triNodes(6,1) =  0.25;  

  // Generic array for the output values; needs to be properly resized depending on the operator type
  FieldContainer<double> vals;
*/

  try{
    ordinal_type nthrow = 0, ncatch = 0;
#ifdef HAVE_INTREPID2_DEBUG

    // Define array containing the 3 vertices of the reference TRI and its 3 edge midpoints.
    DynRankView ConstructWithLabel(triNodes, 7, 2);

    // Generic array for the output values; needs to be properly resized depending on the operator type
    const auto numFields = triBasis.getCardinality();
    const auto numPoints = triNodes.dimension(0);
    const auto spaceDim  = triBasis.getBaseCellTopology().getDimension();

    DynRankView vals;
    vals = DynRankView("vals", numFields, numPoints);

    {
    // exception #1: GRAD cannot be applied to HCURL functions 
    // resize vals to rank-3 container with dimensions (num. basis functions, num. points, arbitrary)
//    vals.resize(triBasis.getCardinality(), triNodes.dimension(0), 4 );
    INTREPID2_TEST_ERROR_EXPECTED( triBasis.getValues(vals, triNodes, OPERATOR_GRAD) );
    }
    {
    // exception #2: DIV cannot be applied to HCURL functions
    // resize vals to rank-2 container with dimensions (num. points, num. basis functions)
//    vals.resize(triBasis.getCardinality(), triNodes.dimension(0) );
    INTREPID2_TEST_ERROR_EXPECTED( triBasis.getValues(vals, triNodes, OPERATOR_DIV) );
    }
    // Exceptions 3-7: all bf tags/bf Ids below are wrong and should cause getDofOrdinal() and 
    // getDofTag() to access invalid array elements thereby causing bounds check exception
    {
    // exception #3
    INTREPID2_TEST_ERROR_EXPECTED( triBasis.getDofOrdinal(3,0,0) );
    // exception #4
    INTREPID2_TEST_ERROR_EXPECTED( triBasis.getDofOrdinal(1,1,1) );
    // exception #5
    INTREPID2_TEST_ERROR_EXPECTED( triBasis.getDofOrdinal(0,4,1) );
    // exception #6
    INTREPID2_TEST_ERROR_EXPECTED( triBasis.getDofTag(12) );
    // exception #7
    INTREPID2_TEST_ERROR_EXPECTED( triBasis.getDofTag(-1) );
    }
    // Exceptions 8-15 test exception handling with incorrectly dimensioned input/output arrays
    {
    // exception #8: input points array must be of rank-2
    DynRankView ConstructWithLabel(badPoints1, 4, 5, 3);
    INTREPID2_TEST_ERROR_EXPECTED( triBasis.getValues(vals, badPoints1, OPERATOR_VALUE) );
    }
    // exception #9 dimension 1 in the input point array must equal space dimension of the cell
    {
    DynRankView ConstructWithLabel(badPoints2, 4, triBasis.getBaseCellTopology().getDimension() + 1);
    INTREPID2_TEST_ERROR_EXPECTED( triBasis.getValues(vals, badPoints2, OPERATOR_VALUE) ); 
    }
    {
    // exception #10 output values must be of rank-3 for OPERATOR_VALUE in 2D
    DynRankView ConstructWithLabel(badVals1, 4, 3);
    INTREPID2_TEST_ERROR_EXPECTED( triBasis.getValues(badVals1, triNodes, OPERATOR_VALUE) ); 
    }
    {
    DynRankView ConstructWithLabel(badCurls1,4,3,2);
    // exception #11 output values must be of rank-2 for OPERATOR_CURL
    INTREPID2_TEST_ERROR_EXPECTED( triBasis.getValues(badCurls1, triNodes, OPERATOR_CURL) ); 
    }
    {
    // exception #12 incorrect 0th dimension of output array (must equal number of basis functions)
    DynRankView ConstructWithLabel(badVals2, triBasis.getCardinality() + 1, triNodes.dimension(0), triBasis.getBaseCellTopology().getDimension());
    INTREPID2_TEST_ERROR_EXPECTED( triBasis.getValues(badVals2, triNodes, OPERATOR_VALUE) ) ;
    }
    {
    // exception #13 incorrect 1st  dimension of output array (must equal number of points)
    DynRankView ConstructWithLabel(badVals3, triBasis.getCardinality(), triNodes.dimension(0) + 1, triBasis.getBaseCellTopology().getDimension() );
    INTREPID2_TEST_ERROR_EXPECTED( triBasis.getValues(badVals3, triNodes, OPERATOR_VALUE) ) ;
    }
    {
    // exception #14: incorrect 2nd dimension of output array for VALUE (must equal the space dimension)
    DynRankView ConstructWithLabel(badVals4, triBasis.getCardinality(), triNodes.dimension(0), triBasis.getBaseCellTopology().getDimension() - 1);
    INTREPID2_TEST_ERROR_EXPECTED( triBasis.getValues(badVals4, triNodes, OPERATOR_VALUE) ) ;
    } 
    // exception #15: D2 cannot be applied to HCURL functions 
    // resize vals to rank-3 container with dimensions (num. basis functions, num. points, arbitrary)
//    vals.resize(triBasis.getCardinality(), 
//                triNodes.dimension(0),  
//                Intrepid2::getDkCardinality(OPERATOR_D2, triBasis.getBaseCellTopology().getDimension()));
//    INTREPID2_TEST_ERROR_EXPECTED( triBasis.getValues(vals, triNodes, OPERATOR_D2) ); 
#endif
  // Check if number of thrown exceptions matches the one we expect
    if (nthrow != ncatch) {
      errorFlag++;
      *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
      *outStream << "# of catch ("<< ncatch << ") is different from # of throw (" << ncatch << ")\n";
    }
  } catch (std::logic_error err) {
    *outStream << "UNEXPECTED ERROR !!! ----------------------------------------------------------\n";
    *outStream << err.what() << '\n';
    *outStream << "-------------------------------------------------------------------------------" << "\n\n";
    errorFlag = -1000;
  }
  
  *outStream
    << "\n"
    << "===============================================================================\n"
    << "| TEST 2: correctness of tag to enum and enum to tag lookups                  |\n"
    << "===============================================================================\n";
  
  // all tags are on host space
  try{
    const auto numFields = triBasis.getCardinality();
    const auto allTags = triBasis.getAllDofTags();
    
    // Loop over all tags, lookup the associated dof enumeration and then lookup the tag again
    const auto dofTagSize = allTags.dimension(0);
    for (auto i = 0; i < dofTagSize; ++i) {
      const auto bfOrd  = triBasis.getDofOrdinal(allTags(i,0), allTags(i,1), allTags(i,2));
      
      const auto myTag = triBasis.getDofTag(bfOrd);
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
    for( auto bfOrd = 0; bfOrd < triBasis.getCardinality(); bfOrd++) {
      const auto myTag  = triBasis.getDofTag(bfOrd);
      const auto myBfOrd = triBasis.getDofOrdinal(myTag(0), myTag(1), myTag(2));
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
  catch (std::logic_error err){
    *outStream << err.what() << "\n\n";
    errorFlag = -1000;
  };
  
  *outStream
    << "\n"
    << "===============================================================================\n"
    << "| TEST 3: correctness of basis function values                                |\n"
    << "===============================================================================\n";
  
  outStream -> precision(20);
  
  // VALUE: correct values in (P,F,D) layout
  const ValueType basisValues[] = {
    1.000, 0, 0, 0, 0, -1.000, 1.000, 1.000, 0, 1.000, 0, 0, 0, 0, \
    -1.000, 0, -1.000, -1.000, 1.000, 0.5000, 0, 0.5000, 0, -0.5000, \
    0.5000, 0.5000, -0.5000, 0.5000, -0.5000, -0.5000, 0.5000, 0, \
    -0.5000, 0, -0.5000, -1.000, 0.7500, 0.2500, -0.2500, 0.2500, \
    -0.2500, -0.7500};
  
  // CURL: correct values in (P,F) layout
  const ValueType basisCurls[] = {   
    2.0,  2.0,  2.0,  
    2.0,  2.0,  2.0,  
    2.0,  2.0,  2.0,  
    2.0,  2.0,  2.0,  
    2.0,  2.0,  2.0,  
    2.0,  2.0,  2.0,  
    2.0,  2.0,  2.0  
  };
  
  try{
    DynRankView ConstructWithLabel(triNodesHost, 7, 2);
    triNodesHost(0,0) =  0.0;  triNodesHost(0,1) =  0.0;  
    triNodesHost(1,0) =  1.0;  triNodesHost(1,1) =  0.0;  
    triNodesHost(2,0) =  0.0;  triNodesHost(2,1) =  1.0;  
    // edge midpoints
    triNodesHost(3,0) =  0.5;  triNodesHost(3,1) =  0.0;  
    triNodesHost(4,0) =  0.5;  triNodesHost(4,1) =  0.5;  
    triNodesHost(5,0) =  0.0;  triNodesHost(5,1) =  0.5;  
    // Inside Triangle
    triNodesHost(6,0) =  0.25; triNodesHost(6,1) =  0.25;  

    auto triNodes = Kokkos::create_mirror_view(typename DeviceSpaceType::memory_space(), triNodesHost);
    Kokkos::deep_copy(triNodes, triNodesHost);

    // Dimensions for the output arrays:
    const auto numFields = triBasis.getCardinality();
    const auto numPoints = triNodes.dimension(0);
    const auto spaceDim  = triBasis.getBaseCellTopology().getDimension();        
    
    {
    // Generic array for values and curls that will be properly sized before each call
    DynRankView ConstructWithLabel(vals, numFields, numPoints, spaceDim);
    // Check VALUE of basis functions: resize vals to rank-3 container:
//    vals.resize(numFields, numPoints, spaceDim);
    triBasis.getValues(vals, triNodes, OPERATOR_VALUE);
    auto vals_host = Kokkos::create_mirror_view(typename HostSpaceType::memory_space(), vals);
    Kokkos::deep_copy(vals_host, vals);
    for (auto i = 0; i < numFields; ++i) {
      for (auto j = 0; j < numPoints; ++j) {
        for (auto k = 0; k < spaceDim; ++k) {
          
          // compute offset for (P,F,D) data layout: indices are P->j, F->i, D->k
           auto l = k + i * spaceDim + j * spaceDim * numFields;
           if (std::abs(vals_host(i,j,k) - basisValues[l]) > tol) {
             errorFlag++;
             *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";

             // Output the multi-index of the value where the error is:
             *outStream << " At multi-index { ";
             *outStream << i << " ";*outStream << j << " ";*outStream << k << " ";
             *outStream << "}  computed value: " << vals(i,j,k)
               << " but reference value: " << basisValues[l] << "\n";
            }
         }
      }
    }
    } 
    
    {
    // Check CURL of basis function: resize vals to rank-2 container
    DynRankView ConstructWithLabel(vals, numFields, numPoints);
//    vals.resize(numFields, numPoints);
    triBasis.getValues(vals, triNodes, OPERATOR_CURL);
    auto vals_host = Kokkos::create_mirror_view(typename HostSpaceType::memory_space(), vals);
    Kokkos::deep_copy(vals_host, vals);
    for (auto i = 0; i < numFields; ++i) {
      for (auto j = 0; j < numPoints; ++j) {
        int l =  i + j * numFields;
        if (std::abs(vals_host(i,j) - basisCurls[l]) > tol) {
          errorFlag++;
          *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
          
          // Output the multi-index of the value where the error is:
          *outStream << " At multi-index { ";
          *outStream << i << " ";*outStream << j << " ";
          *outStream << "}  computed curl component: " << vals(i,j)
            << " but reference curl component: " << basisCurls[l] << "\n";
        }
      }
    }  
    }
  } //end try
   
  // Catch unexpected errors
  catch (std::logic_error err) {
    *outStream << err.what() << "\n\n";
    errorFlag = -1000;
  };
 
  *outStream 
    << "\n"
    << "===============================================================================\n"
    << "| TEST 4: correctness of DoF locations                                        |\n"
    << "===============================================================================\n";

  try{
    const auto numFields = triBasis.getCardinality();
    const auto spaceDim  = triBasis.getBaseCellTopology().getDimension();

    // Check exceptions.
    ordinal_type nthrow = 0, ncatch = 0;
#ifdef HAVE_INTREPID2_DEBUG
    {
      DynRankView ConstructWithLabel(badVals, 1, 2, 3);
      INTREPID2_TEST_ERROR_EXPECTED( triBasis.getDofCoords(badVals) );
    }
    {
      DynRankView ConstructWithLabel(badVals, 4, 2);
      INTREPID2_TEST_ERROR_EXPECTED( triBasis.getDofCoords(badVals) );
    }
    {
      DynRankView ConstructWithLabel(badVals, 4, 3);
      INTREPID2_TEST_ERROR_EXPECTED( triBasis.getDofCoords(badVals) );
    }
#endif
//    cvals.resize(3,spaceDim);
//    INTREPID_TEST_COMMAND( coord_iface->getDofCoords(cvals), throwCounter, nException ); nException--;
    // Check if number of thrown exceptions matches the one we expect
    if (nthrow != ncatch) {
      errorFlag++;
      *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
      *outStream << "# of catch ("<< ncatch << ") is different from # of throw (" << ncatch << ")\n";
    }

    // Check mathematical correctness
    DynRankView ConstructWithLabel(bvals, numFields, numFields, spaceDim);
    DynRankView ConstructWithLabel(cvals, numFields, spaceDim); 
    triBasis.getDofCoords(cvals);
    triBasis.getValues(bvals, cvals, OPERATOR_VALUE);

    auto cvals_host = Kokkos::create_mirror_view(typename HostSpaceType::memory_space(), cvals);
    Kokkos::deep_copy(cvals_host, cvals);
    auto bvals_host = Kokkos::create_mirror_view(typename HostSpaceType::memory_space(), bvals);
    Kokkos::deep_copy(bvals_host, bvals);

    Kokkos::View<ValueType**, typename DynRankView::array_layout, typename HostSpaceType::execution_space> ConstructWithLabel(tangents, numFields, spaceDim); 
    tangents(0,0) =  1.0; tangents(0,1) =  0.0;
    tangents(1,0) = -1.0; tangents(1,1) =  1.0;
    tangents(2,0) =  0.0; tangents(2,1) = -1.0;

    char buffer[120];
    for (auto i=0; i<bvals_host.dimension(0); i++) { 
      for (auto j=0; j<bvals_host.dimension(1); j++) {

        double tangent = 0.0;
        for(auto d=0;d<spaceDim;d++)
           tangent += bvals_host(i,j,d)*tangents(j,d);

        if ((i != j) && (std::abs(tangent - 0.0) > tol )) {
          errorFlag++;
          sprintf(buffer, "\nValue of basis function %d at (%6.4e, %6.4e) is %6.4e but should be %6.4e!\n", i, cvals_host(i,0), cvals_host(i,1), tangent, 0.0);
          *outStream << buffer;
        }
        else if ((i == j) && (std::abs(tangent - 1.0) > tol )) {
          errorFlag++;
          sprintf(buffer, "\nValue of basis function %d at (%6.4e, %6.4e) is %6.4e but should be %6.4e!\n", i, cvals_host(i,0), cvals_host(i,1), tangent, 1.0);
          *outStream << buffer;
        }
      }
    }

  }
  catch (std::logic_error err){
    *outStream << err.what() << "\n\n";
    errorFlag = -1000;
  };

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
