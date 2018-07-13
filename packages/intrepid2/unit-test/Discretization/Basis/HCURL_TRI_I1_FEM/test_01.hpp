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

#include "Intrepid2_Basis.hpp"
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
    catch (std::exception err) {                                        \
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
      typedef Kokkos::DynRankView<ValueType,HostSpaceType> DynRankViewHost;
#define ConstructWithLabel(obj, ...) obj(#obj, __VA_ARGS__)
      const ValueType tol = tolerence();

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

  try{
    ordinal_type nthrow = 0, ncatch = 0;
#ifdef HAVE_INTREPID2_DEBUG

    // Define array containing the 3 vertices of the reference TRI and its 3 edge midpoints.
    DynRankView ConstructWithLabel(triNodes, 7, 2);

    // Generic array for the output values; needs to be properly resized depending on the operator type
    const ordinal_type cardinality = triBasis.getCardinality();
    const ordinal_type numPoints = triNodes.extent(0);

    DynRankView vals;
    vals = DynRankView("vals", cardinality, numPoints);

    {
    // exception #1: GRAD cannot be applied to HCURL functions 
      INTREPID2_TEST_ERROR_EXPECTED( triBasis.getValues(vals, triNodes, OPERATOR_GRAD) );
    }
    {
    // exception #2: DIV cannot be applied to HCURL functions
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
      INTREPID2_TEST_ERROR_EXPECTED( triBasis.getDofTag(cardinality) );
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
    // exception #11 output values must be of rank-2 for OPERATOR_CURL
      DynRankView ConstructWithLabel(badCurls1,4,3,2);
      INTREPID2_TEST_ERROR_EXPECTED( triBasis.getValues(badCurls1, triNodes, OPERATOR_CURL) ); 
    }
    {
    // exception #12 incorrect 0th dimension of output array (must equal number of basis functions)
      DynRankView ConstructWithLabel(badVals2, triBasis.getCardinality() + 1, triNodes.extent(0), triBasis.getBaseCellTopology().getDimension());
      INTREPID2_TEST_ERROR_EXPECTED( triBasis.getValues(badVals2, triNodes, OPERATOR_VALUE) ) ;
    }
    {
    // exception #13 incorrect 1st  dimension of output array (must equal number of points)
      DynRankView ConstructWithLabel(badVals3, triBasis.getCardinality(), triNodes.extent(0) + 1, triBasis.getBaseCellTopology().getDimension() );
      INTREPID2_TEST_ERROR_EXPECTED( triBasis.getValues(badVals3, triNodes, OPERATOR_VALUE) ) ;
    }
    {
    // exception #14: incorrect 2nd dimension of output array for VALUE (must equal the space dimension)
      DynRankView ConstructWithLabel(badVals4, triBasis.getCardinality(), triNodes.extent(0), triBasis.getBaseCellTopology().getDimension() - 1);
      INTREPID2_TEST_ERROR_EXPECTED( triBasis.getValues(badVals4, triNodes, OPERATOR_VALUE) ) ;
    } 
    // exception #15: D2 cannot be applied to HCURL functions 
    // resize vals to rank-3 container with dimensions (num. basis functions, num. points, arbitrary)
//    vals.resize(triBasis.getCardinality(), 
//                triNodes.extent(0),  
//                Intrepid2::getDkCardinality(OPERATOR_D2, triBasis.getBaseCellTopology().getDimension()));
//    INTREPID2_TEST_ERROR_EXPECTED( triBasis.getValues(vals, triNodes, OPERATOR_D2) ); 
#endif
  // Check if number of thrown exceptions matches the one we expect
    if (nthrow != ncatch) {
      errorFlag++;
      *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
      *outStream << "# of catch ("<< ncatch << ") is different from # of throw (" << nthrow << ")\n";
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
    const auto allTags = triBasis.getAllDofTags();
    
    // Loop over all tags, lookup the associated dof enumeration and then lookup the tag again
    const ordinal_type dofTagSize = allTags.extent(0);
    for (ordinal_type i = 0; i < dofTagSize; ++i) {
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
    for( ordinal_type bfOrd = 0; bfOrd < triBasis.getCardinality(); bfOrd++) {
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
    DynRankViewHost ConstructWithLabel(triNodesHost, 7, 2);
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
    const ordinal_type cardinality = triBasis.getCardinality();
    const ordinal_type numPoints = triNodes.extent(0);
    const ordinal_type spaceDim  = triBasis.getBaseCellTopology().getDimension();
    
    {
    // Check VALUE of basis functions: resize vals to rank-3 container:
    DynRankView ConstructWithLabel(vals, cardinality, numPoints, spaceDim);
    triBasis.getValues(vals, triNodes, OPERATOR_VALUE);
    auto vals_host = Kokkos::create_mirror_view(typename HostSpaceType::memory_space(), vals);
    Kokkos::deep_copy(vals_host, vals);
    for (ordinal_type i = 0; i < cardinality; ++i) {
      for (ordinal_type j = 0; j < numPoints; ++j) {
        for (ordinal_type k = 0; k < spaceDim; ++k) {
          
          // compute offset for (P,F,D) data layout: indices are P->j, F->i, D->k
           ordinal_type l = k + i * spaceDim + j * spaceDim * cardinality;
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
    // Check CURL of basis function: resize vals to rank-2 container
    DynRankView ConstructWithLabel(vals, cardinality, numPoints);
    triBasis.getValues(vals, triNodes, OPERATOR_CURL);
    auto vals_host = Kokkos::create_mirror_view(typename HostSpaceType::memory_space(), vals);
    Kokkos::deep_copy(vals_host, vals);
    for (ordinal_type i = 0; i < cardinality; ++i) {
      for (ordinal_type j = 0; j < numPoints; ++j) {
        ordinal_type l =  i + j * cardinality;
        if (std::abs(vals_host(i,j) - basisCurls[l]) > tol) {
          errorFlag++;
          *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
          
          // Output the multi-index of the value where the error is:
          *outStream << " At multi-index { ";
          *outStream << i << " ";*outStream << j << " ";
          *outStream << "}  computed curl component: " << vals_host(i,j)
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
    const ordinal_type cardinality = triBasis.getCardinality();
    const ordinal_type spaceDim  = triBasis.getBaseCellTopology().getDimension();

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
    // Check if number of thrown exceptions matches the one we expect
    if (nthrow != ncatch) {
      errorFlag++;
      *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
      *outStream << "# of catch ("<< ncatch << ") is different from # of throw (" << ncatch << ")\n";
    }

    // Check mathematical correctness
    DynRankView ConstructWithLabel(bvals, cardinality, cardinality, spaceDim);
    DynRankView ConstructWithLabel(cvals, cardinality, spaceDim);
    DynRankView ConstructWithLabel(dofCoeffs, cardinality, spaceDim);
    triBasis.getDofCoeffs(dofCoeffs);
    triBasis.getDofCoords(cvals);
    triBasis.getValues(bvals, cvals, OPERATOR_VALUE);

    auto cvals_host = Kokkos::create_mirror_view(typename HostSpaceType::memory_space(), cvals);
    Kokkos::deep_copy(cvals_host, cvals);
    auto dofCoeffs_host = Kokkos::create_mirror_view(typename HostSpaceType::memory_space(), dofCoeffs);
    Kokkos::deep_copy(dofCoeffs_host, dofCoeffs);
    auto bvals_host = Kokkos::create_mirror_view(typename HostSpaceType::memory_space(), bvals);
    Kokkos::deep_copy(bvals_host, bvals);

    char buffer[120];
    for (ordinal_type i=0; i<cardinality; ++i) {
      for (ordinal_type j=0; j<cardinality; ++j) {

        double dofValue = 0.0;
        for(ordinal_type d=0;d<spaceDim;++d)
          dofValue += bvals_host(i,j,d)*dofCoeffs_host(j,d);

        if ((i != j) && (std::abs(dofValue - 0.0) > tol )) {
          errorFlag++;
          sprintf(buffer, "\nValue of basis function %d at (%6.4e, %6.4e) is %6.4e but should be %6.4e!\n", i, cvals_host(i,0), cvals_host(i,1), dofValue, 0.0);
          *outStream << buffer;
        }
        else if ((i == j) && (std::abs(dofValue - 1.0) > tol )) {
          errorFlag++;
          sprintf(buffer, "\nValue of basis function %d at (%6.4e, %6.4e) is %6.4e but should be %6.4e!\n", i, cvals_host(i,0), cvals_host(i,1), dofValue, 1.0);
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
