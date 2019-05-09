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
    \brief  Unit tests for the Intrepid2::HDIV_TET_I1_FEM class.
    \author Created by P. Bochev, D. Ridzal, and K. Peterson.
*/

#include "Intrepid2_config.h"

#ifdef HAVE_INTREPID2_DEBUG
#define INTREPID2_TEST_FOR_DEBUG_ABORT_OVERRIDE_TO_CONTINUE
#endif

#include "Intrepid2_Basis.hpp"
#include "Intrepid2_HDIV_TET_I1_FEM.hpp"

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
    int HDIV_TET_I1_FEM_Test01(const bool verbose) {

      Teuchos::RCP<std::ostream> outStream;
      Teuchos::oblackholestream bhs; // outputs nothing

      if (verbose)
        outStream = Teuchos::rcp(&std::cout, false);
      else
        outStream = Teuchos::rcp(&bhs,       false);

      Teuchos::oblackholestream oldFormatState;
      oldFormatState.copyfmt(std::cout);

      typedef typename
        Kokkos::Impl::is_space<DeviceSpaceType>::host_mirror_space::execution_space HostSpaceType ;

      *outStream << "DeviceSpace::  "; DeviceSpaceType::print_configuration(*outStream, false);
      *outStream << "HostSpace::    ";   HostSpaceType::print_configuration(*outStream, false);
  
      *outStream                                                       
        << "===============================================================================\n"
        << "|                                                                             |\n"
        << "|                 Unit Test (Basis_HDIV_TET_I1_FEM)                           |\n"
        << "|                                                                             |\n"
        << "|     1) Conversion of Dof tags into Dof ordinals and back                    |\n"
        << "|     2) Basis values for VALUE and HDIV operators                            |\n"
        << "|                                                                             |\n"
        << "|  Questions? Contact  Pavel Bochev  (pbboche@sandia.gov),                    |\n"
        << "|                      Denis Ridzal  (dridzal@sandia.gov),                    |\n"
        << "|                      Kara Peterson (kjpeter@sandia.gov).                    |\n"
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
      Basis_HDIV_TET_I1_FEM<DeviceSpaceType,outputValueType,pointValueType> tetBasis;

      *outStream
      << "\n"
      << "===============================================================================\n"
      << "| TEST 1: constructors and exceptions                                         |\n"
      << "===============================================================================\n";

      try {
#ifdef HAVE_INTREPID2_DEBUG
        ordinal_type nthrow = 0, ncatch = 0;

        // Define array containing the 4 vertices of the reference TET and its center.  
        DynRankView ConstructWithLabel(tetNodes, 10, 3);

        const auto numFields = tetBasis.getCardinality();
        const auto numPoints = tetNodes.extent(0);
        const auto spaceDim  = tetBasis.getBaseCellTopology().getDimension();

        DynRankView vals("vals", numFields, numPoints, spaceDim);

        // exception #1: GRAD cannot be applied to HDIV functions
        INTREPID2_TEST_ERROR_EXPECTED( tetBasis.getValues(vals, tetNodes, OPERATOR_GRAD));

        // exception #2: CURL cannot be applied to HDIV functions
        INTREPID2_TEST_ERROR_EXPECTED( tetBasis.getValues(vals, tetNodes, OPERATOR_CURL));

        // Exceptions 3-7: all bf tags/bf Ids below are wrong and should cause getDofOrdinal() and
        // getDofTag() to access invalid array elements thereby causing bounds check exception
        // exception #3
        INTREPID2_TEST_ERROR_EXPECTED( tetBasis.getDofOrdinal(3,0,0));
        // exception #4
        INTREPID2_TEST_ERROR_EXPECTED( tetBasis.getDofOrdinal(1,1,1));
        // exception #5
        INTREPID2_TEST_ERROR_EXPECTED( tetBasis.getDofOrdinal(0,4,1));
        // exception #6
        INTREPID2_TEST_ERROR_EXPECTED( tetBasis.getDofTag(numFields));
        // exception #7
        INTREPID2_TEST_ERROR_EXPECTED( tetBasis.getDofTag(-1));

        // Exceptions 8-15 test exception handling with incorrectly dimensioned input/output arrays
        // exception #8: input points array must be of rank-2
        DynRankView ConstructWithLabel(badPoints1, 4, 5, 3);
        INTREPID2_TEST_ERROR_EXPECTED( tetBasis.getValues(vals, badPoints1, OPERATOR_VALUE));

        // exception #9 dimension 1 in the input point array must equal space dimension of the cell
        DynRankView ConstructWithLabel(badPoints2, 4, 2);
        INTREPID2_TEST_ERROR_EXPECTED( tetBasis.getValues(vals, badPoints2, OPERATOR_VALUE));
        
        // exception #10 output values must be of rank-3 for OPERATOR_VALUE
        DynRankView ConstructWithLabel(badVals1, 4, 3);
        INTREPID2_TEST_ERROR_EXPECTED( tetBasis.getValues(badVals1, tetNodes, OPERATOR_VALUE));

        // exception #11 output values must be of rank-2 for OPERATOR_DIV
        DynRankView ConstructWithLabel(badVals2, 4, 3, 1);
        INTREPID2_TEST_ERROR_EXPECTED( tetBasis.getValues(badVals2, tetNodes, OPERATOR_VALUE));

        // exception #12 incorrect 0th dimension of output array (must equal number of basis functions)
        DynRankView ConstructWithLabel(badVals3, tetBasis.getCardinality() + 1, tetNodes.extent(0), 3);
        INTREPID2_TEST_ERROR_EXPECTED( tetBasis.getValues(badVals3, tetNodes, OPERATOR_VALUE));

        // exception #13 incorrect 0th dimension of output array (must equal number of basis functions)
        DynRankView ConstructWithLabel(badVals4, tetBasis.getCardinality() + 1, tetNodes.extent(0));
        INTREPID2_TEST_ERROR_EXPECTED( tetBasis.getValues(badVals4, tetNodes, OPERATOR_DIV));

        // exception #14 incorrect 1st dimension of output array (must equal number of points)
        DynRankView ConstructWithLabel(badVals5, tetBasis.getCardinality(), tetNodes.extent(0) + 1, 3);
        INTREPID2_TEST_ERROR_EXPECTED( tetBasis.getValues(badVals5, tetNodes, OPERATOR_VALUE));
    
        // exception #15 incorrect 1st dimension of output array (must equal number of points)
        DynRankView ConstructWithLabel(badVals6, tetBasis.getCardinality(), tetNodes.extent(0) + 1);
        INTREPID2_TEST_ERROR_EXPECTED( tetBasis.getValues(badVals6, tetNodes, OPERATOR_DIV));
    
        // exception #16: incorrect 2nd dimension of output array (must equal the space dimension)
        DynRankView ConstructWithLabel(badVals7, tetBasis.getCardinality(), tetNodes.extent(0), 4);
        INTREPID2_TEST_ERROR_EXPECTED( tetBasis.getValues(badVals7, tetNodes, OPERATOR_VALUE));
#endif
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
      
      try {
        const auto allTags = tetBasis.getAllDofTags();

        // Loop over all tags, lookup the associated dof enumeration and then lookup the tag again
        const auto dofTagSize = allTags.extent(0);
        for (unsigned i = 0; i < dofTagSize; i++) {
          int bfOrd  = tetBasis.getDofOrdinal(allTags(i,0), allTags(i,1), allTags(i,2));

          auto myTag = tetBasis.getDofTag(bfOrd);
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
        for( int bfOrd = 0; bfOrd < tetBasis.getCardinality(); bfOrd++) {
          auto myTag  = tetBasis.getDofTag(bfOrd);
          int myBfOrd = tetBasis.getDofOrdinal(myTag(0), myTag(1), myTag(2));
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
      } catch (std::logic_error err) {
        *outStream << err.what() << "\n\n";
        errorFlag = -1000;
      }
  
      *outStream 
        << "\n"
        << "===============================================================================\n"
        << "| TEST 3: correctness of basis function values                                |\n"
        << "===============================================================================\n";
  
      outStream -> precision(20);
  
      // VALUE: Correct basis values in (P,F,D) format: each row gives the 4x3 correct basis set values
      // at an evaluation point. Note that getValues returns results as an (F,P,D) array.
      double basisValues[] = {
        // 4 vertices
        0.,-2.0,0.,    0.,0.,0.,    -2.0,0.,0.,     0.,0.,-2.0,
        2.0,-2.0,0.,   2.0,0.,0.,    0.,0.,0.,      2.0,0.,-2.0,
        0.,0.,0.,      0.,2.0,0.,   -2.0,2.0,0.,    0,2.0,-2.0,
        0.,-2.0,2.0,   0.,0.,2.0,   -2.0,0.,2.0,    0.,0.,0.,
        // 6 edge midpoints
        1.0,-2.0,0.,   1.0,0.,0.,    -1.0,0.,0.,     1.0,0.,-2.0,
        1.0,-1.0,0.,   1.0,1.0,0.,   -1.0,1.0,0.,    1.0,1.0,-2.0,
        0.,-1.0,0.,    0.,1.0,0.,    -2.0,1.0,0.,    0.,1.0,-2.0,
        0.,-2.0,1.0,   0.,0.,1.0,    -2.0,0.,1.0,    0.,0.,-1.0,
        1.0,-2.0,1.0,  1.0,0.,1.0,   -1.0,0.,1.0,    1.0,0.,-1.0,
        0.,-1.0,1.0,   0.,1.0,1.0,   -2.0,1.0,1.0,   0.,1.0,-1.0
        // bf0         bf1                bf2            bf3
      };
  
      // DIV: each row gives the 4 correct values of the divergence of the 4 basis functions
      double basisDivs[] = {
        // 4 vertices
         6.0, 6.0, 6.0, 6.0,
         6.0, 6.0, 6.0, 6.0,
         6.0, 6.0, 6.0, 6.0,
         6.0, 6.0, 6.0, 6.0,
        // 6 edge midpoints
         6.0, 6.0, 6.0, 6.0,
         6.0, 6.0, 6.0, 6.0,
         6.0, 6.0, 6.0, 6.0,
         6.0, 6.0, 6.0, 6.0,
         6.0, 6.0, 6.0, 6.0,
         6.0, 6.0, 6.0, 6.0
      };
  
      try {
        // Define array containing the 4 vertices of the reference TET and its 6 edge midpoints.
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
  
        const auto tetNodes = Kokkos::create_mirror_view(typename DeviceSpaceType::memory_space(), tetNodesHost);
        Kokkos::deep_copy(tetNodes, tetNodesHost);
  
          // Dimensions for the output arrays:
        const auto numFields = tetBasis.getCardinality();
        const auto numPoints = tetNodes.extent(0);
        const auto spaceDim  = tetBasis.getBaseCellTopology().getDimension();

        {
          // Check VALUE of basis functions:
          DynRankView ConstructWithLabel(vals, numFields, numPoints, spaceDim);
          tetBasis.getValues(vals, tetNodes, OPERATOR_VALUE);
          const auto vals_host = Kokkos::create_mirror_view(typename HostSpaceType::memory_space(), vals);
          Kokkos::deep_copy(vals_host, vals);
          for (int i = 0; i < numFields; i++) {
            for (size_type j = 0; j < numPoints; j++) {
              for (size_type k = 0; k < spaceDim; k++) {
                // basisValues is (P,F,D) array so its multiindex is (j,i,k) and not (i,j,k)!
                 int l = k + i * spaceDim + j * spaceDim * numFields;
                 if (std::abs(vals_host(i,j,k) - basisValues[l]) > tol) {
                   errorFlag++;
                   *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";

                   // Output the multi-index of the value where the error is:
                   *outStream << " At (Field,Point,Dim) multi-index { ";
                   *outStream << i << " ";*outStream << j << " ";*outStream << k << " ";
                   *outStream << "}  computed value: " << vals_host(i,j,k)
                     << " but reference value: " << basisValues[l] << "\n";
                  }
               }
            }
          }
        }

        {
          // Check DIV of basis function:
          DynRankView ConstructWithLabel(vals, numFields, numPoints);
          tetBasis.getValues(vals, tetNodes, OPERATOR_DIV);
          const auto vals_host = Kokkos::create_mirror_view(typename HostSpaceType::memory_space(), vals);
          Kokkos::deep_copy(vals_host, vals);
          for (int i = 0; i < numFields; i++) {
            for (size_type j = 0; j < numPoints; j++) {
                int l =  i + j * numFields;
                 if (std::abs(vals_host(i,j) - basisDivs[l]) > tol) {
                   errorFlag++;
                   *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";

                   // Output the multi-index of the value where the error is:
                   *outStream << " At multi-index { ";
                   *outStream << i << " ";*outStream << j << " ";
                   *outStream << "}  computed divergence component: " << vals_host(i,j)
                     << " but reference divergence component: " << basisDivs[l] << "\n";
               }
            }
          }
        }
      } catch (std::logic_error err) {
        *outStream << err.what() << "\n\n";
        errorFlag = -1000;
      }

     *outStream
       << "\n"
       << "===============================================================================\n"
       << "| TEST 4: DOF correctness (Kronecker property)                                |\n"
       << "===============================================================================\n";

     try {
       const ordinal_type numFields = tetBasis.getCardinality();
       const ordinal_type spaceDim  = tetBasis.getBaseCellTopology().getDimension();

        // Check exceptions.
        ordinal_type nthrow = 0, ncatch = 0;
#ifdef HAVE_INTREPID2_DEBUG
        {
          DynRankView ConstructWithLabel(badVals, 1,2,3);
          INTREPID2_TEST_ERROR_EXPECTED( tetBasis.getDofCoords(badVals));
        }
        {
          DynRankView ConstructWithLabel(badVals, 3,2);
          INTREPID2_TEST_ERROR_EXPECTED( tetBasis.getDofCoords(badVals));
        }
        {
          DynRankView ConstructWithLabel(badVals, 4,2);
          INTREPID2_TEST_ERROR_EXPECTED( tetBasis.getDofCoords(badVals));
        }
#endif
        if (nthrow != ncatch) {
          errorFlag++;
          *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
          *outStream << "# of catch ("<< ncatch << ") is different from # of throw (" << ncatch << ")\n";
        }
       
        DynRankView ConstructWithLabel(bvals, numFields, numFields, spaceDim);
        DynRankView ConstructWithLabel(cvals, numFields, spaceDim);
        DynRankView ConstructWithLabel(dofCoeffs, numFields, spaceDim);
       
        // Check mathematical correctness.
        tetBasis.getDofCoords(cvals);
        tetBasis.getValues(bvals, cvals, OPERATOR_VALUE);
        tetBasis.getDofCoeffs(dofCoeffs);
    
        auto cvals_host = Kokkos::create_mirror_view(typename HostSpaceType::memory_space(), cvals);
        Kokkos::deep_copy(cvals_host, cvals);

        auto bvals_host = Kokkos::create_mirror_view(typename HostSpaceType::memory_space(), bvals);
        Kokkos::deep_copy(bvals_host, bvals);

        auto dofCoeffs_host = Kokkos::create_mirror_view(typename HostSpaceType::memory_space(), dofCoeffs);
        Kokkos::deep_copy(dofCoeffs_host, dofCoeffs);
    
        for (ordinal_type i=0;i<numFields;++i) {
          for (ordinal_type j=0;j<numFields;++j) {

            ValueType dofValue = 0.0;
            for(ordinal_type d=0;d<spaceDim;++d) {
              dofValue += bvals_host(i,j,d)*dofCoeffs_host(j,d);
            }

            const ValueType expected_dofValue = (i == j);
            if (std::abs(dofValue - expected_dofValue) > tol || std::isnan(dofValue)) {
              errorFlag++;
              std::stringstream ss;
              ss << "\nDegree of freedom " << j << " of basis function " << i << " at (" << cvals_host(j,0) << ", " << cvals_host(j,1)<< ", " << cvals_host(j,2) << ") is " << dofValue << " but should be " << expected_dofValue << "\n";
              *outStream << ss.str();
            }
          }
        }
      } catch (std::logic_error err){
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
  }
}

