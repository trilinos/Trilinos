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

/** \file test_01.cpp
\brief  Unit tests for the Intrepid2::HDIV_TRI_I1_FEM class.
\author Created by P. Bochev, D. Ridzal, and K. Peterson.
*/

#include "Intrepid2_config.h"

#ifdef HAVE_INTREPID2_DEBUG
#define INTREPID2_TEST_FOR_DEBUG_ABORT_OVERRIDE_TO_CONTINUE
#endif

#include "Intrepid2_Basis.hpp"
#include "Intrepid2_HDIV_TRI_I1_FEM.hpp"

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
    int HDIV_TRI_I1_FEM_Test01(const bool verbose) {

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
        << "|                 Unit Test (Basis_HDIV_TRI_I1_FEM)                           |\n"
        << "|                                                                             |\n"
        << "|     1) Conversion of Dof tags into Dof ordinals and back                    |\n"
        << "|     2) Basis values for VALUE and DIV operators                             |\n"
        << "|                                                                             |\n"
        << "|  Questions? Contact  Pavel Bochev  (pbboche@sandia.gov),                    |\n"
        << "|                      Denis Ridzal  (dridzal@sandia.gov),                    |\n"
        << "|                      Kara Peterson (dridzal@sandia.gov).                    |\n"
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
      Basis_HDIV_TRI_I1_FEM<DeviceSpaceType,outputValueType,pointValueType> triBasis;
      //typedef typename decltype(triBasis)::outputViewType outputViewType;
      //typedef typename decltype(triBasis)::pointViewType  pointViewType;

      *outStream
        << "\n"
        << "===============================================================================\n"
        << "| TEST 1: Basis creation, exceptions tests                                    |\n"
        << "===============================================================================\n";

      try {
        ordinal_type nthrow = 0, ncatch = 0;
#ifdef HAVE_INTREPID2_DEBUG

        // Define array containing the 3 vertices of the reference TRI and its 3 edge midpoints.
        DynRankView ConstructWithLabel(triNodes, 6, 2);

        // Generic array for the output values; needs to be properly resized depending on the operator type
        const ordinal_type numFields = triBasis.getCardinality();
        const ordinal_type numPoints = triNodes.extent(0);

        DynRankView vals;
        vals = DynRankView("vals", numFields, numPoints);


        // exception #1: GRAD cannot be applied to HDIV functions
        // resize vals to rank-3 container with dimensions (num. basis functions, num. points, arbitrary)
        INTREPID2_TEST_ERROR_EXPECTED( triBasis.getValues(vals, triNodes, OPERATOR_GRAD) );

        // exception #2: CURL cannot be applied to HDIV functions
        INTREPID2_TEST_ERROR_EXPECTED( triBasis.getValues(vals, triNodes, OPERATOR_CURL) );

        // Exceptions 3-7: all bf tags/bf Ids below are wrong and should cause getDofOrdinal() and
        // getDofTag() to access invalid array elements thereby causing bounds check exception
        // exception #3
        INTREPID2_TEST_ERROR_EXPECTED( triBasis.getDofOrdinal(3,0,0) );
        // exception #4
        INTREPID2_TEST_ERROR_EXPECTED( triBasis.getDofOrdinal(1,1,1) );
        // exception #5
        INTREPID2_TEST_ERROR_EXPECTED( triBasis.getDofOrdinal(0,2,1) );
        // exception #6
        INTREPID2_TEST_ERROR_EXPECTED( triBasis.getDofTag(numFields) );
        // exception #7
        INTREPID2_TEST_ERROR_EXPECTED( triBasis.getDofTag(-1) );
        // exception #8
        INTREPID2_TEST_ERROR_EXPECTED( triBasis.getDofOrdinal(2,0,0) );

        // Exceptions 9-16 test exception handling with incorrectly dimensioned input/output arrays
        // exception #9: input points array must be of rank-2
        DynRankView ConstructWithLabel(badPoints1, 4, 5, 3);
        INTREPID2_TEST_ERROR_EXPECTED( triBasis.getValues(vals, badPoints1, OPERATOR_VALUE) );

        // exception #10 dimension 1 in the input point array must equal space dimension of the cell
        DynRankView ConstructWithLabel(badPoints2, 4, triBasis.getBaseCellTopology().getDimension() + 1);
        INTREPID2_TEST_ERROR_EXPECTED( triBasis.getValues(vals, badPoints2, OPERATOR_VALUE) );
        
        // exception #11 output values must be of rank-3 for OPERATOR_VALUE
        DynRankView ConstructWithLabel(badVals1, 4, 3);
        INTREPID2_TEST_ERROR_EXPECTED( triBasis.getValues(badVals1, triNodes, OPERATOR_VALUE) );

        // exception #12 output values must be of rank-2 for OPERATOR_DIV
        DynRankView ConstructWithLabel(badVals2, 4, 3, 1);
        INTREPID2_TEST_ERROR_EXPECTED( triBasis.getValues(badVals2, triNodes, OPERATOR_VALUE) );

        // exception #13 incorrect 0th dimension of output array for OPERATOR_VALUE (must equal number of basis functions)
        DynRankView ConstructWithLabel(badVals3, triBasis.getCardinality() + 1, triNodes.extent(0), triBasis.getBaseCellTopology().getDimension());
        INTREPID2_TEST_ERROR_EXPECTED( triBasis.getValues(badVals3, triNodes, OPERATOR_VALUE) );

        // exception #14 incorrect 0th dimension of output array for OPERATOR_DIV (must equal number of basis functions)
        DynRankView ConstructWithLabel(badVals4, triBasis.getCardinality() + 1, triNodes.extent(0));
        INTREPID2_TEST_ERROR_EXPECTED( triBasis.getValues(badVals4, triNodes, OPERATOR_DIV) );

        // exception #15 incorrect 1st dimension of output array (must equal number of points)
        DynRankView ConstructWithLabel(badVals5, triBasis.getCardinality(), triNodes.extent(0) + 1, triBasis.getBaseCellTopology().getDimension());
        INTREPID2_TEST_ERROR_EXPECTED( triBasis.getValues(badVals5, triNodes, OPERATOR_VALUE) );
    
        // exception #16 incorrect 1st dimension of output array (must equal number of points)
        DynRankView ConstructWithLabel(badVals6, triBasis.getCardinality(), triNodes.extent(0) + 1);
        INTREPID2_TEST_ERROR_EXPECTED( triBasis.getValues(badVals6, triNodes, OPERATOR_DIV) );
    
        // exception #17: incorrect 2nd dimension of output array (must equal the space dimension)
        DynRankView ConstructWithLabel(badVals7, triBasis.getCardinality(), triNodes.extent(0), triBasis.getBaseCellTopology().getDimension() + 1);
        INTREPID2_TEST_ERROR_EXPECTED( triBasis.getValues(badVals7, triNodes, OPERATOR_VALUE) );
    
#endif
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
      try {
        const ordinal_type numFields = triBasis.getCardinality();
        const auto allTags = triBasis.getAllDofTags();
   
        // Loop over all tags, lookup the associated dof enumeration and then lookup the tag again
        const ordinal_type dofTagSize = allTags.extent(0);
        
        // Loop over all tags, lookup the associated dof enumeration and then lookup the tag again
        for (ordinal_type i = 0; i < dofTagSize; i++) {
          const auto bfOrd  = triBasis.getDofOrdinal(allTags(i,0), allTags(i,1), allTags(i,2));

          const auto myTag = triBasis.getDofTag(bfOrd);
          if( !( (myTag(0) == allTags(i, 0)) &&
                  (myTag(1) == allTags(i, 1)) &&
                  (myTag(2) == allTags(i, 2)) &&
                  (myTag(3) == allTags(i, 3)) ) ) {
            errorFlag++;
            *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
            *outStream << " getDofOrdinal( {"
              << allTags(i, 0) << ", "
              << allTags(i, 1) << ", "
              << allTags(i, 2) << ", "
              << allTags(i, 3) << "}) = " << bfOrd <<" but \n";
            *outStream << " getDofTag(" << bfOrd << ") = { "
              << myTag(0) << ", "
              << myTag(1) << ", "
              << myTag(2) << ", "
              << myTag(3) << "}\n";
          }
        }
    
        // Now do the same but loop over basis functions
        for( ordinal_type bfOrd = 0; bfOrd < numFields; bfOrd++) {
          auto myTag  = triBasis.getDofTag(bfOrd);
          auto myBfOrd = triBasis.getDofOrdinal(myTag[0], myTag[1], myTag[2]);
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

      *outStream \
        << "\n"
        << "===============================================================================\n"\
        << "| TEST 3: correctness of basis function values                                |\n"\
        << "===============================================================================\n";

      outStream -> precision(20);

      // VALUE: Correct basis values in (F,P,D) format: each group of two rows gives basis function
      // values at vertices followed by midpoints. This is the same array format as the output from getValues.
      double basisValues[] = {
        // basis function 0 at 3 vertices followed by 3 midpoints
        0.0,-1.0,                  1.0,-1.0,                               0.0, 0.0,
        0.5,-1.0,                  0.5,-0.5,                               0.0,-0.5,
        // basis function 1 at 3 vertices followed by 3 midpoints
        0.0, 0.0,                  1.0, 0.0,                0.0, 1.0,
        0.5, 0.0,                  0.5, 0.5,                0.0, 0.5,
        // basis function 2 at 3 vertices followed by 3 midpoints
        -1.0, 0.0,                 0.0, 0.0,                              -1.0, 1.0,
        -0.5, 0.0,                -0.5, 0.5,                              -1.0, 0.5
      };

      // DIV: each row gives the 3 correct values of the divergence of the 3 basis functions
      double basisDivs[] = {
        // 3 vertices
        2.0,  2.0,   2.0,
        2.0,  2.0,   2.0,
        2.0,  2.0,   2.0,
        // 3 edge centers
        2.0,  2.0,   2.0,
        2.0,  2.0,   2.0,
        2.0,  2.0,   2.0,
      };



      try{

        DynRankViewHost ConstructWithLabel(triNodesHost, 6, 2);

        triNodesHost(0,0) =  0.0;  triNodesHost(0,1) =  0.0;
        triNodesHost(1,0) =  1.0;  triNodesHost(1,1) =  0.0;
        triNodesHost(2,0) =  0.0;  triNodesHost(2,1) =  1.0;
        // edge midpoints
        triNodesHost(3,0) =  0.5;  triNodesHost(3,1) =  0.0;
        triNodesHost(4,0) =  0.5;  triNodesHost(4,1) =  0.5;
        triNodesHost(5,0) =  0.0;  triNodesHost(5,1) =  0.5;

        auto triNodes = Kokkos::create_mirror_view(typename DeviceSpaceType::memory_space(), triNodesHost);
        Kokkos::deep_copy(triNodes, triNodesHost);
        
        // Dimensions for the output arrays:
        const ordinal_type numFields = triBasis.getCardinality();
        const ordinal_type numPoints = triNodes.extent(0);
        const ordinal_type spaceDim  = triBasis.getBaseCellTopology().getDimension();
    
        {
          // Check VALUE of basis functions: resize vals to rank-3 container:
          DynRankView vals = DynRankView("vals", numFields, numPoints, spaceDim);
          triBasis.getValues(vals, triNodes, OPERATOR_VALUE);
          const auto vals_host = Kokkos::create_mirror_view(typename HostSpaceType::memory_space(), vals);
          Kokkos::deep_copy(vals_host, vals);
          for (ordinal_type i = 0; i < numFields; i++) {
            for (ordinal_type j = 0; j < numPoints; j++) {
              for (ordinal_type k = 0; k < spaceDim; k++) {
                // basisValues are in (F,P,D) format and the multiindex is (i,j,k), here's the offset:
                 ordinal_type l = k + j * spaceDim + i * spaceDim * numPoints;

                 if (std::abs(vals_host(i,j,k) - basisValues[l]) > tol) {
                   errorFlag++;
                   *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";

                   // Output the multi-index of the value where the error is:
                   *outStream << " address =  "<< l <<"\n";
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
          // Check DIV of basis function: resize vals to rank-2 container
          DynRankView vals = DynRankView("vals", numFields, numPoints);
          triBasis.getValues(vals, triNodes, OPERATOR_DIV);
          const auto vals_host = Kokkos::create_mirror_view(typename HostSpaceType::memory_space(), vals);
          Kokkos::deep_copy(vals_host, vals);
          for (ordinal_type i = 0; i < numFields; i++) {
            for (ordinal_type j = 0; j < numPoints; j++) {
              ordinal_type l =  i + j * numFields;
              if (std::abs(vals_host(i,j) - basisDivs[l]) > tol) {
                errorFlag++;
                *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";

                // Output the multi-index of the value where the error is:
                *outStream << " At (Field,Point,Dim) multi-index { ";
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
        const ordinal_type numFields = triBasis.getCardinality();
        const ordinal_type spaceDim  = triBasis.getBaseCellTopology().getDimension();

         // Check exceptions.
         ordinal_type nthrow = 0, ncatch = 0;
 #ifdef HAVE_INTREPID2_DEBUG
         {
           DynRankView ConstructWithLabel(badVals, 1,2,3);
           INTREPID2_TEST_ERROR_EXPECTED( triBasis.getDofCoords(badVals));
         }
         {
           DynRankView ConstructWithLabel(badVals, 4,2);
           INTREPID2_TEST_ERROR_EXPECTED( triBasis.getDofCoords(badVals));
         }
         {
           DynRankView ConstructWithLabel(badVals, 3,3);
           INTREPID2_TEST_ERROR_EXPECTED( triBasis.getDofCoords(badVals));
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
         triBasis.getDofCoords(cvals);
         triBasis.getValues(bvals, cvals, OPERATOR_VALUE);
         triBasis.getDofCoeffs(dofCoeffs);


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
