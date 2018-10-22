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
    \brief  Unit tests for the Intrepid2::G_TRIPRISM_C1_FEM class.
    \author Created by P. Bochev, D. Ridzal, K. Peterson, M. Perego and Kyungjoo Kim.
*/

#include "Intrepid2_config.h"

#ifdef HAVE_INTREPID2_DEBUG
#define INTREPID2_TEST_FOR_DEBUG_ABORT_OVERRIDE_TO_CONTINUE
#endif

#include "Intrepid2_HGRAD_PYR_C1_FEM.hpp"

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
    int HGRAD_PYR_C1_FEM_Test01(const bool verbose) {

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
        << "|                 Unit Test (Basis_HGRAD_PYR_C1_FEM)                          |\n"
        << "|                                                                             |\n"
        << "|     1) Conversion of Dof tags into Dof ordinals and back                    |\n"
        << "|     2) Basis values for VALUE, GRAD, CURL, and Dk operators                 |\n"
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
      typedef Kokkos::DynRankView<ValueType,HostSpaceType>   DynRankViewHost;
#define ConstructWithLabel(obj, ...) obj(#obj, __VA_ARGS__)

      const ValueType tol = tolerence();
      int errorFlag = 0;

      // for virtual function, value and point types are declared in the class
      typedef ValueType outputValueType;
      typedef ValueType pointValueType;
      Basis_HGRAD_PYR_C1_FEM<DeviceSpaceType,outputValueType,pointValueType> pyrBasis;

     *outStream
       << "\n"
       << "===============================================================================\n"
       << "| TEST 1: constructors and exceptions                                         |\n"
       << "===============================================================================\n";

     try {
       ordinal_type nthrow = 0, ncatch = 0;
#ifdef HAVE_INTREPID2_DEBUG
       // Define array containing the 4 vertices of the reference PYR and its center.
       DynRankView ConstructWithLabel(pyrNodes, 10, 3);

       // Generic array for the output values; needs to be properly resized depending on the operator type
       const auto numFields = pyrBasis.getCardinality();
       const auto numPoints = pyrNodes.extent(0);
       const auto spaceDim  = pyrBasis.getBaseCellTopology().getDimension();

       DynRankView vals("vals", numFields, numPoints);
       DynRankView vals_vec("vals", numFields, numPoints, spaceDim);

       {
         // exception #1: CURL cannot be applied to scalar functions
         INTREPID2_TEST_ERROR_EXPECTED( pyrBasis.getValues(vals_vec, pyrNodes, OPERATOR_CURL) );
         // exception #2: DIV cannot be applied to scalar functions
         INTREPID2_TEST_ERROR_EXPECTED( pyrBasis.getValues(vals_vec, pyrNodes, OPERATOR_DIV) );
       }

       // Exceptions 3-7: all bf tags/bf Ids below are wrong and should cause getDofOrdinal() and
       // getDofTag() to access invalid array elements thereby causing bounds check exception
       {
         // exception #3
         INTREPID2_TEST_ERROR_EXPECTED( pyrBasis.getDofOrdinal(3,0,0) );
         // exception #4
         INTREPID2_TEST_ERROR_EXPECTED( pyrBasis.getDofOrdinal(1,1,1) );
         // exception #5
         INTREPID2_TEST_ERROR_EXPECTED( pyrBasis.getDofOrdinal(0,6,0) );
         // exception #6
         INTREPID2_TEST_ERROR_EXPECTED( pyrBasis.getDofTag(numFields) );
         // exception #7
         INTREPID2_TEST_ERROR_EXPECTED( pyrBasis.getDofTag(-1) );
       }

       // Exceptions 8-18 test exception handling with incorrectly dimensioned input/output arrays
       {
         // exception #8: input points array must be of rank-2
         DynRankView ConstructWithLabel(badPoints1, 4, 5, 3);
         INTREPID2_TEST_ERROR_EXPECTED( pyrBasis.getValues(vals, badPoints1, OPERATOR_VALUE) );
       }
       {
         // exception #9 dimension 1 in the input point array must equal space dimension of the cell
         DynRankView ConstructWithLabel(badPoints2, 4, 2);
         INTREPID2_TEST_ERROR_EXPECTED( pyrBasis.getValues(vals, badPoints2, OPERATOR_VALUE) );
       }
       {
         // exception #10 output values must be of rank-2 for OPERATOR_VALUE
         DynRankView ConstructWithLabel(badVals1, 4, 3, 1);
         INTREPID2_TEST_ERROR_EXPECTED( pyrBasis.getValues(badVals1, pyrNodes, OPERATOR_VALUE) );
       }
       {
         // exception #11 output values must be of rank-3 for OPERATOR_GRAD
         DynRankView ConstructWithLabel(badVals2, 4, 3);
         INTREPID2_TEST_ERROR_EXPECTED( pyrBasis.getValues(badVals2, pyrNodes, OPERATOR_GRAD) );

         // exception #12 output values must be of rank-3 for OPERATOR_D1
         INTREPID2_TEST_ERROR_EXPECTED( pyrBasis.getValues(badVals2, pyrNodes, OPERATOR_D1) );

         // exception #13 output values must be of rank-3 for OPERATOR_D2
         INTREPID2_TEST_ERROR_EXPECTED( pyrBasis.getValues(badVals2, pyrNodes, OPERATOR_D2) );
       }
       {
         // exception #14 incorrect 0th dimension of output array (must equal number of basis functions)
         DynRankView ConstructWithLabel(badVals3, numFields + 1, numPoints);
         INTREPID2_TEST_ERROR_EXPECTED( pyrBasis.getValues(badVals3, pyrNodes, OPERATOR_VALUE) );
       }
       {
         // exception #15 incorrect 1st dimension of output array (must equal number of points)
         DynRankView ConstructWithLabel(badVals4, numFields, numPoints + 1);
         INTREPID2_TEST_ERROR_EXPECTED( pyrBasis.getValues(badVals4, pyrNodes, OPERATOR_VALUE) );
       }
       {
         // exception #16: incorrect 2nd dimension of output array (must equal the space dimension)
         DynRankView ConstructWithLabel(badVals5, numFields, numPoints, 4);
         INTREPID2_TEST_ERROR_EXPECTED( pyrBasis.getValues(badVals5, pyrNodes, OPERATOR_GRAD) );
       }
       {
         // exception #17: incorrect 2nd dimension of output array (must equal D2 cardinality in 3D)
         DynRankView ConstructWithLabel(badVals6, numFields, numPoints, 40);
         INTREPID2_TEST_ERROR_EXPECTED( pyrBasis.getValues(badVals6, pyrNodes, OPERATOR_D2) );
         //exception #18: incorrect 2nd dimension of output array (must equal D3 cardinality in 3D)
         INTREPID2_TEST_ERROR_EXPECTED( pyrBasis.getValues(badVals6, pyrNodes, OPERATOR_D3) );
       }
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

     *outStream                                 \
       << "\n"
       << "===============================================================================\n" \
       << "| TEST 2: correctness of tag to enum and enum to tag lookups                  |\n" \
       << "===============================================================================\n";

     try {
       const ordinal_type numFields = pyrBasis.getCardinality();
       const auto allTags = pyrBasis.getAllDofTags();

       // Loop over all tags, lookup the associated dof enumeration and then lookup the tag again
       const ordinal_type dofTagSize = allTags.extent(0);
       for (ordinal_type i=0;i<dofTagSize;++i) {
         const auto bfOrd = pyrBasis.getDofOrdinal(allTags(i,0), allTags(i,1), allTags(i,2));

         const auto myTag = pyrBasis.getDofTag(bfOrd);
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
       for (ordinal_type bfOrd=0;bfOrd<numFields;++bfOrd) {
         const auto myTag = pyrBasis.getDofTag(bfOrd);

         const auto myBfOrd = pyrBasis.getDofOrdinal(myTag(0), myTag(1), myTag(2));
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

     *outStream                                 \
       << "\n"
       << "===============================================================================\n" \
       << "| TEST 3: correctness of basis function values                                |\n" \
       << "===============================================================================\n";

     outStream -> precision(20);

     // VALUE: Each row gives the 4 correct basis set values at an evaluation point
     const ValueType basisValues[] = {
       1.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 1.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 1.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 1.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 1.0,
       //
       0.0515625,   0.0984375,  0.4265625,  0.2234375,  0.2,
       0.2,         0,          0.,         0.5,        0.3,
       0.025,       0.025,      0.,         0,          0.95,
       0.18,        0.045,      0.005,      0.02,       0.75,
       0.035,       0.015,      0.285,      0.665,      0.,
     };

     // GRAD and D1: each row gives the 3 x 5 correct values of the gradients of the 5 basis functions
     const ValueType basisGrads[] = {
       -0.5, -0.5,  0.0,  0.5,  0.0, -0.5,  0.0,  0.0,  0.0,  0.0,  0.5, -0.5,  0.0,  0.0,  1.0, \
       -0.5,  0.0, -0.5,  0.5, -0.5,  0.0,  0.0,  0.5, -0.5,  0.0,  0.0,  0.0,  0.0,  0.0,  1.0, \
       0.0,  0.0,  0.0,  0.0, -0.5, -0.5,  0.5,  0.5,  0.0, -0.5,  0.0, -0.5,  0.0,  0.0,  1.0, \
       0.0, -0.5, -0.5,  0.0,  0.0,  0.0,  0.5,  0.0, -0.5, -0.5,  0.5,  0.0,  0.0,  0.0,  1.0, \
       -0.25,-0.25,-0.25, 0.25,-0.25,-0.25, 0.25, 0.25,-0.25,-0.25, 0.25,-0.25, 0.0,  0.0,  1.0, \
       -0.09375, -0.171875, -0.201171875,  0.09375, -0.328125, -0.298828125, 0.40625,  0.328125, -0.201171875, -0.40625,  0.171875, -0.298828125,  0.0,  0.0,  1.0, \
       -0.1428571428571429, -0.5, -0.3571428571428571,  0.1428571428571429,  0.0, -0.1428571428571429,  0.3571428571428571, 0.0, -0.3571428571428571, -0.3571428571428571,  0.5, -0.1428571428571429,  0.0,  0.0,  1.0, \
       -0.5, -0.25, -0.25,  0.5, -0.25, -0.25,  0.0,  0.25, -0.25,  0.,  0.25, -0.25, 0.0,  0.0,  1.0, \
       -0.45, -0.4, -0.13,  0.45, -0.1, -0.37,  0.05,  0.1, -0.13, -0.05,  0.4, -0.37,  0.0,  0.0,  1.0, \
       -0.025, -0.35, -0.34,  0.025, -0.15, -0.16,  0.475,  0.15, -0.34, -0.475,  0.35, -0.16,  0.0,  0.0,  1.0
     };


     //D2: flat array with the values of D2 applied to basis functions. Multi-index is (P,F,K)
     const auto eps = epsilon();
     const ValueType basisD2[] = {
       0, 0.25,-0.25, 0,-0.25, 0.5, 0,-0.25, 0.25, 0, 0.25,-0.5, 0, 0.25,-0.25, 0,-0.25, 0.5, 0,-0.25, 0.25, 0, 0.25,-0.5, 0, 0, 0, 0, 0, 0, \
       0, 0.25,-0.25, 0, 0.25,-0.5, 0,-0.25, 0.25, 0,-0.25, 0.5, 0, 0.25,-0.25, 0, 0.25,-0.5, 0,-0.25, 0.25, 0,-0.25, 0.5, 0, 0, 0, 0, 0, 0, \
       0, 0.25, 0.25, 0, 0.25, 0.5, 0,-0.25,-0.25, 0,-0.25,-0.5, 0, 0.25, 0.25, 0, 0.25, 0.5, 0,-0.25,-0.25, 0,-0.25,-0.5, 0, 0, 0, 0, 0, 0, \
       0, 0.25, 0.25, 0,-0.25,-0.5, 0,-0.25,-0.25, 0, 0.25, 0.5, 0, 0.25, 0.25, 0,-0.25,-0.5, 0,-0.25,-0.25, 0, 0.25, 0.5, 0, 0, 0, 0, 0, 0, \
       0, 0.25/eps, 0, 0, 0, 0, 0,-0.25/eps, 0, 0, 0, 0, 0, 0.25/eps, 0, 0, 0, 0, 0,-0.25/eps, 0, 0, 0, 0, 0, 0, 0, 0, 0,0, \
       0, 0.3125, 0.1953125, 0, 0.09765625, 0.1220703125, 0,-0.3125,-0.1953125, 0,-0.09765625,-0.1220703125, 0, 0.3125, 0.1953125, 0, 0.09765625, 0.1220703125, 0,-0.3125,-0.1953125, 0,-0.09765625,-0.1220703125, 0, 0, 0, 0, 0, 0, \
       0, 0.3571428571428571, 0.1530612244897959, 0,-0.3571428571428571,-0.306122448979592, 0,-0.3571428571428572,-0.1530612244897959, 0, 0.3571428571428571, 0.306122448979592, 0, 0.3571428571428571, 0.1530612244897959, 0,-0.3571428571428571,-0.306122448979592, 0,-0.3571428571428571,-0.1530612244897959, 0, 0.3571428571428571, 0.306122448979592, 0, 0, 0, 0, 0, 0, \
       0, 5,-5, 0, 0, 0, 0,-5, 5, 0, 0, 0, 0, 5,-5, 0, 0, 0, 0, -5, 5, 0, 0, 0, 0, 0, 0, 0, 0, 0, \
       0, 1,-0.8, 0,-0.6, 0.96, 0, -1, 0.8, 0, 0.6,-0.96, 0,1,-0.8, 0, -0.6, 0.96, 0, -1, 0.8, 0, 0.6, -0.96, 0, 0, 0, 0, 0, 0, \
       0, 0.25, 0.225, 0,-0.1,-0.18, 0,-0.25,-0.225, 0, 0.1, 0.18, 0, 0.25, 0.225, 0,-0.1,-0.18, 0,-0.25,-0.225,0,0.1,0.18, 0, 0, 0, 0, 0, 0
     };

     try {
       DynRankViewHost ConstructWithLabel(pyrNodesHost, 10, 3);

       pyrNodesHost(0,0) = -1.0;  pyrNodesHost(0,1) = -1.0;  pyrNodesHost(0,2) =  0;
       pyrNodesHost(1,0) =  1.0;  pyrNodesHost(1,1) = -1.0;  pyrNodesHost(1,2) =  0;
       pyrNodesHost(2,0) =  1.0;  pyrNodesHost(2,1) =  1.0;  pyrNodesHost(2,2) =  0;
       pyrNodesHost(3,0) = -1.0;  pyrNodesHost(3,1) =  1.0;  pyrNodesHost(3,2) =  0;
       pyrNodesHost(4,0) =  0.0;  pyrNodesHost(4,1) =  0.0;  pyrNodesHost(4,2) =  1.0;

       pyrNodesHost(5,0) =  0.25; pyrNodesHost(5,1) =  0.5;  pyrNodesHost(5,2) = 0.2;
       pyrNodesHost(6,0) = -0.7 ; pyrNodesHost(6,1) =  0.3;  pyrNodesHost(6,2) = 0.3;
       pyrNodesHost(7,0) =  0.;   pyrNodesHost(7,1) = -0.05; pyrNodesHost(7,2) = 0.95;
       pyrNodesHost(8,0) = -0.15; pyrNodesHost(8,1) = -0.2;  pyrNodesHost(8,2) = 0.75;
       pyrNodesHost(9,0) = -0.4;  pyrNodesHost(9,1) =  0.9;  pyrNodesHost(9,2) = 0.0;

       auto pyrNodes = Kokkos::create_mirror_view(typename DeviceSpaceType::memory_space(), pyrNodesHost);
       Kokkos::deep_copy(pyrNodes, pyrNodesHost);

       // Dimensions for the output arrays:
       const ordinal_type numFields = pyrBasis.getCardinality();
       const ordinal_type numPoints = pyrNodes.extent(0);
       const ordinal_type spaceDim  = pyrBasis.getBaseCellTopology().getDimension();
       const ordinal_type D2Cardin  = getDkCardinality(OPERATOR_D2, spaceDim);

       // Check VALUE of basis functions: resize vals to rank-2 container:
       {
         DynRankView vals = DynRankView("vals", numFields, numPoints);
         pyrBasis.getValues(vals, pyrNodes, OPERATOR_VALUE);
         auto vals_host = Kokkos::create_mirror_view(typename HostSpaceType::memory_space(), vals);
         Kokkos::deep_copy(vals_host, vals);
         for (ordinal_type i=0;i<numFields;++i) {
           for (ordinal_type j=0;j<numPoints;++j) {
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

       // Check GRAD of basis function: resize vals to rank-3 container
       {
         DynRankView vals = DynRankView("vals", numFields, numPoints, spaceDim);
         pyrBasis.getValues(vals, pyrNodes, OPERATOR_GRAD);
         auto vals_host = Kokkos::create_mirror_view(typename HostSpaceType::memory_space(), vals);
         Kokkos::deep_copy(vals_host, vals);
         for (ordinal_type i=0;i<numFields;++i) {
           for (ordinal_type j=0;j<numPoints;++j) {
             for (ordinal_type k=0;k<spaceDim;++k) {
               const ordinal_type l = k + i * spaceDim + j * spaceDim * numFields;
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
       }

       // Check D1 of basis function
       {
         DynRankView vals = DynRankView("vals", numFields, numPoints, spaceDim);
         pyrBasis.getValues(vals, pyrNodes, OPERATOR_D1);
         auto vals_host = Kokkos::create_mirror_view(typename HostSpaceType::memory_space(), vals);
         Kokkos::deep_copy(vals_host, vals);
         for (ordinal_type i=0;i<numFields;++i) {
           for (ordinal_type j=0;j<numPoints;++j) {
             for (ordinal_type k=0;k<spaceDim;++k) {
               const ordinal_type l = k + i * spaceDim + j * spaceDim * numFields;
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

       // Check D2 of basis function
       {
         DynRankView vals = DynRankView("vals", numFields, numPoints, D2Cardin);
         pyrBasis.getValues(vals, pyrNodes, OPERATOR_D2);
         auto vals_host = Kokkos::create_mirror_view(typename HostSpaceType::memory_space(), vals);
         Kokkos::deep_copy(vals_host, vals);
         for (ordinal_type i=0;i<numFields;++i) {
           for (ordinal_type j=0;j<numPoints;++j) {
             // derivatives are singular when z = 1; using the same eps, it can be comparable
             //if (j == 4) continue; 
             for (ordinal_type k=0;k<D2Cardin;++k) {
               const ordinal_type l = k + i * D2Cardin + j * D2Cardin * numFields;
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

       // Check all higher derivatives - must be zero.
       {
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
           const ordinal_type DkCardin  = getDkCardinality(op, spaceDim);
           DynRankView vals("vals", numFields, numPoints, DkCardin);

           pyrBasis.getValues(vals, pyrNodes, op);
           auto vals_host = Kokkos::create_mirror_view(typename HostSpaceType::memory_space(), vals);
           Kokkos::deep_copy(vals_host, vals);
           for (ordinal_type i1=0;i1<numFields; i1++)
             for (ordinal_type i2=0;i2<numPoints; i2++)
               for (ordinal_type i3=0;i3<DkCardin; i3++) {
                 if (std::abs(vals_host(i1,i2,i3)) > tol) {
                   errorFlag++;
                   *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";

                   // Get the multi-index of the value where the error is and the operator order
                   const ordinal_type ord = Intrepid2::getOperatorOrder(op);
                   *outStream << " At multi-index { "<<i1<<" "<<i2 <<" "<<i3;
                   *outStream << "}  computed D"<< ord <<" component: " << vals(i1,i2,i3)
                              << " but reference D" << ord << " component:  0 \n";
                 }
               }
         }
       }
     } catch (std::logic_error err) {
       *outStream << err.what() << "\n\n";
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
  }
}






