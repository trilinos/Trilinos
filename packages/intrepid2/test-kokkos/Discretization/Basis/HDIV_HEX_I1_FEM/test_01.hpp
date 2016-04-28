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
\brief  Unit tests for the Intrepid2::C_HEX_I1_FEM class.
\author Created by P. Bochev, D. Ridzal, and K. Peterson.
*/
#include "Intrepid2_config.h"

#ifdef HAVE_INTREPID2_DEBUG
#define INTREPID2_TEST_FOR_DEBUG_ABORT_OVERRIDE_TO_CONTINUE
#endif

#include "Intrepid2_HDIV_HEX_I1_FEM.hpp"

#include "Teuchos_oblackholestream.hpp"
#include "Teuchos_RCP.hpp"

namespace Intrepid2 {

  namespace Test {

#define INTREPID2_TEST_ERROR_EXPECTED( S )              \
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
    int HDIV_HEX_I1_FEM_Test01(const bool verbose) {

      typedef ValueType value_type;

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
      << "===============================================================================\n" \
      << "|                                                                             |\n" \
      << "|                 Unit Test (Basis_HDIV_HEX_I1_FEM)                           |\n" \
      << "|                                                                             |\n" \
      << "|     1) Conversion of Dof tags into Dof ordinals and back                    |\n" \
      << "|     2) Basis values for VALUE and DIV operators                             |\n" \
      << "|                                                                             |\n" \
      << "|  Questions? Contact  Pavel Bochev  (pbboche@sandia.gov),                    |\n" \
      << "|                      Denis Ridzal  (dridzal@sandia.gov),                    |\n" \
      << "|                      Kara Peterson (kjpeter@sandia.gov).                    |\n" \
      << "|                                                                             |\n" \
      << "|  Intrepid's website: http://trilinos.sandia.gov/packages/intrepid           |\n" \
      << "|  Trilinos website:   http://trilinos.sandia.gov                             |\n" \
      << "|                                                                             |\n" \
      << "===============================================================================\n";

      typedef Kokkos::DynRankView<value_type,DeviceSpaceType> DynRankView;
#define ConstructWithLabel(obj, ...) obj(#obj, __VA_ARGS__)
      const value_type tol = Parameters::Tolerence;
      int errorFlag = 0;

      *outStream
        << "\n"
        << "===============================================================================\n"
        << "| TEST 1: Basis creation, exceptions tests                                    |\n"
        << "===============================================================================\n";

      try{
        ordinal_type nthrow = 0, ncatch = 0;
#ifdef HAVE_INTREPID2_DEBUG
        Basis_HDIV_HEX_I1_FEM<DeviceSpaceType> hexBasis;

        // Define array containing the 8 vertices of the reference HEX, its center and 6 face centers
        DynRankView ConstructWithLabel(hexNodes, 15, 3);

        int spaceDim  = hexBasis.getBaseCellTopology().getDimension();

        hexNodes(0,0) = -1.0;  hexNodes(0,1) = -1.0;  hexNodes(0,2) = -1.0;
        hexNodes(1,0) =  1.0;  hexNodes(1,1) = -1.0;  hexNodes(1,2) = -1.0;
        hexNodes(2,0) =  1.0;  hexNodes(2,1) =  1.0;  hexNodes(2,2) = -1.0;
        hexNodes(3,0) = -1.0;  hexNodes(3,1) =  1.0;  hexNodes(3,2) = -1.0;

        hexNodes(4,0) = -1.0;  hexNodes(4,1) = -1.0;  hexNodes(4,2) =  1.0;
        hexNodes(5,0) =  1.0;  hexNodes(5,1) = -1.0;  hexNodes(5,2) =  1.0;
        hexNodes(6,0) =  1.0;  hexNodes(6,1) =  1.0;  hexNodes(6,2) =  1.0;
        hexNodes(7,0) = -1.0;  hexNodes(7,1) =  1.0;  hexNodes(7,2) =  1.0;

        hexNodes(8,0) =  0.0;  hexNodes(8,1) =  0.0;  hexNodes(8,2) =  0.0;

        hexNodes(9,0) =  1.0;  hexNodes(9,1) =  0.0;  hexNodes(9,2) =  0.0;
        hexNodes(10,0)= -1.0;  hexNodes(10,1)=  0.0;  hexNodes(10,2)=  0.0;
        
        hexNodes(11,0)=  0.0;  hexNodes(11,1)=  1.0;  hexNodes(11,2)=  0.0;
        hexNodes(12,0)=  0.0;  hexNodes(12,1)= -1.0;  hexNodes(12,2)=  0.0;
        
        hexNodes(13,0)=  0.0;  hexNodes(13,1)=  0.0;  hexNodes(13,2)=  1.0;
        hexNodes(14,0)=  0.0;  hexNodes(14,1)=  0.0;  hexNodes(14,2)= -1.0;

        // exception #1: GRAD cannot be applied to HDIV functions
        // resize vals to rank-3 container with dimensions (num. basis functions, num. points, arbitrary)
        DynRankView ConstructWithLabel(vals, hexBasis.getCardinality(), hexNodes.dimension(0), spaceDim );
        INTREPID2_TEST_ERROR_EXPECTED( hexBasis.getValues(vals, hexNodes, OPERATOR_GRAD) );

        // exception #2: CURL cannot be applied to HDIV functions
        INTREPID2_TEST_ERROR_EXPECTED( hexBasis.getValues(vals, hexNodes, OPERATOR_CURL) );

        // Exceptions 3-7: all bf tags/bf Ids below are wrong and should cause getDofOrdinal() and
        // getDofTag() to access invalid array elements thereby causing bounds check exception
        // exception #3
        INTREPID2_TEST_ERROR_EXPECTED( hexBasis.getDofOrdinal(3,0,0) );
        // exception #4
        INTREPID2_TEST_ERROR_EXPECTED( hexBasis.getDofOrdinal(1,1,1) );
        // exception #5
        INTREPID2_TEST_ERROR_EXPECTED( hexBasis.getDofOrdinal(0,4,1) );
        // exception #6
        INTREPID2_TEST_ERROR_EXPECTED( hexBasis.getDofTag(12) );
        // exception #7
        INTREPID2_TEST_ERROR_EXPECTED( hexBasis.getDofTag(-1) );

        // Exceptions 8- test exception handling with incorrectly dimensioned input/output arrays
        // exception #8: input points array must be of rank-2
        DynRankView ConstructWithLabel(badPoints1, 4, 5, 3);
        INTREPID2_TEST_ERROR_EXPECTED( hexBasis.getValues(vals, badPoints1, OPERATOR_VALUE) );

        // exception #9 dimension 1 in the input point array must equal space dimension of the cell
        DynRankView ConstructWithLabel(badPoints2, 4, 2);
        INTREPID2_TEST_ERROR_EXPECTED( hexBasis.getValues(vals, badPoints2, OPERATOR_VALUE) );

        // exception #10 output values must be of rank-3 for OPERATOR_VALUE
        DynRankView ConstructWithLabel(badVals1, 4, 3);
        INTREPID2_TEST_ERROR_EXPECTED( hexBasis.getValues(badVals1, hexNodes, OPERATOR_VALUE) );

        // exception #11 output values must be of rank-2 for OPERATOR_DIV
        DynRankView ConstructWithLabel(badVals2, 4, 3, 3);
        INTREPID2_TEST_ERROR_EXPECTED( hexBasis.getValues(badVals2, hexNodes, OPERATOR_DIV) );

        // exception #12 incorrect 0th dimension of output array (must equal number of basis functions)
        DynRankView ConstructWithLabel(badVals3, hexBasis.getCardinality() + 1, hexNodes.dimension(0), 3);
        INTREPID2_TEST_ERROR_EXPECTED( hexBasis.getValues(badVals3, hexNodes, OPERATOR_VALUE) );

        // exception #13 incorrect 0th dimension of output array (must equal number of basis functions)
        DynRankView ConstructWithLabel(badVals4, hexBasis.getCardinality() + 1, hexNodes.dimension(0));
        INTREPID2_TEST_ERROR_EXPECTED( hexBasis.getValues(badVals4, hexNodes, OPERATOR_DIV) );
    
        // exception #14 incorrect 1st dimension of output array (must equal number of points)
        DynRankView ConstructWithLabel(badVals5, hexBasis.getCardinality(), hexNodes.dimension(0) + 1, 3);
        INTREPID2_TEST_ERROR_EXPECTED( hexBasis.getValues(badVals5, hexNodes, OPERATOR_VALUE) );
    
        // exception #15 incorrect 1st dimension of output array (must equal number of points)
        DynRankView ConstructWithLabel(badVals6, hexBasis.getCardinality(), hexNodes.dimension(0) + 1);
        INTREPID2_TEST_ERROR_EXPECTED( hexBasis.getValues(badVals6, hexNodes, OPERATOR_DIV) );
    
        // exception #16: incorrect 2nd dimension of output array (must equal the space dimension)
        DynRankView ConstructWithLabel(badVals7, hexBasis.getCardinality(), hexNodes.dimension(0), 4);
        INTREPID2_TEST_ERROR_EXPECTED( hexBasis.getValues(badVals7, hexNodes, OPERATOR_VALUE) );
    
        // Check if number of thrown exceptions matches the one we expect
        if (nthrow != ncatch) {
          errorFlag++;
          *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
        }
#endif
      }
      catch (std::logic_error err) {
        *outStream << "UNEXPECTED ERROR !!! ----------------------------------------------------------\n";
        *outStream << err.what() << '\n';
        *outStream << "-------------------------------------------------------------------------------" << "\n\n";
        errorFlag = -1000;
      };
  

  
      *outStream \
        << "\n"
        << "===============================================================================\n"\
        << "| TEST 2: correctness of tag to enum and enum to tag lookups                  |\n"\
        << "===============================================================================\n";

      try{
        Basis_HDIV_HEX_I1_FEM<DeviceSpaceType> hexBasis;

        const auto numFields = hexBasis.getCardinality();
        const auto allTags = hexBasis.getAllDofTags();

        // Loop over all tags, lookup the associated dof enumeration and then lookup the tag again
        const auto dofTagSize = allTags.dimension(0);
        for (size_type i=0;i<dofTagSize;++i) {
          const auto bfOrd = hexBasis.getDofOrdinal(allTags(i,0), allTags(i,1), allTags(i,2));

          const auto myTag = hexBasis.getDofTag(bfOrd);
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
        for(auto bfOrd=0;bfOrd<numFields;++bfOrd) {
          const auto myTag  = hexBasis.getDofTag(bfOrd);
          const auto myBfOrd = hexBasis.getDofOrdinal(myTag(0), myTag(1), myTag(2));
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
      } catch (std::logic_error err){
        *outStream << "UNEXPECTED ERROR !!! ----------------------------------------------------------\n";
        *outStream << err.what() << '\n';
        *outStream << "-------------------------------------------------------------------------------" << "\n\n";
        errorFlag = -1000;
      };

      *outStream \
        << "\n"
        << "===============================================================================\n"\
        << "| TEST 3: correctness of basis function values                                |\n"\
        << "===============================================================================\n";

      outStream -> precision(20);

      try{

        // VALUE: Each row pair gives the 6x3 correct basis set values at an evaluation point: (P,F,D) layout
        value_type basisValues[] = {
          // bottom 4 vertices
            0.,-0.25,0.,  0.,0.,0.,  0.,0.,0., -0.25,0.,0.,  0.,0.,-0.25,  0.,0.,0.,
            0.,-0.25,0.,  0.25,0.,0.,  0.,0.,0.,  0.,0.,0.,  0.,0.,-0.25,  0.,0.,0.,
            0.,0.,0.,  0.25,0.,0.,  0.,0.25,0.,  0.,0.,0.,  0.,0.,-0.25,  0.,0.,0.,
            0.,0.,0.,  0.,0.,0.,  0.,0.25,0., -0.25,0.,0.,  0.,0.,-0.25,  0.,0.,0.,
          // top 4 vertices
            0.,-0.25,0.,  0.,0.,0.,  0.,0.,0.,  -0.25,0.,0.,  0.,0.,0.,  0.,0.,0.25,
            0.,-0.25,0.,  0.25,0.,0.,  0.,0.,0.,  0.,0.,0.,  0.,0.,0.,  0.,0.,0.25,
            0.,0.,0.,  0.25,0.,0.,  0.,0.25,0.,  0.,0.,0.,  0.,0.,0.,  0.,0.,0.25,
            0.,0.,0.,  0.,0.,0.,  0.,0.25,0.,  -0.25,0.,0.,  0.,0.,0.,  0.,0.,0.25,
          // center {0, 0, 0}
            0.,-0.125,0.,  0.125,0.,0.,  0.,0.125,0.,  -0.125,0.,0.,  0.,0.,-0.125,  0.,0.,0.125,
          // faces { 1, 0, 0} and {-1, 0, 0}
            0.,-0.125,0.,  0.25,0.,0.,  0.,0.125,0.,  0.,0.,0.,  0.,0.,-0.125,  0.,0.,0.125,
            0.,-0.125,0.,  0.,0.,0.,  0.,0.125,0., -0.25,0.,0.,  0.,0.,-0.125,  0.,0.,0.125,
          // faces { 0, 1, 0} and { 0,-1, 0}
            0.,0.,0.,  0.125,0.,0.,  0.,0.25,0.,  -0.125,0.,0.,  0.,0.,-0.125,  0.,0.,0.125,
            0.,-0.25,0.,  0.125,0.,0.,  0.,0.,0.,  -0.125,0.,0.,  0.,0.,-0.125,  0.,0.,0.125,
          // faces {0, 0, 1} and {0, 0, -1}
            0.,-0.125,0.,  0.125,0.,0.,  0.,0.125,0.,  -0.125,0.,0.,  0.,0.,0.,  0.,0.,0.25,
            0.,-0.125,0.,  0.125,0.,0.,  0.,0.125,0.,  -0.125,0.,0.,  0.,0.,-0.25,  0.,0.,0.
        };

        // DIV: each row gives the 6 correct values of the divergence of the 6 basis functions: (P,F) layout
        value_type basisDivs[] = {
          // bottom 4 vertices
            0.125, 0.125, 0.125, 0.125, 0.125, 0.125,
            0.125, 0.125, 0.125, 0.125, 0.125, 0.125,
            0.125, 0.125, 0.125, 0.125, 0.125, 0.125,
            0.125, 0.125, 0.125, 0.125, 0.125, 0.125,
          // top 4 vertices
            0.125, 0.125, 0.125, 0.125, 0.125, 0.125,
            0.125, 0.125, 0.125, 0.125, 0.125, 0.125,
            0.125, 0.125, 0.125, 0.125, 0.125, 0.125,
            0.125, 0.125, 0.125, 0.125, 0.125, 0.125,
          // center {0, 0, 0}
            0.125, 0.125, 0.125, 0.125, 0.125, 0.125,
          // faces { 1, 0, 0} and {-1, 0, 0}
            0.125, 0.125, 0.125, 0.125, 0.125, 0.125,
            0.125, 0.125, 0.125, 0.125, 0.125, 0.125,
          // faces { 0, 1, 0} and { 0,-1, 0}
            0.125, 0.125, 0.125, 0.125, 0.125, 0.125,
            0.125, 0.125, 0.125, 0.125, 0.125, 0.125,
          // faces {0, 0, 1} and {0, 0, -1}
            0.125, 0.125, 0.125, 0.125, 0.125, 0.125,
            0.125, 0.125, 0.125, 0.125, 0.125, 0.125,
        };

        Basis_HDIV_HEX_I1_FEM<DeviceSpaceType> hexBasis;
        // Define array containing the 8 vertices of the reference HEX, its center and 6 face centers
        DynRankView ConstructWithLabel(hexNodes, 15, 3);
  
        hexNodes(0,0) = -1.0;  hexNodes(0,1) = -1.0;  hexNodes(0,2) = -1.0;
        hexNodes(1,0) =  1.0;  hexNodes(1,1) = -1.0;  hexNodes(1,2) = -1.0;
        hexNodes(2,0) =  1.0;  hexNodes(2,1) =  1.0;  hexNodes(2,2) = -1.0;
        hexNodes(3,0) = -1.0;  hexNodes(3,1) =  1.0;  hexNodes(3,2) = -1.0;
  
        hexNodes(4,0) = -1.0;  hexNodes(4,1) = -1.0;  hexNodes(4,2) =  1.0;
        hexNodes(5,0) =  1.0;  hexNodes(5,1) = -1.0;  hexNodes(5,2) =  1.0;
        hexNodes(6,0) =  1.0;  hexNodes(6,1) =  1.0;  hexNodes(6,2) =  1.0;
        hexNodes(7,0) = -1.0;  hexNodes(7,1) =  1.0;  hexNodes(7,2) =  1.0;
  
        hexNodes(8,0) =  0.0;  hexNodes(8,1) =  0.0;  hexNodes(8,2) =  0.0;
  
        hexNodes(9,0) =  1.0;  hexNodes(9,1) =  0.0;  hexNodes(9,2) =  0.0;
        hexNodes(10,0)= -1.0;  hexNodes(10,1)=  0.0;  hexNodes(10,2)=  0.0;
  
        hexNodes(11,0)=  0.0;  hexNodes(11,1)=  1.0;  hexNodes(11,2)=  0.0;
        hexNodes(12,0)=  0.0;  hexNodes(12,1)= -1.0;  hexNodes(12,2)=  0.0;
  
        hexNodes(13,0)=  0.0;  hexNodes(13,1)=  0.0;  hexNodes(13,2)=  1.0;
        hexNodes(14,0)=  0.0;  hexNodes(14,1)=  0.0;  hexNodes(14,2)= -1.0;

        
        // Dimensions for the output arrays:
        int numPoints = hexNodes.dimension(0);
        int numFields = hexBasis.getCardinality();
        int spaceDim  = hexBasis.getBaseCellTopology().getDimension();


        // Generic array for values and curls that will be properly sized before each call
        DynRankView ConstructWithLabel(vals, numFields, numPoints, spaceDim);

        // Check VALUE of basis functions: resize vals to rank-3 container:
        hexBasis.getValues(vals, hexNodes, OPERATOR_VALUE);
        for (int i = 0; i < numFields; i++) {
          for (int j = 0; j < numPoints; j++) {
            for (int k = 0; k < spaceDim; k++) {

              // compute offset for (P,F,D) data layout: indices are P->j, F->i, D->k
               int l = k + i * spaceDim + j * spaceDim * numFields;
               if (std::abs(vals(i,j,k) - basisValues[l]) > tol) {
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

        // Check DIV of basis function: resize vals to rank-2 container
        vals = DynRankView("vals", numFields, numPoints);
        hexBasis.getValues(vals, hexNodes, OPERATOR_DIV);
        for (int i = 0; i < numFields; i++) {
          for (int j = 0; j < numPoints; j++) {
              int l =  i + j * numFields;
               if (std::abs(vals(i,j) - basisDivs[l]) > tol) {
                 errorFlag++;
                 *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";

                 // Output the multi-index of the value where the error is:
                 *outStream << " At multi-index { ";
                 *outStream << i << " ";*outStream << j << " ";
                 *outStream << "}  computed divergence component: " << vals(i,j)
                   << " but reference divergence component: " << basisDivs[l] << "\n";
             }
          }
        }
      }

      // Catch unexpected errors
      catch (std::logic_error err) {
        *outStream << err.what() << "\n\n";
        errorFlag = -1000;
      };
    
      *outStream \
        << "\n"
        << "===============================================================================\n"\
        << "| TEST 4: correctness of DoF locations                                        |\n"\
        << "===============================================================================\n";
    
      try{
        Basis_HDIV_HEX_I1_FEM<DeviceSpaceType> hexBasis;
        const auto numFields = hexBasis.getCardinality();
        const auto spaceDim  = hexBasis.getBaseCellTopology().getDimension();
    
        // Check exceptions.
        ordinal_type nthrow = 0, ncatch = 0;
#ifdef HAVE_INTREPID2_DEBUG
        {
          DynRankView ConstructWithLabel(badVals,1,2,3);
          INTREPID2_TEST_ERROR_EXPECTED( hexBasis.getDofCoords(badVals) );
        }
        {
          DynRankView ConstructWithLabel(badVals, 3,2);
          INTREPID2_TEST_ERROR_EXPECTED( hexBasis.getDofCoords(badVals) );
        }
        {
          DynRankView ConstructWithLabel(badVals, 4,2);
          INTREPID2_TEST_ERROR_EXPECTED( hexBasis.getDofCoords(badVals) );
        }

         // Check if number of thrown exceptions matches the one we expect
         if (nthrow != ncatch) {
           errorFlag++;
           *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
         }
#endif
        // Check mathematical correctness
        DynRankView ConstructWithLabel(normals, hexBasis.getCardinality(),spaceDim); // normals at each point basis point
        normals(0,0)  =  0.0; normals(0,1)  = -4.0; normals(0,2)  =  0.0;
        normals(1,0)  =  4.0; normals(1,1)  =  0.0; normals(1,2)  =  0.0;
        normals(2,0)  =  0.0; normals(2,1)  =  4.0; normals(2,2)  =  0.0;
        normals(3,0)  = -4.0; normals(3,1)  =  0.0; normals(3,2)  =  0.0;
        normals(4,0)  =  0.0; normals(4,1)  =  0.0; normals(4,2)  = -4.0;
        normals(5,0)  =  0.0; normals(5,1)  =  0.0; normals(5,2)  =  4.0;

        DynRankView ConstructWithLabel(cvals, numFields,spaceDim);
        DynRankView ConstructWithLabel(bvals, numFields, numFields, spaceDim); // last dimension is spatial dim


        hexBasis.getDofCoords(cvals);
        hexBasis.getValues(bvals, cvals, OPERATOR_VALUE);

        value_type expected_normal;
        for (size_type i=0; i<bvals.dimension(0); i++) {
          for (size_type j=0; j<bvals.dimension(1); j++) {

            value_type normal = 0.0;
            for(size_type d=0;d<spaceDim;d++)
               normal += bvals(i,j,d)*normals(j,d);

            expected_normal = (i == j);
            if (std::abs(normal - expected_normal) > tol) {
              errorFlag++;
              std::stringstream ss;
              ss << "\nValue of basis function " << i << " at (" << cvals(i,0) << ", " << cvals(i,1)<< ", " << cvals(i,2) << ") is " << normal << " but should be " << expected_normal << "\n";
              *outStream << ss.str();
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
  }
}
