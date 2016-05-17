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

/** \file test_01.hpp
    \brief  Unit tests for the Intrepid2::C_HEX_I1_FEM class.
    \author Created by P. Bochev, D. Ridzal, K. Peterson and Kyungjoo Kim
*/

#include "Intrepid2_config.h"

#ifdef HAVE_INTREPID2_DEBUG
#define INTREPID2_TEST_FOR_DEBUG_ABORT_OVERRIDE_TO_CONTINUE
#endif

#include "Intrepid2_HCURL_HEX_I1_FEM.hpp"

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
    
    template <typename ValueType, typename DeviceSpaceType>
    int HCURL_HEX_I1_FEM_Test01(const bool verbose) {

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
        << "|                 Unit Test (Basis_HCURL_HEX_I1_FEM)                          |\n"
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

      *outStream
        << "\n"
        << "===============================================================================\n"
        << "| TEST 1: Basis creation, exceptions tests                                    |\n"
        << "===============================================================================\n";
  
      try {
        ordinal_type nthrow = 0, ncatch = 0;
#ifdef HAVE_INTREPID2_DEBUG
        Basis_HCURL_HEX_I1_FEM<DeviceSpaceType> hexBasis;

        // Define array containing the 8 vertices of the reference HEX, its center and 6 face centers
        DynRankView ConstructWithLabel(hexNodes, 15, 3);

        const auto numPoints = hexNodes.dimension(0);
        const auto numFields = hexBasis.getCardinality();
        const auto spaceDim  = hexBasis.getBaseCellTopology().getDimension();

        DynRankView ConstructWithLabel(vals, numFields, numPoints, spaceDim );
        
        {
          // exception #1: GRAD cannot be applied to HCURL functions 
          INTREPID2_TEST_ERROR_EXPECTED( hexBasis.getValues(vals, hexNodes, OPERATOR_GRAD) );
          
          // exception #2: DIV cannot be applied to HCURL functions
          INTREPID2_TEST_ERROR_EXPECTED( hexBasis.getValues(vals, hexNodes, OPERATOR_DIV) );
        }
        // Exceptions 3-7: all bf tags/bf Ids below are wrong and should cause getDofOrdinal() and 
        // getDofTag() to access invalid array elements thereby causing bounds check exception
        {
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
        }

        // Exceptions 8-15 test exception handling with incorrectly dimensioned input/output arrays
        // exception #8: input points array must be of rank-2
        {
          DynRankView ConstructWithLabel( badPoints1, 4, 5, 3);
          INTREPID2_TEST_ERROR_EXPECTED( hexBasis.getValues(vals, badPoints1, OPERATOR_VALUE) );
        }
        {
          // exception #9 dimension 1 in the input point array must equal space dimension of the cell
          DynRankView ConstructWithLabel( badPoints2, 4, 2);
          INTREPID2_TEST_ERROR_EXPECTED( hexBasis.getValues(vals, badPoints2, OPERATOR_VALUE) );
        }
        {
          // exception #10 output values must be of rank-3 for OPERATOR_VALUE
          DynRankView ConstructWithLabel( badVals1, 4, 3);
          INTREPID2_TEST_ERROR_EXPECTED( hexBasis.getValues(badVals1, hexNodes, OPERATOR_VALUE) );
        }
        {
          // exception #11 output values must be of rank-3 for OPERATOR_CURL
          DynRankView ConstructWithLabel( badVals1, 4, 3);
          INTREPID2_TEST_ERROR_EXPECTED( hexBasis.getValues(badVals1, hexNodes, OPERATOR_CURL) );
        }
        {
          // exception #12 incorrect 0th dimension of output array (must equal number of basis functions)
          DynRankView ConstructWithLabel( badVals2, numFields + 1, numPoints, spaceDim);
          INTREPID2_TEST_ERROR_EXPECTED( hexBasis.getValues(badVals2, hexNodes, OPERATOR_VALUE) ) ;
        }
        {
          // exception #13 incorrect 1st  dimension of output array (must equal number of points)
          DynRankView ConstructWithLabel( badVals3, numFields, numPoints + 1, spaceDim);
          INTREPID2_TEST_ERROR_EXPECTED( hexBasis.getValues(badVals3, hexNodes, OPERATOR_VALUE) ) ;
        }
        {
          // exception #14: incorrect 2nd dimension of output array (must equal the space dimension)
          DynRankView ConstructWithLabel( badVals4, numFields, numPoints, spaceDim + 1);
          INTREPID2_TEST_ERROR_EXPECTED( hexBasis.getValues(badVals4, hexNodes, OPERATOR_VALUE) ) ;
        }
        {
          // exception #15: incorrect 2nd dimension of output array (must equal the space dimension)
          DynRankView ConstructWithLabel( badVals4, numFields, numPoints, spaceDim + 1);
          INTREPID2_TEST_ERROR_EXPECTED( hexBasis.getValues(badVals4, hexNodes, OPERATOR_CURL) ) ;
        }
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
  
      try {
        Basis_HCURL_HEX_I1_FEM<DeviceSpaceType> hexBasis;

        const auto numFields = hexBasis.getCardinality();
        const auto allTags = hexBasis.getAllDofTags();
        
        // Loop over all tags, lookup the associated dof enumeration and then lookup the tag again
        const auto dofTagSize = allTags.dimension(0);
        for (auto i=0;i<dofTagSize;++i) {
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

      try {

        // VALUE: Each row pair gives the 12x3 correct basis set values at an evaluation point: (P,F,D) layout
        ValueType basisValues[] = {
          // bottom 4 vertices
          0.5,0.,0.,  0.,0.,0.,  0.,0.,0.,  0.,-0.5,0.,  0.,0.,0.,  0.,0.,0.,  
          0.,0.,0.,  0.,0.,0.,  0.,0.,0.5,  0.,0.,0.,  0.,0.,0.,  0.,0.,0.,
          
          0.5,0.,0.,  0.,0.5,0.,  0.,0.,0.,  0.,0.,0.,  0.,0.,0.,  0.,0.,0.,
          0.,0.,0.,  0.,0.,0.,  0.,0.,0.,  0.,0.,0.5,  0.,0.,0.,  0.,0.,0.,
          
          0.,0.,0.,  0.,0.5,0.,  -0.5,0.,0.,  0.,0.,0.,  0.,0.,0.,  0.,0.,0.,
          0.,0.,0.,  0.,0.,0.,  0.,0.,0.,  0.,0.,0.,  0.,0.,0.5,  0.,0.,0.,
          
          0.,0.,0.,  0.,0.,0.,  -0.5,0.,0.,  0.,-0.5,0.,  0.,0.,0.,  0.,0.,0.,
          0.,0.,0.,  0.,0.,0.,  0.,0.,0.,  0.,0.,0.,  0.,0.,0.,  0.,0.,0.5,
          
          // top 4 vertices
          0.,0.,0.,  0.,0.,0.,  0.,0.,0.,  0.,0.,0.,  0.5,0.,0.,  0.,0.,0.,
          0.,0.,0.,  0.,-0.5,0.,  0.,0.,0.5,  0.,0.,0.,  0.,0.,0.,  0.,0.,0.,
          
          0.,0.,0.,  0.,0.,0.,  0.,0.,0.,  0.,0.,0.,  0.5,0.,0.,  0.,0.5,0.,
          0.,0.,0.,  0.,0.,0.,  0.,0.,0.,  0.,0.,0.5,  0.,0.,0.,  0.,0.,0.,
          
          0.,0.,0.,  0.,0.,0.,  0.,0.,0.,  0.,0.,0.,  0.,0.,0.,  0.,0.5,0.,
          -0.5,0.,0.,  0.,0.,0.,  0.,0.,0.,  0.,0.,0.,  0.,0.,0.5,  0.,0.,0.,
          
          0.,0.,0.,  0.,0.,0.,  0.,0.,0.,  0.,0.,0.,  0.,0.,0.,  0.,0.,0.,
          -0.5,0.,0.,  0.,-0.5,0.,  0.,0.,0.,  0.,0.,0.,  0.,0.,0.,  0.,0.,0.5,
          
          // center {0, 0, 0}
          0.125,0.,0.,  0.,0.125,0.,  -0.125,0.,0.,  0.,-0.125,0.,  0.125,0.,0.,  0.,0.125,0.,
          -0.125,0.,0.,  0.,-0.125,0.,  0.,0.,0.125,  0.,0.,0.125,  0.,0.,0.125,  0.,0.,0.125,
          
          // faces { 1, 0, 0} and {-1, 0, 0}
          0.125,0.,0.,  0.,0.25,0.,  -0.125,0.,0.,  0.,0.,0.,  0.125,0.,0.,  0.,0.25,0.,
          -0.125,0.,0.,  0.,0.,0.,  0.,0.,0.,  0.,0.,0.25,  0.,0.,0.25,  0.,0.,0.,  
          
          0.125,0.,0.,  0.,0.,0.,  -0.125,0.,0.,  0.,-0.25,0.,  0.125,0.,0.,  0.,0.,0.,
          -0.125,0.,0.,  0.,-0.25,0.,  0.,0.,0.25,  0.,0.,0.,  0.,0.,0.,  0.,0.,0.25,
          
          // faces { 0, 1, 0} and { 0,-1, 0}
          0.,0.,0.,  0.,0.125,0.,  -0.25,0.,0.,  0.,-0.125,0.,  0.,0.,0.,  0.,0.125,0.,
          -0.25,0.,0.,  0.,-0.125,0.,  0.,0.,0.,  0.,0.,0.,  0.,0.,0.25,  0.,0.,0.25,
          
          0.25,0.,0.,  0.,0.125,0.,  0.,0.,0.,  0.,-0.125,0.,  0.25,0.,0.,  0.,0.125,0.,
          0.,0.,0.,  0.,-0.125,0.,  0.,0.,0.25,  0.,0.,0.25,  0.,0.,0.,  0.,0.,0.,
          
          // faces {0, 0, 1} and {0, 0, -1}
          0.,0.,0.,  0.,0.,0.,  0.,0.,0.,  0.,0.,0.,  0.25,0.,0.,  0.,0.25,0.,
          -0.25,0.,0.,  0.,-0.25,0.,  0.,0.,0.125,  0.,0.,0.125,  0.,0.,0.125,  0.,0.,0.125,
          
          0.25,0.,0.,  0.,0.25,0.,  -0.25,0.,0.,  0.,-0.25,0.,  0.,0.,0.,  0.,0.,0.,  
          0.0,0.,0.,  0.,0.,0.0,  0.,0.,0.125,  0.,0.,0.125,  0.,0.,0.125,  0.,0.,0.125
        };
        
        // CURL: each row pair gives the 3x12 correct values of the curls of the 12 basis functions: (P,F,D) layout
        ValueType basisCurls[] = {   
          // bottom 4 vertices
          0.,-0.25,0.25,  0.,0.,0.25,  0.,0.,0.25,  -0.25,0.,0.25,  0.,0.25,0.,  0.,0.,0.,  
          0.,0.,0.,  0.25,0.,0.,  -0.25,0.25,0.,  0.,-0.25,0.,  0.,0.,0.,  0.25,0.,0.,
          
          0.,-0.25,0.25,  0.25,0.,0.25,  0.,0.,0.25,  0.,0.,0.25,  0.,0.25,0.,  -0.25,0.,0., 
          0.,0.,0.,  0.,0.,0.,  0.,0.25,0.,  -0.25,-0.25,0.,  0.25,0.,0.,  0.,0.,0.,  
          
          0.,0.,0.25,  0.25,0.,0.25,  0.,0.25,0.25,  0.,0.,0.25,  0.,0.,0.,  -0.25,0.,0.,
          0.,-0.25,0.,  0.,0.,0.,  0.,0.,0.,  -0.25,0.,0.,  0.25,-0.25,0.,  0.,0.25,0.,
          
          0.,0.,0.25,  0.,0.,0.25,  0.,0.25,0.25,  -0.25,0.,0.25,  0.,0.,0.,  0.,0.,0.,
          0.,-0.25,0.,  0.25,0.,0.,  -0.25,0.,0.,  0.,0.,0.,  0.,-0.25,0.,  0.25,0.25,0.,
          
          // top 4 vertices
          0.,-0.25,0.,  0.,0.,0.,  0.,0.,0.,  -0.25,0.,0.,  0.,0.25,0.25,  0.,0.,0.25,
          0.,0.,0.25,  0.25,0.,0.25,  -0.25,0.25,0.,  0.,-0.25,0.,  0.,0.,0.,  0.25,0.,0.,
          
          0.,-0.25,0.,  0.25,0.,0.,  0.,0.,0.,  0.,0.,0.,  0.,0.25,0.25,  -0.25,0.,0.25,
          0.,0.,0.25,  0.,0.,0.25,  0.,0.25,0.,  -0.25,-0.25,0.,  0.25,0.,0.,  0.,0.,0.,
          
          0.,0.,0.,  0.25,0.,0.,  0.,0.25,0.,  0.,0.,0.,  0.,0.,0.25,  -0.25,0.,0.25,
          0.,-0.25,0.25,  0.,0.,0.25,  0.,0.,0.,  -0.25,0.,0.,  0.25,-0.25,0.,  0.,0.25,0.,
          
          0.,0.,0.,  0.,0.,0.,  0.,0.25,0.,  -0.25,0.,0.,  0.,0.,0.25,  0.,0.,0.25,
          0.,-0.25,0.25,  0.25,0.,0.25,  -0.25,0.,0., 0.,0.,0.,  0.,-0.25,0.,  0.25,0.25,0.,
          
          // center {0, 0, 0}
          0.,-0.125,0.125,  0.125,0.,0.125,  0.,0.125,0.125,  -0.125,0.,0.125,  0.,0.125,0.125,  -0.125,0.,0.125,
          0.,-0.125,0.125,  0.125,0.,0.125,  -0.125,0.125,0.,  -0.125,-0.125,0.,  0.125,-0.125,0.,  0.125,0.125,0.,
          
          // faces { 1, 0, 0} and {-1, 0, 0}
          0.,-0.125,0.125,  0.25,0.,0.125,  0.,0.125,0.125,  0.,0.,0.125,  0.,0.125,0.125,  -0.25,0.,0.125, 
          0.,-0.125,0.125,  0.,0.,0.125,  0.,0.125,0.,  -0.25,-0.125,0.,  0.25,-0.125,0.,  0.,0.125,0., 
          
          0.,-0.125,0.125,  0.,0.,0.125,  0.,0.125,0.125,  -0.25,0.,0.125,  0.,0.125,0.125,  0.,0.,0.125,
          0.,-0.125,0.125,  0.25,0.,0.125,  -0.25,0.125,0.,  0.,-0.125,0.,  0.,-0.125,0.,  0.25,0.125,0.,
          
          // faces { 0, 1, 0} and { 0,-1, 0}
          0.,0.,0.125,  0.125,0.,0.125,  0.,0.25,0.125,  -0.125,0.,0.125,  0.,0.,0.125,  -0.125,0.,0.125,
          0.,-0.25,0.125,  0.125,0.,0.125,  -0.125,0.,0.,  -0.125,0.,0.,  0.125,-0.25,0.,  0.125,0.25,0.,
          
          0.,-0.25,0.125,  0.125,0.,0.125,  0.,0.,0.125,  -0.125,0.,0.125,  0.,0.25,0.125,  -0.125,0.,0.125,
          0.,0.,0.125,  0.125,0.,0.125,  -0.125,0.25,0.,  -0.125,-0.25,0.,  0.125,0.,0.,  0.125,0.,0.,
          
          // faces {0, 0, 1} and {0, 0, -1}
          0.,-0.125,0.,  0.125,0.,0.,  0.,0.125,0.,  -0.125,0.,0.,  0.,0.125,0.25,  -0.125,0.,0.25,
          0.,-0.125,0.25,  0.125,0.,0.25,  -0.125,0.125,0.,  -0.125,-0.125,0.,  0.125,-0.125,0.,  0.125,0.125,0.,
          
          0.,-0.125,0.25,  0.125,0.,0.25,  0.,0.125,0.25,  -0.125,0.,0.25,  0.,0.125,0.,  -0.125,0.,0.,
          0.,-0.125,0.,  0.125,0.,0.,  -0.125,0.125,0.,  -0.125,-0.125,0.,  0.125,-0.125,0.,  0.125,0.125,0.
        };

        Basis_HCURL_HEX_I1_FEM<DeviceSpaceType> hexBasis;

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
        const auto numPoints = hexNodes.dimension(0);
        const auto numFields = hexBasis.getCardinality();
        const auto spaceDim  = hexBasis.getBaseCellTopology().getDimension();
        
        // Generic array for values and curls that will be properly sized before each call
        DynRankView ConstructWithLabel(vals, numFields, numPoints, spaceDim);
    
        // Check VALUE of basis functions: resize vals to rank-3 container:
        hexBasis.getValues(vals, hexNodes, OPERATOR_VALUE);
        for (auto i=0;i<numFields;++i) 
          for (auto j=0;j<numPoints;++j) 
            for (auto k=0;k<spaceDim;++k) {
          
              // compute offset for (P,F,D) data layout: indices are P->j, F->i, D->k
              const auto l = k + i * spaceDim + j * spaceDim * numFields;
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
        
        // Check CURL of basis function: resize vals to rank-3 container
        hexBasis.getValues(vals, hexNodes, OPERATOR_CURL);

        for (auto i=0;i<numFields;++i)
          for (auto j=0;j<numPoints;++j)
            for (auto k=0;k<spaceDim;++k) {
              
              // compute offset for (P,F,D) data layout: indices are P->j, F->i, D->k
              const auto l = k + i * spaceDim + j * spaceDim * numFields;
              if (std::abs(vals(i,j,k) - basisCurls[l]) > tol) {
                errorFlag++;
                *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
                
                // Output the multi-index of the value where the error is: 
                *outStream << " At multi-index { ";
                *outStream << i << " ";*outStream << j << " ";*outStream << k << " ";
                *outStream << "}  computed curl component: " << vals(i,j,k) 
                           << " but reference curl component: " << basisCurls[l] << "\n";
              }
            }
        
        // Catch unexpected errors
      } catch (std::logic_error err) {
        *outStream << err.what() << "\n\n";
        errorFlag = -1000;
      }
      
      *outStream                               
        << "\n"
        << "===============================================================================\n" 
        << "| TEST 4: correctness of DoF locations                                        |\n"
        << "===============================================================================\n";
      
      try{
        Basis_HCURL_HEX_I1_FEM<DeviceSpaceType> hexBasis;
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
#endif
        // Check if number of thrown exceptions matches the one we expect
        if (nthrow != ncatch) {
          errorFlag++;
          *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
        }

        // Check mathematical correctness
        DynRankView ConstructWithLabel(tangents, numFields,spaceDim); // tangents at each point basis$
        tangents(0,0)  =  2.0; tangents(0,1)  =  0.0; tangents(0,2)  = 0.0;
        tangents(1,0)  =  0.0; tangents(1,1)  =  2.0; tangents(1,2)  = 0.0;
        tangents(2,0)  = -2.0; tangents(2,1)  =  0.0; tangents(2,2)  = 0.0;
        tangents(3,0)  =  0.0; tangents(3,1)  = -2.0; tangents(3,2)  = 0.0;
        tangents(4,0)  =  2.0; tangents(4,1)  =  0.0; tangents(4,2)  = 0.0;
        tangents(5,0)  =  0.0; tangents(5,1)  =  2.0; tangents(5,2)  = 0.0;
        tangents(6,0)  = -2.0; tangents(6,1)  =  0.0; tangents(6,2)  = 0.0;
        tangents(7,0)  =  0.0; tangents(7,1)  = -2.0; tangents(7,2)  = 0.0;
        tangents(8,0)  =  0.0; tangents(8,1)  =  0.0; tangents(8,2)  = 2.0;
        tangents(9,0)  =  0.0; tangents(9,1)  =  0.0; tangents(9,2)  = 2.0;
        tangents(10,0) =  0.0; tangents(10,1) =  0.0; tangents(10,2) = 2.0;
        tangents(11,0) =  0.0; tangents(11,1) =  0.0; tangents(11,2) = 2.0;
        
        DynRankView ConstructWithLabel(cvals, numFields, spaceDim);
        DynRankView ConstructWithLabel(bvals, numFields, numFields, spaceDim); // last dimension is spatial dim

        hexBasis.getDofCoords(cvals);
        hexBasis.getValues(bvals, cvals, OPERATOR_VALUE);

        for (size_type i=0;i<numFields;++i) 
          for (size_type j=0;j<numFields;++j) {
            
            ValueType tangent = 0.0;
            for(size_type d=0;d<spaceDim;++d)
              tangent += bvals(i,j,d)*tangents(j,d);

            const ValueType expected_tangent = (i == j);
            if (std::abs(tangent - expected_tangent) > tol) {
              errorFlag++;
              std::stringstream ss;
              ss << "\nValue of basis function " << i << " at (" << cvals(i,0) << ", " << cvals(i,1)<< ", " << cvals(i,2) << ") is " << tangent << " but should be " << expected_tangent << "\n";
              *outStream << ss.str();
            }
          }
      } catch (std::logic_error err) {
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
