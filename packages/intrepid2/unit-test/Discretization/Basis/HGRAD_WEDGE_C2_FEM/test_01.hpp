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

/**   \file   test_01.cpp
      \brief  Unit tests for the Intrepid2::HGRAD_WEDGE_C2_FEM class.
      \author Created by P. Bochev, D. Ridzal, K. Peterson and Kyungjoo Kim.
*/

#include "Intrepid2_config.h"

#ifdef HAVE_INTREPID2_DEBUG
#define INTREPID2_TEST_FOR_DEBUG_ABORT_OVERRIDE_TO_CONTINUE
#endif

#include "Intrepid2_HGRAD_WEDGE_C2_FEM.hpp"

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
    int HGRAD_WEDGE_C2_FEM_Test01(const bool verbose) {
      
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
        << "|                 Unit Test (Basis_HGRAD_WEDGE_C2_FEM)                        |\n"
        << "|                                                                             |\n"
        << "|     1) Conversion of Dof tags into Dof ordinals and back                    |\n"
        << "|     2) Basis values for VALUE, GRAD, and Dk operators                       |\n"
        << "|                                                                             |\n"
        << "|  Questions? Contact  Pavel Bochev  (pbboche@sandia.gov),                    |\n"
        << "|                      Denis Ridzal  (dridzal@sandia.gov),                    |\n"
        << "|                      Kara Peterson (kjpeter@sandia.gov),                    |\n"
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
      Basis_HGRAD_WEDGE_C2_FEM<DeviceSpaceType,outputValueType,pointValueType> wedgeBasis;

      *outStream
        << "\n"
        << "===============================================================================\n"
        << "| TEST 1: constructors and exceptions                                         |\n"
        << "===============================================================================\n";


      try {
        ordinal_type nthrow = 0, ncatch = 0;
#ifdef HAVE_INTREPID2_DEBUG

        // Define array containing the 4 vertices of the reference WEDGE and its center.
        DynRankView ConstructWithLabel(wedgeNodes, 18, 3);

        // Generic array for the output values; needs to be properly resized depending on the operator type
        const ordinal_type numFields = wedgeBasis.getCardinality();
        const ordinal_type numPoints = wedgeNodes.extent(0);
        const ordinal_type spaceDim  = wedgeBasis.getBaseCellTopology().getDimension();

        DynRankView vals("vals", numFields, numPoints);
        DynRankView vals_vec("vals", numFields, numPoints, spaceDim);

        {
          // exception #1: CURL cannot be applied to scalar functions
          INTREPID2_TEST_ERROR_EXPECTED( wedgeBasis.getValues(vals_vec, wedgeNodes, OPERATOR_DIV) );

          // exception #2: DIV cannot be applied to scalar functions
          INTREPID2_TEST_ERROR_EXPECTED( wedgeBasis.getValues(vals_vec, wedgeNodes, OPERATOR_DIV) );
        }

        // Exceptions 3-7: all bf tags/bf Ids below are wrong and should cause getDofOrdinal() and 
        // getDofTag() to access invalid array elements thereby causing bounds check exception
        {
          // exception #3
          INTREPID2_TEST_ERROR_EXPECTED( wedgeBasis.getDofOrdinal(3,0,0) );
          // exception #4
          INTREPID2_TEST_ERROR_EXPECTED( wedgeBasis.getDofOrdinal(1,1,1) );
          // exception #5
          INTREPID2_TEST_ERROR_EXPECTED( wedgeBasis.getDofOrdinal(0,9,0) );
          // exception #6
          INTREPID2_TEST_ERROR_EXPECTED( wedgeBasis.getDofTag(numFields) );
          // exception #7
          INTREPID2_TEST_ERROR_EXPECTED( wedgeBasis.getDofTag(-1) );
        }
        // Exceptions 8-18 test exception handling with incorrectly dimensioned input/output arrays
        {
          // exception #8: input points array must be of rank-2
          DynRankView ConstructWithLabel( badPoints1, 4, 5, 3);
          INTREPID2_TEST_ERROR_EXPECTED( wedgeBasis.getValues(vals, badPoints1, OPERATOR_VALUE) );
        }
        {
          // exception #9 dimension 1 in the input point array must equal space dimension of the cell
          DynRankView ConstructWithLabel( badPoints2, 4, spaceDim + 1);
          INTREPID2_TEST_ERROR_EXPECTED( wedgeBasis.getValues(vals, badPoints2, OPERATOR_VALUE) );
        }
        {
          // exception #10 output values must be of rank-2 for OPERATOR_VALUE
          DynRankView ConstructWithLabel( badVals1, 4, 3, 1);
          INTREPID2_TEST_ERROR_EXPECTED( wedgeBasis.getValues(badVals1, wedgeNodes, OPERATOR_VALUE) );
        }
        {
          // exception #11 output values must be of rank-3 for OPERATOR_GRAD
          DynRankView ConstructWithLabel( badVals2, 4, 3);
          INTREPID2_TEST_ERROR_EXPECTED( wedgeBasis.getValues(badVals2, wedgeNodes, OPERATOR_GRAD) );

          // exception #12 output values must be of rank-3 for OPERATOR_D1
          INTREPID2_TEST_ERROR_EXPECTED( wedgeBasis.getValues(badVals2, wedgeNodes, OPERATOR_D1) );
    
          // exception #13 output values must be of rank-3 for OPERATOR_D2
          INTREPID2_TEST_ERROR_EXPECTED( wedgeBasis.getValues(badVals2, wedgeNodes, OPERATOR_D2) );
        }
        {
          // exception #14 incorrect 0th dimension of output array (must equal number of basis functions)
          DynRankView ConstructWithLabel( badVals3, numFields + 1, numPoints);
          INTREPID2_TEST_ERROR_EXPECTED( wedgeBasis.getValues(badVals3, wedgeNodes, OPERATOR_VALUE) );
        }
        {
          // exception #15 incorrect 1st dimension of output array (must equal number of points)
          DynRankView ConstructWithLabel( badVals4, numFields, numPoints + 1);
          INTREPID2_TEST_ERROR_EXPECTED( wedgeBasis.getValues(badVals4, wedgeNodes, OPERATOR_VALUE) );
        }
        {
          // exception #16: incorrect 2nd dimension of output array (must equal the space dimension)
          DynRankView ConstructWithLabel( badVals5, numFields, numPoints, spaceDim - 1);
          INTREPID2_TEST_ERROR_EXPECTED( wedgeBasis.getValues(badVals5, wedgeNodes, OPERATOR_GRAD) );
        }
        {
          // exception #17: incorrect 2nd dimension of output array (must equal D2 cardinality in 3D)
          DynRankView ConstructWithLabel( badVals6, numFields, numPoints, 40);
          INTREPID2_TEST_ERROR_EXPECTED( wedgeBasis.getValues(badVals6, wedgeNodes, OPERATOR_D2) );
          // exception #18: incorrect 2nd dimension of output array (must equal D3 cardinality in 3D)
          INTREPID2_TEST_ERROR_EXPECTED( wedgeBasis.getValues(badVals6, wedgeNodes, OPERATOR_D3) );
        }
#endif
        if (nthrow != ncatch) {
          errorFlag++;
          *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
          *outStream << "# of catch ("<< ncatch << ") is different from # of throw (" << ncatch << ")\n";
        }
      } catch (std::exception err) {
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
        const ordinal_type numFields = wedgeBasis.getCardinality();
        const auto allTags = wedgeBasis.getAllDofTags();

        // Loop over all tags, lookup the associated dof enumeration and then lookup the tag again
        const ordinal_type dofTagSize = allTags.extent(0);
        for (ordinal_type i=0;i<dofTagSize;++i) {
          const auto bfOrd = wedgeBasis.getDofOrdinal(allTags(i,0), allTags(i,1), allTags(i,2));
          
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
        for (ordinal_type bfOrd=0;bfOrd<numFields;++bfOrd) {
          const auto myTag = wedgeBasis.getDofTag(bfOrd);
          
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
  
      // VALUE: correct basis function values in (F,P) format
      const ValueType basisValues[] = {
        1.00, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.00, 0, \
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.00, 0, 0, 0, 0, \
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.00, 0, 0, 0, 0, 0, 0, 0, \
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.00, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, \
        0, 0, 0, 0, 0, 0, 0, 0, 1.00, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, \
        0, 0, 0, 0, 0, 1.00, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, \
        0, 0, 1.00, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, \
        1.00, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.00, 0, \
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.00, 0, 0, 0, 0, \
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.00, 0, 0, 0, 0, 0, 0, 0, \
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.00, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, \
        0, 0, 0, 0, 0, 0, 0, 0, 1.00, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, \
        0, 0, 0, 0, 0, 1.00, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, \
        0, 0, 1.00, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, \
        1.00, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.00  };

      // GRAD and D1 values are stored in (F,P,D) format in a data file. Read file and do the test     
      std::vector<ValueType> basisGrads;           // Flat array for the gradient values.
      { 
        std::ifstream dataFile("../testdata/WEDGE_C2_GradVals.dat");
        
        INTREPID2_TEST_FOR_EXCEPTION( !dataFile.good(), std::logic_error,
                                      ">>> ERROR (HGRAD_WEDGE_C2/test01): could not open GRAD values data file, test aborted.");
        while (!dataFile.eof() ){
          ValueType temp;
          std::string line;                            // string for one line of input file
          std::getline(dataFile, line);           // get next line from file
          std::stringstream data_line(line);           // convert to stringstream
          while (data_line >> temp)               // extract value from line
            basisGrads.push_back(temp);           // push into vector
        }
      }
  
      //D2: flat array with the values of D2 applied to basis functions. Multi-index is (F,P,D2cardinality)
      std::vector<ValueType> basisD2;
      { 
        std::ifstream dataFile("../testdata/WEDGE_C2_D2Vals.dat");
        
        INTREPID2_TEST_FOR_EXCEPTION( !dataFile.good(), std::logic_error,
                                      ">>> ERROR (HGRAD_WEDGE_C2/test01): could not open D2 values data file, test aborted.");
        while (!dataFile.eof() ){
          ValueType temp;
          std::string line;                            // string for one line of input file
          std::getline(dataFile, line);           // get next line from file
          std::stringstream data_line(line);           // convert to stringstream
          while (data_line >> temp)               // extract value from line
            basisD2.push_back(temp);           // push into vector
        }
      }

      //D3: flat array with the values of D3 applied to basis functions. Multi-index is (F,P,D3cardinality)
      std::vector<ValueType> basisD3;
      { 
        std::ifstream dataFile("../testdata/WEDGE_C2_D3Vals.dat");
        
        INTREPID2_TEST_FOR_EXCEPTION( !dataFile.good(), std::logic_error,
                                      ">>> ERROR (HGRAD_WEDGE_C2/test01): could not open D3 values data file, test aborted.");
        while (!dataFile.eof() ){
          ValueType temp;
          std::string line;                            // string for one line of input file
          std::getline(dataFile, line);           // get next line from file
          std::stringstream data_line(line);           // convert to stringstream
          while (data_line >> temp)               // extract value from line
            basisD3.push_back(temp);           // push into vector
        }
      }
  
      //D4: flat array with the values of D3 applied to basis functions. Multi-index is (F,P,D4cardinality)
      std::vector<ValueType> basisD4;
      { 
        std::ifstream dataFile("../testdata/WEDGE_C2_D4Vals.dat");
        
        INTREPID2_TEST_FOR_EXCEPTION( !dataFile.good(), std::logic_error,
                                      ">>> ERROR (HGRAD_WEDGE_C2/test01): could not open D4 values data file, test aborted.");
        while (!dataFile.eof() ){
          ValueType temp;
          std::string line;                            // string for one line of input file
          std::getline(dataFile, line);           // get next line from file
          std::stringstream data_line(line);           // convert to stringstream
          while (data_line >> temp)               // extract value from line
            basisD4.push_back(temp);           // push into vector
        }
      }
      
      try {
        DynRankViewHost ConstructWithLabel(wedgeNodesHost, 18, 3);

        wedgeNodesHost(0,0) =  0.0;  wedgeNodesHost(0,1) =  0.0;  wedgeNodesHost(0,2) = -1.0;  
        wedgeNodesHost(1,0) =  1.0;  wedgeNodesHost(1,1) =  0.0;  wedgeNodesHost(1,2) = -1.0;  
        wedgeNodesHost(2,0) =  0.0;  wedgeNodesHost(2,1) =  1.0;  wedgeNodesHost(2,2) = -1.0;
        wedgeNodesHost(3,0) =  0.0;  wedgeNodesHost(3,1) =  0.0;  wedgeNodesHost(3,2) =  1.0;  
        wedgeNodesHost(4,0) =  1.0;  wedgeNodesHost(4,1) =  0.0;  wedgeNodesHost(4,2) =  1.0;  
        wedgeNodesHost(5,0) =  0.0;  wedgeNodesHost(5,1) =  1.0;  wedgeNodesHost(5,2) =  1.0;
        
        wedgeNodesHost(6,0) =  0.5;  wedgeNodesHost(6,1) =  0.0;  wedgeNodesHost(6,2) = -1.0;  
        wedgeNodesHost(7,0) =  0.5;  wedgeNodesHost(7,1) =  0.5;  wedgeNodesHost(7,2) = -1.0;  
        wedgeNodesHost(8,0) =  0.0;  wedgeNodesHost(8,1) =  0.5;  wedgeNodesHost(8,2) = -1.0;
        wedgeNodesHost(9,0) =  0.0;  wedgeNodesHost(9,1) =  0.0;  wedgeNodesHost(9,2) =  0.0;
        wedgeNodesHost(10,0)=  1.0;  wedgeNodesHost(10,1)=  0.0;  wedgeNodesHost(10,2)=  0.0;  
        wedgeNodesHost(11,0)=  0.0;  wedgeNodesHost(11,1)=  1.0;  wedgeNodesHost(11,2)=  0.0;  
        
        wedgeNodesHost(12,0)=  0.5;  wedgeNodesHost(12,1)=  0.0;  wedgeNodesHost(12,2)=  1.0;  
        wedgeNodesHost(13,0)=  0.5;  wedgeNodesHost(13,1)=  0.5;  wedgeNodesHost(13,2)=  1.0;  
        wedgeNodesHost(14,0)=  0.0;  wedgeNodesHost(14,1)=  0.5;  wedgeNodesHost(14,2)=  1.0;  
        wedgeNodesHost(15,0)=  0.5;  wedgeNodesHost(15,1)=  0.0;  wedgeNodesHost(15,2)=  0.0;  
        wedgeNodesHost(16,0)=  0.5;  wedgeNodesHost(16,1)=  0.5;  wedgeNodesHost(16,2)=  0.0;  
        wedgeNodesHost(17,0)=  0.0;  wedgeNodesHost(17,1)=  0.5;  wedgeNodesHost(17,2)=  0.0;  
        
        auto wedgeNodes = Kokkos::create_mirror_view(typename DeviceSpaceType::memory_space(), wedgeNodesHost);
        Kokkos::deep_copy(wedgeNodes, wedgeNodesHost);

        // Dimensions for the output arrays:
        const ordinal_type numFields = wedgeBasis.getCardinality();
        const ordinal_type numPoints = wedgeNodes.extent(0);
        const ordinal_type spaceDim  = wedgeBasis.getBaseCellTopology().getDimension();
        const ordinal_type D2Cardin  = getDkCardinality(OPERATOR_D2, spaceDim);
        const ordinal_type D3Cardin  = getDkCardinality(OPERATOR_D3, spaceDim);
        const ordinal_type D4Cardin  = getDkCardinality(OPERATOR_D4, spaceDim);

        // Check VALUE of basis functions: resize vals to rank-2 container:
        {
          DynRankView vals = DynRankView("vals", numFields, numPoints);
          wedgeBasis.getValues(vals, wedgeNodes, OPERATOR_VALUE);
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
    
    
        // Check GRAD of basis function: resize vals to rank-3 container
        {
          DynRankView vals = DynRankView("vals", numFields, numPoints, spaceDim);
          wedgeBasis.getValues(vals, wedgeNodes, OPERATOR_GRAD);
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
        }

        // Check D1 of basis function (do not resize vals because it has the correct size: D1 = GRAD)
        {
          DynRankView vals = DynRankView("vals", numFields, numPoints, spaceDim);
          wedgeBasis.getValues(vals, wedgeNodes, OPERATOR_D1);
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
          wedgeBasis.getValues(vals, wedgeNodes, OPERATOR_D2);
          auto vals_host = Kokkos::create_mirror_view(typename HostSpaceType::memory_space(), vals);
          Kokkos::deep_copy(vals_host, vals);
          for (ordinal_type i = 0; i < numFields; ++i) {
            for (ordinal_type j = 0; j < numPoints; ++j) {
              for (ordinal_type k = 0; k < D2Cardin; ++k) {

                // basisGrads is (F,P,D), compute offset:
                const ordinal_type l = k + j * D2Cardin + i * D2Cardin * numPoints;
                if (std::abs(vals_host(i,j,k) - basisD2[l]) > tol) {
                  errorFlag++;
                  *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";

                  // Output the multi-index of the value where the error is:
                  *outStream << " At multi-index { ";
                  *outStream << i << " ";*outStream << j << " ";*outStream << k << " ";
                  *outStream << "}  computed D2 component: " << vals_host(i,j,k)
                             << " but reference D2 component: " << basisGrads[l] << "\n";
                }
              }
            }
          }
        }

        // Check D3 of basis function
        {
          DynRankView vals = DynRankView("vals", numFields, numPoints, D3Cardin);
          wedgeBasis.getValues(vals, wedgeNodes, OPERATOR_D3);
          auto vals_host = Kokkos::create_mirror_view(typename HostSpaceType::memory_space(), vals);
          Kokkos::deep_copy(vals_host, vals);
          for (ordinal_type i = 0; i < numFields; ++i) {
            for (ordinal_type j = 0; j < numPoints; ++j) {
              for (ordinal_type k = 0; k < D3Cardin; ++k) {

                // basisGrads is (F,P,D), compute offset:
                const ordinal_type l = k + j * D3Cardin + i * D3Cardin * numPoints;
                if (std::abs(vals_host(i,j,k) - basisD3[l]) > tol) {
                  errorFlag++;
                  *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";

                  // Output the multi-index of the value where the error is:
                  *outStream << " At multi-index { ";
                  *outStream << i << " ";*outStream << j << " ";*outStream << k << " ";
                  *outStream << "}  computed D3 component: " << vals_host(i,j,k)
                             << " but reference D3 component: " << basisD3[l] << "\n";
                }
              }
            }
          }
        }

        // Check D4 of basis function
        {
          DynRankView vals = DynRankView("vals", numFields, numPoints, D4Cardin);
          wedgeBasis.getValues(vals, wedgeNodes, OPERATOR_D4);
          auto vals_host = Kokkos::create_mirror_view(typename HostSpaceType::memory_space(), vals);
          Kokkos::deep_copy(vals_host, vals);
          for (ordinal_type i = 0; i < numFields; ++i) {
            for (ordinal_type j = 0; j < numPoints; ++j) {
              for (ordinal_type k = 0; k < D4Cardin; ++k) {

                // basisGrads is (F,P,D), compute offset:
                const ordinal_type l = k + j * D4Cardin + i * D4Cardin * numPoints;
                if (std::abs(vals_host(i,j,k) - basisD4[l]) > tol) {
                  errorFlag++;
                  *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";

                  // Output the multi-index of the value where the error is:
                  *outStream << " At multi-index { ";
                  *outStream << i << " ";*outStream << j << " ";*outStream << k << " ";
                  *outStream << "}  computed D4 component: " << vals_host(i,j,k)
                             << " but reference D4 component: " << basisD4[l] << "\n";
                }
              }
            }
          }
        }

        // Check all higher derivatives - must be zero.
        {
          const EOperator ops[] = { OPERATOR_D5,
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

            wedgeBasis.getValues(vals, wedgeNodes, op);
            auto vals_host = Kokkos::create_mirror_view(typename HostSpaceType::memory_space(), vals);
            Kokkos::deep_copy(vals_host, vals);
            for (ordinal_type i1=0;i1<numFields; i1++)
              for (ordinal_type i2=0;i2<numPoints; i2++)
                for (ordinal_type i3=0;i3<DkCardin; i3++) {
                  if (std::abs(vals_host(i1,i2,i3)) > tol) {
                    errorFlag++;
                    *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
                    
                    // Get the multi-index of the value where the error is and the operator order
                    const auto ord = Intrepid2::getOperatorOrder(op);
                    *outStream << " At multi-index { "<<i1<<" "<<i2 <<" "<<i3;
                    *outStream << "}  computed D"<< ord <<" component: " << vals_host(i1,i2,i3)
                               << " but reference D" << ord << " component:  0 \n";
                  }
                }
          }
        }
      } catch (std::exception err) {
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
