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
    \brief  Unit tests for the Intrepid2::HGRAD_HEX_C2_FEM class.
    \author Created by P. Bochev, D. Ridzal, K. Peterson and Kyungjoo Kim
*/
#include "Intrepid2_config.h"

#ifdef HAVE_INTREPID2_DEBUG
#define INTREPID2_TEST_FOR_DEBUG_ABORT_OVERRIDE_TO_CONTINUE
#endif

#include "Intrepid2_HGRAD_HEX_C2_FEM.hpp"

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
    int HGRAD_HEX_C2_FEM_Test01(const bool verbose) {

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
        << "\n"
        << "===============================================================================\n"
        << "|                                                                             |\n"
        << "|                 Unit Test (Basis_HGRAD_HEX_C2_FEM)                          |\n"
        << "|                                                                             |\n"
        << "|     1) Conversion of Dof tags into Dof ordinals and back                    |\n"
        << "|     2) Basis values for VALUE, GRAD, and Dk operators                       |\n"
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
      typedef Kokkos::DynRankView<ValueType,HostSpaceType>   DynRankViewHost;
#define ConstructWithLabel(obj, ...) obj(#obj, __VA_ARGS__)

      const ValueType tol = tolerence();
      int errorFlag = 0;

      // for virtual function, value and point types are declared in the class
      typedef ValueType outputValueType;
      typedef ValueType pointValueType;
      Basis_HGRAD_HEX_C2_FEM<DeviceSpaceType,outputValueType,pointValueType> hexBasis;
      //typedef typename decltype(hexBasis)::outputViewType outputViewType;
      //typedef typename decltype(hexBasis)::pointViewType  pointViewType;

      *outStream
        << "\n"
        << "===============================================================================\n"
        << "| TEST 1: Basis creation, exception testing                                   |\n"
        << "===============================================================================\n";

      try{
        ordinal_type nthrow = 0, ncatch = 0;
#ifdef HAVE_INTREPID2_DEBUG

        // Define array containing the 8 vertices of the reference HEX, its center, 12 edge nodes and 6 face centers
        DynRankView ConstructWithLabel( hexNodes, 27, 3);

        // Generic array for the output values; needs to be properly resized depending on the operator type
        const ordinal_type numFields = hexBasis.getCardinality();
        const ordinal_type numPoints = hexNodes.extent(0);
        const ordinal_type spaceDim  = hexBasis.getBaseCellTopology().getDimension();
        const ordinal_type D2Cardin  = getDkCardinality(OPERATOR_D2, spaceDim);

        const ordinal_type workSize  = numFields*numPoints*D2Cardin;
        DynRankView ConstructWithLabel(work, workSize);

        DynRankView vals(work.data(), numFields, numPoints);
        {
          // exception #1: CURL cannot be applied to scalar functions in 3D
          // resize vals to rank-3 container with dimensions (num. basis functions, num. points, arbitrary)
          DynRankView tmpvals = DynRankView(work.data(), numFields, numPoints, 4);
          INTREPID2_TEST_ERROR_EXPECTED( hexBasis.getValues(tmpvals, hexNodes, OPERATOR_CURL) );
        }
        {
          // exception #2: DIV cannot be applied to scalar functions in 3D
          // resize vals to rank-2 container with dimensions (num. basis functions, num. points)
          DynRankView tmpvals = DynRankView(work.data(), numFields, numPoints);
          INTREPID2_TEST_ERROR_EXPECTED( hexBasis.getValues(tmpvals, hexNodes, OPERATOR_DIV) );
        }
        // Exceptions 3-7: all bf tags/bf Ids below are wrong and should cause getDofOrdinal() and
        // getDofTag() to access invalid array elements thereby causing bounds check exception
        {
          // exception #3
          INTREPID2_TEST_ERROR_EXPECTED( hexBasis.getDofOrdinal(3,10,0) );
          // exception #4
          INTREPID2_TEST_ERROR_EXPECTED( hexBasis.getDofOrdinal(1,2,1) );
          // exception #5
          INTREPID2_TEST_ERROR_EXPECTED( hexBasis.getDofOrdinal(0,4,1) );
          // exception #6
          INTREPID2_TEST_ERROR_EXPECTED( hexBasis.getDofTag(numFields) );
          // exception #7
          INTREPID2_TEST_ERROR_EXPECTED( hexBasis.getDofTag(-1) );
        }
        // Exceptions 8-18 test exception handling with incorrectly dimensioned input/output arrays
        {
          // exception #8: input points array must be of rank-2
          DynRankView ConstructWithLabel(badPoints1, 4, 5, 3);
          INTREPID2_TEST_ERROR_EXPECTED( hexBasis.getValues(vals, badPoints1, OPERATOR_VALUE) );
        }
        {
          // exception #9 dimension 1 in the input point array must equal space dimension of the cell
          DynRankView ConstructWithLabel(badPoints2, 4, spaceDim - 1);
          INTREPID2_TEST_ERROR_EXPECTED( hexBasis.getValues(vals, badPoints2, OPERATOR_VALUE) );
        }
        {
          // exception #10 output values must be of rank-2 for OPERATOR_VALUE
          DynRankView ConstructWithLabel(badVals1, 4, 3, 1);
          INTREPID2_TEST_ERROR_EXPECTED( hexBasis.getValues(badVals1, hexNodes, OPERATOR_VALUE) );
        }
        {
          // exception #11 output values must be of rank-3 for OPERATOR_GRAD
          DynRankView ConstructWithLabel(badVals2, 4, 3);
          INTREPID2_TEST_ERROR_EXPECTED( hexBasis.getValues(badVals2, hexNodes, OPERATOR_GRAD) );
          // exception #12 output values must be of rank-3 for OPERATOR_D1
          INTREPID2_TEST_ERROR_EXPECTED( hexBasis.getValues(badVals2, hexNodes, OPERATOR_D1) );
          // exception #13 output values must be of rank-3 for OPERATOR_D2
          INTREPID2_TEST_ERROR_EXPECTED( hexBasis.getValues(badVals2, hexNodes, OPERATOR_D2) );
        }
        {
          // exception #14 incorrect 0th dimension of output array (must equal number of basis functions)
          DynRankView ConstructWithLabel(badVals3, numFields + 1, numPoints);
          INTREPID2_TEST_ERROR_EXPECTED( hexBasis.getValues(badVals3, hexNodes, OPERATOR_VALUE) );
        }
        {
          // exception #15 incorrect 1st dimension of output array (must equal number of points)
          DynRankView ConstructWithLabel(badVals4, numFields, numPoints + 1);
          INTREPID2_TEST_ERROR_EXPECTED( hexBasis.getValues(badVals4, hexNodes, OPERATOR_VALUE) );
        }
        {
          // exception #16: incorrect 2nd dimension of output array (must equal the space dimension)
          DynRankView ConstructWithLabel(badVals5, numFields, numPoints, spaceDim - 1);
          INTREPID2_TEST_ERROR_EXPECTED( hexBasis.getValues(badVals5, hexNodes, OPERATOR_GRAD) );
        }
        {
          // exception #17: incorrect 2nd dimension of output array (must equal D2 cardinality in 3D)
          DynRankView ConstructWithLabel(badVals6, numFields, numPoints, 40);
          INTREPID2_TEST_ERROR_EXPECTED( hexBasis.getValues(badVals6, hexNodes, OPERATOR_D2) );
        }
        {
          // exception #18: incorrect 2nd dimension of output array (must equal D3 cardinality in 3D)
          DynRankView ConstructWithLabel(badVals7, numFields, numPoints, 50);
          INTREPID2_TEST_ERROR_EXPECTED( hexBasis.getValues(badVals7, hexNodes, OPERATOR_D3) );
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
      };

      *outStream
        << "\n"
        << "===============================================================================\n"
        << "| TEST 2: correctness of tag to enum and enum to tag lookups                  |\n"
        << "===============================================================================\n";

      try{

        const ordinal_type numFields = hexBasis.getCardinality();
        const auto allTags = hexBasis.getAllDofTags();

        // Loop over all tags, lookup the associated dof enumeration and then lookup the tag again
        const ordinal_type dofTagSize = allTags.extent(0);
        for (ordinal_type i = 0; i < dofTagSize; ++i) {
          const auto bfOrd  = hexBasis.getDofOrdinal(allTags(i,0), allTags(i,1), allTags(i,2));

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
        for( ordinal_type bfOrd = 0; bfOrd < numFields; ++bfOrd) {
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
        *outStream << err.what() << "\n\n";
        errorFlag = -1000;
      };


      *outStream
        << "\n"
        << "===============================================================================\n"
        << "| TEST 3: correctness of basis function values                                |\n"
        << "===============================================================================\n";

      outStream -> precision(20);

      try{
        // VALUE: Each row gives the 8 correct basis set values at an evaluation point
        const ValueType basisValues[] = {
          1.000, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, \
          0, 0, 0, 0, 0, 0, 1.000, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, \
          0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.000, 0, 0, 0, 0, 0, 0, 0, 0, 0, \
          0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.000, 0, 0, 0, \
          0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, \
          0, 1.000, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, \
          0, 0, 0, 0, 0, 0, 0, 1.000, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, \
          0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.000, 0, 0, 0, 0, 0, 0, 0, 0, \
          0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.000, 0, 0, \
          0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, \
          0, 0, 1.000, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, \
          0, 0, 0, 0, 0, 0, 0, 0, 1.000, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, \
          0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.000, 0, 0, 0, 0, 0, 0, 0, \
          0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.000, 0, \
          0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, \
          0, 0, 0, 1.000, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, \
          0, 0, 0, 0, 0, 0, 0, 0, 0, 1.000, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, \
          0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.000, 0, 0, 0, 0, 0, 0, \
          0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.000, \
          0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, \
          0, 0, 0, 0, 1.000, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, \
          0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.000, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, \
          0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.000, 0, 0, 0, 0, 0, \
          0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, \
          1.000, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, \
          0, 0, 0, 0, 0, 0, 1.000, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, \
          0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.000, 0, 0, 0, 0, 0, 0, 0, 0, 0, \
          0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.000, 0, 0, 0, \
          0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, \
          0, 1.000, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, \
          0, 0, 0, 0, 0, 0, 0, 1.000, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, \
          0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.000, 0, 0, 0, 0, 0, 0, 0, 0, \
          0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.000  };


        // GRAD, D1, D2, D3 and D4 test values are stored in files due to their large size
        std::string     fileName;
        std::ifstream   dataFile;

        // GRAD and D1 values are stored in (F,P,D) format in a data file. Read file and do the test
        std::vector<double> basisGrads;           // Flat array for the gradient values.

        fileName = "../testdata/HEX_C2_GradVals.dat";
        dataFile.open(fileName.c_str());
        INTREPID2_TEST_FOR_EXCEPTION( !dataFile.good(), std::logic_error,
                                      ">>> ERROR (HGRAD_HEX_C2/test01): could not open GRAD values data file, test aborted.");
        while (!dataFile.eof() ){
          double temp;
          std::string line;                            // string for one line of input file
          std::getline(dataFile, line);           // get next line from file
          std::stringstream data_line(line);           // convert to stringstream
          while(data_line >> temp){               // extract value from line
            basisGrads.push_back(temp);           // push into vector
          }
        }
        // It turns out that just closing and then opening the ifstream variable does not reset it
        // and subsequent open() command fails. One fix is to explicitely clear the ifstream, or
        // scope the variables.
        dataFile.close();
        dataFile.clear();


        //D2: flat array with the values of D2 applied to basis functions. Multi-index is (F,P,D2cardinality)
        std::vector<double> basisD2;
        fileName = "../testdata/HEX_C2_D2Vals.dat";
        dataFile.open(fileName.c_str());
        INTREPID2_TEST_FOR_EXCEPTION( !dataFile.good(), std::logic_error,
                                      ">>> ERROR (HGRAD_HEX_C2/test01): could not open D2 values data file, test aborted.");
        while (!dataFile.eof() ){
          double temp;
          std::string line;                            // string for one line of input file
          std::getline(dataFile, line);           // get next line from file
          std::stringstream data_line(line);           // convert to stringstream
          while(data_line >> temp){               // extract value from line
            basisD2.push_back(temp);              // push into vector
          }
        }
        dataFile.close();
        dataFile.clear();


        //D3: flat array with the values of D3 applied to basis functions. Multi-index is (F,P,D3cardinality)
        std::vector<double> basisD3;

        fileName = "../testdata/HEX_C2_D3Vals.dat";
        dataFile.open(fileName.c_str());
        INTREPID2_TEST_FOR_EXCEPTION( !dataFile.good(), std::logic_error,
                                      ">>> ERROR (HGRAD_HEX_C2/test01): could not open D3 values data file, test aborted.");

        while (!dataFile.eof() ){
          double temp;
          std::string line;                            // string for one line of input file
          std::getline(dataFile, line);           // get next line from file
          std::stringstream data_line(line);           // convert to stringstream
          while(data_line >> temp){               // extract value from line
            basisD3.push_back(temp);              // push into vector
          }
        }
        dataFile.close();
        dataFile.clear();


        //D4: flat array with the values of D3 applied to basis functions. Multi-index is (F,P,D4cardinality)
        std::vector<double> basisD4;

        fileName = "../testdata/HEX_C2_D4Vals.dat";
        dataFile.open(fileName.c_str());
        INTREPID2_TEST_FOR_EXCEPTION( !dataFile.good(), std::logic_error,
                                      ">>> ERROR (HGRAD_HEX_C2/test01): could not open D4 values data file, test aborted.");

        while (!dataFile.eof() ){
          double temp;
          std::string line;                            // string for one line of input file
          std::getline(dataFile, line);           // get next line from file
          std::stringstream data_line(line);           // convert to stringstream
          while(data_line >> temp){               // extract value from line
            basisD4.push_back(temp);              // push into vector
          }
        }
        dataFile.close();
        dataFile.clear();

        DynRankViewHost ConstructWithLabel(hexNodesHost, 27, 3);

        // vertices
        hexNodesHost(0, 0) = -1.0;  hexNodesHost(0, 1) = -1.0;  hexNodesHost(0, 2) = -1.0;
        hexNodesHost(1, 0) =  1.0;  hexNodesHost(1, 1) = -1.0;  hexNodesHost(1, 2) = -1.0;
        hexNodesHost(2, 0) =  1.0;  hexNodesHost(2, 1) =  1.0;  hexNodesHost(2, 2) = -1.0;
        hexNodesHost(3, 0) = -1.0;  hexNodesHost(3, 1) =  1.0;  hexNodesHost(3, 2) = -1.0;

        hexNodesHost(4, 0) = -1.0;  hexNodesHost(4, 1) = -1.0;  hexNodesHost(4, 2) =  1.0;
        hexNodesHost(5, 0) =  1.0;  hexNodesHost(5, 1) = -1.0;  hexNodesHost(5, 2) =  1.0;
        hexNodesHost(6, 0) =  1.0;  hexNodesHost(6, 1) =  1.0;  hexNodesHost(6, 2) =  1.0;
        hexNodesHost(7, 0) = -1.0;  hexNodesHost(7, 1) =  1.0;  hexNodesHost(7, 2) =  1.0;

        // nodes on edges
        hexNodesHost(8, 0) =  0.0;   hexNodesHost(8, 1) = -1.0;  hexNodesHost(8, 2) = -1.0;
        hexNodesHost(9, 0) =  1.0;   hexNodesHost(9, 1) =  0.0;  hexNodesHost(9, 2) = -1.0;
        hexNodesHost(10,0) =  0.0;   hexNodesHost(10,1) =  1.0;  hexNodesHost(10,2) = -1.0;
        hexNodesHost(11,0) = -1.0;   hexNodesHost(11,1) =  0.0;  hexNodesHost(11,2) = -1.0;
        hexNodesHost(12,0) = -1.0;   hexNodesHost(12,1) = -1.0;  hexNodesHost(12,2) =  0.0;
        hexNodesHost(13,0) =  1.0;   hexNodesHost(13,1) = -1.0;  hexNodesHost(13,2) =  0.0;
        hexNodesHost(14,0) =  1.0;   hexNodesHost(14,1) =  1.0;  hexNodesHost(14,2) =  0.0;
        hexNodesHost(15,0) = -1.0;   hexNodesHost(15,1) =  1.0;  hexNodesHost(15,2) =  0.0;
        hexNodesHost(16,0) =  0.0;   hexNodesHost(16,1) = -1.0;  hexNodesHost(16,2) =  1.0;
        hexNodesHost(17,0) =  1.0;   hexNodesHost(17,1) =  0.0;  hexNodesHost(17,2) =  1.0;
        hexNodesHost(18,0) =  0.0;   hexNodesHost(18,1) =  1.0;  hexNodesHost(18,2) =  1.0;
        hexNodesHost(19,0) = -1.0;   hexNodesHost(19,1) =  0.0;  hexNodesHost(19,2) =  1.0;

        // center
        hexNodesHost(20,0) =  0.0;  hexNodesHost(20,1) =  0.0;   hexNodesHost(20,2) =  0.0;

        // Face nodes
        hexNodesHost(21,0) =  0.0;   hexNodesHost(21,1) =  0.0;  hexNodesHost(21,2) = -1.0;
        hexNodesHost(22,0) =  0.0;   hexNodesHost(22,1) =  0.0;  hexNodesHost(22,2) =  1.0;
        hexNodesHost(23,0) = -1.0;   hexNodesHost(23,1) =  0.0;  hexNodesHost(23,2) =  0.0;
        hexNodesHost(24,0) =  1.0;   hexNodesHost(24,1) =  0.0;  hexNodesHost(24,2) =  0.0;
        hexNodesHost(25,0) =  0.0;   hexNodesHost(25,1) = -1.0;  hexNodesHost(25,2) =  0.0;
        hexNodesHost(26,0) =  0.0;   hexNodesHost(26,1) =  1.0;  hexNodesHost(26,2) =  0.0;

        auto hexNodes = Kokkos::create_mirror_view(typename DeviceSpaceType::memory_space(), hexNodesHost);
        Kokkos::deep_copy(hexNodes, hexNodesHost);

        // Dimensions for the output arrays:
        const ordinal_type numFields = hexBasis.getCardinality();
        const ordinal_type numPoints = hexNodes.extent(0);
        const ordinal_type spaceDim  = hexBasis.getBaseCellTopology().getDimension();
        const ordinal_type D2Cardin  = getDkCardinality(OPERATOR_D2, spaceDim);
        const ordinal_type D3Cardin  = getDkCardinality(OPERATOR_D3, spaceDim);
        const ordinal_type D4Cardin  = getDkCardinality(OPERATOR_D4, spaceDim);

        {
          // Generic array for values, grads, curls, etc. that will be properly sized before each call
          DynRankView ConstructWithLabel(vals, numFields, numPoints);
          // Check VALUE of basis functions: resize vals to rank-2 container:
          hexBasis.getValues(vals, hexNodes, OPERATOR_VALUE);
          auto vals_host = Kokkos::create_mirror_view(typename HostSpaceType::memory_space(), vals);
          Kokkos::deep_copy(vals_host, vals);
          for (ordinal_type i = 0; i < numFields; ++i) {
            for (ordinal_type j = 0; j < numPoints; ++j) {
              const ordinal_type l =  i + j * numFields;
              if (std::abs(vals_host(i,j) - basisValues[l]) > tol ) {
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
          DynRankView ConstructWithLabel(vals, numFields, numPoints, spaceDim);
          // Check GRAD of basis function: resize vals to rank-3 container
          hexBasis.getValues(vals, hexNodes, OPERATOR_GRAD);
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
          hexBasis.getValues(vals, hexNodes, OPERATOR_D1);
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
          DynRankView ConstructWithLabel(vals, numFields, numPoints, D2Cardin);
          hexBasis.getValues(vals, hexNodes, OPERATOR_D2);
          auto vals_host = Kokkos::create_mirror_view(typename HostSpaceType::memory_space(), vals);
          Kokkos::deep_copy(vals_host, vals);
          for (ordinal_type i = 0; i < numFields; ++i) {
            for (ordinal_type j = 0; j < numPoints; ++j) {
              for (ordinal_type k = 0; k < D2Cardin; ++k) {

                // basisD2 is (F,P,Dk), compute offset:
                const ordinal_type l = k + j * D2Cardin + i * D2Cardin * numPoints;
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
          // Check D3 of basis function
          DynRankView ConstructWithLabel(vals, numFields, numPoints, D3Cardin);
          hexBasis.getValues(vals, hexNodes, OPERATOR_D3);
          auto vals_host = Kokkos::create_mirror_view(typename HostSpaceType::memory_space(), vals);
          Kokkos::deep_copy(vals_host, vals);
          for (ordinal_type i = 0; i < numFields; ++i) {
            for (ordinal_type j = 0; j < numPoints; ++j) {
              for (ordinal_type k = 0; k < D3Cardin; ++k) {

                // basisD3 is (F,P,Dk), compute offset:
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

        {
          // Check D4 of basis function
          DynRankView ConstructWithLabel(vals, numFields, numPoints, D4Cardin);
          hexBasis.getValues(vals, hexNodes, OPERATOR_D4);
          auto vals_host = Kokkos::create_mirror_view(typename HostSpaceType::memory_space(), vals);
          Kokkos::deep_copy(vals_host, vals);
          for (ordinal_type i = 0; i < numFields; i++) {
            for (ordinal_type j = 0; j < numPoints; j++) {
              for (ordinal_type k = 0; k < D4Cardin; k++) {

                // basisD4 is (F,P,Dk), compute offset:
                int l = k + j * D4Cardin + i * D4Cardin * numPoints;
                if (std::abs(vals_host(i,j,k) - basisD4[l]) > tol) {
                  errorFlag++;
                  *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";

                  // Output the multi-index of the value where the error is:
                  *outStream << " At multi-index { ";
                  *outStream << i << " ";*outStream << j << " ";*outStream << k << " ";
                  *outStream << "}  computed D4 component: " << vals_host(i,j,k)
                             << " but reference D4 component: " << basisD2[l] << "\n"; //D2 same as D4?
                }
              }
            }
          }
        }

        {
          // Check D7 to D10 - must be zero. This basis does not support D5 and D6

          const EOperator ops[] = { OPERATOR_D7,
                                    OPERATOR_D8,
                                    OPERATOR_D9,
                                    OPERATOR_D10,
                                    OPERATOR_MAX };
          for (auto h=0;ops[h]!=OPERATOR_MAX;++h) {
            const auto op = ops[h];
            // The last dimension is the number of kth derivatives and needs to be resized for every Dk
            const ordinal_type DkCardin  = getDkCardinality(op, spaceDim);
            DynRankView ConstructWithLabel(vals, numFields, numPoints, DkCardin);
            hexBasis.getValues(vals, hexNodes, op);
            auto vals_host = Kokkos::create_mirror_view(typename HostSpaceType::memory_space(), vals);
            Kokkos::deep_copy(vals_host, vals);

            for (ordinal_type i1 = 0; i1 < numFields; ++i1)
              for (ordinal_type i2 = 0; i2 < numPoints; ++i2)
                for (ordinal_type i3 = 0; i3 < DkCardin; ++i3) {
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
      } catch (std::logic_error err) {
        *outStream << err.what() << "\n\n";
        errorFlag = -1000;
      };

      *outStream
        << "\n"
        << "===============================================================================\n"
        << "| TEST 4: correctness of DoF locations                                        |\n"
        << "===============================================================================\n";

      try{
        const ordinal_type numFields = hexBasis.getCardinality();
        const ordinal_type spaceDim  = hexBasis.getBaseCellTopology().getDimension();

        // Check exceptions.
        ordinal_type nthrow = 0, ncatch = 0;
#ifdef HAVE_INTREPID2_DEBUG
        {
          DynRankView ConstructWithLabel(badVals, 1, 2, 3);
          INTREPID2_TEST_ERROR_EXPECTED( hexBasis.getDofCoords(badVals) );
        }
        {
          DynRankView ConstructWithLabel(badVals, 3, 2);
          INTREPID2_TEST_ERROR_EXPECTED( hexBasis.getDofCoords(badVals) );
        }
        {
          DynRankView ConstructWithLabel(badVals, 27, 2);
          INTREPID2_TEST_ERROR_EXPECTED( hexBasis.getDofCoords(badVals) );
        }
#endif
        // Check if number of thrown exceptions matches the one we expect
        if (nthrow != ncatch) {
          errorFlag++;
          *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
          *outStream << "# of catch ("<< ncatch << ") is different from # of throw (" << ncatch << ")\n";
        }

        DynRankView ConstructWithLabel(bvals, numFields, numFields);
        DynRankView ConstructWithLabel(cvals, numFields, spaceDim);

        // Check mathematical correctness.
        hexBasis.getDofCoords(cvals);
        hexBasis.getValues(bvals, cvals, OPERATOR_VALUE);
        auto cvals_host = Kokkos::create_mirror_view(typename HostSpaceType::memory_space(), cvals);
        Kokkos::deep_copy(cvals_host, cvals);
        auto bvals_host = Kokkos::create_mirror_view(typename HostSpaceType::memory_space(), bvals);
        Kokkos::deep_copy(bvals_host, bvals);
        char buffer[120];
        ordinal_type b_dim0(bvals.extent(0)), b_dim1(bvals.extent(1));
        for (ordinal_type i=0; i<b_dim0; ++i) {
          for (ordinal_type j=0; j<b_dim1; ++j) {
            if ((i != j) && (std::abs(bvals_host(i,j) - 0.0) > tol)) {
              errorFlag++;
              sprintf(buffer, "\nValue of basis function %d at (%6.4e, %6.4e, %6.4e) is %6.4e but should be %6.4e!\n", i, cvals_host(i,0), cvals_host(i,1), cvals_host(i,2), bvals_host(i,j), 0.0);
              *outStream << buffer;
            }
            else if ((i == j) && (std::abs(bvals_host(i,j) - 1.0) > tol)) {
              errorFlag++;
              sprintf(buffer, "\nValue of basis function %d at (%6.4e, %6.4e, %6.4e) is %6.4e but should be %6.4e!\n", i, cvals_host(i,0), cvals_host(i,1), cvals_host(i,2), bvals_host(i,j), 1.0);
              *outStream << buffer;
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

  } //end namespace
} //end namespace
