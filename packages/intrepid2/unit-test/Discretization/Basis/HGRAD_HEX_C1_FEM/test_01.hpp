// @HEADER
// *****************************************************************************
//                           Intrepid2 Package
//
// Copyright 2007 NTESS and the Intrepid2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/** \file test_01.hpp
    \brief  Unit tests for the Intrepid2::HGRAD_HEX_C1_FEM class.
    \author Created by P. Bochev, D. Ridzal, K. Peterson and Kyungjoo Kim
*/
#include "Intrepid2_config.h"

#ifdef HAVE_INTREPID2_DEBUG
#define INTREPID2_TEST_FOR_DEBUG_ABORT_OVERRIDE_TO_CONTINUE
#endif

#include "Intrepid2_HGRAD_HEX_C1_FEM.hpp"

#include "Teuchos_oblackholestream.hpp"
#include "Teuchos_RCP.hpp"

#include "packages/intrepid2/unit-test/Discretization/Basis/Macros.hpp"
#include "packages/intrepid2/unit-test/Discretization/Basis/Setup.hpp"

namespace Intrepid2 {

  namespace Test {

    using HostSpaceType = Kokkos::DefaultHostExecutionSpace;

    template<typename ValueType, typename DeviceType>
    int HGRAD_HEX_C1_FEM_Test01(const bool verbose) {

      //! Create an execution space instance.
      const auto space = Kokkos::Experimental::partition_space(typename DeviceType::execution_space {}, 1)[0];

      //! Setup test output stream.
      Teuchos::RCP<std::ostream> outStream = setup_output_stream<DeviceType>(
        verbose, "Basis_HGRAD_HEX_C1_FEM", {
          "1) Conversion of Dof tags into Dof ordinals and back",
          "2) Basis values for VALUE, GRAD, CURL, and Dk operators"
      });

      Teuchos::oblackholestream oldFormatState;
      oldFormatState.copyfmt(std::cout);

      typedef Kokkos::DynRankView<ValueType,DeviceType> DynRankView;
      typedef Kokkos::DynRankView<ValueType,HostSpaceType>   DynRankViewHost;

      const ValueType tol = tolerence();
      int errorFlag = 0;

      // for virtual function, value and point types are declared in the class
      typedef ValueType outputValueType;
      typedef ValueType pointValueType;
      Basis_HGRAD_HEX_C1_FEM<DeviceType,outputValueType,pointValueType> hexBasis;
      //typedef typename decltype(hexBasis)::OutputViewType OutputViewType;
      //typedef typename decltype(hexBasis)::PointViewType  PointViewType;

      *outStream
        << "\n"
        << "===============================================================================\n"
        << "| TEST 1: Basis creation, exceptions tests                                    |\n"
        << "===============================================================================\n";


      try {
        ordinal_type nthrow = 0, ncatch = 0;
#ifdef HAVE_INTREPID2_DEBUG

        // Define array containing the 8 vertices of the reference HEX, its center and 6 face centers
        DynRankView ConstructWithLabel(hexNodes, 15, 3);

        // Generic array for the output values; needs to be properly resized depending on the operator type
        const auto numFields = hexBasis.getCardinality();
        const auto numPoints = hexNodes.extent(0);
        const auto spaceDim  = hexBasis.getBaseCellTopology().getDimension();
        const auto D2Cardin  = getDkCardinality(OPERATOR_D2, spaceDim);

        const auto workSize  = numFields*numPoints*D2Cardin;
        DynRankView ConstructWithLabel(work, workSize);

        // resize vals to rank-2 container with dimensions
        DynRankView vals = DynRankView(work.data(), numFields, numPoints);
        {
          // exception #1: CURL cannot be applied to scalar functions in 3D
          // resize vals to rank-3 container with dimensions (num. basis functions, num. points, arbitrary)
          DynRankView tmpvals = DynRankView(work.data(), numFields, numPoints, spaceDim);
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
          INTREPID2_TEST_ERROR_EXPECTED( hexBasis.getDofOrdinal(3,0,0) );
          INTREPID2_TEST_ERROR_EXPECTED( hexBasis.getDofOrdinal(1,1,1) );
          INTREPID2_TEST_ERROR_EXPECTED( hexBasis.getDofOrdinal(0,4,1) );
          INTREPID2_TEST_ERROR_EXPECTED( hexBasis.getDofTag(numFields) );
          INTREPID2_TEST_ERROR_EXPECTED( hexBasis.getDofTag(-1)        );
        }

        // Exceptions 8-18 test exception handling with incorrectly dimensioned input/output arrays
        {
          // exception #8: input points array must be of rank-2
          DynRankView ConstructWithLabel(badPoints, 4, 5, 3);
          INTREPID2_TEST_ERROR_EXPECTED( hexBasis.getValues(vals, badPoints, OPERATOR_VALUE) );
        }
        {
          // exception #9 dimension 1 in the input point array must equal space dimension of the cell
          DynRankView ConstructWithLabel(badPoints, 4, 2);
          INTREPID2_TEST_ERROR_EXPECTED( hexBasis.getValues(vals, badPoints, OPERATOR_VALUE) );
        }
        {
          // exception #10 output values must be of rank-2 for OPERATOR_VALUE
          DynRankView ConstructWithLabel(badVals, 4, 3, 1);
          INTREPID2_TEST_ERROR_EXPECTED( hexBasis.getValues(badVals, hexNodes, OPERATOR_VALUE) );
        }
        {
          // exception #11 output values must be of rank-3 for OPERATOR_GRAD
          DynRankView ConstructWithLabel(badVals, 4, 3);
          INTREPID2_TEST_ERROR_EXPECTED( hexBasis.getValues(badVals, hexNodes, OPERATOR_GRAD) );

          // exception #12 output values must be of rank-3 for OPERATOR_D1
          INTREPID2_TEST_ERROR_EXPECTED( hexBasis.getValues(badVals, hexNodes, OPERATOR_D1) );

          // exception #13 output values must be of rank-3 for OPERATOR_D2
          INTREPID2_TEST_ERROR_EXPECTED( hexBasis.getValues(badVals, hexNodes, OPERATOR_D2) );
        }
        {
          // exception #14 incorrect 0th dimension of output array (must equal number of basis functions)
          DynRankView ConstructWithLabel(badVals, numFields + 1, numPoints);
          INTREPID2_TEST_ERROR_EXPECTED( hexBasis.getValues(badVals, hexNodes, OPERATOR_VALUE) );
        }
        {
          // exception #15 incorrect 1st dimension of output array (must equal number of points)
          DynRankView ConstructWithLabel(badVals, numFields, numPoints + 1);
          INTREPID2_TEST_ERROR_EXPECTED( hexBasis.getValues(badVals, hexNodes, OPERATOR_VALUE) );
        }
        {
          // exception #16: incorrect 2nd dimension of output array (must equal the space dimension)
          DynRankView ConstructWithLabel(badVals, numFields, numPoints, 4);
          INTREPID2_TEST_ERROR_EXPECTED( hexBasis.getValues(badVals, hexNodes, OPERATOR_GRAD) );
        }
        {
          // exception #17: incorrect 2nd dimension of output array (must equal D2 cardinality in 3D)
          DynRankView ConstructWithLabel(badVals, numFields, numPoints, 40);
          INTREPID2_TEST_ERROR_EXPECTED( hexBasis.getValues(badVals, hexNodes, OPERATOR_D2) );
        }
        {
          // exception #18: incorrect 2nd dimension of output array (must equal D3 cardinality in 3D)
          DynRankView ConstructWithLabel(badVals, numFields, numPoints, 50);
          INTREPID2_TEST_ERROR_EXPECTED( hexBasis.getValues(badVals, hexNodes, OPERATOR_D3) );
        }
#endif
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
      }

      *outStream
        << "\n"
        << "===============================================================================\n"
        << "| TEST 2: correctness of tag to enum and enum to tag lookups                  |\n"
        << "===============================================================================\n";

      try{

        const auto numFields = hexBasis.getCardinality();
        const auto allTags = hexBasis.getAllDofTags();

        // Loop over all tags, lookup the associated dof enumeration and then lookup the tag again
        const auto dofTagSize = allTags.extent(0);
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
      } catch (std::logic_error &err){
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

      try {
        // VALUE: Each row gives the 8 correct basis set values at an evaluation point
        const ValueType basisValues[][8] = {
          // bottom 4 vertices
          { 1.0, 0.0, 0.0, 0.0,  0.0, 0.0, 0.0, 0.0 },
          { 0.0, 1.0, 0.0, 0.0,  0.0, 0.0, 0.0, 0.0 },
          { 0.0, 0.0, 1.0, 0.0,  0.0, 0.0, 0.0, 0.0 },
          { 0.0, 0.0, 0.0, 1.0,  0.0, 0.0, 0.0, 0.0 },
          // top 4 vertices
          { 0.0, 0.0, 0.0, 0.0,  1.0, 0.0, 0.0, 0.0 },
          { 0.0, 0.0, 0.0, 0.0,  0.0, 1.0, 0.0, 0.0 },
          { 0.0, 0.0, 0.0, 0.0,  0.0, 0.0, 1.0, 0.0 },
          { 0.0, 0.0, 0.0, 0.0,  0.0, 0.0, 0.0, 1.0 },
          // center {0, 0, 0}
          { 0.125, 0.125, 0.125, 0.125,  0.125, 0.125, 0.125, 0.125 },
          // faces { 1, 0, 0} and {-1, 0, 0}
          { 0.0,   0.25,  0.25,  0.0,    0.0,   0.25,  0.25,  0.0 },
          { 0.25,  0.0,   0.0,   0.25,   0.25,  0.0,   0.0,   0.25 },
          // faces { 0, 1, 0} and { 0,-1, 0}
          { 0.0,   0.0,   0.25,  0.25,   0.0,   0.0,   0.25,  0.25 },
          { 0.25,  0.25,  0.0,   0.0,    0.25,  0.25,  0.0,   0.0 },
          // faces {0, 0, 1} and {0, 0, -1}
          { 0.0,   0.0,   0.0,   0.0,    0.25,  0.25,  0.25,  0.25} ,
          { 0.25,  0.25,  0.25,  0.25,   0.0,   0.0,   0.0,   0.0 },
        };

        // GRAD and D1: each row gives the 3x8 correct values of the gradients of the 8 basis functions
        const ValueType basisGrads[][8][3] = {
          // points 0-3
          { { -0.5,-0.5,-0.5 },  { 0.5, 0.0, 0.0 },  { 0.0, 0.0, 0.0 },  { 0.0, 0.5, 0.0 },
            {  0.0, 0.0, 0.5 },  { 0.0, 0.0, 0.0 },  { 0.0, 0.0, 0.0 },  { 0.0, 0.0, 0.0 } },
          { { -0.5, 0.0, 0.0 },  { 0.5,-0.5,-0.5 },  { 0.0, 0.5, 0.0 },  { 0.0, 0.0, 0.0 },
            {  0.0, 0.0, 0.0 },  { 0.0, 0.0, 0.5 },  { 0.0, 0.0, 0.0 },  { 0.0, 0.0, 0.0 } },
          { {  0.0, 0.0, 0.0 },  { 0.0,-0.5, 0.0 },  { 0.5, 0.5,-0.5 },  {-0.5, 0.0, 0.0 },
            {  0.0, 0.0, 0.0 },  { 0.0, 0.0, 0.0 },  { 0.0, 0.0, 0.5 },  { 0.0, 0.0, 0.0 } },
          { {  0.0,-0.5, 0.0 },  { 0.0, 0.0, 0.0 },  { 0.5, 0.0, 0.0 },  {-0.5, 0.5,-0.5 },
            {  0.0, 0.0, 0.0 },  { 0.0, 0.0, 0.0 },  { 0.0, 0.0, 0.0 },  { 0.0, 0.0, 0.5 } },
          // points 4-7
          { {  0.0, 0.0,-0.5 },  { 0.0, 0.0, 0.0 },  { 0.0, 0.0, 0.0 },  { 0.0, 0.0, 0.0 },
            { -0.5,-0.5, 0.5 },  { 0.5, 0.0, 0.0 },  { 0.0, 0.0, 0.0 },  { 0.0, 0.5, 0.0 } },
          { {  0.0, 0.0, 0.0 },  { 0.0, 0.0,-0.5 },  { 0.0, 0.0, 0.0 },  { 0.0, 0.0, 0.0 },
            { -0.5, 0.0, 0.0 },  { 0.5,-0.5, 0.5 },  { 0.0, 0.5, 0.0 },  { 0.0, 0.0, 0.0 } },
          { {  0.0, 0.0, 0.0 },  { 0.0, 0.0, 0.0 },  { 0.0, 0.0,-0.5 },  { 0.0, 0.0, 0.0 },
            {  0.0, 0.0, 0.0 },  { 0.0,-0.5, 0.0 },  { 0.5, 0.5, 0.5 },  {-0.5, 0.0, 0.0 } },
          { {  0.0, 0.0, 0.0 },  { 0.0, 0.0, 0.0 },  { 0.0, 0.0, 0.0 },  { 0.0, 0.0,-0.5 },
            {  0.0,-0.5, 0.0 },  { 0.0, 0.0, 0.0 },  { 0.5, 0.0, 0.0 },  {-0.5, 0.5, 0.5 } },
          // point 8
          { { -0.125,-0.125,-0.125 },  { 0.125,-0.125,-0.125 }, { 0.125, 0.125,-0.125 },  { -0.125, 0.125,-0.125 },
            { -0.125,-0.125, 0.125 },  { 0.125,-0.125, 0.125 }, {  0.125, 0.125, 0.125 }, {-0.125, 0.125, 0.125 } },
          // point 9
          { { -0.125, 0.0,   0.0 },    { 0.125,-0.25, -0.25 },  { 0.125, 0.25, -0.25 },   { -0.125, 0.0, 0.0 },
            { -0.125, 0.0,   0.0 },    { 0.125,-0.25,  0.25 },  { 0.125, 0.25,  0.25 },   { -0.125, 0.0, 0.0 } },
          // point 10
          { { -0.125,-0.25, -0.25 },   { 0.125, 0.0,   0.0 },   { 0.125, 0.0,   0.0 },    { -0.125, 0.25, -0.25 },
            { -0.125,-0.25,  0.25 },   { 0.125, 0.0,   0.0 },   { 0.125, 0.0,   0.0 },    { -0.125, 0.25,  0.25 } },
          // point 11
          { {  0.0,  -0.125, 0.0 },    { 0.0,  -0.125, 0.0 },    { 0.25,  0.125,-0.25 },  { -0.25,  0.125,-0.25 },
            {  0.0,  -0.125, 0.0 },    { 0.0,  -0.125, 0.0 },    { 0.25,  0.125, 0.25 },  { -0.25,  0.125, 0.25 } },
          // point 12
          { { -0.25, -0.125,-0.25 },   { 0.25, -0.125,-0.25 },   { 0.0,   0.125, 0.0 },   { 0.0,   0.125, 0.0 },
            { -0.25, -0.125, 0.25 },   { 0.25, -0.125, 0.25 },   { 0.0,   0.125, 0.0 },   { 0.0,   0.125, 0.0 } },
          // point 13
          { {  0.0,   0.0,  -0.125 },  { 0.0,   0.0,  -0.125 },  { 0.0,   0.0,  -0.125 }, { 0.0,   0.0,  -0.125 },
            { -0.25, -0.25,  0.125 },  { 0.25, -0.25,  0.125 },  { 0.25,  0.25,  0.125 }, { -0.25,  0.25,  0.125 } },
          // point 14
          { { -0.25, -0.25, -0.125 },  { 0.25, -0.25, -0.125 },  { 0.25,  0.25, -0.125 }, { -0.25,  0.25, -0.125 },
            {  0.0,   0.0,   0.125 },  { 0.0,   0.0,   0.125 },  { 0.0,   0.0,   0.125 }, {  0.0,   0.0,   0.125 } }
        };

        //D2: flat array with the values of D2 applied to basis functions. Multi-index is (P,F,K)
        const ValueType basisD2[][8][6] = {
          // pt 0
          { {  0.00000,   0.25000,   0.25000,   0.00000,   0.25000,   0.00000 },
            {  0.00000,  -0.25000,  -0.25000,   0.00000,   0.00000,   0.00000 },
            {  0.00000,   0.25000,   0.00000,   0.00000,   0.00000,   0.00000 },
            {  0.00000,  -0.25000,   0.00000,   0.00000,  -0.25000,   0.00000 },
            {  0.00000,   0.00000,  -0.25000,   0.00000,  -0.25000,   0.00000 },
            {  0.00000,   0.00000,   0.25000,   0.00000,   0.00000,   0.00000 },
            {  0.00000,   0.00000,   0.00000,   0.00000,   0.00000,   0.00000 },
            {  0.00000,   0.00000,   0.00000,   0.00000,   0.25000,   0.00000 } },
          // pt 1
          { {  0.00000,   0.25000,   0.25000,   0.00000,   0.00000,   0.00000 },
            {  0.00000,  -0.25000,  -0.25000,   0.00000,   0.25000,   0.00000 },
            {  0.00000,   0.25000,   0.00000,   0.00000,  -0.25000,   0.00000 },
            {  0.00000,  -0.25000,   0.00000,   0.00000,   0.00000,   0.00000 },
            {  0.00000,   0.00000,  -0.25000,   0.00000,   0.00000,   0.00000 },
            {  0.00000,   0.00000,   0.25000,   0.00000,  -0.25000,   0.00000 },
            {  0.00000,   0.00000,   0.00000,   0.00000,   0.25000,   0.00000 },
            {  0.00000,   0.00000,   0.00000,   0.00000,   0.00000,   0.00000 } },
          // Pt 2
          { {  0.00000,   0.25000,   0.00000,   0.00000,   0.00000,   0.00000 },
            {  0.00000,  -0.25000,   0.00000,   0.00000,   0.25000,   0.00000 },
            {  0.00000,   0.25000,  -0.25000,   0.00000,  -0.25000,   0.00000 },
            {  0.00000,  -0.25000,   0.25000,   0.00000,   0.00000,   0.00000 },
            {  0.00000,   0.00000,   0.00000,   0.00000,   0.00000,   0.00000 },
            {  0.00000,   0.00000,   0.00000,   0.00000,  -0.25000,   0.00000 },
            {  0.00000,   0.00000,   0.25000,   0.00000,   0.25000,   0.00000 },
            {  0.00000,   0.00000,  -0.25000,   0.00000,   0.00000,   0.00000 } },
          // Pt 3
          { {  0.00000,   0.25000,   0.00000,   0.00000,   0.25000,   0.00000 },
            {  0.00000,  -0.25000,   0.00000,   0.00000,   0.00000,   0.00000 },
            {  0.00000,   0.25000,  -0.25000,   0.00000,   0.00000,   0.00000 },
            {  0.00000,  -0.25000,   0.25000,   0.00000,  -0.25000,   0.00000 },
            {  0.00000,   0.00000,   0.00000,   0.00000,  -0.25000,   0.00000 },
            {  0.00000,   0.00000,   0.00000,   0.00000,   0.00000,   0.00000 },
            {  0.00000,   0.00000,   0.25000,   0.00000,   0.00000,   0.00000 },
            {  0.00000,   0.00000,  -0.25000,   0.00000,   0.25000,   0.00000 } },
          // Pt 4
          { {  0.00000,   0.00000,   0.25000,   0.00000,   0.25000,   0.00000 },
            {  0.00000,   0.00000,  -0.25000,   0.00000,   0.00000,   0.00000 },
            {  0.00000,   0.00000,   0.00000,   0.00000,   0.00000,   0.00000 },
            {  0.00000,   0.00000,   0.00000,   0.00000,  -0.25000,   0.00000 },
            {  0.00000,   0.25000,  -0.25000,   0.00000,  -0.25000,   0.00000 },
            {  0.00000,  -0.25000,   0.25000,   0.00000,   0.00000,   0.00000 },
            {  0.00000,   0.25000,   0.00000,   0.00000,   0.00000,   0.00000 },
            {  0.00000,  -0.25000,   0.00000,   0.00000,   0.25000,   0.00000 } },
          // Pt 5
          { {  0.00000,   0.00000,   0.25000,   0.00000,   0.00000,   0.00000 },
            {  0.00000,   0.00000,  -0.25000,   0.00000,   0.25000,   0.00000 },
            {  0.00000,   0.00000,   0.00000,   0.00000,  -0.25000,   0.00000 },
            {  0.00000,   0.00000,   0.00000,   0.00000,   0.00000,   0.00000 },
            {  0.00000,   0.25000,  -0.25000,   0.00000,   0.00000,   0.00000 },
            {  0.00000,  -0.25000,   0.25000,   0.00000,  -0.25000,   0.00000 },
            {  0.00000,   0.25000,   0.00000,   0.00000,   0.25000,   0.00000 },
            {  0.00000,  -0.25000,   0.00000,   0.00000,   0.00000,   0.00000 } },
          // Pt 6
          { {  0.00000,   0.00000,   0.00000,   0.00000,   0.00000,   0.00000 },
            {  0.00000,   0.00000,   0.00000,   0.00000,   0.25000,   0.00000 },
            {  0.00000,   0.00000,  -0.25000,   0.00000,  -0.25000,   0.00000 },
            {  0.00000,   0.00000,   0.25000,   0.00000,   0.00000,   0.00000 },
            {  0.00000,   0.25000,   0.00000,   0.00000,   0.00000,   0.00000 },
            {  0.00000,  -0.25000,   0.00000,   0.00000,  -0.25000,   0.00000 },
            {  0.00000,   0.25000,   0.25000,   0.00000,   0.25000,   0.00000 },
            {  0.00000,  -0.25000,  -0.25000,   0.00000,   0.00000,   0.00000 } },
          // Pt 7
          { {  0.00000,   0.00000,   0.00000,   0.00000,   0.25000,   0.00000 },
            {  0.00000,   0.00000,   0.00000,   0.00000,   0.00000,   0.00000 },
            {  0.00000,   0.00000,  -0.25000,   0.00000,   0.00000,   0.00000 },
            {  0.00000,   0.00000,   0.25000,   0.00000,  -0.25000,   0.00000 },
            {  0.00000,   0.25000,   0.00000,   0.00000,  -0.25000,   0.00000 },
            {  0.00000,  -0.25000,   0.00000,   0.00000,   0.00000,   0.00000 },
            {  0.00000,   0.25000,   0.25000,   0.00000,   0.00000,   0.00000 },
            {  0.00000,  -0.25000,  -0.25000,   0.00000,   0.25000,   0.00000 } },
          // Pt 8
          { {  0.00000,   0.12500,   0.12500,   0.00000,   0.12500,   0.00000 },
            {  0.00000,  -0.12500,  -0.12500,   0.00000,   0.12500,   0.00000 },
            {  0.00000,   0.12500,  -0.12500,   0.00000,  -0.12500,   0.00000 },
            {  0.00000,  -0.12500,   0.12500,   0.00000,  -0.12500,   0.00000 },
            {  0.00000,   0.12500,  -0.12500,   0.00000,  -0.12500,   0.00000 },
            {  0.00000,  -0.12500,   0.12500,   0.00000,  -0.12500,   0.00000 },
            {  0.00000,   0.12500,   0.12500,   0.00000,   0.12500,   0.00000 },
            {  0.00000,  -0.12500,  -0.12500,   0.00000,   0.12500,   0.00000 } },
          // Pt 9
          { {  0.00000,   0.12500,   0.12500,   0.00000,   0.00000,   0.00000 },
            {  0.00000,  -0.12500,  -0.12500,   0.00000,   0.25000,   0.00000 },
            {  0.00000,   0.12500,  -0.12500,   0.00000,  -0.25000,   0.00000 },
            {  0.00000,  -0.12500,   0.12500,   0.00000,   0.00000,   0.00000 },
            {  0.00000,   0.12500,  -0.12500,   0.00000,   0.00000,   0.00000 },
            {  0.00000,  -0.12500,   0.12500,   0.00000,  -0.25000,   0.00000 },
            {  0.00000,   0.12500,   0.12500,   0.00000,   0.25000,   0.00000 },
            {  0.00000,  -0.12500,  -0.12500,   0.00000,   0.00000,   0.00000 } },
          // Pt 10
          { {  0.00000,   0.12500,   0.12500,   0.00000,   0.25000,   0.00000 },
            {  0.00000,  -0.12500,  -0.12500,   0.00000,   0.00000,   0.00000 },
            {  0.00000,   0.12500,  -0.12500,   0.00000,   0.00000,   0.00000 },
            {  0.00000,  -0.12500,   0.12500,   0.00000,  -0.25000,   0.00000 },
            {  0.00000,   0.12500,  -0.12500,   0.00000,  -0.25000,   0.00000 },
            {  0.00000,  -0.12500,   0.12500,   0.00000,   0.00000,   0.00000 },
            {  0.00000,   0.12500,   0.12500,   0.00000,   0.00000,   0.00000 },
            {  0.00000,  -0.12500,  -0.12500,   0.00000,   0.25000,   0.00000 } },
          // Pt 11
          { {  0.00000,   0.12500,   0.00000,   0.00000,   0.12500,   0.00000 },
            {  0.00000,  -0.12500,   0.00000,   0.00000,   0.12500,   0.00000 },
            {  0.00000,   0.12500,  -0.25000,   0.00000,  -0.12500,   0.00000 },
            {  0.00000,  -0.12500,   0.25000,   0.00000,  -0.12500,   0.00000 },
            {  0.00000,   0.12500,   0.00000,   0.00000,  -0.12500,   0.00000 },
            {  0.00000,  -0.12500,   0.00000,   0.00000,  -0.12500,   0.00000 },
            {  0.00000,   0.12500,   0.25000,   0.00000,   0.12500,   0.00000 },
            {  0.00000,  -0.12500,  -0.25000,   0.00000,   0.12500,   0.00000 } },
          // Pt 12
          { {  0.00000,   0.12500,   0.25000,   0.00000,   0.12500,   0.00000 },
            {  0.00000,  -0.12500,  -0.25000,   0.00000,   0.12500,   0.00000 },
            {  0.00000,   0.12500,   0.00000,   0.00000,  -0.12500,   0.00000 },
            {  0.00000,  -0.12500,   0.00000,   0.00000,  -0.12500,   0.00000 },
            {  0.00000,   0.12500,  -0.25000,   0.00000,  -0.12500,   0.00000 },
            {  0.00000,  -0.12500,   0.25000,   0.00000,  -0.12500,   0.00000 },
            {  0.00000,   0.12500,   0.00000,   0.00000,   0.12500,   0.00000 },
            {  0.00000,  -0.12500,   0.00000,   0.00000,   0.12500,   0.00000 } },
          // Pt 13
          { {  0.00000,   0.00000,   0.12500,   0.00000,   0.12500,   0.00000 },
            {  0.00000,   0.00000,  -0.12500,   0.00000,   0.12500,   0.00000 },
            {  0.00000,   0.00000,  -0.12500,   0.00000,  -0.12500,   0.00000 },
            {  0.00000,   0.00000,   0.12500,   0.00000,  -0.12500,   0.00000 },
            {  0.00000,   0.25000,  -0.12500,   0.00000,  -0.12500,   0.00000 },
            {  0.00000,  -0.25000,   0.12500,   0.00000,  -0.12500,   0.00000 },
            {  0.00000,   0.25000,   0.12500,   0.00000,   0.12500,   0.00000 },
            {  0.00000,  -0.25000,  -0.12500,   0.00000,   0.12500,   0.00000 } },
          // Pt 14
          { {  0.00000,   0.25000,   0.12500,   0.00000,   0.12500,   0.00000 },
            {  0.00000,  -0.25000,  -0.12500,   0.00000,   0.12500,   0.00000 },
            {  0.00000,   0.25000,  -0.12500,   0.00000,  -0.12500,   0.00000 },
            {  0.00000,  -0.25000,   0.12500,   0.00000,  -0.12500,   0.00000 },
            {  0.00000,   0.00000,  -0.12500,   0.00000,  -0.12500,   0.00000 },
            {  0.00000,   0.00000,   0.12500,   0.00000,  -0.12500,   0.00000 },
            {  0.00000,   0.00000,   0.12500,   0.00000,   0.12500,   0.00000 },
            {  0.00000,   0.00000,  -0.12500,   0.00000,   0.12500,   0.00000 } }
        };

        // the only nonzeros for D3 are in the "4" slot of the operator column, with values ±1/8, as indicated below:
        const ValueType basisD3Nonzeros[8] = { -0.125, 0.125,-0.125, 0.125, 0.125,-0.125,0.125,-0.125  };


        // Define array containing the 8 vertices of the reference HEX, its center and 6 face centers
        DynRankViewHost ConstructWithLabel(hexNodesHost, 15, 3);

        hexNodesHost(0,0) = -1.0;  hexNodesHost(0,1) = -1.0;  hexNodesHost(0,2) = -1.0;
        hexNodesHost(1,0) =  1.0;  hexNodesHost(1,1) = -1.0;  hexNodesHost(1,2) = -1.0;
        hexNodesHost(2,0) =  1.0;  hexNodesHost(2,1) =  1.0;  hexNodesHost(2,2) = -1.0;
        hexNodesHost(3,0) = -1.0;  hexNodesHost(3,1) =  1.0;  hexNodesHost(3,2) = -1.0;

        hexNodesHost(4,0) = -1.0;  hexNodesHost(4,1) = -1.0;  hexNodesHost(4,2) =  1.0;
        hexNodesHost(5,0) =  1.0;  hexNodesHost(5,1) = -1.0;  hexNodesHost(5,2) =  1.0;
        hexNodesHost(6,0) =  1.0;  hexNodesHost(6,1) =  1.0;  hexNodesHost(6,2) =  1.0;
        hexNodesHost(7,0) = -1.0;  hexNodesHost(7,1) =  1.0;  hexNodesHost(7,2) =  1.0;

        hexNodesHost(8,0) =  0.0;  hexNodesHost(8,1) =  0.0;  hexNodesHost(8,2) =  0.0;

        hexNodesHost(9,0) =  1.0;  hexNodesHost(9,1) =  0.0;  hexNodesHost(9,2) =  0.0;
        hexNodesHost(10,0)= -1.0;  hexNodesHost(10,1)=  0.0;  hexNodesHost(10,2)=  0.0;

        hexNodesHost(11,0)=  0.0;  hexNodesHost(11,1)=  1.0;  hexNodesHost(11,2)=  0.0;
        hexNodesHost(12,0)=  0.0;  hexNodesHost(12,1)= -1.0;  hexNodesHost(12,2)=  0.0;

        hexNodesHost(13,0)=  0.0;  hexNodesHost(13,1)=  0.0;  hexNodesHost(13,2)=  1.0;
        hexNodesHost(14,0)=  0.0;  hexNodesHost(14,1)=  0.0;  hexNodesHost(14,2)= -1.0;

        auto hexNodes = Kokkos::create_mirror_view(typename DeviceType::memory_space(), hexNodesHost);
        Kokkos::deep_copy(hexNodes, hexNodesHost);

        // Generic array for the output values; needs to be properly resized depending on the operator type
        const auto numFields = hexBasis.getCardinality();
        const auto numPoints = hexNodes.extent(0);
        const auto spaceDim  = hexBasis.getBaseCellTopology().getDimension();
        const auto D2Cardin  = getDkCardinality(OPERATOR_D2, spaceDim);
        const auto D3Cardin  = getDkCardinality(OPERATOR_D3, spaceDim);
        const auto D10Cardin = getDkCardinality(OPERATOR_D10, spaceDim);

        const auto workSize  = numFields*numPoints*D10Cardin;
        DynRankView ConstructWithLabel(work, workSize);


        // Check VALUE of basis functions: resize vals to rank-2 container:
        {
          DynRankView vals = DynRankView(work.data(), numFields, numPoints);
          hexBasis.getValues(space, vals, hexNodes, OPERATOR_VALUE);
          auto vals_host = Kokkos::create_mirror_view(typename HostSpaceType::memory_space(), vals);
          Kokkos::deep_copy(vals_host, vals);
          for (auto i=0;i<numFields;++i)
            for (size_type j=0;j<numPoints;++j)
              if (std::abs(vals_host(i,j) - basisValues[j][i]) > tol) {
                errorFlag++;
                *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";

                // Output the multi-index of the value where the error is:
                *outStream << " At multi-index { ";
                *outStream << i << " ";*outStream << j << " ";
                *outStream << "}  computed value: " << vals_host(i,j)
                           << " but reference value: " << basisValues[j][i] << "\n";
              }
        }


        // Check GRAD of basis function: resize vals to rank-3 container
        {
          const EOperator ops[] = { OPERATOR_GRAD,
                                    OPERATOR_D1,
                                    OPERATOR_MAX };
          for (auto h=0;ops[h]!=OPERATOR_MAX;++h) {
            const auto op = ops[h];
            DynRankView vals = DynRankView(work.data(), numFields, numPoints, spaceDim);
            hexBasis.getValues(space, vals, hexNodes, op);
            auto vals_host = Kokkos::create_mirror_view(typename HostSpaceType::memory_space(), vals);
            Kokkos::deep_copy(vals_host, vals);
            for (auto i=0;i<numFields;++i)
              for (size_type j=0;j<numPoints;++j)
                for (size_type k=0;k<spaceDim;++k)
                  if (std::abs(vals_host(i,j,k) - basisGrads[j][i][k]) > tol) {
                    errorFlag++;
                    *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";

                    // Output the multi-index of the value where the error is:
                    *outStream << " At multi-index { ";
                    *outStream << i << " ";*outStream << j << " ";*outStream << k << " ";
                    *outStream << "}  computed grad component: " << vals_host(i,j,k)
                               << " but reference grad component: " << basisGrads[j][i][k] << "\n";
                  }
          }
        }


        // Check D2 of basis function
        {
          DynRankView vals = DynRankView(work.data(), numFields, numPoints, D2Cardin);
          hexBasis.getValues(space, vals, hexNodes, OPERATOR_D2);
          auto vals_host = Kokkos::create_mirror_view(typename HostSpaceType::memory_space(), vals);
          Kokkos::deep_copy(vals_host, vals);
          for (auto i=0;i<numFields;++i)
            for (size_type j=0;j<numPoints;++j)
              for (auto k=0;k<D2Cardin;++k)
                if (std::abs(vals_host(i,j,k) - basisD2[j][i][k]) > tol) {
                  errorFlag++;
                  *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";

                  // Output the multi-index of the value where the error is:
                  *outStream << " At multi-index { ";
                  *outStream << i << " ";*outStream << j << " ";*outStream << k << " ";
                  *outStream << "}  computed D2 component: " << vals_host(i,j,k)
                             << " but reference D2 component: " << basisD2[j][i][k] << "\n";
                }
        }

        // Check D3 of basis function
        {
          DynRankView vals = DynRankView(work.data(), numFields, numPoints, D3Cardin);
          hexBasis.getValues(space, vals, hexNodes, OPERATOR_D3);
          auto vals_host = Kokkos::create_mirror_view(typename HostSpaceType::memory_space(), vals);
          Kokkos::deep_copy(vals_host, vals);
          for (auto i=0;i<numFields;++i)
            for (size_type j=0;j<numPoints;++j)
              for (auto k=0;k<D3Cardin;++k)
              {
                const ValueType expected_value = (k==4) ? basisD3Nonzeros[i] : 0.0;
                if (std::abs(vals_host(i,j,k) - expected_value) > tol) {
                  errorFlag++;
                  *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";

                  // Output the multi-index of the value where the error is:
                  *outStream << " At multi-index { ";
                  *outStream << i << " ";*outStream << j << " ";*outStream << k << " ";
                  *outStream << "}  computed D3 component: " << vals_host(i,j,k)
                             << " but reference D3 component: " << expected_value << "\n";
                }
              }
        }


        // Check all higher derivatives - must be zero.
        {
          const EOperator ops[] = { OPERATOR_D4,
                                    OPERATOR_D5,
                                    OPERATOR_D6,
                                    OPERATOR_D7,
                                    OPERATOR_D8,
                                    OPERATOR_D9,
                                    OPERATOR_D10,
                                    OPERATOR_MAX };
          for (auto h=0;ops[h]!=OPERATOR_MAX;++h) {
            const auto op = ops[h];
            const auto DkCardin  = getDkCardinality(op, spaceDim);
            DynRankView vals = DynRankView(work.data(), numFields, numPoints, DkCardin);

            hexBasis.getValues(space, vals, hexNodes, op);
            auto vals_host = Kokkos::create_mirror_view(typename HostSpaceType::memory_space(), vals);
            Kokkos::deep_copy(vals_host, vals);
            for (auto i1=0;i1<numFields;++i1)
              for (size_type i2=0;i2<numPoints;++i2)
                for (auto i3=0;i3<DkCardin;++i3)
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
      } catch (std::logic_error &err) {
        *outStream << "UNEXPECTED ERROR !!! ----------------------------------------------------------\n";
        *outStream << err.what() << '\n';
        *outStream << "-------------------------------------------------------------------------------" << "\n\n";
        errorFlag = -1000;
      };

      *outStream
        << "\n"
        << "===============================================================================\n"
        << "| TEST 4: correctness of DoF locations                                        |\n"
        << "===============================================================================\n";

      try{
        const auto numFields = hexBasis.getCardinality();
        const auto spaceDim  = hexBasis.getBaseCellTopology().getDimension();

        // Check exceptions.
        ordinal_type nthrow = 0, ncatch = 0;
#ifdef HAVE_INTREPID2_DEBUG
        {
          DynRankView ConstructWithLabel(badVals, 1,2,3);
          INTREPID2_TEST_ERROR_EXPECTED( hexBasis.getDofCoords(badVals) );
        }
        {
          DynRankView ConstructWithLabel(badVals, 3,2);
          INTREPID2_TEST_ERROR_EXPECTED( hexBasis.getDofCoords(badVals) );
        }
        {
          DynRankView ConstructWithLabel(badVals, 8,2);
          INTREPID2_TEST_ERROR_EXPECTED( hexBasis.getDofCoords(badVals) );
        }
#endif
        if (nthrow != ncatch) {
          errorFlag++;
          *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
          *outStream << "# of catch ("<< ncatch << ") is different from # of throw (" << ncatch << ")\n";
        }

        DynRankView ConstructWithLabel(bvals_dev, numFields, numFields);
        DynRankView ConstructWithLabel(cvals_dev, numFields, spaceDim);

        // Check mathematical correctness.
        hexBasis.getDofCoords(cvals_dev);
        hexBasis.getValues(space, bvals_dev, cvals_dev, OPERATOR_VALUE);

        auto bvals = Kokkos::create_mirror_view(typename HostSpaceType::memory_space(), bvals_dev);
        Kokkos::deep_copy(bvals, bvals_dev);
        auto cvals = Kokkos::create_mirror_view(typename HostSpaceType::memory_space(), cvals_dev);
        Kokkos::deep_copy(cvals, cvals_dev);
        for (auto i=0;i<numFields;++i) {
          for (auto j=0;j<numFields;++j) {
            if (i != j && (std::abs(bvals(i,j) - 0.0) > tol)) {
              errorFlag++;
              std::stringstream ss;
              ss << "\n Value of basis function " << i << " at (" << cvals(i,0) << ", " << cvals(i,1) << ") is " << bvals(i,j) << " but should be 0.0\n";
              *outStream << ss.str();
            }
            else if ((i == j) && (std::abs(bvals(i,j) - 1.0) > tol)) {
              errorFlag++;
              std::stringstream ss;
              ss << "\n Value of basis function " << i << " at (" << cvals(i,0) << ", " << cvals(i,1) << ") is " << bvals(i,j) << " but should be 1.0\n";
              *outStream << ss.str();
            }
          }
        }

      } catch (std::logic_error &err){
        *outStream << "UNEXPECTED ERROR !!! ----------------------------------------------------------\n";
        *outStream << err.what() << '\n';
        *outStream << "-------------------------------------------------------------------------------" << "\n\n";
        errorFlag = -1000;
      };

      *outStream
      << "\n"
      << "===============================================================================\n"
      << "| TEST 5: Function Space is Correct                                           |\n"
      << "===============================================================================\n";

      try {
        const EFunctionSpace fs = hexBasis.getFunctionSpace();

        if (fs != FUNCTION_SPACE_HGRAD)
        {
          *outStream << std::setw(70) << "------------- TEST FAILURE! -------------" << "\n";

          // Output the multi-index of the value where the error is:
          *outStream << " Expected a function space of FUNCTION_SPACE_HGRAD (enum value " << FUNCTION_SPACE_HGRAD << "),";
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
  } // end of test
} // end of intrepid2
