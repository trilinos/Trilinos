// @HEADER
// *****************************************************************************
//                           Intrepid2 Package
//
// Copyright 2007 NTESS and the Intrepid2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/** \file test_01.cpp
    \brief  Unit tests for the Intrepid2::Basis_HGRAD_TRI_C2_FEM class.
    \author Created by P. Bochev, D. Ridzal, K. Peterson and Kyungjoo Kim.
*/

#include "Intrepid2_config.h"

#ifdef HAVE_INTREPID2_DEBUG
#define INTREPID2_TEST_FOR_DEBUG_ABORT_OVERRIDE_TO_CONTINUE
#endif

#include "Intrepid2_HGRAD_TRI_C2_FEM.hpp"

#include "Teuchos_oblackholestream.hpp"
#include "Teuchos_RCP.hpp"

#include "packages/intrepid2/unit-test/Discretization/Basis/Macros.hpp"
#include "packages/intrepid2/unit-test/Discretization/Basis/Setup.hpp"

namespace Intrepid2 {

  namespace Test {

    using HostSpaceType = Kokkos::DefaultHostExecutionSpace;

    template<typename ValueType, typename DeviceType>
    int HGRAD_TRI_C2_FEM_Test01(const bool verbose) {

      //! Create an execution space instance.
      const auto space = Kokkos::Experimental::partition_space(typename DeviceType::execution_space {}, 1)[0];

      //! Setup test output stream.
      Teuchos::RCP<std::ostream> outStream = setup_output_stream<DeviceType>(
        verbose, "Basis_HGRAD_TRI_C2_FEM", {
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
      Basis_HGRAD_TRI_C2_FEM<DeviceType,outputValueType,pointValueType> triBasis;
      //typedef typename decltype(triBasis)::OutputViewType OutputViewType;
      //typedef typename decltype(triBasis)::PointViewType  PointViewType;

      *outStream
        << "\n"
        << "===============================================================================\n"
        << "| TEST 1: Basis creation, exception testing                                   |\n"
        << "===============================================================================\n";

      try{
        ordinal_type nthrow = 0, ncatch = 0;
#ifdef HAVE_INTREPID2_DEBUG

        DynRankView ConstructWithLabel(triNodes, 7, 2);

        const ordinal_type numFields = triBasis.getCardinality();
        const ordinal_type numPoints = triNodes.extent(0);
        const ordinal_type spaceDim  = triBasis.getBaseCellTopology().getDimension();
        //const ordinal_type D2Cardin  = getDkCardinality(OPERATOR_D2, spaceDim);

        //const auto workSize = numFields*numPoints*D2Cardin;

        DynRankView ConstructWithLabel(vals, numFields, numPoints);
        {
          // exception #1: DIV cannot be applied to scalar functions
          // resize vals to rank-2 container with dimensions (num. points, num. basis functions)
          INTREPID2_TEST_ERROR_EXPECTED( triBasis.getValues(vals, triNodes, OPERATOR_DIV) );
        }
        // Exceptions 2-6: all bf tags/bf Ids below are wrong and should cause getDofOrdinal() and 
        // getDofTag() to access invalid array elements thereby causing bounds check exception
        {
          // exception #2
          INTREPID2_TEST_ERROR_EXPECTED( triBasis.getDofOrdinal(2,0,0) );
          // exception #3
          INTREPID2_TEST_ERROR_EXPECTED( triBasis.getDofOrdinal(1,1,1) );
          // exception #4
          INTREPID2_TEST_ERROR_EXPECTED( triBasis.getDofOrdinal(0,4,0) );
          // exception #5
          INTREPID2_TEST_ERROR_EXPECTED( triBasis.getDofTag(numFields) );
          // exception #6
          INTREPID2_TEST_ERROR_EXPECTED( triBasis.getDofTag(-1) );
        }
        // Exceptions 7-17 test exception handling with incorrectly dimensioned input/output arrays
        {
          // exception #7: input points array must be of rank-2
          DynRankView ConstructWithLabel(badPoints1, 4, 5, 3);
          INTREPID2_TEST_ERROR_EXPECTED( triBasis.getValues(vals, badPoints1, OPERATOR_VALUE) );
        }
        {
          // exception #8 dimension 1 in the input point array must equal space dimension of the cell
          DynRankView ConstructWithLabel(badPoints2, 4, spaceDim + 1);
          INTREPID2_TEST_ERROR_EXPECTED( triBasis.getValues(vals, badPoints2, OPERATOR_VALUE) );
        }
        { 
          // exception #9 output values must be of rank-2 for OPERATOR_VALUE
          DynRankView ConstructWithLabel(badVals1, 4, 3, 1);
          INTREPID2_TEST_ERROR_EXPECTED( triBasis.getValues(badVals1, triNodes, OPERATOR_VALUE) );
        }
        { 
          // exception #10 output values must be of rank-3 for OPERATOR_GRAD
          DynRankView ConstructWithLabel(badVals2, 4, 3);
          INTREPID2_TEST_ERROR_EXPECTED( triBasis.getValues(badVals2, triNodes, OPERATOR_GRAD) );
          // exception #11 output values must be of rank-3 for OPERATOR_CURL
          INTREPID2_TEST_ERROR_EXPECTED( triBasis.getValues(badVals2, triNodes, OPERATOR_CURL) );
          // exception #12 output values must be of rank-3 for OPERATOR_D2
          INTREPID2_TEST_ERROR_EXPECTED( triBasis.getValues(badVals2, triNodes, OPERATOR_D2) );
        }
        { 
          // exception #13 incorrect 1st dimension of output array (must equal number of basis functions)
          DynRankView ConstructWithLabel(badVals3, numFields + 1, numPoints);
          INTREPID2_TEST_ERROR_EXPECTED( triBasis.getValues(badVals3, triNodes, OPERATOR_VALUE) );
        }
        { 
          // exception #14 incorrect 0th dimension of output array (must equal number of points)
          DynRankView ConstructWithLabel(badVals4, numFields, numPoints + 1);
          INTREPID2_TEST_ERROR_EXPECTED( triBasis.getValues(badVals4, triNodes, OPERATOR_VALUE) );
        }
        { 
          // exception #15: incorrect 2nd dimension of output array (must equal the space dimension)
          DynRankView ConstructWithLabel(badVals5, numFields, numPoints, spaceDim + 1);
          INTREPID2_TEST_ERROR_EXPECTED( triBasis.getValues(badVals5, triNodes, OPERATOR_GRAD) );
        }
        { 
          // exception #16: incorrect 2nd dimension of output array (must equal D2 cardinality in 2D)
          DynRankView ConstructWithLabel(badVals6, numFields, numPoints, 40);
          INTREPID2_TEST_ERROR_EXPECTED( triBasis.getValues(badVals6, triNodes, OPERATOR_D2) );
          // exception #17: incorrect 2nd dimension of output array (must equal D3 cardinality in 2D)
          INTREPID2_TEST_ERROR_EXPECTED( triBasis.getValues(badVals6, triNodes, OPERATOR_D3) );
        }
#endif
        // Check if number of thrown exceptions matches the one we expect
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
      };
  
      *outStream
        << "\n"
        << "===============================================================================\n"
        << "| TEST 2: correctness of tag to enum and enum to tag lookups                  |\n"
        << "===============================================================================\n";

      try{
        const ordinal_type numFields = triBasis.getCardinality();
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
        for( ordinal_type bfOrd = 0; bfOrd < numFields; bfOrd++) {
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
      } catch (std::logic_error &err){
        *outStream << err.what() << "\n\n";
        errorFlag = -1000;
      };
  
      *outStream \
        << "\n"
        << "===============================================================================\n"\
        << "| TEST 3: correctness of basis function values                                |\n"\
        << "===============================================================================\n";

      outStream -> precision(20);

      // VALUE: Correct basis values in (F,P) format: each row gives the 6 correct basis values at
      // the lattice points (3 vertices and 3 edge midpoints)
      const ValueType basisValues[] = {
        1.0, 0.0, 0.0,    0.0, 0.0, 0.0,    -0.10710168,
        0.0, 1.0, 0.0,    0.0, 0.0, 0.0,    -0.11320352,
        0.0, 0.0, 1.0,    0.0, 0.0, 0.0,     0.23015592,
        0.0, 0.0, 0.0,    1.0, 0.0, 0.0,     0.10766112,
        0.0, 0.0, 0.0,    0.0, 1.0, 0.0,     0.46514592,
        0.0, 0.0, 0.0,    0.0, 0.0, 1.0,     0.41734224
      };

      // GRAD and D1: Correct gradients and D1 in (F,P,D) format: each row contains 6x2 values of
      // gradients of basis functions. Same array is used to check correctness of CURL.
      const ValueType basisGrads[] = {
        //   V0            V1           V2                   M0             M1            M2          RP
        -3.0, -3.0,     1.0, 1.0,   1.0, 1.0,           -1.0, -1.0,     1.0, 1.0,     -1.0,-1.0,     0.3784,  0.3784,
        -1.0,  0.0,     3.0, 0.0,  -1.0, 0.0,            1.0,  0.0,     1.0, 0.0,     -1.0, 0.0,    -0.3072,  0.0,
        0.0, -1.0,     0.0,-1.0,   0.0, 3.0,            0.0, -1.0,     0.0, 1.0,      0.0, 1.0,     0.0,     1.6856,
        4.0,  0.0,    -4.0,-4.0,   0.0, 0.0,            0.0, -2.0,    -2.0,-2.0,      2.0, 0.0,    -0.0712, -0.6928, 
        0.0,  0.0,     0.0, 4.0,   4.0, 0.0,            0.0,  2.0,     2.0, 2.0,      2.0, 0.0,     2.6856,  0.6928,
        0.0,  4.0,     0.0, 0.0,  -4.0,-4.0,            0.0,  2.0,    -2.0,-2.0,     -2.0, 0.0,    -2.6856, -2.064
      };
  
      // D2: Correct multiset of second order partials in (F,P,Dk) format. D2 cardinality = 3 for 2D 
      const ValueType basisD2[] = {
        //   V0                  V1              V2                   M0              M1                M2            RP
        4.0, 4.0, 4.0,    4.0, 4.0, 4.0,    4.0, 4.0, 4.0,    4.0, 4.0, 4.0,    4.0, 4.0, 4.0,    4.0, 4.0, 4.0,   4.0, 4.0, 4.0,
        4.0, 0.0, 0.0,    4.0, 0.0, 0.0,    4.0, 0.0, 0.0,    4.0, 0.0, 0.0,    4.0, 0.0, 0.0,    4.0, 0.0, 0.0,   4.0, 0.0, 0.0,
        0.0, 0.0, 4.0,    0.0, 0.0, 4.0,    0.0, 0.0, 4.0,    0.0, 0.0, 4.0,    0.0, 0.0, 4.0,    0.0, 0.0, 4.0,   0.0, 0.0, 4.0,
        -8.0,-4.0, 0.0,   -8.0,-4.0, 0.0,   -8.0,-4.0, 0.0,   -8.0,-4.0, 0.0,   -8.0,-4.0, 0.0,   -8.0,-4.0, 0.0,  -8.0,-4.0, 0.0,
        0.0, 4.0, 0.0,    0.0, 4.0, 0.0,    0.0, 4.0, 0.0,    0.0, 4.0, 0.0,    0.0, 4.0, 0.0,    0.0, 4.0, 0.0,   0.0, 4.0, 0.0,
        0.0,-4.0,-8.0,    0.0,-4.0,-8.0,    0.0,-4.0,-8.0,    0.0,-4.0,-8.0,    0.0,-4.0,-8.0,    0.0,-4.0,-8.0,   0.0,-4.0,-8.0
      };

      try{
        DynRankViewHost ConstructWithLabel(triNodesHost, 7, 2); 

        triNodesHost(0,0) =  0.0;     triNodesHost(0,1) =  0.0;
        triNodesHost(1,0) =  1.0;     triNodesHost(1,1) =  0.0;
        triNodesHost(2,0) =  0.0;     triNodesHost(2,1) =  1.0;
        triNodesHost(3,0) =  0.5;     triNodesHost(3,1) =  0.0;
        triNodesHost(4,0) =  0.5;     triNodesHost(4,1) =  0.5;
        triNodesHost(5,0) =  0.0;     triNodesHost(5,1) =  0.5;
        triNodesHost(6,0) =  0.1732;  triNodesHost(6,1) = 0.6714;

        auto triNodes = Kokkos::create_mirror_view(typename DeviceType::memory_space(), triNodesHost);
        Kokkos::deep_copy(triNodes,triNodesHost);

        // Dimensions for the output arrays:
        const ordinal_type numFields = triBasis.getCardinality();
        const ordinal_type numPoints = triNodes.extent(0);
        const ordinal_type spaceDim  = triBasis.getBaseCellTopology().getDimension();
        const ordinal_type D2cardinality = getDkCardinality(OPERATOR_D2,spaceDim);
    
        // Check VALUE of basis functions: resize vals to rank-2 container:
        {
          DynRankView ConstructWithLabel(vals, numFields, numPoints);
          triBasis.getValues(space, vals, triNodes, OPERATOR_VALUE);
          auto vals_host = Kokkos::create_mirror_view(typename HostSpaceType::memory_space(), vals);
          Kokkos::deep_copy(vals_host, vals);
          for (ordinal_type i = 0; i < numFields; ++i) {
            for (ordinal_type j = 0; j < numPoints; ++j) {
        
              // Compute offset for (F,P) container
              const ordinal_type l =  j + i * numPoints;
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
          DynRankView ConstructWithLabel(vals, numFields, numPoints, spaceDim);
          triBasis.getValues(space, vals, triNodes, OPERATOR_GRAD);
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
          triBasis.getValues(space, vals, triNodes, OPERATOR_D1);
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

        // Check CURL of basis function: resize vals just for illustration! 
        {
          DynRankView ConstructWithLabel(vals, numFields, numPoints, spaceDim);
          triBasis.getValues(space, vals, triNodes, OPERATOR_CURL);
          auto vals_host = Kokkos::create_mirror_view(typename HostSpaceType::memory_space(), vals);
          Kokkos::deep_copy(vals_host, vals);
          for (ordinal_type i = 0; i < numFields; ++i) {
            for (ordinal_type j = 0; j < numPoints; ++j) {
              // We will use "rotated" basisGrads to check CURL: get offsets to extract (u_y, -u_x)
              const ordinal_type curl_0 = 1 + j * spaceDim + i * spaceDim * numPoints;               // position of y-derivative
              const ordinal_type curl_1 = 0 + j * spaceDim + i * spaceDim * numPoints;               // position of x-derivative
        
              const auto curl_value_0 = basisGrads[curl_0];
              const auto curl_value_1 =-basisGrads[curl_1];
              if (std::abs(vals_host(i,j,0) - curl_value_0) > tol) {
                errorFlag++;
                *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
                // Output the multi-index of the value where the error is:
                *outStream << " At multi-index { ";
                *outStream << i << " ";*outStream << j << " ";*outStream << 0 << " ";
                *outStream << "}  computed curl component: " << vals_host(i,j,0)
                           << " but reference curl component: " << curl_value_0 << "\n";
              }
              if (std::abs(vals_host(i,j,1) - curl_value_1) > tol) {
                errorFlag++;
                *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
                // Output the multi-index of the value where the error is:
                *outStream << " At multi-index { ";
                *outStream << i << " ";*outStream << j << " ";*outStream << 1 << " ";
                *outStream << "}  computed curl component: " << vals_host(i,j,1)
                           << " but reference curl component: " << curl_value_1 << "\n";
              }
            }
          }
        }
    
        // Check D2 of basis function
        {
          DynRankView ConstructWithLabel(vals, numFields, numPoints, D2cardinality);
          triBasis.getValues(space, vals, triNodes, OPERATOR_D2);
          auto vals_host = Kokkos::create_mirror_view(typename HostSpaceType::memory_space(), vals);
          Kokkos::deep_copy(vals_host, vals);
          for (ordinal_type i = 0; i < numFields; i++) {
            for (ordinal_type j = 0; j < numPoints; j++) {
              for (ordinal_type k = 0; k < D2cardinality; k++) {
          
                // basisD2 is (F,P,Dk), compute offset:
                const ordinal_type l = k + j * D2cardinality + i * D2cardinality * numPoints;
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

            triBasis.getValues(space, vals, triNodes, op);
            auto vals_host = Kokkos::create_mirror_view(typename HostSpaceType::memory_space(), vals);
            Kokkos::deep_copy(vals_host, vals);
            for (ordinal_type i1 = 0; i1 < numFields; i1++)
              for (ordinal_type i2 = 0; i2 < numPoints; i2++)
                for (ordinal_type i3 = 0; i3 < DkCardin; i3++) {
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
      }

      // Catch unexpected errors
      catch (std::logic_error &err) {
        *outStream << err.what() << "\n\n";
        errorFlag = -1000;
      };

      *outStream
      << "\n"
      << "===============================================================================\n"
      << "| TEST 4: Function Space is Correct                                           |\n"
      << "===============================================================================\n";
      
      try {
        const EFunctionSpace fs = triBasis.getFunctionSpace();
        
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

  } //end namespace
} //end namespace
