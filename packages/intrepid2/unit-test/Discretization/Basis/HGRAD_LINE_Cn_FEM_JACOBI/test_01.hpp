// @HEADER
// *****************************************************************************
//                           Intrepid2 Package
//
// Copyright 2007 NTESS and the Intrepid2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/** \file test_01.cpp
    \brief  Unit tests for the Intrepid2::Basis_HGRAD_LINE_Cn_FEM_JACOBI class.
    \author Created by P. Bochev, D. Ridzal, K. Peterson, Kyungjoo Kim.
*/

#include "Intrepid2_config.h"

#ifdef HAVE_INTREPID2_DEBUG
#define INTREPID2_TEST_FOR_DEBUG_ABORT_OVERRIDE_TO_CONTINUE
#endif

#include "Intrepid2_Types.hpp"
#include "Intrepid2_Utils.hpp"

#include "Intrepid2_CubatureDirectLineGauss.hpp"
#include "Intrepid2_FunctionSpaceTools.hpp"
#include "Intrepid2_HGRAD_LINE_Cn_FEM_JACOBI.hpp"

#include "packages/intrepid2/unit-test/Discretization/Basis/Macros.hpp"
#include "packages/intrepid2/unit-test/Discretization/Basis/Setup.hpp"

namespace Intrepid2 {

  namespace Test {

    using HostSpaceType = Kokkos::DefaultHostExecutionSpace;

    template<typename ValueType, typename DeviceType>
    int HGRAD_LINE_Cn_FEM_JACOBI_Test01(const bool verbose) {

      //! Create an execution space instance.
      const auto space = Kokkos::Experimental::partition_space(typename DeviceType::execution_space {}, 1)[0];

      Teuchos::RCP<std::ostream> outStream = setup_output_stream<DeviceType>(
        verbose, "Basis_HGRAD_LINE_Cn_FEM_JACOBI", {
          "1) Conversion of Dof tags into Dof ordinals and back",
          "2) Basis values for VALUE, GRAD, CURL, and Dk operators"
      });

      Teuchos::oblackholestream oldFormatState;
      oldFormatState.copyfmt(std::cout);

      typedef Kokkos::DynRankView<ValueType,DeviceType> DynRankView;
      typedef Kokkos::DynRankView<ValueType,Kokkos::HostSpace>   DynRankViewHost;
      
      const ValueType tol = tolerence();
      int errorFlag = 0;
      
      // for virtual function, value and point types are declared in the class
      typedef ValueType outputValueType;
      typedef ValueType pointValueType;
      typedef ValueType weightValueType;

      typedef Basis_HGRAD_LINE_Cn_FEM_JACOBI<DeviceType,outputValueType,pointValueType> LineBasisType;
      typedef CubatureDirectLineGauss<DeviceType,pointValueType,weightValueType> CubatureLineType;

      constexpr ordinal_type maxOrder = Parameters::MaxOrder;
      
      typedef FunctionSpaceTools<DeviceType> fst;
      
      *outStream
        << "\n"
        << "===============================================================================\n"
        << "| TEST 1: Basis creation, exceptions tests                                    |\n"
        << "===============================================================================\n";
      
      try{

#ifdef HAVE_INTREPID2_DEBUG
        ordinal_type nthrow = 0, ncatch = 0;
        constexpr ordinal_type order = 5;
        if(order <= maxOrder) {
          const double alpha = 0.0, beta = 0.0;

          LineBasisType lineBasis(order, alpha, beta);

          // Define array containing array of nodes to evaluate
          DynRankView ConstructWithLabel(lineNodes, 10, 1);

          // Generic array for the output values; needs to be properly resized depending on the operator type
          const auto numFields = lineBasis.getCardinality();
          const auto numPoints = lineNodes.extent(0);
          //const auto spaceDim  = lineBasis.getBaseCellTopology().getDimension();


          // Exceptions 1-5: all bf tags/bf Ids below are wrong and should cause getDofOrdinal() and
          // getDofTag() to access invalid array elements thereby causing bounds check exception
          {
            INTREPID2_TEST_ERROR_EXPECTED( lineBasis.getDofOrdinal(2,0,0) );
            INTREPID2_TEST_ERROR_EXPECTED( lineBasis.getDofOrdinal(1,1,1) );
            INTREPID2_TEST_ERROR_EXPECTED( lineBasis.getDofOrdinal(1,0,7) );
            INTREPID2_TEST_ERROR_EXPECTED( lineBasis.getDofTag(numFields) );
            INTREPID2_TEST_ERROR_EXPECTED( lineBasis.getDofTag(-1) );

            // no exception; if happens, it is unexpected;
            lineBasis.getDofOrdinal(1,0,5);
            lineBasis.getDofTag(5);
          }

          // Exceptions 6-16 test exception handling with incorrectly dimensioned input/output arrays
          {
            DynRankView ConstructWithLabel(vals, numFields, numPoints);
            {
              // exception #6: input points array must be of rank-2
              DynRankView ConstructWithLabel(badPoints, 4, 5, 3);
              INTREPID2_TEST_ERROR_EXPECTED( lineBasis.getValues(vals, badPoints, OPERATOR_VALUE) );
            }
            {
              // exception #7: dimension 1 in the input point array must equal space dimension of the cell
              DynRankView ConstructWithLabel(badPoints, 4, 3);
              INTREPID2_TEST_ERROR_EXPECTED( lineBasis.getValues(vals, badPoints, OPERATOR_VALUE) );
            }
            {
              // exception #8: output values must be of rank-2 for OPERATOR_VALUE
              DynRankView ConstructWithLabel(badVals, 4, 3, 1);
              INTREPID2_TEST_ERROR_EXPECTED( lineBasis.getValues(badVals, lineNodes, OPERATOR_VALUE) );
            }
            lineBasis.getValues(vals, lineNodes, OPERATOR_VALUE);
          }
          {
            // exception #9: output values must be of rank-3 for OPERATOR_GRAD
            DynRankView ConstructWithLabel(badVals, 4, 3);
            INTREPID2_TEST_ERROR_EXPECTED( lineBasis.getValues(badVals, lineNodes, OPERATOR_GRAD) );

            // exception #10: output values must be of rank-3 for OPERATOR_CURL
            INTREPID2_TEST_ERROR_EXPECTED( lineBasis.getValues(badVals, lineNodes, OPERATOR_CURL) );

            // exception #11: output values must be of rank-2 for OPERATOR_DIV
            INTREPID2_TEST_ERROR_EXPECTED( lineBasis.getValues(badVals, lineNodes, OPERATOR_DIV) );

            // exception #12: output values must be of rank-2 for OPERATOR_D1
            INTREPID2_TEST_ERROR_EXPECTED( lineBasis.getValues(badVals, lineNodes, OPERATOR_D1) );
          }
          {
            // exception #13: incorrect 0th dimension of output array (must equal number of basis functions)
            DynRankView ConstructWithLabel(badVals, lineBasis.getCardinality() + 1, lineNodes.extent(0));
            INTREPID2_TEST_ERROR_EXPECTED( lineBasis.getValues(badVals, lineNodes, OPERATOR_VALUE) );
          }
          {
            // exception #14: incorrect 1st dimension of output array (must equal number of points)
            DynRankView ConstructWithLabel(badVals, lineBasis.getCardinality(), lineNodes.extent(0) + 1);
            INTREPID2_TEST_ERROR_EXPECTED( lineBasis.getValues(badVals, lineNodes, OPERATOR_VALUE) );
          }
          {
            // exception #15: incorrect 2nd dimension of output array (must equal spatial dimension)
            DynRankView ConstructWithLabel(badVals, lineBasis.getCardinality(), lineNodes.extent(0), 2);
            INTREPID2_TEST_ERROR_EXPECTED( lineBasis.getValues(badVals, lineNodes, OPERATOR_GRAD) );
          }
        }
        if (nthrow != ncatch) {
          errorFlag++;
          *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
          *outStream << "# of catch ("<< ncatch << ") is different from # of throw (" << nthrow << ")\n";
        }
#endif
      } catch (std::exception &err) {
        *outStream << "UNEXPECTED ERROR !!! ----------------------------------------------------------\n";
        *outStream << err.what() << '\n';
        *outStream << "-------------------------------------------------------------------------------" << "\n\n";
        errorFlag = -1000;
      };

      *outStream
        << "\n"
        << "===============================================================================\n"
        << "| TEST 2: orthogonality of basis functions                                    |\n"
        << "===============================================================================\n";

      outStream->precision(20);

      try {
        const double alpha = 0.0, beta = 0.0;
        const shards::CellTopology line(shards::getCellTopologyData< shards::Line<> >()); 

        const auto maxp = Parameters::MaxOrder;
        CubatureLineType lineCub(maxp + maxp);
        
        const auto numCubPoints = lineCub.getNumPoints();
        const auto cubDimension = lineCub.getDimension();

        DynRankView ConstructWithLabel(cubPoints, numCubPoints, cubDimension); 
        DynRankView ConstructWithLabel(cubWeightsCell, 1, numCubPoints);

        auto cubWeights = Kokkos::subdynrankview(cubWeightsCell, 0, Kokkos::ALL());
        lineCub.getCubature(cubPoints, cubWeights);        

        for (auto ip=0;ip<maxp;++ip) {
          LineBasisType lineBasisLeft(ip, alpha, beta);

          const auto leftCard = lineBasisLeft.getCardinality();
          DynRankView ConstructWithLabel(valsLeftCell, 1, leftCard, numCubPoints);

          auto valsLeft = Kokkos::subdynrankview(valsLeftCell, 0, Kokkos::ALL(), Kokkos::ALL());
          lineBasisLeft.getValues(valsLeft, cubPoints, OPERATOR_VALUE);

          for (auto jp=0;jp<maxp;++jp) {
            LineBasisType lineBasisRight(jp, alpha, beta);

            const auto rightCard = lineBasisRight.getCardinality();
            DynRankView ConstructWithLabel(valsRightCell, 1, rightCard, numCubPoints);
            
            auto valsRight = Kokkos::subdynrankview(valsRightCell, 0, Kokkos::ALL(), Kokkos::ALL());
            lineBasisRight.getValues(valsRight, cubPoints, OPERATOR_VALUE);

            DynRankView ConstructWithLabel(massMatrix, 1, leftCard, rightCard);

            fst::scalarMultiplyDataField(valsRightCell, cubWeightsCell, valsRightCell);
            fst::integrate(massMatrix, valsLeftCell, valsRightCell);

            // host mirror for comparison
            auto massMatrixHost = Kokkos::create_mirror_view(massMatrix);
            Kokkos::deep_copy(massMatrixHost, massMatrix);
            
            // check orthogonality property
            for (auto i=0;i<leftCard;++i) 
              for (auto j=0;j<rightCard;++j) {
                const ValueType exactVal = (i == j)*(2.0/(2.0*j+1.0));
                const auto val = massMatrixHost(0,i,j);
                if ( std::isnan(val) || std::abs(val-exactVal) > tol) {
                  *outStream << "Incorrect value for i=" << i << ", j=" << j << ": "
                             << massMatrixHost(0,i,j) << " != " << exactVal << "\n\n";
                  errorFlag++;
                }
              }
          }
        }
      } catch (std::exception &err) {
        *outStream << "UNEXPECTED ERROR !!! ----------------------------------------------------------\n";
        *outStream << err.what() << '\n';
        *outStream << "-------------------------------------------------------------------------------" << "\n\n";
        errorFlag = -1000;
      };

      *outStream 
        << "\n"
        << "===============================================================================\n"
        << "| TEST 4: correctness of basis function derivatives                           |\n"
        << "===============================================================================\n";

      outStream->precision(20);
      
      // function values stored by bf, then pt
      const ValueType basisValues[] = {
        1.000000000000000, 1.000000000000000, 1.000000000000000,	
        1.000000000000000, -1.000000000000000, -0.3333333333333333, 
        0.3333333333333333, 1.000000000000000, 1.000000000000000,	
        -0.3333333333333333, -0.3333333333333333, 1.000000000000000,	
        -1.000000000000000, 0.4074074074074074, -0.4074074074074074,	
        1.000000000000000};

      const ValueType basisD1Values[] = 
        {0, 0, 0, 0, 1.000000000000000, 1.000000000000000, 1.000000000000000, 
         1.000000000000000, -3.000000000000000, -1.000000000000000,		
         1.000000000000000, 3.000000000000000, 6.000000000000000,		
         -0.6666666666666667, -0.6666666666666667, 6.000000000000000};
      
      const ValueType basisD2Values[] = 
        {0, 0, 0, 0, 0, 0, 0, 0, 3.000000000000000, 3.000000000000000,	
         3.000000000000000, 3.000000000000000, -15.00000000000000,		
         -5.000000000000000, 5.000000000000000, 15.00000000000000};
      
      const ValueType basisD3Values[] = 
        {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 15.00000000000000,		
         15.00000000000000, 15.00000000000000, 15.00000000000000};

      try {
        constexpr ordinal_type order = 3;
        if(order <= maxOrder) {
        const double alpha = 0.0, beta = 0.0;

        LineBasisType lineBasis(order, alpha, beta);
        const ordinal_type numFields = lineBasis.getCardinality();
        const ordinal_type spaceDim  = lineBasis.getBaseCellTopology().getDimension();

        DynRankViewHost ConstructWithLabel(lineNodesHost, order+1, spaceDim);
        const ordinal_type numPoints = lineNodesHost.extent(0);

        for (ordinal_type i=0;i<numPoints;++i)
          lineNodesHost(i, 0) = -1.0+(2.0*i)/(numPoints-1);

        const auto lineNodes = Kokkos::create_mirror_view(typename DeviceType::memory_space(), lineNodesHost);
        Kokkos::deep_copy(lineNodes, lineNodesHost);
        
        // test basis values
        {
          *outStream << " -- Comparing OPERATOR_VALUE -- \n\n"; 
          DynRankView ConstructWithLabel(vals, numFields, numPoints);
          lineBasis.getValues(space, vals, lineNodes, OPERATOR_VALUE);

          // host mirror for comparison
          auto valsHost = Kokkos::create_mirror_view(vals);
          Kokkos::deep_copy(valsHost, vals);
          
          for (ordinal_type i=0;i<numFields;++i)
            for (ordinal_type j=0;j<numPoints;++j) {
              // Compute offset for (F,P) container
              const ordinal_type l =  j + i * numPoints;

              const auto val = valsHost(i,j);
              const auto exactVal = basisValues[l];

              if ( std::isnan(val) || std::abs(val - exactVal) > tol ) {
                errorFlag++;
                *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
                
                // Output the multi-index of the value where the error is:
                *outStream << " At multi-index { ";
                *outStream << i << " ";*outStream << j << " ";
                *outStream << "}  computed value: " << val
                           << " but reference value: " << exactVal << "\n";
              }
            }
        }
        
        {
          *outStream << " -- Comparing OPERATOR_D1 -- \n\n"; 
          DynRankView ConstructWithLabel(vals, numFields, numPoints, spaceDim);
          lineBasis.getValues(space, vals, lineNodes, OPERATOR_D1);

          // host mirror for comparison
          auto valsHost = Kokkos::create_mirror_view(vals);
          Kokkos::deep_copy(valsHost, vals);

          for (ordinal_type i=0;i<numFields;++i)
            for (ordinal_type j=0;j<numPoints;++j) {
              // Compute offset for (F,P) container
              const ordinal_type l =  j + i * numPoints;

              const auto val = valsHost(i,j,0);
              const auto exactVal = basisD1Values[l];

              if ( std::isnan(val) || std::abs(val - exactVal) > tol) {
                errorFlag++;
                *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
                
                // Output the multi-index of the value where the error is:
                *outStream << " At multi-index { ";
                *outStream << i << " ";*outStream << j << " ";
                *outStream << "}  computed value: " << val
                           << " but reference value: " << exactVal << "\n";
              }
            }
        }

        {
          *outStream << " -- Comparing OPERATOR_D2 -- \n\n"; 
          DynRankView ConstructWithLabel(vals, numFields, numPoints, spaceDim);
          lineBasis.getValues(space, vals, lineNodes, OPERATOR_D2);

          // host mirror for comparison
          auto valsHost = Kokkos::create_mirror_view(vals);
          Kokkos::deep_copy(valsHost, vals);

          for (ordinal_type i=0;i<numFields;++i)
            for (ordinal_type j=0;j<numPoints;++j) {
              // Compute offset for (F,P) container
              const ordinal_type l =  j + i * numPoints;

              const auto val = valsHost(i,j,0);
              const auto exactVal = basisD2Values[l];

              if ( std::isnan(val) || std::abs(val - exactVal) > tol ) {
                errorFlag++;
                *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
                
                // Output the multi-index of the value where the error is:
                *outStream << " At multi-index { ";
                *outStream << i << " ";*outStream << j << " ";
                *outStream << "}  computed value: " << val
                           << " but reference value: " << exactVal << "\n";
              }
            }
        }

        {
          *outStream << " -- Comparing OPERATOR_D3 -- \n\n"; 
          DynRankView ConstructWithLabel(vals, numFields, numPoints, spaceDim);
          lineBasis.getValues(space, vals, lineNodes, OPERATOR_D3);

          // host mirror for comparison
          auto valsHost = Kokkos::create_mirror_view(vals);
          Kokkos::deep_copy(valsHost, vals);

          for (ordinal_type i=0;i<numFields;++i)
            for (ordinal_type j=0;j<numPoints;++j) {
              // Compute offset for (F,P) container
              const ordinal_type l =  j + i * numPoints;

              const auto val = valsHost(i,j,0);
              const auto exactVal = basisD3Values[l];

              if ( std::isnan(val) || std::abs(val - exactVal) > tol) {
                errorFlag++;
                *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
                
                // Output the multi-index of the value where the error is:
                *outStream << " At multi-index { ";
                *outStream << i << " ";*outStream << j << " ";
                *outStream << "}  computed value: " << val
                           << " but reference value: " << exactVal << "\n";
              }
            }
        }
        }
      } catch (std::exception &err) {
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
        int order = 2;
        const double alpha = 0.0, beta = 0.0;
        LineBasisType lineBasis(order, alpha, beta);
        
        const EFunctionSpace fs = lineBasis.getFunctionSpace();
        
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
  }
}
