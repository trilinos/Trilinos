// @HEADER
// *****************************************************************************
//                           Intrepid2 Package
//
// Copyright 2007 NTESS and the Intrepid2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/** \file   test_01.cpp
    \brief  Unit tests for the Intrepid2::HDIV_QUAD_I1_FEM class.
    \author Created by P. Bochev, D. Ridzal, K. Peterson and Kyungjoo Kim.
*/
#include "Intrepid2_config.h"

#ifdef HAVE_INTREPID2_DEBUG
#define INTREPID2_TEST_FOR_DEBUG_ABORT_OVERRIDE_TO_CONTINUE
#endif

#include "Intrepid2_HDIV_QUAD_I1_FEM.hpp"

#include "Teuchos_oblackholestream.hpp"
#include "Teuchos_RCP.hpp"

#include "packages/intrepid2/unit-test/Discretization/Basis/Macros.hpp"

namespace Intrepid2 {

  namespace Test {

    template<typename ValueType, typename DeviceType>
    int HDIV_QUAD_I1_FEM_Test01(const bool verbose) {

      Teuchos::RCP<std::ostream> outStream;
      Teuchos::oblackholestream bhs; // outputs nothing

      if (verbose)
        outStream = Teuchos::rcp(&std::cout, false);
      else
        outStream = Teuchos::rcp(&bhs,       false);

      Teuchos::oblackholestream oldFormatState;
      oldFormatState.copyfmt(std::cout);
      using DeviceSpaceType = typename DeviceType::execution_space;
      typedef typename
        Kokkos::DefaultHostExecutionSpace HostSpaceType ;

      *outStream << "DeviceSpace::  "; DeviceSpaceType().print_configuration(*outStream, false);
      *outStream << "HostSpace::    ";   HostSpaceType().print_configuration(*outStream, false);

      *outStream
        << "===============================================================================\n"
        << "|                                                                             |\n"
        << "|                 Unit Test (Basis_HDIV_QUAD_I1_FEM)                          |\n"
        << "|                                                                             |\n"
        << "|     1) Conversion of Dof tags into Dof ordinals and back                    |\n"
        << "|     2) Basis values for VALUE and DIV operators                             |\n"
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

      typedef Kokkos::DynRankView<ValueType,DeviceType> DynRankView;
      typedef Kokkos::DynRankView<ValueType,HostSpaceType> DynRankViewHost;

      const ValueType tol = tolerence();
      int errorFlag = 0;

      typedef ValueType outputValueType;
      typedef ValueType pointValueType;
      Basis_HDIV_QUAD_I1_FEM<DeviceType,outputValueType,pointValueType> quadBasis;

      *outStream
        << "\n"
        << "===============================================================================\n"
        << "| TEST 1: Basis creation, exceptions tests                                    |\n"
        << "===============================================================================\n";

      try{
        ordinal_type nthrow = 0, ncatch = 0;
#ifdef HAVE_INTREPID2_DEBUG

        // Array with the 4 vertices of the reference Quadrilateral, its center and 4 more points
        DynRankView ConstructWithLabel(quadNodes, 9, 2);

        const ordinal_type numFields = quadBasis.getCardinality();
        const ordinal_type numPoints = quadNodes.extent(0);
        const ordinal_type spaceDim = quadBasis.getBaseCellTopology().getDimension();
        DynRankView ConstructWithLabel(vals, numFields, numPoints, spaceDim );

        {
          // exception #1: GRAD cannot be applied to HDIV functions
          INTREPID2_TEST_ERROR_EXPECTED( quadBasis.getValues(vals, quadNodes, OPERATOR_GRAD) );
          
          // exception #2: CURL cannot be applied to HDIV functions
          INTREPID2_TEST_ERROR_EXPECTED( quadBasis.getValues(vals, quadNodes, OPERATOR_CURL) );
        }
        {
          // Exceptions 3-7: all bf tags/bf Ids below are wrong and should cause getDofOrdinal() and
          // getDofTag() to access invalid array elements thereby causing bounds check exception
          // exception #3
          INTREPID2_TEST_ERROR_EXPECTED( quadBasis.getDofOrdinal(3,0,0) );
          // exception #4
          INTREPID2_TEST_ERROR_EXPECTED( quadBasis.getDofOrdinal(1,1,1) );
          // exception #5
          INTREPID2_TEST_ERROR_EXPECTED( quadBasis.getDofOrdinal(0,4,1) );
          // exception #6
          INTREPID2_TEST_ERROR_EXPECTED( quadBasis.getDofTag(numFields) );
          // exception #7
          INTREPID2_TEST_ERROR_EXPECTED( quadBasis.getDofTag(-1) );
        }
    
        // Exceptions 8- test exception handling with incorrectly dimensioned input/output arrays
        {
          // exception #8: input points array must be of rank-2
          DynRankView ConstructWithLabel(badPoints1, 4, 5, 3);
          INTREPID2_TEST_ERROR_EXPECTED( quadBasis.getValues(vals, badPoints1, OPERATOR_VALUE) );
        }
        {
          // exception #9 dimension 1 in the input point array must equal space dimension of the cell
          DynRankView ConstructWithLabel(badPoints2, 4, spaceDim + 1);
          INTREPID2_TEST_ERROR_EXPECTED( quadBasis.getValues(vals, badPoints2, OPERATOR_VALUE) );
        }
        {
          // exception #10 output values must be of rank-3 for OPERATOR_VALUE
          DynRankView ConstructWithLabel(badVals1, 4, 5);
          INTREPID2_TEST_ERROR_EXPECTED( quadBasis.getValues(badVals1, quadNodes, OPERATOR_VALUE) );
        }
        {
          // exception #11 output values must be of rank-2 for OPERATOR_DIV
          DynRankView ConstructWithLabel(badVals2, 4, 5, spaceDim);
          INTREPID2_TEST_ERROR_EXPECTED( quadBasis.getValues(badVals2, quadNodes, OPERATOR_DIV) );
        }
        {
          // exception #12 incorrect 0th dimension of output array (must equal number of basis functions)
          DynRankView ConstructWithLabel(badVals3, quadBasis.getCardinality() + 1, quadNodes.extent(0), spaceDim);
          INTREPID2_TEST_ERROR_EXPECTED( quadBasis.getValues(badVals3, quadNodes, OPERATOR_VALUE) );
        }
        {
          // exception #13 incorrect 0th dimension of output array (must equal number of basis functions)
          DynRankView ConstructWithLabel(badVals4, quadBasis.getCardinality() + 1, quadNodes.extent(0));
          INTREPID2_TEST_ERROR_EXPECTED( quadBasis.getValues(badVals4, quadNodes, OPERATOR_DIV) );
        }
        {
          // exception #14 incorrect 1st dimension of output array (must equal number of points)
          DynRankView ConstructWithLabel(badVals5, quadBasis.getCardinality(), quadNodes.extent(0) + 1, spaceDim);
          INTREPID2_TEST_ERROR_EXPECTED( quadBasis.getValues(badVals5, quadNodes, OPERATOR_VALUE) );
        }
        {
          // exception #15 incorrect 1st dimension of output array (must equal number of points)
          DynRankView ConstructWithLabel(badVals6, quadBasis.getCardinality(), quadNodes.extent(0) + 1);
          INTREPID2_TEST_ERROR_EXPECTED( quadBasis.getValues(badVals6, quadNodes, OPERATOR_DIV) );
        }
        {
          // exception #16: incorrect 2nd dimension of output array (must equal the space dimension)
          DynRankView ConstructWithLabel(badVals7, quadBasis.getCardinality(), quadNodes.extent(0), spaceDim + 1);
          INTREPID2_TEST_ERROR_EXPECTED( quadBasis.getValues(badVals7, quadNodes, OPERATOR_VALUE) );
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
  
      try {
    
        const ordinal_type numFields = quadBasis.getCardinality();
        const auto allTags = quadBasis.getAllDofTags();

        // Loop over all tags, lookup the associated dof enumeration and then lookup the tag again
        const ordinal_type dofTagSize = allTags.extent(0);
        for (ordinal_type i=0;i<dofTagSize;++i) {
          const auto bfOrd = quadBasis.getDofOrdinal(allTags(i,0), allTags(i,1), allTags(i,2));
    
          const auto myTag = quadBasis.getDofTag(bfOrd);
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
        for(ordinal_type bfOrd=0;bfOrd<numFields;++bfOrd) {
          const auto myTag  = quadBasis.getDofTag(bfOrd);
          const auto myBfOrd = quadBasis.getDofOrdinal(myTag(0), myTag(1), myTag(2));
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
  
      try{

        // VALUE: Each row pair gives the 6x3 correct basis set values at an evaluation point: (P,F,D) layout
        double basisValues[] = {
          0, -1.0, 0, 0, 0, 0, -1.0, 0, 0, -1.0, 1.0, 0, 0, \
          0, 0, 0, 0, 0, 1.0, 0, 0, 1.0, 0, 0, 0, 0, 0, 0, 0, \
          1.0, -1.0, 0, 0, -0.5, 0.5, 0, 0, 0.5, \
          -0.5, 0, 0, -0.75, 0.5, 0, 0, 0.25, -0.5, 0, \
          0, -0.25, 0.5, 0, 0, 0.75, -0.5, 0, 0, -0.5, \
          0.25, 0, 0, 0.5, -0.75, 0, 0, -0.5, 0.75, 0, 0, \
          0.5, -0.25, 0  };

        // DIV: each row gives the 6 correct values of the divergence of the 6 basis functions: (P,F) layout
        double basisDivs[] = {
            0.5, 0.5, 0.5, 0.5,
            0.5, 0.5, 0.5, 0.5,
            0.5, 0.5, 0.5, 0.5,
            0.5, 0.5, 0.5, 0.5,
            0.5, 0.5, 0.5, 0.5,
            0.5, 0.5, 0.5, 0.5,
            0.5, 0.5, 0.5, 0.5,
            0.5, 0.5, 0.5, 0.5,
            0.5, 0.5, 0.5, 0.5,
        };
  
        DynRankViewHost ConstructWithLabel(quadNodesHost, 9, 2);

        quadNodesHost(0,0) = -1.0;  quadNodesHost(0,1) = -1.0;
        quadNodesHost(1,0) =  1.0;  quadNodesHost(1,1) = -1.0;
        quadNodesHost(2,0) =  1.0;  quadNodesHost(2,1) =  1.0;
        quadNodesHost(3,0) = -1.0;  quadNodesHost(3,1) =  1.0;

        quadNodesHost(4,0) =  0.0;  quadNodesHost(4,1) =  0.0;
        quadNodesHost(5,0) =  0.0;  quadNodesHost(5,1) = -0.5;
        quadNodesHost(6,0) =  0.0;  quadNodesHost(6,1) =  0.5;
        quadNodesHost(7,0) = -0.5;  quadNodesHost(7,1) =  0.0;
        quadNodesHost(8,0) =  0.5;  quadNodesHost(8,1) =  0.0;

        const auto quadNodes = Kokkos::create_mirror_view(typename DeviceType::memory_space(), quadNodesHost);
        Kokkos::deep_copy(quadNodes, quadNodesHost);

        // Dimensions for the output arrays:
        const ordinal_type numPoints = quadNodes.extent(0);
        const ordinal_type numFields = quadBasis.getCardinality();
        const ordinal_type spaceDim  = quadBasis.getBaseCellTopology().getDimension();

        { 
        DynRankView ConstructWithLabel(vals, numFields, numPoints, spaceDim);
        quadBasis.getValues(vals, quadNodes, OPERATOR_VALUE);
        const auto vals_host = Kokkos::create_mirror_view(typename HostSpaceType::memory_space(), vals);
        Kokkos::deep_copy(vals_host, vals);
        for (ordinal_type i=0;i<numFields;++i)
          for (ordinal_type j=0;j<numPoints;++j)
            for (ordinal_type k=0;k<spaceDim;++k) {
              
              // compute offset for (P,F,D) data layout: indices are P->j, F->i, D->k
              int l = k + i * spaceDim + j * spaceDim * numFields;
              if (std::abs(vals_host(i,j,k) - basisValues[l]) > tol) {
                errorFlag++;
                *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
                
                // Output the multi-index of the value where the error is:
                *outStream << " At multi-index { ";
                *outStream << i << " ";*outStream << j << " ";*outStream << k << " ";
                *outStream << "}  computed value: " << vals_host(i,j,k)
                           << " but reference value: " << basisValues[l] << "\n";
              }
            }
        }

        {
        // Check DIV of basis function: resize vals to rank-2 container
        DynRankView ConstructWithLabel(vals, numFields, numPoints);
        quadBasis.getValues(vals, quadNodes, OPERATOR_DIV);
        const auto vals_host = Kokkos::create_mirror_view(typename HostSpaceType::memory_space(), vals);
        Kokkos::deep_copy(vals_host, vals);
        for (int i = 0; i < numFields; i++) {
          for (int j = 0; j < numPoints; j++) {
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

      // Catch unexpected errors
      } catch (std::logic_error &err) {
        *outStream << err.what() << "\n\n";
        errorFlag = -1000;
      };

      *outStream \
        << "\n"
        << "===============================================================================\n"\
        << "| TEST 4: correctness of DoF locations                                        |\n"\
        << "===============================================================================\n";

      try{
        const ordinal_type numFields = quadBasis.getCardinality();
        const ordinal_type spaceDim  = quadBasis.getBaseCellTopology().getDimension();

        // Check exceptions.
        ordinal_type nthrow = 0, ncatch = 0;
#ifdef HAVE_INTREPID2_DEBUG
        {
          DynRankView ConstructWithLabel(badVals,1,2,3);
          INTREPID2_TEST_ERROR_EXPECTED( quadBasis.getDofCoords(badVals) );
        }
        {
          DynRankView ConstructWithLabel(badVals, 3,2);
          INTREPID2_TEST_ERROR_EXPECTED( quadBasis.getDofCoords(badVals) );
        }
        {
          DynRankView ConstructWithLabel(badVals, 4,3);
          INTREPID2_TEST_ERROR_EXPECTED( quadBasis.getDofCoords(badVals) );
        }
#endif
         // Check if number of thrown exceptions matches the one we expect
         if (nthrow != ncatch) {
           errorFlag++;
           *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
         }
         
        // Check mathematical correctness
        DynRankView ConstructWithLabel(cvals_dev, numFields,spaceDim);
        DynRankView ConstructWithLabel(dofCoeffs_dev, numFields, spaceDim);
        DynRankView ConstructWithLabel(bvals_dev, numFields, numFields, spaceDim); // last dimension is spatial dim

        quadBasis.getDofCoords(cvals_dev);
        quadBasis.getDofCoeffs(dofCoeffs_dev);
        quadBasis.getValues(bvals_dev, cvals_dev, OPERATOR_VALUE);

        const auto bvals = Kokkos::create_mirror_view(typename HostSpaceType::memory_space(), bvals_dev);
        Kokkos::deep_copy(bvals, bvals_dev);
        const auto dofCoeffs = Kokkos::create_mirror_view(typename HostSpaceType::memory_space(), dofCoeffs_dev);
        Kokkos::deep_copy(dofCoeffs, dofCoeffs_dev);
        const auto cvals = Kokkos::create_mirror_view(typename HostSpaceType::memory_space(), cvals_dev);
        Kokkos::deep_copy(cvals, cvals_dev);

        for (ordinal_type i=0;i<numFields;++i)
          for (ordinal_type j=0;j<numFields;++j) {

            ValueType dofValue = 0.0;
            for(ordinal_type d=0;d<spaceDim;++d)
              dofValue += bvals(i,j,d)*dofCoeffs(j,d);

            const ValueType expected_dofValue = (i == j);
            if (std::abs(dofValue - expected_dofValue) > tol) {
              errorFlag++;
              std::stringstream ss;
              ss << "\nValue of basis function " << i << " at (" << cvals(i,0) << ", " << cvals(i,1) << ") is " << dofValue << " but should be " << expected_dofValue << "\n";
              *outStream << ss.str();
            }
          }
      } catch (std::logic_error &err) { 
        *outStream << err.what() << "\n\n";
        errorFlag = -1000;
      };
      
      *outStream
      << "\n"
      << "===============================================================================\n"
      << "| TEST 4: Function Space is Correct                                           |\n"
      << "===============================================================================\n";
      
      try {
        const EFunctionSpace fs = quadBasis.getFunctionSpace();
        
        if (fs != FUNCTION_SPACE_HDIV)
        {
          *outStream << std::setw(70) << "------------- TEST FAILURE! -------------" << "\n";
          
          // Output the multi-index of the value where the error is:
          *outStream << " Expected a function space of FUNCTION_SPACE_HDIV (enum value " << FUNCTION_SPACE_HDIV << "),";
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
