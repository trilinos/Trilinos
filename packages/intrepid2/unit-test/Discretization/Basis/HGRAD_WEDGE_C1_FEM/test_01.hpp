// @HEADER
// *****************************************************************************
//                           Intrepid2 Package
//
// Copyright 2007 NTESS and the Intrepid2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/** \file test_01.hpp
    \brief  Unit tests for the Intrepid2::HGRAD_WEDGE_C1_FEM class.
    \author Created by P. Bochev, D. Ridzal, K. Peterson and Kyungjoo Kim.
*/

#include "Intrepid2_config.h"

#ifdef HAVE_INTREPID2_DEBUG
#define INTREPID2_TEST_FOR_DEBUG_ABORT_OVERRIDE_TO_CONTINUE
#endif

#include "Intrepid2_HGRAD_WEDGE_C1_FEM.hpp"

#include "Teuchos_oblackholestream.hpp"
#include "Teuchos_RCP.hpp"

#include "packages/intrepid2/unit-test/Discretization/Basis/Macros.hpp"

namespace Intrepid2 {

  namespace Test {

    template<typename ValueType, typename DeviceType>
    int HGRAD_WEDGE_C1_FEM_Test01(const bool verbose) {

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
        << "|                 Unit Test (Basis_HGRAD_WEDGE_C1_FEM)                        |\n"
        << "|                                                                             |\n"
        << "|     1) Conversion of Dof tags into Dof ordinals and back                    |\n"
        << "|     2) Basis values for VALUE, GRAD, CURL, and Dk operators                 |\n"
        << "|                                                                             |\n"
        << "|  Questions? Contact  Pavel Bochev  (pbboche@sandia.gov),                    |\n"
        << "|                      Denis Ridzal  (dridzal@sandia.gov),                    |\n"
        << "|                      Kara Peterson (kjpeter@sandia.gov).                    |\n"
        << "|                                                                             |\n"
        << "|  Intrepid's website: http://trilinos.sandia.gov/packages/intrepid           |\n"
        << "|  Trilinos website:   http://trilinos.sandia.gov                             |\n"
        << "|                                                                             |\n"
        << "===============================================================================\n";

        typedef Kokkos::DynRankView<ValueType,DeviceType> DynRankView;
        typedef Kokkos::DynRankView<ValueType,HostSpaceType>   DynRankViewHost;

      const ValueType tol = tolerence();
      int errorFlag = 0;

      // for virtual function, value and point types are declared in the class
      typedef ValueType outputValueType;
      typedef ValueType pointValueType;
      Basis_HGRAD_WEDGE_C1_FEM<DeviceType,outputValueType,pointValueType> wedgeBasis;

      *outStream
        << "\n"
        << "===============================================================================\n"
        << "| TEST 1: constructors and exceptions                                         |\n"
        << "===============================================================================\n";


      try {
        ordinal_type nthrow = 0, ncatch = 0;
#ifdef HAVE_INTREPID2_DEBUG

        // Define array containing the 4 vertices of the reference WEDGE and its center.
        DynRankView ConstructWithLabel(wedgeNodes, 12, 3);

        // Generic array for the output values; needs to be properly resized depending on the operator type
        const ordinal_type numFields = wedgeBasis.getCardinality();
        const ordinal_type numPoints = wedgeNodes.extent(0);
        const ordinal_type spaceDim  = wedgeBasis.getBaseCellTopology().getDimension();

        DynRankView vals("vals", numFields, numPoints);
        DynRankView vals_vec("vals", numFields, numPoints, spaceDim);
        {
          // exception #1: CURL cannot be applied to scalar functions
          INTREPID2_TEST_ERROR_EXPECTED( wedgeBasis.getValues(vals_vec, wedgeNodes, OPERATOR_CURL) );

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
          INTREPID2_TEST_ERROR_EXPECTED( wedgeBasis.getDofOrdinal(0,6,0) );
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
          DynRankView ConstructWithLabel( badPoints2, 4, 2);
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
          DynRankView ConstructWithLabel( badVals5, numFields, numPoints, 4);
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
      } catch (std::exception &err) {
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
      } catch (std::logic_error &err) {
        *outStream << err.what() << "\n\n";
        errorFlag = -1000;
      }

      *outStream
        << "\n"
        << "===============================================================================\n"
        << "| TEST 3: correctness of basis function values                                |\n"
        << "===============================================================================\n";
  
      outStream -> precision(20);
  
      // VALUE: Each row gives the 4 correct basis set values at an evaluation point
      const ValueType basisValues[] = {
        1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 
        0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 
        0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 
        0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 
        0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 
        0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 
        //
        0.25,     0.25,     0.5,    0.,     0.,     0.,
        0.125,    0.25,     0.125,  0.125,  0.25,   0.125,
        0.,       0.,       0.,     0.45,   0.25,   0.3,
        0.0875,   0.0375,   0,      0.6125, 0.2625, 0,
        0.3444,   0,        0.2706, 0.2156, 0,      0.1694,
        0.,       0.2,      0.3,    0.,     0.2,    0.3
      };
  
      // GRAD and D1: each row gives the 3 x 4 correct values of the gradients of the 4 basis functions
      const ValueType basisGrads[] = {
        -1., -1., -0.5, 1., 0, 0, 0, 1., 0, 0., 0., 0.5, 0., 0, 0, 0, 0., 0, \
        -1., -1., 0., 1., 0, -0.5, 0, 1., 0, 0., 0., 0., 0., 0, 0.5, 0, 0., \
        0, -1., -1., 0., 1., 0, 0, 0, 1., -0.5, 0., 0., 0., 0., 0, 0, 0, 0., \
        0.5, 0., 0., -0.5, 0., 0, 0, 0, 0., 0, -1., -1., 0.5, 1., 0, 0, 0, \
        1., 0, 0., 0., 0., 0., 0, -0.5, 0, 0., 0, -1., -1., 0., 1., 0, 0.5, \
        0, 1., 0, 0., 0., 0., 0., 0, 0, 0, 0., -0.5, -1., -1., 0., 1., 0, 0, \
        0, 1., 0.5, -1., -1., -0.125, 1., 0, -0.125, 0, 1., -0.25, 0., 0., \
        0.125, 0., 0, 0.125, 0, 0., 0.25, -0.5, -0.5, -0.125, 0.5, 0, -0.25, \
        0, 0.5, -0.125, -0.5, -0.5, 0.125, 0.5, 0, 0.25, 0, 0.5, 0.125, 0., \
        0., -0.225, 0., 0, -0.125, 0, 0., -0.15, -1., -1., 0.225, 1., 0, \
        0.125, 0, 1., 0.15, -0.125, -0.125, -0.35, 0.125, 0, -0.15, 0, 0.125, \
        0, -0.875, -0.875, 0.35, 0.875, 0, 0.15, 0, 0.875, 0, -0.615, -0.615, \
        -0.28, 0.615, 0, 0, 0, 0.615, -0.22, -0.385, -0.385, 0.28, 0.385, 0, \
        0, 0, 0.385, 0.22, -0.5, -0.5, 0., 0.5, 0, -0.2, 0, 0.5, -0.3, -0.5, \
        -0.5, 0., 0.5, 0, 0.2, 0, 0.5, 0.3
      };

      //D2: flat array with the values of D2 applied to basis functions. Multi-index is (P,F,K)
      const ValueType basisD2[] = {
        0, 0, 0.5, 0, 0.5, 0, 0, 0, -0.5, 0, 0, 0, 0, 0, 0, 0, -0.5, 0, 0, 0, \
        -0.5, 0, -0.5, 0, 0, 0, 0.5, 0, 0, 0, 0, 0, 0, 0, 0.5, 0, 0, 0, 0.5, \
        0, 0.5, 0, 0, 0, -0.5, 0, 0, 0, 0, 0, 0, 0, -0.5, 0, 0, 0, -0.5, 0, \
        -0.5, 0, 0, 0, 0.5, 0, 0, 0, 0, 0, 0, 0, 0.5, 0, 0, 0, 0.5, 0, 0.5, \
        0, 0, 0, -0.5, 0, 0, 0, 0, 0, 0, 0, -0.5, 0, 0, 0, -0.5, 0, -0.5, 0, \
        0, 0, 0.5, 0, 0, 0, 0, 0, 0, 0, 0.5, 0, 0, 0, 0.5, 0, 0.5, 0, 0, 0, \
        -0.5, 0, 0, 0, 0, 0, 0, 0, -0.5, 0, 0, 0, -0.5, 0, -0.5, 0, 0, 0, \
        0.5, 0, 0, 0, 0, 0, 0, 0, 0.5, 0, 0, 0, 0.5, 0, 0.5, 0, 0, 0, -0.5, \
        0, 0, 0, 0, 0, 0, 0, -0.5, 0, 0, 0, -0.5, 0, -0.5, 0, 0, 0, 0.5, 0, \
        0, 0, 0, 0, 0, 0, 0.5, 0, 0, 0, 0.5, 0, 0.5, 0, 0, 0, -0.5, 0, 0, 0, \
        0, 0, 0, 0, -0.5, 0, 0, 0, -0.5, 0, -0.5, 0, 0, 0, 0.5, 0, 0, 0, 0, \
        0, 0, 0, 0.5, 0, 0, 0, 0.5, 0, 0.5, 0, 0, 0, -0.5, 0, 0, 0, 0, 0, 0, \
        0, -0.5, 0, 0, 0, -0.5, 0, -0.5, 0, 0, 0, 0.5, 0, 0, 0, 0, 0, 0, 0, \
        0.5, 0, 0, 0, 0.5, 0, 0.5, 0, 0, 0, -0.5, 0, 0, 0, 0, 0, 0, 0, -0.5, \
        0, 0, 0, -0.5, 0, -0.5, 0, 0, 0, 0.5, 0, 0, 0, 0, 0, 0, 0, 0.5, 0, 0, \
        0, 0.5, 0, 0.5, 0, 0, 0, -0.5, 0, 0, 0, 0, 0, 0, 0, -0.5, 0, 0, 0, \
        -0.5, 0, -0.5, 0, 0, 0, 0.5, 0, 0, 0, 0, 0, 0, 0, 0.5, 0, 0, 0, 0.5, \
        0, 0.5, 0, 0, 0, -0.5, 0, 0, 0, 0, 0, 0, 0, -0.5, 0, 0, 0, -0.5, 0, \
        -0.5, 0, 0, 0, 0.5, 0, 0, 0, 0, 0, 0, 0, 0.5, 0, 0, 0, 0.5, 0, 0.5, \
        0, 0, 0, -0.5, 0, 0, 0, 0, 0, 0, 0, -0.5, 0, 0, 0, -0.5, 0, -0.5, 0, \
        0, 0, 0.5, 0, 0, 0, 0, 0, 0, 0, 0.5, 0, 0, 0, 0.5, 0, 0.5, 0, 0, 0, \
        -0.5, 0, 0, 0, 0, 0, 0, 0, -0.5, 0, 0, 0, -0.5, 0, -0.5, 0, 0, 0, \
        0.5, 0, 0, 0, 0, 0, 0, 0, 0.5, 0
      };
      
      try {
        DynRankViewHost ConstructWithLabel(wedgeNodesHost, 12, 3);

        wedgeNodesHost(0,0) =  0.0;  wedgeNodesHost(0,1) =  0.0;  wedgeNodesHost(0,2) = -1.0;  
        wedgeNodesHost(1,0) =  1.0;  wedgeNodesHost(1,1) =  0.0;  wedgeNodesHost(1,2) = -1.0;  
        wedgeNodesHost(2,0) =  0.0;  wedgeNodesHost(2,1) =  1.0;  wedgeNodesHost(2,2) = -1.0;
        wedgeNodesHost(3,0) =  0.0;  wedgeNodesHost(3,1) =  0.0;  wedgeNodesHost(3,2) =  1.0;  
        wedgeNodesHost(4,0) =  1.0;  wedgeNodesHost(4,1) =  0.0;  wedgeNodesHost(4,2) =  1.0;  
        wedgeNodesHost(5,0) =  0.0;  wedgeNodesHost(5,1) =  1.0;  wedgeNodesHost(5,2) =  1.0;
        
        wedgeNodesHost(6,0) =  0.25; wedgeNodesHost(6,1) =  0.5;  wedgeNodesHost(6,2) = -1.0;  
        wedgeNodesHost(7,0) =  0.5;  wedgeNodesHost(7,1) =  0.25; wedgeNodesHost(7,2) =  0.0;  
        wedgeNodesHost(8,0) =  0.25; wedgeNodesHost(8,1) =  0.30; wedgeNodesHost(8,2) =  1.0;
        wedgeNodesHost(9,0) =  0.3;  wedgeNodesHost(9,1) =  0.0;  wedgeNodesHost(9,2) =  0.75;
        wedgeNodesHost(10,0)=  0.0;  wedgeNodesHost(10,1)=  0.44; wedgeNodesHost(10,2)= -0.23;  
        wedgeNodesHost(11,0)=  0.4;  wedgeNodesHost(11,1)=  0.6;  wedgeNodesHost(11,2)=  0.0;  

        auto wedgeNodes = Kokkos::create_mirror_view(typename DeviceType::memory_space(), wedgeNodesHost);
        Kokkos::deep_copy(wedgeNodes, wedgeNodesHost);

        // Dimensions for the output arrays:
        const ordinal_type numFields = wedgeBasis.getCardinality();
        const ordinal_type numPoints = wedgeNodes.extent(0);
        const ordinal_type spaceDim  = wedgeBasis.getBaseCellTopology().getDimension();
        const ordinal_type D2Cardin  = getDkCardinality(OPERATOR_D2, spaceDim);
        
        // Check VALUE of basis functions: resize vals to rank-2 container:
        {
          DynRankView vals = DynRankView("vals", numFields, numPoints);
          wedgeBasis.getValues(vals, wedgeNodes, OPERATOR_VALUE);
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
          wedgeBasis.getValues(vals, wedgeNodes, OPERATOR_GRAD);
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

        // Check D1 of basis function (do not resize vals because it has the correct size: D1 = GRAD)
        {
          DynRankView vals = DynRankView("vals", numFields, numPoints, spaceDim);
          wedgeBasis.getValues(vals, wedgeNodes, OPERATOR_D1);
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
          wedgeBasis.getValues(vals, wedgeNodes, OPERATOR_D2);
          auto vals_host = Kokkos::create_mirror_view(typename HostSpaceType::memory_space(), vals);
          Kokkos::deep_copy(vals_host, vals);
          for (ordinal_type i=0;i<numFields;++i) {
            for (ordinal_type j=0;j<numPoints;++j) {
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
                    *outStream << "}  computed D"<< ord <<" component: " << vals(i1,i2,i3)
                               << " but reference D" << ord << " component:  0 \n";
                  }
                }
          }
        }
      } catch (std::exception &err) {
        *outStream << err.what() << "\n\n";
        errorFlag = -1000;
      }

      *outStream
        << "\n"
        << "===============================================================================\n"
        << "| TEST 4: correctness of DoF locations                                        |\n"
        << "===============================================================================\n";

      try {
        const ordinal_type numFields = wedgeBasis.getCardinality();
        const ordinal_type spaceDim  = wedgeBasis.getBaseCellTopology().getDimension();

        // Check exceptions.
        ordinal_type nthrow = 0, ncatch = 0;
#ifdef HAVE_INTREPID2_DEBUG
        {
          DynRankView ConstructWithLabel(badVals, 1,2,3);
          INTREPID2_TEST_ERROR_EXPECTED( wedgeBasis.getDofCoords(badVals) );
        }
        {
          DynRankView ConstructWithLabel(badVals, 4,2);
          INTREPID2_TEST_ERROR_EXPECTED( wedgeBasis.getDofCoords(badVals) );
        }
        {
          DynRankView ConstructWithLabel(badVals, 3,3);
          INTREPID2_TEST_ERROR_EXPECTED( wedgeBasis.getDofCoords(badVals) );
        }
#endif
        if (nthrow != ncatch) {
          errorFlag++;
          *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
          *outStream << "# of catch ("<< ncatch << ") is different from # of throw (" << ncatch << ")\n";
        }

        DynRankView ConstructWithLabel(bvals, numFields, numFields);
        DynRankView ConstructWithLabel(cvals, numFields, spaceDim);

        // Check mathematical correctness.
        wedgeBasis.getDofCoords(cvals);
        wedgeBasis.getValues(bvals, cvals, OPERATOR_VALUE);

        auto cvals_host = Kokkos::create_mirror_view(typename HostSpaceType::memory_space(), cvals);
        Kokkos::deep_copy(cvals_host, cvals);

        auto bvals_host = Kokkos::create_mirror_view(typename HostSpaceType::memory_space(), bvals);

        Kokkos::deep_copy(bvals_host, bvals);
        for (ordinal_type i=0;i<numFields;++i) {
          for (ordinal_type j=0;j<numFields;++j) {
            const ValueType expected_value = (i == j);
            const ValueType value = bvals_host(i,j);
            if (std::abs(value - expected_value) > tol) {
              errorFlag++;
              std::stringstream ss;
              ss << "\nValue of basis function " << i << " at (" << cvals_host(i,0) << ", " << cvals_host(i,1)<< ") is " << value << " but should be " << expected_value << "\n";
              *outStream << ss.str();
            }
          }
        }
      } catch (std::logic_error &err) {
        *outStream << err.what() << "\n\n";
        errorFlag = -1000;
      }
      
      *outStream
      << "\n"
      << "===============================================================================\n"
      << "| TEST 5: Function Space is Correct                                           |\n"
      << "===============================================================================\n";
      
      try {
        const EFunctionSpace fs = wedgeBasis.getFunctionSpace();
        
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

