// @HEADER
// *****************************************************************************
//                           Intrepid2 Package
//
// Copyright 2007 NTESS and the Intrepid2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/** \file test_01.hpp
    \brief  Unit tests for the Intrepid2::HVOL_QUAD_Cn_FEM class.
    \author Created by P. Bochev, D. Ridzal, K. Peterson, Kyungjoo Kim
*/


#include "Intrepid2_config.h"

#ifdef HAVE_INTREPID2_DEBUG
#define INTREPID2_TEST_FOR_DEBUG_ABORT_OVERRIDE_TO_CONTINUE
#endif

#include "Intrepid2_Types.hpp"
#include "Intrepid2_Utils.hpp"

#include "Intrepid2_PointTools.hpp"
#include "Intrepid2_HVOL_QUAD_Cn_FEM.hpp"

#include "Teuchos_oblackholestream.hpp"
#include "Teuchos_RCP.hpp"

#include "packages/intrepid2/unit-test/Discretization/Basis/Macros.hpp"

namespace Intrepid2 {

  namespace Test {

    template<typename OutValueType, typename PointValueType, typename DeviceType>
    int HVOL_QUAD_Cn_FEM_Test01(const bool verbose) {

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
        << "|                 Unit Test (Basis_HVOL_QUAD_Cn_FEM)                            |\n"
        << "|                                                                             |\n"
        << "|     1) Conversion of Dof tags into Dof ordinals and back                    |\n"
        << "|     2) Basis values for VALUE, GRAD, CURL, and Dk operators                 |\n"
        << "|                                                                             |\n"
        << "|  Questions? Contact  Pavel Bochev  (pbboche@sandia.gov),                    |\n"
        << "|                      Robert Kirby  (robert.c.kirby@ttu.edu),                |\n"
        << "|                      Denis Ridzal  (dridzal@sandia.gov),                    |\n"
        << "|                      Kara Peterson (kjpeter@sandia.gov),                    |\n"
        << "|                      Kyungjoo Kim  (kyukim@sandia.gov).                     |\n"
        << "|                                                                             |\n"
        << "|  Intrepid's website: http://trilinos.sandia.gov/packages/intrepid           |\n"
        << "|  Trilinos website:   http://trilinos.sandia.gov                             |\n"
        << "|                                                                             |\n"
        << "===============================================================================\n";

      typedef Kokkos::DynRankView<PointValueType,DeviceType> DynRankViewPointValueType;
      typedef Kokkos::DynRankView<OutValueType,DeviceType> DynRankViewOutValueType;
      typedef typename ScalarTraits<OutValueType>::scalar_type scalar_type;
      typedef Kokkos::DynRankView<scalar_type, DeviceType> DynRankViewScalarValueType;

      const scalar_type tol = tolerence();
      int errorFlag = 0;

      typedef Basis_HVOL_QUAD_Cn_FEM<DeviceType,OutValueType,PointValueType> QuadBasisType;
      constexpr ordinal_type maxOrder = Parameters::MaxOrder;


      *outStream
        << "\n"
        << "===============================================================================\n"
        << "| TEST 1: Basis creation, exceptions tests                                    |\n"
        << "===============================================================================\n";

      try{
#ifdef HAVE_INTREPID2_DEBUG
        ordinal_type nthrow = 0, ncatch = 0;
        constexpr  ordinal_type order = 5;
        if (order <= maxOrder) {
          QuadBasisType quadBasis(order);

          // Define array containing array of nodes to evaluate
          DynRankViewPointValueType ConstructWithLabelPointView(quadNodes, 10, 2);

          // Generic array for the output values; needs to be properly resized depending on the operator type
          const auto numFields = quadBasis.getCardinality();
          const auto numPoints = quadNodes.extent(0);
          //const auto spaceDim  = quadBasis.getBaseCellTopology().getDimension();

          // Exceptions 1-5: all bf tags/bf Ids below are wrong and should cause getDofOrdinal() and
          // getDofTag() to access invalid array elements thereby causing bounds check exception
          {
            INTREPID2_TEST_ERROR_EXPECTED( quadBasis.getDofOrdinal(3,0,0) );
            INTREPID2_TEST_ERROR_EXPECTED( quadBasis.getDofOrdinal(1,0,4) );
            INTREPID2_TEST_ERROR_EXPECTED( quadBasis.getDofOrdinal(0,4,0) );
            INTREPID2_TEST_ERROR_EXPECTED( quadBasis.getDofTag(numFields) );
            INTREPID2_TEST_ERROR_EXPECTED( quadBasis.getDofTag(-1) );
          }

          // Exceptions 6-16 test exception handling with incorrectly dimensioned input/output arrays
          {
            DynRankViewOutValueType ConstructWithLabelOutView(vals, numFields, numPoints);
            {
              // exception #6: input points array must be of rank-2
              DynRankViewPointValueType ConstructWithLabelPointView(badPoints, 4, 5, 3);
              INTREPID2_TEST_ERROR_EXPECTED( quadBasis.getValues(vals, badPoints, OPERATOR_VALUE) );
            }
            {
              // exception #7: dimension 1 in the input point array must equal space dimension of the cell
              DynRankViewPointValueType ConstructWithLabelPointView(badPoints, 4, 3);
              INTREPID2_TEST_ERROR_EXPECTED( quadBasis.getValues(vals, badPoints, OPERATOR_VALUE) );
            }
            {
              // exception #8: output values must be of rank-2 for OPERATOR_VALUE
              DynRankViewOutValueType ConstructWithLabelOutView(badVals, 4, 3, 1);
              INTREPID2_TEST_ERROR_EXPECTED( quadBasis.getValues(badVals, quadNodes, OPERATOR_VALUE) );
            }
            {
              DynRankViewOutValueType ConstructWithLabelOutView(badVals, 4, 3);

              // exception #9: output values must be of rank-3 for OPERATOR_GRAD
              INTREPID2_TEST_ERROR_EXPECTED( quadBasis.getValues(badVals, quadNodes, OPERATOR_GRAD) );

              // exception #10: output values must be of rank-3 for OPERATOR_CURL
              INTREPID2_TEST_ERROR_EXPECTED( quadBasis.getValues(badVals, quadNodes, OPERATOR_CURL) );

              // exception #11: output values must be of rank-3 for OPERATOR_DIV
              INTREPID2_TEST_ERROR_EXPECTED( quadBasis.getValues(badVals, quadNodes, OPERATOR_CURL) );

              // exception #12: output values must be of rank-3 for OPERATOR_D2
              INTREPID2_TEST_ERROR_EXPECTED( quadBasis.getValues(badVals, quadNodes, OPERATOR_D2) );
            }
          }
          {
            // exception #13: incorrect 0th dimension of output array (must equal number of basis functions)
            DynRankViewOutValueType ConstructWithLabelOutView(badVals, quadBasis.getCardinality() + 1, quadNodes.extent(0));
            INTREPID2_TEST_ERROR_EXPECTED( quadBasis.getValues(badVals, quadNodes, OPERATOR_VALUE) );
          }
          {
            // exception #14: incorrect 1st dimension of output array (must equal number of points)
            DynRankViewOutValueType ConstructWithLabelOutView(badVals, quadBasis.getCardinality(), quadNodes.extent(0) + 1);
            INTREPID2_TEST_ERROR_EXPECTED( quadBasis.getValues(badVals, quadNodes, OPERATOR_VALUE) );
          }
          {
            // exception #15: incorrect 2nd dimension of output array (must equal spatial dimension)
            DynRankViewOutValueType ConstructWithLabelOutView(badVals, quadBasis.getCardinality(), quadNodes.extent(0), 3);
            INTREPID2_TEST_ERROR_EXPECTED( quadBasis.getValues(badVals, quadNodes, OPERATOR_GRAD) );
          }
          {
            DynRankViewOutValueType ConstructWithLabelOutView(badVals, quadBasis.getCardinality(), quadNodes.extent(0), 40);

            // exception #16: incorrect 2nd dimension of output array (must equal spatial dimension)
            INTREPID2_TEST_ERROR_EXPECTED( quadBasis.getValues(badVals, quadNodes, OPERATOR_D2) );

            // exception #17: incorrect 2nd dimension of output array (must equal spatial dimension)
            INTREPID2_TEST_ERROR_EXPECTED( quadBasis.getValues(badVals, quadNodes, OPERATOR_D3) );
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
        << "| TEST 2: correctness of tag to enum and enum to tag lookups                  |\n"
        << "===============================================================================\n";

      try {
        const auto order = std::min(5, maxOrder);
        QuadBasisType quadBasis(order);

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

      *outStream
        << "\n"
        << "===============================================================================\n"
        << "| TEST 3: correctness of basis function values                                |\n"
        << "===============================================================================\n";
      outStream->precision(20);

      // check Kronecker property for Lagrange polynomials.
      try {
        const EPointType pts[3] = { POINTTYPE_EQUISPACED, POINTTYPE_WARPBLEND, POINTTYPE_GAUSS };
        for (auto idx=0;idx<3;++idx) {
          *outStream << " -- Testing " << EPointTypeToString(pts[idx]) << " -- \n";
          for (auto ip=0;ip<std::min(5, maxOrder);++ip) {
            QuadBasisType quadBasis(ip);

            const auto numDofs   = quadBasis.getCardinality();
            const auto numPoints = quadBasis.getCardinality();
            const auto spaceDim  = quadBasis.getBaseCellTopology().getDimension();

            DynRankViewScalarValueType ConstructWithLabel(dofCoords_scalar, numPoints, spaceDim);
            DynRankViewPointValueType ConstructWithLabelPointView(dofCoords, numDofs , spaceDim);

            quadBasis.getDofCoords(dofCoords_scalar);
            RealSpaceTools<DeviceType>::clone(dofCoords,dofCoords_scalar);

            DynRankViewOutValueType ConstructWithLabelOutView(vals, numDofs, numPoints);
            quadBasis.getValues(vals, dofCoords, OPERATOR_VALUE);

            // host mirror for comparison
            auto valsHost = Kokkos::create_mirror_view(vals);
            Kokkos::deep_copy(valsHost, vals);

            for (auto i=0;i<numDofs;++i)
              for (int j=0;j<numPoints;++j) {
                const scalar_type exactVal = (i == j);
                const auto val = get_scalar_value(valsHost(i,j));
                if (std::isnan(val) || std::abs(val-exactVal) > tol) {
                  errorFlag++;
                  *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
                  *outStream << " Basis function at i= " << i << ", j=" << j << ": "
                             << val << " != " << exactVal << "\n\n";
                }
              }
          }
        }
      } catch (std::exception &err) {
        *outStream << "UNEXPECTED ERROR !!! ----------------------------------------------------------\n";
        *outStream << err.what() << "\n\n";
        errorFlag = -1000;
      };

      *outStream
      << "\n"
      << "===============================================================================\n"
      << "| TEST 4: Function Space is Correct                                           |\n"
      << "===============================================================================\n";

      try {
        for (auto ip=0;ip<std::min(5, maxOrder);++ip) {
          QuadBasisType quadBasis(ip);

          const EFunctionSpace fs = quadBasis.getFunctionSpace();

          if (fs != FUNCTION_SPACE_HVOL)
          {
            *outStream << std::setw(70) << "------------- TEST FAILURE! -------------" << "\n";

            // Output the multi-index of the value where the error is:
            *outStream << " Expected a function space of FUNCTION_SPACE_HVOL (enum value " << FUNCTION_SPACE_HVOL << "),";
            *outStream << " but got " << fs << "\n";
            if (fs == FUNCTION_SPACE_MAX)
            {
              *outStream << "Note that this matches the default value defined by superclass, FUNCTION_SPACE_MAX.  Likely the subclass has failed to set the superclass functionSpace_ field.\n";
            }
            errorFlag++;
          }
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
