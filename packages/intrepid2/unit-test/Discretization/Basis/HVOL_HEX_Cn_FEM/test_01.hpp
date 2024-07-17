// @HEADER
// *****************************************************************************
//                           Intrepid2 Package
//
// Copyright 2007 NTESS and the Intrepid2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/** \file test_01.cpp
    \brief  Unit tests for the Intrepid2::G_HEX_C2_FEM class.
    \author Created by P. Bochev, D. Ridzal, K. Peterson and Kyungjoo Kim.
 */

#include "Intrepid2_config.h"

#ifdef HAVE_INTREPID2_DEBUG
#define INTREPID2_TEST_FOR_DEBUG_ABORT_OVERRIDE_TO_CONTINUE
#endif

#include "Intrepid2_Types.hpp"
#include "Intrepid2_Utils.hpp"

#include "Intrepid2_PointTools.hpp"
#include "Intrepid2_HVOL_HEX_Cn_FEM.hpp"

#include "Teuchos_oblackholestream.hpp"
#include "Teuchos_RCP.hpp"

#include "packages/intrepid2/unit-test/Discretization/Basis/Macros.hpp"

namespace Intrepid2 {

namespace Test {

template<typename OutValueType, typename PointValueType, typename DeviceType>
int HVOL_HEX_Cn_FEM_Test01(const bool verbose) {

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
  << "|                 Unit Test (Basis_HVOL_HEX_C2_FEM)                             |\n"
  << "|                                                                             |\n"
  << "|     1) Conversion of Dof tags into Dof ordinals and back                    |\n"
  << "|     2) Basis values for VALUE, GRAD, and Dk operators                       |\n"
  << "|                                                                             |\n"
  << "|  Questions? Contact  Pavel Bochev  (pbboche@sandia.gov),                    |\n"
  << "|                      Robert Kirby  (robert.c.kirby@ttu.edu),                |\n"
  << "|                      Denis Ridzal  (dridzal@sandia.gov),                    |\n"
  << "|                      Kara Peterson (kjpeter@sandia.gov).                    |\n"
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

  typedef Basis_HVOL_HEX_Cn_FEM<DeviceType,OutValueType,PointValueType> HexBasisType;
  constexpr ordinal_type maxOrder = Parameters::MaxOrder;


  *outStream
  << "\n"
  << "===============================================================================\n"
  << "| TEST 1: Basis creation, exceptions tests                                    |\n"
  << "===============================================================================\n";

  try {

#ifdef HAVE_INTREPID2_DEBUG
    ordinal_type nthrow = 0, ncatch = 0;
    constexpr ordinal_type order = 3;
    if(order < maxOrder) {
      HexBasisType hexBasis(order);

      // Define array containing array of nodes to evaluate
      DynRankViewPointValueType ConstructWithLabelPointView(hexNodes, 27, 3);

      // Generic array for the output values; needs to be properly resized depending on the operator type
      const ordinal_type numFields = hexBasis.getCardinality();
      const ordinal_type numPoints = hexNodes.extent(0);
      //const ordinal_type spaceDim  = hexBasis.getBaseCellTopology().getDimension();

      // exception 1 - 2: CURL and DIV is not supported.
      {
        DynRankViewOutValueType ConstructWithLabelOutView(vals, numFields, numPoints, 3);
        INTREPID2_TEST_ERROR_EXPECTED( hexBasis.getValues(vals, hexNodes, OPERATOR_CURL) );
        INTREPID2_TEST_ERROR_EXPECTED( hexBasis.getValues(vals, hexNodes, OPERATOR_DIV) );
      }

      // Exceptions 3-7: all bf tags/bf Ids below are wrong and should cause getDofOrdinal() and
      // getDofTag() to access invalid array elements thereby causing bounds check exception
      {
        INTREPID2_TEST_ERROR_EXPECTED( hexBasis.getDofOrdinal(3,10,0) );
        INTREPID2_TEST_ERROR_EXPECTED( hexBasis.getDofOrdinal(1,2,3) );
        INTREPID2_TEST_ERROR_EXPECTED( hexBasis.getDofOrdinal(0,4,1) );
        INTREPID2_TEST_ERROR_EXPECTED( hexBasis.getDofTag(numFields) );
        INTREPID2_TEST_ERROR_EXPECTED( hexBasis.getDofTag(-1) );
      }

      // Exceptions 8-18 test exception handling with incorrectly dimensioned input/output arrays
      // exception #8: input points array must be of rank-2
      {
        DynRankViewOutValueType ConstructWithLabelOutView(vals, numFields, numPoints);
        {
          // exception #8: input points array must be of rank-2
          DynRankViewPointValueType ConstructWithLabelPointView(badPoints, 4, 5, 3);
          INTREPID2_TEST_ERROR_EXPECTED( hexBasis.getValues(vals, badPoints, OPERATOR_VALUE) );
        }
        {
          // exception #9: dimension 1 in the input point array must equal space dimension of the cell
          DynRankViewPointValueType ConstructWithLabelPointView(badPoints, 4, 4);
          INTREPID2_TEST_ERROR_EXPECTED( hexBasis.getValues(vals, badPoints, OPERATOR_VALUE) );
        }
        {
          // exception #10: output values must be of rank-2 for OPERATOR_VALUE
          DynRankViewOutValueType ConstructWithLabelOutView(badVals, 4, 3, 1);
          INTREPID2_TEST_ERROR_EXPECTED( hexBasis.getValues(badVals, hexNodes, OPERATOR_VALUE) );
        }
        {
          DynRankViewOutValueType ConstructWithLabelOutView(badVals, 4, 3);

          // exception #11: output values must be of rank-3 for OPERATOR_GRAD
          INTREPID2_TEST_ERROR_EXPECTED( hexBasis.getValues(badVals, hexNodes, OPERATOR_GRAD) );

          // exception #12: output values must be of rank-3 for OPERATOR_CURL
          INTREPID2_TEST_ERROR_EXPECTED( hexBasis.getValues(badVals, hexNodes, OPERATOR_D1) );

          // exception #13: output values must be of rank-3 for OPERATOR_D2
          INTREPID2_TEST_ERROR_EXPECTED( hexBasis.getValues(badVals, hexNodes, OPERATOR_D2) );
        }
      }
      {
        // exception #14: incorrect 0th dimension of output array (must equal number of basis functions)
        DynRankViewOutValueType ConstructWithLabelOutView(badVals, hexBasis.getCardinality() + 1, hexNodes.extent(0));
        INTREPID2_TEST_ERROR_EXPECTED( hexBasis.getValues(badVals, hexNodes, OPERATOR_VALUE) );
      }
      {
        // exception #15: incorrect 1st dimension of output array (must equal number of points)
        DynRankViewOutValueType ConstructWithLabelOutView(badVals, hexBasis.getCardinality(), hexNodes.extent(0) + 1);
        INTREPID2_TEST_ERROR_EXPECTED( hexBasis.getValues(badVals, hexNodes, OPERATOR_VALUE) );
      }
      {
        // exception #16: incorrect 2nd dimension of output array (must equal spatial dimension)
        DynRankViewOutValueType ConstructWithLabelOutView(badVals, hexBasis.getCardinality(), hexNodes.extent(0), 2);
        INTREPID2_TEST_ERROR_EXPECTED( hexBasis.getValues(badVals, hexNodes, OPERATOR_GRAD) );
      }
      {
        DynRankViewOutValueType ConstructWithLabelOutView(badVals, hexBasis.getCardinality(), hexNodes.extent(0), 40);

        // exception #17: incorrect 2nd dimension of output array (must equal spatial dimension)
        INTREPID2_TEST_ERROR_EXPECTED( hexBasis.getValues(badVals, hexNodes, OPERATOR_D2) );

        // exception #18: incorrect 2nd dimension of output array (must equal spatial dimension)
        INTREPID2_TEST_ERROR_EXPECTED( hexBasis.getValues(badVals, hexNodes, OPERATOR_D3) );
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

    const ordinal_type order = std::min(5, maxOrder);
    HexBasisType hexBasis(order);

    const ordinal_type numFields = hexBasis.getCardinality();
    const auto allTags = hexBasis.getAllDofTags();

    // Loop over all tags, lookup the associated dof enumeration and then lookup the tag again
    const ordinal_type dofTagSize = allTags.extent(0);
    for (ordinal_type i=0;i<dofTagSize;++i) {
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
    for(ordinal_type bfOrd=0;bfOrd<numFields;++bfOrd) {
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
  } catch (std::exception &err) {
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
        HexBasisType hexBasis(ip);

        const auto numDofs   = hexBasis.getCardinality();
        const auto numPoints = hexBasis.getCardinality();
        const auto dim  = hexBasis.getBaseCellTopology().getDimension();

        DynRankViewScalarValueType ConstructWithLabel(dofCoords_scalar, numPoints , dim);
        DynRankViewPointValueType ConstructWithLabelPointView(dofCoords, numPoints , dim);
        hexBasis.getDofCoords(dofCoords_scalar);
        RealSpaceTools<DeviceType>::clone(dofCoords,dofCoords_scalar);

        DynRankViewOutValueType ConstructWithLabelOutView(vals, numDofs, numPoints);
        hexBasis.getValues(vals, dofCoords, OPERATOR_VALUE);

        // host mirror for comparison
        auto valsHost = Kokkos::create_mirror_view(vals);
        Kokkos::deep_copy(valsHost, vals);

        for (auto i=0;i<numDofs;++i)
          for (int j=0;j<numPoints;++j) {
            const scalar_type exactVal = (i == j);
            const auto val = valsHost(i,j);
            if (std::isnan(get_scalar_value(val)) || std::abs(val-exactVal) > tol) {
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
      HexBasisType hexBasis(ip);

      const EFunctionSpace fs = hexBasis.getFunctionSpace();

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








