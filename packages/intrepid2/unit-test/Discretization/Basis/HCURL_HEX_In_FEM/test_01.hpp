// @HEADER
// *****************************************************************************
//                           Intrepid2 Package
//
// Copyright 2007 NTESS and the Intrepid2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/** \file test_01.hpp
    \brief  Unit tests for the Intrepid2::HCURL_HEX_In_FEM class.
    \author Created by P. Bochev, D. Ridzal, K. Peterson, Kyungjoo Kim and Mauro Perego
 */


#include "Intrepid2_config.h"

#ifdef HAVE_INTREPID2_DEBUG
#define INTREPID2_TEST_FOR_DEBUG_ABORT_OVERRIDE_TO_CONTINUE
#endif

#include "Intrepid2_Types.hpp"
#include "Intrepid2_Utils.hpp"

#include "Intrepid2_PointTools.hpp"
#include "Intrepid2_HCURL_HEX_In_FEM.hpp"


#include "Teuchos_oblackholestream.hpp"
#include "Teuchos_RCP.hpp"

#include "packages/intrepid2/unit-test/Discretization/Basis/Macros.hpp"

namespace Intrepid2 {

namespace Test {

template<typename OutValueType, typename PointValueType, typename DeviceType>
int HCURL_HEX_In_FEM_Test01(const bool verbose) {
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
  << "|                      Unit Test HCURL_HEX_In_FEM                             |\n"
  << "|                                                                             |\n"
  << "|  Questions? Contact  Pavel Bochev  (pbboche@sandia.gov),                    |\n"
  << "|                      Robert Kirby  (robert.c.kirby@ttu.edu),                |\n"
  << "|                      Denis Ridzal  (dridzal@sandia.gov),                    |\n"
  << "|                      Kara Peterson (kjpeter@sandia.gov),                    |\n"
  << "|                      Kyungjoo Kim  (kyukim@sandia.gov),                     |\n"
  << "|                      Mauro Perego  (mperego@sandia.gov).                    |\n"
  << "|                                                                             |\n"
  << "===============================================================================\n";

  typedef Kokkos::DynRankView<PointValueType,DeviceType> DynRankViewPointValueType;
  typedef Kokkos::DynRankView<OutValueType,DeviceType> DynRankViewOutValueType;
  typedef typename ScalarTraits<OutValueType>::scalar_type scalar_type;
  typedef Kokkos::DynRankView<scalar_type, DeviceType> DynRankViewScalarValueType;
  typedef Kokkos::DynRankView<scalar_type, HostSpaceType> DynRankViewHostScalarValueType;

  const scalar_type tol = tolerence();
  int errorFlag = 0;

  typedef Basis_HCURL_HEX_In_FEM<DeviceType,OutValueType,PointValueType> HexBasisType;
  constexpr ordinal_type maxOrder = Parameters::MaxOrder ;
  constexpr ordinal_type dim = 3;

  *outStream
  << "\n"
  << "===============================================================================\n"
  << "| TEST 1: Exception testing                                                   |\n"
  << "===============================================================================\n";

  try {

#ifdef HAVE_INTREPID2_DEBUG

    ordinal_type nthrow = 0, ncatch = 0;

    const ordinal_type order = 1;
    HexBasisType hexBasis(order);

    // Array of reference hex nodes - used for evaluation of basis
    DynRankViewHostScalarValueType ConstructWithLabel(hexNodesHost, 8, 3);

    hexNodesHost(0,0) = -1.0; hexNodesHost(0,1) = -1.0; hexNodesHost(0,2) = -1.0;
    hexNodesHost(1,0) =  1.0; hexNodesHost(1,1) = -1.0; hexNodesHost(1,2) = -1.0;
    hexNodesHost(2,0) = -1.0; hexNodesHost(2,1) =  1.0; hexNodesHost(2,2) = -1.0;
    hexNodesHost(3,0) =  1.0; hexNodesHost(3,1) =  1.0; hexNodesHost(3,2) = -1.0;
    hexNodesHost(4,0) = -1.0; hexNodesHost(4,1) = -1.0; hexNodesHost(4,2) =  1.0;
    hexNodesHost(5,0) =  1.0; hexNodesHost(5,1) = -1.0; hexNodesHost(5,2) =  1.0;
    hexNodesHost(6,0) = -1.0; hexNodesHost(6,1) =  1.0; hexNodesHost(6,2) =  1.0;
    hexNodesHost(7,0) =  1.0; hexNodesHost(7,1) =  1.0; hexNodesHost(7,2) =  1.0;

    auto hexNodes_scalar = Kokkos::create_mirror_view(typename DeviceType::memory_space(), hexNodesHost);
    Kokkos::deep_copy(hexNodes_scalar, hexNodesHost);

    DynRankViewPointValueType ConstructWithLabelPointView(hexNodes, 8, 3);
    RealSpaceTools<DeviceType>::clone(hexNodes, hexNodes_scalar);

    // Array dimensions
    const ordinal_type numFields = hexBasis.getCardinality();
    const ordinal_type numPoints = hexNodes.extent(0);
    const ordinal_type spaceDim  = hexBasis.getBaseCellTopology().getDimension();

    {
      DynRankViewOutValueType ConstructWithLabelOutView(vals, numFields, numPoints, spaceDim);

      // exception #1: GRAD cannot be applied to HCURL functions
      INTREPID2_TEST_ERROR_EXPECTED( hexBasis.getValues(vals, hexNodes, OPERATOR_GRAD) );

      // exception #2: DIV cannot be applied to HCURL functions
      INTREPID2_TEST_ERROR_EXPECTED( hexBasis.getValues(vals, hexNodes, OPERATOR_DIV) );
    }

    {
      // Exceptions 3-7: all bf tags/bf Ids below are wrong and should cause getDofOrdinal() and
      // getDofTag() to access invalid array elements thereby causing bounds check exception
      //
      INTREPID2_TEST_ERROR_EXPECTED( hexBasis.getDofOrdinal(3,0,0) );
      INTREPID2_TEST_ERROR_EXPECTED( hexBasis.getDofOrdinal(1,1,1) );
      INTREPID2_TEST_ERROR_EXPECTED( hexBasis.getDofOrdinal(0,4,1) );
      INTREPID2_TEST_ERROR_EXPECTED( hexBasis.getDofTag(12) );
      INTREPID2_TEST_ERROR_EXPECTED( hexBasis.getDofTag(-1) );
    }

    {
      // Exceptions 8-15 test exception handling with incorrectly dimensioned input/output arrays
      DynRankViewOutValueType ConstructWithLabelOutView(vals, numFields, numPoints, spaceDim);

      {
        // Exception #8: input points array must be of rank-2
        DynRankViewPointValueType ConstructWithLabelPointView(badPoints,4,5,3);
        INTREPID2_TEST_ERROR_EXPECTED( hexBasis.getValues(vals, badPoints, OPERATOR_VALUE));
      }

      {
        // Exception #9: dimension 1 in the input point array must equal space dimension of the cell
        DynRankViewPointValueType ConstructWithLabelPointView(badPoints,4,2);
        INTREPID2_TEST_ERROR_EXPECTED( hexBasis.getValues(vals, badPoints, OPERATOR_VALUE));
      }

      {
        // Exception #10: output values must be of rank-3 for OPERATOR_VALUE
        DynRankViewOutValueType ConstructWithLabelOutView(badVals,4,3);
        INTREPID2_TEST_ERROR_EXPECTED( hexBasis.getValues(badVals, hexNodes, OPERATOR_VALUE));
      }

      {
        // Exception #11: output values must be of rank-3 for OPERATOR_CURL
        DynRankViewOutValueType ConstructWithLabelOutView(badVals,4,3);
        INTREPID2_TEST_ERROR_EXPECTED( hexBasis.getValues(badVals, hexNodes, OPERATOR_CURL));
      }

      {
        // Exception #12: incorrect 0th dimension of output array (must equal number of basis functions)
        DynRankViewOutValueType ConstructWithLabelOutView(badVals,numFields+1,numPoints,3);
        INTREPID2_TEST_ERROR_EXPECTED( hexBasis.getValues(badVals, hexNodes, OPERATOR_VALUE));
      }

      {
        // Exception #13: incorrect 1st dimension of output array (must equal number of points)
        DynRankViewOutValueType ConstructWithLabelOutView(badVals,numFields,numPoints+1,3);
        INTREPID2_TEST_ERROR_EXPECTED( hexBasis.getValues(badVals, hexNodes, OPERATOR_VALUE));
      }

      {
        // Exception #14: incorrect 2nd dimension of output array (must equal the space dimension)
        DynRankViewOutValueType ConstructWithLabelOutView(badVals,numFields,numPoints,4);
        INTREPID2_TEST_ERROR_EXPECTED( hexBasis.getValues(badVals, hexNodes, OPERATOR_VALUE));
      }

      {
        // Exception #15: incorrect 2nd dimension of output array (must equal the space dimension)
        DynRankViewOutValueType ConstructWithLabelOutView(badVals,numFields,numPoints,4);
        INTREPID2_TEST_ERROR_EXPECTED( hexBasis.getValues(badVals, hexNodes, OPERATOR_CURL));
      }
    }

    {
      // Exception #16: D2 cannot be applied to HCURL functions
      DynRankViewOutValueType ConstructWithLabelOutView(badVals,numFields,numPoints,4);
      INTREPID2_TEST_ERROR_EXPECTED( hexBasis.getValues(badVals, hexNodes, OPERATOR_D2));
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
  << "| TEST 2: Testing OPERATOR_VALUE (Kronecker property using dof coefficients)  |\n"
  << "===============================================================================\n";

  try {

    const ordinal_type order = std::min(3, maxOrder);
    HexBasisType hexBasis(order);

    const ordinal_type numFields = hexBasis.getCardinality();
    DynRankViewScalarValueType ConstructWithLabel(dofCoords_scalar, numFields, dim);
    hexBasis.getDofCoords(dofCoords_scalar);

    DynRankViewScalarValueType ConstructWithLabel(dofCoeffs, numFields, dim);
    hexBasis.getDofCoeffs(dofCoeffs);

    DynRankViewPointValueType ConstructWithLabelPointView(dofCoords, numFields , dim);
    RealSpaceTools<DeviceType>::clone(dofCoords,dofCoords_scalar);

    DynRankViewOutValueType ConstructWithLabelOutView(basisAtDofCoords, numFields, numFields, dim);
    hexBasis.getValues(basisAtDofCoords, dofCoords, OPERATOR_VALUE);

    auto h_dofCoords = Kokkos::create_mirror_view(dofCoords);
    Kokkos::deep_copy(h_dofCoords, dofCoords);

    auto h_dofCoeffs = Kokkos::create_mirror_view(dofCoeffs);
    Kokkos::deep_copy(h_dofCoeffs, dofCoeffs);

    auto h_basisAtDofCoords = Kokkos::create_mirror_view(basisAtDofCoords);
    Kokkos::deep_copy(h_basisAtDofCoords, basisAtDofCoords);

    for (int i=0;i<numFields;i++) {
      for (int j=0;j<numFields;j++) {

        OutValueType dofValue = 0.0;
        for(ordinal_type d=0;d<dim;++d)
          dofValue += h_basisAtDofCoords(i,j,d)*h_dofCoeffs(j,d);

        // check values
        const scalar_type expected_dofValue = (i == j);
        if (std::abs(dofValue - expected_dofValue) > tol) {
          errorFlag++;
          std::stringstream ss;
          ss << "\nValue of basis function " << i << " at (" << h_dofCoords(i,0) << ", " << h_dofCoords(i,1) << ", "<< h_dofCoords(i,2) << ") is " << dofValue << " but should be " << expected_dofValue << "\n";
          *outStream << ss.str();
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
  << "| TEST 3: Testing OPERATOR_VALUE (Kronecker property)                         |\n"
  << "===============================================================================\n";

  try {

    const ordinal_type order = std::min(3, maxOrder);
    HexBasisType hexBasis(order);

    shards::CellTopology hex_8(shards::getCellTopologyData<shards::Hexahedron<8> >());
    const ordinal_type numFields = hexBasis.getCardinality();
    DynRankViewScalarValueType ConstructWithLabel(dofCoords_scalar, numFields, dim);
    hexBasis.getDofCoords(dofCoords_scalar);

    DynRankViewPointValueType ConstructWithLabelPointView(dofCoords, numFields , dim);
    RealSpaceTools<DeviceType>::clone(dofCoords,dofCoords_scalar);

    DynRankViewOutValueType ConstructWithLabelOutView(basisAtDofCoords, numFields, numFields, dim);
    hexBasis.getValues(basisAtDofCoords, dofCoords, OPERATOR_VALUE);

    auto h_basisAtDofCoords = Kokkos::create_mirror_view(basisAtDofCoords);
    Kokkos::deep_copy(h_basisAtDofCoords, basisAtDofCoords);

    // These are tensor product basis functions where the local edge/basis orientation
    // is assumed to be in the positive x, y, or z direction.
    // For simplicity we do not need to compute tangents, but can just access values
    // from the correct indices assuming that tangents are (1,0,0), (0,1,0), or (0,0,1)

    // test for Kronecker property
    for (int i=0;i<numFields;i++) {
      for (int j=0;j<numFields;j++) {

        // initialize
        OutValueType dofValue = 0;

        // get index given that dofs are arranged by dimension, all x dofs, all y dofs, all z dofs
        auto dofsPerDim = numFields/dim;
        auto k_ind = j/dofsPerDim;
        dofValue = h_basisAtDofCoords(i,j,k_ind);

        // check values
        if ( i==j && std::abs( dofValue - 1.0 ) > tol ) {
          errorFlag++;
          *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
          *outStream << " Basis function " << i << " does not have unit value at its node (" << dofValue <<")\n";
        }
        if ( i!=j && std::abs( dofValue ) > tol ) {
          errorFlag++;
          *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
          *outStream << " Basis function " << i << " does not vanish at node " << j << "(" << dofValue <<")\n";
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
  << "| TEST 4: Test correctness of basis function values and curls                 |\n"
  << "===============================================================================\n";

  // VALUE: Each row pair gives the 12x3 correct basis set values at an evaluation point: (P,F,D) layout
  const scalar_type basisValues[] = {
      1,0,0, 1,0,0, 0,0,0, 0,0,0, 0,0,0, 0,0,0, 0,0,0, 0,0,0,
      0,0,0, 0,0,0, 1,0,0, 1,0,0, 0,0,0, 0,0,0, 0,0,0, 0,0,0,
      0,0,0, 0,0,0, 0,0,0, 0,0,0, 1,0,0, 1,0,0, 0,0,0, 0,0,0,
      0,0,0, 0,0,0, 0,0,0, 0,0,0, 0,0,0, 0,0,0, 1,0,0, 1,0,0,
      0,1,0, 0,0,0, 0,1,0, 0,0,0, 0,0,0, 0,0,0, 0,0,0, 0,0,0,
      0,0,0, 0,1,0, 0,0,0, 0,1,0, 0,0,0, 0,0,0, 0,0,0, 0,0,0,
      0,0,0, 0,0,0, 0,0,0, 0,0,0, 0,1,0, 0,0,0, 0,1,0, 0,0,0,
      0,0,0, 0,0,0, 0,0,0, 0,0,0, 0,0,0, 0,1,0, 0,0,0, 0,1,0,
      0,0,1, 0,0,0, 0,0,0, 0,0,0, 0,0,1, 0,0,0, 0,0,0, 0,0,0,
      0,0,0, 0,0,1, 0,0,0, 0,0,0, 0,0,0, 0,0,1, 0,0,0, 0,0,0,
      0,0,0, 0,0,0, 0,0,1, 0,0,0, 0,0,0, 0,0,0, 0,0,1, 0,0,0,
      0,0,0, 0,0,0, 0,0,0, 0,0,1, 0,0,0, 0,0,0, 0,0,0, 0,0,1
  };

  // CURL
  const scalar_type basisCurls[] = {   
      0,-0.5,0.5, 0,-0.5,0.5, 0,0,0.5, 0,0,0.5, 0,-0.5,0, 0,-0.5,0, 0,0,0, 0,0,0,
      0,0,-0.5, 0,0,-0.5, 0,-0.5,-0.5, 0,-0.5,-0.5, 0,0,0, 0,0,0, 0,-0.5,0, 0,-0.5,0,
      0,0.5,0, 0,0.5,0, 0,0,0, 0,0,0, 0,0.5,0.5, 0,0.5,0.5, 0,0,0.5, 0,0,0.5,
      0,0,0, 0,0,0, 0,0.5,0, 0,0.5,0, 0,0,-0.5, 0,0,-0.5, 0,0.5,-0.5, 0,0.5,-0.5,
      // y-component basis functions
      // first y-component bf is (0,1/4(1-x)(1-z),0)
      // curl is (1/4(1-x),0,-1/4(1-z))
      0.5,0,-0.5, 0,0,-0.5, 0.5,0,-0.5, 0,0,-0.5, 0.5,0,0, 0,0,0, 0.5,0,0, 0,0,0,
      // second y-component bf is (0,1/4(1+x)(1-z),0)
      // curl is (1/4(1+x),0,1/4(1-z))
      0,0,0.5, 0.5,0,0.5, 0,0,0.5, 0.5,0,0.5, 0,0,0, 0.5,0,0, 0,0,0, 0.5,0,0,
      // third y-component bf is (0,1/4(1-x)(1+z),0)
      // curl is (-1/4(1-x),0,-1/4(1+z))
      -0.5,0,0, 0,0,0, -0.5,0,0, 0,0,0, -0.5,0,-0.5, 0,0,-0.5, -0.5,0,-0.5, 0,0,-0.5,
      // fourth y-component bf is (0,1/4(1+x)(1+z),0)
      // curl is (-1/4(1+x),0,1/4(1+z))
      0.0,0,0, -0.5,0,0, 0.0,0,0, -0.5,0,0, 0.0,0,0.5, -0.5,0,0.5, 0.0,0,0.5, -0.5,0,0.5,
      // first z-component bf is (0,0,1/4(1-x)(1-y))
      // curl is (-1/4(1-x),1/4(1-y),0)
      -0.5,0.5,0, 0,0.5,0, -0.5,0,0, 0,0,0, -0.5,0.5,0, 0,0.5,0, -0.5,0,0, 0,0,0,
      // second z-component bf is (0,0,1/4(1+x)(1-y))
      // curl is (-1/4(1+x),1/4(1-y),0)
      0.0,-0.5,0, -0.5,-0.5,0, 0,0,0, -0.5,0,0, 0,-0.5,0, -0.5,-0.5,0, 0,0,0, -0.5,0,0,
      // third z-component bf is (0,0,1/4(1-x)(1+y))
      // curl is (1/4(1-x),1/4(1+y),0)
      0.5,0,0, 0,0,0, 0.5,0.5,0, 0,0.5,0, 0.5,0,0, 0,0,0, 0.5,0.5,0, 0,0.5,0,
      // fourth z-component bf is (0,0,1/4(1+x)(1+y))
      // curl is (1/4(1+x),-1/4(1+y),0)
      0,0,0, 0.5,0,0, 0,-0.5,0, 0.5,-0.5,0, 0,0,0, 0.5,0,0, 0,-0.5,0, 0.5,-0.5,0
  };


  try {

    const ordinal_type order = 1;
    HexBasisType hexBasis(order);

    // Array of reference hex nodes - used for evaluation of basis
    DynRankViewHostScalarValueType ConstructWithLabel(hexNodesHost, 8, 3);
    DynRankViewPointValueType ConstructWithLabelPointView(hexNodes, 8, 3);

    hexNodesHost(0,0) = -1.0; hexNodesHost(0,1) = -1.0; hexNodesHost(0,2) = -1.0;
    hexNodesHost(1,0) =  1.0; hexNodesHost(1,1) = -1.0; hexNodesHost(1,2) = -1.0;
    hexNodesHost(2,0) = -1.0; hexNodesHost(2,1) =  1.0; hexNodesHost(2,2) = -1.0;
    hexNodesHost(3,0) =  1.0; hexNodesHost(3,1) =  1.0; hexNodesHost(3,2) = -1.0;
    hexNodesHost(4,0) = -1.0; hexNodesHost(4,1) = -1.0; hexNodesHost(4,2) =  1.0;
    hexNodesHost(5,0) =  1.0; hexNodesHost(5,1) = -1.0; hexNodesHost(5,2) =  1.0;
    hexNodesHost(6,0) = -1.0; hexNodesHost(6,1) =  1.0; hexNodesHost(6,2) =  1.0;
    hexNodesHost(7,0) =  1.0; hexNodesHost(7,1) =  1.0; hexNodesHost(7,2) =  1.0;

    auto hexNodes_scalar = Kokkos::create_mirror_view(typename DeviceType::memory_space(), hexNodesHost);
    Kokkos::deep_copy(hexNodes_scalar, hexNodesHost);
    RealSpaceTools<DeviceType>::clone(hexNodes, hexNodes_scalar);

    // Array dimensions
    const ordinal_type numFields = hexBasis.getCardinality();
    const ordinal_type numPoints = hexNodes.extent(0);
    const ordinal_type spaceDim  = hexBasis.getBaseCellTopology().getDimension();

    *outStream << " -- Testing OPERATOR_VALUE \n";
    {
      // Check VALUE of basis functions
      DynRankViewOutValueType ConstructWithLabelOutView(vals, numFields, numPoints, spaceDim);
      hexBasis.getValues(vals, hexNodes, OPERATOR_VALUE);
      auto vals_host = Kokkos::create_mirror_view(vals);
      Kokkos::deep_copy(vals_host, vals);

      for (ordinal_type i = 0; i < numFields; i++) {
        for (ordinal_type j = 0; j < numPoints; j++) {
          for (ordinal_type k = 0; k < spaceDim; k++) {

            // compute offset for (P,F,D) data layout: indices are P->j, F->i, D->k
            const ordinal_type l = k + i * spaceDim * numPoints + j * spaceDim;
            if (std::abs(vals_host(i,j,k) - basisValues[l]) > tol) {
              errorFlag++;
              *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";

              // Output the multi-index of the value where the error is:
              *outStream << " At multi-index { ";
              *outStream << i << " ";*outStream << j << " ";*outStream << k << " ";
              *outStream << "}  computed value: " << vals(i,j,k)
                                    << " but reference value: " << basisValues[l] << "\n";
            }
          }
        }
      }
    }
    *outStream << " -- Testing OPERATOR_CURL \n";
    {
      // Check CURL of basis function:
      DynRankViewOutValueType ConstructWithLabelOutView(curls, numFields, numPoints, spaceDim);
      hexBasis.getValues(curls, hexNodes, OPERATOR_CURL);
      auto curls_host = Kokkos::create_mirror_view(curls);
      Kokkos::deep_copy(curls_host, curls);

      for (ordinal_type i = 0; i < numFields; ++i) {
        for (ordinal_type j = 0; j < numPoints; ++j) {
          for (ordinal_type k = 0; k < spaceDim; k++) {
            const ordinal_type l = k + i * spaceDim * numPoints + j * spaceDim;
            if (std::abs(curls_host(i,j,k) - basisCurls[l]) > tol) {
              errorFlag++;
              *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";

              // Output the multi-index of the curl where the error is:
              *outStream << " At multi-index { ";
              *outStream << i << " ";*outStream << j << " ";*outStream << k << " ";
              *outStream << "}  computed curl component: " << curls(i,j,k)
                                    << " but reference curl component: " << basisCurls[l] << "\n";
            }
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
    const ordinal_type order = std::min(3, maxOrder);
    HexBasisType hexBasis(order);
    
    const EFunctionSpace fs = hexBasis.getFunctionSpace();
    
    if (fs != FUNCTION_SPACE_HCURL)
    {
      *outStream << std::setw(70) << "------------- TEST FAILURE! -------------" << "\n";
      
      // Output the multi-index of the value where the error is:
      *outStream << " Expected a function space of FUNCTION_SPACE_HCURL (enum value " << FUNCTION_SPACE_HCURL << "),";
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
