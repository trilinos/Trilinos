// @HEADER
// *****************************************************************************
//                           Intrepid2 Package
//
// Copyright 2007 NTESS and the Intrepid2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/** \file   test_01.hpp
    \brief  Unit tests for the Intrepid2::HDIV_QUAD_In_FEM class.
    \author Created by P. Bochev, D. Ridzal, K. Peterson, Kyungjoo Kim
 */


#include "Intrepid2_config.h"

#ifdef HAVE_INTREPID2_DEBUG
#define INTREPID2_TEST_FOR_DEBUG_ABORT_OVERRIDE_TO_CONTINUE
#endif

#include "Intrepid2_Types.hpp"
#include "Intrepid2_Utils.hpp"

#include "Intrepid2_PointTools.hpp"
#include "Intrepid2_HDIV_QUAD_In_FEM.hpp"

#include "Teuchos_oblackholestream.hpp"
#include "Teuchos_RCP.hpp"

#include "packages/intrepid2/unit-test/Discretization/Basis/Macros.hpp"

namespace Intrepid2 {

namespace Test {

template<typename OutValueType, typename PointValueType, typename DeviceType>
int HDIV_QUAD_In_FEM_Test01(const bool verbose) {

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
  << "|                 Unit Test (Basis_HDIV_QUAD_In_FEM)                          |\n"
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

  typedef Kokkos::DynRankView<PointValueType,DeviceType> DynRankViewPointValueType;
  typedef Kokkos::DynRankView<OutValueType,DeviceType> DynRankViewOutValueType;
  typedef typename ScalarTraits<OutValueType>::scalar_type scalar_type;
  typedef Kokkos::DynRankView<scalar_type, DeviceType> DynRankViewScalarValueType;
  typedef Kokkos::DynRankView<scalar_type, HostSpaceType> DynRankViewHostScalarValueType;

  const scalar_type tol = tolerence();
  int errorFlag = 0;

  typedef Basis_HDIV_QUAD_In_FEM<DeviceType,OutValueType,PointValueType> QuadBasisType;
  constexpr ordinal_type maxOrder = Parameters::MaxOrder ;

  *outStream
  << "\n"
  << "===============================================================================\n"
  << "| TEST 1: Basis creation, exceptions tests                                    |\n"
  << "===============================================================================\n";

  try{

#ifdef HAVE_INTREPID2_DEBUG
    ordinal_type nthrow = 0, ncatch = 0;
    constexpr  ordinal_type order = 5;
    if(order <= maxOrder) {
      QuadBasisType quadBasis(order);

      // Define array containing array of nodes to evaluate
      DynRankViewPointValueType ConstructWithLabelPointView(quadNodes, 9, 2);

      // Generic array for the output values; needs to be properly resized depending on the operator type
      const ordinal_type numFields = quadBasis.getCardinality();
      const ordinal_type numPoints = quadNodes.extent(0);
      const ordinal_type spaceDim  = quadBasis.getBaseCellTopology().getDimension();

      {
        DynRankViewOutValueType ConstructWithLabelOutView(vals, numFields, numPoints, spaceDim);

        // exception #1: GRAD cannot be applied to HDIV functions
        INTREPID2_TEST_ERROR_EXPECTED( quadBasis.getValues(vals, quadNodes, OPERATOR_GRAD) );

        // exception #2: CURL cannot be applied to HDIV functions
        INTREPID2_TEST_ERROR_EXPECTED( quadBasis.getValues(vals, quadNodes, OPERATOR_CURL) );
      }

      // Exceptions 3-7: all bf tags/bf Ids below are wrong and should cause getDofOrdinal() and
      // getDofTag() to access invalid array elements thereby causing bounds check exception
      {
        INTREPID2_TEST_ERROR_EXPECTED( quadBasis.getDofOrdinal(3,0,0) );
        INTREPID2_TEST_ERROR_EXPECTED( quadBasis.getDofOrdinal(1,0,5) );
        INTREPID2_TEST_ERROR_EXPECTED( quadBasis.getDofOrdinal(0,4,0) );
        INTREPID2_TEST_ERROR_EXPECTED( quadBasis.getDofTag(numFields) );
        INTREPID2_TEST_ERROR_EXPECTED( quadBasis.getDofTag(-1) );
      }

      // Exceptions 8-15 test exception handling with incorrectly dimensioned input/output arrays
      {
        DynRankViewOutValueType ConstructWithLabelOutView(vals, numFields, numPoints, spaceDim);
        {
          // exception #8: input points array must be of rank-2
          DynRankViewPointValueType ConstructWithLabelPointView(badPoints, 4, 5, 3);
          INTREPID2_TEST_ERROR_EXPECTED( quadBasis.getValues(vals, badPoints, OPERATOR_VALUE) );
        }
        {
          // exception #9 dimension 1 in the input point array must equal space dimension of the cell
          DynRankViewPointValueType ConstructWithLabelPointView(badPoints, 4, 3);
          INTREPID2_TEST_ERROR_EXPECTED( quadBasis.getValues(vals, badPoints, OPERATOR_VALUE) );
        }
        {
          // exception #10 output values must be of rank-3 for OPERATOR_VALUE in 2D
          DynRankViewOutValueType ConstructWithLabelOutView(badVals, 4, 3);
          INTREPID2_TEST_ERROR_EXPECTED( quadBasis.getValues(badVals, quadNodes, OPERATOR_VALUE) );
        }
        {
          // exception #11 output values must be of rank-2 for OPERATOR_DIV
          DynRankViewOutValueType ConstructWithLabelOutView(badVals, 4, 3, 2);
          INTREPID2_TEST_ERROR_EXPECTED( quadBasis.getValues(badVals, quadNodes, OPERATOR_DIV) );
        }
        {
          // exception #12 incorrect 0th dimension of output array (must equal number of basis functions)
          DynRankViewOutValueType ConstructWithLabelOutView(badVals, quadBasis.getCardinality() + 1, quadNodes.extent(0), quadBasis.getBaseCellTopology().getDimension());
          INTREPID2_TEST_ERROR_EXPECTED( quadBasis.getValues(badVals, quadNodes, OPERATOR_VALUE) ) ;
        }
        {
          // exception #13 incorrect 1st  dimension of output array (must equal number of points)
          DynRankViewOutValueType ConstructWithLabelOutView(badVals, quadBasis.getCardinality(), quadNodes.extent(0) + 1, quadBasis.getBaseCellTopology().getDimension() );
          INTREPID2_TEST_ERROR_EXPECTED( quadBasis.getValues(badVals, quadNodes, OPERATOR_VALUE) ) ;
        }
        {
          // exception #14: incorrect 2nd dimension of output array for VALUE (must equal the space dimension)
          DynRankViewOutValueType ConstructWithLabelOutView(badVals, quadBasis.getCardinality(), quadNodes.extent(0), quadBasis.getBaseCellTopology().getDimension() - 1);
          INTREPID2_TEST_ERROR_EXPECTED( quadBasis.getValues(badVals, quadNodes, OPERATOR_VALUE) ) ;
        }
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
  << "| TEST 2: Testing OPERATOR_VALUE (Kronecker property using dof coeffs)        |\n"
  << "===============================================================================\n";

  try {

    const ordinal_type order = std::min(3, maxOrder);
    QuadBasisType quadBasis(order);
    const ordinal_type spaceDim  = quadBasis.getBaseCellTopology().getDimension();

    const ordinal_type numFields = quadBasis.getCardinality();
    DynRankViewScalarValueType ConstructWithLabel(dofCoords_scalar, numFields, spaceDim);
    quadBasis.getDofCoords(dofCoords_scalar);

    DynRankViewScalarValueType ConstructWithLabel(dofCoeffs, numFields, spaceDim);
    quadBasis.getDofCoeffs(dofCoeffs);

    DynRankViewPointValueType ConstructWithLabelPointView(dofCoords, numFields , spaceDim);
    RealSpaceTools<DeviceType>::clone(dofCoords,dofCoords_scalar);

    DynRankViewOutValueType ConstructWithLabelOutView(basisAtDofCoords, numFields, numFields, spaceDim);
    quadBasis.getValues(basisAtDofCoords, dofCoords, OPERATOR_VALUE);

    auto h_dofCoords = Kokkos::create_mirror_view(dofCoords);
    Kokkos::deep_copy(h_dofCoords, dofCoords);

    auto h_dofCoeffs = Kokkos::create_mirror_view(dofCoeffs);
    Kokkos::deep_copy(h_dofCoeffs, dofCoeffs);

    auto h_basisAtDofCoords = Kokkos::create_mirror_view(basisAtDofCoords);
    Kokkos::deep_copy(h_basisAtDofCoords, basisAtDofCoords);

    // test for Kronecker property
    for (int i=0;i<numFields;i++) {
      for (int j=0;j<numFields;j++) {

        OutValueType dofValue = 0.0;
        for(ordinal_type d=0;d<spaceDim;++d)
          dofValue += h_basisAtDofCoords(i,j,d)*h_dofCoeffs(j,d);

        // check values
        const scalar_type expected_dofValue = (i == j);
        if (std::abs(dofValue - expected_dofValue) > tol) {
          errorFlag++;
          std::stringstream ss;
          ss << "\nValue of basis function " << i << " at (" << h_dofCoords(j,0) << ", " << h_dofCoords(j,1) << ") is " << dofValue << " but should be " << expected_dofValue << "\n";
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
    QuadBasisType quadBasis(order);
    const ordinal_type spaceDim  = quadBasis.getBaseCellTopology().getDimension();

    const ordinal_type numFields = quadBasis.getCardinality();
    DynRankViewScalarValueType ConstructWithLabel(dofCoords_scalar, numFields, spaceDim);
    quadBasis.getDofCoords(dofCoords_scalar);

    DynRankViewPointValueType ConstructWithLabelPointView(dofCoords, numFields , spaceDim);
    RealSpaceTools<DeviceType>::clone(dofCoords,dofCoords_scalar);

    DynRankViewOutValueType ConstructWithLabelOutView(basisAtDofCoords, numFields, numFields, spaceDim);
    quadBasis.getValues(basisAtDofCoords, dofCoords, OPERATOR_VALUE);

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

        // get index given that dofs are arranged by dimension, all x dofs, all y dofs
        auto dofsPerDim = numFields/spaceDim;
        auto k_ind = spaceDim-1-j/dofsPerDim;
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
  << "| TEST 5: correctness of tag to enum and enum to tag lookups                  |\n"
  << "===============================================================================\n";

  try {
    const ordinal_type order = std::min(5, maxOrder);
    QuadBasisType quadBasis(order);

    const ordinal_type numFields = quadBasis.getCardinality();
    const auto allTags = quadBasis.getAllDofTags();

    // Loop over all tags, lookup the associated dof enumeration and then lookup the tag again
    const ordinal_type dofTagSize = allTags.extent(0);
    for (ordinal_type i=0;i<dofTagSize;++i) {
      const ordinal_type bfOrd = quadBasis.getDofOrdinal(allTags(i,0), allTags(i,1), allTags(i,2));

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
            << allTags(i,3) << "}) = " << bfOrd <<" but ";
        *outStream << " getDofTag(" << bfOrd << ") = {"
            << myTag(0) << ", "
            << myTag(1) << ", "
            << myTag(2) << ", "
            << myTag(3) << "}\n";
      }
    }

    // Now do the same but loop over basis functions
    for(ordinal_type bfOrd=0;bfOrd<numFields;++bfOrd) {
      const auto myTag  = quadBasis.getDofTag(bfOrd);
      const ordinal_type myBfOrd = quadBasis.getDofOrdinal(myTag(0), myTag(1), myTag(2));
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
  << "| TEST 6: correctness of basis function values                                |\n"
  << "===============================================================================\n";

  outStream->precision(20);

  // VALUE: Each row pair gives the correct basis set values at an evaluation point: (F, P, D) layout
  const scalar_type basisValues[] = {
      0.0, 1.0, 0.0, 1.0, 0.0, 1.0,   0.0, 0.5, 0.0, 0.5, 0.0, 0.5,   0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0,   0.0, 0.5, 0.0, 0.5, 0.0, 0.5,   0.0, 1.0, 0.0, 1.0, 0.0, 1.0,
      1.0, 0.0, 0.5, 0.0, 0.0, 0.0,   1.0, 0.0, 0.5, 0.0, 0.0, 0.0,   1.0, 0.0, 0.5, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.5, 0.0, 1.0, 0.0,   0.0, 0.0, 0.5, 0.0, 1.0, 0.0,   0.0, 0.0, 0.5, 0.0, 1.0, 0.0
  };

  // DIV: correct values in (P,F) format
  const scalar_type basisDivs[] = {
      -0.5, -0.5, -0.5, -0.5, -0.5, -0.5, -0.5, -0.5, -0.5,
      0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5,
      -0.5, -0.5, -0.5, -0.5, -0.5, -0.5, -0.5, -0.5, -0.5,
      0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5,
  };

  try{
    const ordinal_type order = 1;
    QuadBasisType quadBasis(order);

    // Define array containing array of nodes to evaluate
    DynRankViewHostScalarValueType ConstructWithLabel(quadNodesHost, 9, 2);
    DynRankViewPointValueType ConstructWithLabelPointView(quadNodes, 9, 2);

    quadNodesHost(0,0) = -1.0;  quadNodesHost(0,1) = -1.0;
    quadNodesHost(1,0) =  0.0;  quadNodesHost(1,1) = -1.0;
    quadNodesHost(2,0) =  1.0;  quadNodesHost(2,1) = -1.0;
    quadNodesHost(3,0) = -1.0;  quadNodesHost(3,1) =  0.0;
    quadNodesHost(4,0) =  0.0;  quadNodesHost(4,1) =  0.0;
    quadNodesHost(5,0) =  1.0;  quadNodesHost(5,1) =  0.0;
    quadNodesHost(6,0) = -1.0;  quadNodesHost(6,1) =  1.0;
    quadNodesHost(7,0) =  0.0;  quadNodesHost(7,1) =  1.0;
    quadNodesHost(8,0) =  1.0;  quadNodesHost(8,1) =  1.0;

    auto quadNodes_scalar = Kokkos::create_mirror_view(typename DeviceType::memory_space(), quadNodesHost);
    Kokkos::deep_copy(quadNodes_scalar, quadNodesHost);
    RealSpaceTools<DeviceType>::clone(quadNodes, quadNodes_scalar);

    // Generic array for the output values; needs to be properly resized depending on the operator type
    const ordinal_type numFields = quadBasis.getCardinality();
    const ordinal_type numPoints = quadNodes.extent(0);
    const ordinal_type spaceDim  = quadBasis.getBaseCellTopology().getDimension();

    *outStream << " -- Testing OPERATOR_VALUE \n";
    {
      // Check VALUE of basis functions: resize vals to rank-3 container:
      DynRankViewOutValueType ConstructWithLabelOutView(vals, numFields, numPoints, spaceDim);
      quadBasis.getValues(vals, quadNodes, OPERATOR_VALUE);
      auto vals_host = Kokkos::create_mirror_view(vals);
      Kokkos::deep_copy(vals_host, vals);


      for (ordinal_type i = 0; i < numFields; ++i) {
        for (ordinal_type j = 0; j < numPoints; ++j) {
          for (ordinal_type k = 0; k < spaceDim; ++k) {

            // compute offset for (P,F,D) data layout: indices are P->j, F->i, D->k
            const ordinal_type l = i * spaceDim * numPoints + j * spaceDim + k;
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
    *outStream << " -- Testing OPERATOR_DIV \n";
    {
      // Check DIV of basis function: resize vals to rank-2 container
      DynRankViewOutValueType ConstructWithLabelOutView(vals, numFields, numPoints);
      quadBasis.getValues(vals, quadNodes, OPERATOR_DIV);

      auto vals_host = Kokkos::create_mirror_view(vals);
      Kokkos::deep_copy(vals_host, vals);

      for (ordinal_type i = 0; i < numFields; ++i) {
        for (ordinal_type j = 0; j < numPoints; ++j) {
          const ordinal_type l =  j + i * numPoints;
          if (std::abs(vals_host(i,j) - basisDivs[l]) > tol) {
            errorFlag++;
            *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";

            // Output the multi-index of the value where the error is:
            *outStream << " At multi-index { ";
            *outStream << i << " ";*outStream << j << " ";
            *outStream << "}  computed div component: " << vals(i,j)
                                   << " but reference div component: " << basisDivs[l] << "\n";
          }
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
  << "| TEST 7: Function Space is Correct                                           |\n"
  << "===============================================================================\n";
  
  try {
    const int order = 2;
    QuadBasisType quadBasis(order);
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
