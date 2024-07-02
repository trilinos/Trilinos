// @HEADER
// *****************************************************************************
//                           Intrepid2 Package
//
// Copyright 2007 NTESS and the Intrepid2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/** \file test_01.hpp
    \brief  Unit tests for the Intrepid2::HGRAD_QUAD_Cn_FEM class.
    \author Created by P. Bochev, D. Ridzal, K. Peterson, Kyungjoo Kim
 */


#include "Intrepid2_config.h"

#ifdef HAVE_INTREPID2_DEBUG
#define INTREPID2_TEST_FOR_DEBUG_ABORT_OVERRIDE_TO_CONTINUE
#endif

#include "Intrepid2_Types.hpp"
#include "Intrepid2_Utils.hpp"

#include "Intrepid2_PointTools.hpp"
#include "Intrepid2_HGRAD_QUAD_Cn_FEM.hpp"

#include "Teuchos_oblackholestream.hpp"
#include "Teuchos_RCP.hpp"

#include "packages/intrepid2/unit-test/Discretization/Basis/Macros.hpp"
#include "packages/intrepid2/unit-test/Discretization/Basis/Setup.hpp"

namespace Intrepid2 {

namespace Test {

using HostSpaceType = Kokkos::DefaultHostExecutionSpace;

template<typename OutValueType, typename PointValueType, typename DeviceType>
int HGRAD_QUAD_Cn_FEM_Test01(const bool verbose) {

  //! Create an execution space instance.
  const auto space = Kokkos::Experimental::partition_space(typename DeviceType::execution_space {}, 1)[0];

  //! Setup test output stream.
  Teuchos::RCP<std::ostream> outStream = setup_output_stream<DeviceType>(
    verbose, "Basis_HGRAD_QUAD_Cn_FEM", {
      "1) Conversion of Dof tags into Dof ordinals and back",
      "2) Basis values for VALUE, GRAD, CURL, and Dk operators"
  });

  Teuchos::oblackholestream oldFormatState;
  oldFormatState.copyfmt(std::cout);

  typedef Kokkos::DynRankView<PointValueType,DeviceType> DynRankViewPointValueType;
  typedef Kokkos::DynRankView<OutValueType,DeviceType> DynRankViewOutValueType;
  typedef typename ScalarTraits<OutValueType>::scalar_type scalar_type;
  typedef Kokkos::DynRankView<scalar_type, DeviceType> DynRankViewScalarValueType;      
  typedef Kokkos::DynRankView<scalar_type, HostSpaceType> DynRankViewHostScalarValueType;      

  const scalar_type tol = tolerence();
  int errorFlag = 0;

  typedef Basis_HGRAD_QUAD_Cn_FEM<DeviceType,OutValueType,PointValueType> QuadBasisType;
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


  try {

    *outStream
    << "\n"
    << "===============================================================================\n"
    << "| TEST 3: Testing Kronecker property of basis functions                       |\n"
    << "===============================================================================\n";


    const ordinal_type order = std::min(2,maxOrder);
    QuadBasisType quadBasis(order, POINTTYPE_WARPBLEND);
    constexpr ordinal_type dim=2;
    const ordinal_type basisCardinality = quadBasis.getCardinality();
    DynRankViewScalarValueType ConstructWithLabel(dofCoords, basisCardinality , dim);
    DynRankViewPointValueType ConstructWithLabelPointView(lattice, basisCardinality , dim);

    quadBasis.getDofCoords(dofCoords);
    RealSpaceTools<DeviceType>::clone(lattice,dofCoords);

    auto lattice_host = Kokkos::create_mirror_view(lattice);

    DynRankViewOutValueType ConstructWithLabelOutView(basisAtLattice, basisCardinality, basisCardinality);
    quadBasis.getValues(space, basisAtLattice, lattice, OPERATOR_VALUE);

    auto h_basisAtLattice = Kokkos::create_mirror_view(basisAtLattice);
    Kokkos::deep_copy(h_basisAtLattice, basisAtLattice);

    auto numInternalDofs = quadBasis.getDofCount(2,0);
    for(ordinal_type i=0; i<numInternalDofs; i++) {
      auto idof = quadBasis.getDofOrdinal(2,0,i);
      for(ordinal_type j=0; j<numInternalDofs; j++) {
        auto jdof = quadBasis.getDofOrdinal(2,0,j);
        if ( idof==jdof && std::abs( h_basisAtLattice(idof,jdof) - 1.0 ) > tol ) {
          errorFlag++;
          *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
          *outStream << " Basis function " << idof << " does not have unit value at its node (" << h_basisAtLattice(idof,jdof) <<")\n";
        }
        if ( i!=j && std::abs( h_basisAtLattice(idof,jdof) ) > tol ) {
          errorFlag++;
          *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
          *outStream << " Basis function " << idof << " does not vanish at node " << jdof << "\n";
          *outStream << " Basis function value is " << h_basisAtLattice(idof,jdof) << "\n";
        }
      }
    }


    // test for Kronecker property
    for (int i=0;i<basisCardinality;i++) {
      for (int j=0;j<basisCardinality;j++) {
        if ( i==j && std::abs( h_basisAtLattice(i,j) - 1.0 ) > tol ) {
          errorFlag++;
          *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
          *outStream << " Basis function " << i << " does not have unit value at its node (" << h_basisAtLattice(i,j) <<")\n";
        }
        if ( i!=j && std::abs( h_basisAtLattice(i,j) ) > tol ) {
          errorFlag++;
          *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
          *outStream << " Basis function " << i << " does not vanish at node " << j << "\n";
          *outStream << " Basis function value is " << h_basisAtLattice(i,j) << "\n";
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
  << "| TEST 3: correctness of basis function values                                |\n"
  << "===============================================================================\n";

  outStream->precision(20);

  // VALUE: Correct basis values in (F,P) format:
  const scalar_type basisValues[] = {
      1, 0, 0, 0, 0, 0, 0, 0, 0, -0.05333333333333334, \
      0, 0, 0, 0, 1, 0, 0, 0, 0, 0.4266666666666667, \
      0, 1, 0, 0, 0, 0, 0, 0, 0, 0.1066666666666667, \
      0, 0, 0, 0, 0, 0, 0, 1, 0, -0.07111111111111112 , \
      0, 0, 0, 0, 0, 0, 0, 0, 1, 0.5688888888888890, \
      0, 0, 0, 0, 0, 1, 0, 0, 0, 0.1422222222222222 ,\
      0, 0, 0, 1, 0, 0, 0, 0, 0, 0.01333333333333333, \
      0, 0, 0, 0, 0, 0, 1, 0, 0, -0.1066666666666667, \
      0, 0, 1, 0, 0, 0, 0, 0, 0, -0.02666666666666666 };

  // GRAD and D1: Correct gradients and D1 in (F,P,D) format
  // 9 basis functions, each evaluated at 10 points, with two
  // components at each point.
  // that looks like 10 per to me.
  const scalar_type basisGrads[] = {
      //
      -1.500000000000000, -1.500000000000000, 0.5000000000000000, 0, 0, 0, 0, 0.5000000000000000, -0.5000000000000000, 0, \
      0, 0, 0, 0, 0, -0.5000000000000000, 0, 0, -0.08000000000000002, 0.1222222222222222, \
      //
      2.000000000000000, 0, -2.000000000000000, 0, 0, 0, 0, 0, 0, -1.500000000000000, \
      0, 0, 0, 0.5000000000000000, 0, 0, 0, -0.5000000000000000,  -0.3199999999999999, -0.9777777777777779, \
      //
      -0.5000000000000000, 0, 1.500000000000000, -1.500000000000000, 0, 0.5000000000000000, 0, 0, 0.5000000000000000, 0, \
      0, -0.5000000000000000, 0, 0, 0, 0, 0, 0, 0.3999999999999999, -0.2444444444444444, \
      //
      0, 2.0, 0, 0, 0, 0, 0, -2.000000000000000, 0, 0, \
      0.5000000000000000, 0, 0, 0, -1.50, 0, -0.50, 0, -0.1066666666666667, -0.1333333333333333, \
      //
      0, 0, 0, 0, 0, 0, 0, 0, 0, 2.0,\
      -2.00, 0, 0, -2.0, 2.0, 0, 0, 0, -0.4266666666666667, 1.066666666666667 ,  \
      //
      0, 0, 0, 2.000000000000000, 0, -2.000000000000000, 0, 0, 0, 0, \
      1.5, 0, 0, 0, -0.5, 0, 0.5000000000000000, 0, 0.5333333333333334, 0.2666666666666666 ,  \
      //
      0, -0.5000000000000000, 0, 0, 0.5000000000000000, 0, -1.500000000000000, 1.500000000000000, 0, 0, \
      0, 0, -0.5000000000000000, 0, 0, 0.5000000000000000, 0, 0, 0.02000000000000000, 0.01111111111111112 , \
      //
      0, 0, 0, 0, -2.0, 0, 2.0, 0, 0, -0.50,				\
      0, 0, 0, 1.5, 0, 0, 0, 0.5000000000000000, 0.07999999999999997, -0.08888888888888888, \
      //
      0, 0, 0, -0.5000000000000000, 1.500000000000000, 1.500000000000000, -0.5000000000000000, 0, 0, 0, \
      0, 0.5000000000000000, 0.5000000000000000, 0, 0, 0, 0, 0, -0.09999999999999998, -0.02222222222222221 \
      //
  };

  // D2: Correct multiset of second order partials in (F,P,Dk) format. D2 cardinality = 3 for 2D
  // 10 quad points and 3 values per point, so
  // each bf consists of 30 values.
  const scalar_type basisD2[] = {
      1.0, 2.25, 1.0, 1.0, -0.75, 0, 0, 0.25, 0,  0, -0.75, 1.0, 1.0,  0.75, 0, 0, -0.25, 0, 0,  -0.25, 0, 0, 0.75, 1.0, 0,  0.25, 0, 0.48, 0.1833333333333334, -0.1111111111111111,
      //
      -2.0, -3.0, 0, -2.0, 3.0, 0, 0, -1.0, 0, \
      0, 1.0, 0, -2.0, 0, 1.0, 0, \
      1.0, 0, 0, 0, 1.0, 0, -1.0, \
      0, 0, 0, 1.0, -0.96, 0.7333333333333332, \
      0.8888888888888890, \
      //
      1.0, 0.75, 0, 1.0, -2.25, 1.0, 0, 0.75, 1.0, 0, -0.25, 0, \
      1.0, -0.75, 0, 0, -0.75, 1.0, 0, 0.25, 0, 0, 0.25, \
      0, 0, -0.25, 0, 0.48, -0.9166666666666666, 0.2222222222222222,
      //
      0, -3.0, -2.0, 0, 1.0, 0, 0, -1.0, 0, 0, 3.0, \
      -2.0, 0, -1.0, 0, 1.0, 0, 0, 0, 1.0, 0, 1.0, 0, -2.0, \
      1.0, 0, 0, 0.6400000000000001, -0.2000000000000001,  0.2222222222222222, \
      //
      0, 4.0, 0, 0, -4.0,	0, 0, 4.0, 0, 0, -4.0, 0, 0, 0, \
      -2.0, -2.0, 0, 0, 0, 0, -2.0, -2.0, 0, 0, -2.0, 0, \
      -2.0, -1.280000000000000, -0.7999999999999998, -1.777777777777778 ,
      //
      0, -1.0, 0, 0, 3.0,	-2.0, 0, -3.0, -2.0, 0, \
      1.0, 0, 0, 1.0, 0, 1.0, 0, -2.0, 0, -1.0, 0, 1.0, 0, \
      0, 1.0, 0, 0, 0.6400000000000001, 1.0, -0.4444444444444444, \
      //
      0, 0.75, 1.0, 0, -0.25, 0, 1.0, 0.75, 0, 1.0, -2.25, 1.0, 0, \
      0.25, 0, 0, 0.25, 0, 1.0, -0.75, 0, 0, -0.75, 1.0, 0, \
      -0.25, 0, -0.12, 0.01666666666666666,  -0.1111111111111111, \
      //
      0, -1.0, 0, 0, 1.0,	0, -2.0, -3.0, 0, -2.0, 3.0, 0, 0, 0, 1.0, 0, -1.0, \
      0, -2.0, 0, 1.0, 0, 1.0, 0, \
      0, 0, 1.0, 0.24, 0.06666666666666665,  0.8888888888888890,	\
      //
      0, 0.25, 0, 0, -0.75, 1.0, 1.0, 2.25, 1.0, 1.0, \
      -0.75, 0, 0, -0.25, 0, 0, 0.75, 1.0, 1.0, \
      0.75, 0, 0, -0.25, 0, 0, 0.25, 0, -0.12, -0.08333333333333331, 0.2222222222222222 \
  };

  //D3: Correct multiset of second order partials in (F,P,Dk) format. D3 cardinality = 4 for 2D
  const scalar_type basisD3[] = {
      0, -1.5, -1.5, 0, 0, -1.5, 0.5, 0, 0, 0.5,
      0.5, 0, 0, 0.5, -1.5, 0, 0, -1.5, -0.5, 0,
      0, -0.5, 0.5, 0, 0, 0.5, -0.5, 0, 0, -0.5,
      -1.5, 0, 0, -0.5, -0.5, 0, 0, -1.1, -0.1666666666666667, 0,
      //
      0, 3.0, 2.0, 0, 0, 3.0, -2.0, 0, 0, -1.0,
      -2.0, 0, 0, -1.0, 2.0, 0, 0, 3.0, 0, 0,
      0, 1.0, -2.0, 0, 0,	-1.0, 0, 0, 0, 1.0,
      2.0, 0, 0, 1.0, 0, 0, 0, 2.2, -0.6666666666666665, 0,
      //
      0, -1.5, -0.5, 0, 0, -1.5, 1.5, 0, 0, 0.5,
      1.5, 0, 0, 0.5, -0.5, 0, 0, -1.5, 0.5, 0,
      0, -0.5, 1.5, 0, 0, 0.5, 0.5, 0, 0, -0.5,
      -0.5, 0, 0, -0.5, 0.5, 0, 0, -1.1, 0.8333333333333333, 0,
      //
      0, 2.0, 3.0, 0, 0, 2.0, -1.0, 0, 0, -2.0,
      -1.0, 0, 0, -2.0, 3.0, 0, 0, 2.0,  1.0, 0,
      0, 0, -1.0, 0, 0, -2.0, 1.0, 0, 0, 0,
      3.0, 0, 0, 0, 1.0, 0, 0, 1.2, 0.3333333333333334, 0,
      //
      0, -4.0, -4.0, 0, 0, -4.0, 4.0, 0, 0, 4.0,
      4.0, 0, 0, 4.0, -4.0, 0, 0,	-4.0, 0, 0,
      0, 0, 4.0, 0, 0, 4.0, 0, 0, 0, 0,
      -4.0, 0, 0, 0, 0, 0, 0, -2.40, 1.333333333333333, 0,
      //
      0, 2.0, 1.0, 0, 0, 2.0, -3.0, 0, 0, -2.0,
      -3.0, 0, 0, -2.0, 1.0, 0, 0, 2.0, -1.0, 0,
      0, 0, -3.0, 0, 0, -2.0, -1.0, 0, 0, 0,
      1.0, 0, 0, 0, -1.0, 0, 0, 1.2, -1.666666666666667, 0 ,
      //
      0, -0.5, -1.5, 0, 0, -0.5, 0.5, 0, 0, 1.5,
      0.5, 0, 0, 1.5, -1.5, 0, 0, -0.5, -0.5, 0,
      0, 0.5, 0.5, 0, 0, 1.5, -0.5, 0, 0, 0.5,
      -1.5, 0, 0, 0.5, -0.5, 0,  0, -0.09999999999999998, -0.1666666666666667, 0,
      //
      0, 1.0, 2.0, 0, 0, 1.0, -2.0, 0, 0,	-3.0,
      -2.0, 0, 0, -3.0, 2.0, 0, 0, 1.0, 0, 0,
      0, -1.0, -2.0, 0, 0, -3.0, 0,  0, 0, -1.0,
      2.0, 0, 0, -1.0, 0, 0, 0, 0.2, -0.6666666666666665, 0,
      //
      0, -0.5, -0.5, 0, 0, -0.5, 1.5, 0, 0, 1.5,
      1.5, 0, 0, 1.5, -0.5, 0, 0, -0.5, 0.5, 0,
      0, 0.5, 1.5, 0, 0, 1.5, 0.5, 0, 0, 0.5,
      -0.5, 0, 0, 0.5, 0.5, 0, 0, -0.09999999999999998, 0.8333333333333333, 0
  };
  //D4: Correct multiset of second order partials in (F,P,Dk) format. D4 cardinality = 5 for 2D
  const scalar_type basisD4[] = {
      0, 0, 1.0, 0, 0, 0, 0, 1.0, 0, 0,
      0, 0, 1.0, 0, 0, 0, 0, 1.0, 0, 0,
      0, 0, 1.0, 0, 0, 0, 0, 1.0, 0, 0,
      0, 0, 1.0, 0, 0, 0, 0, 1.0, 0, 0,
      0, 0, 1.0, 0, 0, 0, 0, 1.0, 0, 0,
      //
      0, 0, -2.0, 0, 0, 0, 0, -2.0, 0, 0,
      0, 0, -2.0, 0, 0, 0, 0, -2.0, 0, 0,
      0, 0, -2.0, 0, 0, 0, 0, -2.0, 0, 0,
      0, 0, -2.0, 0, 0, 0, 0, -2.0, 0, 0,
      0, 0, -2.0, 0, 0, 0, 0, -2.0, 0, 0,
      //
      0, 0, 1.0, 0, 0, 0, 0, 1.0, 0, 0,
      0, 0, 1.0, 0, 0, 0, 0, 1.0, 0, 0,
      0, 0, 1.0, 0, 0, 0, 0, 1.0, 0, 0,
      0, 0, 1.0, 0, 0, 0, 0, 1.0, 0, 0,
      0, 0, 1.0, 0, 0, 0, 0, 1.0, 0, 0,
      //
      0, 0, -2.0, 0, 0, 0, 0, -2.0, 0, 0,
      0, 0, -2.0, 0, 0, 0, 0, -2.0, 0, 0,
      0, 0, -2.0, 0, 0, 0, 0, -2.0, 0, 0,
      0, 0, -2.0, 0, 0, 0, 0, -2.0, 0, 0,
      0, 0, -2.0, 0, 0, 0, 0, -2.0, 0, 0,
      //
      0, 0, 4.0, 0, 0, 0, 0, 4.0, 0, 0,
      0, 0, 4.0, 0, 0, 0, 0, 4.0, 0, 0,
      0, 0, 4.0, 0, 0, 0, 0, 4.0, 0, 0,
      0, 0, 4.0, 0, 0, 0, 0, 4.0, 0, 0,
      0, 0, 4.0, 0, 0, 0, 0, 4.0, 0, 0,
      //
      0, 0, -2.0, 0, 0, 0, 0, -2.0, 0, 0,
      0, 0, -2.0, 0, 0, 0, 0, -2.0, 0, 0,
      0, 0, -2.0, 0, 0, 0, 0, -2.0, 0, 0,
      0, 0, -2.0, 0, 0, 0, 0, -2.0, 0, 0,
      0, 0, -2.0, 0, 0, 0, 0, -2.0, 0, 0,
      //
      0, 0, 1.0, 0, 0, 0, 0, 1.0, 0, 0,
      0, 0, 1.0, 0, 0, 0, 0, 1.0, 0, 0,
      0, 0, 1.0, 0, 0, 0, 0, 1.0, 0, 0,
      0, 0, 1.0, 0, 0, 0, 0, 1.0, 0, 0,
      0, 0, 1.0, 0, 0, 0, 0, 1.0, 0, 0,
      //
      0, 0, -2.0, 0, 0, 0, 0, -2.0, 0, 0,
      0, 0, -2.0, 0, 0, 0, 0, -2.0, 0, 0,
      0, 0, -2.0, 0, 0, 0, 0, -2.0, 0, 0,
      0, 0, -2.0, 0, 0, 0, 0, -2.0, 0, 0,
      0, 0, -2.0, 0, 0, 0, 0, -2.0, 0, 0,
      //
      0, 0, 1.0, 0, 0, 0, 0, 1.0, 0, 0,
      0, 0, 1.0, 0, 0, 0, 0, 1.0, 0, 0,
      0, 0, 1.0, 0, 0, 0, 0, 1.0, 0, 0,
      0, 0, 1.0, 0, 0, 0, 0, 1.0, 0, 0,
      0, 0, 1.0, 0, 0, 0, 0, 1.0, 0, 0
  };

  try {
    constexpr ordinal_type order = 2;
    if(order <= maxOrder) {
      QuadBasisType quadBasis(order);

      DynRankViewHostScalarValueType ConstructWithLabel(quadNodesHost, 10, 2);
      DynRankViewPointValueType ConstructWithLabelPointView(quadNodes, 10, 2);

      quadNodesHost(0,0) = -1.0;  quadNodesHost(0,1) = -1.0;
      quadNodesHost(1,0) =  1.0;  quadNodesHost(1,1) = -1.0;
      quadNodesHost(2,0) =  1.0;  quadNodesHost(2,1) =  1.0;
      quadNodesHost(3,0) = -1.0;  quadNodesHost(3,1) =  1.0;
      // edge midpoints
      quadNodesHost(4,0) =  0.0;  quadNodesHost(4,1) = -1.0;
      quadNodesHost(5,0) =  1.0;  quadNodesHost(5,1) =  0.0;
      quadNodesHost(6,0) =  0.0;  quadNodesHost(6,1) =  1.0;
      quadNodesHost(7,0) = -1.0;  quadNodesHost(7,1) =  0.0;
      // center & random point
      quadNodesHost(8,0) =  0.0;  quadNodesHost(8,1) =  0.0;
      quadNodesHost(9,0) =1./3.;  quadNodesHost(9,1) =-3./5.;


      auto quadNodes_scalar = Kokkos::create_mirror_view(typename DeviceType::memory_space(), quadNodesHost);
      Kokkos::deep_copy(quadNodes_scalar, quadNodesHost);
      RealSpaceTools<DeviceType>::clone(quadNodes, quadNodes_scalar);
      // auto quadNodesHost2 = Kokkos::create_mirror_view(quadNodes);
      // Kokkos::deep_copy(quadNodesHost2, quadNodes);

      // for(int i=0; i< quadNodesHost.extent(0); ++i)
      // std::cout <<  quadNodesHost(i,0) << " "<<  quadNodesHost(i,1) <<std::endl;
      // for(int i=0; i< quadNodesHost2.extent(0); ++i)
      // std::cout <<  quadNodesHost2(i,0) << " "<<  quadNodesHost2(i,1) <<std::endl;

      // Dimensions for the output arrays:
      const ordinal_type numFields = quadBasis.getCardinality();
      const ordinal_type numPoints = quadNodes.extent(0);
      const ordinal_type spaceDim  = quadBasis.getBaseCellTopology().getDimension();
      const ordinal_type D2Cardin  = getDkCardinality(OPERATOR_D2, spaceDim);
      const ordinal_type D3Cardin  = getDkCardinality(OPERATOR_D3, spaceDim);
      const ordinal_type D4Cardin  = getDkCardinality(OPERATOR_D4, spaceDim);

      *outStream << " -- Testing OPERATOR_VALUE \n";
      {
        // Check VALUE of basis functions: resize vals to rank-2 container:
        DynRankViewOutValueType ConstructWithLabelOutView(vals, numFields, numPoints);
        quadBasis.getValues(space, vals, quadNodes, OPERATOR_VALUE);
        auto vals_host = Kokkos::create_mirror_view(vals);
        Kokkos::deep_copy(vals_host, vals);
        for (ordinal_type i = 0; i < numFields; ++i) {
          for (ordinal_type j = 0; j < numPoints; ++j) {

            // Compute offset for (F,P) container
            ordinal_type l =  j + i * numPoints;
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
      *outStream << " -- Testing OPERATOR_GRAD \n";
      {
        // Check GRAD of basis function: resize vals to rank-3 container
        DynRankViewOutValueType ConstructWithLabelOutView(vals, numFields, numPoints, spaceDim);
        quadBasis.getValues(space, vals, quadNodes, OPERATOR_GRAD);
        auto vals_host = Kokkos::create_mirror_view(vals);
        Kokkos::deep_copy(vals_host, vals);
        for (ordinal_type i = 0; i < numFields; ++i) {
          for (ordinal_type j = 0; j < numPoints; ++j) {
            for (ordinal_type k = 0; k < spaceDim; ++k) {

              // basisGrads is (F,P,D), compute offset:
              ordinal_type const l = k + j * spaceDim + i * spaceDim * numPoints;
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
      *outStream << " -- Testing OPERATOR_D1 \n";
      {
        // Check D1 of basis function (do not resize vals because it has the correct size: D1 = GRAD)
        DynRankViewOutValueType ConstructWithLabelOutView(vals, numFields, numPoints, spaceDim);
        quadBasis.getValues(space, vals, quadNodes, OPERATOR_D1);
        auto vals_host = Kokkos::create_mirror_view(vals);
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
      *outStream << " -- Testing OPERATOR_CURL \n";
      {
        // Check CURL of basis function: resize vals just for illustration!
        DynRankViewOutValueType ConstructWithLabelOutView(vals, numFields, numPoints, spaceDim);
        quadBasis.getValues(space, vals, quadNodes, OPERATOR_CURL);
        auto vals_host = Kokkos::create_mirror_view(vals);
        Kokkos::deep_copy(vals_host, vals);
        for (ordinal_type i = 0; i < numFields; ++i) {
          for (ordinal_type j = 0; j < numPoints; ++j) {
            // We will use "rotated" basisGrads to check CURL: get offsets to extract (u_y, -u_x)
            const ordinal_type curl_0 = 1 + j * spaceDim + i * spaceDim * numPoints;               // position of y-derivat$
            const ordinal_type curl_1 = 0 + j * spaceDim + i * spaceDim * numPoints;               // position of x-derivat$

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
      *outStream << " -- Testing OPERATOR_D2 \n";
      {
        // Check D2 of basis function
        DynRankViewOutValueType ConstructWithLabelOutView(vals, numFields, numPoints, D2Cardin);
        quadBasis.getValues(space, vals, quadNodes, OPERATOR_D2);
        auto vals_host = Kokkos::create_mirror_view( vals);
        Kokkos::deep_copy(vals_host, vals);
        for (ordinal_type i = 0; i < numFields; ++i) {
          for (ordinal_type j = 0; j < numPoints; ++j) {
            for (ordinal_type k = 0; k < 1; ++k) { //D2Cardin
              const auto val = get_scalar_value(vals_host(i,j,k)) ;
              // basisD2 is (F,P,Dk), compute offset:
              const ordinal_type l = k + j * D2Cardin + i * D2Cardin * numPoints;
              if (std::isnan(val) || std::abs(val - basisD2[l]) > tol) {
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
      *outStream << " -- Testing OPERATOR_D3 \n";
      {
        // Check D3 of basis function
        DynRankViewOutValueType ConstructWithLabelOutView(vals, numFields, numPoints, D3Cardin);
        quadBasis.getValues(space, vals, quadNodes, OPERATOR_D3);
        auto vals_host = Kokkos::create_mirror_view(vals);
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
                                 << " but reference D3 component: " << basisD2[l] << "\n"; //not D3?
              }
            }
          }
        }
      }
      *outStream << " -- Testing OPERATOR_D4 \n";
      {
        // Check D4 of basis function
        DynRankViewOutValueType ConstructWithLabelOutView(vals, numFields, numPoints, D4Cardin);
        quadBasis.getValues(space, vals, quadNodes, OPERATOR_D4);
        auto vals_host = Kokkos::create_mirror_view(vals);
        Kokkos::deep_copy(vals_host, vals);

        for (ordinal_type i = 0; i < numFields; ++i) {
          for (ordinal_type j = 0; j < numPoints; ++j) {
            for (ordinal_type k = 0; k < D4Cardin; ++k) {

              // basisD4 is (F,P,Dk), compute offset:
              ordinal_type l = k + j * D4Cardin + i * D4Cardin * numPoints;
              if (std::abs(vals_host(i,j,k) - basisD4[l]) > tol) {
                errorFlag++;
                *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";

                // Output the multi-index of the value where the error is:
                *outStream << " At multi-index { ";
                *outStream << i << " ";*outStream << j << " ";*outStream << k << " ";
                *outStream << "}  computed D4 component: " << vals_host(i,j,k)
                                 << " but reference D4 component: " << basisD2[l] << "\n"; //not D4?
              }
            }
          }
        }
      }
      // Check all higher derivatives - must be zero.
      if (0) {
        const EOperator ops[] = { OPERATOR_D5,
            OPERATOR_D6,
            OPERATOR_D7,
            OPERATOR_D8,
            OPERATOR_D9,
            OPERATOR_D10,
            OPERATOR_MAX };
        for (auto h=0;ops[h]!=OPERATOR_MAX;++h) {
          //    for(EOperator op = OPERATOR_D5; op < OPERATOR_MAX; op++) {
            const auto op = ops[h];
            // The last dimension is the number of kth derivatives and needs to be resized for every Dk
            const ordinal_type DkCardin  = getDkCardinality(op, spaceDim);
            DynRankViewOutValueType ConstructWithLabelOutView(vals, numFields, numPoints, DkCardin);
            quadBasis.getValues(space, vals, quadNodes, op);
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
  << "| TEST 4: Function Space is Correct                                           |\n"
  << "===============================================================================\n";
  
  try {
    const auto order = std::min(5, maxOrder);
    QuadBasisType quadBasis(order);
    
    const EFunctionSpace fs = quadBasis.getFunctionSpace();
    
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
