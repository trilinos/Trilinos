// @HEADER
// *****************************************************************************
//                           Intrepid2 Package
//
// Copyright 2007 NTESS and the Intrepid2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/** \file test_01.cpp
    \brief  Unit tests for the Intrepid2::HCURL_TET_In_FEM class.
    \author Created by P. Bochev, D. Ridzal, K. Peterson, Kyungjoo Kim and Mauro Perego
 */


#include "Intrepid2_config.h"

#ifdef HAVE_INTREPID2_DEBUG
#define INTREPID2_TEST_FOR_DEBUG_ABORT_OVERRIDE_TO_CONTINUE
#endif

#include "Intrepid2_Types.hpp"
#include "Intrepid2_Utils.hpp"

#include "Intrepid2_PointTools.hpp"
#include "Intrepid2_HCURL_TET_In_FEM.hpp"


#include "Teuchos_oblackholestream.hpp"
#include "Teuchos_RCP.hpp"

#include "packages/intrepid2/unit-test/Discretization/Basis/Macros.hpp"
#include "packages/intrepid2/unit-test/Discretization/Basis/Setup.hpp"

namespace Intrepid2 {

namespace Test {

using HostSpaceType = Kokkos::DefaultHostExecutionSpace;

template<typename OutValueType, typename PointValueType, typename DeviceType>
int HCURL_TET_In_FEM_Test01(const bool verbose) {

  //! Create an execution space instance.
  const auto space = Kokkos::Experimental::partition_space(typename DeviceType::execution_space {}, 1)[0];

  //! Setup test output stream.
  Teuchos::RCP<std::ostream> outStream = setup_output_stream<DeviceType>(
    verbose, "HCURL_TET_In_FEM", {
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

  // for virtual function, value and point types are declared in the class
  typedef OutValueType outputValueType;
  typedef PointValueType pointValueType;

  typedef Basis_HCURL_TET_In_FEM<DeviceType,outputValueType,pointValueType> TetBasisType;

  constexpr ordinal_type maxOrder = Parameters::MaxOrder;
  constexpr ordinal_type dim = 3;

  try {

    *outStream
    << "\n"
    << "===============================================================================\n"
    << "| TEST 1: Testing OPERATOR_VALUE (Kronecker property using getDofCoeffs)      |\n"
    << "===============================================================================\n";


    const ordinal_type order = std::min(3, maxOrder);
    TetBasisType tetBasis(order, POINTTYPE_WARPBLEND);

    const ordinal_type cardinality = tetBasis.getCardinality();
    DynRankViewScalarValueType ConstructWithLabel(dofCoords_scalar, cardinality , dim);
    tetBasis.getDofCoords(dofCoords_scalar);

    DynRankViewPointValueType ConstructWithLabelPointView(dofCoords, cardinality , dim);
    RealSpaceTools<DeviceType>::clone(dofCoords, dofCoords_scalar);

    DynRankViewScalarValueType ConstructWithLabel(dofCoeffs, cardinality , dim);
    tetBasis.getDofCoeffs(dofCoeffs);

    DynRankViewOutValueType ConstructWithLabelOutView(basisAtDofCoords, cardinality , cardinality, dim);
    tetBasis.getValues(space, basisAtDofCoords, dofCoords, OPERATOR_VALUE);

    auto h_basisAtDofCoords = Kokkos::create_mirror_view(basisAtDofCoords);
    Kokkos::deep_copy(h_basisAtDofCoords, basisAtDofCoords);

    auto h_dofCoeffs = Kokkos::create_mirror_view(dofCoeffs);
    Kokkos::deep_copy(h_dofCoeffs, dofCoeffs);

    // test for Kronecker property
    for (int i=0;i<cardinality;i++) {
      for (int j=0;j<cardinality;j++) {
        outputValueType dofValue = 0;
        for (ordinal_type k=0;k<dim; k++)
          dofValue += h_basisAtDofCoords(i,j,k)*h_dofCoeffs(j,k);

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

  try {

    *outStream
    << "\n"
    << "=======================================================================================\n"
    << "| TEST 2: Testing OPERATOR_VALUE (Kronecker property, reconstructing dofs using tags) |\n"
    << "=======================================================================================\n";


    const ordinal_type order = std::min(3, maxOrder);
    TetBasisType tetBasis(order, POINTTYPE_WARPBLEND);
    shards::CellTopology tet_4(shards::getCellTopologyData<shards::Tetrahedron<4> >());


    const ordinal_type cardinality = tetBasis.getCardinality();
    DynRankViewScalarValueType ConstructWithLabel(dofCoords_scalar, cardinality , dim);
    tetBasis.getDofCoords(dofCoords_scalar);

    DynRankViewPointValueType ConstructWithLabelPointView(dofCoords, cardinality , dim);
    RealSpaceTools<DeviceType>::clone(dofCoords, dofCoords_scalar);

    DynRankViewOutValueType ConstructWithLabelOutView(basisAtDofCoords, cardinality , cardinality, dim);
    tetBasis.getValues(space, basisAtDofCoords, dofCoords, OPERATOR_VALUE);

    auto h_basisAtDofCoords = Kokkos::create_mirror_view(basisAtDofCoords);
    Kokkos::deep_copy(h_basisAtDofCoords, basisAtDofCoords);

    // Dimensions for the output arrays:
    const ordinal_type numFields = tetBasis.getCardinality();

    //Normals at each edge
    DynRankViewHostScalarValueType ConstructWithLabel(tangents, numFields,dim); // normals at each point basis point
    DynRankViewHostScalarValueType ConstructWithLabel(edgeTan, dim );
    DynRankViewHostScalarValueType ConstructWithLabel(faceTan1, dim );
    DynRankViewHostScalarValueType ConstructWithLabel(faceTan2, dim );

    const auto allTags = tetBasis.getAllDofTags();

    // test for Kronecker property
    for (int i=0;i<numFields;i++) {
      for (int j=0;j<numFields;j++) {
        outputValueType dofValue = 0;
        if(allTags(j,0) == 1) { //edge
          auto edgeId = allTags(j,1);
          CellTools<Kokkos::HostSpace::execution_space>::getReferenceEdgeTangent( edgeTan ,
              edgeId ,
              tet_4 );

          for (ordinal_type k=0;k<dim; k++)
            dofValue += h_basisAtDofCoords(i,j,k)*edgeTan(k);
        }
        else if(allTags(j,0) == 2) { //face
          auto faceId = allTags(j,1);
          CellTools<Kokkos::HostSpace::execution_space>::getReferenceFaceTangents( faceTan1 ,
              faceTan2,
              faceId ,
              tet_4  );

          auto faceDofId =  allTags(j,2);
          dofValue = 0;
          if(faceDofId%2 ==0 )
            for (ordinal_type k=0;k<dim; k++)
              dofValue += h_basisAtDofCoords(i,j,k)*faceTan1(k);
          else
            for (ordinal_type k=0;k<dim; k++)
              dofValue += h_basisAtDofCoords(i,j,k)*faceTan2(k);
        }
        else { //elem
          auto elemDofId =  allTags(j,2);
          dofValue = h_basisAtDofCoords(i,j,elemDofId%dim);
        }

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

  // Intrepid2 basis have been redefined and they are no longer
  // equivalent to FIAT basis. However they are proportional to the FIAT basis
  // with a scaling that depends on the geometric entity (points, edge, face, cell)
  // associated to the basis
  scalar_type scaling_factor[4] = {0,2,1,1};

  try {

    *outStream
    << "\n"
    << "===============================================================================\n"
    << "| TEST 3: Testing OPERATOR_VALUE (FIAT Values)                                |\n"
    << "===============================================================================\n";


    const ordinal_type order = 1;
    TetBasisType tetBasis(order, POINTTYPE_EQUISPACED);

    shards::CellTopology tet_4(shards::getCellTopologyData<shards::Tetrahedron<4> >());
    const ordinal_type np_lattice = PointTools::getLatticeSize(tet_4, order,0);
    const ordinal_type cardinality = tetBasis.getCardinality();

    //Need to use Scalar type for lattice because PointTools dont's work with FAD types
    DynRankViewScalarValueType ConstructWithLabel(lattice_scalar, np_lattice , dim);
    PointTools::getLattice(lattice_scalar, tet_4, order, 0, POINTTYPE_EQUISPACED);
    DynRankViewPointValueType ConstructWithLabelPointView(lattice, np_lattice , dim);
    RealSpaceTools<DeviceType>::clone(lattice,lattice_scalar);

    DynRankViewOutValueType ConstructWithLabelOutView(basisAtLattice, cardinality , np_lattice, dim);
    tetBasis.getValues(space, basisAtLattice, lattice, OPERATOR_VALUE);

    auto h_basisAtLattice = Kokkos::create_mirror_view(basisAtLattice);
    Kokkos::deep_copy(h_basisAtLattice, basisAtLattice);

    // Dimensions for the output arrays:
    const ordinal_type numFields = tetBasis.getCardinality();

    const scalar_type fiat_vals[] = {
        1.000000000000001e+00, -2.498001805406602e-16, -1.665334536937735e-16,
        9.999999999999998e-01, 1.000000000000000e+00, 1.000000000000000e+00,
        5.828670879282072e-16, 1.110223024625157e-16, 2.498001805406602e-16,
        7.771561172376096e-16, 8.326672684688674e-17, 1.110223024625157e-16,
        2.081668171172169e-16, -2.914335439641036e-16, 1.280865063236792e-16,
        -3.191891195797325e-16, 1.000000000000000e+00, -4.293998586504916e-17,
        -9.999999999999994e-01, 2.081668171172169e-16, 2.400576428367544e-16,
        2.220446049250313e-16, -5.551115123125783e-17, 1.084013877651281e-16,
        3.469446951953614e-16, -1.000000000000000e+00, 1.387778780781446e-16,
        -1.804112415015879e-16, 1.942890293094024e-16, -1.387778780781446e-16,
        -9.999999999999993e-01, -9.999999999999996e-01, -9.999999999999998e-01,
        5.551115123125783e-17, -2.220446049250313e-16, -8.326672684688674e-17,
        -2.220446049250313e-16, -5.551115123125783e-17, 9.999999999999999e-01,
        1.665334536937735e-16, 1.110223024625157e-16, -6.383782391594650e-16,
        1.110223024625157e-16, 1.110223024625157e-16, -1.110223024625157e-16,
        9.999999999999990e-01, 9.999999999999994e-01, 9.999999999999996e-01,
        1.387778780781446e-16, -2.496931404305374e-17, -1.665334536937735e-16,
        -2.498001805406602e-16, -2.149987498083074e-16, 1.000000000000000e+00,
        8.326672684688674e-17, -3.769887250591415e-17, 8.326672684688674e-17,
        -9.999999999999994e-01, 1.556977698723022e-16, 2.220446049250313e-16,
        -9.422703950001342e-18, 1.665334536937735e-16, -2.359223927328458e-16,
        -9.422703950001268e-18, -8.326672684688674e-17, 1.387778780781446e-17,
        -7.525083148581445e-17, 2.775557561562891e-17, 1.000000000000000e+00,
        2.789513560273035e-16, -9.999999999999998e-01, -5.551115123125783e-17
    };

    const auto allTags = tetBasis.getAllDofTags();
    ordinal_type cur=0;
    for (ordinal_type i=0;i<numFields;i++) {
      auto scaling = scaling_factor[allTags(i,0)];
      for (ordinal_type j=0;j<np_lattice;j++) {
        for (ordinal_type k=0;k<dim; k++) {
          if (std::fabs( h_basisAtLattice(i,j,k) - scaling*fiat_vals[cur] ) > tol ) {
            errorFlag++;
            *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";

            // Output the multi-index of the value where the error is:
            *outStream << " At multi-index { ";
            *outStream << i << " " << j << " " << k;
            *outStream << "}  computed value: " <<  h_basisAtLattice(i,j,k)
                              << " but correct value: " << scaling*fiat_vals[cur] << "\n";
            *outStream << "Difference: " << std::fabs(  h_basisAtLattice(i,j,k) - scaling*fiat_vals[cur] ) << "\n";
          }
          cur++;
        }
      }
    }

  } catch (std::exception &err) {
    *outStream << "UNEXPECTED ERROR !!! ----------------------------------------------------------\n";
    *outStream << err.what() << '\n';
    *outStream << "-------------------------------------------------------------------------------" << "\n\n";
    errorFlag = -1000;
  };

  try {

    *outStream
    << "\n"
    << "===============================================================================\n"
    << "| TEST 4: Testing Operator CURL                                               |\n"
    << "===============================================================================\n";


    const ordinal_type order = 1;
    TetBasisType tetBasis(order, POINTTYPE_EQUISPACED);

    shards::CellTopology tet_4(shards::getCellTopologyData<shards::Tetrahedron<4> >());
    const ordinal_type np_lattice = PointTools::getLatticeSize(tet_4, order,0);
    const ordinal_type cardinality = tetBasis.getCardinality();

    DynRankViewScalarValueType ConstructWithLabel(lattice_scalar, np_lattice , dim);
    PointTools::getLattice(lattice_scalar, tet_4, order, 0, POINTTYPE_EQUISPACED);
    DynRankViewPointValueType ConstructWithLabelPointView(lattice, np_lattice , dim);
    RealSpaceTools<DeviceType>::clone(lattice,lattice_scalar);

    DynRankViewOutValueType ConstructWithLabelOutView(curlBasisAtLattice, cardinality , np_lattice, dim);
    tetBasis.getValues(space, curlBasisAtLattice, lattice, OPERATOR_CURL);

    auto h_curlBasisAtLattice = Kokkos::create_mirror_view(curlBasisAtLattice);
    Kokkos::deep_copy(h_curlBasisAtLattice, curlBasisAtLattice);

    // Dimensions for the output arrays:
    const ordinal_type numFields = tetBasis.getCardinality();

    const scalar_type fiat_curls[] = {
        -5.551115123125783e-16, -2.000000000000000e+00, 2.000000000000000e+00,
        -5.551115123125783e-16, -2.000000000000000e+00, 2.000000000000000e+00,
        -5.551115123125783e-16, -2.000000000000000e+00, 2.000000000000000e+00,
        -5.551115123125783e-16, -2.000000000000000e+00, 2.000000000000000e+00,
        -4.440892098500626e-16, -2.775557561562891e-16, 2.000000000000000e+00,
        -4.440892098500626e-16, -2.775557561562891e-16, 2.000000000000000e+00,
        -4.440892098500626e-16, -2.775557561562891e-16, 2.000000000000000e+00,
        -4.440892098500626e-16, -2.775557561562891e-16, 2.000000000000000e+00,
        -2.000000000000000e+00, 5.551115123125783e-17, 2.000000000000000e+00,
        -2.000000000000000e+00, 5.551115123125783e-17, 2.000000000000000e+00,
        -2.000000000000000e+00, 5.551115123125783e-17, 2.000000000000000e+00,
        -2.000000000000000e+00, 5.551115123125783e-17, 2.000000000000000e+00,
        -2.000000000000000e+00, 2.000000000000000e+00, 9.861075762086680e-17,
        -2.000000000000000e+00, 2.000000000000000e+00, 9.861075762086680e-17,
        -2.000000000000000e+00, 2.000000000000000e+00, 9.861075762086680e-17,
        -2.000000000000000e+00, 2.000000000000000e+00, 9.861075762086680e-17,
        -2.775557561562891e-17, -2.000000000000000e+00, 4.287451790760826e-16,
        -2.775557561562891e-17, -2.000000000000000e+00, 4.287451790760826e-16,
        -2.775557561562891e-17, -2.000000000000000e+00, 4.287451790760826e-16,
        -2.775557561562891e-17, -2.000000000000000e+00, 4.287451790760826e-16,
        2.000000000000000e+00, -2.185751579730777e-16, 1.526556658859590e-16,
        2.000000000000000e+00, -2.185751579730777e-16, 1.526556658859590e-16,
        2.000000000000000e+00, -2.185751579730777e-16, 1.526556658859590e-16,
        2.000000000000000e+00, -2.185751579730777e-16, 1.526556658859590e-16
    };

    ordinal_type cur=0;
    const auto allTags = tetBasis.getAllDofTags();
    for (ordinal_type i=0;i<numFields;i++) {
      auto scaling = scaling_factor[allTags(i,0)];
      for (ordinal_type j=0;j<np_lattice;j++) {
        for (ordinal_type k=0;k<dim; k++) {
          if (std::abs( h_curlBasisAtLattice(i,j,k) - scaling*fiat_curls[cur] ) > tol ) {
            errorFlag++;
            *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";

            // Output the multi-index of the value where the error is:
            *outStream << " At multi-index { ";
            *outStream << i << " " << j;
            *outStream << "}  computed value: " <<  h_curlBasisAtLattice(i,j,k)
                                << " but correct value: " << scaling*fiat_curls[cur] << "\n";
            *outStream << "Difference: " << std::fabs(  h_curlBasisAtLattice(i,j,k) - scaling*fiat_curls[cur] ) << "\n";
          }
          cur++;
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
    TetBasisType tetBasis(order, POINTTYPE_WARPBLEND);

    const EFunctionSpace fs = tetBasis.getFunctionSpace();

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
