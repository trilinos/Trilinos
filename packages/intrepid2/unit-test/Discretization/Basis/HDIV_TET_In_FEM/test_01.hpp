// @HEADER
// *****************************************************************************
//                           Intrepid2 Package
//
// Copyright 2007 NTESS and the Intrepid2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER


/** \file
    \brief  Unit tests for the Intrepid2::HDIV_TET_In_FEM class. 
    \author Created by R. Kirby
 */


#include "Intrepid2_config.h"

#ifdef HAVE_INTREPID2_DEBUG
#define INTREPID2_TEST_FOR_DEBUG_ABORT_OVERRIDE_TO_CONTINUE
#endif

#include "Intrepid2_Types.hpp"
#include "Intrepid2_Utils.hpp"

#include "Intrepid2_PointTools.hpp"
#include "Intrepid2_HDIV_TET_In_FEM.hpp"

#include "Teuchos_oblackholestream.hpp"
#include "Teuchos_RCP.hpp"

#include "packages/intrepid2/unit-test/Discretization/Basis/Macros.hpp"

namespace Intrepid2 {

namespace Test {

template<typename OutValueType, typename PointValueType, typename DeviceType>
int HDIV_TET_In_FEM_Test01(const bool verbose) {

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
  <<"\n"
  << "===============================================================================\n"
  << "|                                                                             |\n"
  << "|                 Unit Test HDIV_TET_In_FEM                                   |\n"
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

  // for virtual function, value and point types are declared in the class
  typedef OutValueType outputValueType;
  typedef PointValueType pointValueType;

  typedef Basis_HDIV_TET_In_FEM<DeviceType,outputValueType,pointValueType> TetBasisType;

  constexpr ordinal_type maxOrder = Parameters::MaxOrder;
  constexpr ordinal_type dim = 3;

  try {

    *outStream
    << "\n"
    << "===============================================================================\n"
    << "| TEST 1: Testing OPERATOR_VALUE (Kronecker property using getDofCoeffs)      |\n"
    << "===============================================================================\n";

    const ordinal_type order = std::min(4, maxOrder);
    TetBasisType tetBasis(order, POINTTYPE_EQUISPACED);

    const ordinal_type cardinality = tetBasis.getCardinality();
    DynRankViewScalarValueType ConstructWithLabel(dofCoords_scalar, cardinality , dim);
    tetBasis.getDofCoords(dofCoords_scalar);

    DynRankViewPointValueType ConstructWithLabelPointView(dofCoords, cardinality , dim);
    RealSpaceTools<DeviceType>::clone(dofCoords, dofCoords_scalar);

    DynRankViewScalarValueType ConstructWithLabel(dofCoeffs, cardinality , dim);
    tetBasis.getDofCoeffs(dofCoeffs);

    DynRankViewOutValueType ConstructWithLabelOutView(basisAtDofCoords, cardinality , cardinality, dim);
    tetBasis.getValues(basisAtDofCoords, dofCoords, OPERATOR_VALUE);

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

    const ordinal_type order = std::min(4, maxOrder);
    TetBasisType tetBasis(order, POINTTYPE_EQUISPACED);
    const ordinal_type cardinality = tetBasis.getCardinality();
    DynRankViewScalarValueType ConstructWithLabel(dofCoords_scalar, cardinality , dim);
    tetBasis.getDofCoords(dofCoords_scalar);

    DynRankViewPointValueType ConstructWithLabelPointView(dofCoords, cardinality , dim);
    RealSpaceTools<DeviceType>::clone(dofCoords, dofCoords_scalar);

    DynRankViewOutValueType ConstructWithLabelOutView(basisAtDofCoords, cardinality , cardinality, dim);
    tetBasis.getValues(basisAtDofCoords, dofCoords, OPERATOR_VALUE);

    auto h_basisAtDofCoords = Kokkos::create_mirror_view(basisAtDofCoords);
    Kokkos::deep_copy(h_basisAtDofCoords, basisAtDofCoords);

    //Normals at each face
    DynRankViewHostScalarValueType ConstructWithLabel(normals, cardinality, dim); // normals at each point basis point
    shards::CellTopology tet_4(shards::getCellTopologyData<shards::Tetrahedron<4> >());
    for (int sideId = 0; sideId < 4; ++sideId) {
      auto normal = Kokkos::subview(normals, sideId, Kokkos::ALL());
    CellTools<Kokkos::HostSpace::execution_space>::getReferenceSideNormal( normal ,
        sideId ,
        tet_4 );
    }
    const auto allTags = tetBasis.getAllDofTags();

    for (int i=0;i<cardinality;i++) {
      for (int j=0;j<cardinality;j++) {
        outputValueType dofValue = 0;
        if(allTags(j,0) == dim-1) { //face
          auto faceId = allTags(j,1);
          for (ordinal_type k=0;k<dim; k++) {
            dofValue += h_basisAtDofCoords(i,j,k)*normals(faceId,k);
          }
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
  scalar_type scaling_factor[4] = {0,0,0.5,1};

  try {

    *outStream
    << "\n"
    << "===============================================================================\n"
    << "| TEST 3: Testing OPERATOR_VALUE (FIAT Values)                                |\n"
    << "===============================================================================\n";

    constexpr ordinal_type order = 2;
    if(order <= maxOrder) {
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
      tetBasis.getValues(basisAtLattice, lattice, OPERATOR_VALUE);

      auto h_basisAtLattice = Kokkos::create_mirror_view(basisAtLattice);
      Kokkos::deep_copy(h_basisAtLattice, basisAtLattice);

      const scalar_type fiat_vals[] = {
          0.000000000000000e+00, -6.000000000000000e+00, 0.000000000000000e+00,
          1.000000000000000e+00, -2.000000000000000e+00, 0.000000000000000e+00,
          -2.000000000000000e+00, 2.000000000000000e+00, 0.000000000000000e+00,
          0.000000000000000e+00, -1.000000000000000e+00, 0.000000000000000e+00,
          -1.000000000000000e+00, 1.000000000000000e+00, 0.000000000000000e+00,
          0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00,
          0.000000000000000e+00, -2.000000000000000e+00, 1.000000000000000e+00,
          -1.000000000000000e+00, 2.000000000000000e+00, -1.000000000000000e+00,
          0.000000000000000e+00, 1.000000000000000e+00, -1.000000000000000e+00,
          0.000000000000000e+00, 2.000000000000000e+00, -2.000000000000000e+00,
          0.000000000000000e+00, 2.000000000000000e+00, 0.000000000000000e+00,
          1.000000000000000e+00, -2.000000000000000e+00, 0.000000000000000e+00,
          6.000000000000000e+00, -6.000000000000000e+00, 0.000000000000000e+00,
          0.000000000000000e+00, 1.000000000000000e+00, 0.000000000000000e+00,
          1.000000000000000e+00, -1.000000000000000e+00, 0.000000000000000e+00,
          0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00,
          0.000000000000000e+00, 2.000000000000000e+00, -1.000000000000000e+00,
          1.000000000000000e+00, -2.000000000000000e+00, 1.000000000000000e+00,
          0.000000000000000e+00, 1.000000000000000e+00, -1.000000000000000e+00,
          0.000000000000000e+00, 2.000000000000000e+00, -2.000000000000000e+00,
          0.000000000000000e+00, 2.000000000000000e+00, 0.000000000000000e+00,
          -1.000000000000000e+00, 2.000000000000000e+00, 0.000000000000000e+00,
          -2.000000000000000e+00, 2.000000000000000e+00, 0.000000000000000e+00,
          0.000000000000000e+00, 1.000000000000000e+00, 0.000000000000000e+00,
          -1.000000000000000e+00, 1.000000000000000e+00, 0.000000000000000e+00,
          0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00,
          0.000000000000000e+00, -2.000000000000000e+00, 1.000000000000000e+00,
          1.000000000000000e+00, -2.000000000000000e+00, 1.000000000000000e+00,
          0.000000000000000e+00, -1.000000000000000e+00, 1.000000000000000e+00,
          0.000000000000000e+00, -6.000000000000000e+00, 6.000000000000000e+00,
          0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00,
          1.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00,
          6.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00,
          0.000000000000000e+00, -1.000000000000000e+00, 0.000000000000000e+00,
          1.000000000000000e+00, 1.000000000000000e+00, 0.000000000000000e+00,
          0.000000000000000e+00, -2.000000000000000e+00, 0.000000000000000e+00,
          0.000000000000000e+00, 0.000000000000000e+00, -1.000000000000000e+00,
          1.000000000000000e+00, 0.000000000000000e+00, 1.000000000000000e+00,
          0.000000000000000e+00, -1.000000000000000e+00, -1.000000000000000e+00,
          0.000000000000000e+00, 0.000000000000000e+00, -2.000000000000000e+00,
          0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00,
          -1.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00,
          -2.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00,
          0.000000000000000e+00, 1.000000000000000e+00, 0.000000000000000e+00,
          1.000000000000000e+00, 1.000000000000000e+00, 0.000000000000000e+00,
          0.000000000000000e+00, 6.000000000000000e+00, 0.000000000000000e+00,
          0.000000000000000e+00, 0.000000000000000e+00, -1.000000000000000e+00,
          -1.000000000000000e+00, 0.000000000000000e+00, -1.000000000000000e+00,
          0.000000000000000e+00, 1.000000000000000e+00, 1.000000000000000e+00,
          0.000000000000000e+00, 0.000000000000000e+00, -2.000000000000000e+00,
          0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00,
          -1.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00,
          -2.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00,
          0.000000000000000e+00, -1.000000000000000e+00, 0.000000000000000e+00,
          -1.000000000000000e+00, -1.000000000000000e+00, 0.000000000000000e+00,
          0.000000000000000e+00, -2.000000000000000e+00, 0.000000000000000e+00,
          0.000000000000000e+00, 0.000000000000000e+00, 1.000000000000000e+00,
          1.000000000000000e+00, 0.000000000000000e+00, 1.000000000000000e+00,
          0.000000000000000e+00, 1.000000000000000e+00, 1.000000000000000e+00,
          0.000000000000000e+00, 0.000000000000000e+00, 6.000000000000000e+00,
          -6.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00,
          -1.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00,
          0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00,
          -2.000000000000000e+00, 1.000000000000000e+00, 0.000000000000000e+00,
          1.000000000000000e+00, -1.000000000000000e+00, 0.000000000000000e+00,
          2.000000000000000e+00, -2.000000000000000e+00, 0.000000000000000e+00,
          -2.000000000000000e+00, 0.000000000000000e+00, 1.000000000000000e+00,
          1.000000000000000e+00, 0.000000000000000e+00, -1.000000000000000e+00,
          2.000000000000000e+00, -1.000000000000000e+00, -1.000000000000000e+00,
          2.000000000000000e+00, 0.000000000000000e+00, -2.000000000000000e+00,
          2.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00,
          1.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00,
          0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00,
          2.000000000000000e+00, -1.000000000000000e+00, 0.000000000000000e+00,
          1.000000000000000e+00, -1.000000000000000e+00, 0.000000000000000e+00,
          2.000000000000000e+00, -2.000000000000000e+00, 0.000000000000000e+00,
          -2.000000000000000e+00, 0.000000000000000e+00, 1.000000000000000e+00,
          -1.000000000000000e+00, 0.000000000000000e+00, 1.000000000000000e+00,
          -2.000000000000000e+00, 1.000000000000000e+00, 1.000000000000000e+00,
          -6.000000000000000e+00, 0.000000000000000e+00, 6.000000000000000e+00,
          2.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00,
          1.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00,
          0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00,
          -2.000000000000000e+00, 1.000000000000000e+00, 0.000000000000000e+00,
          -1.000000000000000e+00, 1.000000000000000e+00, 0.000000000000000e+00,
          -6.000000000000000e+00, 6.000000000000000e+00, 0.000000000000000e+00,
          2.000000000000000e+00, 0.000000000000000e+00, -1.000000000000000e+00,
          1.000000000000000e+00, 0.000000000000000e+00, -1.000000000000000e+00,
          -2.000000000000000e+00, 1.000000000000000e+00, 1.000000000000000e+00,
          2.000000000000000e+00, 0.000000000000000e+00, -2.000000000000000e+00,
          0.000000000000000e+00, 0.000000000000000e+00, -6.000000000000000e+00,
          1.000000000000000e+00, 0.000000000000000e+00, -2.000000000000000e+00,
          -2.000000000000000e+00, 0.000000000000000e+00, 2.000000000000000e+00,
          0.000000000000000e+00, 1.000000000000000e+00, -2.000000000000000e+00,
          -1.000000000000000e+00, -1.000000000000000e+00, 2.000000000000000e+00,
          0.000000000000000e+00, -2.000000000000000e+00, 2.000000000000000e+00,
          0.000000000000000e+00, 0.000000000000000e+00, -1.000000000000000e+00,
          -1.000000000000000e+00, 0.000000000000000e+00, 1.000000000000000e+00,
          0.000000000000000e+00, -1.000000000000000e+00, 1.000000000000000e+00,
          0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00,
          0.000000000000000e+00, 0.000000000000000e+00, 2.000000000000000e+00,
          -1.000000000000000e+00, 0.000000000000000e+00, 2.000000000000000e+00,
          -2.000000000000000e+00, 0.000000000000000e+00, 2.000000000000000e+00,
          0.000000000000000e+00, 1.000000000000000e+00, -2.000000000000000e+00,
          1.000000000000000e+00, 1.000000000000000e+00, -2.000000000000000e+00,
          0.000000000000000e+00, 6.000000000000000e+00, -6.000000000000000e+00,
          0.000000000000000e+00, 0.000000000000000e+00, 1.000000000000000e+00,
          -1.000000000000000e+00, 0.000000000000000e+00, 1.000000000000000e+00,
          0.000000000000000e+00, 1.000000000000000e+00, -1.000000000000000e+00,
          0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00,
          0.000000000000000e+00, 0.000000000000000e+00, 2.000000000000000e+00,
          1.000000000000000e+00, 0.000000000000000e+00, -2.000000000000000e+00,
          6.000000000000000e+00, 0.000000000000000e+00, -6.000000000000000e+00,
          0.000000000000000e+00, -1.000000000000000e+00, 2.000000000000000e+00,
          1.000000000000000e+00, 1.000000000000000e+00, -2.000000000000000e+00,
          0.000000000000000e+00, -2.000000000000000e+00, 2.000000000000000e+00,
          0.000000000000000e+00, 0.000000000000000e+00, 1.000000000000000e+00,
          1.000000000000000e+00, 0.000000000000000e+00, -1.000000000000000e+00,
          0.000000000000000e+00, -1.000000000000000e+00, 1.000000000000000e+00,
          0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00,
          0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00,
          2.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00,
          0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00,
          0.000000000000000e+00, 1.000000000000000e+00, 0.000000000000000e+00,
          1.000000000000000e+00, -1.000000000000000e+00, 0.000000000000000e+00,
          0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00,
          0.000000000000000e+00, 0.000000000000000e+00, 1.000000000000000e+00,
          1.000000000000000e+00, 0.000000000000000e+00, -1.000000000000000e+00,
          0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00,
          0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00,
          0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00,
          1.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00,
          0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00,
          0.000000000000000e+00, 2.000000000000000e+00, 0.000000000000000e+00,
          -1.000000000000000e+00, 1.000000000000000e+00, 0.000000000000000e+00,
          0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00,
          0.000000000000000e+00, 0.000000000000000e+00, 1.000000000000000e+00,
          0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00,
          0.000000000000000e+00, 1.000000000000000e+00, -1.000000000000000e+00,
          0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00,
          0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00,
          1.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00,
          0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00,
          0.000000000000000e+00, 1.000000000000000e+00, 0.000000000000000e+00,
          0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00,
          0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00,
          0.000000000000000e+00, 0.000000000000000e+00, 2.000000000000000e+00,
          -1.000000000000000e+00, 0.000000000000000e+00, 1.000000000000000e+00,
          0.000000000000000e+00, -1.000000000000000e+00, 1.000000000000000e+00,
          0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00

      };

      ordinal_type cur=0;
      const auto allTags = tetBasis.getAllDofTags();
      for (ordinal_type i=0;i<cardinality;i++) {
        auto scaling = scaling_factor[allTags(i,0)];
        for (ordinal_type j=0;j<np_lattice;j++) {
          for (ordinal_type k=0;k<dim; k++) {
            if (std::abs( h_basisAtLattice(i,j,k) - scaling*fiat_vals[cur] ) > tol ) {
              errorFlag++;
              *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";

              // Output the multi-index of the value where the error is:
              *outStream << " At multi-index { ";
              *outStream << i << " " << j << " " << k;
              *outStream << "}  computed value: " <<  h_basisAtLattice(i,j,k)
                                          << " but correct value: " << scaling*fiat_vals[cur] << "\n";
              *outStream << "Difference: " << std::abs(  h_basisAtLattice(i,j,k) - scaling*fiat_vals[cur] ) << "\n";
            }
            cur++;
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

  try {

    *outStream
    << "\n"
    << "===============================================================================\n"
    << "| TEST 4: Testing Operator DIV                                                |\n"
    << "===============================================================================\n";


    constexpr ordinal_type order = 2;
    if(order <= maxOrder) {
      TetBasisType tetBasis(order, POINTTYPE_EQUISPACED);

      shards::CellTopology tet_4(shards::getCellTopologyData<shards::Tetrahedron<4> >());
      const ordinal_type np_lattice = PointTools::getLatticeSize(tet_4, order,0);
      const ordinal_type cardinality = tetBasis.getCardinality();
      DynRankViewScalarValueType ConstructWithLabel(lattice_scalar, np_lattice , dim);
      PointTools::getLattice(lattice_scalar, tet_4, order, 0, POINTTYPE_EQUISPACED);
      DynRankViewPointValueType ConstructWithLabelPointView(lattice, np_lattice , dim);
      RealSpaceTools<DeviceType>::clone(lattice,lattice_scalar);

      DynRankViewOutValueType ConstructWithLabelOutView(basisDivAtLattice, cardinality , np_lattice);
      tetBasis.getValues(basisDivAtLattice, lattice, OPERATOR_DIV);

      auto h_basisDivAtLattice = Kokkos::create_mirror_view(basisDivAtLattice);
      Kokkos::deep_copy(h_basisDivAtLattice, basisDivAtLattice);

      const scalar_type fiat_divs[] = {
          2.600000000000000e+01,
          1.000000000000000e+01,
          -6.000000000000000e+00,
          1.000000000000000e+01,
          -6.000000000000000e+00,
          -6.000000000000000e+00,
          1.000000000000000e+01,
          -6.000000000000000e+00,
          -6.000000000000000e+00,
          -6.000000000000000e+00,
          -6.000000000000000e+00,
          1.000000000000000e+01,
          2.600000000000000e+01,
          -6.000000000000000e+00,
          1.000000000000000e+01,
          -6.000000000000000e+00,
          -6.000000000000000e+00,
          1.000000000000000e+01,
          -6.000000000000000e+00,
          -6.000000000000000e+00,
          -6.000000000000000e+00,
          -6.000000000000000e+00,
          -6.000000000000000e+00,
          -6.000000000000000e+00,
          -6.000000000000000e+00,
          -6.000000000000000e+00,
          1.000000000000000e+01,
          1.000000000000000e+01,
          1.000000000000000e+01,
          2.600000000000000e+01,
          -6.000000000000000e+00,
          1.000000000000000e+01,
          2.600000000000000e+01,
          -6.000000000000000e+00,
          1.000000000000000e+01,
          -6.000000000000000e+00,
          -6.000000000000000e+00,
          1.000000000000000e+01,
          -6.000000000000000e+00,
          -6.000000000000000e+00,
          -6.000000000000000e+00,
          -6.000000000000000e+00,
          -6.000000000000000e+00,
          1.000000000000000e+01,
          1.000000000000000e+01,
          2.600000000000000e+01,
          -6.000000000000000e+00,
          -6.000000000000000e+00,
          1.000000000000000e+01,
          -6.000000000000000e+00,
          -6.000000000000000e+00,
          -6.000000000000000e+00,
          -6.000000000000000e+00,
          -6.000000000000000e+00,
          -6.000000000000000e+00,
          -6.000000000000000e+00,
          1.000000000000000e+01,
          1.000000000000000e+01,
          1.000000000000000e+01,
          2.600000000000000e+01,
          2.600000000000000e+01,
          1.000000000000000e+01,
          -6.000000000000000e+00,
          1.000000000000000e+01,
          -6.000000000000000e+00,
          -6.000000000000000e+00,
          1.000000000000000e+01,
          -6.000000000000000e+00,
          -6.000000000000000e+00,
          -6.000000000000000e+00,
          -6.000000000000000e+00,
          -6.000000000000000e+00,
          -6.000000000000000e+00,
          -6.000000000000000e+00,
          -6.000000000000000e+00,
          -6.000000000000000e+00,
          1.000000000000000e+01,
          1.000000000000000e+01,
          1.000000000000000e+01,
          2.600000000000000e+01,
          -6.000000000000000e+00,
          -6.000000000000000e+00,
          -6.000000000000000e+00,
          1.000000000000000e+01,
          1.000000000000000e+01,
          2.600000000000000e+01,
          -6.000000000000000e+00,
          -6.000000000000000e+00,
          1.000000000000000e+01,
          -6.000000000000000e+00,
          2.600000000000000e+01,
          1.000000000000000e+01,
          -6.000000000000000e+00,
          1.000000000000000e+01,
          -6.000000000000000e+00,
          -6.000000000000000e+00,
          1.000000000000000e+01,
          -6.000000000000000e+00,
          -6.000000000000000e+00,
          -6.000000000000000e+00,
          -6.000000000000000e+00,
          -6.000000000000000e+00,
          -6.000000000000000e+00,
          1.000000000000000e+01,
          1.000000000000000e+01,
          2.600000000000000e+01,
          -6.000000000000000e+00,
          -6.000000000000000e+00,
          1.000000000000000e+01,
          -6.000000000000000e+00,
          -6.000000000000000e+00,
          1.000000000000000e+01,
          2.600000000000000e+01,
          -6.000000000000000e+00,
          1.000000000000000e+01,
          -6.000000000000000e+00,
          -6.000000000000000e+00,
          1.000000000000000e+01,
          -6.000000000000000e+00,
          -6.000000000000000e+00,
          1.600000000000000e+01,
          0.000000000000000e+00,
          -1.600000000000000e+01,
          8.000000000000000e+00,
          -8.000000000000000e+00,
          0.000000000000000e+00,
          8.000000000000000e+00,
          -8.000000000000000e+00,
          0.000000000000000e+00,
          0.000000000000000e+00,
          1.600000000000000e+01,
          8.000000000000000e+00,
          0.000000000000000e+00,
          0.000000000000000e+00,
          -8.000000000000000e+00,
          -1.600000000000000e+01,
          8.000000000000000e+00,
          0.000000000000000e+00,
          -8.000000000000000e+00,
          0.000000000000000e+00,
          1.600000000000000e+01,
          8.000000000000000e+00,
          0.000000000000000e+00,
          8.000000000000000e+00,
          0.000000000000000e+00,
          0.000000000000000e+00,
          0.000000000000000e+00,
          -8.000000000000000e+00,
          -8.000000000000000e+00,
          -1.600000000000000e+01

      };

      const auto allTags = tetBasis.getAllDofTags();
      ordinal_type cur=0;
      for (ordinal_type i=0;i<cardinality;i++) {
        auto scaling = scaling_factor[allTags(i,0)];
        for (ordinal_type j=0;j<np_lattice;j++) {
          if (std::abs( h_basisDivAtLattice(i,j) - scaling*fiat_divs[cur] ) > 10*tol ) {
            errorFlag++;
            *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";

            // Output the multi-index of the value where the error is:
            *outStream << " At multi-index { ";
            *outStream << i << " " << j;
            *outStream << "}  computed value: " <<  h_basisDivAtLattice(i,j)
                                          << " but correct value: " << scaling*fiat_divs[cur] << "\n";
            *outStream << "Difference: " << std::abs(  h_basisDivAtLattice(i,j) - scaling*fiat_divs[cur] ) << "\n";
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
    const ordinal_type order = std::min(4, maxOrder);
    TetBasisType tetBasis(order, POINTTYPE_EQUISPACED);
    
    const EFunctionSpace fs = tetBasis.getFunctionSpace();
    
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
