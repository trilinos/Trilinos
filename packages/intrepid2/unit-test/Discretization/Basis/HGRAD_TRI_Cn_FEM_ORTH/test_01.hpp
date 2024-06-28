// @HEADER
// *****************************************************************************
//                           Intrepid2 Package
//
// Copyright 2007 NTESS and the Intrepid2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/** \file test_01.hpp
 \brief  Unit tests for the Intrepid2::HGRAD_TRI_CN_FEM_ORTH class.
 \author Created by R. Kirby, M. Perego and K. Kim.
 */

#include "Intrepid2_config.h"

#ifdef HAVE_INTREPID2_DEBUG
#define INTREPID2_TEST_FOR_DEBUG_ABORT_OVERRIDE_TO_CONTINUE
#endif

#include "Intrepid2_HGRAD_TRI_Cn_FEM_ORTH.hpp"
#include "Intrepid2_PointTools.hpp"


#include "Teuchos_oblackholestream.hpp"
#include "Teuchos_RCP.hpp"
#include "Intrepid2_CubatureDirectTriDefault.hpp"

#include "packages/intrepid2/unit-test/Discretization/Basis/Macros.hpp"
#include "packages/intrepid2/unit-test/Discretization/Basis/Setup.hpp"

namespace Intrepid2 {

namespace Test {

template<typename ValueType, typename DeviceType>
int HGRAD_TRI_Cn_FEM_ORTH_Test01(const bool verbose) {

  //! Create an execution space instance.
  const auto space = Kokkos::Experimental::partition_space(typename DeviceType::execution_space {}, 1)[0];

  Teuchos::RCP<std::ostream> outStream = setup_output_stream<DeviceType>(
    verbose, "OrthogonalBases", {
        "1) Tests orthogonality of triangular orthogonal basis"
  });

  Teuchos::oblackholestream oldFormatState;
  oldFormatState.copyfmt(std::cout);

  typedef Kokkos::DynRankView<ValueType, DeviceType> DynRankView;

  const ValueType tol = tolerence();
  int errorFlag = 0;

  // for virtual function, value and point types are declared in the class
  typedef ValueType outputValueType;
  typedef ValueType pointValueType;
  typedef ValueType weightValueType;
  typedef Basis_HGRAD_TRI_Cn_FEM_ORTH<DeviceType, outputValueType,
      pointValueType> triBasisType;
  typedef CubatureDirectTriDefault<DeviceType, pointValueType,
      weightValueType> cubatureTriType;
  *outStream << "\n"
      << "===============================================================================\n"
      << "| TEST 1: constructors and exceptions                                         |\n"
      << "===============================================================================\n";


  const ordinal_type dim = 2;
  constexpr ordinal_type maxOrder = Parameters::MaxOrder ;

  try {
    const ordinal_type order = std::min(3, maxOrder);
    triBasisType triBasis(order);

    const ordinal_type polydim = triBasis.getCardinality();

    // First, get a reference quadrature rule
    cubatureTriType cub(2 * order);
    const ordinal_type npts = cub.getNumPoints();
    DynRankView ConstructWithLabel(cubPoints, npts, dim);
    DynRankView ConstructWithLabel(cubWeights, npts);

    cub.getCubature(cubPoints, cubWeights);
    auto h_cubWeights = Kokkos::create_mirror_view(cubWeights);
    Kokkos::deep_copy(h_cubWeights, cubWeights);

    // Tabulate the basis functions at the cubature points
    DynRankView ConstructWithLabel(basisAtCubPts, polydim, npts);

    triBasis.getValues(space, basisAtCubPts, cubPoints, OPERATOR_VALUE);
    auto h_basisAtCubPts = Kokkos::create_mirror_view(basisAtCubPts);
    Kokkos::deep_copy(h_basisAtCubPts, basisAtCubPts);

    // Now let's compute the mass matrix
    for (ordinal_type i = 0; i < polydim; i++) {
      for (ordinal_type j = i; j < polydim; j++) {
        ValueType cur = 0.0;
        for (ordinal_type k = 0; k < cub.getNumPoints(); k++) {
          cur += h_cubWeights(k) * h_basisAtCubPts(i, k)* h_basisAtCubPts(j, k);
        }
        if (i != j && fabs(cur) > tol) {
          errorFlag++;
        }
        if (i == j && fabs(cur) < tol) {
          errorFlag++;
        }
      }
    }
  }
  catch (std::exception &err) {
    *outStream << err.what() << "\n\n";
    errorFlag = -1000;
  }

  *outStream << "\n"
      << "===============================================================================\n"
      << "| TEST 2: Operator D1                                                         |\n"
      << "===============================================================================\n";

  try {
    constexpr ordinal_type order = 3;
    if(order <= maxOrder) {
      triBasisType triBasis(order);
      shards::CellTopology tri_3(
          shards::getCellTopologyData<shards::Triangle<3> >());
      const ordinal_type np_lattice = PointTools::getLatticeSize(tri_3, order,
          0);
      const ordinal_type polydim = triBasis.getCardinality();
      DynRankView ConstructWithLabel(lattice, np_lattice , dim);
      PointTools::getLattice(lattice, tri_3, order, 0, POINTTYPE_EQUISPACED);

      DynRankView ConstructWithLabel(dBasisAtLattice, polydim , np_lattice , dim);
      triBasis.getValues(space, dBasisAtLattice, lattice, OPERATOR_D1);

      auto h_dBasisAtLattice = Kokkos::create_mirror_view(dBasisAtLattice);
      Kokkos::deep_copy(h_dBasisAtLattice, dBasisAtLattice);

      const double fiat_vals[] = {
          0.000000000000000e+00, 0.000000000000000e+00,
          0.000000000000000e+00, 0.000000000000000e+00,
          0.000000000000000e+00, 0.000000000000000e+00,
          0.000000000000000e+00, 0.000000000000000e+00,
          0.000000000000000e+00, 0.000000000000000e+00,
          0.000000000000000e+00, 0.000000000000000e+00,
          0.000000000000000e+00, 0.000000000000000e+00,
          0.000000000000000e+00, 0.000000000000000e+00,
          0.000000000000000e+00, 0.000000000000000e+00,
          0.000000000000000e+00, 0.000000000000000e+00,
          3.464101615137754e+00, 1.732050807568877e+00,
          3.464101615137754e+00, 1.732050807568877e+00,
          3.464101615137754e+00, 1.732050807568877e+00,
          3.464101615137754e+00, 1.732050807568877e+00,
          3.464101615137754e+00, 1.732050807568877e+00,
          3.464101615137754e+00, 1.732050807568877e+00,
          3.464101615137754e+00, 1.732050807568877e+00,
          3.464101615137754e+00, 1.732050807568877e+00,
          3.464101615137754e+00, 1.732050807568877e+00,
          3.464101615137754e+00, 1.732050807568877e+00,
          0.000000000000000e+00, 3.000000000000000e+00,
          0.000000000000000e+00, 3.000000000000000e+00,
          0.000000000000000e+00, 3.000000000000000e+00,
          0.000000000000000e+00, 3.000000000000000e+00,
          0.000000000000000e+00, 3.000000000000000e+00,
          0.000000000000000e+00, 3.000000000000000e+00,
          0.000000000000000e+00, 3.000000000000000e+00,
          0.000000000000000e+00, 3.000000000000000e+00,
          0.000000000000000e+00, 3.000000000000000e+00,
          0.000000000000000e+00, 3.000000000000000e+00,
          -1.643167672515498e+01, -5.477225575051661e+00,
          -5.477225575051661e+00, 0.000000000000000e+00,
          5.477225575051660e+00, 5.477225575051660e+00,
          1.643167672515498e+01, 1.095445115010332e+01,
          -1.095445115010332e+01, -3.651483716701107e+00,
          -9.121412916732176e-16, 1.825741858350553e+00,
          1.095445115010332e+01, 7.302967433402213e+00,
          -5.477225575051661e+00, -1.825741858350554e+00,
          5.477225575051660e+00, 3.651483716701107e+00,
          0.000000000000000e+00, 0.000000000000000e+00,
          -4.242640687119285e+00, -1.272792206135786e+01,
          -4.242640687119285e+00, -5.656854249492381e+00,
          -4.242640687119285e+00, 1.414213562373094e+00,
          -4.242640687119285e+00, 8.485281374238570e+00,
          2.828427124746189e+00, -5.656854249492381e+00,
          2.828427124746189e+00, 1.414213562373094e+00,
          2.828427124746189e+00, 8.485281374238568e+00,
          9.899494936611664e+00, 1.414213562373094e+00,
          9.899494936611664e+00, 8.485281374238568e+00,
          1.697056274847714e+01, 8.485281374238570e+00,
          0.000000000000000e+00, -9.797958971132712e+00,
          0.000000000000000e+00, -9.797958971132712e+00,
          0.000000000000000e+00, -9.797958971132710e+00,
          0.000000000000000e+00, -9.797958971132712e+00,
          0.000000000000000e+00, -1.632993161855452e+00,
          0.000000000000000e+00, -1.632993161855452e+00,
          0.000000000000000e+00, -1.632993161855452e+00,
          0.000000000000000e+00, 6.531972647421806e+00,
          0.000000000000000e+00, 6.531972647421806e+00,
          0.000000000000000e+00, 1.469693845669907e+01,
          4.489988864128730e+01, 1.122497216032182e+01,
          -4.988876515698587e+00, -6.236095644623235e+00,
          -4.988876515698591e+00, 1.247219128924645e+00,
          4.489988864128730e+01, 3.367491648096547e+01,
          1.995550606279436e+01, 4.988876515698590e+00,
          -4.988876515698589e+00, -2.494438257849295e+00,
          1.995550606279435e+01, 1.496662954709576e+01,
          4.988876515698590e+00, 1.247219128924648e+00,
          4.988876515698586e+00, 3.741657386773940e+00,
          0.000000000000000e+00, 0.000000000000000e+00,
          1.897366596101028e+01, 2.846049894151541e+01,
          6.324555320336759e+00, -7.378647873726218e+00,
          -6.324555320336757e+00, -1.370320319406298e+01,
          -1.897366596101028e+01, 9.486832980505138e+00,
          -1.686548085423136e+01, 4.216370213557840e+00,
          -1.404333387430680e-15, -2.108185106778921e+00,
          1.686548085423135e+01, 2.108185106778919e+01,
          -2.319003617456811e+01, -5.270462766947300e+00,
          2.319003617456811e+01, 1.791957340762081e+01,
          0.000000000000000e+00, 0.000000000000000e+00,
          4.898979485566356e+00, 3.184336665618131e+01,
          4.898979485566356e+00, 1.224744871391589e+01,
          4.898979485566356e+00, -7.348469228349532e+00,
          4.898979485566356e+00, -2.694438717061496e+01,
          -3.265986323710904e+00, -4.898979485566356e+00,
          -3.265986323710904e+00, -1.632993161855452e+00,
          -3.265986323710904e+00, 1.632993161855451e+00,
          1.143095213298816e+01, -7.348469228349537e+00,
          1.143095213298816e+01, 1.877942136133769e+01,
          4.898979485566356e+01, 2.449489742783178e+01,
          0.000000000000000e+00, 2.121320343559643e+01,
          0.000000000000000e+00, 2.121320343559643e+01,
          0.000000000000000e+00, 2.121320343559642e+01,
          0.000000000000000e+00, 2.121320343559643e+01,
          0.000000000000000e+00, -4.714045207910317e+00,
          0.000000000000000e+00, -4.714045207910317e+00,
          0.000000000000000e+00, -4.714045207910317e+00,
          0.000000000000000e+00, 2.357022603955157e+00,
          0.000000000000000e+00, 2.357022603955157e+00,
          0.000000000000000e+00, 4.242640687119285e+01
      };

      int fiat_index_cur = 0;
      for (ordinal_type i=0;i<polydim;i++) {
        for (ordinal_type j=0;j<np_lattice;j++) {
          for (ordinal_type k=0;k<dim;k++) {
            if (std::abs( h_dBasisAtLattice(i,j,k) - fiat_vals[fiat_index_cur] ) > 10.0*tol ) {
              errorFlag++;
              *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";

              // Output the multi-index of the value where the error is:
              *outStream << " At multi-index { ";
              *outStream << i << " " << j << " " << k;
              *outStream << "}  computed value: " << h_dBasisAtLattice(i,j,k)
                                             << " but correct value: " << fiat_vals[fiat_index_cur] << "\n";
              *outStream << "Difference: " << std::abs( h_dBasisAtLattice(i,j,k) - fiat_vals[fiat_index_cur] ) << "\n";
            }
            fiat_index_cur++;
          }
        }
      }
    }
  }  catch (std::exception &err) {
    *outStream << err.what() << "\n\n";
    errorFlag = -1000;
  }

  // do second order derivatives
#if 0 //#ifdef HAVE_INTREPID2_SACADO
  try {
    *outStream << "\n"
        << "===============================================================================\n"
        << "| TEST 3: Operator D2                                                         |\n"
        << "===============================================================================\n";

    constexpr ordinal_type order = 3;
    if(order <= maxOrder) {
      triBasisType triBasis(order);
      shards::CellTopology tri_3(
          shards::getCellTopologyData<shards::Triangle<3> >());
      const ordinal_type np_lattice = PointTools::getLatticeSize(tri_3, order,
          0);
      const ordinal_type polydim = triBasis.getCardinality();
      DynRankView ConstructWithLabel(lattice, np_lattice , dim);
      PointTools::getLattice(lattice, tri_3, order, 0, POINTTYPE_EQUISPACED);

      DynRankView ConstructWithLabel(dBasisAtLattice, polydim , np_lattice , 3);
      triBasis.getValues(space, dBasisAtLattice, lattice, OPERATOR_D2);

      auto h_dBasisAtLattice = Kokkos::create_mirror_view(dBasisAtLattice);
      Kokkos::deep_copy(h_dBasisAtLattice, dBasisAtLattice);

      const double fiat_vals[] = {
          0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00,
          0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00,
          0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00,
          0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00,
          0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00,
          0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00,
          0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00,
          0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00,
          0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00,
          0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00,
          0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00,
          0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00,
          0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00,
          0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00,
          0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00,
          0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00,
          0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00,
          0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00,
          0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00,
          0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00,
          0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00,
          0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00,
          0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00,
          0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00,
          0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00,
          0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00,
          0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00,
          0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00,
          0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00,
          0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00,
          3.286335345030997e+01, 1.643167672515498e+01, 5.477225575051661e+00,
          3.286335345030997e+01, 1.643167672515498e+01, 5.477225575051661e+00,
          3.286335345030997e+01, 1.643167672515498e+01, 5.477225575051661e+00,
          3.286335345030997e+01, 1.643167672515498e+01, 5.477225575051661e+00,
          3.286335345030997e+01, 1.643167672515498e+01, 5.477225575051661e+00,
          3.286335345030997e+01, 1.643167672515498e+01, 5.477225575051661e+00,
          3.286335345030997e+01, 1.643167672515498e+01, 5.477225575051661e+00,
          3.286335345030997e+01, 1.643167672515498e+01, 5.477225575051661e+00,
          3.286335345030997e+01, 1.643167672515498e+01, 5.477225575051661e+00,
          3.286335345030997e+01, 1.643167672515498e+01, 5.477225575051661e+00,
          0.000000000000000e+00, 2.121320343559642e+01, 2.121320343559642e+01,
          0.000000000000000e+00, 2.121320343559642e+01, 2.121320343559642e+01,
          0.000000000000000e+00, 2.121320343559642e+01, 2.121320343559642e+01,
          0.000000000000000e+00, 2.121320343559642e+01, 2.121320343559642e+01,
          0.000000000000000e+00, 2.121320343559642e+01, 2.121320343559642e+01,
          0.000000000000000e+00, 2.121320343559642e+01, 2.121320343559642e+01,
          0.000000000000000e+00, 2.121320343559642e+01, 2.121320343559642e+01,
          0.000000000000000e+00, 2.121320343559642e+01, 2.121320343559642e+01,
          0.000000000000000e+00, 2.121320343559642e+01, 2.121320343559642e+01,
          0.000000000000000e+00, 2.121320343559642e+01, 2.121320343559642e+01,
          0.000000000000000e+00, 0.000000000000000e+00, 2.449489742783178e+01,
          0.000000000000000e+00, 0.000000000000000e+00, 2.449489742783178e+01,
          0.000000000000000e+00, 0.000000000000000e+00, 2.449489742783177e+01,
          0.000000000000000e+00, 0.000000000000000e+00, 2.449489742783178e+01,
          0.000000000000000e+00, 0.000000000000000e+00, 2.449489742783178e+01,
          0.000000000000000e+00, 0.000000000000000e+00, 2.449489742783178e+01,
          0.000000000000000e+00, 0.000000000000000e+00, 2.449489742783178e+01,
          0.000000000000000e+00, 0.000000000000000e+00, 2.449489742783177e+01,
          0.000000000000000e+00, 0.000000000000000e+00, 2.449489742783177e+01,
          0.000000000000000e+00, 0.000000000000000e+00, 2.449489742783178e+01,
          -2.244994432064365e+02, -8.979977728257460e+01, -2.244994432064365e+01,
          -7.483314773547885e+01, -1.496662954709577e+01, 7.483314773547880e+00,
          7.483314773547882e+01, 5.986651818838305e+01, 3.741657386773942e+01,
          2.244994432064365e+02, 1.346996659238619e+02, 6.734983296193094e+01,
          -1.496662954709577e+02, -5.986651818838306e+01, -1.496662954709577e+01,
          -1.246222254316567e-14, 1.496662954709576e+01, 1.496662954709576e+01,
          1.496662954709576e+02, 8.979977728257458e+01, 4.489988864128730e+01,
          -7.483314773547885e+01, -2.993325909419154e+01, -7.483314773547884e+00,
          7.483314773547882e+01, 4.489988864128729e+01, 2.244994432064365e+01,
          0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00,
          -3.794733192202055e+01, -1.517893276880822e+02, -9.486832980505139e+01,
          -3.794733192202055e+01, -6.324555320336759e+01, -6.324555320336759e+00,
          -3.794733192202055e+01, 2.529822128134703e+01, 8.221921916437785e+01,
          -3.794733192202055e+01, 1.138419957660617e+02, 1.707629936490925e+02,
          5.059644256269407e+01, -6.324555320336759e+01, -5.059644256269407e+01,
          5.059644256269407e+01, 2.529822128134703e+01, 3.794733192202055e+01,
          5.059644256269407e+01, 1.138419957660616e+02, 1.264911064067352e+02,
          1.391402170474087e+02, 2.529822128134704e+01, -6.324555320336762e+00,
          1.391402170474087e+02, 1.138419957660617e+02, 8.221921916437786e+01,
          2.276839915321233e+02, 1.138419957660617e+02, 3.794733192202055e+01,
          0.000000000000000e+00, -5.878775382679627e+01, -1.616663230236897e+02,
          0.000000000000000e+00, -5.878775382679627e+01, -9.308061022576078e+01,
          0.000000000000000e+00, -5.878775382679627e+01, -2.449489742783179e+01,
          0.000000000000000e+00, -5.878775382679627e+01, 4.409081537009720e+01,
          0.000000000000000e+00, 9.797958971132708e+00, -5.878775382679628e+01,
          0.000000000000000e+00, 9.797958971132708e+00, 9.797958971132703e+00,
          0.000000000000000e+00, 9.797958971132708e+00, 7.838367176906168e+01,
          0.000000000000000e+00, 7.838367176906168e+01, 4.409081537009718e+01,
          0.000000000000000e+00, 7.838367176906168e+01, 1.126765281680262e+02,
          0.000000000000000e+00, 1.469693845669907e+02, 1.469693845669907e+02,
          0.000000000000000e+00, 0.000000000000000e+00, -1.272792206135786e+02,
          0.000000000000000e+00, 0.000000000000000e+00, -1.272792206135786e+02,
          0.000000000000000e+00, 0.000000000000000e+00, -1.272792206135785e+02,
          0.000000000000000e+00, 0.000000000000000e+00, -1.272792206135786e+02,
          0.000000000000000e+00, 0.000000000000000e+00, -2.828427124746190e+01,
          0.000000000000000e+00, 0.000000000000000e+00, -2.828427124746190e+01,
          0.000000000000000e+00, 0.000000000000000e+00, -2.828427124746190e+01,
          0.000000000000000e+00, 0.000000000000000e+00, 7.071067811865474e+01,
          0.000000000000000e+00, 0.000000000000000e+00, 7.071067811865474e+01,
          0.000000000000000e+00, 0.000000000000000e+00, 1.697056274847714e+02

      };
      int fiat_index_cur = 0;
      for (int i=0;i<polydim;i++) {
        for (int j=0;j<np_lattice;j++) {
          for (int k=0;k<3;k++) {
            if (std::abs( h_dBasisAtLattice(i,j,k) - fiat_vals[fiat_index_cur] ) > 10.0*tol ) {
              errorFlag++;
              *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";

              // Output the multi-index of the value where the error is:
              *outStream << " At multi-index { ";
              *outStream << i << " " << j << " " << k;
              *outStream << "}  computed value: " << h_dBasisAtLattice(i,j,k)
                                          << " but correct value: " << fiat_vals[fiat_index_cur] << "\n";
              *outStream << "Difference: " << std::abs( h_dBasisAtLattice(i,j,k) - fiat_vals[fiat_index_cur] ) << "\n";
            }
            fiat_index_cur++;
          }
        }
      }
    }
  }
  catch (std::exception &err) {
    *outStream << err.what() << "\n\n";
    errorFlag = -1000;
  }
#endif
  
  *outStream
  << "\n"
  << "===============================================================================\n"
  << "| TEST 4: Function Space is Correct                                           |\n"
  << "===============================================================================\n";
  
  try {
    for (auto order=0;order<std::min(3, maxOrder);++order) {
      triBasisType triBasis(order);
      
      const EFunctionSpace fs = triBasis.getFunctionSpace();
      
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




