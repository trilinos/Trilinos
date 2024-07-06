// @HEADER
// *****************************************************************************
//                           Intrepid2 Package
//
// Copyright 2007 NTESS and the Intrepid2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/** \file test_04.cpp
    \brief  Unit test for the FunctionSpaceTools class, testing volume
    \author Created by D. Ridzal, P. Bochev, K. Peterson and Kyungjoo Kim.
*/
#include "Intrepid2_config.h"

#ifdef HAVE_INTREPID2_DEBUG
#define INTREPID2_TEST_FOR_DEBUG_ABORT_OVERRIDE_TO_CONTINUE
#endif

#include "Intrepid2_Types.hpp"
#include "Intrepid2_Utils.hpp"
#include "Intrepid2_Utils_ExtData.hpp"

#include "Intrepid2_CellTools.hpp"
#include "Intrepid2_FunctionSpaceTools.hpp"

#include "Intrepid2_DefaultCubatureFactory.hpp"

#include "Teuchos_oblackholestream.hpp"
#include "Teuchos_RCP.hpp"

namespace Intrepid2 {

  namespace Test {
#define INTREPID2_TEST_ERROR_EXPECTED( S )                              \
    try {                                                               \
      ++nthrow;                                                         \
      S ;                                                               \
    } catch (std::logic_error &err) {                                    \
      ++ncatch;                                                         \
      *outStream << "Expected Error ----------------------------------------------------------------\n"; \
      *outStream << err.what() << '\n';                                 \
      *outStream << "-------------------------------------------------------------------------------" << "\n\n"; \
    };
#define ConstructWithLabel(obj, ...) obj(#obj, __VA_ARGS__)
    
    template<typename ValueType, typename DeviceType>
    int FunctionSpaceTools_Test04(const bool verbose) {
      using ExecSpaceType = typename DeviceType::execution_space;

      Teuchos::RCP<std::ostream> outStream;
      Teuchos::oblackholestream bhs; // outputs nothing

      if (verbose)
        outStream = Teuchos::rcp(&std::cout, false);
      else
        outStream = Teuchos::rcp(&bhs, false);

      Teuchos::oblackholestream oldFormatState;
      oldFormatState.copyfmt(std::cout);

      typedef typename
        Kokkos::DefaultHostExecutionSpace HostSpaceType ;

      *outStream << "DeviceSpace::  ";   ExecSpaceType().print_configuration(*outStream, false);
      *outStream << "HostSpace::    ";   HostSpaceType().print_configuration(*outStream, false);

      *outStream                                                       
        << "===============================================================================\n"
        << "|                                                                             |\n"
        << "|                      Unit Test (FunctionSpaceTools)                         |\n"
        << "|                                                                             |\n"
        << "|     1) volume integration on tetrahedra, testing dataIntegral               |\n"
        << "|                                                                             |\n"
        << "|  Questions? Contact  Pavel Bochev (pbboche@sandia.gov) or                   |\n"
        << "|                      Denis Ridzal (dridzal@sandia.gov) or                   |\n"
        << "|                      Kyungjoo Kim (kyukim@sandia.gov).                      |\n"
        << "|                                                                             |\n"
        << "|  Intrepid's website: http://trilinos.sandia.gov/packages/intrepid           |\n"
        << "|  Trilinos website:   http://trilinos.sandia.gov                             |\n"
        << "|                                                                             |\n"
        << "===============================================================================\n";

      typedef CellTools<DeviceType> ct;
      typedef FunctionSpaceTools<DeviceType> fst;
      typedef Kokkos::DynRankView<ValueType,DeviceType> DynRankView;

      const auto tol = tolerence();

      int errorFlag = 0;

      *outStream
        << "\n"
        << "===============================================================================\n"
        << "| TEST 1: correctness of cell volumes                                         |\n"
        << "===============================================================================\n";

      outStream->precision(20);

      try {
        DefaultCubatureFactory cub_factory;

        shards::CellTopology cell_topo = shards::getCellTopologyData< shards::Tetrahedron<4> >();

        const auto cub_degree = 0;
        auto cub = cub_factory.create<DeviceType,ValueType,ValueType>(cell_topo, cub_degree);

        const auto space_dim = cub->getDimension();
        const auto num_cub_points = cub->getNumPoints();

        /* Cell geometries. */
        const auto num_cells = 4;
        const auto num_nodes = 4;

        const ValueType tetnodes[] = {
          // tet 0
          0.0, 0.0, 0.0,
          1.0, 0.0, 0.0,
          0.0, 1.0, 0.0,
          0.0, 0.0, 1.0,
          // tet 1
          4.0, 5.0, 1.0,
          -6.0, 2.0, 0.0,
          4.0, -3.0, -1.0,
          0.0, 2.0, 5.0,
          // tet 2
          -6.0, -3.0, 1.0,
          9.0, 2.0, 1.0,
          8.9, 2.1, 0.9,
          8.9, 2.1, 1.1,
          // tet 3
          -6.0, -3.0, 1.0,
          12.0, 3.0, 1.0,
          2.9, 0.1, 0.9,
          2.9, 0.1, 1.1
        };
        
        /* Analytic volumes. */
        const ValueType tetvols[] = {1.0/6.0, 194.0/3.0, 1.0/15.0, 2.0/25.0};

      /* Computational arrays. */
        DynRankView ConstructWithLabel( cub_points,  num_cub_points, space_dim);
        DynRankView ConstructWithLabel( cub_weights, num_cub_points);

        DynRankView ConstructWithLabel( cell_nodes,       num_cells, num_nodes, space_dim);

        DynRankView ConstructWithLabel( jacobian,         num_cells, num_cub_points, space_dim, space_dim);
        DynRankView ConstructWithLabel( jacobian_det,     num_cells, num_cub_points);
        DynRankView ConstructWithLabel( weighted_measure, num_cells, num_cub_points);

        DynRankView ConstructWithLabel( data_one, num_cells, num_cub_points);
        DynRankView ConstructWithLabel( volumes, num_cells);

        /******************* START COMPUTATION ***********************/

        // get cubature points and weights
        cub->getCubature(cub_points, cub_weights);

        const Kokkos::DynRankView<ValueType,Kokkos::LayoutRight,Kokkos::HostSpace> cell_nodes_host (const_cast<ValueType*>(&tetnodes[0]),  num_cells, num_nodes, space_dim);

        // 1. mirror allocation
        // 2. deep copy preserving layout
        // 3. remap to native layout of the device
        auto cell_nodes_device  = create_mirror_view(typename DeviceType::memory_space(), cell_nodes_host);
        Kokkos::deep_copy( cell_nodes_device, cell_nodes_host );
        Kokkos::deep_copy( cell_nodes , cell_nodes_device );

        // compute geometric cell information
        ct::setJacobian(jacobian, cub_points, cell_nodes, cell_topo);
        ct::setJacobianDet(jacobian_det, jacobian);

        // compute weighted measure
        fst::computeCellMeasure(weighted_measure, jacobian_det, cub_weights);

        Kokkos::deep_copy(data_one, 1.0);

        // compute volumes
        fst::integrate(volumes, data_one, weighted_measure);
        
        /******************* STOP COMPUTATION ***********************/
        
        // memcpy do implicit sync
        const auto volumes_host = Kokkos::create_mirror_view(typename HostSpaceType::memory_space(), volumes);
        Kokkos::deep_copy(volumes_host, volumes);

        /******************* START COMPARISON ***********************/
        for (auto cid=0;cid<num_cells-1;++cid) {
          *outStream << "Volume of cell " << cid 
                     << " = " << std::setw(24) << volumes_host(cid) 
                     << "    vs.    Analytic value =  " << std::setw(24) << tetvols[cid] << "\n";
          errorFlag += (std::fabs(volumes_host(cid)-tetvols[cid]) > tol);
          
        }
        /******************* STOP COMPARISON ***********************/
        
        *outStream << "\n";
      } catch (std::logic_error &err) {
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
