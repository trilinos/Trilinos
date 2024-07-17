// @HEADER
// *****************************************************************************
//                           Intrepid2 Package
//
// Copyright 2007 NTESS and the Intrepid2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/** \file test_02.cpp
    \brief  Unit test for the FunctionSpaceTools class, testing H-curl.
    \author Created by D. Ridzal, P. Bochev, K. Peterson and Kyungjoo Kim
*/

#include "Intrepid2_config.h"

#ifdef HAVE_INTREPID2_DEBUG
#define INTREPID2_TEST_FOR_DEBUG_ABORT_OVERRIDE_TO_CONTINUE
#endif

#include "Intrepid2_Types.hpp"
#include "Intrepid2_Utils.hpp"
#include "Intrepid2_Utils_ExtData.hpp"

#include "Intrepid2_CellTools.hpp"
#include "Intrepid2_HGRAD_TET_C1_FEM.hpp"
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
    int FunctionSpaceTools_Test01(const bool verbose) {
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
        << "|     1) basic operator transformations and integration in HCURL              |\n"
        << "|                                                                             |\n"
        << "|  Questions? Contact  Pavel Bochev (pbboche@sandia.gov) or                   |\n"
        << "|                      Denis Ridzal (dridzal@sandia.gov).                     |\n"
        << "|                      Kyungjoo Kim (kyukim@sandia.gov).                      |\n"
        << "|                                                                             |\n"
        << "|  Intrepid's website: http://trilinos.sandia.gov/packages/intrepid           |\n"
        << "|  Trilinos website:   http://trilinos.sandia.gov                             |\n"
        << "|                                                                             |\n"
        << "===============================================================================\n";

      typedef CellTools<DeviceType> ct;
      typedef FunctionSpaceTools<DeviceType> fst;
      typedef Kokkos::DynRankView<ValueType,DeviceType> DynRankView;
      
      int errorFlag = 0;
      
      *outStream
        << "\n"
        << "===============================================================================\n"
        << "| TEST 1: correctness of math operations                                      |\n"
        << "===============================================================================\n";
      
      outStream->precision(20);
      
      try {

        DefaultCubatureFactory cub_factory;

        shards::CellTopology cell_topo = shards::getCellTopologyData< shards::Tetrahedron<4> >();

        const auto cub_degree = 2;
        auto cub = cub_factory.create<DeviceType,ValueType,ValueType>(cell_topo, cub_degree);

        const auto space_dim = cub->getDimension();
        const auto num_cub_points = cub->getNumPoints();

        Basis_HGRAD_TET_C1_FEM<DeviceType> tetBasis;
        const auto num_fields = tetBasis.getCardinality();

        /* Cell geometries and orientations. */
        const auto num_cells = 4;
        const auto num_nodes = 4;

        ValueType tetnodes[] = {
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

        /* Computational arrays. */
        DynRankView ConstructWithLabel( cub_points,  num_cub_points, space_dim);
        DynRankView ConstructWithLabel( cub_weights, num_cub_points);

        DynRankView ConstructWithLabel( cell_nodes,       num_cells, num_nodes, space_dim);
        DynRankView ConstructWithLabel( field_signs,      num_cells, num_fields);

        DynRankView ConstructWithLabel( jacobian,         num_cells, num_cub_points, space_dim, space_dim);
        DynRankView ConstructWithLabel( jacobian_inv,     num_cells, num_cub_points, space_dim, space_dim);
        DynRankView ConstructWithLabel( jacobian_det,     num_cells, num_cub_points);
        DynRankView ConstructWithLabel( weighted_measure, num_cells, num_cub_points);

        DynRankView ConstructWithLabel( grad_of_basis_at_cub_points,                                 num_fields, num_cub_points, space_dim);
        DynRankView ConstructWithLabel( transformed_grad_of_basis_at_cub_points,          num_cells, num_fields, num_cub_points, space_dim);
        DynRankView ConstructWithLabel( weighted_transformed_grad_of_basis_at_cub_points, num_cells, num_fields, num_cub_points, space_dim);
        DynRankView ConstructWithLabel( stiffness_matrices,                               num_cells, num_fields, num_fields);
        
        DynRankView ConstructWithLabel( value_of_basis_at_cub_points,                                 num_fields, num_cub_points);
        DynRankView ConstructWithLabel( transformed_value_of_basis_at_cub_points,          num_cells, num_fields, num_cub_points);
        DynRankView ConstructWithLabel( weighted_transformed_value_of_basis_at_cub_points, num_cells, num_fields, num_cub_points);
        DynRankView ConstructWithLabel( mass_matrices,                                     num_cells, num_fields, num_fields);

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
        ct::setJacobianInv(jacobian_inv, jacobian);
        ct::setJacobianDet(jacobian_det, jacobian);

        // compute weighted measure
        fst::computeCellMeasure(weighted_measure, jacobian_det, cub_weights);

        // ** Computing stiffness matrices:
        tetBasis.getValues(grad_of_basis_at_cub_points, cub_points, OPERATOR_GRAD);

        // transform grads of basis functions 
        fst::HGRADtransformGRAD(transformed_grad_of_basis_at_cub_points,
                                jacobian_inv,
                                grad_of_basis_at_cub_points);

        // multiply with weighted measure
        fst::multiplyMeasure(weighted_transformed_grad_of_basis_at_cub_points,
                             weighted_measure,
                             transformed_grad_of_basis_at_cub_points);
        
        // compute stiffness matrices
        fst::integrate(stiffness_matrices,
                       transformed_grad_of_basis_at_cub_points,
                       weighted_transformed_grad_of_basis_at_cub_points);

        // ** Computing mass matrices:
        tetBasis.getValues(value_of_basis_at_cub_points, cub_points, OPERATOR_VALUE);

        // transform values of basis functions 
        fst::HGRADtransformVALUE(transformed_value_of_basis_at_cub_points,
                                 value_of_basis_at_cub_points);
        
        // multiply with weighted measure
        fst::multiplyMeasure(weighted_transformed_value_of_basis_at_cub_points,
                             weighted_measure,
                             transformed_value_of_basis_at_cub_points);
        
        // compute mass matrices
        fst::integrate(mass_matrices,
                       transformed_value_of_basis_at_cub_points,
                       weighted_transformed_value_of_basis_at_cub_points);

        ExecSpaceType().fence();

        /*******************  STOP COMPUTATION ***********************/


        /******************* START COMPARISON ***********************/
        auto mass_matrices_host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), mass_matrices);
        auto stiffness_matrices_host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), stiffness_matrices);

        std::string basedir = "./testdata";
        for (auto cid=0;cid<num_cells-1;++cid) {
          std::stringstream namestream;
          std::string filename;
          namestream <<  basedir << "/mass_TET_FEM_P1" << "_" << "0" << cid+1 << ".dat";
          namestream >> filename;

          *outStream << "\nCell ID : " << cid << "  mass matrix comparing with " << filename << "\n\n";

          std::ifstream massfile(&filename[0]);
          if (massfile.is_open()) {
            const auto mass_matrix_cell = Kokkos::subdynrankview(mass_matrices_host, cid, Kokkos::ALL(), Kokkos::ALL());
            errorFlag += compareToAnalytic(massfile,
                                           mass_matrix_cell, 
                                           1e-10, 
                                           verbose);
            massfile.close();
          } else {
            errorFlag = -1;
            INTREPID2_TEST_FOR_EXCEPTION( true, std::runtime_error,
                                          "Failed to open a file" );
          }

          namestream.clear();
          namestream << basedir << "/stiff_TET_FEM_P1" << "_" << "0" << cid+1 << ".dat";
          namestream >> filename;

          *outStream << "\nCell ID : " << cid << "  stiffness matrix comparing with " << filename << "\n\n";

          std::ifstream stifffile(&filename[0]);
          if (stifffile.is_open()) {
            const auto stiffness_matrix_cell = Kokkos::subdynrankview(stiffness_matrices_host, cid, Kokkos::ALL(), Kokkos::ALL());
            errorFlag += compareToAnalytic(stifffile,
                                           stiffness_matrix_cell,
                                           1e-10,
                                           verbose);
          } else {
            errorFlag = -1;
            INTREPID2_TEST_FOR_EXCEPTION( true, std::runtime_error,
                                          "Failed to open a file" );
          }
        }

        /******************* STOP COMPARISON ***********************/

        *outStream << "\n";
      }
      catch (std::logic_error &err) {
        *outStream << "UNEXPECTED ERROR !!! ----------------------------------------------------------\n";
        *outStream << err.what() << '\n';
        *outStream << "-------------------------------------------------------------------------------" << "\n\n";
        errorFlag = -1000;
      };
      
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
