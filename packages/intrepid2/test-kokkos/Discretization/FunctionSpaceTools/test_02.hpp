// @HEADER
// ************************************************************************
//
//                           Intrepid2 Package
//                 Copyright (2007) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Kyungjoo Kim  (kyukim@sandia.gov), or
//                    Mauro Perego  (mperego@sandia.gov)
//
// ************************************************************************
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
#include "Intrepid2_HCURL_HEX_I1_FEM.hpp"
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
    } catch (std::logic_error err) {                                    \
      ++ncatch;                                                         \
      *outStream << "Expected Error ----------------------------------------------------------------\n"; \
      *outStream << err.what() << '\n';                                 \
      *outStream << "-------------------------------------------------------------------------------" << "\n\n"; \
    };
#define ConstructWithLabel(obj, ...) obj(#obj, __VA_ARGS__)
    
    template<typename ValueType, typename DeviceSpaceType>
    int FunctionSpaceTools_Test02(const bool verbose) {

      Teuchos::RCP<std::ostream> outStream;
      Teuchos::oblackholestream bhs; // outputs nothing

      if (verbose)
        outStream = Teuchos::rcp(&std::cout, false);
      else
        outStream = Teuchos::rcp(&bhs, false);

      Teuchos::oblackholestream oldFormatState;
      oldFormatState.copyfmt(std::cout);

      typedef typename
        Kokkos::Impl::is_space<DeviceSpaceType>::host_mirror_space::execution_space HostSpaceType ;
      
      *outStream << "DeviceSpace::  "; DeviceSpaceType::print_configuration(*outStream, false);
      *outStream << "HostSpace::    ";   HostSpaceType::print_configuration(*outStream, false);
      
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

      typedef CellTools<DeviceSpaceType> ct;
      typedef FunctionSpaceTools<DeviceSpaceType> fst;
      typedef Kokkos::DynRankView<ValueType,DeviceSpaceType> DynRankView;
      
      int errorFlag = 0;
      
      *outStream
        << "\n"
        << "===============================================================================\n"
        << "| TEST 1: correctness of math operations                                      |\n"
        << "===============================================================================\n";
      
      outStream->precision(20);
      
      try {

        DefaultCubatureFactory cub_factory;

        shards::CellTopology cell_topo = shards::getCellTopologyData< shards::Hexahedron<8> >();

        const auto cub_degree = 20;
        auto cub = cub_factory.create<DeviceSpaceType,ValueType,ValueType>(cell_topo, cub_degree);

        const auto space_dim = cub->getDimension();
        const auto num_cub_points = cub->getNumPoints();

        Basis_HCURL_HEX_I1_FEM<DeviceSpaceType> hexBasis;
        const auto num_fields = hexBasis.getCardinality();

        /* Cell geometries and orientations. */
        const auto num_cells = 4;
        const auto num_nodes = 8;

        ValueType hexnodes[] = {
          // hex 0  -- affine
          -1.0, -1.0, -1.0,
          1.0, -1.0, -1.0,
          1.0, 1.0, -1.0,
          -1.0, 1.0, -1.0,
          -1.0, -1.0, 1.0,
          1.0, -1.0, 1.0,
          1.0, 1.0, 1.0,
          -1.0, 1.0, 1.0,
          // hex 1  -- affine
          -3.0, -3.0, 1.0,
          6.0, 3.0, 1.0,
          7.0, 8.0, 0.0,
          -2.0, 2.0, 0.0,
          -3.0, -3.0, 4.0,
          6.0, 3.0, 4.0,
          7.0, 8.0, 3.0,
          -2.0, 2.0, 3.0,
          // hex 2  -- affine
          -3.0, -3.0, 0.0,
          9.0, 3.0, 0.0,
          15.0, 6.1, 0.0,
          3.0, 0.1, 0.0,
          9.0, 3.0, 0.1,
          21.0, 9.0, 0.1,
          27.0, 12.1, 0.1,
          15.0, 6.1, 0.1,
          // hex 3  -- nonaffine
          -2.0, -2.0, 0.0,
          2.0, -1.0, 0.0,
          1.0, 6.0, 0.0,
          -1.0, 1.0, 0.0,
          0.0, 0.0, 1.0,
          1.0, 0.0, 1.0,
          1.0, 1.0, 1.0,
          0.0, 1.0, 1.0
        };

        ValueType edgesigns[] = {
          1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
          1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1,
          -1, -1, -1, 1, 1, 1, -1, -1, 1, 1, -1, 1,
          1, 1, 1, 1, -1, -1, -1, -1, 1, 1, 1, 1
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

        DynRankView ConstructWithLabel( curl_of_basis_at_cub_points,                                 num_fields, num_cub_points, space_dim);
        DynRankView ConstructWithLabel( transformed_curl_of_basis_at_cub_points,          num_cells, num_fields, num_cub_points, space_dim);
        DynRankView ConstructWithLabel( weighted_transformed_curl_of_basis_at_cub_points, num_cells, num_fields, num_cub_points, space_dim);
        DynRankView ConstructWithLabel( stiffness_matrices,                               num_cells, num_fields, num_fields);
        
        DynRankView ConstructWithLabel( value_of_basis_at_cub_points,                                 num_fields, num_cub_points, space_dim);
        DynRankView ConstructWithLabel( transformed_value_of_basis_at_cub_points,          num_cells, num_fields, num_cub_points, space_dim);
        DynRankView ConstructWithLabel( weighted_transformed_value_of_basis_at_cub_points, num_cells, num_fields, num_cub_points, space_dim);
        DynRankView ConstructWithLabel( mass_matrices,                                     num_cells, num_fields, num_fields);

        /******************* START COMPUTATION ***********************/

        // get cubature points and weights
        cub->getCubature(cub_points, cub_weights);

        const Kokkos::DynRankView<const ValueType,Kokkos::LayoutRight,Kokkos::HostSpace> cell_nodes_host (&hexnodes[0],  num_cells, num_nodes, space_dim);
        const Kokkos::DynRankView<const ValueType,Kokkos::LayoutRight,Kokkos::HostSpace> field_signs_host(&edgesigns[0], num_cells, num_fields);

        Kokkos::deep_copy( cell_nodes,  cell_nodes_host  );
        Kokkos::deep_copy( field_signs, field_signs_host );

        // compute geometric cell information
        ct::setJacobian(jacobian, cub_points, cell_nodes, cell_topo);
        ct::setJacobianInv(jacobian_inv, jacobian);
        ct::setJacobianDet(jacobian_det, jacobian);

        // compute weighted measure
        fst::computeCellMeasure(weighted_measure, jacobian_det, cub_weights);

        // ** Computing stiffness matrices:
        hexBasis.getValues(curl_of_basis_at_cub_points, cub_points, OPERATOR_CURL);

        // transform curls of basis functions 
        fst::HCURLtransformCURL(transformed_curl_of_basis_at_cub_points,
                                jacobian,
                                jacobian_det,
                                curl_of_basis_at_cub_points);

        // multiply with weighted measure
        fst::multiplyMeasure(weighted_transformed_curl_of_basis_at_cub_points,
                             weighted_measure,
                             transformed_curl_of_basis_at_cub_points);
        
        // we can apply the field signs to the basis function arrays, or after the fact, see below
        fst::applyFieldSigns(transformed_curl_of_basis_at_cub_points, field_signs);
        fst::applyFieldSigns(weighted_transformed_curl_of_basis_at_cub_points, field_signs);
        
        // compute stiffness matrices
        fst::integrate(stiffness_matrices,
                       transformed_curl_of_basis_at_cub_points,
                       weighted_transformed_curl_of_basis_at_cub_points);

        // ** Computing mass matrices:
        hexBasis.getValues(value_of_basis_at_cub_points, cub_points, OPERATOR_VALUE);

        // transform values of basis functions 
        fst::HCURLtransformVALUE(transformed_value_of_basis_at_cub_points,
                                 jacobian_inv,
                                 value_of_basis_at_cub_points);
        
        // multiply with weighted measure
        fst::multiplyMeasure(weighted_transformed_value_of_basis_at_cub_points,
                             weighted_measure,
                             transformed_value_of_basis_at_cub_points);
        
        // compute mass matrices
        fst::integrate(mass_matrices,
                       transformed_value_of_basis_at_cub_points,
                       weighted_transformed_value_of_basis_at_cub_points);

        // apply field signs (after the fact, as a post-processing step)
        fst::applyLeftFieldSigns(mass_matrices, field_signs);
        fst::applyRightFieldSigns(mass_matrices, field_signs);

        /*******************  STOP COMPUTATION ***********************/


        /******************* START COMPARISON ***********************/
        std::string basedir = "../testdata";
        for (auto cid=0;cid<num_cells-1;++cid) {
          std::stringstream namestream;
          std::string filename;
          namestream <<  basedir << "/mass_HCURL_HEX_I1_FEM" << "_" << "0" << cid+1 << ".dat";
          namestream >> filename;

          *outStream << "\nCell ID : " << cid << "  mass matrix comparing with " << filename << "\n\n";

          std::ifstream massfile(&filename[0]);
          if (massfile.is_open()) {
            const auto mass_matrix_cell = Kokkos::subdynrankview(mass_matrices, cid, Kokkos::ALL(), Kokkos::ALL());
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
          namestream << basedir << "/stiff_HCURL_HEX_I1_FEM" << "_" << "0" << cid+1 << ".dat";
          namestream >> filename;

          *outStream << "\nCell ID : " << cid << "  stiffness matrix comparing with " << filename << "\n\n";

          std::ifstream stifffile(&filename[0]);
          if (stifffile.is_open()) {
            const auto stiffness_matrix_cell = Kokkos::subdynrankview(stiffness_matrices, cid, Kokkos::ALL(), Kokkos::ALL());
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

        for (auto cid=3;cid<num_cells;++cid) {
          std::stringstream namestream;
          std::string filename;
          namestream <<  basedir << "/mass_fp_HCURL_HEX_I1_FEM" << "_" << "0" << cid+1 << ".dat";
          namestream >> filename;

          *outStream << "\nCell ID : " << cid << "  mass matrix comparing with " << filename << "\n\n";

          std::ifstream massfile(&filename[0]);
          if (massfile.is_open()) {
            const auto mass_matrix_cell = Kokkos::subdynrankview(mass_matrices, cid, Kokkos::ALL(), Kokkos::ALL());
            errorFlag += compareToAnalytic(massfile,
                                           mass_matrix_cell,
                                           1e-4,
                                           verbose,
                                           INTREPID2_UTILS_SCALAR);
            massfile.close();
          } else {
            errorFlag = -1;
            INTREPID2_TEST_FOR_EXCEPTION( true, std::runtime_error,
                                          "Failed to open a file" );
          }
          
          namestream.clear();
          namestream << basedir << "/stiff_fp_HCURL_HEX_I1_FEM" << "_" << "0" << cid+1 << ".dat";
          namestream >> filename;

          *outStream << "\nCell ID : " << cid << "  stiffness matrix comparing with " << filename << "\n\n";

          std::ifstream stifffile(&filename[0]);
          if (stifffile.is_open()) {
            const auto stiffness_matrix_cell = Kokkos::subdynrankview(stiffness_matrices, cid, Kokkos::ALL(), Kokkos::ALL());
            errorFlag += compareToAnalytic(stifffile,
                                           stiffness_matrix_cell,
                                           1e-4,
                                           verbose,
                                           INTREPID2_UTILS_SCALAR);
            stifffile.close();
          } else {
            errorFlag = -1;
            INTREPID2_TEST_FOR_EXCEPTION( true, std::runtime_error,
                                          "Failed to open a file" );
          }
        }

        /******************* STOP COMPARISON ***********************/

        *outStream << "\n";
      }
      catch (std::logic_error err) {
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
