// @HEADER
// *****************************************************************************
//                           Intrepid2 Package
//
// Copyright 2007 NTESS and the Intrepid2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER


/** \file
    \brief  Test of the CellTools class.
    \author Created by P. Bochev, D. Ridzal, K. Peterson and Kyungjoo Kim
*/

#include "Intrepid2_config.h"

#ifdef HAVE_INTREPID2_DEBUG
#define INTREPID2_TEST_FOR_DEBUG_ABORT_OVERRIDE_TO_CONTINUE
#endif

#include "Intrepid2_DefaultCubatureFactory.hpp"

#include "Intrepid2_CellTools.hpp"
#include "Intrepid2_CellTools_Serial.hpp"

#include "Teuchos_oblackholestream.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_ScalarTraits.hpp"


namespace Intrepid2 {

  namespace Test {

    template<typename OutputViewType,
             typename InputViewType,
             typename WorksetViewType,
             typename ImplBasis>      
    struct F_mapToReferenceFrame {
      OutputViewType output_;
      const InputViewType input_;
      const WorksetViewType workset_;
      const typename InputViewType::value_type shellThickness_;


      KOKKOS_INLINE_FUNCTION
      F_mapToReferenceFrame(OutputViewType output,
                            const InputViewType input,
                            const WorksetViewType workset,
                            const typename InputViewType::value_type shellThickness) 
        : output_(output), 
          input_(input), 
          workset_(workset),
          shellThickness_(shellThickness) {}
      
      KOKKOS_INLINE_FUNCTION
      void
      operator()(const ordinal_type cl, const ordinal_type pt) const {
        auto nodes = Kokkos::subview(workset_, cl, Kokkos::ALL(), Kokkos::ALL()); // N,D
        auto in  = Kokkos::subview(input_,  cl, pt, Kokkos::ALL()); // D
        auto out = Kokkos::subview(output_, cl, pt, Kokkos::ALL()); // D

        const auto tol = tolerence<typename OutputViewType::value_type>();        
          
        Impl::CellTools::Serial::mapToReferenceFrame<ImplBasis>(out, in, nodes, tol, shellThickness_);
        }
      };

    template<typename DeviceType,
             typename OutputViewType,
             typename InputViewType,
             typename WorksetViewType>
      static void
      mapToReferenceFrameSerialImpl(      OutputViewType output,
                    const InputViewType input,
                    const WorksetViewType workset,
                    const shards::CellTopology cellTopo,
                    typename InputViewType::value_type shellThickness) {

        using range_policy_type = Kokkos::MDRangePolicy<typename DeviceType::execution_space, Kokkos::Rank<2>, Kokkos::IndexType<ordinal_type> >;
        range_policy_type policy( { 0, 0 }, { input.extent(0), input.extent(1) } );

        OutputViewType cellCenter("cellCenter", cellTopo.getDimension());
        CellTools<DeviceType>::getReferenceCellCenter(cellCenter, cellTopo);
        RealSpaceTools<DeviceType>::clone(output, cellCenter);

        #define KOKKOS_MAP_TO_REF_PARALLEL_FOR(BasisName)      \
          do {                                                                                                                      \
            using FunctorType = F_mapToReferenceFrame<OutputViewType, InputViewType, WorksetViewType, Intrepid2::Impl::BasisName>; \
            Kokkos::parallel_for(policy, FunctorType(output, input, workset, shellThickness));                                                               \
          } while (0)

        switch (cellTopo.getKey()) {
        case shards::Line<2>::key:                KOKKOS_MAP_TO_REF_PARALLEL_FOR(Basis_HGRAD_LINE_C1_FEM); break; 
        case shards::Triangle<3>::key:            KOKKOS_MAP_TO_REF_PARALLEL_FOR(Basis_HGRAD_TRI_C1_FEM); break;
        case shards::Quadrilateral<4>::key:       KOKKOS_MAP_TO_REF_PARALLEL_FOR(Basis_HGRAD_QUAD_C1_FEM); break;
        case shards::Tetrahedron<4>::key:         KOKKOS_MAP_TO_REF_PARALLEL_FOR(Basis_HGRAD_TET_C1_FEM); break;
        case shards::Hexahedron<8>::key:          KOKKOS_MAP_TO_REF_PARALLEL_FOR(Basis_HGRAD_HEX_C1_FEM); break;
        case shards::Wedge<6>::key:               KOKKOS_MAP_TO_REF_PARALLEL_FOR(Basis_HGRAD_WEDGE_C1_FEM); break;
        case shards::Pyramid<5>::key:             KOKKOS_MAP_TO_REF_PARALLEL_FOR(Basis_HGRAD_PYR_C1_FEM); break;

        case shards::Line<3>::key:                KOKKOS_MAP_TO_REF_PARALLEL_FOR(Basis_HGRAD_LINE_C2_FEM); break;   
        case shards::Triangle<6>::key:            KOKKOS_MAP_TO_REF_PARALLEL_FOR(Basis_HGRAD_TRI_C2_FEM); break;
        case shards::Quadrilateral<8>::key:       KOKKOS_MAP_TO_REF_PARALLEL_FOR(Basis_HGRAD_QUAD_DEG2_FEM<true>); break;
        case shards::Quadrilateral<9>::key:       KOKKOS_MAP_TO_REF_PARALLEL_FOR(Basis_HGRAD_QUAD_DEG2_FEM<false>); break;
        case shards::Tetrahedron<10>::key:        KOKKOS_MAP_TO_REF_PARALLEL_FOR(Basis_HGRAD_TET_C2_FEM); break;
        case shards::Tetrahedron<11>::key:        KOKKOS_MAP_TO_REF_PARALLEL_FOR(Basis_HGRAD_TET_COMP12_FEM); break;
        case shards::Hexahedron<20>::key:         KOKKOS_MAP_TO_REF_PARALLEL_FOR(Basis_HGRAD_HEX_DEG2_FEM<true>); break;
        case shards::Hexahedron<27>::key:         KOKKOS_MAP_TO_REF_PARALLEL_FOR(Basis_HGRAD_HEX_DEG2_FEM<false>); break;
        case shards::Wedge<15>::key:              KOKKOS_MAP_TO_REF_PARALLEL_FOR(Basis_HGRAD_WEDGE_DEG2_FEM<true>); break;
        case shards::Wedge<18>::key:              KOKKOS_MAP_TO_REF_PARALLEL_FOR(Basis_HGRAD_WEDGE_DEG2_FEM<false>); break;
        case shards::Pyramid<13>::key:            KOKKOS_MAP_TO_REF_PARALLEL_FOR(Basis_HGRAD_PYR_I2_FEM); break;

        case shards::Beam<2>::key:                KOKKOS_MAP_TO_REF_PARALLEL_FOR(Basis_HGRAD_LINE_C1_FEM); break;
        case shards::Beam<3>::key:                KOKKOS_MAP_TO_REF_PARALLEL_FOR(Basis_HGRAD_LINE_C2_FEM); break;
        case shards::ShellLine<2>::key:           KOKKOS_MAP_TO_REF_PARALLEL_FOR(Basis_HGRAD_LINE_C1_FEM); break;
        case shards::ShellLine<3>::key:           KOKKOS_MAP_TO_REF_PARALLEL_FOR(Basis_HGRAD_LINE_C2_FEM); break;
        case shards::ShellTriangle<3>::key:       KOKKOS_MAP_TO_REF_PARALLEL_FOR(Basis_HGRAD_TRI_C1_FEM); break;
        case shards::ShellTriangle<6>::key:       KOKKOS_MAP_TO_REF_PARALLEL_FOR(Basis_HGRAD_TRI_C2_FEM); break;
        case shards::ShellQuadrilateral<4>::key:  KOKKOS_MAP_TO_REF_PARALLEL_FOR(Basis_HGRAD_QUAD_C1_FEM); break;
        case shards::ShellQuadrilateral<8>::key:  KOKKOS_MAP_TO_REF_PARALLEL_FOR(Basis_HGRAD_QUAD_DEG2_FEM<true>); break;
        case shards::ShellQuadrilateral<9>::key:  KOKKOS_MAP_TO_REF_PARALLEL_FOR(Basis_HGRAD_QUAD_DEG2_FEM<false>); break;
        
        default: {
          INTREPID2_TEST_FOR_EXCEPTION( true, std::invalid_argument,
                                        ">>> ERROR (Intrepid2::CellTools::createHGradBasis): Cell topology not supported.");
        }
        }
      }



#define INTREPID2_TEST_ERROR_EXPECTED( S )                              \
    try {                                                               \
      S ;                                                               \
    }                                                                   \
    catch (std::logic_error &err) {                                      \
      *outStream << "Expected Error ----------------------------------------------------------------\n"; \
      *outStream << err.what() << '\n';                                 \
      *outStream << "-------------------------------------------------------------------------------" << "\n\n"; \
    };
#define ConstructWithLabel(obj, ...) obj(#obj, __VA_ARGS__)
    
    template<typename ValueType, typename DeviceType>
    int CellTools_Test04(const bool verbose) {

      Teuchos::RCP<std::ostream> outStream;
      Teuchos::oblackholestream bhs; // outputs nothing

      if (verbose)
        outStream = Teuchos::rcp(&std::cout, false);
      else
        outStream = Teuchos::rcp(&bhs,       false);

      Teuchos::oblackholestream oldFormatState;
      oldFormatState.copyfmt(std::cout);

      *outStream                                                       
        << "===============================================================================\n"
        << "|                                                                             |\n"
        << "|                              Unit Test CellTools                            |\n"
        << "|                                                                             |\n"
        << "|     1) Mapping to and from reference cells with base and extended topologies|\n"
        << "|        using default initial guesses when computing the inverse F^{-1}      |\n"
        << "|     2) Repeat all tests from 1) using user-defined initial guess for F^{-1} |\n"
        << "|     3) Exception testing                                                    |\n"
        << "|                                                                             |\n"
        << "|  Questions? Contact  Pavel Bochev (pbboche@sandia.gov)                      |\n"
        << "|                      Denis Ridzal (dridzal@sandia.gov), or                  |\n"
        << "|                      Kara Peterson(kjpeter@sandia.gov), or                  |\n"
        << "|                      Kyungjoo Kim (kyukim@sandia.gov)                       |\n"
        << "|                                                                             |\n"
        << "|  Intrepid's website: http://trilinos.sandia.gov/packages/intrepid           |\n"
        << "|  Trilinos website:   http://trilinos.sandia.gov                             |\n"
        << "|                                                                             |\n"
        << "===============================================================================\n";

      using ct = CellTools<DeviceType>;
      using DynRankView = Kokkos::DynRankView<ValueType,DeviceType>;

      const ValueType tol = tolerence<ValueType>()*100.0;

      int errorFlag  = 0;

      // Collect all supported cell topologies
      std::vector<shards::CellTopology> standardBaseTopologies;
      shards::getTopologies(standardBaseTopologies, 4, shards::ALL_CELLS, shards::ALL_TOPOLOGIES);
      const auto topoSize = standardBaseTopologies.size();

      for (auto testCase=0;testCase<3;++testCase) {
        try {
          switch (testCase) {
          case 0: {
            *outStream
              << "\n"
              << "===============================================================================\n" 
              << "| Test 1: computing F(x) and F^{-1}(x) using default initial guesses.         |\n" 
              << "===============================================================================\n\n";
            break;
          }
          case 1: {
            *outStream
              << "\n"
              << "===============================================================================\n" 
              << "| Test 2: computing F(x) and F^{-1}(x) using user-defined initial guess.      |\n" 
              << "===============================================================================\n\n";
            break;
          }
          case 2: {
            *outStream
              << "\n"
              << "===============================================================================\n" 
              << "| Test 3: computing F(x) and F^{-1}(x) using serial implementation.            |\n" 
              << "===============================================================================\n\n";
            break;
          }
          } 

          /*
           *  Test summary:
           *
           *    A reference point set is mapped to physical frame and then back to reference frame.
           *    Test passes if the final set of points matches the first set of points. The cell workset
           *    is generated by perturbing randomly the worksetCell of a reference cell with the specified 
           *    cell topology. 
           *
           */
          DefaultCubatureFactory cubFactory;   
          
          // Initialize testing env.
          const auto numCells = 10;
          const auto testAccuracy = 4;

          // Loop over cell topologies, make cell workset for each one by perturbing the worksetCell & test methods
          for (size_type topoOrd=0;topoOrd<topoSize;++topoOrd) {
            const auto cell = standardBaseTopologies[topoOrd];
            
            if (!ct::hasReferenceCell(cell))
              continue;
      
            // 1.   Define a single reference point set using cubature factory with order 4 cubature
            const auto cellCubature = cubFactory.create<DeviceType,ValueType,ValueType>(cell, testAccuracy);
            const auto cubDim = cellCubature->getDimension();
            const auto cubNumPoints = cellCubature->getNumPoints();
            const auto cellDim  = cell.getDimension();
            
            DynRankView ConstructWithLabel(cubPoints,  cubNumPoints, cellDim );
            auto cubDimRange = std::pair<int,int>(0,cubDim);
            DynRankView ConstructWithLabel(cubWeights, cubNumPoints );
            cellCubature->getCubature(Kokkos::subview(cubPoints, Kokkos::ALL(), cubDimRange), cubWeights);

            // 2.   Define a cell workset by perturbing the worksetCell of the reference cell with the specified topology
            // 2.1  Resize dimensions of the rank-3 (C,N,D) cell workset array for the current topology
            const auto numNodes = cell.getNodeCount();

            DynRankView ConstructWithLabel(worksetCell, numCells, numNodes, cellDim);
            
            // 2.2  Copy worksetCell of the reference cell with the same topology to temp rank-2 (N,D) array
            DynRankView ConstructWithLabel(refCellNodes, numNodes, cellDim);
            ct::getReferenceSubcellNodes(refCellNodes, cellDim, 0, cell);
            
            // create mirror host views
            auto hRefCellNodes = Kokkos::create_mirror_view_and_copy(
                Kokkos::HostSpace(), refCellNodes);

            auto hWorksetCell = Kokkos::create_mirror_view(worksetCell);

            // 2.3  Create randomly perturbed version of the reference cell and save in the cell workset array
            for (auto cellOrd=0;cellOrd<numCells;++cellOrd) {
              // Move vertices +/-0.03125 along their axes. Gives nontrivial cells for base and extended topologies 
              for (size_type nodeOrd=0;nodeOrd<numNodes;++nodeOrd) 
                for(size_type d=0;d<cellDim;++d) {
                  const auto delta = Teuchos::ScalarTraits<double>::random()/32.0;
                  hWorksetCell(cellOrd, nodeOrd, d) = hRefCellNodes(nodeOrd, d) + delta;
                } 
            }
            Kokkos::deep_copy(worksetCell,hWorksetCell);
            
            /* 
             * 3.1 Test 1: single point set to single physical cell: map ref. point set in rank-2 (P,D) array
             *      to a physical point set in rank-2 (P,D) array for a specified cell ordinal. Use the cub.
             *      points array for this test. Resize physPoints and controlPoints to rank-2 (P,D) arrays.
             */
            DynRankView ConstructWithLabel(physPoints,    numCells, cubNumPoints, cellDim );
            DynRankView ConstructWithLabel(controlPoints, numCells, cubNumPoints, cellDim );
            
            *outStream 
              << " Mapping a set of " << cubNumPoints << " points to one cell in a workset of " << numCells << " " 
              << cell.getName() << " cells. \n";

            // Forward map:: requires cell ordinal; input (P,D) -> output (C,P,D)
            ct::mapToPhysicalFrame(physPoints, cubPoints, worksetCell, cell);
            
            const auto key = cell.getBaseKey();
            const bool isShellorBeam = (key == shards::Beam<>::key) || (key == shards::ShellLine<>::key) || (key == shards::ShellTriangle<>::key) || (key == shards::ShellQuadrilateral<>::key);
            ValueType shellThickness = (isShellorBeam) ? 0.01 : -1.0;
            
            // Inverse map: requires cell ordinal; input (C,P,D) -> output (C,P,D)
            switch (testCase) {
            case 0: {
              ct::mapToReferenceFrame(controlPoints, physPoints, worksetCell, cell, shellThickness);
              break;
            }
            case 1: {
              DynRankView ConstructWithLabel(initGuess, numCells, cubNumPoints, cellDim );
              // cubPoints become the initial guess
              RealSpaceTools<DeviceType>::clone(initGuess, cubPoints);
              ct::mapToReferenceFrameInitGuess(controlPoints, initGuess, physPoints, worksetCell, cell, shellThickness);
              break;
            }
            case 2: {
              mapToReferenceFrameSerialImpl<DeviceType>(controlPoints, physPoints, worksetCell, cell, shellThickness);
              break;
            }
            }
            
            // copy to mirror host views
            auto hControlPoints = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), controlPoints);
            auto hCubPoints = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), cubPoints);

            // Points in controlPoints should match the originals in cubPoints up to a tolerance
            for (auto cellOrd=0;cellOrd<numCells;++cellOrd) 
              for (auto pt=0;pt<cubNumPoints;++pt) 
                for (size_type d=0;d<cellDim;++d) 
                  if( std::abs( hControlPoints(cellOrd, pt, d) - hCubPoints(pt, d) ) > tol ) {
                    errorFlag++;
                    *outStream
                      << std::setw(70) << "^^^^----FAILURE!" << "\n"
                      << " Mapping a single point set to a single physical cell in a workset failed for: \n"
                      << "                    Cell Topology = " << cell.getName() << "\n"
                      << " Physical cell ordinal in workset = " << cellOrd << "\n"
                      << "          Reference point ordinal = " << std::setprecision(12) << pt << "\n"
                      << "    At reference point coordinate = " << std::setprecision(12) << d << "\n"
                      << "                   Original value = " << hCubPoints(pt, d) << "\n"
                      << "                     F^{-1}F(P_d) = " << hControlPoints(cellOrd, pt, d) <<"\n";
                  }
          }
        } catch (std::logic_error &err) {
          *outStream << err.what() << "\n";
          errorFlag = -1000;
        }
      }

      if (errorFlag != 0)
        std::cout << "End Result: TEST FAILED\n";
      else
        std::cout << "End Result: TEST PASSED\n";

      // reset format state of std::cout
      std::cout.copyfmt(oldFormatState);

      return errorFlag;
    }
  } // end test
} // end intrepid2
