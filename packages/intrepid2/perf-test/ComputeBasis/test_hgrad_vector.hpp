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

/** \file test_01.cpp
\brief  Performance test comparing dynrankview overhead
\author Created by Kyungjoo Kim.
*/

#include "Intrepid2_config.h"
#include "Intrepid2_Types.hpp"
#include "Intrepid2_Utils.hpp"


#include "Intrepid2_HGRAD_QUAD_C1_FEM.hpp"
#include "Intrepid2_Cubature.hpp"
#include "Intrepid2_DefaultCubatureFactory.hpp"

#include "Intrepid2_CellTools.hpp"
#include "Intrepid2_FunctionSpaceTools.hpp"


#include "Teuchos_oblackholestream.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_ScalarTraits.hpp"

#include "KokkosKernels_Vector.hpp"

using namespace KokkosKernels;

namespace Intrepid2 {
  
  namespace Test {

    template<typename VectorTagType, typename DeviceSpaceType>
    int ComputeBasis_HGRAD_Vector(const ordinal_type nworkset,
                                  const ordinal_type C,
                                  const ordinal_type order,
                                  const bool verbose) {
      typedef Vector<VectorTagType> VectorType;

      typedef typename VectorTagType::value_type ValueType;
      constexpr int VectorLength = VectorTagType::length;

      Teuchos::RCP<std::ostream> verboseStream;
      Teuchos::oblackholestream bhs; // outputs nothing

      if (verbose) 
        verboseStream = Teuchos::rcp(&std::cout, false);
      else
        verboseStream = Teuchos::rcp(&bhs,       false);

      Teuchos::oblackholestream oldFormatState;
      oldFormatState.copyfmt(std::cout);

      typedef typename
        Kokkos::Impl::is_space<DeviceSpaceType>::host_mirror_space::execution_space HostSpaceType ;

      *verboseStream << "DeviceSpace::  "; DeviceSpaceType::print_configuration(*verboseStream, false);
      *verboseStream << "HostSpace::    ";   HostSpaceType::print_configuration(*verboseStream, false);
      *verboseStream << "VectorLength::  " << (VectorLength) << "\n";

      using BasisTypeHost = Basis_HGRAD_HEX_C1_FEM<HostSpaceType,ValueType,ValueType>;
      using ImplBasisType = Impl::Basis_HGRAD_HEX_C1_FEM;
      using range_type = Kokkos::pair<ordinal_type,ordinal_type>;

      constexpr size_t LLC_CAPACITY = 32*1024*1024;
      Intrepid2::Test::Flush<LLC_CAPACITY,DeviceSpaceType> flush;
      
      Kokkos::Impl::Timer timer;
      double t_vectorize = 0;
      int errorFlag = 0;

      BasisTypeHost hostBasis;
      const auto cellTopo = hostBasis.getBaseCellTopology();

      auto cubature = DefaultCubatureFactory::create<DeviceSpaceType,ValueType,ValueType>(cellTopo, order);

      const ordinal_type 
        numCells = C,
        numCellsAdjusted = C/VectorLength + (C%VectorLength > 0),
        numVerts = cellTopo.getVertexCount(),
        numDofs = hostBasis.getCardinality(),
        numPoints = cubature->getNumPoints(), 
        spaceDim = cubature->getDimension();

      Kokkos::DynRankView<ValueType,HostSpaceType> dofCoordsHost("dofCoordsHost", numDofs, spaceDim);
      hostBasis.getDofCoords(dofCoordsHost);
      const auto refNodesHost = Kokkos::subview(dofCoordsHost, range_type(0, numVerts), Kokkos::ALL());
      
      // pertub nodes
      Kokkos::DynRankView<VectorType,HostSpaceType> worksetCellsHost("worksetCellsHost", numCellsAdjusted, numVerts, spaceDim);
      for (ordinal_type cell=0;cell<numCells;++cell) {
        for (ordinal_type i=0;i<numVerts;++i)
          for (ordinal_type j=0;j<spaceDim;++j) {
            ValueType val = (rand()/(RAND_MAX + 1.0))*0.2 -0.1;
            worksetCellsHost(cell/VectorLength, i, j)[cell%VectorLength] = refNodesHost(i, j) + val; 
          }
      }

      auto worksetCells = Kokkos::create_mirror_view(typename DeviceSpaceType::memory_space(), worksetCellsHost);
      Kokkos::deep_copy(worksetCells, worksetCellsHost);

      Kokkos::DynRankView<ValueType,DeviceSpaceType> refPoints("refPoints", numPoints, spaceDim), refWeights("refWeights", numPoints);
      cubature->getCubature(refPoints, refWeights);

      std::cout
        << "===============================================================================\n" 
        << " Performance Test evaluating ComputeBasis \n"
        << " # of workset = " << nworkset << "\n" 
        << " Test Array Structure (C,F,P,D) = " << numCells << ", " << numDofs << ", " << numPoints << ", " << spaceDim << "\n"
        << "===============================================================================\n";

      *verboseStream
        << "\n"
        << "===============================================================================\n"
        << "TEST 1: evaluateFields vector version\n"
        << "===============================================================================\n";
      
      try {
        Kokkos::DynRankView<ValueType,DeviceSpaceType> 
          refBasisValues("refBasisValues", numDofs, numPoints),
          refBasisGrads ("refBasisGrads",  numDofs, numPoints, spaceDim);
        
        ImplBasisType::getValues<DeviceSpaceType>(refBasisValues, refPoints, OPERATOR_VALUE);
        ImplBasisType::getValues<DeviceSpaceType>(refBasisGrads,  refPoints, OPERATOR_GRAD);
        
        const ordinal_type ibegin = -3;

        // testing vertical approach
        {          
          Kokkos::DynRankView<VectorType,DeviceSpaceType> 
            weightedBasisValues("weightedBasisValues", numCellsAdjusted, numDofs, numPoints),
            weightedBasisGrads ("weightedBasisGrads",  numCellsAdjusted, numDofs, numPoints, spaceDim);

          typedef F_hgrad_eval<VectorType,ValueType,DeviceSpaceType> FunctorType;

          using range_policy_type = Kokkos::Experimental::MDRangePolicy
            < DeviceSpaceType, Kokkos::Experimental::Rank<2>, Kokkos::IndexType<ordinal_type> >;
          range_policy_type policy( {                0,         0 },
                                    { numCellsAdjusted, numPoints } );

          FunctorType functor(weightedBasisValues,
                              weightedBasisGrads,
                              refBasisGrads,
                              worksetCells,
                              refWeights,
                              refBasisValues,
                              refBasisGrads);
          
          for (ordinal_type iwork=ibegin;iwork<nworkset;++iwork) {
            flush.run();

            DeviceSpaceType::fence();
            timer.reset();
            
            Kokkos::parallel_for(policy, functor);

            DeviceSpaceType::fence();
            t_vectorize += (iwork >= 0)*timer.seconds();
          }

        }

      } catch (std::exception err) {
        *verboseStream << "UNEXPECTED ERROR !!! ----------------------------------------------------------\n";
        *verboseStream << err.what() << '\n';
        *verboseStream << "-------------------------------------------------------------------------------" << "\n\n";
        errorFlag = -1000;
      }

      std::cout 
        << "TEST HGRAD " 
        << " t_vectorize = " << (t_vectorize/nworkset)
        << std::endl;
      
      if (errorFlag != 0)
        std::cout << "End Result: TEST FAILED\n";
      else
        std::cout << "End Result: TEST PASSED\n";
      
      // reset format state of std::cout
      std::cout.copyfmt(oldFormatState);
      
      return errorFlag;
    }
  } // end of namespace TEST
} // end of namespace Intrepid2
