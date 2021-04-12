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


/** \file
    \brief  Test of the CellTools class.
    \author Kyungjoo Kim
*/

#include "Intrepid2_config.h"
#include "Intrepid2_CellTopologyTags.hpp"

#include "Intrepid2_DefaultCubatureFactory.hpp"

#include "Teuchos_oblackholestream.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_ScalarTraits.hpp"

namespace Intrepid2 {
  
  namespace Test {
#define INTREPID2_TEST_ERROR_EXPECTED( S )                              \
    try {                                                               \
      S ;                                                               \
    }                                                                   \
    catch (std::logic_error &err) {                                      \
      *outStream << "Expected Error ----------------------------------------------------------------\n"; \
      *outStream << err.what() << '\n';                                 \
      *outStream << "-------------------------------------------------------------------------------" << "\n\n"; \
    };                                                                  


#define INTREPID2_TEST_CHECK_POINT_INCLUSION(offset, expected, shtopo, celltag) \
    {                                                                   \
      typedef shtopo shardsTopology;                                    \
      typedef celltag CellTopologyTag;                                  \
                                                                        \
      const auto order = 3;                                             \
      const auto cellTopo = shards::CellTopology(shards::getCellTopologyData< shardsTopology >()); \
      auto cub = DefaultCubatureFactory::create<DeviceType,ValueType,ValueType>(cellTopo, order); \
      const auto P = cub->getNumPoints();                               \
      const auto D = 3;                                                 \
                                                                        \
      Kokkos::DynRankView<ValueType,DeviceType> pts("pts", P, D);  \
      Kokkos::DynRankView<ValueType,DeviceType> wts("wts", P);     \
      Kokkos::DynRankView<int,DeviceType> check("check", P);       \
                                                                        \
      cub->getCubature(pts, wts);                                       \
                                                                        \
      Kokkos::RangePolicy<DeviceType> policy(0, P);                \
      typedef F_checkPointInclusion<CellTopologyTag,decltype(check),decltype(pts)> FunctorType; \
      Kokkos::parallel_for(policy, FunctorType(offset, check, pts));    \
                                                                        \
      auto check_host = Kokkos::create_mirror_view(typename HostSpaceType::memory_space(), check); \
      Kokkos::deep_copy(check_host, check);                             \
                                                                        \
      for (ordinal_type i=0;i<P;++i) {                                  \
        const double diff = std::abs(check_host(i) - expected);         \
        if (diff > tol) {                                               \
          *outStream << "Error : checkPointInclusion at ("              \
                     << i                                               \
                     << ") with diff = " << diff << "\n";               \
          errorFlag++;                                                  \
        }                                                               \
      }                                                                 \
    }                                                                   \
      
    template<typename cellTopologyTagType,
             typename OutputViewType,
             typename inputViewType>
    struct F_checkPointInclusion {
      double _offset;
      OutputViewType _output;
      inputViewType _input;

      KOKKOS_INLINE_FUNCTION
      F_checkPointInclusion(const double offset_, 
                            OutputViewType output_,
                            inputViewType input_)
        : _offset(offset_), 
          _output(output_), 
          _input(input_) {}

      KOKKOS_INLINE_FUNCTION
      void
      operator()(const ordinal_type i) const {
        const auto in = Kokkos::subview(_input,i,Kokkos::ALL());
        for (int k=0;k<3;++k) in(k)+=_offset;
        const auto check = cellTopologyTagType::checkPointInclusion(in, 0.0);
        
        _output(i) = check;        
      }
    };
        
    template<typename ValueType, typename DeviceType>
    int CellTools_Test07(const bool verbose) {

      using ExecSpaceType = typename DeviceType::execution_space;

      Teuchos::RCP<std::ostream> outStream;
      Teuchos::oblackholestream bhs; // outputs nothing

      if (verbose)
        outStream = Teuchos::rcp(&std::cout, false);
      else
        outStream = Teuchos::rcp(&bhs,       false);

      Teuchos::oblackholestream oldFormatState;
      oldFormatState.copyfmt(std::cout);

      typedef typename
        Kokkos::Impl::is_space<DeviceType>::host_mirror_space::execution_space HostSpaceType ;

      *outStream << "DeviceSpace::  ";   ExecSpaceType::print_configuration(*outStream, false);
      *outStream << "HostSpace::    ";   HostSpaceType::print_configuration(*outStream, false);
      
      *outStream
        << "===============================================================================\n"
        << "|                                                                             |\n"
        << "|                              Unit Test CellTools                            |\n"
        << "|                                                                             |\n"
        << "|     1) check point inclusion and cell topology tag tests                    |\n"
        << "|                                                                             |\n"
        << "===============================================================================\n";
  
      const ValueType tol = tolerence()*100.0;

      int errorFlag = 0;
      
      try {
        
        *outStream
          << "\n"
          << "===============================================================================\n" 
          << "| Test 1: test cubature points\n"
          << "===============================================================================\n\n";

        {
          double offset = 0.0;
          INTREPID2_TEST_CHECK_POINT_INCLUSION(offset, true, shards::Line<>,          Impl::Line<2>);
          
          INTREPID2_TEST_CHECK_POINT_INCLUSION(offset, true, shards::Triangle<>,      Impl::Triangle<3>);
          INTREPID2_TEST_CHECK_POINT_INCLUSION(offset, true, shards::Quadrilateral<>, Impl::Quadrilateral<4>);
          
          INTREPID2_TEST_CHECK_POINT_INCLUSION(offset, true, shards::Tetrahedron<>,   Impl::Tetrahedron<4>);
          INTREPID2_TEST_CHECK_POINT_INCLUSION(offset, true, shards::Hexahedron<>,    Impl::Hexahedron<8>);
          
          INTREPID2_TEST_CHECK_POINT_INCLUSION(offset, true, shards::Pyramid<>,       Impl::Pyramid<5>);
          INTREPID2_TEST_CHECK_POINT_INCLUSION(offset, true, shards::Wedge<>,         Impl::Wedge<6>);
        }
        {
          double offset = 3.0;
          INTREPID2_TEST_CHECK_POINT_INCLUSION(offset, false, shards::Line<>,          Impl::Line<2>);
          
          INTREPID2_TEST_CHECK_POINT_INCLUSION(offset, false, shards::Triangle<>,      Impl::Triangle<3>);
          INTREPID2_TEST_CHECK_POINT_INCLUSION(offset, false, shards::Quadrilateral<>, Impl::Quadrilateral<4>);
          
          INTREPID2_TEST_CHECK_POINT_INCLUSION(offset, false, shards::Tetrahedron<>,   Impl::Tetrahedron<4>);
          INTREPID2_TEST_CHECK_POINT_INCLUSION(offset, false, shards::Hexahedron<>,    Impl::Hexahedron<8>);
          
          INTREPID2_TEST_CHECK_POINT_INCLUSION(offset, false, shards::Pyramid<>,       Impl::Pyramid<5>);
          INTREPID2_TEST_CHECK_POINT_INCLUSION(offset, false, shards::Wedge<>,         Impl::Wedge<6>);
        }

      } catch (std::logic_error &err) {
        //============================================================================================//
        // Wrap up test: check if the test broke down unexpectedly due to an exception                //
        //============================================================================================//
        *outStream << err.what() << "\n";
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
    

