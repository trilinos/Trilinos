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

#include "Intrepid2_config.h"

#ifdef HAVE_INTREPID2_DEBUG
#define INTREPID2_TEST_FOR_DEBUG_ABORT_OVERRIDE_TO_CONTINUE
#endif

#include "Intrepid2_Basis.hpp"
#include "Intrepid2_HVOL_C0_FEM.hpp"

#include "Teuchos_oblackholestream.hpp"
#include "Teuchos_RCP.hpp"

namespace Intrepid2 {

  namespace Test {
    
#define INTREPID2_TEST_ERROR_EXPECTED( S )                              \
    try {                                                               \
      ++nthrow;                                                         \
      S ;                                                               \
    }                                                                   \
    catch (std::exception err) {                                        \
      ++ncatch;                                                         \
      *outStream << "Expected Error ----------------------------------------------------------------\n"; \
      *outStream << err.what() << '\n';                                 \
      *outStream << "-------------------------------------------------------------------------------" << "\n\n"; \
    }
    
    template<typename ValueType, typename DeviceSpaceType>
    int BasisConst_Test01(const bool verbose) {

      Teuchos::RCP<std::ostream> outStream;
      Teuchos::oblackholestream bhs; // outputs nothing

      if (verbose)
        outStream = Teuchos::rcp(&std::cout, false);
      else
        outStream = Teuchos::rcp(&bhs,       false);

      Teuchos::oblackholestream oldFormatState;
      oldFormatState.copyfmt(std::cout);

      typedef typename
        Kokkos::Impl::is_space<DeviceSpaceType>::host_mirror_space::execution_space HostSpaceType ;

      *outStream << "DeviceSpace::  "; DeviceSpaceType::print_configuration(*outStream, false);
      *outStream << "HostSpace::    ";   HostSpaceType::print_configuration(*outStream, false);

      *outStream
        << "===============================================================================\n"
        << "|                                                                             |\n"
        << "|                 Unit Test (Basis_HVOL_C0_FEM)                              |\n"
        << "|                                                                             |\n"
        << "|  Intrepid's website: http://trilinos.sandia.gov/packages/intrepid           |\n"
        << "|  Trilinos website:   http://trilinos.sandia.gov                             |\n"
        << "|                                                                             |\n"
        << "===============================================================================\n";

      typedef Kokkos::DynRankView<ValueType,DeviceSpaceType> DynRankView;
#define ConstructWithLabel(obj, ...) obj(#obj, __VA_ARGS__)

      const ValueType tol = tolerence();
      int errorFlag = 0;

      // for virtual function, value and point types are declared in the class
      typedef ValueType outputValueType;
      typedef ValueType pointValueType;

      const shards::CellTopology cells[6] = {
        shards::CellTopology(shards::getCellTopologyData<shards::Quadrilateral<4> >() ),
        shards::CellTopology(shards::getCellTopologyData<shards::Hexahedron<8> >() ),
        shards::CellTopology(shards::getCellTopologyData<shards::Pyramid<5> >() ),
        shards::CellTopology(shards::getCellTopologyData<shards::Tetrahedron<4> >() ),
        shards::CellTopology(shards::getCellTopologyData<shards::Triangle<3> >() ),
        shards::CellTopology(shards::getCellTopologyData<shards::Wedge<6> >() )
      };
      
      *outStream
        << "\n"
        << "===============================================================================\n"
        << "| TEST 1: Basis creation, exceptions tests                                    |\n"
        << "===============================================================================\n";

      try {
        for (ordinal_type id=0;id<6;++id) {
          Basis_HVOL_C0_FEM<DeviceSpaceType,outputValueType,pointValueType> basis(cells[id]);
        }

      } catch (std::logic_error err){
        *outStream << "UNEXPECTED ERROR !!! ----------------------------------------------------------\n";
        *outStream << err.what() << '\n';
        *outStream << "-------------------------------------------------------------------------------" << "\n\n";
        errorFlag = -1000;
      } 

      *outStream
        << "\n"
        << "===============================================================================\n"
        << "| TEST 2: correctness of tag to enum and enum to tag lookups                  |\n"
        << "===============================================================================\n";

      try {
        for (ordinal_type id=0;id<6;++id) {
          Basis_HVOL_C0_FEM<DeviceSpaceType,outputValueType,pointValueType> basis(cells[id]);

          const ordinal_type numFields = basis.getCardinality();
          const auto allTags = basis.getAllDofTags();
          
          // Loop over all tags, lookup the associated dof enumeration and then lookup the tag again
          const ordinal_type dofTagSize = allTags.extent(0);
          for (ordinal_type i=0;i<dofTagSize;++i) {
            const ordinal_type bfOrd = basis.getDofOrdinal(allTags(i,0), allTags(i,1), allTags(i,2));
            
            const auto myTag = basis.getDofTag(bfOrd);
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
            const auto myTag  = basis.getDofTag(bfOrd);
            const ordinal_type myBfOrd = basis.getDofOrdinal(myTag(0), myTag(1), myTag(2));
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
        }
      } catch (std::logic_error err){
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

      outStream -> precision(20);

      try {
        
        for (ordinal_type id=0;id<6;++id) {
          Basis_HVOL_C0_FEM<DeviceSpaceType,outputValueType,pointValueType> basis(cells[id]);

          const ordinal_type numFields = basis.getCardinality();
          const ordinal_type numPoints = 1;
          const ordinal_type spaceDim  = basis.getBaseCellTopology().getDimension();

          
          DynRankView cellCenter("cellCenter", 1, spaceDim), cellVert("cellVert", spaceDim);
          CellTools<DeviceSpaceType>::getReferenceCellCenter(Kokkos::subview(cellCenter, 0, Kokkos::ALL()),
                                                             cellVert, cells[id]);

          // Check VALUE of basis functions: resize vals to rank-2 container:
          {
            DynRankView ConstructWithLabel(vals, numFields, numPoints);
            basis.getValues(vals, cellCenter, OPERATOR_VALUE);
            auto vals_host = Kokkos::create_mirror_view(typename HostSpaceType::memory_space(), vals);
            Kokkos::deep_copy(vals_host, vals);
            for (ordinal_type i=0;i<numFields;++i)
              for (ordinal_type j=0;j<numPoints;++j)
                if (std::abs(vals_host(i,j) - 1.0) > tol) {
                  errorFlag++;
                  *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
                  
                  // Output the multi-index of the value where the error is:
                  *outStream << " At multi-index { ";
                  *outStream << i << " ";*outStream << j << " ";
                  *outStream << "}  computed value: " << vals_host(i,j)
                             << " but reference value: " << 1.0 << "\n";
                }
          }
          
        }
      } catch (std::logic_error err) {
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
      Kokkos::finalize();
      return errorFlag;
    }
  }
}
