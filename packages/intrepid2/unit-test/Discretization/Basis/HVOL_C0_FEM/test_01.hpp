// @HEADER
// *****************************************************************************
//                           Intrepid2 Package
//
// Copyright 2007 NTESS and the Intrepid2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Intrepid2_config.h"

#ifdef HAVE_INTREPID2_DEBUG
#define INTREPID2_TEST_FOR_DEBUG_ABORT_OVERRIDE_TO_CONTINUE
#endif

#include "Intrepid2_Basis.hpp"
#include "Intrepid2_HVOL_C0_FEM.hpp"

#include "Teuchos_oblackholestream.hpp"
#include "Teuchos_RCP.hpp"

#include "packages/intrepid2/unit-test/Discretization/Basis/Macros.hpp"

namespace Intrepid2 {

  namespace Test {

    template<typename ValueType, typename DeviceType>
    int BasisConst_Test01(const bool verbose) {

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
        << "===============================================================================\n"
        << "|                                                                             |\n"
        << "|                 Unit Test (Basis_HVOL_C0_FEM)                              |\n"
        << "|                                                                             |\n"
        << "|  Intrepid's website: http://trilinos.sandia.gov/packages/intrepid           |\n"
        << "|  Trilinos website:   http://trilinos.sandia.gov                             |\n"
        << "|                                                                             |\n"
        << "===============================================================================\n";

      typedef Kokkos::DynRankView<ValueType,DeviceType> DynRankView;

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
          Basis_HVOL_C0_FEM<DeviceType,outputValueType,pointValueType> basis(cells[id]);
        }

      } catch (std::logic_error &err){
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
          Basis_HVOL_C0_FEM<DeviceType,outputValueType,pointValueType> basis(cells[id]);

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
      } catch (std::logic_error &err){
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
          Basis_HVOL_C0_FEM<DeviceType,outputValueType,pointValueType> basis(cells[id]);

          const ordinal_type numFields = basis.getCardinality();
          const ordinal_type numPoints = 1;
          const ordinal_type spaceDim  = basis.getBaseCellTopology().getDimension();

          
          DynRankView cellCenter("cellCenter", 1, spaceDim);
          auto cellCenter_host = Kokkos::create_mirror_view(cellCenter);
          CellTools<typename HostSpaceType::device_type>
            ::getReferenceCellCenter(Kokkos::subview(cellCenter_host, 0, Kokkos::ALL()),
                                     cells[id]);
          Kokkos::deep_copy(cellCenter, cellCenter_host);

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
      } catch (std::logic_error &err) {
        *outStream << "UNEXPECTED ERROR !!! ----------------------------------------------------------\n";
        *outStream << err.what() << '\n';
        *outStream << "-------------------------------------------------------------------------------" << "\n\n";
        errorFlag = -1000;
      };

      *outStream
      << "\n"
      << "===============================================================================\n"
      << "| TEST 4: Function Space is Correct                                           |\n"
      << "===============================================================================\n";
      
      try {
        for (ordinal_type id=0;id<6;++id) {
          Basis_HVOL_C0_FEM<DeviceType,outputValueType,pointValueType> basis(cells[id]);
          
          const EFunctionSpace fs      = basis.getFunctionSpace();
          
          if (fs != FUNCTION_SPACE_HVOL)
          {
            *outStream << std::setw(70) << "------------- TEST FAILURE! -------------" << "\n";
            
            // Output the multi-index of the value where the error is:
            *outStream << " Expected a function space of FUNCTION_SPACE_HVOL (enum value " << FUNCTION_SPACE_HVOL << "),";
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
