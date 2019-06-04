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
    \brief  Unit tests for the Intrepid2::HVOL_TRI_Cn_FEM class.
    \author Created by P. Bochev, D. Ridzal, K. Peterson, Kyungjoo Kim
*/


#include "Intrepid2_config.h"

#ifdef HAVE_INTREPID2_DEBUG
#define INTREPID2_TEST_FOR_DEBUG_ABORT_OVERRIDE_TO_CONTINUE
#endif

#include "Intrepid2_Types.hpp"
#include "Intrepid2_Utils.hpp"

#include "Intrepid2_PointTools.hpp"
#include "Intrepid2_HVOL_TRI_Cn_FEM.hpp"

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


    template<typename OutValueType, typename PointValueType, typename DeviceSpaceType>
    int HVOL_TRI_Cn_FEM_Test01(const bool verbose) {

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
        << "|                 Unit Test (HVOL_TRI_Cn_FEM)                                   |\n"
        << "|                                                                             |\n"
        << "|  Questions? Contact  Pavel Bochev  (pbboche@sandia.gov),                    |\n"
        << "|                      Robert Kirby  (robert.c.kirby@ttu.edu),                |\n"
        << "|                      Denis Ridzal  (dridzal@sandia.gov),                    |\n"
        << "|                      Kara Peterson (kjpeter@sandia.gov),                    |\n"
        << "|                      Kyungjoo Kim  (kyukim@sandia.gov).                     |\n"
        << "|                                                                             |\n"
        << "===============================================================================\n";

      typedef Kokkos::DynRankView<PointValueType,DeviceSpaceType> DynRankViewPointValueType;
      typedef Kokkos::DynRankView<OutValueType,DeviceSpaceType> DynRankViewOutValueType;
      typedef typename ScalarTraits<OutValueType>::scalar_type scalar_type;
      typedef Kokkos::DynRankView<scalar_type, DeviceSpaceType> DynRankViewScalarValueType;

#define ConstructWithLabelScalar(obj, ...) obj(#obj, __VA_ARGS__)

      const scalar_type tol = tolerence();
      int errorFlag = 0;

      typedef Basis_HVOL_TRI_Cn_FEM<DeviceSpaceType,OutValueType,PointValueType> TriBasisType;
      constexpr ordinal_type maxOrder = Parameters::MaxOrder;

      const ordinal_type spaceDim = 2;

      try {

        *outStream
          << "\n"
          << "===============================================================================\n"
          << "| TEST 1: Testing OPERATOR_VALUE                                              |\n"
          << "===============================================================================\n";


        const ordinal_type order = maxOrder;
        TriBasisType triBasis(order, POINTTYPE_WARPBLEND);

        const ordinal_type polydim = triBasis.getCardinality();
        const ordinal_type numPoints = triBasis.getCardinality();
        DynRankViewScalarValueType ConstructWithLabelScalar(lattice_scalar, numPoints, spaceDim);
        DynRankViewPointValueType ConstructWithLabelPointView(lattice, numPoints , spaceDim);

        triBasis.getDofCoords(lattice_scalar);
        RealSpaceTools<DeviceSpaceType>::clone(lattice, lattice_scalar);        

        DynRankViewOutValueType ConstructWithLabelOutView(basisAtLattice, polydim , numPoints);
        triBasis.getValues(basisAtLattice, lattice, OPERATOR_VALUE);

        auto h_basisAtLattice = Kokkos::create_mirror_view(basisAtLattice);
        Kokkos::deep_copy(h_basisAtLattice, basisAtLattice);

            // test for Kronecker property
        for (int i=0;i<polydim;i++) {
          for (int j=0;j<numPoints;j++) {
            if ( i==j && std::abs( h_basisAtLattice(i,j) - 1.0 ) > tol ) {
              errorFlag++;
              *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
              *outStream << " Basis function " << i << " does not have unit value at its node (" << h_basisAtLattice(i,j) <<")\n";
            }
            if ( i!=j && std::abs( h_basisAtLattice(i,j) ) > tol ) {
              errorFlag++;
              *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
              *outStream << " Basis function " << i << " does not vanish at node " << j << "\n";
              *outStream << " Basis function value is " << h_basisAtLattice(i,j) << "\n";
            }
          }
        }
      } catch (std::exception err) {
        *outStream << "UNEXPECTED ERROR !!! ----------------------------------------------------------\n";
        *outStream << err.what() << '\n';
        *outStream << "-------------------------------------------------------------------------------" << "\n\n";
        errorFlag = -1000;
      };

      try {

        *outStream
          << "\n"
          << "===============================================================================\n"
          << "| TEST 2: Testing DOF Data                                              |\n"
          << "===============================================================================\n";


        const ordinal_type order = std::min(3, maxOrder);
        TriBasisType triBasis(order, POINTTYPE_WARPBLEND);
        auto dofData = triBasis.getAllDofOrdinal();
        
      for (unsigned d=0;d<dofData.extent(0);d++) {
      std::cout << "Dimension " << d << "\n";
      for (unsigned f=0;f<dofData.extent(1);f++) {
        int print=-1;
        for (unsigned n=0;n<dofData.extent(2);n++)
          print = std::max(print,dofData(d,f,n));
        if(print == -1) continue;
        std::cout << "\tFacet number " << f << "\n";
        std::cout << "\t\tDOFS:\n";
        for (unsigned n=0;n<dofData.extent(2);n++) {
          if(dofData(d,f,n)>=0)
          std::cout << "\t\t\t" << dofData(d,f,n) << "\n";
        }
      }
    }

      } catch (std::exception err) {
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
