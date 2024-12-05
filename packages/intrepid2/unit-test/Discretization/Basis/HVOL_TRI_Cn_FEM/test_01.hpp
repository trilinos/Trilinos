// @HEADER
// *****************************************************************************
//                           Intrepid2 Package
//
// Copyright 2007 NTESS and the Intrepid2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
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

#include "packages/intrepid2/unit-test/Discretization/Basis/Macros.hpp"

namespace Intrepid2 {

  namespace Test {

    template<typename OutValueType, typename PointValueType, typename DeviceType>
    int HVOL_TRI_Cn_FEM_Test01(const bool verbose) {

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
        << "|                 Unit Test (HVOL_TRI_Cn_FEM)                                   |\n"
        << "|                                                                             |\n"
        << "|  Questions? Contact  Pavel Bochev  (pbboche@sandia.gov),                    |\n"
        << "|                      Robert Kirby  (robert.c.kirby@ttu.edu),                |\n"
        << "|                      Denis Ridzal  (dridzal@sandia.gov),                    |\n"
        << "|                      Kara Peterson (kjpeter@sandia.gov),                    |\n"
        << "|                      Kyungjoo Kim  (kyukim@sandia.gov).                     |\n"
        << "|                                                                             |\n"
        << "===============================================================================\n";

      typedef Kokkos::DynRankView<PointValueType,DeviceType> DynRankViewPointValueType;
      typedef Kokkos::DynRankView<OutValueType,DeviceType> DynRankViewOutValueType;
      typedef typename ScalarTraits<OutValueType>::scalar_type scalar_type;
      typedef Kokkos::DynRankView<scalar_type, DeviceType> DynRankViewScalarValueType;

      const scalar_type tol = tolerence();
      int errorFlag = 0;

      typedef Basis_HVOL_TRI_Cn_FEM<DeviceType,OutValueType,PointValueType> TriBasisType;
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
        DynRankViewScalarValueType ConstructWithLabel(lattice_scalar, numPoints, spaceDim);
        DynRankViewPointValueType ConstructWithLabelPointView(lattice, numPoints , spaceDim);

        triBasis.getDofCoords(lattice_scalar);
        RealSpaceTools<DeviceType>::clone(lattice, lattice_scalar);        

        DynRankViewOutValueType ConstructWithLabelOutView(basisAtLattice, polydim , numPoints);
        triBasis.getValues(basisAtLattice, lattice, OPERATOR_VALUE);

        auto h_basisAtLattice = Kokkos::create_mirror_view(basisAtLattice);
        Kokkos::deep_copy(h_basisAtLattice, basisAtLattice);

            // test for Kronecker property
        for (int i=0;i<polydim;i++) {
          for (int j=0;j<numPoints;j++) {
            if ( i==j && std::abs( h_basisAtLattice(i,j) - 1.0 ) > tol * 10 ) { // relax tolerance now that we support orders up to 20
              errorFlag++;
              *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
              *outStream << " Basis function " << i << " does not have unit value at its node (" << h_basisAtLattice(i,j) <<")\n";
            }
            if ( i!=j && std::abs( h_basisAtLattice(i,j) ) > tol * 10 ) { // relax tolerance now that we support orders up to 20
              errorFlag++;
              *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
              *outStream << " Basis function " << i << " does not vanish at node " << j << "\n";
              *outStream << " Basis function value is " << h_basisAtLattice(i,j) << "\n";
            }
          }
        }
      } catch (std::exception &err) {
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

      } catch (std::exception &err) {
        *outStream << "UNEXPECTED ERROR !!! ----------------------------------------------------------\n";
        *outStream << err.what() << '\n';
        *outStream << "-------------------------------------------------------------------------------" << "\n\n";
        errorFlag = -1000;
      };


      *outStream
      << "\n"
      << "===============================================================================\n"
      << "| TEST 3: Function Space is Correct                                           |\n"
      << "===============================================================================\n";

      try {
        for (auto order=0;order<std::min(3, maxOrder);++order) {
          TriBasisType triBasis(order, POINTTYPE_WARPBLEND);

          const EFunctionSpace fs = triBasis.getFunctionSpace();

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
