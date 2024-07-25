// @HEADER
// *****************************************************************************
//                           Intrepid2 Package
//
// Copyright 2007 NTESS and the Intrepid2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/** \file test_01.cpp
    \brief  Unit tests for the Intrepid2::HGRAD_TRI_Cn_FEM class.
    \author Created by P. Bochev, D. Ridzal, K. Peterson, Kyungjoo Kim
 */


#include "Intrepid2_config.h"

#ifdef HAVE_INTREPID2_DEBUG
#define INTREPID2_TEST_FOR_DEBUG_ABORT_OVERRIDE_TO_CONTINUE
#endif

#include "Intrepid2_Types.hpp"
#include "Intrepid2_Utils.hpp"

#include "Intrepid2_PointTools.hpp"
#include "Intrepid2_HGRAD_TET_Cn_FEM.hpp"

#include "Teuchos_oblackholestream.hpp"
#include "Teuchos_RCP.hpp"

#include "packages/intrepid2/unit-test/Discretization/Basis/Macros.hpp"
#include "packages/intrepid2/unit-test/Discretization/Basis/Setup.hpp"

namespace Intrepid2 {

namespace Test {

template<typename OutValueType, typename PointValueType, typename DeviceType>
int HGRAD_TET_Cn_FEM_Test01(const bool verbose) {

  //! Create an execution space instance.
  const auto space = Kokkos::Experimental::partition_space(typename DeviceType::execution_space {}, 1)[0];

  //! Setup test output stream.
  Teuchos::RCP<std::ostream> outStream = setup_output_stream<DeviceType>(
    verbose, "HGRAD_TET_Cn_FEM", {}
  );

  Teuchos::oblackholestream oldFormatState;
  oldFormatState.copyfmt(std::cout);

  typedef Kokkos::DynRankView<PointValueType,DeviceType> DynRankViewPointValueType;
  typedef Kokkos::DynRankView<OutValueType,DeviceType> DynRankViewOutValueType;
  typedef typename ScalarTraits<OutValueType>::scalar_type scalar_type;
  typedef Kokkos::DynRankView<scalar_type, DeviceType> DynRankViewScalarValueType;

  const scalar_type tol = tolerence();
  int errorFlag = 0;

  // for virtual function, value and point types are declared in the class
  typedef OutValueType outputValueType;
  typedef PointValueType pointValueType;

  typedef Basis_HGRAD_TET_Cn_FEM<DeviceType,outputValueType,pointValueType> TetBasisType;

  constexpr ordinal_type maxOrder = Parameters::MaxOrder ;
  const ordinal_type dim = 3;

  try {

    *outStream
    << "\n"
    << "===============================================================================\n"
    << "| TEST 1: Testing Kronecker property of basis functions                                              |\n"
    << "===============================================================================\n";


    const ordinal_type order = std::min(2,maxOrder);
    TetBasisType tetBasis(order, POINTTYPE_WARPBLEND);

    shards::CellTopology tet_4(shards::getCellTopologyData<shards::Tetrahedron<4> >());
    const ordinal_type np_lattice = PointTools::getLatticeSize(tet_4, order,0);
    const ordinal_type polydim = tetBasis.getCardinality();

    DynRankViewScalarValueType ConstructWithLabel(lattice_scalar, np_lattice , dim);
    tetBasis.getDofCoords(lattice_scalar);

    DynRankViewPointValueType ConstructWithLabelPointView(lattice, np_lattice , dim);

    RealSpaceTools<DeviceType>::clone(lattice,lattice_scalar);
    DynRankViewOutValueType ConstructWithLabelOutView(basisAtLattice, polydim , np_lattice);
    tetBasis.getValues(space, basisAtLattice, lattice, OPERATOR_VALUE);

    auto h_basisAtLattice = Kokkos::create_mirror_view(basisAtLattice);
    Kokkos::deep_copy(h_basisAtLattice, basisAtLattice);

    // Dimensions for the output arrays:
    const ordinal_type numFields = tetBasis.getCardinality();

    // test for Kronecker property
    for (int i=0;i<numFields;i++) {
      for (int j=0;j<np_lattice;j++) {
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


    const ordinal_type order = std::min(4,maxOrder);
    TetBasisType tetBasis(order, POINTTYPE_WARPBLEND);
    auto dofData = tetBasis.getAllDofOrdinal();

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
  try {

    *outStream
    << "\n"
    << "===============================================================================\n"
    << "| TEST 3: Testing OPERATOR_GRAD                                              |\n"
    << "===============================================================================\n";


    const ordinal_type order = std::min(2,maxOrder);
    TetBasisType tetBasis(order, POINTTYPE_EQUISPACED);

    shards::CellTopology tet_4(shards::getCellTopologyData<shards::Tetrahedron<4> >());
    const ordinal_type np_lattice = PointTools::getLatticeSize(tet_4, order,0);
    const ordinal_type polydim = tetBasis.getCardinality();

    //Need to use Scalar type for lattice because PointTools dont's work with FAD types
    DynRankViewScalarValueType ConstructWithLabel(lattice_scalar, np_lattice , dim);
    PointTools::getLattice(lattice_scalar, tet_4, order, 0, POINTTYPE_WARPBLEND);
    DynRankViewPointValueType ConstructWithLabelPointView(lattice, np_lattice , dim);
    RealSpaceTools<DeviceType>::clone(lattice,lattice_scalar);

    DynRankViewOutValueType ConstructWithLabelOutView(dbasisAtLattice, polydim , np_lattice , dim);
    tetBasis.getValues(space, dbasisAtLattice, lattice, OPERATOR_GRAD);

    auto h_dbasisAtLattice = Kokkos::create_mirror_view(dbasisAtLattice);
    Kokkos::deep_copy(h_dbasisAtLattice, dbasisAtLattice);

    std::cout << "GRAD: [";
    for(ordinal_type i=0; i<polydim; i++)
      for(ordinal_type j=0; j<np_lattice; j++)
        for(ordinal_type k=0; k<dim; k++)
          std::cout << h_dbasisAtLattice(i,j,k) << " ";
    std::cout << "]"<< std::endl;


  } catch (std::exception &err) {
    *outStream << "UNEXPECTED ERROR !!! ----------------------------------------------------------\n";
    *outStream << err.what() << '\n';
    *outStream << "-------------------------------------------------------------------------------" << "\n\n";
    errorFlag = -1000;
  };

#if 0 //#ifdef HAVE_INTREPID2_SACADO
  try {

    *outStream
    << "\n"
    << "===============================================================================\n"
    << "| TEST 4: Testing OPERATOR_D2                                                 |\n"
    << "===============================================================================\n";


    const ordinal_type order = 1;
    TetBasisType tetBasis(order, POINTTYPE_EQUISPACED);

    shards::CellTopology tet_4(shards::getCellTopologyData<shards::Tetrahedron<4> >());
    const ordinal_type np_lattice = PointTools::getLatticeSize(tet_4, order,0);
    const ordinal_type polydim = tetBasis.getCardinality();

    DynRankViewScalarValueType ConstructWithLabel(lattice_scalar, np_lattice , dim);
    tetBasis.getDofCoords(lattice_scalar);
    DynRankViewPointValueType ConstructWithLabelPointView(lattice, np_lattice , dim);
    RealSpaceTools<DeviceType>::clone(lattice,lattice_scalar);
    int deriv_order = 2;
    DynRankViewOutValueType ConstructWithLabelOutView(dbasisAtLattice, polydim , np_lattice , (deriv_order+1)*(deriv_order+2)/2);
    tetBasis.getValues(dbasisAtLattice, lattice, OPERATOR_D2);

    auto h_dbasisAtLattice = Kokkos::create_mirror_view(dbasisAtLattice);
    Kokkos::deep_copy(h_dbasisAtLattice, dbasisAtLattice);

    for (int i=0;i<polydim;i++) {
      for (int j=0;j<np_lattice;j++)
        for(ordinal_type k=0; k< (ordinal_type) h_dbasisAtLattice.extent(3); k++){
          if ( std::abs( h_dbasisAtLattice(i,j,k)) > tol ) {
            errorFlag++;
            *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
            *outStream << "Component " << k <<  " of second order derivative of first order polynomial basis function " << i << " does not vanish at node " << j << "\n";
            *outStream << " derivative value is " << h_dbasisAtLattice(i,j,k) << "\n";
          }
        }
    }
  } catch (std::exception &err) {
    *outStream << "UNEXPECTED ERROR !!! ----------------------------------------------------------\n";
    *outStream << err.what() << '\n';
    *outStream << "-------------------------------------------------------------------------------" << "\n\n";
    errorFlag = -1000;
  };
#endif

  *outStream
  << "\n"
  << "===============================================================================\n"
  << "| TEST 5: Function Space is Correct                                           |\n"
  << "===============================================================================\n";

  try {
    const ordinal_type order = std::min(2,maxOrder);
    TetBasisType tetBasis(order, POINTTYPE_WARPBLEND);

    const EFunctionSpace fs = tetBasis.getFunctionSpace();

    if (fs != FUNCTION_SPACE_HGRAD)
    {
      *outStream << std::setw(70) << "------------- TEST FAILURE! -------------" << "\n";

      // Output the multi-index of the value where the error is:
      *outStream << " Expected a function space of FUNCTION_SPACE_HGRAD (enum value " << FUNCTION_SPACE_HGRAD << "),";
      *outStream << " but got " << fs << "\n";
      if (fs == FUNCTION_SPACE_MAX)
      {
        *outStream << "Note that this matches the default value defined by superclass, FUNCTION_SPACE_MAX.  Likely the subclass has failed to set the superclass functionSpace_ field.\n";
      }
      errorFlag++;
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
