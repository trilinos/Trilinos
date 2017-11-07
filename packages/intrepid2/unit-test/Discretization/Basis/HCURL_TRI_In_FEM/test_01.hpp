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
    \brief  Unit tests for the Intrepid2::HCURL_TRI_In_FEM class.
    \author Created by P. Bochev, D. Ridzal, K. Peterson, Kyungjoo Kim and Mauro Perego
 */


#include "Intrepid2_config.h"

#ifdef HAVE_INTREPID2_DEBUG
#define INTREPID2_TEST_FOR_DEBUG_ABORT_OVERRIDE_TO_CONTINUE
#endif

#include "Intrepid2_Types.hpp"
#include "Intrepid2_Utils.hpp"

#include "Intrepid2_PointTools.hpp"
#include "Intrepid2_HCURL_TRI_In_FEM.hpp"
#include "Intrepid2_HDIV_TRI_In_FEM.hpp"


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
int HCURL_TRI_In_FEM_Test01(const bool verbose) {

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

//  *outStream << "DeviceSpace::  "; DeviceSpaceType::print_configuration(*outStream, false);
//  *outStream << "HostSpace::    ";   HostSpaceType::print_configuration(*outStream, false);

  *outStream
  << "===============================================================================\n"
  << "|                                                                             |\n"
  << "|                      Unit Test HCURL_TRI_In_FEM                             |\n"
  << "|                                                                             |\n"
  << "|  Questions? Contact  Pavel Bochev  (pbboche@sandia.gov),                    |\n"
  << "|                      Robert Kirby  (robert.c.kirby@ttu.edu),                |\n"
  << "|                      Denis Ridzal  (dridzal@sandia.gov),                    |\n"
  << "|                      Kara Peterson (kjpeter@sandia.gov),                    |\n"
  << "|                      Kyungjoo Kim  (kyukim@sandia.gov),                     |\n"
  << "|                      Mauro Perego  (mperego@sandia.gov).                    |\n"
  << "|                                                                             |\n"
  << "===============================================================================\n";

  typedef Kokkos::DynRankView<PointValueType,DeviceSpaceType> DynRankViewPointValueType;
  typedef Kokkos::DynRankView<OutValueType,DeviceSpaceType> DynRankViewOutValueType;
  typedef typename ScalarTraits<OutValueType>::scalar_type scalar_type;
  typedef Kokkos::DynRankView<scalar_type, DeviceSpaceType> DynRankViewScalarValueType;
  typedef Kokkos::DynRankView<scalar_type, HostSpaceType> DynRankViewHostScalarValueType;

#define ConstructWithLabelScalar(obj, ...) obj(#obj, __VA_ARGS__)

  const scalar_type tol = tolerence();
  int errorFlag = 0;
  constexpr ordinal_type dim =2;

  typedef Basis_HCURL_TRI_In_FEM<DeviceSpaceType,OutValueType,PointValueType> TriBasisType;
  typedef Basis_HDIV_TRI_In_FEM<DeviceSpaceType,OutValueType,PointValueType> TriBasisHDivType;
  constexpr ordinal_type maxOrder = Parameters::MaxOrder ;

//  constexpr ordinal_type dim = 2;

  try {

    *outStream
    << "\n"
    << "===============================================================================\n"
    << "| TEST 1: Testing OPERATOR_VALUE (Kronecker property using getDofCoeffs)      |\n"
    << "===============================================================================\n";


    const ordinal_type order = std::min(4, maxOrder);
    TriBasisType triBasis(order, POINTTYPE_WARPBLEND);

    const ordinal_type cardinality = triBasis.getCardinality();
    DynRankViewScalarValueType ConstructWithLabelScalar(dofCoords_scalar, cardinality, dim);
    triBasis.getDofCoords(dofCoords_scalar);

    DynRankViewPointValueType ConstructWithLabelPointView(dofCoords, cardinality , dim);
    RealSpaceTools<DeviceSpaceType>::clone(dofCoords,dofCoords_scalar);

    DynRankViewScalarValueType ConstructWithLabelScalar(dofCoeffs, cardinality , dim);
    triBasis.getDofCoeffs(dofCoeffs);

    DynRankViewOutValueType ConstructWithLabelOutView(basisAtDofCoords, cardinality , cardinality, dim);
    triBasis.getValues(basisAtDofCoords, dofCoords, OPERATOR_VALUE);

    auto h_basisAtDofCoords = Kokkos::create_mirror_view(basisAtDofCoords);
    Kokkos::deep_copy(h_basisAtDofCoords, basisAtDofCoords);

    auto h_dofCoeffs = Kokkos::create_mirror_view(dofCoeffs);
    Kokkos::deep_copy(h_dofCoeffs, dofCoeffs);

    // test for Kronecker property
    for (int i=0;i<cardinality;i++) {
      for (int j=0;j<cardinality;j++) {
        OutValueType dofValue = 0;
        for (ordinal_type k=0;k<dim; k++)
          dofValue += h_basisAtDofCoords(i,j,k)*h_dofCoeffs(j,k);

        if ( i==j && std::abs( dofValue - 1.0 ) > tol ) {
          errorFlag++;
          *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
          *outStream << " Basis function " << i << " does not have unit value at its node (" << dofValue <<")\n";
        }
        if ( i!=j && std::abs( dofValue ) > tol ) {
          errorFlag++;
          *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
          *outStream << " Basis function " << i << " does not vanish at node " << j << "(" << dofValue <<")\n";
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
    << "=======================================================================================\n"
    << "| TEST 2: Testing OPERATOR_VALUE (Kronecker property, reconstructing dofs using tags) |\n"
    << "=======================================================================================\n";


    const ordinal_type order = std::min(4, maxOrder);
    TriBasisType triBasis(order, POINTTYPE_WARPBLEND);

    const ordinal_type cardinality = triBasis.getCardinality();
    DynRankViewScalarValueType ConstructWithLabelScalar(dofCoords_scalar, cardinality , dim);
    triBasis.getDofCoords(dofCoords_scalar);

    DynRankViewPointValueType ConstructWithLabelPointView(dofCoords, cardinality , dim);
    RealSpaceTools<DeviceSpaceType>::clone(dofCoords,dofCoords_scalar);

    DynRankViewOutValueType ConstructWithLabelOutView(basisAtDofCoords, cardinality , cardinality, dim);
    triBasis.getValues(basisAtDofCoords, dofCoords, OPERATOR_VALUE);

    auto h_basisAtDofCoords = Kokkos::create_mirror_view(basisAtDofCoords);
    Kokkos::deep_copy(h_basisAtDofCoords, basisAtDofCoords);

    // Dimensions for the output arrays:
    const ordinal_type numFields = triBasis.getCardinality();

    //Normals at each edge
    DynRankViewHostScalarValueType ConstructWithLabelScalar(tangents, numFields,dim); // normals at each point basis point
    tangents(0,0)  =  1.0; tangents(0,1)  = 0.0;
    tangents(1,0)  =  -1.0; tangents(1,1)  =  1.0;
    tangents(2,0)  =  0.0; tangents(2,1)  =  -1.0;

    const auto allTags = triBasis.getAllDofTags();

    // test for Kronecker property
    for (int i=0;i<numFields;i++) {
      for (int j=0;j<numFields;j++) {
        OutValueType dofValue = 0;
        if(allTags(j,0) == dim-1) { //edge
          auto edgeId = allTags(j,1);
          for (ordinal_type k=0;k<dim; k++)
            dofValue += h_basisAtDofCoords(i,j,k)*tangents(edgeId,k);
        }
        else { //elem
          auto elemDofId =  allTags(j,2);
          dofValue = h_basisAtDofCoords(i,j,elemDofId%dim);
        }

        if ( i==j && std::abs( dofValue - 1.0 ) > tol ) {
          errorFlag++;
          *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
          *outStream << " Basis function " << i << " does not have unit value at its node (" << dofValue <<")\n";
        }
        if ( i!=j && std::abs( dofValue ) > tol ) {
          errorFlag++;
          *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
          *outStream << " Basis function " << i << " does not vanish at node " << j << "(" << dofValue <<")\n";
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
    << "| TEST 3: Testing OPERATOR_VALUE (orthogonality with HDiv)                    |\n"
    << "===============================================================================\n";


    const ordinal_type order = 1;

    TriBasisType triBasis(order, POINTTYPE_EQUISPACED);
    TriBasisHDivType triBasisHDiv(order, POINTTYPE_EQUISPACED);

    shards::CellTopology tri_3(shards::getCellTopologyData<shards::Triangle<3> >());
    const ordinal_type np_lattice = PointTools::getLatticeSize(tri_3, order,0);
    const ordinal_type cardinality = triBasis.getCardinality();
    DynRankViewHostScalarValueType ConstructWithLabelScalar(lattice_host_scalar, np_lattice , dim);
    PointTools::getLattice(lattice_host_scalar, tri_3, order, 0, POINTTYPE_EQUISPACED);

    auto lattice_scalar = Kokkos::create_mirror_view(typename DeviceSpaceType::memory_space(), lattice_host_scalar);
    deep_copy(lattice_scalar, lattice_host_scalar);

    DynRankViewPointValueType ConstructWithLabelPointView(lattice, np_lattice , dim);
    RealSpaceTools<DeviceSpaceType>::clone(lattice,lattice_scalar);

    DynRankViewOutValueType ConstructWithLabelOutView(basisAtLattice, cardinality , np_lattice, dim);
    triBasis.getValues(basisAtLattice, lattice, OPERATOR_VALUE);
    DynRankViewOutValueType ConstructWithLabelOutView(basisHDivAtLattice, cardinality , np_lattice, dim);
    triBasisHDiv.getValues(basisHDivAtLattice, lattice, OPERATOR_VALUE);

    auto h_basisAtLattice = Kokkos::create_mirror_view(basisAtLattice);
    Kokkos::deep_copy(h_basisAtLattice, basisAtLattice);
    auto h_basisHDivAtLattice = Kokkos::create_mirror_view(basisHDivAtLattice);
    Kokkos::deep_copy(h_basisHDivAtLattice, basisHDivAtLattice);
    const ordinal_type numFields = triBasis.getCardinality();

    for (int i=0;i<numFields;i++) {
      for (int j=0;j<np_lattice;j++) {
        if (std::abs( h_basisHDivAtLattice(i,j,1) + h_basisAtLattice(i,j,0) ) > tol ) {
          errorFlag++;
          *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
          // Output the multi-index of the value where the error is:
          *outStream << " At multi-index { ";
          *outStream << i << " " << j << " and component 0";
          *outStream << "}  computed value: " << h_basisAtLattice(i,j,0)
                            << " but correct value: " << -h_basisHDivAtLattice(i,j,1) << "\n";
          *outStream << "Difference: " << std::abs( h_basisAtLattice(i,j,0) + h_basisHDivAtLattice(i,j,1) ) << "\n";
        }
        if (std::abs( h_basisHDivAtLattice(i,j,0) - h_basisAtLattice(i,j,1) ) > tol ) {
          errorFlag++;
          *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
          // Output the multi-index of the value where the error is:
          *outStream << " At multi-index { ";
          *outStream << i << " " << j << " and component 1";
          *outStream << "}  computed value: " << h_basisAtLattice(i,j,1)
                            << " but correct value: " << h_basisHDivAtLattice(i,j,0) << "\n";
          *outStream << "Difference: " << std::abs( h_basisAtLattice(i,j,1) - h_basisHDivAtLattice(i,j,1) ) << "\n";
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
    << "| TEST 4: Testing OPERATOR_VALUE (FIAT Values)                                |\n"
    << "===============================================================================\n";


    constexpr ordinal_type order = 2;
    if(order <= maxOrder) {
      TriBasisType triBasis(order, POINTTYPE_EQUISPACED);

      shards::CellTopology tri_3(shards::getCellTopologyData<shards::Triangle<3> >());
      const ordinal_type np_lattice = PointTools::getLatticeSize(tri_3, order,0);
      const ordinal_type cardinality = triBasis.getCardinality();
      DynRankViewHostScalarValueType ConstructWithLabelScalar(lattice_host_scalar, np_lattice , dim);
      PointTools::getLattice(lattice_host_scalar, tri_3, order, 0, POINTTYPE_EQUISPACED);

      auto lattice_scalar = Kokkos::create_mirror_view(typename DeviceSpaceType::memory_space(), lattice_host_scalar);
      deep_copy(lattice_scalar, lattice_host_scalar);

      DynRankViewPointValueType ConstructWithLabelPointView(lattice, np_lattice , dim);
      RealSpaceTools<DeviceSpaceType>::clone(lattice,lattice_scalar);

      DynRankViewOutValueType ConstructWithLabelOutView(basisAtLattice, cardinality , np_lattice, dim);
      triBasis.getValues(basisAtLattice, lattice, OPERATOR_VALUE);

      auto h_basisAtLattice = Kokkos::create_mirror_view(basisAtLattice);
      Kokkos::deep_copy(h_basisAtLattice, basisAtLattice);

      const scalar_type fiat_vals[] = {
          2.000000000000000e+00, 0.000000000000000e+00,
          5.000000000000000e-01, 2.500000000000000e-01,
          -1.000000000000000e+00, -1.000000000000000e+00,
          2.500000000000000e-01, 0.000000000000000e+00,
          -5.000000000000000e-01, -5.000000000000000e-01,
          0.000000000000000e+00, 0.000000000000000e+00,
          -1.000000000000000e+00, 0.000000000000000e+00,
          5.000000000000000e-01, 2.500000000000000e-01,
          2.000000000000000e+00, 2.000000000000000e+00,
          -5.000000000000000e-01, 0.000000000000000e+00,
          2.500000000000000e-01, 2.500000000000000e-01,
          0.000000000000000e+00, 0.000000000000000e+00,
          0.000000000000000e+00, 0.000000000000000e+00,
          0.000000000000000e+00, 2.500000000000000e-01,
          0.000000000000000e+00, 2.000000000000000e+00,
          5.000000000000000e-01, 0.000000000000000e+00,
          -2.500000000000000e-01, 2.500000000000000e-01,
          1.000000000000000e+00, 0.000000000000000e+00,
          0.000000000000000e+00, 0.000000000000000e+00,
          0.000000000000000e+00, -5.000000000000000e-01,
          0.000000000000000e+00, -1.000000000000000e+00,
          -2.500000000000000e-01, 0.000000000000000e+00,
          -2.500000000000000e-01, 2.500000000000000e-01,
          -2.000000000000000e+00, 0.000000000000000e+00,
          0.000000000000000e+00, 1.000000000000000e+00,
          0.000000000000000e+00, 5.000000000000000e-01,
          0.000000000000000e+00, 0.000000000000000e+00,
          -2.500000000000000e-01, -5.000000000000000e-01,
          -2.500000000000000e-01, -2.500000000000000e-01,
          -2.000000000000000e+00, -2.000000000000000e+00,
          0.000000000000000e+00, -2.000000000000000e+00,
          0.000000000000000e+00, -2.500000000000000e-01,
          0.000000000000000e+00, 0.000000000000000e+00,
          -2.500000000000000e-01, -5.000000000000000e-01,
          5.000000000000000e-01, 5.000000000000000e-01,
          1.000000000000000e+00, 1.000000000000000e+00,
          0.000000000000000e+00, 0.000000000000000e+00,
          0.000000000000000e+00, -7.500000000000000e-01,
          0.000000000000000e+00, 0.000000000000000e+00,
          1.500000000000000e+00, 0.000000000000000e+00,
          7.500000000000000e-01, 7.500000000000000e-01,
          0.000000000000000e+00, 0.000000000000000e+00,
          0.000000000000000e+00, 0.000000000000000e+00,
          0.000000000000000e+00, 1.500000000000000e+00,
          0.000000000000000e+00, 0.000000000000000e+00,
          -7.500000000000000e-01, 0.000000000000000e+00,
          7.500000000000000e-01, 7.500000000000000e-01,
          0.000000000000000e+00, 0.000000000000000e+00
      };

      ordinal_type cur=0;
      for (ordinal_type i=0;i<cardinality;i++) {
        for (ordinal_type j=0;j<np_lattice;j++) {
          for (ordinal_type k=0;k<dim; k++) {
            if (std::abs( h_basisAtLattice(i,j,k) - fiat_vals[cur] ) > tol ) {
              errorFlag++;
              *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";

              // Output the multi-index of the value where the error is:
              *outStream << " At multi-index { ";
              *outStream << i << " " << j << " " << k;
              *outStream << "}  computed value: " <<  h_basisAtLattice(i,j,k)
                            << " but correct value: " << fiat_vals[cur] << "\n";
              *outStream << "Difference: " << std::fabs(  h_basisAtLattice(i,j,k) - fiat_vals[cur] ) << "\n";
            }
            cur++;
          }
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
    << "| TEST 5: Testing Operator CURL                                               |\n"
    << "===============================================================================\n";


    constexpr ordinal_type order = 2;
    if(order <= maxOrder) {
      shards::CellTopology tri_3(shards::getCellTopologyData<shards::Triangle<3> >());
      TriBasisType triBasis(order, POINTTYPE_EQUISPACED);
      const ordinal_type cardinality = triBasis.getCardinality();
      const ordinal_type np_lattice = PointTools::getLatticeSize(tri_3, order,0);
      DynRankViewHostScalarValueType ConstructWithLabelScalar(lattice_host_scalar, np_lattice , dim);
      PointTools::getLattice(lattice_host_scalar, tri_3, order, 0, POINTTYPE_EQUISPACED);

      auto lattice_scalar = Kokkos::create_mirror_view(typename DeviceSpaceType::memory_space(), lattice_host_scalar);
      deep_copy(lattice_scalar, lattice_host_scalar);

      DynRankViewPointValueType ConstructWithLabelPointView(lattice, np_lattice , dim);
      RealSpaceTools<DeviceSpaceType>::clone(lattice,lattice_scalar);

      DynRankViewOutValueType ConstructWithLabelOutView(basisDivAtLattice, cardinality , np_lattice);
      triBasis.getValues(basisDivAtLattice, lattice, OPERATOR_CURL);


      auto h_basisDivAtLattice = Kokkos::create_mirror_view(basisDivAtLattice);
      Kokkos::deep_copy(h_basisDivAtLattice, basisDivAtLattice);

      const scalar_type fiat_divs[] = {
        7.000000000000000e+00,
        2.500000000000000e+00,
        -2.000000000000000e+00,
        2.500000000000000e+00,
        -2.000000000000000e+00,
        -2.000000000000000e+00,
        -2.000000000000000e+00,
        2.500000000000000e+00,
        7.000000000000000e+00,
        -2.000000000000000e+00,
        2.500000000000000e+00,
        -2.000000000000000e+00,
        -2.000000000000000e+00,
        2.500000000000000e+00,
        7.000000000000000e+00,
        -2.000000000000000e+00,
        2.500000000000000e+00,
        -2.000000000000000e+00,
        -2.000000000000000e+00,
        -2.000000000000000e+00,
        -2.000000000000000e+00,
        2.500000000000000e+00,
        2.500000000000000e+00,
        7.000000000000000e+00,
        -2.000000000000000e+00,
        -2.000000000000000e+00,
        -2.000000000000000e+00,
        2.500000000000000e+00,
        2.500000000000000e+00,
        7.000000000000000e+00,
        7.000000000000000e+00,
        2.500000000000000e+00,
        -2.000000000000000e+00,
        2.500000000000000e+00,
        -2.000000000000000e+00,
        -2.000000000000000e+00,
        -9.000000000000000e+00,
        -4.500000000000000e+00,
        0.000000000000000e+00,
        0.000000000000000e+00,
        4.500000000000000e+00,
        9.000000000000000e+00,
        9.000000000000000e+00,
        0.000000000000000e+00,
        -9.000000000000000e+00,
        4.500000000000000e+00,
        -4.500000000000000e+00,
        0.000000000000000e+00
    };


    ordinal_type cur=0;
    for (ordinal_type i=0;i<cardinality;i++) {
      for (ordinal_type j=0;j<np_lattice;j++) {
        if (std::abs( h_basisDivAtLattice(i,j) - fiat_divs[cur] ) > tol ) {
          errorFlag++;
          *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";

          // Output the multi-index of the value where the error is:
          *outStream << " At multi-index { ";
          *outStream << i << " " << j;
          *outStream << "}  computed value: " <<  h_basisDivAtLattice(i,j)
                          << " but correct value: " << fiat_divs[cur] << "\n";
          *outStream << "Difference: " << std::fabs(  h_basisDivAtLattice(i,j) - fiat_divs[cur] ) << "\n";
        }
        cur++;
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
