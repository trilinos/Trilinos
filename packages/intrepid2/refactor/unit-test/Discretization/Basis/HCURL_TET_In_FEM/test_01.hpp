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
    \brief  Unit tests for the Intrepid2::HCURL_TET_In_FEM class.
    \author Created by P. Bochev, D. Ridzal, K. Peterson, Kyungjoo Kim and Mauro Perego
 */


#include "Intrepid2_config.h"

#ifdef HAVE_INTREPID2_DEBUG
#define INTREPID2_TEST_FOR_DEBUG_ABORT_OVERRIDE_TO_CONTINUE
#endif

#include "Intrepid2_Types.hpp"
#include "Intrepid2_Utils.hpp"

#include "Intrepid2_PointTools.hpp"
#include "Intrepid2_HCURL_TET_In_FEM.hpp"


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
int HCURL_TET_In_FEM_Test01(const bool verbose) {

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
  << "|                      Unit Test HCURL_TET_In_FEM                             |\n"
  << "|                                                                             |\n"
  << "|  Questions? Contact  Pavel Bochev  (pbboche@sandia.gov),                    |\n"
  << "|                      Robert Kirby  (robert.c.kirby@ttu.edu),                |\n"
  << "|                      Denis Ridzal  (dridzal@sandia.gov),                    |\n"
  << "|                      Kara Peterson (kjpeter@sandia.gov),                    |\n"
  << "|                      Kyungjoo Kim  (kyukim@sandia.gov),                     |\n"
  << "|                      Mauro Perego  (mperego@sandia.gov).                    |\n"
  << "|                                                                             |\n"
  << "===============================================================================\n";

  typedef Kokkos::DynRankView<ValueType,DeviceSpaceType> DynRankView;
  typedef Kokkos::DynRankView<ValueType,HostSpaceType>   DynRankViewHost;
#define ConstructWithLabel(obj, ...) obj(#obj, __VA_ARGS__)

  const ValueType tol = tolerence();
  int errorFlag = 0;

  // for virtual function, value and point types are declared in the class
  typedef ValueType outputValueType;
  typedef ValueType pointValueType;

  typedef Basis_HCURL_TET_In_FEM<DeviceSpaceType,outputValueType,pointValueType> TetBasisType;

  constexpr ordinal_type maxOrder = Parameters::MaxOrder;
  constexpr ordinal_type dim = 3;

  try {

    *outStream
    << "\n"
    << "===============================================================================\n"
    << "| TEST 1: Testing OPERATOR_VALUE (Kronecker property using getDofCoeffs)      |\n"
    << "===============================================================================\n";


    const ordinal_type order = std::min(3, maxOrder);
    TetBasisType tetBasis(order, POINTTYPE_WARPBLEND);

    const ordinal_type cardinality = tetBasis.getCardinality();
    DynRankView ConstructWithLabel(dofCoords, cardinality , dim);
    tetBasis.getDofCoords(dofCoords);

    DynRankView ConstructWithLabel(dofCoeffs, cardinality , dim);
    tetBasis.getDofCoeffs(dofCoeffs);

    DynRankView ConstructWithLabel(basisAtDofCoords, cardinality , cardinality, dim);
    tetBasis.getValues(basisAtDofCoords, dofCoords, OPERATOR_VALUE);

    auto h_basisAtDofCoords = Kokkos::create_mirror_view(basisAtDofCoords);
    Kokkos::deep_copy(h_basisAtDofCoords, basisAtDofCoords);

    auto h_dofCoeffs = Kokkos::create_mirror_view(dofCoeffs);
    Kokkos::deep_copy(h_dofCoeffs, dofCoeffs);

    // test for Kronecker property
    for (int i=0;i<cardinality;i++) {
      for (int j=0;j<cardinality;j++) {
        outputValueType dofValue = 0;
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


    const ordinal_type order = std::min(3, maxOrder);
    TetBasisType tetBasis(order, POINTTYPE_WARPBLEND);

    shards::CellTopology tet_4(shards::getCellTopologyData<shards::Tetrahedron<4> >());
    const ordinal_type polydim = tetBasis.getCardinality();
    DynRankView ConstructWithLabel(dofCoords, polydim , dim);
    tetBasis.getDofCoords(dofCoords);

    DynRankView ConstructWithLabel(basisAtDofCoords, polydim , polydim, dim);
    tetBasis.getValues(basisAtDofCoords, dofCoords, OPERATOR_VALUE);

    auto h_basisAtDofCoords = Kokkos::create_mirror_view(basisAtDofCoords);
    Kokkos::deep_copy(h_basisAtDofCoords, basisAtDofCoords);

    // Dimensions for the output arrays:
    const ordinal_type numFields = tetBasis.getCardinality();

    //Normals at each edge
    DynRankViewHost ConstructWithLabel(tangents, numFields,dim); // normals at each point basis point


    DynRankViewHost ConstructWithLabel(edgeTan, dim );
    DynRankViewHost ConstructWithLabel(faceTan1, dim );
    DynRankViewHost ConstructWithLabel(faceTan2, dim );

    const auto allTags = tetBasis.getAllDofTags();

    // test for Kronecker property
    for (int i=0;i<numFields;i++) {
      for (int j=0;j<numFields;j++) {
        outputValueType dofValue = 0;
        if(allTags(j,0) == 1) { //edge
          auto edgeId = allTags(j,1);
          CellTools<Kokkos::HostSpace::execution_space>::getReferenceEdgeTangent( edgeTan ,
              edgeId ,
              tet_4 );



          for (ordinal_type jj=0;jj<dim;jj++)
                edgeTan(jj) *= 2.0;
          for (ordinal_type k=0;k<dim; k++)
            dofValue += h_basisAtDofCoords(i,j,k)*edgeTan(k);
        }
        else if(allTags(j,0) == 2) { //face
          auto faceId = allTags(j,1);
          CellTools<Kokkos::HostSpace::execution_space>::getReferenceFaceTangents( faceTan1 ,
                    faceTan2,
                    faceId ,
                    tet_4  );

          auto faceDofId =  allTags(j,2);
          dofValue = 0;
          if(faceDofId%2 ==0 )
            for (ordinal_type k=0;k<dim; k++)
              dofValue += h_basisAtDofCoords(i,j,k)*faceTan1(k);
          else
            for (ordinal_type k=0;k<dim; k++)
              dofValue += h_basisAtDofCoords(i,j,k)*faceTan2(k);
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
    << "| TEST 3: Testing OPERATOR_VALUE (FIAT Values)                                |\n"
    << "===============================================================================\n";


    const ordinal_type order = 1;
    TetBasisType tetBasis(order, POINTTYPE_EQUISPACED);

    shards::CellTopology tet_4(shards::getCellTopologyData<shards::Tetrahedron<4> >());
    const ordinal_type np_lattice = PointTools::getLatticeSize(tet_4, order,0);
    const ordinal_type polydim = tetBasis.getCardinality();
    DynRankView ConstructWithLabel(lattice, np_lattice , dim);
    PointTools::getLattice(lattice, tet_4, order, 0, POINTTYPE_EQUISPACED);

    DynRankView ConstructWithLabel(basisAtLattice, polydim , np_lattice, dim);
    tetBasis.getValues(basisAtLattice, lattice, OPERATOR_VALUE);

    auto h_basisAtLattice = Kokkos::create_mirror_view(basisAtLattice);
    Kokkos::deep_copy(h_basisAtLattice, basisAtLattice);

    // Dimensions for the output arrays:
    const ordinal_type numFields = tetBasis.getCardinality();

    const double fiat_vals[] = {
1.000000000000001e+00, -2.498001805406602e-16, -1.665334536937735e-16,
9.999999999999998e-01, 1.000000000000000e+00, 1.000000000000000e+00,
5.828670879282072e-16, 1.110223024625157e-16, 2.498001805406602e-16,
7.771561172376096e-16, 8.326672684688674e-17, 1.110223024625157e-16,
2.081668171172169e-16, -2.914335439641036e-16, 1.280865063236792e-16,
-3.191891195797325e-16, 1.000000000000000e+00, -4.293998586504916e-17,
-9.999999999999994e-01, 2.081668171172169e-16, 2.400576428367544e-16,
2.220446049250313e-16, -5.551115123125783e-17, 1.084013877651281e-16,
3.469446951953614e-16, -1.000000000000000e+00, 1.387778780781446e-16,
-1.804112415015879e-16, 1.942890293094024e-16, -1.387778780781446e-16,
-9.999999999999993e-01, -9.999999999999996e-01, -9.999999999999998e-01,
5.551115123125783e-17, -2.220446049250313e-16, -8.326672684688674e-17,
-2.220446049250313e-16, -5.551115123125783e-17, 9.999999999999999e-01,
1.665334536937735e-16, 1.110223024625157e-16, -6.383782391594650e-16,
1.110223024625157e-16, 1.110223024625157e-16, -1.110223024625157e-16,
9.999999999999990e-01, 9.999999999999994e-01, 9.999999999999996e-01,
1.387778780781446e-16, -2.496931404305374e-17, -1.665334536937735e-16,
-2.498001805406602e-16, -2.149987498083074e-16, 1.000000000000000e+00,
8.326672684688674e-17, -3.769887250591415e-17, 8.326672684688674e-17,
-9.999999999999994e-01, 1.556977698723022e-16, 2.220446049250313e-16,
-9.422703950001342e-18, 1.665334536937735e-16, -2.359223927328458e-16,
-9.422703950001268e-18, -8.326672684688674e-17, 1.387778780781446e-17,
-7.525083148581445e-17, 2.775557561562891e-17, 1.000000000000000e+00,
2.789513560273035e-16, -9.999999999999998e-01, -5.551115123125783e-17
    };

    ordinal_type cur=0;
    for (ordinal_type i=0;i<numFields;i++) {
      for (ordinal_type j=0;j<np_lattice;j++) {
        for (ordinal_type k=0;k<dim; k++) {
          if (std::fabs( h_basisAtLattice(i,j,k) - fiat_vals[cur] ) > tol ) {
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
    << "| TEST 4: Testing Operator CURL                                               |\n"
    << "===============================================================================\n";


    const ordinal_type order = 1;
    TetBasisType tetBasis(order, POINTTYPE_EQUISPACED);

    shards::CellTopology tet_4(shards::getCellTopologyData<shards::Tetrahedron<4> >());
    const ordinal_type np_lattice = PointTools::getLatticeSize(tet_4, order,0);
    const ordinal_type polydim = tetBasis.getCardinality();
    DynRankView ConstructWithLabel(lattice, np_lattice , dim);
    PointTools::getLattice(lattice, tet_4, order, 0, POINTTYPE_EQUISPACED);

    DynRankView ConstructWithLabel(curlBasisAtLattice, polydim , np_lattice, dim);
    tetBasis.getValues(curlBasisAtLattice, lattice, OPERATOR_CURL);

    auto h_curlBasisAtLattice = Kokkos::create_mirror_view(curlBasisAtLattice);
    Kokkos::deep_copy(h_curlBasisAtLattice, curlBasisAtLattice);

    // Dimensions for the output arrays:
    const ordinal_type numFields = tetBasis.getCardinality();

    const double fiat_curls[] = {
      -5.551115123125783e-16, -2.000000000000000e+00, 2.000000000000000e+00,
      -5.551115123125783e-16, -2.000000000000000e+00, 2.000000000000000e+00,
      -5.551115123125783e-16, -2.000000000000000e+00, 2.000000000000000e+00,
      -5.551115123125783e-16, -2.000000000000000e+00, 2.000000000000000e+00,
      -4.440892098500626e-16, -2.775557561562891e-16, 2.000000000000000e+00,
      -4.440892098500626e-16, -2.775557561562891e-16, 2.000000000000000e+00,
      -4.440892098500626e-16, -2.775557561562891e-16, 2.000000000000000e+00,
      -4.440892098500626e-16, -2.775557561562891e-16, 2.000000000000000e+00,
      -2.000000000000000e+00, 5.551115123125783e-17, 2.000000000000000e+00,
      -2.000000000000000e+00, 5.551115123125783e-17, 2.000000000000000e+00,
      -2.000000000000000e+00, 5.551115123125783e-17, 2.000000000000000e+00,
      -2.000000000000000e+00, 5.551115123125783e-17, 2.000000000000000e+00,
      -2.000000000000000e+00, 2.000000000000000e+00, 9.861075762086680e-17,
      -2.000000000000000e+00, 2.000000000000000e+00, 9.861075762086680e-17,
      -2.000000000000000e+00, 2.000000000000000e+00, 9.861075762086680e-17,
      -2.000000000000000e+00, 2.000000000000000e+00, 9.861075762086680e-17,
      -2.775557561562891e-17, -2.000000000000000e+00, 4.287451790760826e-16,
      -2.775557561562891e-17, -2.000000000000000e+00, 4.287451790760826e-16,
      -2.775557561562891e-17, -2.000000000000000e+00, 4.287451790760826e-16,
      -2.775557561562891e-17, -2.000000000000000e+00, 4.287451790760826e-16,
      2.000000000000000e+00, -2.185751579730777e-16, 1.526556658859590e-16,
      2.000000000000000e+00, -2.185751579730777e-16, 1.526556658859590e-16,
      2.000000000000000e+00, -2.185751579730777e-16, 1.526556658859590e-16,
      2.000000000000000e+00, -2.185751579730777e-16, 1.526556658859590e-16
    };


    ordinal_type cur=0;
    for (ordinal_type i=0;i<numFields;i++) {
      for (ordinal_type j=0;j<np_lattice;j++) {
        for (ordinal_type k=0;k<dim; k++) {
          if (std::fabs( h_curlBasisAtLattice(i,j,k) - fiat_curls[cur] ) > tol ) {
            errorFlag++;
            *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";

            // Output the multi-index of the value where the error is:
            *outStream << " At multi-index { ";
            *outStream << i << " " << j;
            *outStream << "}  computed value: " <<  h_curlBasisAtLattice(i,j,k)
                            << " but correct value: " << fiat_curls[cur] << "\n";
            *outStream << "Difference: " << std::fabs(  h_curlBasisAtLattice(i,j,k) - fiat_curls[cur] ) << "\n";
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
