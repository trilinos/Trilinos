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

/** \file test_01.hpp
    \brief  Unit tests for the Intrepid2::HDIV_HEX_In_FEM class.
    \author Created by P. Bochev, D. Ridzal, K. Peterson, Kyungjoo Kim and Mauro Perego
 */


#include "Intrepid2_config.h"

#ifdef HAVE_INTREPID2_DEBUG
#define INTREPID2_TEST_FOR_DEBUG_ABORT_OVERRIDE_TO_CONTINUE
#endif

#include "Intrepid2_Types.hpp"
#include "Intrepid2_Utils.hpp"

#include "Intrepid2_PointTools.hpp"
#include "Intrepid2_HDIV_HEX_In_FEM.hpp"


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
int HDIV_HEX_In_FEM_Test01(const bool verbose) {
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

  *outStream << "\n"
      << "===============================================================================\n"
      << "|                                                                             |\n"
      << "|                      Unit Test HDIV_HEX_In_FEM                              |\n"
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

  typedef Basis_HDIV_HEX_In_FEM<DeviceSpaceType,OutValueType,PointValueType> HexBasisType;
  constexpr ordinal_type maxOrder = Parameters::MaxOrder ;
  constexpr ordinal_type dim = 3;

  *outStream
  << "\n"
  << "===============================================================================\n"
  << "| TEST 1: Exception testing                                                   |\n"
  << "===============================================================================\n";

  try {

#ifdef HAVE_INTREPID2_DEBUG

    ordinal_type nthrow = 0, ncatch = 0;

    const ordinal_type order = 1;
    HexBasisType hexBasis(order);

    // Array of reference hex nodes - used for evaluation of basis
    DynRankViewHostScalarValueType ConstructWithLabelScalar(hexNodesHost, 8, 3);

    hexNodesHost(0,0) = -1.0; hexNodesHost(0,1) = -1.0; hexNodesHost(0,2) = -1.0;
    hexNodesHost(1,0) =  1.0; hexNodesHost(1,1) = -1.0; hexNodesHost(1,2) = -1.0;
    hexNodesHost(2,0) = -1.0; hexNodesHost(2,1) =  1.0; hexNodesHost(2,2) = -1.0;
    hexNodesHost(3,0) =  1.0; hexNodesHost(3,1) =  1.0; hexNodesHost(3,2) = -1.0;
    hexNodesHost(4,0) = -1.0; hexNodesHost(4,1) = -1.0; hexNodesHost(4,2) =  1.0;
    hexNodesHost(5,0) =  1.0; hexNodesHost(5,1) = -1.0; hexNodesHost(5,2) =  1.0;
    hexNodesHost(6,0) = -1.0; hexNodesHost(6,1) =  1.0; hexNodesHost(6,2) =  1.0;
    hexNodesHost(7,0) =  1.0; hexNodesHost(7,1) =  1.0; hexNodesHost(7,2) =  1.0;

    auto hexNodes_scalar = Kokkos::create_mirror_view(typename DeviceSpaceType::memory_space(), hexNodesHost);
    Kokkos::deep_copy(hexNodes_scalar, hexNodesHost);

    DynRankViewPointValueType ConstructWithLabelPointView(hexNodes, 8, 3);
    RealSpaceTools<DeviceSpaceType>::clone(hexNodes, hexNodes_scalar);

    // Array dimensions
    const ordinal_type numFields = hexBasis.getCardinality();
    const ordinal_type numPoints = hexNodes.extent(0);
    const ordinal_type spaceDim  = hexBasis.getBaseCellTopology().getDimension();

    {
      DynRankViewOutValueType ConstructWithLabelOutView(vals, numFields, numPoints, spaceDim);

      // exception #1: GRAD cannot be applied to HDIV functions
      INTREPID2_TEST_ERROR_EXPECTED( hexBasis.getValues(vals, hexNodes, OPERATOR_GRAD) );

      // exception #2: CURL cannot be applied to HDIV functions
      INTREPID2_TEST_ERROR_EXPECTED( hexBasis.getValues(vals, hexNodes, OPERATOR_CURL) );
    }

    {
      // Exceptions 3-7: all bf tags/bf Ids below are wrong and should cause getDofOrdinal() and
      // getDofTag() to access invalid array elements thereby causing bounds check exception
      //
      INTREPID2_TEST_ERROR_EXPECTED( hexBasis.getDofOrdinal(3,0,0) );
      INTREPID2_TEST_ERROR_EXPECTED( hexBasis.getDofOrdinal(1,1,1) );
      INTREPID2_TEST_ERROR_EXPECTED( hexBasis.getDofOrdinal(0,4,1) );
      INTREPID2_TEST_ERROR_EXPECTED( hexBasis.getDofTag(12) );
      INTREPID2_TEST_ERROR_EXPECTED( hexBasis.getDofTag(-1) );
    }

    {
      // Exceptions 8-15 test exception handling with incorrectly dimensioned input/output arrays
      DynRankViewOutValueType ConstructWithLabelOutView(vals, numFields, numPoints, spaceDim);

      {
        // Exception #8: input points array must be of rank-2
        DynRankViewPointValueType ConstructWithLabelPointView(badPoints,4,5,3);
        INTREPID2_TEST_ERROR_EXPECTED( hexBasis.getValues(vals, badPoints, OPERATOR_VALUE));
      }

      {
        // Exception #9: dimension 1 in the input point array must equal space dimension of the cell
        DynRankViewPointValueType ConstructWithLabelPointView(badPoints,4,2);
        INTREPID2_TEST_ERROR_EXPECTED( hexBasis.getValues(vals, badPoints, OPERATOR_VALUE));
      }

      {
        // Exception #10: output values must be of rank-3 for OPERATOR_VALUE
        DynRankViewOutValueType ConstructWithLabelOutView(badVals,4,3);
        INTREPID2_TEST_ERROR_EXPECTED( hexBasis.getValues(badVals, hexNodes, OPERATOR_VALUE));
      }

      {
        // Exception #11: output values must be of rank-2 for OPERATOR_DIV
        DynRankViewOutValueType ConstructWithLabelOutView(badVals,4,3,3);
        INTREPID2_TEST_ERROR_EXPECTED( hexBasis.getValues(badVals, hexNodes, OPERATOR_DIV));
      }

      {
        // Exception #12: incorrect 0th dimension of output array (must equal number of basis functions)
        DynRankViewOutValueType ConstructWithLabelOutView(badVals,numFields+1,numPoints,3);
        INTREPID2_TEST_ERROR_EXPECTED( hexBasis.getValues(badVals, hexNodes, OPERATOR_VALUE));
      }

      {
        // Exception #13: incorrect 0th dimension of output array (must equal number of basis functions)
        DynRankViewOutValueType ConstructWithLabelOutView(badVals,numFields+1,numPoints);
        INTREPID2_TEST_ERROR_EXPECTED( hexBasis.getValues(badVals, hexNodes, OPERATOR_DIV));
      }

      {
        // Exception #14: incorrect 1st dimension of output array (must equal number of points)
        DynRankViewOutValueType ConstructWithLabelOutView(badVals,numFields,numPoints+1,3);
        INTREPID2_TEST_ERROR_EXPECTED( hexBasis.getValues(badVals, hexNodes, OPERATOR_VALUE));
      }

      {
        // Exception #15: incorrect 1st dimension of output array (must equal number of points)
        DynRankViewOutValueType ConstructWithLabelOutView(badVals,numFields,numPoints + 1);
        INTREPID2_TEST_ERROR_EXPECTED( hexBasis.getValues(badVals, hexNodes, OPERATOR_DIV));
      }
      {
        // Exception #16: incorrect 2nd dimension of output array (must equal the space dimension)
        DynRankViewOutValueType ConstructWithLabelOutView(badVals,numFields,numPoints,4);
        INTREPID2_TEST_ERROR_EXPECTED( hexBasis.getValues(badVals, hexNodes, OPERATOR_VALUE));
      }
    }

    if (nthrow != ncatch) {
      errorFlag++;
      *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
      *outStream << "# of catch ("<< ncatch << ") is different from # of throw (" << nthrow << ")\n";
    }
#endif
  } catch (std::exception err) {
    *outStream << "UNEXPECTED ERROR !!! ----------------------------------------------------------\n";
    *outStream << err.what() << '\n';
    *outStream << "-------------------------------------------------------------------------------" << "\n\n";
    errorFlag = -1000;
  };

  *outStream
  << "\n"
  << "===============================================================================\n"
  << "| TEST 2: Testing OPERATOR_VALUE (Kronecker property using dof coefficients)  |\n"
  << "===============================================================================\n";

  try {

    const ordinal_type order = std::min(3, maxOrder);
    HexBasisType hexBasis(order);

    const ordinal_type numFields = hexBasis.getCardinality();
    DynRankViewScalarValueType ConstructWithLabelScalar(dofCoords_scalar, numFields, dim);
    hexBasis.getDofCoords(dofCoords_scalar);

    DynRankViewScalarValueType ConstructWithLabelScalar(dofCoeffs, numFields, dim);
    hexBasis.getDofCoeffs(dofCoeffs);

    DynRankViewPointValueType ConstructWithLabelPointView(dofCoords, numFields , dim);
    RealSpaceTools<DeviceSpaceType>::clone(dofCoords,dofCoords_scalar);

    DynRankViewOutValueType ConstructWithLabelOutView(basisAtDofCoords, numFields, numFields, dim);
    hexBasis.getValues(basisAtDofCoords, dofCoords, OPERATOR_VALUE);

    auto h_dofCoords = Kokkos::create_mirror_view(dofCoords);
    Kokkos::deep_copy(h_dofCoords, dofCoords);

    auto h_dofCoeffs = Kokkos::create_mirror_view(dofCoeffs);
    Kokkos::deep_copy(h_dofCoeffs, dofCoeffs);

    auto h_basisAtDofCoords = Kokkos::create_mirror_view(basisAtDofCoords);
    Kokkos::deep_copy(h_basisAtDofCoords, basisAtDofCoords);

    for (int i=0;i<numFields;i++) {
      for (int j=0;j<numFields;j++) {

        OutValueType dofValue = 0.0;
        for(ordinal_type d=0;d<dim;++d)
          dofValue += h_basisAtDofCoords(i,j,d)*h_dofCoeffs(j,d);

        // check values
        const scalar_type expected_dofValue = (i == j);
        if (std::abs(dofValue - expected_dofValue) > tol) {
          errorFlag++;
          std::stringstream ss;
          ss << "\nValue of basis function " << i << " at (" << h_dofCoords(i,0) << ", " << h_dofCoords(i,1) << ", "<< h_dofCoords(i,2) << ") is " << dofValue << " but should be " << expected_dofValue << "\n";
          *outStream << ss.str();
        }
      }
    }
  } catch (std::exception err) {
    *outStream << "UNEXPECTED ERROR !!! ----------------------------------------------------------\n";
    *outStream << err.what() << '\n';
    *outStream << "-------------------------------------------------------------------------------" << "\n\n";
    errorFlag = -1000;
  };

  *outStream
  << "\n"
  << "===============================================================================\n"
  << "| TEST 3: Testing OPERATOR_VALUE (Kronecker property using tags)              |\n"
  << "===============================================================================\n";

  try {

    const ordinal_type order = std::min(3, maxOrder);
    HexBasisType hexBasis(order);

    shards::CellTopology hex_8(shards::getCellTopologyData<shards::Hexahedron<8> >());
    const ordinal_type numFields = hexBasis.getCardinality();
    DynRankViewScalarValueType ConstructWithLabelScalar(dofCoords_scalar, numFields, dim);
    hexBasis.getDofCoords(dofCoords_scalar);

    DynRankViewPointValueType ConstructWithLabelPointView(dofCoords, numFields , dim);
    RealSpaceTools<DeviceSpaceType>::clone(dofCoords,dofCoords_scalar);

    DynRankViewOutValueType ConstructWithLabelOutView(basisAtDofCoords, numFields, numFields, dim);
    hexBasis.getValues(basisAtDofCoords, dofCoords, OPERATOR_VALUE);

    auto h_basisAtDofCoords = Kokkos::create_mirror_view(basisAtDofCoords);
    Kokkos::deep_copy(h_basisAtDofCoords, basisAtDofCoords);

    // Face normals
    //   Note normals are indexed by face following Shards numbering, but
    //   directions are consistent with basis defintion and not Shards orientation
    DynRankViewHostScalarValueType ConstructWithLabelScalar(normals, numFields,dim); // normals at each point basis point
    normals(0,0)  =  0.0; normals(0,1)  =  1.0; normals(0,2)  =  0.0;
    normals(1,0)  =  1.0; normals(1,1)  =  0.0; normals(1,2)  =  0.0;
    normals(2,0)  =  0.0; normals(2,1)  =  1.0; normals(2,2)  =  0.0;
    normals(3,0)  =  1.0; normals(3,1)  =  0.0; normals(3,2)  =  0.0;
    normals(4,0)  =  0.0; normals(4,1)  =  0.0; normals(4,2)  =  1.0;
    normals(5,0)  =  0.0; normals(5,1)  =  0.0; normals(5,2)  =  1.0;

    const auto allTags = hexBasis.getAllDofTags();

    // test for Kronecker property
    for (int i=0;i<numFields;i++) {
      for (int j=0;j<numFields;j++) {

        // initialize
        OutValueType dofValue = 0;

        if(allTags(j,0) == dim-1) { //face
          auto faceId = allTags(j,1);
          for (ordinal_type k=0;k<dim; k++)
            dofValue += h_basisAtDofCoords(i,j,k)*normals(faceId,k);
        }
        else { //elem
          auto dofsPerDim = numFields/dim;
          auto k_ind = j/dofsPerDim;
          dofValue = h_basisAtDofCoords(i,j,k_ind);
        }

        // check values
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

  *outStream
  << "\n"
  << "===============================================================================\n"
  << "| TEST 4: Test correctness of basis function values and divs                  |\n"
  << "===============================================================================\n";

  outStream -> precision(20);


  // VALUE: Each row pair gives the 12x3 correct basis set values at an evaluation point: (P,F,D) layout
  const scalar_type basisValues[] = {
      // basis function 0 (in from x==-1 plane, y and z are constant functions)
      1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0,
      1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0,
      // basis function 1 (out from x==1 plane, y and z are constant functions
      0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0,
      0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0,
      // basis function 2 (in from y==-1 plane, x and z are constant functions
      0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0,
      0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0,
      // basis function 3 (out from y == 1 plane, x and z are constant function
      0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0,
      0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0,
      // basis function 4 (in from z == -1 plane, x and y are constant function
      0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 1,
      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
      // basis function 4 (out from z == 1 plane, x and y are constant function
      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 1
  };

  // DIV: each row gives the 6 correct values of the divergence of the 6 basis functions: (P,F) layout
  const scalar_type basisDivs[] = {   
      -0.5, -0.5, -0.5, -0.5, -0.5, -0.5, -0.5, -0.5,
      0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5,
      -0.5, -0.5, -0.5, -0.5, -0.5, -0.5, -0.5, -0.5,
      0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5,
      -0.5, -0.5, -0.5, -0.5, -0.5, -0.5, -0.5, -0.5,
      0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5
  };


  try {

    const ordinal_type order = 1;
    HexBasisType hexBasis(order);

    // Array of reference hex nodes - used for evaluation of basis
    DynRankViewHostScalarValueType ConstructWithLabelScalar(hexNodesHost, 8, 3);
    DynRankViewPointValueType ConstructWithLabelPointView(hexNodes, 8, 3);

    hexNodesHost(0,0) = -1.0; hexNodesHost(0,1) = -1.0; hexNodesHost(0,2) = -1.0;
    hexNodesHost(1,0) =  1.0; hexNodesHost(1,1) = -1.0; hexNodesHost(1,2) = -1.0;
    hexNodesHost(2,0) = -1.0; hexNodesHost(2,1) =  1.0; hexNodesHost(2,2) = -1.0;
    hexNodesHost(3,0) =  1.0; hexNodesHost(3,1) =  1.0; hexNodesHost(3,2) = -1.0;
    hexNodesHost(4,0) = -1.0; hexNodesHost(4,1) = -1.0; hexNodesHost(4,2) =  1.0;
    hexNodesHost(5,0) =  1.0; hexNodesHost(5,1) = -1.0; hexNodesHost(5,2) =  1.0;
    hexNodesHost(6,0) = -1.0; hexNodesHost(6,1) =  1.0; hexNodesHost(6,2) =  1.0;
    hexNodesHost(7,0) =  1.0; hexNodesHost(7,1) =  1.0; hexNodesHost(7,2) =  1.0;

    auto hexNodes_scalar = Kokkos::create_mirror_view(typename DeviceSpaceType::memory_space(), hexNodesHost);
    Kokkos::deep_copy(hexNodes_scalar, hexNodesHost);
    RealSpaceTools<DeviceSpaceType>::clone(hexNodes, hexNodes_scalar);

    // Array dimensions
    const ordinal_type numFields = hexBasis.getCardinality();
    const ordinal_type numPoints = hexNodes.extent(0);
    const ordinal_type spaceDim  = hexBasis.getBaseCellTopology().getDimension();

    *outStream << " -- Testing OPERATOR_VALUE \n";
    {
      // Check VALUE of basis functions
      DynRankViewOutValueType ConstructWithLabelOutView(vals, numFields, numPoints, spaceDim);
      hexBasis.getValues(vals, hexNodes, OPERATOR_VALUE);
      auto vals_host = Kokkos::create_mirror_view(vals);
      Kokkos::deep_copy(vals_host, vals);

      for (ordinal_type i = 0; i < numFields; i++) {
        for (ordinal_type j = 0; j < numPoints; j++) {
          for (ordinal_type k = 0; k < spaceDim; k++) {

            // compute offset for (P,F,D) data layout: indices are P->j, F->i, D->k
            const ordinal_type l = k + i * spaceDim * numPoints + j * spaceDim;
            if (std::abs(vals_host(i,j,k) - basisValues[l]) > tol) {
              errorFlag++;
              *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";

              // Output the multi-index of the value where the error is:
              *outStream << " At multi-index { ";
              *outStream << i << " ";*outStream << j << " ";*outStream << k << " ";
              *outStream << "}  computed value: " << vals(i,j,k)
                                << " but reference value: " << basisValues[l] << "\n";
            }
          }
        }
      }
    }
    *outStream << " -- Testing OPERATOR_DIV \n";
    {
      // Check DIV of basis function:
      DynRankViewOutValueType ConstructWithLabelOutView(divs, numFields, numPoints);
      hexBasis.getValues(divs, hexNodes, OPERATOR_DIV);
      auto divs_host = Kokkos::create_mirror_view(divs);
      Kokkos::deep_copy(divs_host, divs);

      for (ordinal_type i = 0; i < numFields; ++i) {
        for (ordinal_type j = 0; j < numPoints; ++j) {
          const ordinal_type l =  i * numPoints + j;
          if (std::abs(divs_host(i,j) - basisDivs[l]) > tol) {
            errorFlag++;
            *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";

            // Output the multi-index of the div where the error is:
            *outStream << " At multi-index { ";
            *outStream << i << " ";*outStream << j << " ";
            *outStream << "}  computed div component: " << divs(i,j)
                                << " but reference div component: " << basisDivs[l] << "\n";
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
