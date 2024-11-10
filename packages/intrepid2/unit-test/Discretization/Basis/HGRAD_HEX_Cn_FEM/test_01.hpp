// @HEADER
// *****************************************************************************
//                           Intrepid2 Package
//
// Copyright 2007 NTESS and the Intrepid2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/** \file test_01.cpp
    \brief  Unit tests for the Intrepid2::G_HEX_C2_FEM class.
    \author Created by P. Bochev, D. Ridzal, K. Peterson and Kyungjoo Kim.
 */

#include "Intrepid2_config.h"

#ifdef HAVE_INTREPID2_DEBUG
#define INTREPID2_TEST_FOR_DEBUG_ABORT_OVERRIDE_TO_CONTINUE
#endif

#include "Intrepid2_Types.hpp"
#include "Intrepid2_Utils.hpp"

#include "Intrepid2_PointTools.hpp"
#include "Intrepid2_HGRAD_HEX_Cn_FEM.hpp"

#include "Teuchos_oblackholestream.hpp"
#include "Teuchos_RCP.hpp"

#include "packages/intrepid2/unit-test/Discretization/Basis/Macros.hpp"
#include "packages/intrepid2/unit-test/Discretization/Basis/Setup.hpp"

namespace Intrepid2 {

namespace Test {

using HostSpaceType = Kokkos::DefaultHostExecutionSpace;

template<typename OutValueType, typename PointValueType, typename DeviceType>
int HGRAD_HEX_Cn_FEM_Test01(const bool verbose) {

  //! Create an execution space instance.
  const auto space = Kokkos::Experimental::partition_space(typename DeviceType::execution_space {}, 1)[0];

  //! Setup test output stream.
  Teuchos::RCP<std::ostream> outStream = setup_output_stream<DeviceType>(
    verbose, "Basis_HGRAD_HEX_Cn_FEM FEM", {
      "1) Conversion of Dof tags into Dof ordinals and back",
      "2) Basis values for VALUE, GRAD, and Dk operators"
  });

  Teuchos::oblackholestream oldFormatState;
  oldFormatState.copyfmt(std::cout);

  typedef Kokkos::DynRankView<PointValueType,DeviceType> DynRankViewPointValueType;
  typedef Kokkos::DynRankView<OutValueType,DeviceType> DynRankViewOutValueType;
  typedef typename ScalarTraits<OutValueType>::scalar_type scalar_type;
  typedef Kokkos::DynRankView<scalar_type, DeviceType> DynRankViewScalarValueType;
  typedef Kokkos::DynRankView<scalar_type, HostSpaceType> DynRankViewHostScalarValueType;

  const scalar_type tol = tolerence();
  int errorFlag = 0;

  typedef Basis_HGRAD_HEX_Cn_FEM<DeviceType,OutValueType,PointValueType> HexBasisType;
  constexpr ordinal_type maxOrder = Parameters::MaxOrder;


  *outStream
  << "\n"
  << "===============================================================================\n"
  << "| TEST 1: Basis creation, exceptions tests                                    |\n"
  << "===============================================================================\n";

  try {

#ifdef HAVE_INTREPID2_DEBUG
    ordinal_type nthrow = 0, ncatch = 0;
    constexpr ordinal_type order = 3;
    if(order < maxOrder) {
      HexBasisType hexBasis(order);

      // Define array containing array of nodes to evaluate
      DynRankViewPointValueType ConstructWithLabelPointView(hexNodes, 27, 3);

      // Generic array for the output values; needs to be properly resized depending on the operator type
      const ordinal_type numFields = hexBasis.getCardinality();
      const ordinal_type numPoints = hexNodes.extent(0);
      //const ordinal_type spaceDim  = hexBasis.getBaseCellTopology().getDimension();

      // exception 1 - 2: CURL and DIV is not supported.
      {
        DynRankViewOutValueType ConstructWithLabelOutView(vals, numFields, numPoints, 3);
        INTREPID2_TEST_ERROR_EXPECTED( hexBasis.getValues(vals, hexNodes, OPERATOR_CURL) );
        INTREPID2_TEST_ERROR_EXPECTED( hexBasis.getValues(vals, hexNodes, OPERATOR_DIV) );
      }

      // Exceptions 3-7: all bf tags/bf Ids below are wrong and should cause getDofOrdinal() and
      // getDofTag() to access invalid array elements thereby causing bounds check exception
      {
        INTREPID2_TEST_ERROR_EXPECTED( hexBasis.getDofOrdinal(3,10,0) );
        INTREPID2_TEST_ERROR_EXPECTED( hexBasis.getDofOrdinal(1,2,3) );
        INTREPID2_TEST_ERROR_EXPECTED( hexBasis.getDofOrdinal(0,4,1) );
        INTREPID2_TEST_ERROR_EXPECTED( hexBasis.getDofTag(numFields) );
        INTREPID2_TEST_ERROR_EXPECTED( hexBasis.getDofTag(-1) );
      }

      // Exceptions 8-18 test exception handling with incorrectly dimensioned input/output arrays
      // exception #8: input points array must be of rank-2
      {
        DynRankViewOutValueType ConstructWithLabelOutView(vals, numFields, numPoints);
        {
          // exception #8: input points array must be of rank-2
          DynRankViewPointValueType ConstructWithLabelPointView(badPoints, 4, 5, 3);
          INTREPID2_TEST_ERROR_EXPECTED( hexBasis.getValues(vals, badPoints, OPERATOR_VALUE) );
        }
        {
          // exception #9: dimension 1 in the input point array must equal space dimension of the cell
          DynRankViewPointValueType ConstructWithLabelPointView(badPoints, 4, 4);
          INTREPID2_TEST_ERROR_EXPECTED( hexBasis.getValues(vals, badPoints, OPERATOR_VALUE) );
        }
        {
          // exception #10: output values must be of rank-2 for OPERATOR_VALUE
          DynRankViewOutValueType ConstructWithLabelOutView(badVals, 4, 3, 1);
          INTREPID2_TEST_ERROR_EXPECTED( hexBasis.getValues(badVals, hexNodes, OPERATOR_VALUE) );
        }
        {
          DynRankViewOutValueType ConstructWithLabelOutView(badVals, 4, 3);

          // exception #11: output values must be of rank-3 for OPERATOR_GRAD
          INTREPID2_TEST_ERROR_EXPECTED( hexBasis.getValues(badVals, hexNodes, OPERATOR_GRAD) );

          // exception #12: output values must be of rank-3 for OPERATOR_CURL
          INTREPID2_TEST_ERROR_EXPECTED( hexBasis.getValues(badVals, hexNodes, OPERATOR_D1) );

          // exception #13: output values must be of rank-3 for OPERATOR_D2
          INTREPID2_TEST_ERROR_EXPECTED( hexBasis.getValues(badVals, hexNodes, OPERATOR_D2) );
        }
      }
      {
        // exception #14: incorrect 0th dimension of output array (must equal number of basis functions)
        DynRankViewOutValueType ConstructWithLabelOutView(badVals, hexBasis.getCardinality() + 1, hexNodes.extent(0));
        INTREPID2_TEST_ERROR_EXPECTED( hexBasis.getValues(badVals, hexNodes, OPERATOR_VALUE) );
      }
      {
        // exception #15: incorrect 1st dimension of output array (must equal number of points)
        DynRankViewOutValueType ConstructWithLabelOutView(badVals, hexBasis.getCardinality(), hexNodes.extent(0) + 1);
        INTREPID2_TEST_ERROR_EXPECTED( hexBasis.getValues(badVals, hexNodes, OPERATOR_VALUE) );
      }
      {
        // exception #16: incorrect 2nd dimension of output array (must equal spatial dimension)
        DynRankViewOutValueType ConstructWithLabelOutView(badVals, hexBasis.getCardinality(), hexNodes.extent(0), 2);
        INTREPID2_TEST_ERROR_EXPECTED( hexBasis.getValues(badVals, hexNodes, OPERATOR_GRAD) );
      }
      {
        DynRankViewOutValueType ConstructWithLabelOutView(badVals, hexBasis.getCardinality(), hexNodes.extent(0), 40);

        // exception #17: incorrect 2nd dimension of output array (must equal spatial dimension)
        INTREPID2_TEST_ERROR_EXPECTED( hexBasis.getValues(badVals, hexNodes, OPERATOR_D2) );

        // exception #18: incorrect 2nd dimension of output array (must equal spatial dimension)
        INTREPID2_TEST_ERROR_EXPECTED( hexBasis.getValues(badVals, hexNodes, OPERATOR_D3) );
      }
    }
    if (nthrow != ncatch) {
      errorFlag++;
      *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
      *outStream << "# of catch ("<< ncatch << ") is different from # of throw (" << nthrow << ")\n";
    }
#endif
  } catch (std::exception &err) {
    *outStream << "UNEXPECTED ERROR !!! ----------------------------------------------------------\n";
    *outStream << err.what() << '\n';
    *outStream << "-------------------------------------------------------------------------------" << "\n\n";
    errorFlag = -1000;
  };

  *outStream
  << "\n"
  << "===============================================================================\n"
  << "| TEST 2: correctness of tag to enum and enum to tag lookups                  |\n"
  << "===============================================================================\n";

  try {

    const ordinal_type order = std::min(5, maxOrder);
    HexBasisType hexBasis(order);

    const ordinal_type numFields = hexBasis.getCardinality();
    const auto allTags = hexBasis.getAllDofTags();

    // Loop over all tags, lookup the associated dof enumeration and then lookup the tag again
    const ordinal_type dofTagSize = allTags.extent(0);
    for (ordinal_type i=0;i<dofTagSize;++i) {
      const auto bfOrd = hexBasis.getDofOrdinal(allTags(i,0), allTags(i,1), allTags(i,2));

      const auto myTag = hexBasis.getDofTag(bfOrd);
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
      const auto myTag  = hexBasis.getDofTag(bfOrd);
      const auto myBfOrd = hexBasis.getDofOrdinal(myTag(0), myTag(1), myTag(2));
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
    << "| TEST 3: Testing Kronecker property of basis functions                                              |\n"
    << "===============================================================================\n";


    const ordinal_type order = std::min(3,maxOrder);
    HexBasisType hexBasis(order, POINTTYPE_WARPBLEND);
    constexpr ordinal_type dim=3;
    const ordinal_type basisCardinality = hexBasis.getCardinality();
    DynRankViewScalarValueType ConstructWithLabel(lattice_scalar, basisCardinality , dim);
    DynRankViewPointValueType ConstructWithLabelPointView(lattice, basisCardinality , dim);

    hexBasis.getDofCoords(lattice_scalar);
    RealSpaceTools<DeviceType>::clone(lattice,lattice_scalar);

    auto lattice_host = Kokkos::create_mirror_view(lattice);

    DynRankViewOutValueType ConstructWithLabelOutView(basisAtLattice, basisCardinality, basisCardinality);
    hexBasis.getValues(space, basisAtLattice, lattice, OPERATOR_VALUE);

    auto h_basisAtLattice = Kokkos::create_mirror_view(basisAtLattice);
    Kokkos::deep_copy(h_basisAtLattice, basisAtLattice);

    for(ordinal_type iface =0; iface<6; iface++) {
      auto numFaceDofs = hexBasis.getDofCount(2,iface);
      for(ordinal_type i=0; i<numFaceDofs; i++) {
        auto idof = hexBasis.getDofOrdinal(2,iface,i);
        for(ordinal_type j=0; j<numFaceDofs; j++) {
          auto jdof = hexBasis.getDofOrdinal(2,iface,j);
          if ( idof==jdof && std::abs( h_basisAtLattice(idof,jdof) - 1.0 ) > tol ) {
            errorFlag++;
            *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
            *outStream << " Basis function " << idof << " does not have unit value at its node (" << h_basisAtLattice(idof,jdof) <<")\n";
          }
          if ( i!=j && std::abs( h_basisAtLattice(idof,jdof) ) > tol ) {
            errorFlag++;
            *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
            *outStream << " Basis function " << idof << " does not vanish at node " << jdof << "\n";
            *outStream << " Basis function value is " << h_basisAtLattice(idof,jdof) << "\n";
          }
        }
      }
    }


    // test for Kronecker property
    for (int i=0;i<basisCardinality;i++) {
      for (int j=0;j<basisCardinality;j++) {
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

  *outStream
  << "\n"
  << "===============================================================================\n"
  << "| TEST 4: correctness of basis function values                                |\n"
  << "===============================================================================\n";

  outStream -> precision(20);

  // VALUE: Each row gives the 27 correct basis set values at an evaluation point
  const scalar_type basisValues[] = {
      1.000, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 1.000, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 1.000, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 1.000, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 1.000, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 1.000, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 1.000, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 1.000, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 1.000, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0, 1.000, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.000, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.000, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.000, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.000, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.000, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.000, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.000, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.000, 0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.000, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.000, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.000, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.000, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.000, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.000, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.000, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.000, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.000  };

  // GRAD, D1, D2, D3 and D4 test values are stored in files due to their large size
  std::string     fileName;
  std::ifstream   dataFile;

  // GRAD and D1 values are stored in (F,P,D) format in a data file. Read file and do the test
  std::vector<scalar_type> basisGrads;           // Flat array for the gradient values.
  {
    fileName = "./testdata/HEX_C2_GradVals.dat";
    dataFile.open(fileName.c_str());
    INTREPID2_TEST_FOR_EXCEPTION( !dataFile.good(), std::logic_error,
        ">>> ERROR (HGRAD_HEX_C2/test01): could not open GRAD values data file, test aborted.");
    while (!dataFile.eof() ){
      double temp;
      std::string line;                       // string for one line of input file
      std::getline(dataFile, line);           // get next line from file
      std::stringstream data_line(line);      // convert to stringstream
      while(data_line >> temp){               // extract value from line
        basisGrads.push_back(temp);           // push into vector
      }
    }
    dataFile.close();
    dataFile.clear();
  }

  //D2: flat array with the values of D2 applied to basis functions. Multi-index is (F,P,D2cardinality)
  std::vector<scalar_type> basisD2;
  {
    fileName = "./testdata/HEX_C2_D2Vals.dat";
    dataFile.open(fileName.c_str());
    INTREPID2_TEST_FOR_EXCEPTION( !dataFile.good(), std::logic_error,
        ">>> ERROR (HGRAD_HEX_C2/test01): could not open D2 values data file, test aborted.");
    while (!dataFile.eof() ){
      double temp;
      std::string line;                       // string for one line of input file
      std::getline(dataFile, line);           // get next line from file
      std::stringstream data_line(line);      // convert to stringstream
      while(data_line >> temp){               // extract value from line
        basisD2.push_back(temp);              // push into vector
      }
    }
    dataFile.close();
    dataFile.clear();
  }

  //D3: flat array with the values of D3 applied to basis functions. Multi-index is (F,P,D3cardinality)
  std::vector<scalar_type> basisD3;
  {
    fileName = "./testdata/HEX_C2_D3Vals.dat";
    dataFile.open(fileName.c_str());
    TEUCHOS_TEST_FOR_EXCEPTION( !dataFile.good(), std::logic_error,
        ">>> ERROR (HGRAD_HEX_C2/test01): could not open D3 values data file, test aborted.");

    while (!dataFile.eof() ){
      double temp;
      std::string line;                            // string for one line of input file
      std::getline(dataFile, line);           // get next line from file
      std::stringstream data_line(line);           // convert to stringstream
      while(data_line >> temp){               // extract value from line
        basisD3.push_back(temp);              // push into vector
      }
    }
    dataFile.close();
    dataFile.clear();
  }

  //D4: flat array with the values of D applied to basis functions. Multi-index is (F,P,D4cardinality)
  std::vector<scalar_type> basisD4;
  {
    fileName = "./testdata/HEX_C2_D4Vals.dat";
    dataFile.open(fileName.c_str());
    TEUCHOS_TEST_FOR_EXCEPTION( !dataFile.good(), std::logic_error,
        ">>> ERROR (HGRAD_HEX_C2/test01): could not open D4 values data file, test aborted.");

    while (!dataFile.eof() ){
      double temp;
      std::string line;                            // string for one line of input file
      std::getline(dataFile, line);           // get next line from file
      std::stringstream data_line(line);           // convert to stringstream
      while(data_line >> temp){               // extract value from line
        basisD4.push_back(temp);              // push into vector
      }
    }
    dataFile.close();
    dataFile.clear();
  }

  try{
    constexpr ordinal_type order = 2;
    if(order < maxOrder) {
      HexBasisType hexBasis(order);

      DynRankViewHostScalarValueType ConstructWithLabel(hexNodesHost, 27, 3);
      DynRankViewPointValueType ConstructWithLabelPointView(hexNodes, 27, 3);


      // do it lexicographically as a lattice
      hexNodesHost(0, 0) = -1.0;   hexNodesHost(0, 1) = -1.0;  hexNodesHost(0, 2) = -1.0;
      hexNodesHost(1, 0) =  0.0;   hexNodesHost(1, 1) = -1.0;  hexNodesHost(1, 2) = -1.0;
      hexNodesHost(2, 0) =  1.0;   hexNodesHost(2, 1) = -1.0;  hexNodesHost(2, 2) = -1.0;
      hexNodesHost(3, 0) = -1.0;   hexNodesHost(3, 1) =  0.0;  hexNodesHost(3, 2) = -1.0;
      hexNodesHost(4, 0) =  0.0;   hexNodesHost(4, 1) =  0.0;  hexNodesHost(4, 2) = -1.0;
      hexNodesHost(5, 0) =  1.0;   hexNodesHost(5, 1) =  0.0;  hexNodesHost(5, 2) = -1.0;
      hexNodesHost(6, 0) = -1.0;   hexNodesHost(6, 1) =  1.0;  hexNodesHost(6, 2) = -1.0;
      hexNodesHost(7, 0) = 0.0;    hexNodesHost(7, 1) =  1.0;  hexNodesHost(7, 2) = -1.0;
      hexNodesHost(8, 0) = 1.0;    hexNodesHost(8, 1) =  1.0;  hexNodesHost(8, 2) = -1.0;
      hexNodesHost(9, 0) = -1.0;   hexNodesHost(9, 1) = -1.0;  hexNodesHost(9, 2) = 0.0;
      hexNodesHost(10, 0) =  0.0;   hexNodesHost(10, 1) = -1.0;  hexNodesHost(10, 2) = 0.0;
      hexNodesHost(11, 0) =  1.0;   hexNodesHost(11, 1) = -1.0;  hexNodesHost(11, 2) = 0.0;
      hexNodesHost(12, 0) = -1.0;   hexNodesHost(12, 1) =  0.0;  hexNodesHost(12, 2) = 0.0;
      hexNodesHost(13, 0) =  0.0;   hexNodesHost(13, 1) =  0.0;  hexNodesHost(13, 2) = 0.0;
      hexNodesHost(14, 0) =  1.0;   hexNodesHost(14, 1) =  0.0;  hexNodesHost(14, 2) = 0.0;
      hexNodesHost(15, 0) = -1.0;   hexNodesHost(15, 1) =  1.0;  hexNodesHost(15, 2) = 0.0;
      hexNodesHost(16, 0) = 0.0;    hexNodesHost(16, 1) =  1.0;  hexNodesHost(16, 2) = 0.0;
      hexNodesHost(17, 0) = 1.0;    hexNodesHost(17, 1) =  1.0;  hexNodesHost(17, 2) = 0.0;
      hexNodesHost(18, 0) = -1.0;   hexNodesHost(18, 1) = -1.0;  hexNodesHost(18, 2) = 1.0;
      hexNodesHost(19, 0) =  0.0;   hexNodesHost(19, 1) = -1.0;  hexNodesHost(19, 2) = 1.0;
      hexNodesHost(20, 0) =  1.0;   hexNodesHost(20, 1) = -1.0;  hexNodesHost(20, 2) = 1.0;
      hexNodesHost(21, 0) = -1.0;   hexNodesHost(21, 1) =  0.0;  hexNodesHost(21, 2) = 1.0;
      hexNodesHost(22, 0) =  0.0;   hexNodesHost(22, 1) =  0.0;  hexNodesHost(22, 2) = 1.0;
      hexNodesHost(23, 0) =  1.0;   hexNodesHost(23, 1) =  0.0;  hexNodesHost(23, 2) = 1.0;
      hexNodesHost(24, 0) = -1.0;   hexNodesHost(24, 1) =  1.0;  hexNodesHost(24, 2) = 1.0;
      hexNodesHost(25, 0) = 0.0;    hexNodesHost(25, 1) =  1.0;  hexNodesHost(25, 2) = 1.0;
      hexNodesHost(26, 0) = 1.0;    hexNodesHost(26, 1) =  1.0;  hexNodesHost(26, 2) = 1.0;

      auto hexNodes_scalar = Kokkos::create_mirror_view(typename DeviceType::memory_space(), hexNodesHost);
      Kokkos::deep_copy(hexNodes_scalar, hexNodesHost);

      RealSpaceTools<DeviceType>::clone(hexNodes, hexNodes_scalar);

      // Dimensions for the output arrays:
      const ordinal_type numFields = hexBasis.getCardinality();
      const ordinal_type numPoints = hexNodes.extent(0);
      const ordinal_type spaceDim  = hexBasis.getBaseCellTopology().getDimension();
      const ordinal_type D2Cardin  = getDkCardinality(OPERATOR_D2, spaceDim);
      const ordinal_type D3Cardin  = getDkCardinality(OPERATOR_D3, spaceDim);
      const ordinal_type D4Cardin  = getDkCardinality(OPERATOR_D4, spaceDim);

      *outStream << " -- Testing OPERATOR_VALUE \n";
      {
        // Check VALUE of basis functions: resize vals to rank-2 container:
        DynRankViewOutValueType ConstructWithLabelOutView(vals, numFields, numPoints);
        hexBasis.getValues(space, vals, hexNodes, OPERATOR_VALUE);
        auto vals_host = Kokkos::create_mirror_view(vals);
        Kokkos::deep_copy(vals_host, vals);
        for (ordinal_type i = 0; i < numFields; ++i) {
          for (ordinal_type j = 0; j < numPoints; ++j) {

            // Compute offset for (F,P) container
            const ordinal_type l = j + i * numPoints;
            if (std::abs(vals_host(i,j) - basisValues[l]) > tol) {
              errorFlag++;
              *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";

              // Output the multi-index of the value where the error is:
              *outStream << " At multi-index { ";
              *outStream << i << " ";*outStream << j << " ";
              *outStream << "}  computed value: " << vals_host(i,j)
                                   << " but reference value: " << basisValues[l] << "\n";
            }
          }
        }
      }

      *outStream << " -- Testing OPERATOR_GRAD \n";
      {
        // Check GRAD of basis function: resize vals to rank-3 container
        DynRankViewOutValueType ConstructWithLabelOutView(vals, numFields, numPoints, spaceDim);
        hexBasis.getValues(space, vals, hexNodes, OPERATOR_GRAD);
        auto vals_host = Kokkos::create_mirror_view(vals);
        Kokkos::deep_copy(vals_host, vals);
        for (ordinal_type i = 0; i < numFields; ++i) {
          for (ordinal_type j = 0; j < numPoints; ++j) {
            for (ordinal_type k = 0; k < spaceDim; ++k) {

              // basisGrads is (F,P,D), compute offset:
              const ordinal_type l = k + j * spaceDim + i * spaceDim * numPoints;
              if (std::abs(vals_host(i,j,k) - basisGrads[l]) > tol) {
                errorFlag++;
                *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";

                // Output the multi-index of the value where the error is:
                *outStream << " At multi-index { ";
                *outStream << i << " ";*outStream << j << " ";*outStream << k << " ";
                *outStream << "}  computed grad component: " << vals_host(i,j,k)
                                     << " but reference grad component: " << basisGrads[l] << "\n";
              }
            }
          }
        }
      }

      *outStream << " -- Testing OPERATOR_D1 \n";
      {
        // Check GRAD of basis function: resize vals to rank-3 container
        DynRankViewOutValueType ConstructWithLabelOutView(vals, numFields, numPoints, spaceDim);
        hexBasis.getValues(space, vals, hexNodes, OPERATOR_D1);
        auto vals_host = Kokkos::create_mirror_view(vals);
        Kokkos::deep_copy(vals_host, vals);
        for (ordinal_type i = 0; i < numFields; ++i) {
          for (ordinal_type j = 0; j < numPoints; ++j) {
            for (ordinal_type k = 0; k < spaceDim; ++k) {

              // basisGrads is (F,P,D), compute offset:
              const ordinal_type l = k + j * spaceDim + i * spaceDim * numPoints;
              if (std::abs(vals_host(i,j,k) - basisGrads[l]) > tol) {
                errorFlag++;
                *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";

                // Output the multi-index of the value where the error is:
                *outStream << " At multi-index { ";
                *outStream << i << " ";*outStream << j << " ";*outStream << k << " ";
                *outStream << "}  computed grad component: " << vals_host(i,j,k)
                                     << " but reference grad component: " << basisGrads[l] << "\n";
              }
            }
          }
        }
      }

      *outStream << " -- Testing OPERATOR_D2 \n";
      {
        // Check GRAD of basis function: resize vals to rank-3 container
        DynRankViewOutValueType ConstructWithLabelOutView(vals, numFields, numPoints, D2Cardin);
        hexBasis.getValues(space, vals, hexNodes, OPERATOR_D2);
        auto vals_host = Kokkos::create_mirror_view(vals);
        Kokkos::deep_copy(vals_host, vals);
        for (ordinal_type i = 0; i < numFields; ++i) {
          for (ordinal_type j = 0; j < numPoints; ++j) {
            for (ordinal_type k = 0; k < D2Cardin; ++k) {

              // basisGrads is (F,P,D), compute offset:
              const ordinal_type l = k + j * D2Cardin + i * D2Cardin * numPoints;
              if (std::abs(vals_host(i,j,k) - basisD2[l]) > tol) {
                errorFlag++;
                *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";

                // Output the multi-index of the value where the error is:
                *outStream << " At multi-index { ";
                *outStream << i << " ";*outStream << j << " ";*outStream << k << " ";
                *outStream << "}  computed grad component: " << vals_host(i,j,k)
                                     << " but reference grad component: " << basisD2[l] << "\n";
              }
            }
          }
        }
      }

      *outStream << " -- Testing OPERATOR_D3 \n";
      {
        // Check GRAD of basis function: resize vals to rank-3 container
        DynRankViewOutValueType ConstructWithLabelOutView(vals, numFields, numPoints, D3Cardin);
        hexBasis.getValues(space, vals, hexNodes, OPERATOR_D3);
        auto vals_host = Kokkos::create_mirror_view(vals);
        Kokkos::deep_copy(vals_host, vals);
        for (ordinal_type i = 0; i < numFields; ++i) {
          for (ordinal_type j = 0; j < numPoints; ++j) {
            for (ordinal_type k = 0; k < D3Cardin; ++k) {

              // basisGrads is (F,P,D), compute offset:
              const ordinal_type l = k + j * D3Cardin + i * D3Cardin * numPoints;
              if (std::abs(vals_host(i,j,k) - basisD3[l]) > tol) {
                errorFlag++;
                *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";

                // Output the multi-index of the value where the error is:
                *outStream << " At multi-index { ";
                *outStream << i << " ";*outStream << j << " ";*outStream << k << " ";
                *outStream << "}  computed grad component: " << vals_host(i,j,k)
                                     << " but reference grad component: " << basisD3[l] << "\n";
              }
            }
          }
        }
      }

      *outStream << " -- Testing OPERATOR_D4 \n";
      {
        // Check GRAD of basis function: resize vals to rank-3 container
        DynRankViewOutValueType ConstructWithLabelOutView(vals, numFields, numPoints, D4Cardin);
        hexBasis.getValues(space, vals, hexNodes, OPERATOR_D4);
        auto vals_host = Kokkos::create_mirror_view(vals);
        Kokkos::deep_copy(vals_host, vals);
        for (ordinal_type i = 0; i < numFields; ++i) {
          for (ordinal_type j = 0; j < numPoints; ++j) {
            for (ordinal_type k = 0; k < D4Cardin; ++k) {

              // basisGrads is (F,P,D), compute offset:
              const ordinal_type l = k + j * D4Cardin + i * D4Cardin * numPoints;
              if (std::abs(vals_host(i,j,k) - basisD4[l]) > tol) {
                errorFlag++;
                *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";

                // Output the multi-index of the value where the error is:
                *outStream << " At multi-index { ";
                *outStream << i << " ";*outStream << j << " ";*outStream << k << " ";
                *outStream << "}  computed grad component: " << vals_host(i,j,k)
                                     << " but reference grad component: " << basisD4[l] << "\n";
              }
            }
          }
        }
      }

      // // Check D7 to D10 - must be zero. This basis does not cover D5 and D6
      const EOperator ops[4] = { OPERATOR_D7,
          OPERATOR_D8,
          OPERATOR_D9,
          OPERATOR_D10 };

      for (ordinal_type oid=0;oid<4;++oid) {
        *outStream << " -- Testing OPERATOR_D" << (oid+7) << "\n";
        const auto op = ops[oid];
        const ordinal_type DkCardin = Intrepid2::getDkCardinality(op, spaceDim);

        DynRankViewOutValueType ConstructWithLabelOutView(vals, numFields, numPoints, DkCardin);
        hexBasis.getValues(space, vals, hexNodes, op);
        auto vals_host = Kokkos::create_mirror_view(vals);
        Kokkos::deep_copy(vals_host, vals);
        for (ordinal_type i = 0; i < numFields; ++i) {
          for (ordinal_type j = 0; j < numPoints; ++j) {
            for (ordinal_type k = 0; k < DkCardin; ++k) {
              // basisGrads is (F,P,D), compute offset:
              if (std::abs(vals_host(i,j,k)) > tol) {
                errorFlag++;
                *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";

                // Output the multi-index of the value where the error is:
                *outStream << " At multi-index { ";
                *outStream << i << " ";*outStream << j << " ";*outStream << k << " ";
                *outStream << "}  computed grad component: " << vals_host(i,j,k)
                                     << " but reference grad component: 0.0\n";
              }
            }
          }
        }
      }
    }
  } catch (std::exception &err) {
    *outStream << err.what() << "\n\n";
    errorFlag = -1000;
  };

  *outStream
  << "\n"
  << "===============================================================================\n"
  << "| TEST 5: Function Space is Correct                                           |\n"
  << "===============================================================================\n";

  try {
    constexpr ordinal_type order = 2;
    HexBasisType hexBasis(order);

    const EFunctionSpace fs = hexBasis.getFunctionSpace();

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








