// ************************************************************************
// 
//          Kokkos: A Fast Kernel Package
//              Copyright (2004) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ************************************************************************
//@HEADER

#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_TimeMonitor.hpp>
#include <Teuchos_Time.hpp>
#include <Teuchos_TypeNameTraits.hpp>
#include <Teuchos_ScalarTraits.hpp>

#include "Kokkos_ConfigDefs.hpp"
#include "Kokkos_DefaultArithmetic.hpp"
#include "Kokkos_DefaultSparseOps.hpp"
#include "Kokkos_CUSPSparseOps.hpp"
#include "Kokkos_LinAlgVersion.hpp"

#include "Kokkos_ThrustGPUNode.hpp"

#include <cuda.h>
#include <cusp/version.h>
#include <cusp/io/matrix_market.h>
#include <cusp/csr_matrix.h>

namespace {

  using Kokkos::MultiVector;
  using Kokkos::CrsMatrix;
  using Kokkos::CrsGraph;
  using Kokkos::DefaultArithmetic;
  using Kokkos::DefaultKernels;
  using Kokkos::CUSPSparseOps;

  using Kokkos::SerialNode;
  using Kokkos::ThrustGPUNode;

  using Teuchos::ArrayRCP;
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::null;

  using std::endl;

  RCP<SerialNode   > serialnode;
  RCP<ThrustGPUNode> thrustnode;
  cusp::crs_matrix<int, double, cusp::host_memory> cuspHostCSR;
  int deviceNumber;
  int verbose;
  std::string matrixFile;
  int printMatrix;  
  int numIters;

  TEUCHOS_STATIC_SETUP()
  {
    Teuchos::CommandLineProcessor &clp = Teuchos::UnitTestRepository::getCLP();
    verbose = 0;
    deviceNumber = 0;
    printMatrix = 0;
    numIters = 10;
    clp.addOutputSetupOptions(true);
    clp.setOption("device-num",&deviceNumber,"Device number for CUDA device.");
    clp.setOption("matrix-file",&matrixFile,"Filename for test matrix.");
    clp.setOption("verbose",&verbose,"Verbosity flag, zero for quiet.");
    clp.setOption("print-matrix",&printMatrix,"Print matrix.");
    clp.setOption("num-iters",&numIters,"Number of power method iterations.");
  }

  template <class Node>
  RCP<Node> getNode() {
    TEST_FOR_EXCEPT(true);
  }

  template <>
  RCP<SerialNode> getNode<SerialNode>() {
    if (serialnode == null) {
      Teuchos::ParameterList pl;
      serialnode = rcp(new SerialNode(pl));
    }
    return serialnode;
  }

  template <>
  RCP<ThrustGPUNode> getNode<ThrustGPUNode>() {
    if (thrustnode == null) {
      Teuchos::ParameterList pl;
      pl.set<int>("Device Number",deviceNumber);
      pl.set<int>("Verbose",verbose);
      thrustnode = rcp(new ThrustGPUNode(pl));
    }
    return thrustnode;
  }

  //
  // UNIT TESTS
  // 

  TEUCHOS_UNIT_TEST( AAAAA_Runs_First, SayHelloAndReadMatrix )
  {
    out << Kokkos::LinAlgVersion() << std::endl
        << "CUDA v"   << (CUDA_VERSION/1000)  << "." << (CUDA_VERSION / 10 % 100) << "." << (CUDA_VERSION % 10)     << std::endl
        << "Thrust v" << THRUST_MAJOR_VERSION << "." << THRUST_MINOR_VERSION      << "." << THRUST_SUBMINOR_VERSION << std::endl
        << "CUSP v"   << CUSP_MAJOR_VERSION   << "." << CUSP_MINOR_VERSION        << "." << CUSP_SUBMINOR_VERSION   << std::endl;
    // 
    out << "\nReading file: " << matrixFile << "..." << std::endl; 
    cusp::io::read_matrix_market_file(cuspHostCSR,matrixFile);
    if (printMatrix) {
      cusp::print_matrix(cuspHostCSR);
    }
  }

  TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( SparseOps, SimplePowerMethod, SparseOps, Node )
  {
    RCP<Node> node = getNode<Node>();
    typedef typename SparseOps<void,int,Node>  OPS;
    typedef CrsGraph<int,Node,OPS>            GRPH;
    typedef CrsMatrix<double,int,Node,OPS>     MAT;
    typedef MultiVector<double,Node>            MV;
    typedef Teuchos::ScalarTraits<double>       ST;
    typedef DefaultArithmetic<MV>             BLAS;
    const double ONE = ST::one(),
                ZERO = ST::zero();
    typename OPS::template rebind<double>::other dsm(node);
    out << "Testing with sparse ops: " << Teuchos::typeName(dsm) << std::endl;
    {
      // grab data from CUSP
      const int numRows = cuspHostCSR.row_offsets.size()-1;
      const int numNZ   = cuspHostCSR.column_indices.size();
      ArrayRCP<size_t>    rowOffsets = arcp<size_t>(numRows+1);
      // one is int, the other size_t. Therefore, we need to copy.
      std::copy( cuspHostCSR.row_offsets.begin(), cuspHostCSR.row_offsets.end(), rowOffsets.begin() );
      ArrayRCP<const int>    colIndices = arcp<const int>(cuspHostCSR.column_indices.data() );
      ArrayRCP<const double> values     = arcp<const int>(cuspHostCSR.values.data() );
      // create Kokkos CrsGraph and CrsMatrix objects, put data there
      GRPH G(numRows,node);
      MAT  A(G);
      G.set1DStructure(colIndices, rowOffsets(0,numRows), rowOffsets(1,numRows));
      A.set1DValues(values);
      A.finalize(true);
      // init sparse ops class
      dsm.initializeStructure(G);
      dsm.initializeValues(A);
    }
    MV z(node), q(node);
    z.initializeValues( numRows,1,node->template allocBuffer<double>(numRows),numRows);
    q.initializeValues(numRows,1,node->template allocBuffer<double>(numRows),numRows);
    BLAS::Init( z, ONE);
    BLAS::Init( q, ZERO);
    // power method, ten iterations
    for (int iter = 0; iter < numIters; ++iter) {
      double z_z = BLAS::Norm2Squared(z);
      dsm.multiply(Teuchos::NO_TRANS, ONE/ST::squareroot(z_z), q, z);
    }
    double lambda = BLAS::Norm2Squared(z);
    out << "After " << numIters << "power iterations, lambda == " << lambda << std::endl;
  }

  TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( SparseOps, SimplePowerMethod, DefaultDeviceSparseOps, SerialNode )
  TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( SparseOps, SimplePowerMethod, DefaultDeviceSparseOps, ThrustGPUNode )
  TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( SparseOps, SimplePowerMethod, CUSPSparseOps         , ThrustGPUNode )

}
