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
#include "Kokkos_DefaultNode.hpp"
#include "Kokkos_CrsMatrix.hpp"
#include "Kokkos_DefaultArithmetic.hpp"
#include "Kokkos_DefaultSparseMultiply.hpp"
#include "Kokkos_Version.hpp"

namespace {

  using std::endl;
  using Kokkos::MultiVector;
  using Kokkos::CrsMatrix;
  using Kokkos::CrsGraph;
  using Kokkos::DefaultArithmetic;
  using Kokkos::DefaultSparseMultiply;
  using Teuchos::ArrayRCP;
  using Teuchos::RCP;

  int N = 1000;

  typedef Kokkos::DefaultNode::DefaultNodeType Node;

  TEUCHOS_STATIC_SETUP()
  {
    Teuchos::CommandLineProcessor &clp = Teuchos::UnitTestRepository::getCLP();
    clp.addOutputSetupOptions(true);
    clp.setOption("test-size",&N,"Vector length for tests.");
  }

  //
  // UNIT TESTS
  // 

  TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( CrsMatrix, SparseMultiply, Scalar, Ordinal )
  {
    typedef CrsGraph<Ordinal,Node>  GRPH;
    typedef CrsMatrix<Scalar,Node>  MAT;
    typedef MultiVector<Scalar,Node> MV;
    typedef typename Node::size_t size_t;
    // generate tridiagonal matrix:
    // [ 1 -1                   ]
    // [-1  2  -1               ]
    // [   -1   2  -1           ]
    // [                        ]
    // [                -1  2 -1]
    // [                   -1  1]
    if (N<2) return;
    RCP<Node> node = Kokkos::DefaultNode::getDefaultNode();
    GRPH G(N,node);
    MAT  A(N,node);
    // allocate buffers for offsets, indices and values
    const size_t totalNNZ = 3*N - 2;
    ArrayRCP<size_t> offsets = node->template allocBuffer<size_t> (N+1);
    ArrayRCP<Ordinal>   inds = node->template allocBuffer<Ordinal>(totalNNZ);
    ArrayRCP<Scalar>    vals = node->template allocBuffer<Scalar >(totalNNZ);
    // fill the buffers on the host
    {
      ArrayRCP<size_t>  offsets_h = node->template viewBufferNonConst<size_t> (Kokkos::WriteOnly,N+1,offsets);
      ArrayRCP<Ordinal>    inds_h = node->template viewBufferNonConst<Ordinal>(Kokkos::WriteOnly,totalNNZ,inds);
      ArrayRCP<Scalar>     vals_h = node->template viewBufferNonConst<Scalar >(Kokkos::WriteOnly,totalNNZ,vals);
      int NNZsofar = 0;
      offsets_h[0] = NNZsofar;
      inds_h[NNZsofar] = 0; inds_h[NNZsofar+1] =  1;
      vals_h[NNZsofar] = 1; vals_h[NNZsofar+1] = -1;
      NNZsofar += 2;
      for (int i=1; i != N-1; ++i) {
        offsets_h[i] = NNZsofar;
        inds_h[NNZsofar] = i-1; inds_h[NNZsofar+1] = i; inds_h[NNZsofar+2] = i+1;
        vals_h[NNZsofar] =  -1; vals_h[NNZsofar+1] = 2; vals_h[NNZsofar+2] =  -1;
        NNZsofar += 3;
      }
      offsets_h[N-1] = NNZsofar;
      inds_h[NNZsofar] = N-2; inds_h[NNZsofar+1] = N-1;
      vals_h[NNZsofar] =  -1; vals_h[NNZsofar+1] = 1;
      NNZsofar += 2;
      offsets_h[N]   = NNZsofar;
      TEST_FOR_EXCEPT(NNZsofar != totalNNZ);
      inds_h    = Teuchos::null;
      vals_h    = Teuchos::null;
      offsets_h = Teuchos::null;
    }
    G.setPackedStructure(offsets, inds);
    A.setPackedValues(vals);
    DefaultSparseMultiply<Scalar,Ordinal,Node> dsm(node);
    Teuchos::DataAccess cv;
    cv = dsm.initializeStructure(G,Teuchos::View);
    out << "DefaultSparseMultiply::initializeStructure<CrsGraph>(G,View) returned "
        << (cv == Teuchos::View ? "View" : "Copy") << endl;
    cv = dsm.initializeValues(A,Teuchos::View);
    out << "DefaultSparseMultiply::initializeValues<CrsMatrix>(A,View) returned "
        << (cv == Teuchos::View ? "View" : "Copy") << endl;

    ArrayRCP<Scalar> xdat, axdat;
    xdat  = node->template allocBuffer<Scalar>(N);
    axdat = node->template allocBuffer<Scalar>(N);
    MV X(node), AX(node);
    X.initializeValues( N,1, xdat,N);
    AX.initializeValues(N,1,axdat,N);
    DefaultArithmetic<MV>::Init(X,1);
    dsm.multiply(Teuchos::NO_TRANS,1.0,X,0.0,AX);
    ArrayRCP<const Scalar> axview = node->template viewBuffer<Scalar>(N,axdat);
    Scalar err = 0.0;
    for (int i=0; i<N; ++i) {
      err = axview[i] * axview[i];
    }
    axview = Teuchos::null;
    err = Teuchos::ScalarTraits<Scalar>::squareroot(err);
    TEST_EQUALITY_CONST(err, 0.0);
    xdat = Teuchos::null;
    axdat = Teuchos::null;
  }

#define UNIT_TEST_GROUP_ORDINAL_SCALAR( ORDINAL, SCALAR ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( CrsMatrix, SparseMultiply, SCALAR, ORDINAL ) \

#    define UNIT_TEST_GROUP_ORDINAL( ORDINAL ) \
         UNIT_TEST_GROUP_ORDINAL_SCALAR(ORDINAL, int) \
         UNIT_TEST_GROUP_ORDINAL_SCALAR(ORDINAL, float)
     UNIT_TEST_GROUP_ORDINAL(int)
     typedef short int ShortInt; UNIT_TEST_GROUP_ORDINAL(ShortInt)

}

