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

#include "Kokkos_ConfigDefs.hpp"
#include "Kokkos_DefaultNode.hpp"
#include "Kokkos_CrsMatrix.hpp"
#include "Kokkos_DefaultArithmetic.hpp"
#include "Kokkos_DefaultSparseMultiply.hpp"
#include "Kokkos_Version.hpp"

namespace {

  using Kokkos::MultiVector;
  using Kokkos::CrsMatrix;
  using Kokkos::DefaultArithmetic;
  using Kokkos::DefaultSparseMultiply;
  using Kokkos::size_type;

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
    typedef CrsMatrix<Scalar,Ordinal,Node>  MAT;
    typedef MultiVector<Scalar,Ordinal,Node> MV;
    // generate tridiagonal matrix:
    // [ 1 -1                   ]
    // [-1  2  -1               ]
    // [   -1   2  -1           ]
    // [                        ]
    // [                -1  2 -1]
    // [                   -1  1]
    if (N<2) return;
    MAT A;
    Node &node = A.getNode();
    TEST_EQUALITY_CONST(A.getNumRows(), 0);
    TEST_EQUALITY_CONST(A.getNumEntries(), 0);
    std::vector<size_type> NNZperRow(N);
    NNZperRow[0] = 2;
    for (int i=1; i<N-1; ++i) NNZperRow[i] = 3;
    NNZperRow[N-1] = 2;
    // fill matrix
    {
      Ordinal inds[3];
      Scalar vals[] = {-1,2,-1};
      A.initializeProfile(N,&NNZperRow[0]);
      for (int i=1; i<N-1; ++i) {
        inds[0] = i-1; inds[1] = i; inds[2] = i+1;
        A.insertEntries(i,3,inds,vals);
      }
      vals[1] = 1;
      inds[0] = 0;   inds[1] = 1;   A.insertEntries(0  ,2,inds,vals+1);
      inds[0] = N-2; inds[1] = N-1; A.insertEntries(N-1,2,inds,vals  );
    }
    DefaultSparseMultiply<MAT,MV> dsm(node);
    dsm.initializeStructure(A);
    dsm.initializeValues(A);

    typename Node::template buffer<Scalar>::buffer_t xdat, axdat;
    xdat  = node.template allocBuffer<Scalar>(N);
    axdat = node.template allocBuffer<Scalar>(N);
    MV X, AX;
    X.initializeValues(N,1,xdat,N);
    AX.initializeValues(N,1,axdat,N);
    DefaultArithmetic<MV>::Init(1.0,X);
    dsm.Apply(X,AX);
    node.template freeBuffer<Scalar>(xdat);
    node.template freeBuffer<Scalar>(axdat);
  }

#define UNIT_TEST_GROUP_ORDINAL_SCALAR( ORDINAL, SCALAR ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( CrsMatrix, SparseMultiply, SCALAR, ORDINAL ) \

#    define UNIT_TEST_GROUP_ORDINAL( ORDINAL ) \
         UNIT_TEST_GROUP_ORDINAL_SCALAR(ORDINAL, int) \
         UNIT_TEST_GROUP_ORDINAL_SCALAR(ORDINAL, float)
     UNIT_TEST_GROUP_ORDINAL(int)
     typedef short int ShortInt; UNIT_TEST_GROUP_ORDINAL(ShortInt)

}

