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
#include "Kokkos_Version.hpp"

#include <vector>

namespace {

  using Kokkos::DefaultNode;
  using Kokkos::CrsMatrix;
  using Kokkos::size_type;
  using Teuchos::ArrayRCP;

  typedef Kokkos::DefaultNode::DefaultNodeType Node;

  int N = 100;

  TEUCHOS_STATIC_SETUP()
  {
    Teuchos::CommandLineProcessor &clp = Teuchos::UnitTestRepository::getCLP();
    clp.addOutputSetupOptions(true);
    clp.setOption("test-size",&N,"Vector length for tests.");
  }

  //
  // UNIT TESTS
  // 

  // check that default constructor zeros out
  TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( CrsMatrix, DefaultConstructor, Scalar, Ordinal )
  {
    CrsMatrix<Scalar,Ordinal,Node> A;
    TEST_EQUALITY_CONST(A.getNumRows(), 0);
    TEST_EQUALITY_CONST(A.getNumEntries(), 0);
  }

  // intialize the structure data
  TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( CrsMatrix, InitializeProfile, Scalar, Ordinal )
  {
    // constant nnz per row
    {
      const size_type NNZperRow = 5;
      CrsMatrix<Scalar,Ordinal,Node> A;
      TEST_EQUALITY_CONST(A.getNumRows(), 0);
      TEST_EQUALITY_CONST(A.getNumEntries(), 0);
      A.initializeProfile(N,NNZperRow);
      TEST_EQUALITY(A.getNumRows(), N);
      TEST_EQUALITY(A.getNumEntries(), N*NNZperRow);
    }
    // variable nnz per row
    {
      CrsMatrix<Scalar,Ordinal,Node> A;
      TEST_EQUALITY_CONST(A.getNumRows(), 0);
      TEST_EQUALITY_CONST(A.getNumEntries(), 0);
      // something interesting...
      // NNZperRow = {0, 1, 2, 3, 4, 5, 0, 1, 2, 3, 4, 5, ...}
      std::vector<size_type> NNZperRow(N);
      size_type expNNZ = 0;
      for (int i=0; i < N; ++i) {NNZperRow[i] = i%6; expNNZ += NNZperRow[i];}
      A.initializeProfile(N,&NNZperRow[0]);
      TEST_EQUALITY(A.getNumRows(), N);
      TEST_EQUALITY(A.getNumEntries(), expNNZ);
    }
  }

  // tridiagonal matrix
  TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( CrsMatrix, TridiagonalMatrix, Scalar, Ordinal )
  {
    if (N<2) return;
    CrsMatrix<Scalar,Ordinal,Node> A;
    TEST_EQUALITY_CONST(A.getNumRows(), 0);
    TEST_EQUALITY_CONST(A.getNumEntries(), 0);
    Node &node = A.getNode();
    std::vector<size_type> NNZperRow(N);
    NNZperRow[0] = 2;
    for (int i=1; i<N-1; ++i) NNZperRow[i] = 3;
    NNZperRow[N-1] = 2;
    size_type expNNZ = 4 + (N-2)*3;
    std::vector<Ordinal> expInds;
    std::vector<Scalar>  expVals;
    {
      expInds.reserve(expNNZ);
      expVals.reserve(expNNZ);
      Scalar vals[] = {-1,1,-1};
      Ordinal inds[3];
      A.initializeProfile(N,&NNZperRow[0]);
      inds[0] = 0; inds[1] = 1; A.insertEntries(0,2,inds,vals+1);
      expInds.insert(expInds.end(), inds, inds+2);
      expVals.insert(expVals.end(), vals+1, vals+3);
      for (int i=1; i<N-1; ++i) {
        inds[0] = i-1; inds[1] = i; inds[2] = i+1;
        A.insertEntries(i,3,inds,vals);
        expInds.insert(expInds.end(), inds, inds+3);
        expVals.insert(expVals.end(), vals, vals+3);
      }
      inds[0] = N-2; inds[1] = N-1; A.insertEntries(N-1,2,inds,vals);
      expInds.insert(expInds.end(), inds, inds+2);
      expVals.insert(expVals.end(), vals, vals+2);
    }
    TEST_EQUALITY(A.getNumRows(), N);
    TEST_EQUALITY(A.getNumEntries(), expNNZ);
    ArrayRCP<const Scalar> actVals = node.template viewBufferConst<Scalar >(expNNZ, A.const_values() ,0);
    ArrayRCP<const Ordinal> actInds = node.template viewBufferConst<Ordinal>(expNNZ, A.const_indices(),0);
    for (size_type i=0; i<expNNZ; ++i) {
      TEST_EQUALITY(expInds[i], actInds[i]);
      TEST_EQUALITY(expVals[i], actVals[i]);
    }
    actVals = Teuchos::null;
    actInds = Teuchos::null;
  }

 #define UNIT_TEST_GROUP_ORDINAL_SCALAR( ORDINAL, SCALAR ) \
       TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( CrsMatrix, DefaultConstructor,  SCALAR, ORDINAL ) \
       TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( CrsMatrix, InitializeProfile , SCALAR, ORDINAL ) \
       TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( CrsMatrix, TridiagonalMatrix , SCALAR, ORDINAL )
 
 #define UNIT_TEST_GROUP_ORDINAL( ORDINAL ) \
          UNIT_TEST_GROUP_ORDINAL_SCALAR(ORDINAL, int) \
          UNIT_TEST_GROUP_ORDINAL_SCALAR(ORDINAL, float)
      UNIT_TEST_GROUP_ORDINAL(int)
      typedef short int ShortInt; UNIT_TEST_GROUP_ORDINAL(ShortInt)
 
}
 
