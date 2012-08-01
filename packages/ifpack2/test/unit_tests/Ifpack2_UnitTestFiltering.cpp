// @HEADER
// ***********************************************************************
// 
//                IFPACK
//                 Copyright (2004) Sandia Corporation
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
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
// @HEADER
#include <Teuchos_ConfigDefs.hpp>
#include <Ifpack2_ConfigDefs.hpp>
#include <Teuchos_UnitTestHarness.hpp>
#include <Ifpack2_Version.hpp>
#include <iostream>

#ifdef HAVE_IFPACK2_QD
#include <qd/dd_real.h>
#endif

#include <Ifpack2_UnitTestHelpers.hpp>
#ifdef HAVE_MPI
#include <Teuchos_DefaultMpiComm.hpp>
#else
#include <Teuchos_DefaultSerialComm.hpp>
#endif
#include <Tpetra_Map.hpp>
#include <Tpetra_CrsMatrix.hpp>
#include <Ifpack2_DiagonalFilter.hpp>
#include <Ifpack2_LocalFilter.hpp>
//#include "Ifpack2_DropFilter.h"
//#include "Ifpack2_SparsityFilter.h"
//#include "Ifpack2_SingletonFilter.h"
#include <Teuchos_RefCountPtr.hpp>

#include <Teuchos_FancyOStream.hpp>



using Tpetra::global_size_t;
typedef tif_utest::Node Node;
using namespace std;
using Teuchos::rcp;

//this macro declares the unit-test-class:
TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL(Ifpack2Filtering, Test0, Scalar, LocalOrdinal, GlobalOrdinal)
{
  std::string version = Ifpack2::Version();
  out << "Ifpack2::Version(): " << version << std::endl;
  
  // Useful matrices and such (tridiagonal test)
  global_size_t num_rows_per_proc = 5;
  const Teuchos::RCP<const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> > rowmap = tif_utest::create_tpetra_map<LocalOrdinal,GlobalOrdinal,Node>(num_rows_per_proc); 
  Teuchos::RCP<const Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> > Matrix = tif_utest::create_test_matrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>(rowmap);
  Tpetra::Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node> x(rowmap), y(rowmap), z(rowmap), a(rowmap),b(rowmap);


  // ====================================== //
  // create a new matrix, diagonally filtered //
  // ====================================== //
  Scalar alpha = 100.0;
  Scalar beta = 1.0;
  Ifpack2::DiagonalFilter<Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> > DropA(Matrix,alpha,beta);
  Teuchos::ScalarTraits<double>::seedrandom(24601);
  x.randomize();

  // Apply w/ Filter
  DropA.apply(x,y);

  // Apply manually
  Matrix->getLocalDiagCopy(b); 
  Matrix->apply(x,z);
  z.update(alpha,x,beta-1,b,1.0);

  // Diff
  TEST_COMPARE_FLOATING_ARRAYS(y.get1dView(), z.get1dView(), 1e4*Teuchos::ScalarTraits<Scalar>::eps());


  // ====================================== //
  // create a new matrix, locally filtered  //
  // ====================================== //
  Ifpack2::LocalFilter<Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> > LocalA(Matrix);
  Teuchos::ScalarTraits<double>::seedrandom(24601);
  x.randomize();


  // Apply w/ filter
  LocalA.apply(x,y);

  // Apply w/ GetRow
  size_t max_nz_per_row=LocalA.getNodeMaxNumRowEntries();
  Teuchos::Array<LocalOrdinal> Indices(max_nz_per_row);
  Teuchos::Array<Scalar> Values(max_nz_per_row);
  Teuchos::ArrayRCP<const Scalar> xview=x.get1dView();

  for(LocalOrdinal i=0; i < (LocalOrdinal)num_rows_per_proc; i++){
    size_t NumEntries;
    LocalA.getLocalRowCopy(i,Indices(),Values(),NumEntries);
    Scalar sum=0;
    for(LocalOrdinal j=0; (size_t) j < NumEntries; j++){
      sum+=Values[j] * xview[Indices[j]];
    }
    z.replaceLocalValue(i,sum);
  }
    
  // Diff
  TEST_COMPARE_FLOATING_ARRAYS(y.get1dView(), z.get1dView(), 1e4*Teuchos::ScalarTraits<Scalar>::eps());




  // ====================================== //
  // create a new matrix, dropping by value //
  // ====================================== //
  //
  // drop all elements below 4.0. Only the upper-right element
  // is maintained, plus all the diagonals that are not
  // considering in dropping.
  //  Ifpack_DropFilter DropA(Matrix,4.0);
  //  assert (DropA.MaxNumEntries() == 2);

  //  cout << "Sparsity, dropping by value" << endl;
  //  Ifpack_PrintSparsity_Simple(DropA);

  // ========================================= //
  // create a new matrix, dropping by sparsity //
  // ========================================= //
  //
  // Mantain 2 off-diagonal elements.
  //  Ifpack_SparsityFilter SparsityA(Matrix,2);

  //  cout << "Sparsity, dropping by sparsity" << endl;
  //  Ifpack_PrintSparsity_Simple(SparsityA);
  //  assert (SparsityA.MaxNumEntries() == 3);

  // ======================================== //
  // create new matrices, dropping singletons //
  // ======================================== //
  //
  // If we apply this filter NumPoints - 1 times, 
  // we end up with a one-row matrix
  //  Ifpack_SingletonFilter Filter1(Matrix);
  //  Ifpack_SingletonFilter Filter2(Teuchos::rcp(&Filter1, false));
  //  Ifpack_SingletonFilter Filter3(Teuchos::rcp(&Filter2, false));
  //  Ifpack_SingletonFilter Filter4(Teuchos::rcp(&Filter3, false));

  //  cout << "Sparsity, dropping singletons 4 times" << endl;
  //  Ifpack_PrintSparsity_Simple(Filter4);
  //  assert (Filter4.NumMyRows() == 1);
}



#define UNIT_TEST_GROUP_SCALAR_ORDINAL(Scalar,LocalOrdinal,GlobalOrdinal) \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( Ifpack2Filtering, Test0, Scalar, LocalOrdinal,GlobalOrdinal)

UNIT_TEST_GROUP_SCALAR_ORDINAL(double, int, int)
#ifndef HAVE_IFPACK2_EXPLICIT_INSTANTIATION
UNIT_TEST_GROUP_SCALAR_ORDINAL(float, short, int)
#endif
