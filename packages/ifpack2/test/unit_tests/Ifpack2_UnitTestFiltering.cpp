// @HEADER
// *****************************************************************************
//       Ifpack2: Templated Object-Oriented Algebraic Preconditioner Package
//
// Copyright 2009 NTESS and the Ifpack2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include <Teuchos_ConfigDefs.hpp>
#include <Ifpack2_ConfigDefs.hpp>
#include <Teuchos_UnitTestHarness.hpp>
#include <Ifpack2_Version.hpp>
#include <iostream>

#include <Ifpack2_UnitTestHelpers.hpp>
#ifdef HAVE_MPI
#include <Teuchos_DefaultMpiComm.hpp>
#else
#include <Teuchos_DefaultSerialComm.hpp>
#endif
#include <Teuchos_RefCountPtr.hpp>
#include <Teuchos_FancyOStream.hpp>

#include <Tpetra_Map.hpp>
#include <Tpetra_CrsMatrix.hpp>
#include <Ifpack2_DiagonalFilter.hpp>
#include <Ifpack2_LocalFilter.hpp>
#include <Ifpack2_DropFilter.hpp>
#include <Ifpack2_SingletonFilter.hpp>
#include <Ifpack2_SparsityFilter.hpp>
#include <Ifpack2_ReorderFilter.hpp>

using Tpetra::global_size_t;
typedef tif_utest::Node Node;

using Teuchos::ArrayRCP;
using Teuchos::RCP;
using Teuchos::rcp;
using Teuchos::rcpFromRef;

//this macro declares the unit-test-class:
TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL(Ifpack2Filtering, Test0, Scalar, LocalOrdinal, GlobalOrdinal)
{
  std::string version = Ifpack2::Version();
  out << "Ifpack2::Version(): " << version << std::endl;
  typedef Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> CRS;
  typedef Tpetra::RowMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> ROW;

  // Useful matrices and such (tridiagonal test)
  global_size_t num_rows_per_proc = 5;
  RCP<const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> > rowmap =
    tif_utest::create_tpetra_map<LocalOrdinal,GlobalOrdinal,Node>(num_rows_per_proc);
  RCP<const CRS> Matrix =
    tif_utest::create_test_matrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>(rowmap);
  Tpetra::Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node> x(rowmap), y(rowmap), z(rowmap), b(rowmap);


  Matrix->describe(out,Teuchos::VERB_EXTREME);

  // Fill x for all time
  Teuchos::ScalarTraits<double>::seedrandom(24601);
  x.randomize();

  // ====================================== //
  // create a new matrix, diagonally filtered //
  // ====================================== //
  Scalar alpha = 100.0;
  Scalar beta = 1.0;
  Ifpack2::DiagonalFilter<ROW> DiagA(Matrix,alpha,beta);

  // Apply w/ Filter
  DiagA.apply(x,y);

  // Apply manually
  Matrix->getLocalDiagCopy(b);
  Matrix->apply(x,z);
  z.update(alpha,x,beta-1,b,1.0);

  // Diff
  TEST_COMPARE_FLOATING_ARRAYS(y.get1dView(), z.get1dView(), 1e4*Teuchos::ScalarTraits<Scalar>::eps());

  // ====================================== //
  // create a new matrix, locally filtered  //
  // ====================================== //
  Ifpack2::LocalFilter<ROW> LocalA(Matrix);
  RCP<const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> > localrowmap =
    LocalA.getRowMap ();
  Tpetra::Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node> lx(rowmap),
    ly(rowmap), lz(rowmap), la(rowmap),lb(rowmap);
  lx.randomize ();

  // Apply w/ filter
  LocalA.apply(lx,ly);

  // Apply w/ GetRow
  size_t max_nz_per_row=LocalA.getLocalMaxNumRowEntries();
  using lids_type = typename Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>::nonconst_local_inds_host_view_type;
  using vals_type = typename Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>::nonconst_values_host_view_type;
  lids_type Indices("Indices",max_nz_per_row);
  vals_type Values("Values",max_nz_per_row);
  {// Host view needs to be scoped
    Teuchos::ArrayRCP<const Scalar> xview=lx.get1dView();

    for(LocalOrdinal i=0; i < (LocalOrdinal)num_rows_per_proc; i++){
      size_t NumEntries;
      LocalA.getLocalRowCopy(i,Indices,Values,NumEntries);
      Scalar sum=0;
      for(LocalOrdinal j=0; (size_t) j < NumEntries; j++){
        sum+=Values[j] * xview[Indices[j]];
      }
      lz.replaceLocalValue(i,sum);
    }
  }

  // Diff
  TEST_COMPARE_FLOATING_ARRAYS(ly.get1dView(), lz.get1dView(), 1e4*Teuchos::ScalarTraits<Scalar>::eps());

  // ====================================== //
  // create a new matrix, dropping by value //
  // ====================================== //
  // drop all elements below 1.5. The matrix then becomes diagonal
  Ifpack2::DropFilter<ROW> DropA(RCP<ROW >(&LocalA,false),1.5);

  // Apply w/ filter
  DropA.apply(lx,ly);

  // Apply via diagonal
  LocalA.getLocalDiagCopy(la);
  lz.elementWiseMultiply(1.0,lx,la,0.0);

  // Diff
  TEST_COMPARE_FLOATING_ARRAYS(ly.get1dView(), lz.get1dView(), 1e4*Teuchos::ScalarTraits<Scalar>::eps());

  // ========================================= //
  // create a new matrix, dropping by sparsity //
  // ========================================= //
  //
  // Mantain 3 off-diagonal elements & 3 bandwidth, which is pretty boring since this is a noop.
  Ifpack2::SparsityFilter<ROW> SparsityA(RCP<ROW >(&LocalA,false),3,3);

  // Apply w/ filter
  SparsityA.apply(lx,ly);

  // Apply via local matrix
  LocalA.apply(lx,lz);

  // Diff
  TEST_COMPARE_FLOATING_ARRAYS(ly.get1dView(), lz.get1dView(), 1e4*Teuchos::ScalarTraits<Scalar>::eps());


  // ======================================== //
  // create new matrices, dropping singletons //
  // ======================================== //

  // This matrix should be the same after the singleton filter since it doesn't have singletons.
  Ifpack2::SingletonFilter<ROW> SingletonA(RCP<ROW >(&LocalA,false));

  // Apply w/ filter
  SingletonA.apply(lx,ly);

  // Apply via local matrix
  LocalA.apply(lx,lz);

  // Diff
  TEST_COMPARE_FLOATING_ARRAYS(ly.get1dView(), lz.get1dView(), 1e4*Teuchos::ScalarTraits<Scalar>::eps());


  // ======================================== //
  // create new matrices, with reordering
  // ======================================== //

  // Fill the permutation AND its inverse with a local reversal.
  // Zoltan2 would normally do this, but for unit testing we do it
  // without Zoltan2.

  // Zoltan2::LocalOrderingSolution<LocalOrdinal>
  //   Ordering ((size_t)num_rows_per_proc);
  // ArrayRCP<LocalOrdinal> l_perm = Ordering.getPermutationRCP ();
  // ArrayRCP<LocalOrdinal> l_invperm = Ordering.getPermutationRCP (true);

  ArrayRCP<LocalOrdinal> l_perm (num_rows_per_proc);
  ArrayRCP<LocalOrdinal> l_invperm (num_rows_per_proc);

  for (LocalOrdinal i = 0; i < static_cast<LocalOrdinal> (num_rows_per_proc); ++i) {
    l_perm[i] = static_cast<LocalOrdinal> (num_rows_per_proc - i - 1);
    l_invperm[i] = static_cast<LocalOrdinal> (num_rows_per_proc - i - 1);
  }

  // Now, build a reordering and a reverse reordering
  //
  // ReorderFilter has the opposite naming convention of Zoltan2, in
  // terms of which permutation is the inverse.
  Ifpack2::ReorderFilter<ROW> Reorder1 (rcpFromRef (LocalA), l_invperm, l_perm);
  Ifpack2::ReorderFilter<ROW> Reorder2 (rcpFromRef (Reorder1), l_invperm, l_perm);

  // Apply w/ double-reversed reordering
  Reorder2.apply(lx,ly);

  // Apply via local matrix
  LocalA.apply(lx,lz);

  // Diff
  TEST_COMPARE_FLOATING_ARRAYS(ly.get1dView(), lz.get1dView(), 1e4*Teuchos::ScalarTraits<Scalar>::eps());
}


#define UNIT_TEST_GROUP_SC_LO_GO(Scalar,LocalOrdinal,GlobalOrdinal) \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( Ifpack2Filtering, Test0, Scalar, LocalOrdinal,GlobalOrdinal)

#include "Ifpack2_ETIHelperMacros.h"

IFPACK2_ETI_MANGLING_TYPEDEFS()

// Test all enabled combinations of Scalar (SC), LocalOrdinal (LO),
// and GlobalOrdinal (GO) types, where Scalar is real.

IFPACK2_INSTANTIATE_SLG_REAL( UNIT_TEST_GROUP_SC_LO_GO )
