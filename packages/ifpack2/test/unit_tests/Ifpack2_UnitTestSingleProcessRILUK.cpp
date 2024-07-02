// @HEADER
// *****************************************************************************
//       Ifpack2: Templated Object-Oriented Algebraic Preconditioner Package
//
// Copyright 2009 NTESS and the Ifpack2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER


/*! \file Ifpack2_UnitTestSingleProcessRILUK.cpp

\brief Ifpack2 single-process unit tests for the RILUK template.
*/

#include "Teuchos_ConfigDefs.hpp"
#include "Ifpack2_ConfigDefs.hpp"
#include "Teuchos_UnitTestHarness.hpp"
#include <iostream>

#include "Tpetra_Core.hpp"
#include "Tpetra_MatrixIO.hpp"
#include "MatrixMarket_Tpetra.hpp"
#include "TpetraExt_MatrixMatrix.hpp"

#include "Ifpack2_UnitTestHelpers.hpp"
#include "Ifpack2_RILUK.hpp"

namespace { // (anonymous)

using Tpetra::global_size_t;
typedef tif_utest::Node Node;

struct IlukImplTypeDetails {
  enum Enum { Serial, KSPILUK };
};

template<class MatrixType, class VectorType>
void remove_diags_and_scale(const MatrixType& L, const MatrixType& U,
                            Teuchos::RCP<MatrixType>& Ln, Teuchos::RCP<MatrixType>& Un, Teuchos::RCP<VectorType>& Dn) {

  typedef typename MatrixType::local_matrix_device_type local_matrix_type;
  typedef typename std::remove_const<typename local_matrix_type::size_type>::type    size_type;
  typedef typename std::remove_const<typename local_matrix_type::ordinal_type>::type ordinal_type;
  typedef typename std::remove_const<typename local_matrix_type::value_type>::type   value_type;
  typedef typename local_matrix_type::device_type device_type;
  typedef typename device_type::execution_space   execution_space;
  
  typedef typename Kokkos::View<size_type*, Kokkos::LayoutLeft, device_type> rowmap_type;
  typedef typename Kokkos::View<ordinal_type*, Kokkos::LayoutLeft, device_type> entries_type;
  typedef typename Kokkos::View<value_type*, Kokkos::LayoutRight, device_type> values_type;

  typedef Kokkos::TeamPolicy<execution_space> team_policy;
  typedef typename Kokkos::TeamPolicy<execution_space>::member_type member_type;

  auto lclL = L.getLocalMatrixDevice();
  auto L_rowmap  = lclL.graph.row_map;
  auto L_entries = lclL.graph.entries;
  auto L_values  = lclL.values;
  auto lclU = U.getLocalMatrixDevice();
  auto U_rowmap  = lclU.graph.row_map;
  auto U_entries = lclU.graph.entries;
  auto U_values  = lclU.values;

  rowmap_type  Ln_rowmap ("Ln_rowmap",  L_rowmap.extent(0));
  entries_type Ln_entries("Ln_entries", L_entries.extent(0) - (L_rowmap.extent(0) - 1));
  values_type  Ln_values ("Ln_values",  L_values.extent(0)  - (L_rowmap.extent(0) - 1));
  rowmap_type  Un_rowmap ("Un_rowmap",  U_rowmap.extent(0));
  entries_type Un_entries("Un_entries", U_entries.extent(0) - (U_rowmap.extent(0) - 1));
  values_type  Un_values ("Un_values",  U_values.extent(0)  - (U_rowmap.extent(0) - 1));

  values_type  Dn_values ("Dn_values",  U_rowmap.extent(0) - 1);

  Kokkos::parallel_for( Kokkos::RangePolicy<execution_space>(0, U_rowmap.extent(0)-1), KOKKOS_LAMBDA(const int& i) {
    Ln_rowmap(i+1) = L_rowmap(i+1) - (i+1);
    Un_rowmap(i+1) = U_rowmap(i+1) - (i+1);
    Dn_values(i)   = 1.0/U_values(U_rowmap(i));
  });

  const team_policy policy( U_rowmap.extent(0)-1, Kokkos::AUTO );

  Kokkos::parallel_for( policy, KOKKOS_LAMBDA(const member_type& teamMember) {
    const int rowid = teamMember.league_rank();

    auto Lentries_src = subview(L_entries,  Kokkos::make_pair(L_rowmap(rowid),  L_rowmap(rowid+1) - 1));
    auto Lvalues_src  = subview(L_values,   Kokkos::make_pair(L_rowmap(rowid),  L_rowmap(rowid+1) - 1));
    auto Lentries_dst = subview(Ln_entries, Kokkos::make_pair(Ln_rowmap(rowid), Ln_rowmap(rowid+1)));
    auto Lvalues_dst  = subview(Ln_values,  Kokkos::make_pair(Ln_rowmap(rowid), Ln_rowmap(rowid+1)));

    auto Uentries_src = subview(U_entries,  Kokkos::make_pair(U_rowmap(rowid)+1, U_rowmap(rowid+1)));
    auto Uvalues_src  = subview(U_values,   Kokkos::make_pair(U_rowmap(rowid)+1, U_rowmap(rowid+1)));
    auto Uentries_dst = subview(Un_entries, Kokkos::make_pair(Un_rowmap(rowid),  Un_rowmap(rowid+1)));
    auto Uvalues_dst  = subview(Un_values,  Kokkos::make_pair(Un_rowmap(rowid),  Un_rowmap(rowid+1)));
	
    Kokkos::Experimental::local_deep_copy(teamMember, Lentries_dst, Lentries_src);
    Kokkos::Experimental::local_deep_copy(teamMember, Uentries_dst, Uentries_src);
    Kokkos::Experimental::local_deep_copy(teamMember, Lvalues_dst, Lvalues_src);
    Kokkos::Experimental::local_deep_copy(teamMember, Uvalues_dst, Uvalues_src);

    teamMember.team_barrier();

    Kokkos::parallel_for(Kokkos::TeamThreadRange(teamMember, Uvalues_dst.extent(0)), [&](const int& i) {
      Uvalues_dst(i) = Uvalues_dst(i)*Dn_values(rowid); 
    });

    teamMember.team_barrier();
  });
  
  Ln = Teuchos::rcp (new MatrixType (L.getRowMap(), L.getColMap(), 
                                     Ln_rowmap, Ln_entries, Ln_values));
  Un = Teuchos::rcp (new MatrixType (U.getRowMap(), U.getColMap(), 
                                     Un_rowmap, Un_entries, Un_values));
  auto Dn_view = Dn->getLocalViewDevice(Tpetra::Access::OverwriteAll);
  Kokkos::deep_copy(subview(Dn_view,Kokkos::ALL(), 0),Dn_values);

  Ln->fillComplete();
  Un->fillComplete();
}

template<typename Scalar, typename LO, typename GO>
void Ifpack2RILUKSingleProcess_test0 (bool& success, Teuchos::FancyOStream& out, const IlukImplTypeDetails::Enum ilukimplType) {
  using Teuchos::RCP;
  using std::endl;

  if (ilukimplType == IlukImplTypeDetails::Serial)
    out << "Ifpack2::RILUK: Test0 -- Serial" << endl;
  else
    out << "Ifpack2::RILUK: Test0 -- Kokkos Kernels SPILUK" << endl;
  
  Teuchos::OSTab tab0 (out);

  global_size_t num_rows_per_proc = 5;
  const RCP<const Tpetra::Map<LO,GO,Node> > rowmap =
    tif_utest::create_tpetra_map<LO,GO,Node>(num_rows_per_proc);

  if (rowmap->getComm()->getSize() > 1) {
    out << endl << "This test may only be run in serial." << endl;
    return;
  }

  RCP<const Tpetra::CrsMatrix<Scalar,LO,GO,Node> > crsmatrix =
    tif_utest::create_test_matrix<Scalar,LO,GO,Node> (rowmap);

  //----------------Default trisolver----------------//
  {
    Ifpack2::RILUK<Tpetra::RowMatrix<Scalar,LO,GO,Node> > prec (crsmatrix);
    
    Teuchos::ParameterList params;
    int fill_level = 1;
    params.set("fact: iluk level-of-fill", fill_level);
    params.set("fact: iluk level-of-overlap", 0);

    if (ilukimplType == IlukImplTypeDetails::KSPILUK)
      params.set("fact: type", "KSPILUK");
	
    TEST_NOTHROW(prec.setParameters(params));
    
    TEST_EQUALITY( prec.getLevelOfFill(), fill_level);
    
    prec.initialize();
    //trivial tests to insist that the preconditioner's domain/range maps are
    //the same as those of the matrix:
    const Tpetra::Map<LO,GO,Node>& mtx_dom_map = *crsmatrix->getDomainMap();
    const Tpetra::Map<LO,GO,Node>& mtx_rng_map = *crsmatrix->getRangeMap();
    
    const Tpetra::Map<LO,GO,Node>& prec_dom_map = *prec.getDomainMap();
    const Tpetra::Map<LO,GO,Node>& prec_rng_map = *prec.getRangeMap();
    
    TEST_ASSERT( prec_dom_map.isSameAs(mtx_dom_map) );
    TEST_ASSERT( prec_rng_map.isSameAs(mtx_rng_map) );
    
    prec.compute();
    
    Tpetra::MultiVector<Scalar,LO,GO,Node> x(rowmap,2), y(rowmap,2);
    x.putScalar (Teuchos::ScalarTraits<Scalar>::one ());
    
    prec.apply(x, y);
    
    Teuchos::ArrayRCP<const Scalar> yview = y.get1dView();
    
    //y should be full of 0.5's now.
    
    Teuchos::ArrayRCP<Scalar> halfs(num_rows_per_proc*2, 0.5);
    
    TEST_COMPARE_FLOATING_ARRAYS(yview, halfs(), Teuchos::ScalarTraits<Scalar>::eps());
  }
  //----------------Kokkos Kernels SPTRSV----------------//
  {
    Ifpack2::RILUK<Tpetra::RowMatrix<Scalar,LO,GO,Node> > prec (crsmatrix);
    
    Teuchos::ParameterList params;
    int fill_level = 1;
    params.set("fact: iluk level-of-fill", fill_level);
    params.set("fact: iluk level-of-overlap", 0);

    if (ilukimplType == IlukImplTypeDetails::KSPILUK)
      params.set("fact: type", "KSPILUK");

    params.set("trisolver: type", "KSPTRSV");
    
    TEST_NOTHROW(prec.setParameters(params));
    
    TEST_EQUALITY( prec.getLevelOfFill(), fill_level);
    
    prec.initialize();
    //trivial tests to insist that the preconditioner's domain/range maps are
    //the same as those of the matrix:
    const Tpetra::Map<LO,GO,Node>& mtx_dom_map = *crsmatrix->getDomainMap();
    const Tpetra::Map<LO,GO,Node>& mtx_rng_map = *crsmatrix->getRangeMap();
    
    const Tpetra::Map<LO,GO,Node>& prec_dom_map = *prec.getDomainMap();
    const Tpetra::Map<LO,GO,Node>& prec_rng_map = *prec.getRangeMap();
    
    TEST_ASSERT( prec_dom_map.isSameAs(mtx_dom_map) );
    TEST_ASSERT( prec_rng_map.isSameAs(mtx_rng_map) );
    
    prec.compute();
    
    Tpetra::MultiVector<Scalar,LO,GO,Node> x(rowmap,2), y(rowmap,2);
    x.putScalar (Teuchos::ScalarTraits<Scalar>::one ());
    
    prec.apply(x, y);
    
    Teuchos::ArrayRCP<const Scalar> yview = y.get1dView();
    
    //y should be full of 0.5's now.
    
    Teuchos::ArrayRCP<Scalar> halfs(num_rows_per_proc*2, 0.5);
    
    TEST_COMPARE_FLOATING_ARRAYS(yview, halfs(), Teuchos::ScalarTraits<Scalar>::eps());
  }
}

template<typename Scalar, typename LO, typename GO>
void Ifpack2RILUKSingleProcess_test1 (bool& success, Teuchos::FancyOStream& out, const IlukImplTypeDetails::Enum ilukimplType) {
  using Kokkos::ArithTraits;
  using Teuchos::RCP;
  using std::endl;
  typedef Tpetra::Map<LO, GO, Node> map_type;
  typedef Tpetra::CrsMatrix<Scalar,LO,GO,Node> crs_matrix_type;
  typedef Tpetra::MultiVector<Scalar,LO,GO,Node> MV;
  typedef Tpetra::RowMatrix<Scalar,LO,GO,Node> row_matrix_type;
  typedef Teuchos::ScalarTraits<Scalar> STS;
  typedef typename MV::impl_scalar_type val_type;
  typedef typename Kokkos::ArithTraits<val_type>::mag_type mag_type;
  typedef typename map_type::device_type device_type;
  const mag_type oneMag = ArithTraits<mag_type>::one ();
  const mag_type twoMag = oneMag + oneMag;

  if (ilukimplType == IlukImplTypeDetails::Serial)
    out << "Ifpack2::RILUK: Test1 -- Serial" << endl;
  else
    out << "Ifpack2::RILUK: Test1 -- Kokkos Kernels SPILUK" << endl;

  Teuchos::OSTab tab1 (out);

  const global_size_t num_rows_per_proc = 5;
  RCP<const map_type> rowmap =
    tif_utest::create_tpetra_map<LO, GO, Node> (num_rows_per_proc);

  // Matrix
  // [ 2 .1  0  0  0]
  // [.1  2  0  0  0]
  // [ 0 .1  2 .1  0]
  // [ 0  0 .1  2 .1]
  // [ 0  0  0 .1  2]

  // Matlab's Factors
  // L
  // Diagonal = 1 (implied)
  // Subdiagonal (approx) = .05, .0501, .0501 .0501

  // U
  // Diagonal (approx)      = 2 1.995 1.995 1.995 1.995 
  // Superdiagonal (approx) = .1 .1 .1 .1


  if (rowmap->getComm ()->getSize () > 1) {
    out << "This test may only be run in serial "
      "or with a single MPI process." << endl;
    return;
  }

  out << "Creating matrix" << endl;
  RCP<const crs_matrix_type> crsmatrix =
    tif_utest::create_test_matrix2<Scalar,LO,GO,Node>(rowmap);

  {//CMS
    out<<"***** A *****"<<std::endl;
    crsmatrix->describe(out,Teuchos::VERB_EXTREME);
  }

  //----------------Default trisolver----------------//
  {
    out << "Creating preconditioner" << endl;
    Ifpack2::RILUK<row_matrix_type> prec (crsmatrix);
    
    out << "Setting preconditioner's parameters" << endl;
    Teuchos::ParameterList params;
    params.set ("fact: iluk level-of-fill", 1);
    params.set ("fact: iluk level-of-overlap", 0);
    if (ilukimplType == IlukImplTypeDetails::KSPILUK)
      params.set("fact: type", "KSPILUK");
    TEST_NOTHROW(prec.setParameters(params));
    
    out << "Calling initialize() and compute()" << endl;
    prec.initialize();
    prec.compute();
   
  {//CMS
    out<<"***** Test L *****"<<std::endl;
    prec.getL().describe(out,Teuchos::VERB_EXTREME);
    out<<"***** Test U *****"<<std::endl;
    prec.getU().describe(out,Teuchos::VERB_EXTREME);
    out<<"***** Test D *****"<<std::endl;
    prec.getD().describe(out,Teuchos::VERB_EXTREME);
  }

 
    out << "Creating test problem" << endl;
    MV x (rowmap, 2);
    MV y (rowmap, 2);
    x.putScalar (STS::one ());
    
    out << "Calling crsmatrix->apply(x, y)" << endl;
    crsmatrix->apply(x,y);
    
    out << "Calling prec.apply(y, x)" << endl;
    //apply the preconditioner to y, putting ~A^-1*y in x
    //(this should set x back to 1's)
    prec.apply(y, x);
    
    //x should be full of 1's now.
    out << "Checking result" << endl;
    
    // Useful things for comparing results.
    Kokkos::View<mag_type*, device_type> norms ("norms", x.getNumVectors ());
    auto norms_h = Kokkos::create_mirror_view (norms);
    MV diff (x.getMap (), x.getNumVectors ());
    
    {
      // FIXME (mfh 04 Oct 2016) This is the bound that I found here
      // when I fixed this test.  Not sure if it's sensible.
      const mag_type bound = twoMag * ArithTraits<val_type>::eps ();
    
      diff.putScalar (STS::one ());
      diff.update (STS::one (), x, -STS::one ());
      diff.normInf (norms);
      Kokkos::deep_copy (norms_h, norms);
      for (LO j = 0; j < static_cast<LO> (norms_h.extent (0)); ++j) {
        const mag_type absVal = ArithTraits<mag_type>::abs (norms_h(j));
        TEST_ASSERT( absVal <= bound );
        if (absVal > bound) {
          out << "\\| x(:," << j << ") - 1 \\|_{\\infty} = "
              << absVal << " > " << bound
              << "; the norm equals " << norms_h(j)
              << endl;
        }
      }
    }
    
    auto test_alpha_beta = [&] (const Scalar& alpha, const Scalar& beta,
                                const Teuchos::ETransp& mode) {
      out << "Testing apply() for alpha = " << alpha
          << " and beta = " << beta << endl;
      const Scalar x_magic_number = -0.42;
      x.putScalar (x_magic_number);
      crsmatrix->apply (x, y, mode);
      MV z = Tpetra::createCopy (x);
      const Scalar z_magic_number = 2.1;
      z.putScalar (z_magic_number);
      // z = beta z + alpha inv(A) y
      //   = beta z + alpha x
      //   = (beta z_magic_number + alpha x_magic_number) 1s
      prec.apply(y, z, mode, alpha, beta);
    
      MV z_true = Tpetra::createCopy (x);
      const Scalar z_true_scalar = beta*z_magic_number + alpha*x_magic_number;
      z_true.putScalar(z_true_scalar);
    
      {
        // FIXME (mfh 04 Oct 2016) This is the bound that I found here
        // when I fixed this test.  Not sure if it's sensible.
        const mag_type bound = 10.0 * ArithTraits<val_type>::eps ();
    
        diff.putScalar (z_true_scalar);
        diff.update (STS::one (), z, -STS::one ());
        diff.normInf (norms);
        Kokkos::deep_copy (norms_h, norms);
    
        for (LO j = 0; j < static_cast<LO> (norms_h.extent (0)); ++j) {
          const mag_type absVal = ArithTraits<mag_type>::abs (norms_h(j));
          TEST_ASSERT( absVal <= bound );
          if (absVal > bound) {
            out << "\\| x(:," << j << ") - 1 \\|_{\\infty} = "
                << absVal << " > " << bound
                << "; the norm equals " << norms_h(j)
                << endl;
          }
        }
      }
    };
    
    for (const auto mode : {Teuchos::NO_TRANS, Teuchos::TRANS}) {
      test_alpha_beta(0.0, 0.0, mode);
      test_alpha_beta(2.0, 0.0, mode);
      test_alpha_beta(0.0, -1.5, mode);
      test_alpha_beta(-0.42, 4.2, mode);
    }
  }

  return;//CMS
  //----------------Kokkos Kernels SPTRSV----------------//
  {
    out << "Creating preconditioner" << endl;
    Ifpack2::RILUK<row_matrix_type> prec (crsmatrix);
    
    out << "Setting preconditioner's parameters" << endl;
    Teuchos::ParameterList params;
    params.set ("fact: iluk level-of-fill", 1);
    params.set ("fact: iluk level-of-overlap", 0);
    if (ilukimplType == IlukImplTypeDetails::KSPILUK)
      params.set("fact: type", "KSPILUK");
    params.set("trisolver: type", "KSPTRSV");
    TEST_NOTHROW(prec.setParameters(params));
    
    out << "Calling initialize() and compute()" << endl;
    prec.initialize();
    prec.compute();
    
    out << "Creating test problem" << endl;
    MV x (rowmap, 2);
    MV y (rowmap, 2);
    x.putScalar (STS::one ());
    
    out << "Calling crsmatrix->apply(x, y)" << endl;
    crsmatrix->apply(x,y);
    
    out << "Calling prec.apply(y, x)" << endl;
    //apply the preconditioner to y, putting ~A^-1*y in x
    //(this should set x back to 1's)
    prec.apply(y, x);
    
    //x should be full of 1's now.
    out << "Checking result" << endl;
    
    // Useful things for comparing results.
    Kokkos::View<mag_type*, device_type> norms ("norms", x.getNumVectors ());
    auto norms_h = Kokkos::create_mirror_view (norms);
    MV diff (x.getMap (), x.getNumVectors ());
    
    {
      // FIXME (mfh 04 Oct 2016) This is the bound that I found here
      // when I fixed this test.  Not sure if it's sensible.
      const mag_type bound = twoMag * ArithTraits<val_type>::eps ();
    
      diff.putScalar (STS::one ());
      diff.update (STS::one (), x, -STS::one ());
      diff.normInf (norms);
      Kokkos::deep_copy (norms_h, norms);
      for (LO j = 0; j < static_cast<LO> (norms_h.extent (0)); ++j) {
        const mag_type absVal = ArithTraits<mag_type>::abs (norms_h(j));
        TEST_ASSERT( absVal <= bound );
        if (absVal > bound) {
          out << "\\| x(:," << j << ") - 1 \\|_{\\infty} = "
              << absVal << " > " << bound
              << "; the norm equals " << norms_h(j)
              << endl;
        }
      }
    }
    
    auto test_alpha_beta = [&] (const Scalar& alpha, const Scalar& beta,
                                const Teuchos::ETransp& mode) {
      out << "Testing apply() for alpha = " << alpha
          << " and beta = " << beta << endl;
      const Scalar x_magic_number = -0.42;
      x.putScalar (x_magic_number);
      crsmatrix->apply (x, y, mode);
      MV z = Tpetra::createCopy (x);
      const Scalar z_magic_number = 2.1;
      z.putScalar (z_magic_number);
      // z = beta z + alpha inv(A) y
      //   = beta z + alpha x
      //   = (beta z_magic_number + alpha x_magic_number) 1s
      prec.apply(y, z, mode, alpha, beta);
    
      MV z_true = Tpetra::createCopy (x);
      const Scalar z_true_scalar = beta*z_magic_number + alpha*x_magic_number;
      z_true.putScalar(z_true_scalar);
    
      {
        // FIXME (mfh 04 Oct 2016) This is the bound that I found here
        // when I fixed this test.  Not sure if it's sensible.
        const mag_type bound = 10.0 * ArithTraits<val_type>::eps ();
    
        diff.putScalar (z_true_scalar);
        diff.update (STS::one (), z, -STS::one ());
        diff.normInf (norms);
        Kokkos::deep_copy (norms_h, norms);
    
        for (LO j = 0; j < static_cast<LO> (norms_h.extent (0)); ++j) {
          const mag_type absVal = ArithTraits<mag_type>::abs (norms_h(j));
          TEST_ASSERT( absVal <= bound );
          if (absVal > bound) {
            out << "\\| x(:," << j << ") - 1 \\|_{\\infty} = "
                << absVal << " > " << bound
                << "; the norm equals " << norms_h(j)
                << endl;
          }
        }
      }
    };
    
    for (const auto mode : {Teuchos::NO_TRANS, Teuchos::TRANS}) {
      test_alpha_beta(0.0, 0.0, mode);
      test_alpha_beta(2.0, 0.0, mode);
      test_alpha_beta(0.0, -1.5, mode);
      test_alpha_beta(-0.42, 4.2, mode);
    }
  }
  out << "Done with test" << endl;
}
	
TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL(Ifpack2RILUKSingleProcess, Test0, Scalar, LO, GO)
{
  Ifpack2RILUKSingleProcess_test0<Scalar, LO, GO> (success, out, IlukImplTypeDetails::Serial);
  Ifpack2RILUKSingleProcess_test0<Scalar, LO, GO> (success, out, IlukImplTypeDetails::KSPILUK);
}

TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL(Ifpack2RILUKSingleProcess, Test1, Scalar, LO, GO)
{
  Ifpack2RILUKSingleProcess_test1<Scalar, LO, GO> (success, out, IlukImplTypeDetails::Serial);
  Ifpack2RILUKSingleProcess_test1<Scalar, LO, GO> (success, out, IlukImplTypeDetails::KSPILUK);
}

TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL(Ifpack2RILUKSingleProcess, FillLevel, Scalar, LO, GO)
{
  // Test that ILU(k) computes correct factors in serial for fill levels 0 to 5.
  // This test does nothing in parallel.
  // 1) Create a banded matrix A with bandwidth equal to lof+2.
  //    Nonzero off-diagonals are subdiagonals +/-1 and +/-(lof+2).
  //    The matrix has 4's on the main diagonal, -1's on off-diagonals.
  // 2) Compute ILU(lof) of A.
  // 3) The incomplete factors should be equal to those of an exact LU decomposition without pivoting.
  //    Note that Ifpack2 stores the inverse of the diagonal separately and scales U.

  using Teuchos::RCP;
  using std::endl;
  typedef Tpetra::MultiVector<Scalar,LO,GO,Node>                multivector_type;
  typedef Tpetra::Vector<Scalar,LO,GO,Node>                     vector_type;
  typedef Tpetra::CrsMatrix<Scalar,LO,GO,Node>                  crs_matrix_type;
  typedef Tpetra::RowMatrix<Scalar,LO,GO,Node>                  row_matrix_type;
  typedef Tpetra::MatrixMarket::Reader<crs_matrix_type>         reader_type;
  typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType magnitudeType;
  typedef Tpetra::Map<LO,GO,Node>                               map_type;
  typedef Teuchos::ScalarTraits<Scalar>                         TST;

  out << "Ifpack2::RILUK: FillLevel" << endl;

  RCP<const Teuchos::Comm<int> > comm = Tpetra::getDefaultComm();

  if (comm->getSize() > 1) {
    out << endl << "This test is only meaningful in serial." << endl;
    return;
  }

  global_size_t num_rows_per_proc = 10;
  const RCP<const Tpetra::Map<LO,GO,Node> > rowmap =
    tif_utest::create_tpetra_map<LO,GO,Node>(num_rows_per_proc);

  for (GO impltype=0; impltype<2; ++impltype) {
    for (GO lof=0; lof<6; ++lof) {
      RCP<const crs_matrix_type > crsmatrix = tif_utest::create_banded_matrix<Scalar,LO,GO,Node>(rowmap,lof+2);
      //std::string aFile = "A_bw=" + Teuchos::toString(lof+2) + ".mm";
      //RCP<crs_matrix_type> crsmatrix = reader_type::readSparseFile (aFile, comm);
      //crsmatrix->describe(out,Teuchos::VERB_EXTREME);
      Ifpack2::RILUK<row_matrix_type> prec(crsmatrix);
    
      Teuchos::ParameterList params;
      params.set("fact: iluk level-of-fill", lof);
      params.set("fact: iluk level-of-overlap", 0);
      params.set("fact: iluk overalloc", 1.01); // Use non-default only to test
                                                // matrix resizing in IlukGraph;
                                                // Usually, 1.01 is too small;
                                                // default value should be used.
      if (impltype == 1)
        params.set("fact: type", "KSPILUK");
      prec.setParameters(params);
      prec.initialize();
      prec.compute();
      //extract incomplete factors
      const crs_matrix_type &iL = prec.getL();
      const crs_matrix_type &iU = prec.getU();
    
      Teuchos::RCP<crs_matrix_type> iLn;
      Teuchos::RCP<crs_matrix_type> iUn;
      Teuchos::RCP<vector_type> iDn;
  
      if (impltype == 0) {
        const vector_type &iD = prec.getD();
        iLn = rcp (new crs_matrix_type (iL));
        iUn = rcp (new crs_matrix_type (iU));
        iDn = rcp (new vector_type (iD));
      }
      else {
        iDn = Teuchos::rcp (new vector_type (rowmap));
        remove_diags_and_scale(iL, iU, iLn, iUn, iDn);
      }

      //read L,U, and D factors from file
      std::string lFile = "Lfactor_bw" + Teuchos::toString(lof+2) + ".mm";
      std::string uFile = "Ufactor_bw" + Teuchos::toString(lof+2) + ".mm";
      std::string dFile = "Dfactor_bw" + Teuchos::toString(lof+2) + ".mm";
      out << "reading " << lFile << ", " << uFile << ", " << dFile << std::endl;
      RCP<crs_matrix_type> L = reader_type::readSparseFile (lFile, comm);
      RCP<crs_matrix_type> U = reader_type::readSparseFile (uFile, comm);
      RCP<const map_type> rm = U->getRowMap();
      RCP<multivector_type> D = reader_type::readVectorFile (dFile, comm, rm);
      
      //compare factors
      out << "bandwidth = " << lof+2 << ", lof = " << lof << std::endl;
      D->update(TST::one(),*iDn,-TST::one());
      RCP<crs_matrix_type> matdiff = Tpetra::MatrixMatrix::add(1.,false,*L,-1.,false,*iLn);
      magnitudeType mag = matdiff->getFrobeniusNorm();
      out << "||L - iL||_fro = " << mag << std::endl;
      TEST_EQUALITY(mag < 1e-12, true);
      out << std::endl;
      matdiff = Tpetra::MatrixMatrix::add(1.,false,*U,-1.,false,*iUn);
      mag = matdiff->getFrobeniusNorm();
      out << "||U - iU||_fro = " << mag << std::endl;
      TEST_EQUALITY(mag < 1e-12, true);
      out << std::endl;
      Teuchos::Array<magnitudeType> norms(1);
      D->norm2(norms);
      out << "||inverse(D) - inverse(iD)||_2 = " << norms[0] << std::endl;
      TEST_EQUALITY(norms[0] < 1e-12, true);
      out << std::endl;
    } //for (GO lof=0; lof<6; ++lof)
  }
} // unit test FillLevel()

TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL(Ifpack2RILUKSingleProcess, IgnoreRowMapGIDs, Scalar, LO, GO)
{
  // Test that ILU(k) ignores ordering of GIDs in the matrix rowmap.  This test is a virtual duplicate
  // of the test "FillLevel", with the exception that the row map GIDs are permuted.
  // This test is associated with bug#6033.
  //
  // This test does nothing in parallel.
  //
  // 1) Create a banded matrix A with bandwidth equal to lof+2.
  //    Nonzero off-diagonals are subdiagonals +/-1 and +/-(lof+2).
  //    The matrix has 4's on the main diagonal, -1's on off-diagonals.
  // 2) Compute ILU(lof) of A.
  // 3) The incomplete factors should be equal to those of an exact LU decomposition without pivoting.
  //    Note that Ifpack2 stores the inverse of the diagonal separately and scales U.

  using Teuchos::RCP;
  using std::endl;
  typedef Tpetra::CrsMatrix<Scalar,LO,GO,Node> crs_matrix_type;
  typedef Tpetra::MultiVector<Scalar,LO,GO,Node> multivector_type;
  typedef Tpetra::Vector<Scalar,LO,GO,Node> vector_type;
  typedef Tpetra::RowMatrix<Scalar,LO,GO,Node> row_matrix_type;
  typedef Tpetra::MatrixMarket::Reader<crs_matrix_type> reader_type;
  typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType magnitudeType;
  typedef Tpetra::Map<LO,GO,Node>                               map_type;
  typedef Teuchos::ScalarTraits<Scalar>                         TST;
  const global_size_t INVALID = Teuchos::OrdinalTraits<global_size_t>::invalid ();

  out << "Ifpack2::RILUK: IgnoreRowMapGIDs" << endl;

  RCP<const Teuchos::Comm<int> > comm = Tpetra::getDefaultComm();

  if (comm->getSize() > 1) {
    out << endl << "This test is only meaningful in serial." << endl;
    return;
  }

  global_size_t num_rows_per_proc = 10;
  const RCP<const Tpetra::Map<LO,GO,Node> > rowMap =
    tif_utest::create_tpetra_map<LO,GO,Node>(num_rows_per_proc);

  //Create a permuted row map.  The first entry is the same as the original row map,
  //the remainder are in descending order.
  Teuchos::ArrayView<const GO> GIDs = rowMap->getLocalElementList();
  Teuchos::Array<GO> permutedGIDs(GIDs.size());
  Teuchos::Array<GO> origToPerm(GIDs.size());
  permutedGIDs[0] = GIDs[0];
  origToPerm[0] = 0;
  for (GO i=1; i<GIDs.size(); ++i) {
    permutedGIDs[i] = GIDs[GIDs.size()-i];
    origToPerm[GIDs[GIDs.size()-i]] = i;
  }
  const LO indexBase = 0;
  Teuchos::RCP<const map_type> permRowMap = Teuchos::rcp(new map_type(INVALID, permutedGIDs(), indexBase, comm));

  for (GO impltype=0; impltype<2; ++impltype) {  
    for (GO lof=0; lof<6; ++lof) {
    
      RCP<const crs_matrix_type > crsmatrix = tif_utest::create_banded_matrix<Scalar,LO,GO,Node>(rowMap,lof+2);
    
      //Copy the banded matrix into a new matrix with permuted row map GIDs.
      //This matrix will have the sparsity pattern as the original matrix.
      RCP<crs_matrix_type> permutedMatrix = Teuchos::rcp(new crs_matrix_type(permRowMap, 5));
      typename crs_matrix_type::nonconst_global_inds_host_view_type Inds("Inds",5), pInds("pInds",5);
      typename crs_matrix_type::nonconst_values_host_view_type Vals("Vals",5), pVals("pVals",5);
      size_t numEntries;
      for (global_size_t i=0; i<num_rows_per_proc; ++i) {
        crsmatrix->getGlobalRowCopy(i,Inds,Vals,numEntries);
        Kokkos::resize(pInds,numEntries);
        Kokkos::resize(pVals,numEntries);
        for (size_t j=0; j<numEntries; ++j) {
          pInds[j] = origToPerm[Inds[j]];
          pVals[j] = Vals[j];
        }
        permutedMatrix->insertGlobalValues(origToPerm[i],numEntries,pVals.data(),pInds.data());
      }
      permutedMatrix->fillComplete();
    
      Ifpack2::RILUK<row_matrix_type> prec(Teuchos::as< RCP<const crs_matrix_type> >(permutedMatrix));
    
      Teuchos::ParameterList params;
      params.set("fact: iluk level-of-fill", lof);
      params.set("fact: iluk level-of-overlap", 0);
      if (impltype == 1)
        params.set("fact: type", "KSPILUK");
    
      prec.setParameters(params);
      prec.initialize();
      prec.compute();
      //extract incomplete factors
      const crs_matrix_type &iL = prec.getL();
      const crs_matrix_type &iU = prec.getU();
    
      Teuchos::RCP<crs_matrix_type> iLn;
      Teuchos::RCP<crs_matrix_type> iUn;
      Teuchos::RCP<vector_type> iDn;
  
      if (impltype == 0) {
        const vector_type &iD = prec.getD();
        iLn = rcp (new crs_matrix_type (iL));
        iUn = rcp (new crs_matrix_type (iU));
        iDn = rcp (new vector_type (iD));
      }
      else {
        iDn = Teuchos::rcp (new vector_type (rowMap));
        remove_diags_and_scale(iL, iU, iLn, iUn, iDn);
      }

      //read L,U, and D factors from file
      std::string lFile = "Lfactor_bw" + Teuchos::toString(lof+2) + ".mm";
      std::string uFile = "Ufactor_bw" + Teuchos::toString(lof+2) + ".mm";
      std::string dFile = "Dfactor_bw" + Teuchos::toString(lof+2) + ".mm";
      out << "reading " << lFile << ", " << uFile << ", " << dFile << std::endl;
      RCP<crs_matrix_type> L = reader_type::readSparseFile (lFile, comm);
      RCP<crs_matrix_type> U = reader_type::readSparseFile (uFile, comm);
      RCP<const map_type> rm = U->getRowMap();
      
      //Compare factors.  We can't use the Frobenius norm, as it uses GIDs.
      //Instead, we use the trick of multiply by the same random vector and comparing the
      //two norm of the results.  One of the random vectors is based on the original rowmap,
      //the other on the permuted row map.  Both contain the same random entries.
      multivector_type randVec(rowMap,1);
      randVec.randomize();
      multivector_type permRandVec(permRowMap,1);
      Teuchos::ArrayRCP<const Scalar> data  = randVec.getData(0);
      Teuchos::ArrayRCP<Scalar> pdata = permRandVec.getDataNonConst(0);
      for (global_size_t i=0; i<num_rows_per_proc; ++i)
        pdata[i] = data[i];
      data = pdata = Teuchos::null;
      
      out << "bandwidth = " << lof+2 << ", lof = " << lof << std::endl;
      multivector_type permResult(permRowMap,1);
      iLn->apply(permRandVec,permResult);
      Teuchos::Array<magnitudeType> n1(1);
      permResult.norm2(n1);
      
      multivector_type result(rowMap,1);
      L->apply(randVec,result);
      Teuchos::Array<magnitudeType> n2(1);
      result.norm2(n2);
      
      out << "||L*randvec||_2 - ||iL*randvec||_2 = " << n1[0]-n2[0] << std::endl;
      TEST_EQUALITY(n1[0]-n2[0] < 1e-6, true);
      out << std::endl;

      iUn->apply(permRandVec,permResult);
      permResult.norm2(n1);
      
      U->apply(randVec,result);
      result.norm2(n2);
      
      out << "||U*randvec||_2 - ||iU*randvec||_2 = " << n1[0]-n2[0] << std::endl;
      {
        typedef typename TST::magnitudeType MT;
        // Heuristic for rounding error: machine epsilon times square
        // root of the number of floating-point operations.
        const MT tol = TST::eps () * TST::squareroot (static_cast<MT> (result.getGlobalLength ()));
        TEST_ASSERT( n1[0]-n2[0] < tol );
      }
      out << std::endl;
      
      RCP<multivector_type> D = reader_type::readVectorFile (dFile, comm, rm);
	  D->update(TST::one(),*iDn,-TST::one());
      Teuchos::Array<magnitudeType> norms(1);
      D->norm2(norms);
      out << "||inverse(D) - inverse(iD)||_2 = " << norms[0] << std::endl;
      TEST_EQUALITY(norms[0] < 1e-7, true);
      out << std::endl;

    } //for (GO lof=0; lof<6; ++lof)
  }
} //unit test IgnoreRowMapGIDs()

TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL(Ifpack2RILUKSingleProcess, TestGIDConsistency, Scalar, LO, GO)
{
  // Test that ILU(k) throws an exception if the ordering of the GIDs
  // in the row Map is not the same as the ordering of the local GIDs
  // in the column Map.  The ILU(k) setup and algorithm assumes this
  // for the moment.

  // 25April2014 JJH: The local filter appears to fix the column Map in parallel so that it's
  //                  consistently ordered with the row Map.  In otherwords, I can't get this
  //                  test to fail in parallel.  So this check is only necessary in serial.

  using Teuchos::RCP;
  using Teuchos::rcp;
  using std::endl;
  typedef Tpetra::CrsMatrix<Scalar, LO, GO, Node> crs_matrix_type;
  typedef Tpetra::RowMatrix<Scalar, LO, GO, Node> row_matrix_type;
  typedef Tpetra::Map<LO, GO, Node> map_type;
  typedef Tpetra::global_size_t GST;

  out << "Ifpack2::RILUK: TestGIDConsistency" << endl;

  const GST INVALID = Teuchos::OrdinalTraits<GST>::invalid ();
  RCP<const Teuchos::Comm<int> > comm = Tpetra::getDefaultComm ();

  if (comm->getSize () > 1) {
    out << endl << "This test only runs in serial." << endl;
    return;
  }

  // Create matrix: 5 rows per process, 3 entries per row
  const LO indexBase = 0;
  GST numRowsPerProc = 5;
  RCP<map_type> rowMap =
    rcp (new map_type (INVALID, numRowsPerProc, indexBase, comm));

  // Create a column Map with the same GIDs at the row Map, but in
  // permuted order.  The first entry is the same as the row Map, the
  // remainder are in descending order.
  Teuchos::ArrayView<const GO> rowGIDs = rowMap->getLocalElementList ();
  Teuchos::Array<GO> colElements (rowGIDs.size ());
  colElements[0] = rowGIDs[0];
  for (GO i = 1; i < rowGIDs.size (); ++i) {
    colElements[i] = rowGIDs[rowGIDs.size () - i];
  }

  RCP<const map_type> colMap =
    rcp (new map_type (INVALID, colElements (), indexBase, comm));
  RCP<crs_matrix_type> A = rcp (new crs_matrix_type (rowMap, colMap, 3));

  // Construct a nondiagonal matrix.  It's not tridiagonal because of
  // process boundaries.
  const Scalar one = Teuchos::ScalarTraits<Scalar>::one ();
  const Scalar two = one + one;
  Teuchos::Array<GO> col (3);
  Teuchos::Array<Scalar> val (3);
  size_t numLocalElts = rowMap->getLocalNumElements ();
  for (LO l_row = 0; static_cast<size_t> (l_row) < numLocalElts; ++l_row) {
    const GO g_row = rowMap->getGlobalElement (l_row);
    size_t i = 0;
    col[i] = g_row;
    val[i++] = two;
    if (l_row>0) {
      col[i] = rowMap->getGlobalElement (l_row - 1);
      val[i++] = -one;
    }
    if (static_cast<size_t> (l_row) < numLocalElts - 1) {
      col[i] = rowMap->getGlobalElement (l_row + 1);
      val[i++] = -one;
    }
    A->insertGlobalValues (g_row, col (0, i), val (0, i));
  }
  A->fillComplete ();

  RCP<const crs_matrix_type> constA = A;
  
  {
    Ifpack2::RILUK<row_matrix_type> prec (constA);
    
    Teuchos::ParameterList params;
    GO lof = 1;
    params.set ("fact: iluk level-of-fill", lof);
    params.set ("fact: iluk level-of-overlap", 0);
      
    prec.setParameters (params);
    TEST_THROW( prec.initialize (), std::runtime_error);
  }
  {
    Ifpack2::RILUK<row_matrix_type> prec (constA);
    
    Teuchos::ParameterList params;
    GO lof = 1;
    params.set ("fact: iluk level-of-fill", lof);
    params.set ("fact: iluk level-of-overlap", 0);
    params.set ("fact: type", "KSPILUK");
      
    prec.setParameters (params);
    TEST_THROW( prec.initialize (), std::runtime_error);
  }

} // unit test TestGIDConsistency()

//
// Instantiate and run unit tests
//

#define UNIT_TEST_GROUP_SC_LO_GO( SC, LO, GO ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( Ifpack2RILUKSingleProcess, Test0, SC, LO, GO ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( Ifpack2RILUKSingleProcess, Test1, SC, LO, GO ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( Ifpack2RILUKSingleProcess, FillLevel, SC, LO, GO ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( Ifpack2RILUKSingleProcess, IgnoreRowMapGIDs, SC, LO, GO ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( Ifpack2RILUKSingleProcess, TestGIDConsistency, SC, LO, GO )

// FIXME (21 Oct 2015) There was a FIXME here a while back about
// matrix-matrix add not getting instantiated for Scalar != double.
// Need to fix that to make this test work for Scalar != double.

#ifdef HAVE_TPETRA_INST_DOUBLE
#  define UNIT_TEST_GROUP_LO_GO( LO, GO ) \
     UNIT_TEST_GROUP_SC_LO_GO( double, LO, GO )
#else // NOT HAVE_TPETRA_INST_DOUBLE
#  define UNIT_TEST_GROUP_LO_GO( LO, GO )
#endif // HAVE_TPETRA_INST_DOUBLE

#include "Ifpack2_ETIHelperMacros.h"

IFPACK2_ETI_MANGLING_TYPEDEFS()

IFPACK2_INSTANTIATE_LG( UNIT_TEST_GROUP_LO_GO )

} // namespace (anonymous)

