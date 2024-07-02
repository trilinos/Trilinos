// @HEADER
// *****************************************************************************
//       Ifpack2: Templated Object-Oriented Algebraic Preconditioner Package
//
// Copyright 2009 NTESS and the Ifpack2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER


/// \file Ifpack2_UnitTestAmesos2solver.cpp
/// \brief Ifpack2 unit test for the Amesos2 wrapper.
///
/// This test only builds nontrivially and executes if Amesos2 is enabled.

#include <Teuchos_ConfigDefs.hpp>
#include <Ifpack2_ConfigDefs.hpp>
#include <Teuchos_UnitTestHarness.hpp>

//#if defined(HAVE_IFPACK2_AMESOS2) && defined(HAVE_AMESOS2_SUPERLU)
#if defined(HAVE_IFPACK2_AMESOS2)

#include <Amesos2_config.h>
#include <Ifpack2_Details_Amesos2Wrapper.hpp>
#include <iostream>

#include <Ifpack2_UnitTestHelpers.hpp>

namespace {

using Teuchos::RCP;
using std::endl;
typedef Tpetra::global_size_t GST;
typedef tif_utest::Node Node;

//this macro declares the unit-test-class:
TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL(Ifpack2Amesos2Wrapper, Test0, Scalar, LO, GO)
{
//we are now in a class method declared by the above macro, and
//that method has these input arguments:
//Teuchos::FancyOStream& out, bool& success

  using Kokkos::ArithTraits;
  typedef Tpetra::Map<LO,GO,Node> map_type;
  typedef Tpetra::CrsMatrix<Scalar,LO,GO,Node> crs_matrix_type;
  typedef Tpetra::RowMatrix<Scalar,LO,GO,Node> row_matrix_type;
  #ifdef HAVE_AMESOS2_SUPERLU
  typedef Tpetra::MultiVector<Scalar,LO,GO,Node> mv_type;
  typedef typename mv_type::impl_scalar_type val_type;
  typedef typename Kokkos::ArithTraits<val_type>::mag_type mag_type;
  typedef typename map_type::device_type device_type;
  //const mag_type oneMag = ArithTraits<mag_type>::one ();
  #endif

  out << "Ifpack2 Amesos2 wrapper: Test0" << endl;
  Teuchos::OSTab tab1 (out);

  GST num_rows_per_proc = 5;
  RCP<const map_type> rowmap =
    tif_utest::create_tpetra_map<LO,GO,Node> (num_rows_per_proc);

  RCP<const crs_matrix_type> crsmatrix = tif_utest::create_test_matrix<Scalar,LO,GO,Node> (rowmap);

  out << "Create Ifpack2 Amesos2 wrapper instance" << endl;
  Ifpack2::Details::Amesos2Wrapper<row_matrix_type> prec (crsmatrix);

  out << "Set solver parameters" << endl;
  Teuchos::ParameterList params;
  params.set("Amesos2 solver name","superlu");
  Teuchos::ParameterList &sublist = params.sublist("Amesos2");
  (sublist.sublist("SuperLU")).set("ILU_Flag",false); //create SuperLU sublist to get rid of unused variable warnings
  TEST_NOTHROW(prec.setParameters(params));

  //trivial tests to insist that the preconditioner's domain/range maps are
  //identically those of the matrix:
  const map_type* mtx_dom_map_ptr = &*crsmatrix->getDomainMap();
  const map_type* mtx_rng_map_ptr = &*crsmatrix->getRangeMap();

  const map_type* prec_dom_map_ptr = &*prec.getDomainMap();
  const map_type* prec_rng_map_ptr = &*prec.getRangeMap();

  TEST_EQUALITY( prec_dom_map_ptr, mtx_dom_map_ptr );
  TEST_EQUALITY( prec_rng_map_ptr, mtx_rng_map_ptr );

# if !defined(HAVE_AMESOS2_SUPERLU)
  out << "NOT TESTING SUPERLU!!!  SUPERLU is NOT enabled!" << endl;
  TEST_THROW(prec.initialize(),std::invalid_argument);
# else
  out << "SuperLU IS enabled." << endl;

  out << "Test prec.initialize()" << endl;
  TEST_NOTHROW(prec.initialize());
  out << "Test prec.compute()" << endl;
  prec.compute();

  mv_type x (rowmap, 2);
  mv_type y (rowmap, 2);
  x.putScalar (Teuchos::ScalarTraits<Scalar>::one ());
  out << "Test prec.apply(x,y)" << endl;
  prec.apply (x, y);

  out << "y should be full of 0.5's now" << endl;
  const Scalar one = Teuchos::ScalarTraits<Scalar>::one ();
  const Scalar two = one + one;
  const Scalar oneHalf = one / two;

  mv_type diff (y.getMap (), y.getNumVectors ());
  diff.putScalar (oneHalf);
  diff.update (one, y, -one);
  Kokkos::View<mag_type*, device_type> norms ("norms", diff.getNumVectors ());
  diff.normInf (norms);
  auto norms_h = Kokkos::create_mirror_view (norms);
  Kokkos::deep_copy (norms_h, norms);

  {
    typedef Teuchos::ScalarTraits<Scalar> STS;
    const mag_type bound = ArithTraits<mag_type>::zero ();

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
# endif
}

TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL(Ifpack2Amesos2Wrapper, Test1, Scalar, LocalOrdinal, GlobalOrdinal)
{
//we are now in a class method declared by the above macro, and
//that method has these input arguments:
//Teuchos::FancyOStream& out, bool& success
  using Teuchos::RCP;
  typedef Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> map_type;
  typedef Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> crs_matrix_type;
  typedef Tpetra::RowMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> row_matrix_type;

  out << "Ifpack2 Amesos2 wrapper: Test1" << endl;
  Teuchos::OSTab tab1 (out);

  GST num_rows_per_proc = 5;

  RCP<const map_type> rowmap = tif_utest::create_tpetra_map<LocalOrdinal,GlobalOrdinal,Node> (num_rows_per_proc);

  RCP<const crs_matrix_type> crsmatrix = tif_utest::create_test_matrix3<Scalar,LocalOrdinal,GlobalOrdinal,Node>(rowmap);

  Ifpack2::Details::Amesos2Wrapper<row_matrix_type> prec (crsmatrix);

  Teuchos::ParameterList params;

  params.set("Amesos2 solver name","superlu");
  Teuchos::ParameterList &sublist = params.sublist("Amesos2");
  (sublist.sublist("SuperLU")).set("ILU_Flag",false); //create SuperLU sublist to get rid of unused variable warnings
  TEST_NOTHROW(prec.setParameters(params));

# if !defined(HAVE_AMESOS2_SUPERLU)
  TEST_THROW(prec.initialize(),std::invalid_argument);
# else
  TEST_NOTHROW(prec.initialize());
  prec.compute();

  typedef Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> mv_type;
  mv_type x(rowmap,2), y(rowmap,2);
  x.putScalar (Teuchos::ScalarTraits<Scalar>::one ());

  crsmatrix->apply (x, y);
  prec.apply (y, x);

  Teuchos::ArrayRCP<const Scalar> xview = x.get1dView();

  //x should be full of 1's now.

  Teuchos::ArrayRCP<Scalar> ones(num_rows_per_proc*2, Teuchos::ScalarTraits<Scalar>::one ());

  TEST_COMPARE_FLOATING_ARRAYS(xview, ones(), 2*Teuchos::ScalarTraits<Scalar>::eps());
# endif
}

#define UNIT_TEST_GROUP_SCALAR_ORDINAL(Scalar,LocalOrdinal,GlobalOrdinal) \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( Ifpack2Amesos2Wrapper, Test0, Scalar, LocalOrdinal, GlobalOrdinal) \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( Ifpack2Amesos2Wrapper, Test1, Scalar, LocalOrdinal, GlobalOrdinal)

// FIXME (mfh 11 Apr 2018) We should test this at least for all
// enabled real Scalar types, if not complex Scalar types as well.

typedef Tpetra::MultiVector<>::scalar_type default_scalar_type;
typedef Tpetra::MultiVector<>::local_ordinal_type default_local_ordinal_type;
typedef Tpetra::MultiVector<>::global_ordinal_type default_global_ordinal_type;

UNIT_TEST_GROUP_SCALAR_ORDINAL( default_scalar_type, default_local_ordinal_type, default_global_ordinal_type )

} // namespace (anonymous)

#endif // HAVE_IFPACK2_AMESOS2
