/*
//@HEADER
// ***********************************************************************
//
//       Ifpack2: Templated Object-Oriented Algebraic Preconditioner Package
//                 Copyright (2009) Sandia Corporation
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
//@HEADER
*/


/*! \file Ifpack2_UnitTestAdditiveSchwarz.cpp

\brief Ifpack2 Unit test for AdditiveSchwarz.
*/

#include "Ifpack2_config.h"
#ifdef HAVE_IFPACK2_AMESOS2
#  include "Amesos2_config.h"
#endif // HAVE_IFPACK2_AMESOS2
#include "Teuchos_UnitTestHarness.hpp"

// Xpetra / Galeri
#ifdef HAVE_IFPACK2_XPETRA
#include "Xpetra_ConfigDefs.hpp"
#include "Xpetra_DefaultPlatform.hpp"
#include "Xpetra_Parameters.hpp"
#include "Xpetra_MapFactory.hpp"
#include "Xpetra_TpetraMap.hpp"
#include "Xpetra_CrsMatrix.hpp"
#include "Xpetra_TpetraCrsMatrix.hpp"
#include "Galeri_XpetraProblemFactory.hpp"
#include "Galeri_XpetraMatrixTypes.hpp"
#include "Galeri_XpetraParameters.hpp"
#include "Galeri_XpetraUtils.hpp"
#include "Galeri_XpetraMaps.hpp"
#endif // HAVE_IFPACK2_XPETRA

#include "Ifpack2_Relaxation.hpp"
#include "Ifpack2_UnitTestHelpers.hpp"
#include "Ifpack2_AdditiveSchwarz.hpp"
#include "Tpetra_RowMatrix.hpp"
#include "Teuchos_CommHelpers.hpp"
#include "MatrixMarket_Tpetra.hpp"

namespace { // (anonymous)

using Tpetra::global_size_t;
using Teuchos::RCP;
using std::endl;
typedef tif_utest::Node Node;

//this macro declares the unit-test-class:
TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL(Ifpack2AdditiveSchwarz, Test0, Scalar, LocalOrdinal, GlobalOrdinal)
{
//we are now in a class method declared by the above macro, and
//that method has these input arguments:
//Teuchos::FancyOStream& out, bool& success

  typedef LocalOrdinal LO;
  typedef GlobalOrdinal GO;

  typedef Tpetra::CrsMatrix<Scalar,LO,GO,Node> crs_matrix_type;
  typedef Tpetra::Map<LO,GO,Node> map_type;
  typedef Tpetra::MultiVector<Scalar,LO,GO,Node> MV;
  typedef Tpetra::RowMatrix<Scalar,LO,GO,Node> row_matrix_type;
  typedef Teuchos::ScalarTraits<Scalar> STS;
  const Scalar one = STS::one ();
  const Scalar two = STS::one () + STS::one ();

  out << "Ifpack2::AdditiveSchwarz: Test0" << endl;
  Teuchos::OSTab tab1 (out);

  global_size_t num_rows_per_proc = 5;

  out << "Creating row Map and CrsMatrix" << endl;

  RCP<const map_type> rowmap =
    tif_utest::create_tpetra_map<LO, GO, Node> (num_rows_per_proc);
  RCP<const crs_matrix_type> crsmatrix =
    tif_utest::create_test_matrix<Scalar, LO, GO, Node> (rowmap);

  out << "Creating AdditiveSchwarz instance" << endl;

  Ifpack2::AdditiveSchwarz<row_matrix_type> prec (crsmatrix);
  Teuchos::ParameterList params, zlist;

  out << "Filling in ParameterList for AdditiveSchwarz" << endl;

  zlist.set ("order_method", "rcm");
  zlist.set ("order_method_type", "local");
  params.set ("inner preconditioner name", "ILUT");
  {
    Teuchos::ParameterList innerParams;
    innerParams.set ("fact: ilut level-of-fill", 1.0);
    innerParams.set ("fact: drop tolerance", 0.0);
    params.set ("inner preconditioner parameters", innerParams);
  }
  params.set ("schwarz: overlap level", static_cast<int> (0));
  params.set ("schwarz: combine mode", "Zero");

#if defined(HAVE_IFPACK2_XPETRA) && defined(HAVE_IFPACK2_ZOLTAN2)
  params.set ("schwarz: use reordering", true);
  params.set ("schwarz: reordering list", zlist);
#else
  params.set ("schwarz: use reordering", false);
#endif

  out << "Setting AdditiveSchwarz's parameters" << endl;

  TEST_NOTHROW(prec.setParameters(params));

  out << "Testing domain and range Maps of AdditiveSchwarz" << endl;

  // FIXME (mfh 26 Jul 2015) The domain and range Maps of the
  // preconditioner don't have to be the same object; they only need
  // to be the same in the sense of Tpetra::Map::isSameAs().
  //
  //trivial tests to insist that the preconditioner's domain/range maps are
  //identically those of the matrix:
  const map_type* mtx_dom_map_ptr = &*crsmatrix->getDomainMap();
  const map_type* mtx_rng_map_ptr = &*crsmatrix->getRangeMap();
  const map_type* prec_dom_map_ptr = &*prec.getDomainMap();
  const map_type* prec_rng_map_ptr = &*prec.getRangeMap();
  TEST_EQUALITY( prec_dom_map_ptr, mtx_dom_map_ptr );
  TEST_EQUALITY( prec_rng_map_ptr, mtx_rng_map_ptr );

  out << "Calling AdditiveSchwarz's initialize()" << endl;
  prec.initialize();

  out << "Calling AdditiveSchwarz's compute()" << endl;
  prec.compute();

  MV x (rowmap, 2), y (rowmap, 2), z (rowmap, 2);
  x.putScalar (one);

  out << "Applying AdditiveSchwarz to a multivector" << endl;
  prec.apply(x, y);

  // The solution should now be full of 1/2s
  z.putScalar (one / two);

  Teuchos::ArrayRCP<const Scalar> yview = y.get1dView();
  Teuchos::ArrayRCP<const Scalar> zview = z.get1dView();

  TEST_COMPARE_FLOATING_ARRAYS(yview, zview, 4*STS::eps ());
}

// ///////////////////////////////////////////////////////////////////// //

TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL(Ifpack2AdditiveSchwarz, Test1, Scalar, LO, GO)
{
  typedef Tpetra::CrsMatrix<Scalar,LO,GO,Node> crs_matrix_type;
  typedef Tpetra::Map<LO,GO,Node> map_type;
  typedef Tpetra::RowMatrix<Scalar,LO,GO,Node> row_matrix_type;

  Teuchos::OSTab tab0 (out);
  out << "AdditiveSchwarz Test1" << endl;

  global_size_t num_rows_per_proc = 5;

  RCP<const map_type> rowmap = tif_utest::create_tpetra_map<LO,GO,Node>(num_rows_per_proc);

  // Don't run in serial.
  if(rowmap->getComm()->getSize()==1) {
    out << "It only makes sense to run this test with 1 MPI process." << endl;
    return;
  }

  RCP<const crs_matrix_type> crsmatrix = tif_utest::create_test_matrix<Scalar,LO,GO,Node>(rowmap);

  out << "Create Ifpack2::AdditiveSchwarz instance" << endl;
  Ifpack2::AdditiveSchwarz<row_matrix_type> prec (crsmatrix);
  Teuchos::ParameterList params, zlist;
  zlist.set ("order_method", "rcm");
  zlist.set ("order_method_type", "local");

  const int overlapLevel=3;
  params.set ("schwarz: overlap level", overlapLevel);
  params.set ("schwarz: combine mode", "Add");
  params.set ("inner preconditioner name", "ILUT");
  {
    Teuchos::ParameterList innerParams;
    innerParams.set ("fact: ilut level-of-fill", 1.0);
    innerParams.set ("fact: drop tolerance", 0.0);
    params.set ("inner preconditioner parameters", innerParams);
  }
#if defined(HAVE_IFPACK2_XPETRA) && defined(HAVE_IFPACK2_ZOLTAN2)
  params.set ("schwarz: use reordering", true);
  params.set ("schwarz: reordering list", zlist);
#else
  params.set ("schwarz: use reordering", false);
#endif

  TEST_NOTHROW(prec.setParameters(params));

  //trivial tests to insist that the preconditioner's domain/range maps are
  //identically those of the matrix:
  const map_type* mtx_dom_map_ptr = &*crsmatrix->getDomainMap();
  const map_type* mtx_rng_map_ptr = &*crsmatrix->getRangeMap();
  const map_type* prec_dom_map_ptr = &*prec.getDomainMap();
  const map_type* prec_rng_map_ptr = &*prec.getRangeMap();
  TEST_EQUALITY( prec_dom_map_ptr, mtx_dom_map_ptr );
  TEST_EQUALITY( prec_rng_map_ptr, mtx_rng_map_ptr );

  prec.initialize();

  prec.compute();

  //prec.describe (* Teuchos::getFancyOStream (Teuchos::rcpFromRef (std::cerr)), Teuchos::VERB_EXTREME);

  Tpetra::MultiVector<Scalar,LO,GO,Node> x(rowmap,1), y(rowmap,1), z(rowmap,1);
  x.putScalar(1);
  prec.apply(x, y);

  // The solution should now be full of 1/2s, except at processor boundaries.  If z[i] is on a processor
  // boundary, the solution will be 0.5 * K, where K is the number of processors that include i in their
  // subdomain. The matrix is a 1D operator and we know the amount of overlap, so this is easy to calculate.
  z.putScalar(0.5);
  Teuchos::ArrayRCP<Scalar> zdata = z.getDataNonConst(0);
  int mypid = rowmap->getComm()->getRank();
  if ( mypid == 0 ) {
    for (int i=0; i<overlapLevel; ++i)
      zdata[num_rows_per_proc-i-1] += 0.5;
  }
  else if (mypid == rowmap->getComm()->getSize()-1) {
    for (GO i=0; i<overlapLevel; ++i)
      zdata[i] += 0.5;
  }
  else {
    for (GO i=0; i<overlapLevel; ++i)
      zdata[i] += 0.5;
    for (GO i=0; i<overlapLevel; ++i)
      zdata[num_rows_per_proc-i-1] += 0.5;
  }
  zdata = Teuchos::null;


  Teuchos::ArrayRCP<const Scalar> yview = y.get1dView();
  Teuchos::ArrayRCP<const Scalar> zview = z.get1dView();

  TEST_COMPARE_FLOATING_ARRAYS(yview, zview, 4*Teuchos::ScalarTraits<Scalar>::eps());
}

// ///////////////////////////////////////////////////////////////////// //

TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL(Ifpack2AdditiveSchwarz, Test2, Scalar, LO, GO)
{
  typedef Tpetra::CrsMatrix<Scalar,LO,GO,Node> crs_matrix_type;
  typedef Tpetra::RowMatrix<Scalar,LO,GO,Node> row_matrix_type;

  out << "Ifpack2::AdditiveSchwarz: Test2" << endl;
  Teuchos::OSTab tab1 (out);

  global_size_t num_rows_per_proc = 5;

  const RCP<const Tpetra::Map<LO,GO,Node> > rowmap =
    tif_utest::create_tpetra_map<LO,GO,Node>(num_rows_per_proc);

  // Don't run in serial.
  if(rowmap->getComm()->getSize()==1) {
    out << endl << "This test must be run on more than one process." << endl << endl;
    return;
  }

  RCP<const crs_matrix_type> crsmatrix =
    tif_utest::create_test_matrix3<Scalar,LO,GO,Node>(rowmap);

  Ifpack2::AdditiveSchwarz<row_matrix_type> prec (crsmatrix);
  Teuchos::ParameterList params, zlist;

  const int overlapLevel=0;
  params.set ("schwarz: overlap level", overlapLevel);
  params.set ("inner preconditioner name", "ILUT");
  {
    Teuchos::ParameterList innerParams;
    innerParams.set ("fact: ilut level-of-fill", 6.0);
    innerParams.set ("fact: drop tolerance", 0.0);
    params.set ("inner preconditioner parameters", innerParams);
  }
  params.set("schwarz: use reordering",false);
  TEST_NOTHROW(prec.setParameters(params));


  prec.initialize();
  prec.compute();


  Tpetra::MultiVector<Scalar,LO,GO,Node> x(rowmap,1), y(rowmap,1), z(rowmap,1);
  x.putScalar(1);
  crsmatrix->apply(x,z);
  prec.apply(z, y);

  // The solution should now be full of 1s
  z.putScalar(1.0);

  Teuchos::ArrayRCP<const Scalar> yview = y.get1dView();
  Teuchos::ArrayRCP<const Scalar> zview = z.get1dView();

  TEST_COMPARE_FLOATING_ARRAYS(yview, zview, 4*Teuchos::ScalarTraits<Scalar>::eps());
}


// ///////////////////////////////////////////////////////////////////// //
// Test RILUK as subdomain solver for AdditiveSchwarz.
TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL(Ifpack2AdditiveSchwarz, RILUK, Scalar, LO, GO)
{
//we are now in a class method declared by the above macro, and
//that method has these input arguments:
//Teuchos::FancyOStream& out, bool& success

  typedef Tpetra::CrsMatrix<Scalar,LO,GO,Node> crs_matrix_type;
  typedef Tpetra::MultiVector<Scalar,LO,GO,Node> MV;
  typedef Tpetra::RowMatrix<Scalar,LO,GO,Node> row_matrix_type;
  typedef Tpetra::Map<LO,GO,Node> map_type;
  typedef Teuchos::ScalarTraits<Scalar> STS;

  Teuchos::OSTab tab0 (out);
  out << "Test RILUK as a subdomain solver for AdditiveSchwarz" << endl;
  Teuchos::OSTab tab1 (out);

  global_size_t num_rows_per_proc = 5;

  out << "Creating row Map and CrsMatrix" << endl;

  RCP<const map_type> rowmap =
    tif_utest::create_tpetra_map<LO,GO,Node>(num_rows_per_proc);
  RCP<const crs_matrix_type> crsmatrix =
    tif_utest::create_test_matrix<Scalar,LO,GO,Node>(rowmap);

  out << "Creating AdditiveSchwarz instance" << endl;

  Ifpack2::AdditiveSchwarz<row_matrix_type> prec (crsmatrix);
  Teuchos::ParameterList params, zlist;

  out << "Filling in ParameterList for AdditiveSchwarz" << endl;

#if defined(HAVE_IFPACK2_XPETRA) && defined(HAVE_IFPACK2_ZOLTAN2)
  params.set ("schwarz: use reordering", true);
#else
  params.set ("schwarz: use reordering", false);
#endif
  params.set ("inner preconditioner name", "RILUK");

  out << "Setting AdditiveSchwarz's parameters" << endl;

  TEST_NOTHROW(prec.setParameters(params));

  out << "Testing domain and range Maps of AdditiveSchwarz" << endl;

  //trivial tests to insist that the preconditioner's domain/range maps are
  //identically those of the matrix:
  const map_type* mtx_dom_map_ptr = &*crsmatrix->getDomainMap();
  const map_type* mtx_rng_map_ptr = &*crsmatrix->getRangeMap();
  const map_type* prec_dom_map_ptr = &*prec.getDomainMap();
  const map_type* prec_rng_map_ptr = &*prec.getRangeMap();
  TEST_EQUALITY( prec_dom_map_ptr, mtx_dom_map_ptr );
  TEST_EQUALITY( prec_rng_map_ptr, mtx_rng_map_ptr );

  out << "Calling AdditiveSchwarz's initialize()" << endl;
  prec.initialize();

  out << "Calling AdditiveSchwarz's compute()" << endl;
  prec.compute();

  MV x(rowmap,2), y(rowmap,2), z(rowmap,2);
  x.putScalar (STS::one ());

  out << "Applying AdditiveSchwarz to a multivector" << endl;
  prec.apply (x, y);

  out << "Testing result of AdditiveSchwarz's apply" << endl;

  // The solution should now be full of 1/2s
  const Scalar one = STS::one ();
  const Scalar two = one + one;
  z.putScalar (one / two);

  Teuchos::ArrayRCP<const Scalar> yview = y.get1dView();
  Teuchos::ArrayRCP<const Scalar> zview = z.get1dView();

  TEST_COMPARE_FLOATING_ARRAYS(yview, zview, 4 * STS::eps ());
}




// ///////////////////////////////////////////////////////////////////// //
// Test RILUK as subdomain solver for AdditiveSchwarz.
TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL(Ifpack2AdditiveSchwarz, RILUK_UserOrdering, Scalar, LO, GO)
{
//we are now in a class method declared by the above macro, and
//that method has these input arguments:
//Teuchos::FancyOStream& out, bool& success

  typedef Tpetra::CrsMatrix<Scalar,LO,GO,Node> crs_matrix_type;
  typedef Tpetra::MultiVector<Scalar,LO,GO,Node> MV;
  typedef Tpetra::RowMatrix<Scalar,LO,GO,Node> row_matrix_type;
  typedef Tpetra::Map<LO,GO,Node> map_type;
  typedef Teuchos::ScalarTraits<Scalar> STS;

  Teuchos::OSTab tab0 (out);
  out << "Test RILUK as a subdomain solver for AdditiveSchwarz" << endl;
  Teuchos::OSTab tab1 (out);

  global_size_t num_rows_per_proc = 5;

  out << "Creating row Map and CrsMatrix" << endl;

  RCP<const map_type> rowmap =
    tif_utest::create_tpetra_map<LO,GO,Node>(num_rows_per_proc);
  RCP<const crs_matrix_type> crsmatrix =
    tif_utest::create_test_matrix<Scalar,LO,GO,Node>(rowmap);
  size_t N = rowmap->getNodeNumElements();

  out << "Creating AdditiveSchwarz instance" << endl;

  Ifpack2::AdditiveSchwarz<row_matrix_type> prec (crsmatrix);
  Teuchos::ParameterList params;

  out << "Filling in ParameterList for AdditiveSchwarz" << endl;

#if defined(HAVE_IFPACK2_XPETRA) && defined(HAVE_IFPACK2_ZOLTAN2)
  params.set ("schwarz: use reordering", true);

  // Now let's do some user ordering.  Reverse ordering, because, why not?
  Teuchos::ArrayRCP<LO> perm(N), revperm(N);
  for(size_t i=0; i<N; i++) {
    perm[i]=(N-1)-i;
    revperm[i]=(N-1)-i;
  }
  Teuchos::ParameterList & zlist = params.sublist("schwarz: reordering list");  
  zlist.set("order_method","user");
  zlist.set("user ordering",perm);
  zlist.set("user reverse ordering",revperm);

#else
  params.set ("schwarz: use reordering", false);
#endif
  params.set ("inner preconditioner name", "RILUK");

  out << "Setting AdditiveSchwarz's parameters" << endl;

  TEST_NOTHROW(prec.setParameters(params));

  out << "Testing domain and range Maps of AdditiveSchwarz" << endl;

  //trivial tests to insist that the preconditioner's domain/range maps are
  //identically those of the matrix:
  const map_type* mtx_dom_map_ptr = &*crsmatrix->getDomainMap();
  const map_type* mtx_rng_map_ptr = &*crsmatrix->getRangeMap();
  const map_type* prec_dom_map_ptr = &*prec.getDomainMap();
  const map_type* prec_rng_map_ptr = &*prec.getRangeMap();
  TEST_EQUALITY( prec_dom_map_ptr, mtx_dom_map_ptr );
  TEST_EQUALITY( prec_rng_map_ptr, mtx_rng_map_ptr );

  out << "Calling AdditiveSchwarz's initialize()" << endl;
  prec.initialize();

  out << "Calling AdditiveSchwarz's compute()" << endl;
  prec.compute();
  prec.describe(out,Teuchos::VERB_EXTREME);

  MV x(rowmap,2), y(rowmap,2), z(rowmap,2);
  x.putScalar (STS::one ());

  out << "Applying AdditiveSchwarz to a multivector" << endl;
  prec.apply (x, y);

  out << "Testing result of AdditiveSchwarz's apply" << endl;

  // The solution should now be full of 1/2s
  const Scalar one = STS::one ();
  const Scalar two = one + one;
  z.putScalar (one / two);

  Teuchos::ArrayRCP<const Scalar> yview = y.get1dView();
  Teuchos::ArrayRCP<const Scalar> zview = z.get1dView();

  TEST_COMPARE_FLOATING_ARRAYS(yview, zview, 4 * STS::eps ());
}


// ///////////////////////////////////////////////////////////////////// //
TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL(Ifpack2AdditiveSchwarz, TestOverlap, Scalar, LO, GO)
{
  // Test that AdditiveSchwarz transfer patterns are correct.
  // A vector v such that v(i) = GID(i) is passed in to AS.  Using the IdentitySolver with any amount
  // of overlap, any subdomain reordering, and combine mode Zero, the solution should be v.

  typedef Tpetra::Map<LO, GO, Node> map_type;
  typedef Tpetra::CrsMatrix<Scalar,LO,GO,Node> crs_matrix_type;
  typedef Tpetra::RowMatrix<Scalar,LO,GO,Node> row_matrix_type;
  typedef Tpetra::MultiVector<Scalar,LO,GO,Node> MV;

  out << "Ifpack2::AdditiveSchwarz: TestOverlap" << endl;
  Teuchos::OSTab tab1 (out);

  out << "Test purpose: verify that AdditiveSchwarz transfer patterns are correct." << endl;

  global_size_t num_rows_per_proc = 10;
  RCP<const map_type> rowmap = tif_utest::create_tpetra_map<LO,GO,Node>(num_rows_per_proc);

  // Don't run in serial.
  if(rowmap->getComm()->getSize()==1) {
    out << endl << "This test must be run on more than one process." << endl << endl;
    return;
  }

  RCP<const crs_matrix_type> crsmatrix =
    tif_utest::create_test_matrix2<Scalar, LO, GO, Node> (rowmap);

  for (int i=0; i<2; ++i) {
    bool reorderSubdomains;
    if (i==0) reorderSubdomains=false;
    else      reorderSubdomains=true;

#if ! defined(HAVE_IFPACK2_XPETRA) || ! defined(HAVE_IFPACK2_ZOLTAN2)
    // mfh 19 Nov 2013: Reordering won't work (will throw an exception
    // in Ifpack2::AdditiveSchwarz) if Trilinos was not built with
    // Xpetra and Zoltan2 enabled.  Don't even bother running the test
    // in that case.
    if (reorderSubdomains) {
      continue;
    }
#endif

    for (int overlapLevel=0; overlapLevel<4; ++overlapLevel) {

      Ifpack2::AdditiveSchwarz<row_matrix_type> prec (crsmatrix);

      Teuchos::ParameterList params, zlist;
      params.set ("schwarz: overlap level", overlapLevel);
      params.set ("schwarz: use reordering",reorderSubdomains);
      zlist.set ("order_method", "rcm");
      zlist.set ("order_method_type", "local");
      params.set ("schwarz: reordering list", zlist);
      params.set("schwarz: combine mode", "Zero");
      params.set ("inner preconditioner name", "IDENTITY");

      prec.setParameters(params);

      prec.initialize();
      prec.compute();


      MV x (rowmap, 1), y (rowmap, 1);
      Teuchos::ArrayRCP<Scalar> xData = x.getDataNonConst(0);
      for (GO j=0; j<xData.size(); ++j)
        xData[j] = rowmap->getGlobalElement(j);
      xData = Teuchos::null;
      prec.apply(x, y);

      out << prec.description() << endl;

      // The solution should now be full of GIDs.
      Teuchos::ArrayRCP<const Scalar> xview = x.get1dView();
      Teuchos::ArrayRCP<const Scalar> yview = y.get1dView();

      //out.setOutputToRootOnly(-1);
      //x.describe(out,Teuchos::VERB_EXTREME);
      //y.describe(out,Teuchos::VERB_EXTREME);

      TEST_COMPARE_FLOATING_ARRAYS(xview, yview, 4*Teuchos::ScalarTraits<Scalar>::eps());
    }
  }
}

// ///////////////////////////////////////////////////////////////////// //
#if defined(HAVE_IFPACK2_AMESOS2) && defined(HAVE_IFPACK2_XPETRA) && (defined(HAVE_AMESOS2_SUPERLU) || defined(HAVE_AMESOS2_KLU2))
// Test sparse direct solver as subdomain solver for AdditiveSchwarz.
TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL(Ifpack2AdditiveSchwarz, SparseDirectSolver, SC, LO, GO)
{
  using Teuchos::ParameterList;
  typedef Node NT;
  typedef typename Tpetra::Vector<SC,LO,GO,NT>::mag_type mag_type;
  typedef Teuchos::ScalarTraits<SC> STS;

  typedef Tpetra::CrsMatrix<SC,LO,GO,NT>   crs_matrix_type;
  typedef Tpetra::RowMatrix<SC,LO,GO,NT>   row_matrix_type;
  typedef Tpetra::MultiVector<SC,LO,GO,NT> MV;
  typedef Tpetra::Map<LO,GO,NT>            map_type;

  typedef Xpetra::TpetraCrsMatrix<SC,LO,GO,NT> XCrsType;
  typedef Xpetra::Map<LO,GO,NT>                XMapType;
  typedef Xpetra::MultiVector<SC,LO,GO,NT>     XMVectorType;

  using Teuchos::REDUCE_MIN;
  using Teuchos::outArg;
  using Teuchos::reduceAll;
  int lclSuccess = 1; // will set again below
  int gblSuccess = 0; // output argument


  out << "Test Amesos2 sparse direct solver"
#if defined(HAVE_AMESOS2_SUPERLU)
  << " (SuperLU)"
#elif defined(HAVE_AMESOS2_KLU2)
  << " (KLU2)"
#endif
  << " as subdomain solver for AdditiveSchwarz" << endl;

  // Generate the matrix using Galeri.  Galeri wraps it in an Xpetra
  // matrix, so after it finishes, ask it for the Tpetra matrix.
  RCP<const Teuchos::Comm<int> > comm = Tpetra::getDefaultComm ();
  Teuchos::CommandLineProcessor clp;
  GO nx = 10, ny=10, nz=10;
  Galeri::Xpetra::Parameters<GO> GaleriParameters (clp, nx, ny, nz, "Laplace2D");
  Xpetra::Parameters xpetraParameters (clp);
  ParameterList GaleriList = GaleriParameters.GetParameterList ();

  RCP<XMapType> xmap =
    Galeri::Xpetra::CreateMap<LO, GO, Node> (xpetraParameters.GetLib (),
                                             "Cartesian2D", comm, GaleriList);
  RCP<Galeri::Xpetra::Problem<XMapType,XCrsType,XMVectorType> > Pr =
    Galeri::Xpetra::BuildProblem<SC,LO,GO,XMapType,XCrsType,XMVectorType> (std::string("Laplace2D"),
                                                                           xmap, GaleriList);

  RCP<XCrsType> XA = Pr->BuildMatrix ();
  RCP<crs_matrix_type> A = XA->getTpetra_CrsMatrixNonConst ();
  TEST_INEQUALITY(A, Teuchos::null);

  RCP<const map_type> rowmap = A->getRowMap ();

  out << "Creating AdditiveSchwarz instance" << endl;

  Ifpack2::AdditiveSchwarz<row_matrix_type> prec (A);
  ParameterList params, zlist;

  out << "Filling in ParameterList for AdditiveSchwarz" << endl;

  params.set ("schwarz: overlap level", 2);

#if defined(HAVE_IFPACK2_XPETRA) && defined(HAVE_IFPACK2_ZOLTAN2)
  params.set ("schwarz: use reordering", true);
#else
  params.set ("schwarz: use reordering", false);
#endif
  //params.set ("inner preconditioner name", "AMESOS2");
  params.set ("subdomain solver name", "AMESOS2");

  /*
    Here is how to set SuperLU options in the Amesos2 subdomain solve:

    <Parameter name="subdomain solver name"  type="string"   value="AMESOS2"/>
                  OR
    <Parameter name="inner preconditioner name"  type="string"   value="AMESOS2"/>
    <ParameterList name="subdomain solver parameters">
      <Parameter name="Amesos2 solver name"   type="string"   value="klu2"/>  //or superlu, superludist, etc.
                                                                              //if omitted, defaults to superlu
      <ParameterList name="Amesos2">
        <ParameterList name="SuperLU">
          <Parameter name="ILU_Flag"           type="bool"   value="false"/>  //set to true to use SuperLU's ILUTP
        </ParameterList>
      </ParameterList>
    </ParameterList>

  */

  ParameterList &subdomainList = params.sublist("subdomain solver parameters");
#if defined(HAVE_AMESOS2_SUPERLU)
  subdomainList.set("Amesos2 solver name","superlu");
  ParameterList &amesos2List = subdomainList.sublist("Amesos2"); // sublist required by Amesos2
  ParameterList &superluList = amesos2List.sublist("SuperLU");
  superluList.set("ILU_Flag",false);
#elif defined(HAVE_AMESOS2_KLU2)
  subdomainList.set("Amesos2 solver name","klu");
  subdomainList.sublist("Amesos2"); // sublist required by Amesos2
#else
#  error "This test requires Amesos2 to have support for either SuperLU or KLU."
#endif

  out << "Setting AdditiveSchwarz's parameters" << endl;

  std::ostringstream ps;
  int indent = 4;
  params.print(ps, indent);
  out << ps.str() << endl;

  TEST_NOTHROW(prec.setParameters(params));
  lclSuccess = success ? 1 : 0;
  gblSuccess = 0; // output argument
  reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
  TEST_EQUALITY_CONST( gblSuccess, 1 );
  if (gblSuccess != 1) {
    out << "prec.setParameters(params) threw an exception on some process!" << endl;
    return; // no sense in continuing
  }

  out << "Testing domain and range Maps of AdditiveSchwarz" << endl;

  //trivial tests to insist that the preconditioner's domain/range maps are
  //identically those of the matrix:
  const map_type* mtx_dom_map_ptr = &*A->getDomainMap();
  const map_type* mtx_rng_map_ptr = &*A->getRangeMap();
  const map_type* prec_dom_map_ptr = &*prec.getDomainMap();
  const map_type* prec_rng_map_ptr = &*prec.getRangeMap();
  TEST_EQUALITY( prec_dom_map_ptr, mtx_dom_map_ptr );
  TEST_EQUALITY( prec_rng_map_ptr, mtx_rng_map_ptr );

  out << endl << "solve using AdditiveSchwarz with sparse direct method as the subdomain solve" << endl;
  out << "Calling AdditiveSchwarz::initialize()" << endl;
  TEST_NOTHROW( prec.initialize() );
  lclSuccess = success ? 1 : 0;
  gblSuccess = 0; // output argument
  reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
  TEST_EQUALITY_CONST( gblSuccess, 1 );
  if (gblSuccess != 1) {
    out << "prec.initialize() threw an exception on some process!" << endl;
    return; // no sense in continuing
  }

  out << "Calling AdditiveSchwarz::compute()" << endl;
  TEST_NOTHROW( prec.compute() );
  lclSuccess = success ? 1 : 0;
  gblSuccess = 0; // output argument
  reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
  TEST_EQUALITY_CONST( gblSuccess, 1 );
  if (gblSuccess != 1) {
    out << "prec.compute() threw an exception on some process!" << endl;
    //return; // no sense in continuing
  }

  out << "Generating input and output (multi)vectors" << endl;
  const int numCols = 2;
  MV x(rowmap, numCols), y(rowmap, numCols);
  x.randomize();

  if (gblSuccess == 1) {
    out << "Calling AdditiveSchwarz::apply()" << endl;
    TEST_NOTHROW( prec.apply (x, y) );
    lclSuccess = success ? 1 : 0;
    gblSuccess = 0; // output argument
    reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
    TEST_EQUALITY_CONST( gblSuccess, 1 );
    if (gblSuccess != 1) {
      out << "prec.apply(x,y) threw an exception on some process!" << endl;
      return; // no sense in continuing
    }
  }

  // Now switch to dense direct solves on the subdomains.
  out << endl << "solve using AdditiveSchwarz with dense direct method as the subdomain solve" << endl;
  params.set ("subdomain solver name", "DENSE");
  prec.setParameters(params);
  out << "Calling AdditiveSchwarz::initialize()" << endl;
  prec.initialize();
  out << "Calling AdditiveSchwarz::compute()" << endl;
  prec.compute();
  MV z(rowmap,numCols);
  out << "Calling AdditiveSchwarz::apply()" << endl;

  prec.apply (x, z);

  if (gblSuccess == 1) {
    out << "Comparing results of two solves" << endl;
  }
  Teuchos::Array<mag_type> ynorms (y.getNumVectors ());
  Teuchos::Array<mag_type> znorms (z.getNumVectors ());
  if (gblSuccess == 1) {
    y.norm2 (ynorms ());
  }
  z.norm2 (znorms ());
  if (gblSuccess == 1) {
    out << "solution norm, sparse direct solve: " << std::setprecision(7) << ynorms[0] << endl;
  }
  out << "solution norm,  dense direct solve: " << std::setprecision(7) << znorms[0] << endl;
  if (gblSuccess == 1) {
    TEST_FLOATING_EQUALITY(ynorms[0], znorms[0], 10* STS::eps ());
  }
}
#endif

// ///////////////////////////////////////////////////////////////////// //
// Test multiple sweeps of additive Schwarz.
TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL(Ifpack2AdditiveSchwarz, MultipleSweeps, Scalar, LocalOrdinal, GlobalOrdinal)
{
  using Teuchos::ArrayRCP;

  typedef LocalOrdinal LO;
  typedef GlobalOrdinal GO;

  typedef Tpetra::CrsMatrix<Scalar,LO,GO,Node>   crs_matrix_type;
  typedef Tpetra::RowMatrix<Scalar,LO,GO,Node>   row_matrix_type;
  typedef Tpetra::Map<LO,GO,Node>                map_type;
  typedef Tpetra::MultiVector<Scalar,LO,GO,Node> multivector_type;
  typedef Teuchos::ScalarTraits<Scalar>          STS;
  typedef typename multivector_type::mag_type    mag_type;

  out << "Test multiple sweeps of AdditiveSchwarz" << endl;
  Teuchos::OSTab tab1 (out);

  global_size_t num_rows_per_proc = 100;

  out << "Creating row Map and CrsMatrix" << endl;

  RCP<const map_type> rowmap = tif_utest::create_tpetra_map<LO,GO,Node>(num_rows_per_proc);
  RCP<const crs_matrix_type> crsmatrix = tif_utest::create_banded_matrix<Scalar,LO,GO,Node>(rowmap,5);

  out << "Creating AdditiveSchwarz instance" << endl;

  Ifpack2::AdditiveSchwarz<row_matrix_type> prec1(crsmatrix), prec2(crsmatrix);
  Teuchos::ParameterList params1, params2, zlist;

  out << "Setting AdditiveSchwarz's parameters" << endl;

  // prec1 assumes initial guess is zero
# if defined(HAVE_IFPACK2_XPETRA) && defined(HAVE_IFPACK2_ZOLTAN2)
  params1.set ("schwarz: use reordering",         true);
# else
  params1.set ("schwarz: use reordering",         false);
# endif
  params1.set ("inner preconditioner name",       "RILUK");
  params1.set ("schwarz: zero starting solution", true);
  TEST_NOTHROW(prec1.setParameters(params1));

  // prec2 assumes initial guess is nonzero
# if defined(HAVE_IFPACK2_XPETRA) && defined(HAVE_IFPACK2_ZOLTAN2)
  params2.set ("schwarz: use reordering",         true);
# else
  params2.set ("schwarz: use reordering",         false);
# endif
  params2.set ("inner preconditioner name",       "RILUK");
  params2.set ("schwarz: zero starting solution", false);
  TEST_NOTHROW(prec2.setParameters(params2));

  out << "Testing domain and range Maps of AdditiveSchwarz" << endl;

  //trivial tests to insist that the preconditioner's domain/range maps are
  //identically those of the matrix:
  const map_type* mtx_dom_map_ptr = &*crsmatrix->getDomainMap();
  const map_type* mtx_rng_map_ptr = &*crsmatrix->getRangeMap();
  const map_type* prec_dom_map_ptr = &*prec1.getDomainMap();
  const map_type* prec_rng_map_ptr = &*prec1.getRangeMap();
  TEST_EQUALITY( prec_dom_map_ptr, mtx_dom_map_ptr );
  TEST_EQUALITY( prec_rng_map_ptr, mtx_rng_map_ptr );

  out << "Calling AdditiveSchwarz's initialize()" << endl;
  prec1.initialize();
  prec2.initialize();

  out << "Calling AdditiveSchwarz's compute()" << endl;
  prec1.compute();
  prec2.compute();

  multivector_type b(rowmap,2),
                   y1(rowmap,2),
                   y2(rowmap,2),
                   z1(rowmap,2),
                   z2(rowmap,2);
  //b.putScalar(1);
  b.randomize();
  y1.putScalar(5);
  y2.putScalar(5);

  // Test that using a nonzero initial guess works.
  out << "Applying AdditiveSchwarz to a multivector" << endl;
  prec1.apply (b, y1);
  prec2.apply (b, y2);

  Teuchos::Array<mag_type> y1norms(2), y2norms(2);
  y1.norm2(y1norms());
  y2.norm2(y2norms());

  out << "schwarz: zero starting solution = "
      << params1.get<bool>("schwarz: zero starting solution")
      << ", ||soln||=" << y1norms[0] << endl;
  out << "schwarz: zero starting solution = "
      << params2.get<bool>("schwarz: zero starting solution")
      << ", ||soln||=" << y2norms[0] << endl;
  //norm(y1) should be different than norm(y2)
  TEST_EQUALITY( (y1norms[0] != y2norms[0]), true);
  TEST_EQUALITY( (y1norms[1] != y2norms[1]), true);


  // Test multiple sweeps.

  const int numSweeps = 5;
  params1.set ("inner preconditioner name",       "RILUK");
  params1.set ("schwarz: zero starting solution", true);
  params1.set ("schwarz: overlap level", 1);
  params1.set ("schwarz: num iterations", 1);
  prec1.setParameters(params1);
  prec1.initialize();
  prec1.compute();

  params2.set ("inner preconditioner name",       "RILUK");
  params2.set ("schwarz: zero starting solution", true);
  params2.set ("schwarz: overlap level", 1);
  params2.set ("schwarz: num iterations", numSweeps);
  prec2.setParameters(params2);
  prec2.initialize();
  prec2.compute();

  y2.putScalar(0);
  prec2.apply(b, y2);

  // Compare against Richardson iterations.
  multivector_type r(rowmap,2), c(rowmap,2);
  Teuchos::Array<mag_type> cnorms(2), rnorms(2);
  y1.putScalar(0);
  for (int i=0; i < numSweeps; ++i) {
    //r=b-A*y1
    Tpetra::deep_copy(r, b);
    crsmatrix->apply(y1, r, Teuchos::NO_TRANS, -STS::one(), STS::one());
    r.norm2(rnorms());
    //solve Ac=r
    c.putScalar(0);
    prec1.apply(r, c);
    //y1 = y1 + c
    y1.update(STS::one(), c, STS::one());
    y1.norm2(y1norms());
    c.norm2(cnorms());
    out << "iter " << i+1
        << "||res||=" << rnorms[0]
        << ", ||correction||=" << cnorms[0]
        << ", ||soln||=" << y1norms[0] << endl;
  }
  y1.norm2(y1norms());
  y2.norm2(y2norms());
  TEST_FLOATING_EQUALITY(y1norms[0], y2norms[0], 10*STS::eps());
  TEST_FLOATING_EQUALITY(y1norms[1], y2norms[1], 10*STS::eps());

# ifdef TODO_TESTING
  //TODO test multiple sweeps vs. handrolled Richardson iteration
  out << "Testing result of AdditiveSchwarz's apply" << endl;

  // The solution should now be full of 1/2s
  z.putScalar(0.5);

  ArrayRCP<const Scalar> yview = y.get1dView();
  ArrayRCP<const Scalar> zview = z.get1dView();

  TEST_COMPARE_FLOATING_ARRAYS(yview, zview, 4*STS::eps());
# endif
}

//test correctness with ILU subdomain solver
//use MATLAB implementation as gold standard
TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL(Ifpack2AdditiveSchwarz, ILU_Overlap, Scalar, LocalOrdinal, GlobalOrdinal)
{ 
  using std::string;
  typedef LocalOrdinal LO;
  typedef GlobalOrdinal GO;

  typedef Tpetra::CrsMatrix<Scalar,LO,GO,Node>   crs_matrix_type;
  typedef Tpetra::RowMatrix<Scalar,LO,GO,Node>   row_matrix_type;
  typedef Tpetra::Map<LO,GO,Node>                map_type;
  typedef Tpetra::MultiVector<Scalar,LO,GO,Node> multivector_type;
  typedef Teuchos::ScalarTraits<Scalar>          STS;
  typedef Tpetra::MatrixMarket::Reader<crs_matrix_type> Reader;

  RCP<const Teuchos::Comm<int> > comm = Tpetra::getDefaultComm ();

  //this test can only be run on 4 procs
  if(comm->getSize() != 4)
  {
    out << "AdditiveSchwarz/ILU correctness test must be run on 4 procs.\n";
    return;
  }

  string rowMapFile = "AdditiveSchwarzMatlab_rowmap.mm";
  RCP<const map_type> rowMap = Reader::readMapFile(rowMapFile, comm);
  string matrixFile = "AdditiveSchwarzMatlab_A.mm";
  RCP<crs_matrix_type> A = Reader::readSparseFile(matrixFile, comm);
  string rhsFile = "AdditiveSchwarzMatlab_rhs.mm";
  RCP<multivector_type> rhs = Reader::readDenseFile(rhsFile, comm, rowMap);
  string solFile = "AdditiveSchwarzILU_O_sol.mm";
  RCP<multivector_type> sol = Reader::readDenseFile(solFile, comm, rowMap);
  RCP<multivector_type> x = rcp(new multivector_type(rowMap, 1));

  //set up additive schwarz
  Ifpack2::AdditiveSchwarz<row_matrix_type> prec (A);
  Teuchos::ParameterList params;

  out << "Filling in ParameterList for AdditiveSchwarz with ILU" << endl;

  //ILU(0) by using RILUK with fill level = 0 and drop tol = 0
  params.set ("inner preconditioner name", "RILUK");
  {
    Teuchos::ParameterList innerParams;
    innerParams.set ("fact: ilut level-of-fill", 0.0);
    innerParams.set ("fact: drop tolerance", 0.0);
    params.set ("inner preconditioner parameters", innerParams);
  }
  params.set ("schwarz: zero starting solution", true);
  params.set ("schwarz: overlap level", 1);
  params.set ("schwarz: combine mode", "ZERO");
  params.set ("schwarz: num iterations", 5);
  params.set ("schwarz: use reordering", false);

  TEST_NOTHROW(prec.setParameters(params));

  out << "Calling AdditiveSchwarz's initialize()" << endl;
  prec.initialize();

  out << "Calling AdditiveSchwarz's compute()" << endl;
  prec.compute();

  out << "Applying AdditiveSchwarz to a multivector" << endl;
  prec.apply(*rhs, *x);

  TEST_COMPARE_FLOATING_ARRAYS(sol->get1dView(), x->get1dView(), 4*STS::eps ());
}


TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL(Ifpack2AdditiveSchwarz, ILU_NonOverlap, Scalar, LocalOrdinal, GlobalOrdinal)
{ 
  using std::string;
  typedef LocalOrdinal LO;
  typedef GlobalOrdinal GO;

  typedef Tpetra::CrsMatrix<Scalar,LO,GO,Node>   crs_matrix_type;
  typedef Tpetra::RowMatrix<Scalar,LO,GO,Node>   row_matrix_type;
  typedef Tpetra::Map<LO,GO,Node>                map_type;
  typedef Tpetra::MultiVector<Scalar,LO,GO,Node> multivector_type;
  typedef Teuchos::ScalarTraits<Scalar>          STS;
  typedef Tpetra::MatrixMarket::Reader<crs_matrix_type> Reader;

  RCP<const Teuchos::Comm<int> > comm = Tpetra::getDefaultComm ();

  //this test can only be run on 4 procs
  if(comm->getSize() != 4)
  {
    out << "AdditiveSchwarz/ILU correctness test must be run on 4 procs.\n";
    return;
  }

  string mapFile = "AdditiveSchwarzMatlab_rowmap.mm";
  RCP<const map_type> rowMap = Reader::readMapFile(mapFile, comm);
  string matrixFile = "AdditiveSchwarzMatlab_A.mm";
  RCP<crs_matrix_type> A = Reader::readSparseFile(matrixFile, comm);
  string rhsFile = "AdditiveSchwarzMatlab_rhs.mm";
  RCP<multivector_type> rhs = Reader::readDenseFile(rhsFile, comm, rowMap);
  string solFile = "AdditiveSchwarzILU_NO_sol.mm";
  RCP<multivector_type> sol = Reader::readDenseFile(solFile, comm, rowMap);
  RCP<multivector_type> x = rcp(new multivector_type(rowMap, 1));

  //set up additive schwarz
  Ifpack2::AdditiveSchwarz<row_matrix_type> prec (A);
  Teuchos::ParameterList params;

  out << "Filling in ParameterList for AdditiveSchwarz with ILU" << endl;

  //RILUK with fill = 0 is ILU(0)
  params.set ("inner preconditioner name", "RILUK");
  {
    Teuchos::ParameterList innerParams;
    innerParams.set ("fact: ilut level-of-fill", 0.0);
    innerParams.set ("fact: drop tolerance", 0.0);
    params.set ("inner preconditioner parameters", innerParams);
  }
  params.set ("schwarz: zero starting solution", true);
  params.set ("schwarz: overlap level", 0);
  params.set ("schwarz: combine mode", "ZERO");
  params.set ("schwarz: num iterations", 5);
  params.set ("schwarz: use reordering", false);

  TEST_NOTHROW(prec.setParameters(params));

  out << "Calling AdditiveSchwarz's initialize()" << endl;
  prec.initialize();

  out << "Calling AdditiveSchwarz's compute()" << endl;
  prec.compute();

  out << "Applying AdditiveSchwarz to a multivector" << endl;
  prec.apply(*rhs, *x);

  auto solView = sol->get1dView();
  auto xView = x->get1dView();

  TEST_COMPARE_FLOATING_ARRAYS(solView, xView, 4 * STS::eps ());
}

//test correctness with SGS subdomain solver (overlap = 0)
//use MATLAB implementation as gold standard
TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL(Ifpack2AdditiveSchwarz, SGS_NonOverlap, Scalar, LocalOrdinal, GlobalOrdinal)
{ 
  using std::string;
  typedef LocalOrdinal LO;
  typedef GlobalOrdinal GO;

  typedef Tpetra::CrsMatrix<Scalar,LO,GO,Node>   crs_matrix_type;
  typedef Tpetra::RowMatrix<Scalar,LO,GO,Node>   row_matrix_type;
  typedef Tpetra::Map<LO,GO,Node>                map_type;
  typedef Tpetra::MultiVector<Scalar,LO,GO,Node> multivector_type;
  typedef Teuchos::ScalarTraits<Scalar>          STS;
  typedef Tpetra::MatrixMarket::Reader<crs_matrix_type> Reader;

  RCP<const Teuchos::Comm<int> > comm = Tpetra::getDefaultComm ();

  //this test can only be run on 4 procs
  if(comm->getSize() != 4)
  {
    out << "AdditiveSchwarz/SGS correctness test must be run on 4 procs.\n";
    return;
  }

  string rowMapFile = "AdditiveSchwarzMatlab_rowmap.mm";
  RCP<const map_type> rowMap = Reader::readMapFile(rowMapFile, comm);
  string matrixFile = "AdditiveSchwarzMatlab_A.mm";
  RCP<crs_matrix_type> A = Reader::readSparseFile(matrixFile, comm);
  string rhsFile = "AdditiveSchwarzMatlab_rhs.mm";
  RCP<multivector_type> rhs = Reader::readDenseFile(rhsFile, comm, rowMap);
  string solFile = "AdditiveSchwarzSGS_NO_sol.mm";
  RCP<multivector_type> sol = Reader::readDenseFile(solFile, comm, rowMap);
  RCP<multivector_type> x = rcp(new multivector_type(rowMap, 1));

  //set up additive schwarz
  Ifpack2::AdditiveSchwarz<row_matrix_type> prec (A);
  Teuchos::ParameterList params;

  out << "Filling in ParameterList for AdditiveSchwarz with SGS" << endl;

  //Inner solver is 1 SGS sweep 
  params.set ("inner preconditioner name", "RELAXATION");
  {
    Teuchos::ParameterList innerParams;
    innerParams.set ("relaxation: type", "Symmetric Gauss-Seidel");
    innerParams.set ("relaxation: sweeps", 1);
    innerParams.set ("relaxation: zero starting solution", true);
    params.set ("inner preconditioner parameters", innerParams);
  }
  params.set ("schwarz: zero starting solution", false);
  params.set ("schwarz: overlap level", 0);
  params.set ("schwarz: combine mode", "ZERO");
  params.set ("schwarz: num iterations", 5);
  params.set ("schwarz: use reordering", false);

  TEST_NOTHROW(prec.setParameters(params));

  out << "Calling AdditiveSchwarz's initialize()" << endl;
  prec.initialize();

  out << "Calling AdditiveSchwarz's compute()" << endl;
  prec.compute();

  out << "Applying AdditiveSchwarz to a multivector" << endl;
  prec.apply(*rhs, *x);

  TEST_COMPARE_FLOATING_ARRAYS(sol->get1dView(), x->get1dView(), 4*STS::eps ());
}

//test correctness with SGS subdomain solver (overlap = 1)
//use MATLAB implementation as gold standard
TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL(Ifpack2AdditiveSchwarz, SGS_Overlap, Scalar, LocalOrdinal, GlobalOrdinal)
{ 
  using std::string;
  typedef LocalOrdinal LO;
  typedef GlobalOrdinal GO;

  typedef Tpetra::CrsMatrix<Scalar,LO,GO,Node>   crs_matrix_type;
  typedef Tpetra::RowMatrix<Scalar,LO,GO,Node>   row_matrix_type;
  typedef Tpetra::Map<LO,GO,Node>                map_type;
  typedef Tpetra::MultiVector<Scalar,LO,GO,Node> multivector_type;
  typedef Teuchos::ScalarTraits<Scalar>          STS;
  typedef Tpetra::MatrixMarket::Reader<crs_matrix_type> Reader;

  RCP<const Teuchos::Comm<int> > comm = Tpetra::getDefaultComm ();

  //this test can only be run on 4 procs
  if(comm->getSize() != 4)
  {
    out << "AdditiveSchwarz/SGS correctness test must be run on 4 procs.\n";
    return;
  }

  string rowMapFile = "AdditiveSchwarzMatlab_rowmap.mm";
  RCP<const map_type> rowMap = Reader::readMapFile(rowMapFile, comm);
  string matrixFile = "AdditiveSchwarzMatlab_A.mm";
  RCP<crs_matrix_type> A = Reader::readSparseFile(matrixFile, comm);
  string rhsFile = "AdditiveSchwarzMatlab_rhs.mm";
  RCP<multivector_type> rhs = Reader::readDenseFile(rhsFile, comm, rowMap);
  string solFile = "AdditiveSchwarzSGS_O_sol.mm";
  RCP<multivector_type> sol = Reader::readDenseFile(solFile, comm, rowMap);
  RCP<multivector_type> x = rcp(new multivector_type(rowMap, 1));

  //set up additive schwarz
  Ifpack2::AdditiveSchwarz<row_matrix_type> prec (A);
  Teuchos::ParameterList params;

  out << "Filling in ParameterList for AdditiveSchwarz with SGS" << endl;

  //Inner solver is 1 SGS sweep 
  params.set ("inner preconditioner name", "RELAXATION");
  {
    Teuchos::ParameterList innerParams;
    innerParams.set ("relaxation: type", "Symmetric Gauss-Seidel");
    innerParams.set ("relaxation: sweeps", 1);
    innerParams.set ("relaxation: zero starting solution", true);
    params.set ("inner preconditioner parameters", innerParams);
  }
  params.set ("schwarz: zero starting solution", false);
  params.set ("schwarz: overlap level", 1);
  params.set ("schwarz: combine mode", "ZERO");
  params.set ("schwarz: num iterations", 5);
  params.set ("schwarz: use reordering", false);

  TEST_NOTHROW(prec.setParameters(params));

  out << "Calling AdditiveSchwarz's initialize()" << endl;
  prec.initialize();

  out << "Calling AdditiveSchwarz's compute()" << endl;
  prec.compute();

  out << "Applying AdditiveSchwarz to a multivector" << endl;
  prec.apply(*rhs, *x);

  TEST_COMPARE_FLOATING_ARRAYS(sol->get1dView(), x->get1dView(), 4*STS::eps ());
}

TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL(Ifpack2AdditiveSchwarz, FastILU, Scalar, LocalOrdinal, GlobalOrdinal)
{
  typedef LocalOrdinal LO;
  typedef GlobalOrdinal GO;

  typedef Tpetra::CrsMatrix<Scalar,LO,GO,Node> crs_matrix_type;
  typedef Tpetra::Map<LO,GO,Node> map_type;
  typedef Tpetra::MultiVector<Scalar,LO,GO,Node> MV;
  typedef Tpetra::RowMatrix<Scalar,LO,GO,Node> row_matrix_type;
  typedef Teuchos::ScalarTraits<Scalar> STS;
  const Scalar one = STS::one ();
  const Scalar two = STS::one () + STS::one ();

  out << "Ifpack2::AdditiveSchwarz: FastILU" << endl;
  Teuchos::OSTab tab1 (out);

  global_size_t num_rows_per_proc = 500;

  out << "Creating row Map and CrsMatrix" << endl;

  RCP<const map_type> rowmap =
    tif_utest::create_tpetra_map<LO, GO, Node> (num_rows_per_proc);
  RCP<const crs_matrix_type> crsmatrix =
    tif_utest::create_test_matrix<Scalar, LO, GO, Node> (rowmap);

  out << "Creating AdditiveSchwarz instance" << endl;

  Ifpack2::AdditiveSchwarz<row_matrix_type> prec (crsmatrix);
  Teuchos::ParameterList params, zlist;

  out << "Filling in ParameterList for AdditiveSchwarz" << endl;

  zlist.set ("order_method", "rcm");
  zlist.set ("order_method_type", "local");
  params.set ("inner preconditioner name", "FAST_ILU");
  {
    Teuchos::ParameterList innerParams;
    params.set ("inner preconditioner parameters", innerParams);
  }
  params.set ("schwarz: overlap level", 0);
  params.set ("schwarz: combine mode", "Zero");

#if defined(HAVE_IFPACK2_XPETRA) && defined(HAVE_IFPACK2_ZOLTAN2)
  params.set ("schwarz: use reordering", true);
  params.set ("schwarz: reordering list", zlist);
#else
  params.set ("schwarz: use reordering", false);
#endif

  out << "Setting AdditiveSchwarz's parameters" << endl;

  TEST_NOTHROW(prec.setParameters(params));

  out << "Testing domain and range Maps of AdditiveSchwarz" << endl;

  // FIXME (mfh 26 Jul 2015) The domain and range Maps of the
  // preconditioner don't have to be the same object; they only need
  // to be the same in the sense of Tpetra::Map::isSameAs().
  //
  //trivial tests to insist that the preconditioner's domain/range maps are
  //identically those of the matrix:
  const map_type* mtx_dom_map_ptr = &*crsmatrix->getDomainMap();
  const map_type* mtx_rng_map_ptr = &*crsmatrix->getRangeMap();
  const map_type* prec_dom_map_ptr = &*prec.getDomainMap();
  const map_type* prec_rng_map_ptr = &*prec.getRangeMap();
  TEST_EQUALITY( prec_dom_map_ptr, mtx_dom_map_ptr );
  TEST_EQUALITY( prec_rng_map_ptr, mtx_rng_map_ptr );

  out << "Calling AdditiveSchwarz's initialize()" << endl;
  prec.initialize();

  out << "Calling AdditiveSchwarz's compute()" << endl;
  prec.compute();

  MV x (rowmap, 2), y (rowmap, 2), z (rowmap, 2);
  x.putScalar (one);

  out << "Applying AdditiveSchwarz to a multivector" << endl;
  prec.apply(x, y);

  // The solution should now be full of 1/2s
  z.putScalar (one / two);

  Teuchos::ArrayRCP<const Scalar> yview = y.get1dView();
  Teuchos::ArrayRCP<const Scalar> zview = z.get1dView();

  TEST_COMPARE_FLOATING_ARRAYS(yview, zview, 4*STS::eps ());
}

TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL(Ifpack2AdditiveSchwarz, FastIC, Scalar, LocalOrdinal, GlobalOrdinal)
{
  typedef LocalOrdinal LO;
  typedef GlobalOrdinal GO;

  typedef Tpetra::CrsMatrix<Scalar,LO,GO,Node> crs_matrix_type;
  typedef Tpetra::Map<LO,GO,Node> map_type;
  typedef Tpetra::MultiVector<Scalar,LO,GO,Node> MV;
  typedef Tpetra::RowMatrix<Scalar,LO,GO,Node> row_matrix_type;
  typedef Teuchos::ScalarTraits<Scalar> STS;
  const Scalar one = STS::one ();
  const Scalar two = STS::one () + STS::one ();

  out << "Ifpack2::AdditiveSchwarz: FastIC" << endl;
  Teuchos::OSTab tab1 (out);

  global_size_t num_rows_per_proc = 500;

  out << "Creating row Map and CrsMatrix" << endl;

  RCP<const map_type> rowmap =
    tif_utest::create_tpetra_map<LO, GO, Node> (num_rows_per_proc);
  RCP<const crs_matrix_type> crsmatrix =
    tif_utest::create_test_matrix<Scalar, LO, GO, Node> (rowmap);

  out << "Creating AdditiveSchwarz instance" << endl;

  Ifpack2::AdditiveSchwarz<row_matrix_type> prec (crsmatrix);
  Teuchos::ParameterList params, zlist;

  out << "Filling in ParameterList for AdditiveSchwarz" << endl;

  zlist.set ("order_method", "rcm");
  zlist.set ("order_method_type", "local");
  params.set ("inner preconditioner name", "FAST_IC");
  {
    Teuchos::ParameterList innerParams;
    params.set ("inner preconditioner parameters", innerParams);
  }
  params.set ("schwarz: overlap level", static_cast<int> (0));
  params.set ("schwarz: combine mode", "Zero");

#if defined(HAVE_IFPACK2_XPETRA) && defined(HAVE_IFPACK2_ZOLTAN2)
  params.set ("schwarz: use reordering", true);
  params.set ("schwarz: reordering list", zlist);
#else
  params.set ("schwarz: use reordering", false);
#endif

  out << "Setting AdditiveSchwarz's parameters" << endl;

  TEST_NOTHROW(prec.setParameters(params));

  out << "Testing domain and range Maps of AdditiveSchwarz" << endl;

  // FIXME (mfh 26 Jul 2015) The domain and range Maps of the
  // preconditioner don't have to be the same object; they only need
  // to be the same in the sense of Tpetra::Map::isSameAs().
  //
  //trivial tests to insist that the preconditioner's domain/range maps are
  //identically those of the matrix:
  const map_type* mtx_dom_map_ptr = &*crsmatrix->getDomainMap();
  const map_type* mtx_rng_map_ptr = &*crsmatrix->getRangeMap();
  const map_type* prec_dom_map_ptr = &*prec.getDomainMap();
  const map_type* prec_rng_map_ptr = &*prec.getRangeMap();
  TEST_EQUALITY( prec_dom_map_ptr, mtx_dom_map_ptr );
  TEST_EQUALITY( prec_rng_map_ptr, mtx_rng_map_ptr );

  out << "Calling AdditiveSchwarz's initialize()" << endl;
  prec.initialize();

  out << "Calling AdditiveSchwarz's compute()" << endl;
  prec.compute();

  MV x (rowmap, 2), y (rowmap, 2), z (rowmap, 2);
  x.putScalar (one);

  out << "Applying AdditiveSchwarz to a multivector" << endl;
  prec.apply(x, y);

  // The solution should now be full of 1/2s
  z.putScalar (one / two);

  Teuchos::ArrayRCP<const Scalar> yview = y.get1dView();
  Teuchos::ArrayRCP<const Scalar> zview = z.get1dView();

  TEST_COMPARE_FLOATING_ARRAYS(yview, zview, 4*STS::eps ());
}

TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL(Ifpack2AdditiveSchwarz, FastILDL, Scalar, LocalOrdinal, GlobalOrdinal)
{
  typedef LocalOrdinal LO;
  typedef GlobalOrdinal GO;

  typedef Tpetra::CrsMatrix<Scalar,LO,GO,Node> crs_matrix_type;
  typedef Tpetra::Map<LO,GO,Node> map_type;
  typedef Tpetra::MultiVector<Scalar,LO,GO,Node> MV;
  typedef Tpetra::RowMatrix<Scalar,LO,GO,Node> row_matrix_type;
  typedef Teuchos::ScalarTraits<Scalar> STS;
  const Scalar one = STS::one ();
  const Scalar two = STS::one () + STS::one ();

  out << "Ifpack2::AdditiveSchwarz: FastILDL" << endl;
  Teuchos::OSTab tab1 (out);

  global_size_t num_rows_per_proc = 500;

  out << "Creating row Map and CrsMatrix" << endl;

  RCP<const map_type> rowmap =
    tif_utest::create_tpetra_map<LO, GO, Node> (num_rows_per_proc);
  RCP<const crs_matrix_type> crsmatrix =
    tif_utest::create_test_matrix<Scalar, LO, GO, Node> (rowmap);

  out << "Creating AdditiveSchwarz instance" << endl;

  Ifpack2::AdditiveSchwarz<row_matrix_type> prec (crsmatrix);
  Teuchos::ParameterList params, zlist;

  out << "Filling in ParameterList for AdditiveSchwarz" << endl;

  zlist.set ("order_method", "rcm");
  zlist.set ("order_method_type", "local");
  params.set ("inner preconditioner name", "FAST_ILDL");
  {
    Teuchos::ParameterList innerParams;
    params.set ("inner preconditioner parameters", innerParams);
  }
  params.set ("schwarz: overlap level", static_cast<int> (0));
  params.set ("schwarz: combine mode", "Zero");

#if defined(HAVE_IFPACK2_XPETRA) && defined(HAVE_IFPACK2_ZOLTAN2)
  params.set ("schwarz: use reordering", true);
  params.set ("schwarz: reordering list", zlist);
#else
  params.set ("schwarz: use reordering", false);
#endif

  out << "Setting AdditiveSchwarz's parameters" << endl;

  TEST_NOTHROW(prec.setParameters(params));

  out << "Testing domain and range Maps of AdditiveSchwarz" << endl;

  // FIXME (mfh 26 Jul 2015) The domain and range Maps of the
  // preconditioner don't have to be the same object; they only need
  // to be the same in the sense of Tpetra::Map::isSameAs().
  //
  //trivial tests to insist that the preconditioner's domain/range maps are
  //identically those of the matrix:
  const map_type* mtx_dom_map_ptr = &*crsmatrix->getDomainMap();
  const map_type* mtx_rng_map_ptr = &*crsmatrix->getRangeMap();
  const map_type* prec_dom_map_ptr = &*prec.getDomainMap();
  const map_type* prec_rng_map_ptr = &*prec.getRangeMap();
  TEST_EQUALITY( prec_dom_map_ptr, mtx_dom_map_ptr );
  TEST_EQUALITY( prec_rng_map_ptr, mtx_rng_map_ptr );

  out << "Calling AdditiveSchwarz's initialize()" << endl;
  prec.initialize();

  out << "Calling AdditiveSchwarz's compute()" << endl;
  prec.compute();

  MV x (rowmap, 2), y (rowmap, 2), z (rowmap, 2);
  x.putScalar (one);

  out << "Applying AdditiveSchwarz to a multivector" << endl;
  prec.apply(x, y);

  // The solution should now be full of 1/2s
  z.putScalar (one / two);

  Teuchos::ArrayRCP<const Scalar> yview = y.get1dView();
  Teuchos::ArrayRCP<const Scalar> zview = z.get1dView();

  TEST_COMPARE_FLOATING_ARRAYS(yview, zview, 4*STS::eps ());
}

#if defined(HAVE_IFPACK2_AMESOS2) and defined(HAVE_IFPACK2_XPETRA) and (defined(HAVE_AMESOS2_SUPERLU) || defined(HAVE_AMESOS2_KLU2))

#  define IFPACK2_AMESOS2_SUPERLU_SCALAR_ORDINAL(Scalar,LocalOrdinal,GlobalOrdinal) \
     TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( Ifpack2AdditiveSchwarz, SparseDirectSolver, Scalar, LocalOrdinal, GlobalOrdinal)
#else
#  define IFPACK2_AMESOS2_SUPERLU_SCALAR_ORDINAL(Scalar,LocalOrdinal,GlobalOrdinal)
#endif

#if defined(HAVE_IFPACK2_SHYLU_NODEFASTILU)
  #define IFPACK2_FASTILU_SCALAR_ORDINAL(Scalar, LocalOrdinal, GlobalOrdinal) \
    TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( Ifpack2AdditiveSchwarz, FastILU, Scalar, LocalOrdinal, GlobalOrdinal) \
    TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( Ifpack2AdditiveSchwarz, FastIC, Scalar, LocalOrdinal, GlobalOrdinal) \
    TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( Ifpack2AdditiveSchwarz, FastILDL, Scalar, LocalOrdinal, GlobalOrdinal)
#else
  #define IFPACK2_FASTILU_SCALAR_ORDINAL(Scalar, LocalOrdinal, GlobalOrdinal)
#endif

#  define UNIT_TEST_GROUP_SCALAR_ORDINAL(Scalar,LocalOrdinal,GlobalOrdinal) \
     TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( Ifpack2AdditiveSchwarz, Test0, Scalar, LocalOrdinal,GlobalOrdinal)  \
     TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( Ifpack2AdditiveSchwarz, Test1, Scalar, LocalOrdinal,GlobalOrdinal)  \
     TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( Ifpack2AdditiveSchwarz, Test2, Scalar, LocalOrdinal,GlobalOrdinal)  \
     TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( Ifpack2AdditiveSchwarz, RILUK, Scalar, LocalOrdinal,GlobalOrdinal) \
     TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( Ifpack2AdditiveSchwarz, RILUK_UserOrdering, Scalar, LocalOrdinal,GlobalOrdinal) \
     TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( Ifpack2AdditiveSchwarz, TestOverlap, Scalar, LocalOrdinal, GlobalOrdinal) \
     TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( Ifpack2AdditiveSchwarz, MultipleSweeps, Scalar, LocalOrdinal, GlobalOrdinal) \
     TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( Ifpack2AdditiveSchwarz, ILU_Overlap, Scalar, LocalOrdinal, GlobalOrdinal) \
     TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( Ifpack2AdditiveSchwarz, ILU_NonOverlap, Scalar, LocalOrdinal, GlobalOrdinal) \
     TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( Ifpack2AdditiveSchwarz, SGS_NonOverlap, Scalar, LocalOrdinal, GlobalOrdinal) \
     TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( Ifpack2AdditiveSchwarz, SGS_Overlap, Scalar, LocalOrdinal, GlobalOrdinal) \
     IFPACK2_AMESOS2_SUPERLU_SCALAR_ORDINAL(Scalar, LocalOrdinal, GlobalOrdinal) \
     IFPACK2_FASTILU_SCALAR_ORDINAL(Scalar, LocalOrdinal, GlobalOrdinal)

// mfh 26 Aug 2015: Ifpack2::AdditiveSchwarz was only getting tested
// for Scalar = double, LocalOrdinal = int, GlobalOrdinal = int, and
// the default Node type.  As part of the fix for Bug 6358, I'm
// removing the assumption that GlobalOrdinal = int exists.

typedef Tpetra::MultiVector<>::scalar_type default_scalar_type;
typedef Tpetra::MultiVector<>::local_ordinal_type default_local_ordinal_type;
typedef Tpetra::MultiVector<>::global_ordinal_type default_global_ordinal_type;

UNIT_TEST_GROUP_SCALAR_ORDINAL(default_scalar_type, default_local_ordinal_type, default_global_ordinal_type)

} // namespace (anonymous)

