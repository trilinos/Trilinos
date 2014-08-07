/*
//@HEADER
// ***********************************************************************
//
//       Ifpack2: Tempated Object-Oriented Algebraic Preconditioner Package
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

// ***********************************************************************
//
//      Ifpack2: Tempated Object-Oriented Algebraic Preconditioner Package
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
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************


/*! \file Ifpack2_UnitTestAdditiveSchwarz.cpp

\brief Ifpack2 Unit test for AdditiveSchwarz.
*/

#include <iostream>
#include <Teuchos_ConfigDefs.hpp>

#include <Ifpack2_ConfigDefs.hpp>
#include <Teuchos_UnitTestHarness.hpp>
#include <Ifpack2_Version.hpp>

#if defined(HAVE_IFPACK2_QD) && !defined(HAVE_TPETRA_EXPLICIT_INSTANTIATION)
#include <qd/dd_real.h>
#endif

// Xpetra / Galeri
#ifdef HAVE_IFPACK2_XPETRA
#include <Xpetra_ConfigDefs.hpp>
#include <Xpetra_DefaultPlatform.hpp>
#include <Xpetra_Parameters.hpp>
#include <Xpetra_MapFactory.hpp>
#include <Xpetra_TpetraMap.hpp>
#include <Xpetra_CrsMatrix.hpp>
#include <Xpetra_TpetraCrsMatrix.hpp>
#include <Galeri_XpetraProblemFactory.hpp>
#include <Galeri_XpetraMatrixTypes.hpp>
#include <Galeri_XpetraParameters.hpp>
#include <Galeri_XpetraUtils.hpp>
#include <Galeri_XpetraMaps.hpp>
#endif

#include <Ifpack2_UnitTestHelpers.hpp>
#include <Ifpack2_AdditiveSchwarz.hpp>

namespace {
using Tpetra::global_size_t;
typedef tif_utest::Node Node;

//this macro declares the unit-test-class:
TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL(Ifpack2AdditiveSchwarz, Test0, Scalar, LocalOrdinal, GlobalOrdinal)
{
  using std::endl;

//we are now in a class method declared by the above macro, and
//that method has these input arguments:
//Teuchos::FancyOStream& out, bool& success

  typedef Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> CrsType;

  std::string version = Ifpack2::Version();
  out << "Ifpack2::Version(): " << version << endl;

  global_size_t num_rows_per_proc = 5;

  out << "Creating row Map and CrsMatrix" << endl;

  const Teuchos::RCP<const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> > rowmap = tif_utest::create_tpetra_map<LocalOrdinal,GlobalOrdinal,Node>(num_rows_per_proc);
  Teuchos::RCP<const Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> > crsmatrix = tif_utest::create_test_matrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>(rowmap);

  out << "Creating AdditiveSchwarz instance" << endl;

  Ifpack2::AdditiveSchwarz<CrsType> prec (crsmatrix);
  Teuchos::ParameterList params, zlist;

  out << "Filling in ParameterList for AdditiveSchwarz" << endl;

  zlist.set ("order_method", "rcm");
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

  //trivial tests to insist that the preconditioner's domain/range maps are
  //identically those of the matrix:
  const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node>* mtx_dom_map_ptr = &*crsmatrix->getDomainMap();
  const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node>* mtx_rng_map_ptr = &*crsmatrix->getRangeMap();
  const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node>* prec_dom_map_ptr = &*prec.getDomainMap();
  const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node>* prec_rng_map_ptr = &*prec.getRangeMap();
  TEST_EQUALITY( prec_dom_map_ptr, mtx_dom_map_ptr );
  TEST_EQUALITY( prec_rng_map_ptr, mtx_rng_map_ptr );

  out << "Calling AdditiveSchwarz's initialize()" << endl;
  prec.initialize();

  out << "Calling AdditiveSchwarz's compute()" << endl;
  prec.compute();

  Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> x(rowmap,2), y(rowmap,2), z(rowmap,2);
  x.putScalar(1);

  out << "Applying AdditiveSchwarz to a multivector" << endl;
  prec.apply(x, y);

  // The solution should now be full of 1/2s
  z.putScalar(0.5);

  Teuchos::ArrayRCP<const Scalar> yview = y.get1dView();
  Teuchos::ArrayRCP<const Scalar> zview = z.get1dView();

  TEST_COMPARE_FLOATING_ARRAYS(yview, zview, 4*Teuchos::ScalarTraits<Scalar>::eps());
}

// ///////////////////////////////////////////////////////////////////// //

TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL(Ifpack2AdditiveSchwarz, Test1, Scalar, LocalOrdinal, GlobalOrdinal)
{
//we are now in a class method declared by the above macro, and
//that method has these input arguments:
//Teuchos::FancyOStream& out, bool& success

  typedef Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> CrsType;

  std::string version = Ifpack2::Version();
  out << "Ifpack2::Version(): " << version << std::endl;

  global_size_t num_rows_per_proc = 5;

  const Teuchos::RCP<const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> > rowmap = tif_utest::create_tpetra_map<LocalOrdinal,GlobalOrdinal,Node>(num_rows_per_proc);

  // Don't run in serial.
  if(rowmap->getComm()->getSize()==1) return;

  Teuchos::RCP<const Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> > crsmatrix = tif_utest::create_test_matrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>(rowmap);

  Ifpack2::AdditiveSchwarz<CrsType> prec (crsmatrix);

  Teuchos::ParameterList params, zlist;
  zlist.set ("order_method", "rcm");

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
  const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node>* mtx_dom_map_ptr = &*crsmatrix->getDomainMap();
  const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node>* mtx_rng_map_ptr = &*crsmatrix->getRangeMap();
  const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node>* prec_dom_map_ptr = &*prec.getDomainMap();
  const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node>* prec_rng_map_ptr = &*prec.getRangeMap();
  TEST_EQUALITY( prec_dom_map_ptr, mtx_dom_map_ptr );
  TEST_EQUALITY( prec_rng_map_ptr, mtx_rng_map_ptr );

  prec.initialize();

  prec.compute();

  //prec.describe (* Teuchos::getFancyOStream (Teuchos::rcpFromRef (std::cerr)), Teuchos::VERB_EXTREME);

  Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> x(rowmap,1), y(rowmap,1), z(rowmap,1);
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
    for (GlobalOrdinal i=0; i<overlapLevel; ++i)
      zdata[i] += 0.5;
  }
  else {
    for (GlobalOrdinal i=0; i<overlapLevel; ++i)
      zdata[i] += 0.5;
    for (GlobalOrdinal i=0; i<overlapLevel; ++i)
      zdata[num_rows_per_proc-i-1] += 0.5;
  }
  zdata = Teuchos::null;


  Teuchos::ArrayRCP<const Scalar> yview = y.get1dView();
  Teuchos::ArrayRCP<const Scalar> zview = z.get1dView();

  TEST_COMPARE_FLOATING_ARRAYS(yview, zview, 4*Teuchos::ScalarTraits<Scalar>::eps());
}

// ///////////////////////////////////////////////////////////////////// //

TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL(Ifpack2AdditiveSchwarz, Test2, Scalar, LocalOrdinal, GlobalOrdinal)
{

  typedef Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> CrsType;

  std::string version = Ifpack2::Version();
  out << "Ifpack2::Version(): " << version << std::endl;

  global_size_t num_rows_per_proc = 5;

  const Teuchos::RCP<const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> > rowmap = tif_utest::create_tpetra_map<LocalOrdinal,GlobalOrdinal,Node>(num_rows_per_proc);

  // Don't run in serial.
  if(rowmap->getComm()->getSize()==1) {
    out << std::endl << "This test must be run on more than one process." << std::endl << std::endl;
    return;
  }

  Teuchos::RCP<const Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> > crsmatrix = tif_utest::create_test_matrix3<Scalar,LocalOrdinal,GlobalOrdinal,Node>(rowmap);

  Ifpack2::AdditiveSchwarz<CrsType> prec (crsmatrix);
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


  Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> x(rowmap,1), y(rowmap,1), z(rowmap,1);
  x.putScalar(1);
  crsmatrix->apply(x,z);
  prec.apply(z, y);

  // The solution should now be full of 1s
  z.putScalar(1.0);

  Teuchos::ArrayRCP<const Scalar> yview = y.get1dView();
  Teuchos::ArrayRCP<const Scalar> zview = z.get1dView();

  TEST_COMPARE_FLOATING_ARRAYS(yview, zview, 4*Teuchos::ScalarTraits<Scalar>::eps());
}


// Test RILUK as subdomain solver for AdditiveSchwarz.
TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL(Ifpack2AdditiveSchwarz, Test3, Scalar, LocalOrdinal, GlobalOrdinal)
{
  using std::endl;

//we are now in a class method declared by the above macro, and
//that method has these input arguments:
//Teuchos::FancyOStream& out, bool& success

  typedef Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> CrsType;

  std::string version = Ifpack2::Version();
  out << "Ifpack2::Version(): " << version << endl;

  global_size_t num_rows_per_proc = 5;

  out << "Creating row Map and CrsMatrix" << endl;

  Teuchos::RCP<const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> > rowmap =
    tif_utest::create_tpetra_map<LocalOrdinal,GlobalOrdinal,Node>(num_rows_per_proc);
  Teuchos::RCP<const Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> > crsmatrix =
    tif_utest::create_test_matrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>(rowmap);

  out << "Creating AdditiveSchwarz instance" << endl;

  Ifpack2::AdditiveSchwarz<CrsType> prec (crsmatrix);
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
  const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node>* mtx_dom_map_ptr = &*crsmatrix->getDomainMap();
  const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node>* mtx_rng_map_ptr = &*crsmatrix->getRangeMap();
  const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node>* prec_dom_map_ptr = &*prec.getDomainMap();
  const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node>* prec_rng_map_ptr = &*prec.getRangeMap();
  TEST_EQUALITY( prec_dom_map_ptr, mtx_dom_map_ptr );
  TEST_EQUALITY( prec_rng_map_ptr, mtx_rng_map_ptr );

  out << "Calling AdditiveSchwarz's initialize()" << endl;
  prec.initialize();

  out << "Calling AdditiveSchwarz's compute()" << endl;
  prec.compute();

  Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> x(rowmap,2), y(rowmap,2), z(rowmap,2);
  x.putScalar(1);

  out << "Applying AdditiveSchwarz to a multivector" << endl;
  prec.apply (x, y);

  out << "Testing result of AdditiveSchwarz's apply" << endl;

  // The solution should now be full of 1/2s
  z.putScalar(0.5);

  Teuchos::ArrayRCP<const Scalar> yview = y.get1dView();
  Teuchos::ArrayRCP<const Scalar> zview = z.get1dView();

  TEST_COMPARE_FLOATING_ARRAYS(yview, zview, 4*Teuchos::ScalarTraits<Scalar>::eps());
}


// ///////////////////////////////////////////////////////////////////// //

TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL(Ifpack2AdditiveSchwarz, TestGIDs, Scalar, LocalOrdinal, GlobalOrdinal)
{
  // Test that AdditiveSchwarz transfer patterns are correct.
  // A vector v such that v(i) = GID(i) is passed in to AS.  Using the IdentitySolver with any amount
  // of overlap, any subdomain reordering, and combine mode Zero, the solution should be v.

  typedef Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> CrsType;

  std::string version = Ifpack2::Version();
  out << "Ifpack2::Version(): " << version << std::endl;

  global_size_t num_rows_per_proc = 10;
  const Teuchos::RCP<const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> > rowmap = tif_utest::create_tpetra_map<LocalOrdinal,GlobalOrdinal,Node>(num_rows_per_proc);

  // Don't run in serial.
  if(rowmap->getComm()->getSize()==1) {
    out << std::endl << "This test must be run on more than one process." << std::endl << std::endl;
    return;
  }

  Teuchos::RCP<const Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> > crsmatrix = tif_utest::create_test_matrix2<Scalar,LocalOrdinal,GlobalOrdinal,Node>(rowmap);


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

      Ifpack2::AdditiveSchwarz<CrsType> prec (crsmatrix);

      Teuchos::ParameterList params, zlist;
      params.set ("schwarz: overlap level", overlapLevel);
      params.set ("schwarz: use reordering",reorderSubdomains);
      zlist.set ("order_method", "rcm");
      params.set ("schwarz: reordering list", zlist);
      params.set("schwarz: combine mode", "Zero");
      params.set ("inner preconditioner name", "IDENTITY");

      prec.setParameters(params);

      prec.initialize();
      prec.compute();


      Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> x(rowmap,1), y(rowmap,1);
      Teuchos::ArrayRCP<Scalar> xData = x.getDataNonConst(0);
      for (GlobalOrdinal j=0; j<xData.size(); ++j)
        xData[j] = rowmap->getGlobalElement(j);
      xData = Teuchos::null;
      prec.apply(x, y);

      out << prec.description() << std::endl;

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

#if defined(HAVE_IFPACK2_AMESOS2) and defined(HAVE_IFPACK2_XPETRA) and defined(HAVE_AMESOS2_SUPERLU)
// Test SuperLU sparse direct solver as subdomain solver for AdditiveSchwarz.
TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL(Ifpack2AdditiveSchwarz, SuperLU, Scalar, LocalOrdinal, GlobalOrdinal)
{
  typedef Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>       CrsType;
  typedef Xpetra::TpetraCrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> XCrsType;
  typedef Xpetra::Map<LocalOrdinal,GlobalOrdinal,Node>                    XMapType;
  typedef Xpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>     XMVectorType;
  typedef Tpetra::Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node>          MultiVectorType;
  typedef Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node>                    MapType;

  using Teuchos::RCP;
  using Teuchos::ParameterList;

  std::string version = Ifpack2::Version();
  out << "Ifpack2::Version(): " << version << std::endl;

  // Generate the matrix using Galeri.
  RCP<const Teuchos::Comm<int> > comm = Tpetra::DefaultPlatform::getDefaultPlatform().getComm();
  Teuchos::CommandLineProcessor clp;
  GlobalOrdinal nx = 10, ny=10, nz=10;
  Galeri::Xpetra::Parameters<GlobalOrdinal> GaleriParameters(clp, nx, ny, nz, "Laplace2D");
  Xpetra::Parameters xpetraParameters(clp);
  ParameterList GaleriList = GaleriParameters.GetParameterList();

  RCP<XMapType > xmap = Galeri::Xpetra::CreateMap<LocalOrdinal, GlobalOrdinal, Node>(xpetraParameters.GetLib(), "Cartesian2D", comm, GaleriList);
  RCP<Galeri::Xpetra::Problem<XMapType,XCrsType,XMVectorType> > Pr = Galeri::Xpetra::BuildProblem<Scalar,LocalOrdinal,GlobalOrdinal,XMapType,XCrsType,XMVectorType>
      (string("Laplace2D"),xmap,GaleriList);

  RCP<XCrsType> XA = Pr->BuildMatrix();
  RCP<CrsType> A = XA->getTpetra_CrsMatrixNonConst();
  TEST_INEQUALITY(A,Teuchos::null);
  
  RCP<const MapType > rowmap = A->getRowMap();

  out << "Creating AdditiveSchwarz instance" << std::endl;

  Ifpack2::AdditiveSchwarz<CrsType> prec (A,1);
  ParameterList params, zlist;

  out << "Filling in ParameterList for AdditiveSchwarz" << std::endl;

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
  subdomainList.set("Amesos2 solver name","superlu");
  ParameterList &amesos2List = subdomainList.sublist("Amesos2");
  ParameterList &superluList = amesos2List.sublist("SuperLU");
  superluList.set("ILU_Flag",false);

  out << "Setting AdditiveSchwarz's parameters" << std::endl;

  std::ostringstream ps;
  int indent = 4;
  params.print(ps, indent);
  out << ps.str() << std::endl;

  TEST_NOTHROW(prec.setParameters(params));

  out << "Testing domain and range Maps of AdditiveSchwarz" << std::endl;

  //trivial tests to insist that the preconditioner's domain/range maps are
  //identically those of the matrix:
  const MapType* mtx_dom_map_ptr = &*A->getDomainMap();
  const MapType* mtx_rng_map_ptr = &*A->getRangeMap();
  const MapType* prec_dom_map_ptr = &*prec.getDomainMap();
  const MapType* prec_rng_map_ptr = &*prec.getRangeMap();
  TEST_EQUALITY( prec_dom_map_ptr, mtx_dom_map_ptr );
  TEST_EQUALITY( prec_rng_map_ptr, mtx_rng_map_ptr );

  out << std::endl << "solve using AdditiveSchwarz with sparse direct method as the subdomain solve" << std::endl;
  out << "Calling AdditiveSchwarz::initialize()" << std::endl;
  prec.initialize();

  out << "Calling AdditiveSchwarz::compute()" << std::endl;
  prec.compute();

  MultiVectorType x(rowmap,2), y(rowmap,2);
  x.randomize();
  out << "Calling AdditiveSchwarz::apply()" << std::endl;
  prec.apply (x, y);

  // Now switch to dense direct solves on the subdomains.
  out << std::endl << "solve using AdditiveSchwarz with dense direct method as the subdomain solve" << std::endl;
  params.set ("subdomain solver name", "DENSE");
  prec.setParameters(params);
  out << "Calling AdditiveSchwarz::initialize()" << std::endl;
  prec.initialize();
  out << "Calling AdditiveSchwarz::compute()" << std::endl;
  prec.compute();
  out << "Calling AdditiveSchwarz::apply()" << std::endl;
  MultiVectorType z(rowmap,2);
  prec.apply (x, z);

  out << "Comparing results of two solves" << std::endl;
  Teuchos::Array<typename Teuchos::ScalarTraits<Scalar>::magnitudeType> ynorms(1),znorms(1);
  y.norm2(ynorms());
  z.norm2(znorms());
  out << "solution norm, sparse direct solve: " << std::setprecision(7) << ynorms[0] << std::endl;
  out << "solution norm,  dense direct solve: " << std::setprecision(7) << znorms[0] << std::endl;
  TEST_FLOATING_EQUALITY(ynorms[0], znorms[0], 10*Teuchos::ScalarTraits<Scalar>::eps());

}
#endif

#if defined(HAVE_IFPACK2_AMESOS2) and defined(HAVE_IFPACK2_XPETRA) and defined(HAVE_AMESOS2_SUPERLU)

#  define UNIT_TEST_GROUP_SCALAR_ORDINAL(Scalar,LocalOrdinal,GlobalOrdinal) \
     TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( Ifpack2AdditiveSchwarz, Test0, Scalar, LocalOrdinal,GlobalOrdinal)  \
     TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( Ifpack2AdditiveSchwarz, Test1, Scalar, LocalOrdinal,GlobalOrdinal)  \
     TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( Ifpack2AdditiveSchwarz, Test2, Scalar, LocalOrdinal,GlobalOrdinal) \
     TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( Ifpack2AdditiveSchwarz, TestGIDs, Scalar, LocalOrdinal, GlobalOrdinal) \
     TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( Ifpack2AdditiveSchwarz, SuperLU, Scalar, LocalOrdinal, GlobalOrdinal)

#else

#  define UNIT_TEST_GROUP_SCALAR_ORDINAL(Scalar,LocalOrdinal,GlobalOrdinal) \
     TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( Ifpack2AdditiveSchwarz, Test0, Scalar, LocalOrdinal,GlobalOrdinal)  \
     TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( Ifpack2AdditiveSchwarz, Test1, Scalar, LocalOrdinal,GlobalOrdinal)  \
     TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( Ifpack2AdditiveSchwarz, Test2, Scalar, LocalOrdinal,GlobalOrdinal) \
     TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( Ifpack2AdditiveSchwarz, TestGIDs, Scalar, LocalOrdinal, GlobalOrdinal)

#endif

//TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( Ifpack2AdditiveSchwarz, TestGIDs, double, int, int)

UNIT_TEST_GROUP_SCALAR_ORDINAL(double, int, int)
#ifndef HAVE_IFPACK2_EXPLICIT_INSTANTIATION
UNIT_TEST_GROUP_SCALAR_ORDINAL(float, short, int)
#endif

#if defined(HAVE_IFPACK2_QD) && !defined(HAVE_TPETRA_EXPLICIT_INSTANTIATION)
UNIT_TEST_GROUP_SCALAR_ORDINAL(dd_real, int, int)
#endif

}//namespace <anonymous>

