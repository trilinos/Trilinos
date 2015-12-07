// @HEADER
//
// ***********************************************************************
//
//             Xpetra: A linear algebra interface package
//                  Copyright 2012 Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
// Questions? Contact
//                    Jonathan Hu       (jhu@sandia.gov)
//                    Andrey Prokopenko (aprokop@sandia.gov)
//                    Ray Tuminaro      (rstumin@sandia.gov)
//
// ***********************************************************************
//
// @HEADER
#include <Teuchos_UnitTestHarness.hpp>
#include <Xpetra_UnitTestHelpers.hpp>
#include <Teuchos_Comm.hpp>

#include "Xpetra_ConfigDefs.hpp"
#include "Xpetra_DefaultPlatform.hpp"

#include "Xpetra_MultiVector.hpp"
#include "Xpetra_Vector.hpp"

#include "Xpetra_MapFactory.hpp"
#include "Xpetra_MultiVectorFactory.hpp"

#ifdef HAVE_XPETRA_TPETRA
#include "Xpetra_TpetraMultiVector.hpp"
#include "Xpetra_TpetraVector.hpp"
#endif

#ifdef HAVE_XPETRA_EPETRA
#include "Xpetra_EpetraMultiVector.hpp"
#include "Xpetra_EpetraVector.hpp"
#endif

namespace {

  TEUCHOS_STATIC_SETUP()
  {
    // Teuchos::CommandLineProcessor &clp = Teuchos::UnitTestRepository::getCLP();
    // clp.addOutputSetupOptions(true);
    // clp.setOption(
    //               "test-mpi", "test-serial", &testMpi,
    //               "Test MPI (if available) or force test of serial.  In a serial build,"
    //               " this option is ignored and a serial comm is always used." );
  }

  // Test getVector() / getVectorNonConst()
  // More specifically, this test verifies that the newly created vector will remain valid after the disappearance of the references to the multivector in user code.
  void XpetraSpecific_GetVector(Xpetra::UnderlyingLib lib, Teuchos::FancyOStream & out, bool & success) {

    using Teuchos::RCP;
    using Teuchos::rcp;

    RCP<const Teuchos::Comm<int> > comm = Xpetra::DefaultPlatform::getDefaultPlatform().getComm();
    const Xpetra::global_size_t INVALID = Teuchos::OrdinalTraits<Xpetra::global_size_t>::invalid();

    const size_t numLocal = 4;
#ifndef XPETRA_TEST_USE_LONGLONG_GO
	typedef int GO;
#else
	typedef long long GO;
#endif

    RCP<const Xpetra::Map<int, GO> > map = Xpetra::MapFactory<int, GO>::createContigMap(lib, INVALID, numLocal, comm);

    RCP< Xpetra::MultiVector<double, int, GO> > mv = Xpetra::MultiVectorFactory<double, int, GO>::Build(map, 3, false);
    for(size_t k=0; k < 3; k++) {
      Teuchos::ArrayRCP<double> mvData = mv->getDataNonConst(k);

      for(size_t i=0; i < numLocal; i++) {
        mvData[i] = i*(k+1);
      }
    }

    Teuchos::RCP< const Xpetra::Vector<double, int, GO> > v         = mv->getVector(1);         // second vector
    Teuchos::RCP<       Xpetra::Vector<double, int, GO> > vNonConst = mv->getVectorNonConst(2); // third vector

    mv = Teuchos::null;

    {
      Teuchos::ArrayRCP<const double> vData = v->getData(0);
      for(size_t i=0; i< numLocal; i++) {
        TEST_EQUALITY(vData[i], i*2);
      }
    }

    {
      Teuchos::ArrayRCP<double> vData = vNonConst->getDataNonConst(0);
      for(size_t i=0; i< numLocal; i++) {
        TEST_EQUALITY(vData[i], i*3);
      }
    }

  }


  // Test getVector() / getVectorNonConst()
  TEUCHOS_UNIT_TEST(MultiVector, XpetraSpecific_GetVector)
  {
#ifdef HAVE_XPETRA_TPETRA
    XpetraSpecific_GetVector(Xpetra::UseTpetra, out, success);
#endif
#ifdef HAVE_XPETRA_EPETRA
    //XpetraSpecific_GetVector(Xpetra::UseEpetra, out, success);
#endif
  }

#ifdef HAVE_XPETRA_KOKKOS_REFACTOR
  // Test getLocalMV()
  void XpetraSpecific_GetHostLocalView(Xpetra::UnderlyingLib lib, Teuchos::FancyOStream & out, bool & success) {

    using Teuchos::RCP;
    using Teuchos::rcp;

#ifndef XPETRA_TEST_USE_LONGLONG_GO
    typedef int GO;
#else
    typedef long long GO;
#endif

    typedef Xpetra::MultiVector<double, int, GO> mv_type;
    typedef typename mv_type::dual_view_type dual_view_type;

    RCP<const Teuchos::Comm<int> > comm = Xpetra::DefaultPlatform::getDefaultPlatform().getComm();
    const Xpetra::global_size_t INVALID = Teuchos::OrdinalTraits<Xpetra::global_size_t>::invalid();

    const size_t numLocal = 4;

    RCP<const Xpetra::Map<int, GO> > map = Xpetra::MapFactory<int, GO>::createContigMap(lib, INVALID, numLocal, comm);

    // create new vector and fill it with data
    RCP< mv_type > mv = Xpetra::MultiVectorFactory<double, int, GO>::Build(map, 3, false);
    for(size_t k=0; k < 3; k++) {
      Teuchos::ArrayRCP<double> mvData = mv->getDataNonConst(k);

      for(size_t i=0; i < numLocal; i++) {
        mvData[i] = i*(k+1) + comm->getRank();
      }
    }

    // get a view of the multivector data on the host memory
    dual_view_type::t_host_um hostView = mv->getHostLocalView ();

    TEST_EQUALITY(hostView.dimension_0(), numLocal);
    TEST_EQUALITY(hostView.dimension_1(), 3);
    TEST_EQUALITY(hostView.size(), numLocal * 3);
    for(size_t k=0; k < 3; k++) {
      for(size_t i = 0; i < hostView.dimension_0(); i++) {
        TEST_EQUALITY(hostView(i,k), i*(k+1) + comm->getRank());
      }
    }

    // overwrite data in hostView
    for(size_t r = 0; r < hostView.dimension_0(); r++) {
      for(size_t c = 0; c < hostView.dimension_1(); c++) {
        hostView(r,c) = comm->getRank() + c*hostView.dimension_1() + r + 42.0;
      }
    }

    // check data in multivector

    for(size_t k=0; k < 3; k++) {
      Teuchos::ArrayRCP<const double> vData = mv->getData(k);
      for(size_t i=0; i< numLocal; i++) {
        TEST_EQUALITY(vData[i], comm->getRank() + k*hostView.dimension_1() + i + 42.0);
      }
    }

    // change data in multivector
    for(size_t k=0; k < 3; k++) {
      Teuchos::ArrayRCP<double> vData = mv->getDataNonConst(k);
      for(size_t i=0; i< numLocal; i++) {
        vData[i] = k * numLocal + i;
      }
    }

    // check updated data in view
    for(size_t r = 0; r < hostView.dimension_0(); r++) {
      for(size_t c = 0; c < hostView.dimension_1(); c++) {
        TEST_EQUALITY(hostView(r,c), c * numLocal + r);
      }
    }

    // delete vector
    mv = Teuchos::null;
  }

  // TODO: enable test for different SC, LO, GO, NO...

  // Test getVector() / getVectorNonConst()
  TEUCHOS_UNIT_TEST(MultiVector, XpetraSpecific_GetHostLocalView)
  {
#ifdef HAVE_XPETRA_TPETRA
    XpetraSpecific_GetHostLocalView(Xpetra::UseTpetra, out, success);
#endif
#ifdef HAVE_XPETRA_EPETRA
    //XpetraSpecific_GetHostLocalView(Xpetra::UseEpetra, out, success);
#endif
  }
#endif // HAVE_XPETRA_KOKKOS_REFACTOR
} // namespace
