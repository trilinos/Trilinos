/*
// @HEADER
// ***********************************************************************
//
//          Tpetra: Templated Linear Algebra Services Package
//                 Copyright (2008) Sandia Corporation
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ************************************************************************
// @HEADER
*/

#include "TpetraCore_config.h"
#ifdef HAVE_TPETRACORE_MPI

#include "Tpetra_TestingUtilities.hpp"
#include "Tpetra_Details_extractMpiCommFromTeuchos.hpp"
#include "Teuchos_DefaultMpiComm.hpp"
#include "Teuchos_DefaultSerialComm.hpp"

namespace { // (anonymous)

  using Tpetra::Details::extractMpiCommFromTeuchos;
  using Teuchos::Comm;
  using Teuchos::MpiComm;
  using Teuchos::RCP;
  using Teuchos::rcp_dynamic_cast;
  using Teuchos::SerialComm;
  using std::endl;

  //
  // UNIT TESTS
  //

  TEUCHOS_UNIT_TEST( extractMpiCommFromTeuchos, basic )
  {
    out << "Testing Tpetra::Details::extractMpiCommFromTeuchos" << endl;
    Teuchos::OSTab tab1 (out);

    int myRank;
    (void) MPI_Comm_rank (MPI_COMM_WORLD, &myRank);
    int numProcs;
    (void) MPI_Comm_size (MPI_COMM_WORLD, &numProcs);

    out << "Test MPI_COMM_WORLD wrapper" << endl;
    // not in an inner scope, because we want to use it later
    RCP<const MpiComm<int> > mpiCommWorld (new MpiComm<int> (MPI_COMM_WORLD));
    {
      TEST_ASSERT( ! mpiCommWorld->getRawMpiComm ().is_null () );
      if (! mpiCommWorld->getRawMpiComm ().is_null ()) {
        TEST_EQUALITY( myRank, mpiCommWorld->getRank () );
        TEST_EQUALITY( numProcs, mpiCommWorld->getSize () );

        MPI_Comm mpiComm = * (mpiCommWorld->getRawMpiComm ());
        int result;
        (void) MPI_Comm_compare (mpiComm, MPI_COMM_WORLD, &result);
        TEST_EQUALITY( result, MPI_IDENT );

        MPI_Comm mpiComm2 = extractMpiCommFromTeuchos (*mpiCommWorld);
        (void) MPI_Comm_compare (mpiComm2, MPI_COMM_WORLD, &result);
        TEST_EQUALITY( result, MPI_IDENT );
      }
    }

    out << "Test MPI_COMM_SELF wrapper" << endl;
    {
      RCP<const MpiComm<int> > mpiCommSelf (new MpiComm<int> (MPI_COMM_SELF));
      TEST_ASSERT( ! mpiCommSelf->getRawMpiComm ().is_null () );
      if (! mpiCommSelf->getRawMpiComm ().is_null ()) {
        MPI_Comm mpiComm = * (mpiCommSelf->getRawMpiComm ());

        int result;
        (void) MPI_Comm_compare (mpiComm, MPI_COMM_SELF, &result);
        TEST_EQUALITY( result, MPI_IDENT );

        MPI_Comm mpiComm2 = extractMpiCommFromTeuchos (*mpiCommSelf);
        (void) MPI_Comm_compare (mpiComm2, MPI_COMM_SELF, &result);
        TEST_EQUALITY( result, MPI_IDENT );
      }
    }

    out << "Test Teuchos::SerialComm" << endl;
    {
      RCP<const SerialComm<int> > teuchosSerialComm (new SerialComm<int> ());
      MPI_Comm serialComm = extractMpiCommFromTeuchos (*teuchosSerialComm);

      int result;
      (void) MPI_Comm_compare (serialComm, MPI_COMM_SELF, &result);
      TEST_EQUALITY( result, MPI_IDENT );
    }

    if (mpiCommWorld->getSize () > 1) {
      out << "Test result of splitting MPI_COMM_WORLD" << endl;

      const int color = (myRank % 2);
      const int key = 0;
      MPI_Comm newMpiComm;
      (void) MPI_Comm_split (MPI_COMM_WORLD, color, key, &newMpiComm);

      RCP<const Comm<int> > newTeuchosComm = mpiCommWorld->split (color, key);
      TEST_ASSERT( ! newTeuchosComm.is_null () );
      if (! newTeuchosComm.is_null ()) {
        int r, p;
        (void) MPI_Comm_rank (newMpiComm, &r);
        (void) MPI_Comm_size (newMpiComm, &p);
        TEST_EQUALITY( r, newTeuchosComm->getRank () );
        TEST_EQUALITY( p, newTeuchosComm->getSize () );

        RCP<const MpiComm<int> > newTeuchosMpiComm =
          rcp_dynamic_cast<const MpiComm<int> > (newTeuchosComm);
        TEST_ASSERT( ! newTeuchosMpiComm.is_null () );

        if (! newTeuchosMpiComm.is_null ()) {
          TEST_ASSERT( ! newTeuchosMpiComm->getRawMpiComm ().is_null () );
          if (! newTeuchosMpiComm->getRawMpiComm ().is_null ()) {
            MPI_Comm newTeuchosMpiCommRaw = * (newTeuchosMpiComm->getRawMpiComm ());

            // We called split twice, so we won't necessarily get the
            // same (MPI_IDENT) communicators.  We must, however, get
            // congruent (MPI_CONGRUENT) communicators.

            int result;
            (void) MPI_Comm_compare (newMpiComm, newTeuchosMpiCommRaw, &result);
            TEST_EQUALITY( result, MPI_CONGRUENT );

            MPI_Comm mpiComm2 = extractMpiCommFromTeuchos (*newTeuchosComm);
            (void) MPI_Comm_compare (mpiComm2, newMpiComm, &result);
            TEST_EQUALITY( result, MPI_CONGRUENT );
          }
        }
      }
    }

    int lclSuccess = success ? 1 : 0;
    int gblSuccess = 0; // output argument
    (void) MPI_Allreduce (&lclSuccess, &gblSuccess, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
    TEST_EQUALITY( gblSuccess, 1 );
  }

} // namespace (anonymous)

#endif // HAVE_TPETRACORE_MPI
