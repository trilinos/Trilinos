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

//
// Test for Github Issue #114.
//

#include "Tpetra_Core.hpp"
#include "Tpetra_Map.hpp"
#include "Tpetra_Import.hpp"
#include "Tpetra_Vector.hpp"
#include "Tpetra_Details_gathervPrint.hpp"

#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_CommHelpers.hpp"
#include "Teuchos_FancyOStream.hpp"
#include "Teuchos_RCP.hpp"


int
main (int argc, char *argv[])
{
  using Teuchos::outArg;
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::REDUCE_MIN;
  using Teuchos::reduceAll;
  using std::endl;
  using map_type = Tpetra::Map<>;
  using import_type = Tpetra::Import<>;
  using GO = Tpetra::Map<>::global_ordinal_type;
  using GST = Tpetra::global_size_t;

  // mfh 09 Aug 2017: Tpetra instantiates Vector with Scalar=int even
  // if GO=int is disabled, because users (e.g., in MueLu and Zoltan2)
  // want to be able to communicate MPI process ranks.
  using IntVector = Tpetra::Vector<int>;

  Tpetra::ScopeGuard tpetraScope (&argc, &argv);
  auto comm = Tpetra::getDefaultComm ();
  auto outPtr = Teuchos::getFancyOStream (Teuchos::rcpFromRef (std::cout));
  auto& out = *outPtr;
  out.setOutputToRootOnly (0);

  Teuchos::OSTab tab0 (out);
  out << "Test for Github Issue #114" << endl;
  Teuchos::OSTab tab1 (out);

  const GST INVALID = Teuchos::OrdinalTraits<GST>::invalid ();
  const GO indexBase = 0;
  const int myRank = comm->getRank ();
  const int numProcs = comm->getSize ();

  if (numProcs != 4) {
    out << "Test FAILED; must be run on exactly 4 MPI processes, "
      "but numProcs = " << numProcs << endl;
    return EXIT_FAILURE;
  }

  // Problem setup
  int num_per_proc;
  if (myRank == 0) {
    num_per_proc=7;
  }
  else {
    num_per_proc=6;
  }

  GO from_gids_p0[7] = {0,1,2,3,4,5,6};
  GO to_gids_p0[7]   = {0,4,8,12,16,20,24};

  GO from_gids_p1[6] = {7,8,9,10,11,12};
  GO to_gids_p1[6]   = {1,5,9,13,17,21};

  GO from_gids_p2[6] = {13,14,15,16,17,18};
  GO to_gids_p2[6]   = {2,6,10,14,18,22};

  GO from_gids_p3[6] = {19,20,21,22,23,24};
  GO to_gids_p3[6]   = {3,7,11,15,19,23};

  // Correctness check array
  int who_owns[25];
  for (int i = 0; i < 7; ++i) {
    who_owns[to_gids_p0[i]] = 0;
  }
  for (int i = 0; i < 6; ++i) {
    who_owns[to_gids_p1[i]] = 1;
    who_owns[to_gids_p2[i]] = 2;
    who_owns[to_gids_p3[i]] = 3;
  }

  GO* from_ptr = NULL;
  GO* to_ptr = NULL;
  if (myRank == 0) {
    from_ptr = &from_gids_p0[0];
    to_ptr = &to_gids_p0[0];
  }
  else if (myRank == 1) {
    from_ptr = &from_gids_p1[0];
    to_ptr=&to_gids_p1[0];
  }
  else if (myRank == 2) {
    from_ptr = &from_gids_p2[0];
    to_ptr=&to_gids_p2[0];
  }
  else if (myRank==3) {
    from_ptr = &from_gids_p3[0];
    to_ptr=&to_gids_p3[0];
  }

  Teuchos::ArrayView<GO> myfromgids (from_ptr, num_per_proc);
  Teuchos::ArrayView<GO> mytogids (to_ptr, num_per_proc);

  // FromMap (from.getRowMap() from Zoltan2)
  RCP<map_type> FromMap = rcp (new map_type (INVALID, myfromgids, indexBase, comm));

  // ToMap (tmap from Zoltan2)
  RCP<map_type> ToMap = rcp (new map_type (INVALID, mytogids, indexBase, comm));

  import_type Importer (FromMap, ToMap);

  // Duplicating what Zoltan2/Tpetra Does
  IntVector FromVector(FromMap);
  IntVector ToVector(ToMap);
  ToVector.putScalar(myRank);
  FromVector.putScalar(-666);

  FromVector.doExport(ToVector,Importer,Tpetra::REPLACE);

  Teuchos::ArrayRCP<const int> f_rcp = FromVector.getData();
  Teuchos::ArrayView<const int> f_view = f_rcp();
  Teuchos::ArrayRCP<const int> t_rcp = ToVector.getData();
  Teuchos::ArrayView<const int> t_view = t_rcp();

  ToVector.setObjectLabel ("ToVector");
  ToVector.describe (out, Teuchos::VERB_EXTREME);
  FromVector.setObjectLabel ("FromVector");
  FromVector.describe (out, Teuchos::VERB_EXTREME);

  // Check the "FromAnswer" answer against who_owns

  int lclSuccess = 1;
  int gblSuccess = 1;

  std::ostringstream lclErrStrm;
  for (size_t i = 0; i < FromMap->getNodeNumElements (); ++i) {
    if (f_view[i] != who_owns[FromMap->getGlobalElement(i)]) {
      lclErrStrm << "[" << myRank << "] ERROR: Ownership of GID "
                 << FromMap->getGlobalElement(i) << " is incorrect!" << endl;
      lclSuccess = 0;
    }
  }

  reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));

  if (gblSuccess == 1) {
    out << "End Result: TEST PASSED" << endl;
    return EXIT_SUCCESS;
  }
  else {
    out << "End Result: TEST FAILED" << endl;
    Teuchos::OSTab tab2 (out);
    Tpetra::Details::gathervPrint (out, lclErrStrm.str (), *comm);
    return EXIT_FAILURE;
  }
}
