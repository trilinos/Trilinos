// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

// Test for github issue #3101 -- test Tpetra::FEMultiVector::doOwnedToOwnedPlusShared

#include "Tpetra_Core.hpp"
#include "Tpetra_FEMultiVector.hpp"
#include "Teuchos_CommHelpers.hpp"

#include <string>
#include <sstream>
#include <iostream>

/////////////////////////////////////////////////////////////////////////////
// Each processor owns ten vertices and copies five from another processor.
// Each processor contributes values to all fifteen vertices it stores.
// Tpetra::FEMultiVector pushes the copies to their owning processors, and
// ADDs them to the owned values.
// Tpetra::FEMultiVector updates copies from their owning processors.

class FEMultiVectorTest {
public:

  typedef Tpetra::Map<> map_t;
  typedef map_t::local_ordinal_type lno_t;
  typedef map_t::global_ordinal_type gno_t;
  typedef int scalar_t;
  typedef Tpetra::FEMultiVector<scalar_t, lno_t, gno_t> femv_t;

  FEMultiVectorTest(Teuchos::RCP<const Teuchos::Comm<int> > &comm_);

  int intTest();

  void printFEMV(const char *msg) ;

private:
  int me;           // my processor rank
  int np;           // number of processors
  int nLocalOwned;  // number of vertices owned by this processor
  int nLocalCopy;   // number of copies of off-processor vertices on this proc
  int nVec;         // number of vectors in multivector

  Teuchos::RCP<const Teuchos::Comm<int> > comm;  // MPI communicator

  Teuchos::RCP<const map_t> mapWithCopies;  // Tpetra::Map including owned
                                            // vertices and copies
  Teuchos::RCP<const map_t> mapOwned;       // Tpetra::Map including only owned
  Teuchos::RCP<femv_t> femv;                // Tpetra::FEMultiVector for test
};

//////////////////////////////////////////////////////////////////////////////
// Constructor
// -  assigns vertices to processors
// -  builds map with owned vertices and
//           map with owned+ghosted (copied) vertices
// -  then builds a Tpetra::FEMultiVector with two vectors using the maps
FEMultiVectorTest::FEMultiVectorTest(
  Teuchos::RCP<const Teuchos::Comm<int> > &comm_
) :
  me(comm_->getRank()),
  np(comm_->getSize()),
  nLocalOwned(10),
  nLocalCopy( np > 1 ? 5 : 0),
  nVec(2),
  comm(comm_)
{
  // Each rank has 15 IDs, the last five of which overlap with the next rank.
  // (IDs and owning processors wrap-around from processor np-1 to 0.)
  const Tpetra::global_size_t nGlobal = np * nLocalOwned;
  lno_t offset = me * nLocalOwned;

  Teuchos::Array<gno_t> gids(nLocalOwned + nLocalCopy);
  for (lno_t i = 0 ; i < nLocalOwned + nLocalCopy; i++)
    gids[i] = static_cast<gno_t> ((offset + i) % nGlobal);

  // Create Map of owned + copies (a.k.a. overlap map); analagous to ColumnMap
  Tpetra::global_size_t dummy =
          Teuchos::OrdinalTraits<Tpetra::global_size_t>::invalid();
  mapWithCopies = rcp(new map_t(dummy, gids(), 0, comm));

  // Create Map of owned only (a.k.a. one-to-one map); analagous to RowMap
  mapOwned = rcp(new map_t(dummy, gids(0, nLocalOwned), 0, comm));

  // Print the entries of each map
  std::cout << me << " MAP WITH COPIES ("
                  << mapWithCopies->getGlobalNumElements() << "):  ";
  lno_t nlocal = lno_t(mapWithCopies->getLocalNumElements());
  for (lno_t i = 0; i < nlocal; i++)
    std::cout << mapWithCopies->getGlobalElement(i) << " ";
  std::cout << std::endl;

  std::cout << me << " ONE TO ONE MAP  ("
                  << mapOwned->getGlobalNumElements() << "):  ";
  nlocal = lno_t(mapOwned->getLocalNumElements());
  for (lno_t i = 0; i < nlocal; i++)
    std::cout << mapOwned->getGlobalElement(i) << " ";
  std::cout << std::endl;

  // Create FEMultiVector
  typedef Tpetra::Import<lno_t, gno_t> import_t;
  Teuchos::RCP<import_t> importer = rcp(new import_t(mapOwned,
                                                     mapWithCopies));
  femv = rcp(new femv_t(mapOwned, importer, nVec, true));
}

//////////////////////////////////////////////////////////////////////////////
// Test FEMultiVector
// Fill first vector with GIDs of vertices
// Fill second vector with processor rank (me)
// Perform communication to add copies' contributions to their owned vertices
// Perform communication to send owned vertices' values back to copies.
int FEMultiVectorTest::intTest()
{
  int ierr = 0;


  // Add contributions to owned vertices and copies of off-processor vertices
  try {
    femv->beginAssembly();
    for (lno_t i = 0; i < nLocalOwned + nLocalCopy; i++) {
      gno_t gid = mapWithCopies->getGlobalElement(i);
      femv->replaceGlobalValue(gid, 0, gid);
      femv->replaceGlobalValue(gid, 1, me);
    }
    femv->endAssembly();
  }
  catch (std::exception &e) {
    std::cout << "FAIL:  Exception thrown in Fill:  " << e.what() << std::endl;
    throw e;
  }

  printFEMV("After Fill ");

  // Update copied vertices from their owners
  try {
    femv->doOwnedToOwnedPlusShared(Tpetra::REPLACE);
  }
  catch (std::exception &e) {
    std::cout << "FAIL:  Exception thrown in doOwnedToOwnedPlusShared:  "
              << e.what() << std::endl;
    throw e;
  }

  printFEMV("After doOwnedToOwnedPlusShared ");

  // Check results:  after ADD in endAssembly,
  // -  overlapping entries of vec 0 should be 2 * gid
  //    nonoverlapping entries of vec 0 should be gid
  // -  overlapping entries of vec 1 should be me + (np + me-1) % np;
  //    nonoverlapping entries of vec 1 should be me

  // Check vector 0:  SOURCE
  auto value = femv->getData(0);
  for (lno_t i = 0; i < nLocalCopy; i++){
    gno_t gid = femv->getMap()->getGlobalElement(i);
    if (value[i] != 2*gid) {
      std::cout << me << " Error in SOURCE vec 0 overlap: gid=" << gid
                      << " value= " << value[i] << " should be " << 2*gid
                      << std::endl;
      ierr++;
    }
  }
  for (lno_t i = nLocalCopy; i < nLocalOwned; i++) {
    gno_t gid = femv->getMap()->getGlobalElement(i);
    if (value[i] != gid) {
      std::cout << me << " Error in SOURCE vec 0:  gid=" << gid
                      << " value= " << value[i] << " should be " << gid
                      << std::endl;
      ierr++;
    }
  }

  // Check vector 0:  TARGET
  femv->switchActiveMultiVector();  // needed to allow testing copies
  value = femv->getData(0);
  for (lno_t i = 0; i < nLocalCopy; i++){
    gno_t gid = femv->getMap()->getGlobalElement(i);
    if (value[i] != 2*gid) {
      std::cout << me << " Error in TARGET vec 0 overlap: gid=" << gid
                      << " value= " << value[i] << " should be " << 2*gid
                      << std::endl;
      ierr++;
    }
  }
  for (lno_t i = nLocalCopy; i < nLocalOwned; i++) {
    gno_t gid = femv->getMap()->getGlobalElement(i);
    if (value[i] != gid) {
      std::cout << me << " Error in TARGET vec 0:  gid=" << gid
                      << " value= " << value[i] << " should be " << gid
                      << std::endl;
      ierr++;
    }
  }
  for (lno_t i = nLocalOwned; i < nLocalCopy + nLocalOwned; i++) {
    gno_t gid = femv->getMap()->getGlobalElement(i);
    if (value[i] != 2*gid) {
      std::cout << me << " Error in TARGET vec 0 copies:  gid=" << gid
                      << " value= " << value[i] << " should be " << gid
                      << std::endl;
      ierr++;
    }
  }
  femv->switchActiveMultiVector();  // restore state

  // Check vector 1:  SOURCE
  value = femv->getData(1);
  int tmp = me + (np + me - 1) % np;
  for (lno_t i = 0; i < nLocalCopy; i++) {
    gno_t gid = mapWithCopies->getGlobalElement(i);
    if (value[i] != tmp) {
      std::cout << me << " Error in SOURCE vec 1 overlap:  gid=" << gid
                      << " value= " << value[i] << " should be " << tmp
                      << std::endl;
      ierr++;
    }
  }
  for (lno_t i = nLocalCopy; i < nLocalOwned; i++) {
    gno_t gid = mapWithCopies->getGlobalElement(i);
    if (value[i] != me) {
      std::cout << me << " Error in SOURCE vec 1:  gid=" << gid
                      << " value= " << value[i] << " should be " << me
                      << std::endl;
      ierr++;
    }
  }

  // Check vector 1:  TARGET
  femv->switchActiveMultiVector();  // needed to allow testing copies
  value = femv->getData(1);
  for (lno_t i = 0; i < nLocalCopy; i++) {
    gno_t gid = mapWithCopies->getGlobalElement(i);
    if (value[i] != tmp) {
      std::cout << me << " Error in TARGET vec 1 overlap:  gid=" << gid
                      << " value= " << value[i] << " should be " << tmp
                      << std::endl;
      ierr++;
    }
  }
  for (lno_t i = nLocalCopy; i < nLocalOwned; i++) {
    gno_t gid = mapWithCopies->getGlobalElement(i);
    if (value[i] != me) {
      std::cout << me << " Error in TARGET vec 1:  gid=" << gid
                      << " value= " << value[i] << " should be " << me
                      << std::endl;
      ierr++;
    }
  }
  tmp = me + (me + 1) % np;
  for (lno_t i = nLocalOwned; i < nLocalCopy + nLocalOwned; i++) {
    gno_t gid = mapWithCopies->getGlobalElement(i);
    if (value[i] != tmp)  {
      std::cout << me << " Error in TARGET vec 1 copies:  gid=" << gid
                      << " value= " << value[i] << " should be " << tmp
                      << std::endl;
      ierr++;
    }
  }
  femv->switchActiveMultiVector();  // restore state

  return ierr;
}

//////////////////////////////////////////////////////////////////////////////
// Print the SOURCE and TARGET multivector entries
void FEMultiVectorTest::printFEMV(const char *msg)
{
  // print Source MV (Owned only)
  for (int v = 0; v < nVec; v++) {
    std::cout << me << " SOURCE " << msg << " FEMV[" << v << "] Owned: ";
    auto value = femv->getData(v);
    for (lno_t i = 0; i < nLocalOwned; i++) std::cout << value[i] << " ";
    std::cout << std::endl;
  }

  // print Target MV (Owned + copies)
  femv->switchActiveMultiVector();  // needed to allow printing copies
  for (int v = 0; v < nVec; v++) {
    std::cout << me << " TARGET " << msg << " FEMV[" << v << "] Owned: ";
    auto value = femv->getData(v);
    for (lno_t i = 0; i < nLocalOwned; i++) std::cout << value[i] << " ";
    std::cout << " Copies: ";
    for (lno_t i = nLocalOwned; i < nLocalOwned+nLocalCopy; i++)
      std::cout << value[i] << " ";
    std::cout << std::endl;
  }
  femv->switchActiveMultiVector();  // restore state
}

///////////////////////////////////////////////////////////////////////////////

int main(int narg, char **arg)
{
  Tpetra::ScopeGuard scope(&narg, &arg);
  Teuchos::RCP<const Teuchos::Comm<int> > comm = Tpetra::getDefaultComm();
  int me = comm->getRank();
  int ierr = 0;

  FEMultiVectorTest femvTest(comm);

  try {
    ierr = femvTest.intTest();
  }
  catch (std::exception &e) {
    std::cout << "FAIL: Exception thrown in intTest " << e.what() << std::endl;
    return 0;
  }

  int gerr = 0;
  Teuchos::reduceAll<int, int>(*comm, Teuchos::REDUCE_SUM, 1, &ierr, &gerr);

  if (me == 0) {
    if (gerr == 0) std::cout << "PASS" << std::endl;
    else std::cout << "FAIL:  " << gerr << " failures" << std::endl;
  }

  return 0;
}
