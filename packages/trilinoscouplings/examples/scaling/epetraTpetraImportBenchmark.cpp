// @HEADER
// ************************************************************************
//
//                           Intrepid Package
//                 Copyright (2007) Sandia Corporation
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
// Questions? Contact Pavel Bochev  (pbboche@sandia.gov),
//                    Denis Ridzal  (dridzal@sandia.gov),
//                    Kara Peterson (kjpeter@sandia.gov).
//
// ************************************************************************
// @HEADER

// TrilinosCouplings includes
#include <TrilinosCouplings_config.h>

#include <Epetra_Import.h>
#include <Epetra_Map.h>
#include <Epetra_Vector.h>
#ifdef EPETRA_MPI
#  include <Epetra_MpiComm.h>
#  include <mpi.h>
#else
#  include <Epetra_SerialComm.h>
#endif

#include <Tpetra_Import.hpp>
#include <Tpetra_Vector.hpp>
#include <Teuchos_TimeMonitor.hpp>

using Tpetra::global_size_t;
using Teuchos::ArrayView;
using Teuchos::as;
using Teuchos::Comm;
using Teuchos::RCP;
using Teuchos::rcp;
using Teuchos::Time;
using Teuchos::TimeMonitor;

typedef double ST;
typedef int LO;
typedef int GO; // So that Epetra and Tpetra can use the same GID lists
typedef Kokkos::SerialNode NT;

// Create a new timer with the given name if it hasn't already been
// created, else get the previously created timer with that name.
RCP<Time> getTimer (const std::string& timerName) {
  RCP<Time> timer = 
    TimeMonitor::lookupCounter (timerName);
  if (timer.is_null ()) {
    timer = TimeMonitor::getNewCounter (timerName);
  }
  return timer;
}

void 
benchmarkTpetraImport (ArrayView<const GO> srcGlobalElts,
		       ArrayView<const GO> destGlobalElts,
		       const GO indexBase,
		       RCP<const Comm<int> > comm,
		       RCP<NT> node,
		       const int numMapCreateTrials,
		       const int numImportCreateTrials,
		       const int numVectorCreateTrials,
		       const int numExecTrials)
{
  typedef Tpetra::Import<LO, GO, NT> import_type;
  typedef Tpetra::Map<LO, GO, NT> map_type;
  typedef Tpetra::Vector<ST, LO, GO, NT> vector_type;
  const global_size_t INVALID = Teuchos::OrdinalTraits<global_size_t>::invalid ();

  TEUCHOS_TEST_FOR_EXCEPTION(
    numMapCreateTrials < 1 && numImportCreateTrials > 0, std::invalid_argument,
    "numMapCreateTrials must be > 0 if numImportCreateTrials > 0.");
  TEUCHOS_TEST_FOR_EXCEPTION(
    numImportCreateTrials < 1 && numExecTrials > 0, std::invalid_argument, 
    "numImportCreateTrials must be > 0 if numExecTrials > 0.");
  TEUCHOS_TEST_FOR_EXCEPTION(
    numVectorCreateTrials < 1 && numMapCreateTrials > 0, std::invalid_argument, 
    "numVectorCreateTrials must be > 0 if numExecTrials > 0.");
  TEUCHOS_TEST_FOR_EXCEPTION(
    numVectorCreateTrials < 1 && numExecTrials > 0, std::invalid_argument, 
    "numVectorCreateTrials must be > 0 if numExecTrials > 0.");

  RCP<Time> mapCreateTimer = getTimer ("Tpetra: Map: Create");
  RCP<Time> importCreateTimer = getTimer ("Tpetra: Import: Create");
  RCP<Time> vectorCreateTimer = getTimer ("Tpetra: Vector: Create");
  RCP<Time> importExecTimer = getTimer ("Tpetra: Import: Execute");

  RCP<map_type> srcMap, destMap;
  {
    TimeMonitor timeMon (*mapCreateTimer);
    for (int k = 0; k < numMapCreateTrials; ++k) {
      srcMap = rcp (new map_type (invalid, srcGlobalElts, indexBase, comm, node));
      destMap = rcp (new map_type (invalid, destGlobalElts, indexBase, comm, node));
    }
  }
  RCP<import_type> import;
  {
    TimeMonitor timeMon (*importCreateTimer);
    for (int k = 0; k < numCreateTrials; ++k) {
      import = rcp (new import_type (srcMap, destMap));
    }
  }
  RCP<vector_type> srcVec, destVec;
  {
    TimeMonitor timeMon (*vectorCreateTimer);
    for (int k = 0; k < numExecTrials; ++k) {
      srcVec = rcp (new vector_type (srcMap));
      destVec = rcp (new vector_type (destMap));
    }
  }
  {
    TimeMonitor timeMon (*importExecTimer);
    for (int k = 0; k < numExecTrials; ++k) {
      destVec->doImport (*import, *srcVec, Tpetra::ADD);
    }
  }
}

void 
benchmarkEpetraImport (ArrayView<const int> srcGlobalElts,
		       ArrayView<const int> destGlobalElts,
		       const int indexBase,
		       const Epetra_Comm& comm,
		       const int numMapCreateTrials,
		       const int numImportCreateTrials,
		       const int numVectorCreateTrials,
		       const int numExecTrials)
{
  const int INVALID = -1;

  TEUCHOS_TEST_FOR_EXCEPTION(
    numMapCreateTrials < 1 && numImportCreateTrials > 0, std::invalid_argument,
    "numMapCreateTrials must be > 0 if numImportCreateTrials > 0.");
  TEUCHOS_TEST_FOR_EXCEPTION(
    numImportCreateTrials < 1 && numExecTrials > 0, std::invalid_argument, 
    "numImportCreateTrials must be > 0 if numExecTrials > 0.");
  TEUCHOS_TEST_FOR_EXCEPTION(
    numVectorCreateTrials < 1 && numMapCreateTrials > 0, std::invalid_argument, 
    "numVectorCreateTrials must be > 0 if numExecTrials > 0.");
  TEUCHOS_TEST_FOR_EXCEPTION(
    numVectorCreateTrials < 1 && numExecTrials > 0, std::invalid_argument, 
    "numVectorCreateTrials must be > 0 if numExecTrials > 0.");

  RCP<Time> mapCreateTimer = getTimer ("Epetra: Map: Create");
  RCP<Time> importCreateTimer = getTimer ("Epetra: Import: Create");
  RCP<Time> vectorCreateTimer = getTimer ("Epetra: Vector: Create");
  RCP<Time> importExecTimer = getTimer ("Epetra: Import: Execute");

  RCP<Epetra_Map> srcMap, destMap;
  {
    TimeMonitor timeMon (*mapCreateTimer);
    for (int k = 0; k < numMapCreateTrials; ++k) {
      const int srcNumElts = as<int> (srcGlobalElts.size ());
      const int* const srcElts = srcGlobalElts.getRawPtr ();
      srcMap = rcp (new Epetra_Map (INVALID, srcNumElts, srcElts, indexBase, comm));

      const int destNumElts = as<int> (destGlobalElts.size ());
      const int* const destElts = destGlobalElts.getRawPtr ();
      destMap = rcp (new Epetra_Map (INVALID, destNumElts, destElts, indexBase, comm));
    }
  }
  RCP<Epetra_Import> import;
  {
    TimeMonitor timeMon (*mapCreateTimer);
    for (int k = 0; k < numImportCreateTrials; ++k) {
      import = rcp (new Epetra_Import (*srcMap, *destMap));
    }
  }
  RCP<Epetra_Vector> srcVec, destVec;
  {
    TimeMonitor timeMon (*vectorCreateTimer);
    for (int k = 0; k < numVectorCreateTrials; ++k) {
      srcVec = rcp (new Epetra_Vector (*srcMap));
      destVec = rcp (new Epetra_Vector (*destMap));
    }
  }
  {
    TimeMonitor timeMon (*importExecTimer);
    for (int k = 0; k < numExecTrials; ++k) {
      (void) destVec->Import (*srcVec, *import, Add);
    }
  }
}

void 
createGidLists (Array<GO>& srcGlobalElts,
		Array<GO>& destGlobalElts,
		const int numProcs,
		const int myRank,
		const int numEltsPerProc)
{
  const int overlap = (myRank == 0 || myRank == numProcs-1) ? 1 : 2;

  srcGlobalElts.resize (numEltsPerProc);
  destGlobalElts.resize (numEltsPerProc + overlap);

  const GO myStartSrcGid = indexBase + myRank * numEltsPerProc;
  const GO myEndSrcGid = myStartSrcGid + numEltsPerProc - 1; // inclusive

  // Put the GIDs in reverse order, to simulate a noncontiguously
  // ordered Map.  Neither Epetra nor Tpetra optimize for this case.
  for (int k = 0; k < numEltsPerProc; ++k) {
    srcGlobalElts[numEltsPerProc - k - 1] = myStartSrcGid + as<GO> (k);
  }

  // 1-D Poisson overlap.
  if (myRank == 0) {
    destGlobalElts[numEltsPerProc] = myEndSrcGid + 1;
  } else if (myRank == numProcs-1) {
    destGlobalElts[numEltsPerProc] = myStartSrcGid - 1; 
  } else {
    destGlobalElts[numEltsPerProc] = myStartSrcGid - 1;
    destGlobalElts[numEltsPerProc+1] = myEndSrcGid + 1;
  }
}


int main (int argc, char* argv[]) {
  Teuchos::oblackholestream blackHole;
  Teuchos::GlobalMPISession mpiSession (&argc, &argv, &blackHole);
  RCP<const Teuchos::Comm<int> > tpetraComm = 
    Tpetra::DefaultPlatform::getDefaultPlatform ().getComm ();
  RCP<NT> node = Kokkos::Details::getNode<NT> ();

#ifdef EPETRA_MPI
  Epetra_MpiComm epetraComm (MPI_COMM_WORLD);
#else
  Epetra_SerialComm epetraComm;
#endif EPETRA_MPI

  const int numProcs = tpetraComm->getSize ();
  const int myRank = tpetraComm->getRank ();
  const int indexBase = 1; // of interest to Sierra

  // Benchmark parameters
  const int numEltsPerProc = 10000;
  const int numMapCreateTrials = 100;
  const int numImportCreateTrials = 100;
  const int numVectorCreateTrials = 100;
  const int numExecTrials = 100;

  // Run the benchmark
  Array<GO> srcGlobalElts, destGlobalElts;
  createGidLists (srcGlobalElts, destGlobalElts, numProcs, myRank, 
		  numEltsPerProc, indexBase);
  benchmarkEpetraImport (srcGlobalElts, destGlobalElts, indexBase, epetraComm,
			 numMapCreateTrials, numImportCreateTrials,
			 numVectorCreateTrials, numExecTrials);
  benchmarkTpetraImport (srcGlobalElts, destGlobalElts, indexBase, comm, node, 
			 numMapCreateTrials, numImportCreateTrials,
			 numVectorCreateTrials, numExecTrials);
  TimeMonitor::report (tpetraComm.ptr ());
  return EXIT_SUCCESS;
}
