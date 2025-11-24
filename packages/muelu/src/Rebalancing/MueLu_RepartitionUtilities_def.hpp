// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef MUELU_REPARTITIONUTILITIES_DEF_HPP
#define MUELU_REPARTITIONUTILITIES_DEF_HPP

#include <algorithm>
#include <iostream>
#include <sstream>

#include "MueLu_RepartitionUtilities_decl.hpp"  // TMP JG NOTE: before other includes, otherwise I cannot test the fwd declaration in _def

#ifdef HAVE_MPI
#include <Teuchos_DefaultMpiComm.hpp>
#include <Teuchos_CommHelpers.hpp>
#include <Teuchos_Details_MpiTypeTraits.hpp>

#include <Xpetra_Map.hpp>
#include <Xpetra_MapFactory.hpp>
#include <Xpetra_MultiVectorFactory.hpp>
#include <Xpetra_VectorFactory.hpp>
#include <Xpetra_Import.hpp>
#include <Xpetra_ImportFactory.hpp>
#include <Xpetra_Export.hpp>
#include <Xpetra_ExportFactory.hpp>

#include "MueLu_Utilities.hpp"

#include "MueLu_CloneRepartitionInterface.hpp"

#include "MueLu_Level.hpp"
#include "MueLu_MasterList.hpp"
#include "MueLu_Monitor.hpp"
#include "MueLu_PerfUtils.hpp"

namespace MueLu {

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
std::tuple<Teuchos::Array<GlobalOrdinal>, GlobalOrdinal>
RepartitionUtilities<Scalar, LocalOrdinal, GlobalOrdinal, Node>::ConstructGIDs(RCP<GOVector> decomposition) {
  ArrayRCP<const GO> decompEntries;
  if (decomposition->getLocalLength() > 0)
    decompEntries = decomposition->getData(0);

  auto rowMap               = decomposition->getMap();
  auto comm                 = rowMap->getComm();
  const auto myRank         = comm->getRank();
  const auto numProcs       = comm->getSize();
  Xpetra::UnderlyingLib lib = rowMap->lib();

  Array<GO> myGIDs;
  myGIDs.reserve(decomposition->getLocalLength());

  // Step 0: Construct mapping
  //    part number -> GIDs I own which belong to this part
  // NOTE: my own part GIDs are not part of the map
  typedef std::map<GO, Array<GO> > map_type;
  map_type sendMap;
  for (LO i = 0; i < decompEntries.size(); i++) {
    GO id  = decompEntries[i];
    GO GID = rowMap->getGlobalElement(i);

    if (id == myRank)
      myGIDs.push_back(GID);
    else
      sendMap[id].push_back(GID);
  }
  decompEntries = Teuchos::null;

  int numSend = sendMap.size(), numRecv;

  // Arrayify map keys
  Array<GO> myParts(numSend), myPart(1);
  int cnt   = 0;
  myPart[0] = myRank;
  for (typename map_type::const_iterator it = sendMap.begin(); it != sendMap.end(); it++)
    myParts[cnt++] = it->first;

  const GO numLocalKept = myGIDs.size();

  // Step 1: Find out how many processors send me data
  // partsIndexBase starts from zero, as the processors ids start from zero
  {
    GO partsIndexBase       = 0;
    RCP<Map> partsIHave     = MapFactory ::Build(lib, Teuchos::OrdinalTraits<Xpetra::global_size_t>::invalid(), myParts(), partsIndexBase, comm);
    RCP<Map> partsIOwn      = MapFactory ::Build(lib, numProcs, myPart(), partsIndexBase, comm);
    RCP<Export> partsExport = ExportFactory::Build(partsIHave, partsIOwn);

    RCP<GOVector> partsISend    = Xpetra::VectorFactory<GO, LO, GO, NO>::Build(partsIHave);
    RCP<GOVector> numPartsIRecv = Xpetra::VectorFactory<GO, LO, GO, NO>::Build(partsIOwn);
    if (numSend) {
      ArrayRCP<GO> partsISendData = partsISend->getDataNonConst(0);
      for (int i = 0; i < numSend; i++)
        partsISendData[i] = 1;
    }
    (numPartsIRecv->getDataNonConst(0))[0] = 0;

    numPartsIRecv->doExport(*partsISend, *partsExport, Xpetra::ADD);
    numRecv = (numPartsIRecv->getData(0))[0];
  }

  RCP<const Teuchos::MpiComm<int> > tmpic = rcp_dynamic_cast<const Teuchos::MpiComm<int> >(comm);
  TEUCHOS_TEST_FOR_EXCEPTION(tmpic == Teuchos::null, Exceptions::RuntimeError, "Cannot cast base Teuchos::Comm to Teuchos::MpiComm object.");
  RCP<const Teuchos::OpaqueWrapper<MPI_Comm> > rawMpiComm = tmpic->getRawMpiComm();

  // Step 2: Get my GIDs from everybody else
  MPI_Datatype MpiType = Teuchos::Details::MpiTypeTraits<GO>::getType();
  int msgTag           = 12345;  // TODO: use Comm::dup for all internal messaging

  // Post sends
  Array<MPI_Request> sendReqs(numSend);
  cnt = 0;
  for (typename map_type::iterator it = sendMap.begin(); it != sendMap.end(); it++)
    MPI_Isend(static_cast<void*>(it->second.getRawPtr()), it->second.size(), MpiType, Teuchos::as<GO>(it->first), msgTag, *rawMpiComm, &sendReqs[cnt++]);

  map_type recvMap;
  size_t totalGIDs = myGIDs.size();
  for (int i = 0; i < numRecv; i++) {
    MPI_Status status;
    MPI_Probe(MPI_ANY_SOURCE, msgTag, *rawMpiComm, &status);

    // Get rank and number of elements from status
    int fromRank = status.MPI_SOURCE, count;
    MPI_Get_count(&status, MpiType, &count);

    recvMap[fromRank].resize(count);
    MPI_Recv(static_cast<void*>(recvMap[fromRank].getRawPtr()), count, MpiType, fromRank, msgTag, *rawMpiComm, &status);

    totalGIDs += count;
  }

  // Do waits on send requests
  if (numSend) {
    Array<MPI_Status> sendStatuses(numSend);
    MPI_Waitall(numSend, sendReqs.getRawPtr(), sendStatuses.getRawPtr());
  }

  // Merge GIDs
  myGIDs.reserve(totalGIDs);
  for (typename map_type::const_iterator it = recvMap.begin(); it != recvMap.end(); it++) {
    int offset = myGIDs.size(), len = it->second.size();
    if (len) {
      myGIDs.resize(offset + len);
      memcpy(myGIDs.getRawPtr() + offset, it->second.getRawPtr(), len * sizeof(GO));
    }
  }
  // NOTE 2: The general sorting algorithm could be sped up by using the knowledge that original myGIDs and all received chunks
  // (i.e. it->second) are sorted. Therefore, a merge sort would work well in this situation.
  std::sort(myGIDs.begin(), myGIDs.end());

  return std::make_tuple(myGIDs, numLocalKept);
}

}  // namespace MueLu

#endif  // ifdef HAVE_MPI

#endif  // MUELU_REPARTITIONUTILITIES_DEF_HPP
