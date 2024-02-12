// @HEADER
//
// ***********************************************************************
//
//        MueLu: A package for multigrid based preconditioning
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
#ifndef MUELU_ZOLTANINTERFACE_DEF_HPP
#define MUELU_ZOLTANINTERFACE_DEF_HPP

#include "MueLu_ZoltanInterface_decl.hpp"
#if defined(HAVE_MUELU_ZOLTAN) && defined(HAVE_MPI)

#include <Teuchos_Utils.hpp>
#include <Teuchos_DefaultMpiComm.hpp>  //TODO: fwd decl.
#include <Teuchos_OpaqueWrapper.hpp>   //TODO: fwd decl.

#include "MueLu_Level.hpp"
#include "MueLu_Exceptions.hpp"
#include "MueLu_Monitor.hpp"

namespace MueLu {

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<const ParameterList> ZoltanInterface<Scalar, LocalOrdinal, GlobalOrdinal, Node>::GetValidParameterList() const {
  RCP<ParameterList> validParamList = rcp(new ParameterList());

  validParamList->set<RCP<const FactoryBase> >("A", Teuchos::null, "Factory of the matrix A");
  validParamList->set<RCP<const FactoryBase> >("number of partitions", Teuchos::null, "Instance of RepartitionHeuristicFactory.");
  validParamList->set<RCP<const FactoryBase> >("Coordinates", Teuchos::null, "Factory of the coordinates");

  return validParamList;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void ZoltanInterface<Scalar, LocalOrdinal, GlobalOrdinal, Node>::DeclareInput(Level &currentLevel) const {
  Input(currentLevel, "A");
  Input(currentLevel, "number of partitions");
  Input(currentLevel, "Coordinates");
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void ZoltanInterface<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build(Level &level) const {
  FactoryMonitor m(*this, "Build", level);

  RCP<Matrix> A            = Get<RCP<Matrix> >(level, "A");
  RCP<BlockedCrsMatrix> bA = Teuchos::rcp_dynamic_cast<BlockedCrsMatrix>(A);
  RCP<const Map> rowMap;
  if (bA != Teuchos::null) {
    // Extracting the full the row map here...
    RCP<const Map> bArowMap       = bA->getRowMap();
    RCP<const BlockedMap> bRowMap = Teuchos::rcp_dynamic_cast<const BlockedMap>(bArowMap);
    rowMap                        = bRowMap->getFullMap();
  } else {
    rowMap = A->getRowMap();
  }

  typedef Xpetra::MultiVector<typename Teuchos::ScalarTraits<Scalar>::magnitudeType, LocalOrdinal, GlobalOrdinal, Node> double_multivector_type;
  RCP<double_multivector_type> Coords = Get<RCP<double_multivector_type> >(level, "Coordinates");
  size_t dim                          = Coords->getNumVectors();
  int numParts                        = Get<int>(level, "number of partitions");

  if (numParts == 1 || numParts == -1) {
    // Running on one processor, so decomposition is the trivial one, all zeros.
    RCP<Xpetra::Vector<GO, LO, GO, NO> > decomposition = Xpetra::VectorFactory<GO, LO, GO, NO>::Build(rowMap, true);
    Set(level, "Partition", decomposition);
    return;
  } else if (numParts == -1) {
    // No repartitioning
    RCP<Xpetra::Vector<GO, LO, GO, NO> > decomposition = Teuchos::null;
    Set(level, "Partition", decomposition);
    return;
  }

  float zoltanVersion_;
  Zoltan_Initialize(0, NULL, &zoltanVersion_);

  RCP<const Teuchos::MpiComm<int> > dupMpiComm            = rcp_dynamic_cast<const Teuchos::MpiComm<int> >(rowMap->getComm()->duplicate());
  RCP<const Teuchos::OpaqueWrapper<MPI_Comm> > zoltanComm = dupMpiComm->getRawMpiComm();

  RCP<Zoltan> zoltanObj_ = rcp(new Zoltan((*zoltanComm)()));  // extract the underlying MPI_Comm handle and create a Zoltan object
  if (zoltanObj_ == Teuchos::null)
    throw Exceptions::RuntimeError("MueLu::Zoltan : Unable to create Zoltan data structure");

  // Tell Zoltan what kind of local/global IDs we will use.
  // In our case, each GID is two ints and there are no local ids.
  // One can skip this step if the IDs are just single ints.
  int rv;
  if ((rv = zoltanObj_->Set_Param("num_gid_entries", "1")) != ZOLTAN_OK)
    throw Exceptions::RuntimeError("MueLu::Zoltan::Setup : setting parameter 'num_gid_entries' returned error code " + Teuchos::toString(rv));
  if ((rv = zoltanObj_->Set_Param("num_lid_entries", "0")) != ZOLTAN_OK)
    throw Exceptions::RuntimeError("MueLu::Zoltan::Setup : setting parameter 'num_lid_entries' returned error code " + Teuchos::toString(rv));
  if ((rv = zoltanObj_->Set_Param("obj_weight_dim", "1")) != ZOLTAN_OK)
    throw Exceptions::RuntimeError("MueLu::Zoltan::Setup : setting parameter 'obj_weight_dim' returned error code " + Teuchos::toString(rv));

  if (GetVerbLevel() & Statistics1)
    zoltanObj_->Set_Param("debug_level", "1");
  else
    zoltanObj_->Set_Param("debug_level", "0");

  zoltanObj_->Set_Param("num_global_partitions", toString(numParts));

  zoltanObj_->Set_Num_Obj_Fn(GetLocalNumberOfRows, (void *)A.getRawPtr());
  zoltanObj_->Set_Obj_List_Fn(GetLocalNumberOfNonzeros, (void *)A.getRawPtr());
  zoltanObj_->Set_Num_Geom_Fn(GetProblemDimension, (void *)&dim);
  zoltanObj_->Set_Geom_Multi_Fn(GetProblemGeometry, (void *)Coords.get());

  // Data pointers that Zoltan requires.
  ZOLTAN_ID_PTR import_gids = NULL;  // Global nums of objs to be imported
  ZOLTAN_ID_PTR import_lids = NULL;  // Local indices to objs to be imported
  int *import_procs         = NULL;  // Proc IDs of procs owning objs to be imported.
  int *import_to_part       = NULL;  // Partition #s to which imported objs should be assigned.
  ZOLTAN_ID_PTR export_gids = NULL;  // Global nums of objs to be exported
  ZOLTAN_ID_PTR export_lids = NULL;  // local indices to objs to be exported
  int *export_procs         = NULL;  // Proc IDs of destination procs for objs to be exported.
  int *export_to_part       = NULL;  // Partition #s for objs to be exported.
  int num_imported;                  // Number of objs to be imported.
  int num_exported;                  // Number of objs to be exported.
  int newDecomp;                     // Flag indicating whether the decomposition has changed
  int num_gid_entries;               // Number of array entries in a global ID.
  int num_lid_entries;

  {
    SubFactoryMonitor m1(*this, "Zoltan RCB", level);
    rv = zoltanObj_->LB_Partition(newDecomp, num_gid_entries, num_lid_entries,
                                  num_imported, import_gids, import_lids, import_procs, import_to_part,
                                  num_exported, export_gids, export_lids, export_procs, export_to_part);
    if (rv == ZOLTAN_FATAL)
      throw Exceptions::RuntimeError("Zoltan::LB_Partition() returned error code");
  }

  // TODO check that A's row map is 1-1.  Zoltan requires this.

  RCP<Xpetra::Vector<GO, LO, GO, NO> > decomposition;
  if (newDecomp) {
    decomposition              = Xpetra::VectorFactory<GO, LO, GO, NO>::Build(rowMap, false);  // Don't initialize, will be overwritten
    ArrayRCP<GO> decompEntries = decomposition->getDataNonConst(0);

    int mypid = rowMap->getComm()->getRank();
    for (typename ArrayRCP<GO>::iterator i = decompEntries.begin(); i != decompEntries.end(); ++i)
      *i = mypid;

    LO blockSize = A->GetFixedBlockSize();
    for (int i = 0; i < num_exported; ++i) {
      // We have assigned Zoltan gids to first row GID in the block
      // NOTE: Zoltan GIDs are different from GIDs in the Coordinates vector
      LO localEl  = rowMap->getLocalElement(export_gids[i]);
      int partNum = export_to_part[i];
      for (LO j = 0; j < blockSize; ++j)
        decompEntries[localEl + j] = partNum;
    }
  }

  Set(level, "Partition", decomposition);

  zoltanObj_->LB_Free_Part(&import_gids, &import_lids, &import_procs, &import_to_part);
  zoltanObj_->LB_Free_Part(&export_gids, &export_lids, &export_procs, &export_to_part);

}  // Build()

//-------------------------------------------------------------------------------------------------------------
// GetLocalNumberOfRows
//-------------------------------------------------------------------------------------------------------------

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
int ZoltanInterface<Scalar, LocalOrdinal, GlobalOrdinal, Node>::GetLocalNumberOfRows(void *data, int *ierr) {
  if (data == NULL) {
    *ierr = ZOLTAN_FATAL;
    return -1;
  }
  Matrix *A = (Matrix *)data;
  *ierr     = ZOLTAN_OK;

  LO blockSize = A->GetFixedBlockSize();
  TEUCHOS_TEST_FOR_EXCEPTION(blockSize == 0, Exceptions::RuntimeError, "MueLu::Zoltan : Matrix has block size 0.");

  return A->getRowMap()->getLocalNumElements() / blockSize;
}  // GetLocalNumberOfRows()

//-------------------------------------------------------------------------------------------------------------
// GetLocalNumberOfNonzeros
//-------------------------------------------------------------------------------------------------------------

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void ZoltanInterface<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    GetLocalNumberOfNonzeros(void *data, int NumGidEntries, int /* NumLidEntries */, ZOLTAN_ID_PTR gids,
                             ZOLTAN_ID_PTR /* lids */, int /* wgtDim */, float *weights, int *ierr) {
  if (data == NULL || NumGidEntries < 1) {
    *ierr = ZOLTAN_FATAL;
    return;
  } else {
    *ierr = ZOLTAN_OK;
  }

  Matrix *A          = (Matrix *)data;
  RCP<const Map> map = A->getRowMap();

  LO blockSize = A->GetFixedBlockSize();
  TEUCHOS_TEST_FOR_EXCEPTION(blockSize == 0, Exceptions::RuntimeError, "MueLu::Zoltan : Matrix has block size 0.");

  size_t numElements          = map->getLocalNumElements();
  ArrayView<const GO> mapGIDs = map->getLocalElementList();

  if (blockSize == 1) {
    for (size_t i = 0; i < numElements; i++) {
      gids[i]    = as<ZOLTAN_ID_TYPE>(mapGIDs[i]);
      weights[i] = A->getNumEntriesInLocalRow(i);
    }

  } else {
    LO numBlockElements = numElements / blockSize;

    for (LO i = 0; i < numBlockElements; i++) {
      // Assign zoltan GID to the first row GID in the block
      // NOTE: Zoltan GIDs are different from GIDs in the Coordinates vector
      gids[i]    = as<ZOLTAN_ID_TYPE>(mapGIDs[i * blockSize]);
      weights[i] = 0.0;
      for (LO j = 0; j < blockSize; j++)
        weights[i] += A->getNumEntriesInLocalRow(i * blockSize + j);
    }
  }
}

//-------------------------------------------------------------------------------------------------------------
// GetProblemDimension
//-------------------------------------------------------------------------------------------------------------

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
int ZoltanInterface<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    GetProblemDimension(void *data, int *ierr) {
  int dim = *((int *)data);
  *ierr   = ZOLTAN_OK;

  return dim;
}  // GetProblemDimension

//-------------------------------------------------------------------------------------------------------------
// GetProblemGeometry
//-------------------------------------------------------------------------------------------------------------

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void ZoltanInterface<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    GetProblemGeometry(void *data, int /* numGIDEntries */, int /* numLIDEntries */, int numObjectIDs,
                       ZOLTAN_ID_PTR /* gids */, ZOLTAN_ID_PTR /* lids */, int dim, double *coordinates, int *ierr) {
  if (data == NULL) {
    *ierr = ZOLTAN_FATAL;
    return;
  }

  typedef Xpetra::MultiVector<typename Teuchos::ScalarTraits<Scalar>::coordinateType, LocalOrdinal, GlobalOrdinal, Node> double_multivector_type;
  double_multivector_type *Coords = (double_multivector_type *)data;

  if (dim != Teuchos::as<int>(Coords->getNumVectors())) {
    // FIXME I'm assuming dim should be 1, 2, or 3 coming in?!
    *ierr = ZOLTAN_FATAL;
    return;
  }

  TEUCHOS_TEST_FOR_EXCEPTION(numObjectIDs != Teuchos::as<int>(Coords->getLocalLength()), Exceptions::Incompatible, "Length of coordinates must be the same as the number of objects");

  ArrayRCP<ArrayRCP<const typename Teuchos::ScalarTraits<Scalar>::coordinateType> > CoordsData(dim);
  for (int j = 0; j < dim; ++j)
    CoordsData[j] = Coords->getData(j);

  size_t numElements = Coords->getLocalLength();
  for (size_t i = 0; i < numElements; ++i)
    for (int j = 0; j < dim; ++j)
      coordinates[i * dim + j] = (double)CoordsData[j][i];

  *ierr = ZOLTAN_OK;

}  // GetProblemGeometry

}  // namespace MueLu

#endif  // if defined(HAVE_MUELU_ZOLTAN) && defined(HAVE_MPI)

#endif  // MUELU_ZOLTANINTERFACE_DEF_HPP
