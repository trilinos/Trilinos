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
//                    Jeremie Gaidamour (jngaida@sandia.gov)
//                    Jonathan Hu       (jhu@sandia.gov)
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
#include <Teuchos_DefaultMpiComm.hpp> //TODO: fwd decl.
#include <Teuchos_OpaqueWrapper.hpp>  //TODO: fwd decl.

#include "MueLu_Level.hpp"
#include "MueLu_Exceptions.hpp"
#include "MueLu_Monitor.hpp"

namespace MueLu {


  //-------------------------------------------------------------------------------------------------------------
  // DeclareInput
  //-------------------------------------------------------------------------------------------------------------

  template <class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void ZoltanInterface<LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::
  DeclareInput(Level & currentLevel) const
  {
    Input(currentLevel, "number of partitions");

    Input(currentLevel, "A");
    Input(currentLevel, "Coordinates");

  } //DeclareInput()

  //-------------------------------------------------------------------------------------------------------------
  // Build
  //-------------------------------------------------------------------------------------------------------------

  template <class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void ZoltanInterface<LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::
  Build(Level &level) const
  {
    FactoryMonitor m(*this, "Build", level);

    RCP<Matrix> A = Get< RCP<Matrix> >(level, "A");
    RCP<const Map> rowMap = A->getRowMap();
    GO numPartitions = Get<GO>(level, "number of partitions");

    if (numPartitions == 1) {
      // Running on one processor, so decomposition is the trivial one, all zeros.
      RCP<Xpetra::Vector<GO, LO, GO, NO> > decomposition = Xpetra::VectorFactory<GO, LO, GO, NO>::Build(rowMap, true);
      Set(level, "Partition", decomposition);
      return;
    }

    // Tell Zoltan what kind of local/global IDs we will use.
    // In our case, each GID is two ints and there are no local ids.
    // One can skip this step if the IDs are just single ints.
    RCP<const Teuchos::MpiComm<int> > mpiComm = rcp_dynamic_cast<const Teuchos::MpiComm<int> >(rowMap->getComm());
    float zoltanVersion_;
    Zoltan_Initialize(0, NULL, &zoltanVersion_);
    //TODO define zoltanComm_ as a subcommunicator?!;
    RCP<const Teuchos::OpaqueWrapper<MPI_Comm> > zoltanComm_ = mpiComm->getRawMpiComm();
    RCP<Zoltan> zoltanObj_ = rcp( new Zoltan( (*zoltanComm_)() ) );  //extract the underlying MPI_Comm handle and create a Zoltan object
    if (zoltanObj_==Teuchos::null) throw(Exceptions::RuntimeError("MueLu::Zoltan : Unable to create Zoltan data structure"));
    int rv;
    if ((rv=zoltanObj_->Set_Param("num_gid_entries", "1")) != ZOLTAN_OK )
      throw(Exceptions::RuntimeError("MueLu::Zoltan::Setup : setting parameter 'num_gid_entries' returned error code " + Teuchos::toString(rv)));
    if ( (rv=zoltanObj_->Set_Param("num_lid_entries", "0") ) != ZOLTAN_OK )
      throw(Exceptions::RuntimeError("MueLu::Zoltan::Setup : setting parameter 'num_lid_entries' returned error code " + Teuchos::toString(rv)));
    if ( (rv=zoltanObj_->Set_Param("obj_weight_dim", "1") ) != ZOLTAN_OK )
      throw(Exceptions::RuntimeError("MueLu::Zoltan::Setup : setting parameter 'obj_weight_dim' returned error code " + Teuchos::toString(rv)));

    zoltanObj_->Set_Param("debug_level", "0");

    std::stringstream ss;
    ss << numPartitions;
    zoltanObj_->Set_Param("num_global_partitions", ss.str());

    //TODO: coordinates should be const

    Array<ArrayRCP<SC> > XYZ; // Using this format because no communications needed here. No need for a map and a Xpetra::MultiVector

//     // Build XYZ from XCoordinates, YCoordinates and ZCoordinates
//     if (level.IsAvailable("XCoordinates")) {

//       {
//         XYZ.push_back(level.Get< ArrayRCP<SC> >("XCoordinates"));
//       }

//       if (level.IsAvailable("YCoordinates")) {
//         XYZ.push_back(level.Get< ArrayRCP<SC> >("YCoordinates"));
//       }

//       if (level.IsAvailable("ZCoordinates")) {
//         TEUCHOS_TEST_FOR_EXCEPTION(!level.IsAvailable("YCoordinates"), Exceptions::RuntimeError, "ZCoordinates specified but no YCoordinates");
//         XYZ.push_back(level.Get< ArrayRCP<SC> >("ZCoordinates"));
//       }

//     } else
    if (IsAvailable(level, "Coordinates")) {

      RCP<Matrix> Aloc = Get<RCP<Matrix> >(level, "A");
      LocalOrdinal blksize = Aloc->GetFixedBlockSize();

      RCP<MultiVector> multiVectorXYZ = Get< RCP<MultiVector> >(level, "Coordinates");
      for (int i=0; i< (int)multiVectorXYZ->getNumVectors(); i++) { //FIXME cast
        XYZ.push_back(coalesceCoordinates(multiVectorXYZ->getDataNonConst(i), blksize)); // If blksize == 1, not copy but it's OK to leave 'open' the MultiVector until the destruction of XYZ because no communications using Xpetra
      }

      // TODO: level.Set(XCoordinates / YCoordinates / ZCoordinates as it is computed and might be needed somewhere else. But can wait for now. This code have to be moved anyway.

    } else {
      throw(Exceptions::RuntimeError("MueLu::ZoltanInterface::Build(): no coordinates available"));
    }

    //~~ size_t problemDimension_ = XYZ->getNumVectors();
    size_t problemDimension_ = XYZ.size();

    zoltanObj_->Set_Num_Obj_Fn(GetLocalNumberOfRows, (void *) &*A);
    zoltanObj_->Set_Obj_List_Fn(GetLocalNumberOfNonzeros, (void *) &*A);
    zoltanObj_->Set_Num_Geom_Fn(GetProblemDimension, (void *) &problemDimension_);
    //~~ zoltanObj_->Set_Geom_Multi_Fn(GetProblemGeometry, (void *) &*XYZ);
    zoltanObj_->Set_Geom_Multi_Fn(GetProblemGeometry, (void *) &XYZ);

    // Data pointers that Zoltan requires.
    ZOLTAN_ID_PTR import_gids = NULL;  // Global nums of objs to be imported
    ZOLTAN_ID_PTR import_lids = NULL;  // Local indices to objs to be imported
    int   *import_procs = NULL;        // Proc IDs of procs owning objs to be imported.
    int   *import_to_part = NULL;      // Partition #s to which imported objs should be assigned.
    ZOLTAN_ID_PTR export_gids = NULL;  // Global nums of objs to be exported
    ZOLTAN_ID_PTR export_lids = NULL;  // local indices to objs to be exported
    int   *export_procs = NULL;        // Proc IDs of destination procs for objs to be exported.
    int   *export_to_part = NULL;      // Partition #s for objs to be exported.
    int   num_imported;                // Number of objs to be imported.
    int   num_exported;                // Number of objs to be exported.
    int   newDecomp;                   // Flag indicating whether the decomposition has changed
    int   num_gid_entries;             // Number of array entries in a global ID.
    int   num_lid_entries;

    {
      SubFactoryMonitor m1(*this, "Zoltan RCB", level);
      rv = zoltanObj_->LB_Partition(newDecomp, num_gid_entries, num_lid_entries,
                                    num_imported, import_gids, import_lids, import_procs, import_to_part,
                                    num_exported, export_gids, export_lids, export_procs, export_to_part);
      if (rv == ZOLTAN_FATAL) {
        throw(Exceptions::RuntimeError("Zoltan::LB_Partition() returned error code"));
      }
    }

    //TODO check that A's row map is 1-1.  Zoltan requires this.
    int mypid = rowMap->getComm()->getRank();
    RCP<Xpetra::Vector<GO, LO, GO, NO> > decomposition;

    // Note: if newDecomp == false, decomposition == Teuchos::null
    if (newDecomp) {
      decomposition = Xpetra::VectorFactory<GO, LO, GO, NO>::Build(rowMap, false); //Don't bother initializing, as this will just be overwritten.
      ArrayRCP<GO> decompEntries = decomposition->getDataNonConst(0);
      for (typename ArrayRCP<GO>::iterator i = decompEntries.begin(); i != decompEntries.end(); ++i)
        *i = mypid;
      LocalOrdinal blockSize = A->GetFixedBlockSize();
      for (int i=0; i< num_exported; ++i) {
        LO localEl = rowMap->getLocalElement(export_gids[i]);
        int partNum = export_to_part[i];
        for (LO j=0; j<blockSize; ++j)
          decompEntries[ localEl + j ] = partNum;
          //decompEntries[ rowMap->getLocalElement(export_gids[i]) + j ] = export_to_part[i];
      }
    }

    Set(level, "Partition", decomposition);

    zoltanObj_->LB_Free_Part(&import_gids, &import_lids, &import_procs, &import_to_part);
    zoltanObj_->LB_Free_Part(&export_gids, &export_lids, &export_procs, &export_to_part);

  } //Build()

  //-------------------------------------------------------------------------------------------------------------
  // GetLocalNumberOfRows
  //-------------------------------------------------------------------------------------------------------------

  template <class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  int ZoltanInterface<LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::
  GetLocalNumberOfRows(void *data, int *ierr)
  {
    if (data == NULL) {
      *ierr = ZOLTAN_FATAL;
      return -1;
    }
    *ierr = ZOLTAN_OK;
    //TODO is there a safer way to cast?
    //Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> *A = (Matrix*) data;
    Matrix *A = (Matrix*) data;
    LocalOrdinal blockSize = A->GetFixedBlockSize(); //FIXME
    if (blockSize==0) throw(Exceptions::RuntimeError("MueLu::Zoltan : Matrix has block size 0."));
    return (A->getRowMap()->getNodeNumElements() / blockSize); //FIXME
  } //GetLocalNumberOfRows()

  //-------------------------------------------------------------------------------------------------------------
  // GetLocalNumberOfNonzeros
  //-------------------------------------------------------------------------------------------------------------

  template <class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void ZoltanInterface<LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::
  GetLocalNumberOfNonzeros(void *data, int NumGidEntries, int NumLidEntries, ZOLTAN_ID_PTR gids,
                           ZOLTAN_ID_PTR lids, int wgtDim, float *weights, int *ierr)
  {
    if (data == NULL || NumGidEntries < 1) {
      *ierr = ZOLTAN_FATAL;
      return;
    } else {
      *ierr = ZOLTAN_OK;
    }

    //TODO is there a safer way to cast?
    Matrix *A = (Matrix*) data;
    RCP<const Map> map = A->getRowMap();
    Teuchos::ArrayView<const LO> cols;
    Teuchos::ArrayView<const SC> vals;
    LocalOrdinal blockSize = A->GetFixedBlockSize(); //FIXME
    if (blockSize==0) throw(Exceptions::RuntimeError("MueLu::Zoltan : Matrix has block size 0."));
    if (blockSize == 1) {
      for (size_t i=0; i<map->getNodeNumElements(); ++i) {
        gids[i] = (ZOLTAN_ID_TYPE) map->getGlobalElement(i);
        A->getLocalRowView(i, cols, vals);
        weights[i] = cols.size();
      }
    } else {
      LocalOrdinal numBlocks = A->getRowMap()->getNodeNumElements() / blockSize;
      for (LocalOrdinal i=0; i<numBlocks; ++i) {
        gids[i] = (ZOLTAN_ID_TYPE) map->getGlobalElement(i*blockSize);
        LO nnz=0;
        for (LocalOrdinal j=i*blockSize; j<(i+1)*blockSize; ++j) {
          A->getLocalRowView(j, cols, vals);
          nnz += vals.size();
        }
        weights[i] = nnz;
      } //for (LocalOrdinal i=0; i<numBlocks; ++i)
    }

  } //GetLocalNumberOfNonzeros()

  //-------------------------------------------------------------------------------------------------------------
  // GetProblemDimension
  //-------------------------------------------------------------------------------------------------------------

  template <class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  int ZoltanInterface<LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::
  GetProblemDimension(void *data, int *ierr)
  {
    //TODO is there a safer way to cast?
    int dim = *((int*)data);
    *ierr = ZOLTAN_OK; /* set error flag */
    return(dim);
  } //GetProblemDimension

  //-------------------------------------------------------------------------------------------------------------
  // GetProblemGeometry
  //-------------------------------------------------------------------------------------------------------------

  template <class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void ZoltanInterface<LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::
  GetProblemGeometry(void *data, int numGIDEntries, int numLIDEntries, int numObjectIDs,
                     ZOLTAN_ID_PTR gids, ZOLTAN_ID_PTR lids, int dim, double *coordinates, int *ierr)
  {
    if (data == NULL) {
      *ierr = ZOLTAN_FATAL;
      return;
    }

    //TODO is there a safer way to cast?
    //~~ MultiVector *XYZ = (MultiVector*) data;
    Array<ArrayRCP<SC> > * XYZpt = (Array<ArrayRCP<SC> > *)data;
    Array<ArrayRCP<SC> > & XYZ   = *XYZpt;

    //~~ if (dim != (int) XYZ->getNumVectors()) {
    if (dim != (int) XYZ.size()) { //FIXME: cast to size_t instead?
      //FIXME I'm assuming dim should be 1, 2, or 3 coming in?!
      *ierr = ZOLTAN_FATAL;
      return;
    }

    //~ assert(numObjectIDs == XYZ->getLocalLength());
    for(int j=0; j<dim; j++) {
      assert(numObjectIDs == XYZ[j].size()); //FIXME: TEST_FOR_EXCEPTION instead?
    }

    /*~~
    ArrayRCP<ArrayRCP<const SC> > XYZdata(dim);
    for (int j=0; j<dim; ++j) XYZdata[j] = XYZ->getData(j);
    for (size_t i=0; i<XYZ->getLocalLength(); ++i) {
      for (int j=0; j<dim; ++j) {
        coordinates[i*dim+j] = (double) XYZdata[j][i];
      }
    }
    */

    for (size_t i=0; i<(size_t)XYZ[0].size(); ++i) { //FIXME cast OK?
      for (int j=0; j<dim; ++j) {
        coordinates[i*dim+j] = (double) XYZ[j][i];
      }
    }

    *ierr = ZOLTAN_OK;

  } //GetProblemGeometry


  template <class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  ArrayRCP<double> ZoltanInterface<LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::coalesceCoordinates(ArrayRCP<double> coord, LocalOrdinal blksize) {
    if (blksize == 1)
      return coord;

    ArrayRCP<double> coalesceCoord(coord.size()/blksize); //TODO: how to avoid automatic initialization of the vector? using arcp()?

    for(int i=0; i<coord.size(); i++) {
#define myDEBUG
#ifdef myDEBUG //FIXME-> HAVE_MUELU_DEBUG
      for(int j=1; j < blksize; j++) {
        TEUCHOS_TEST_FOR_EXCEPTION(coord[i*blksize + j] != coord[i*blksize], Exceptions::RuntimeError, "MueLu::ZoltanInterface: coalesceCoord problem");
      }
#endif
      coalesceCoord[i] = coalesceCoord[i*blksize];
    }

    //std::cout << coord << std::endl;
    //std::cout << coalesceCoord << std::endl;

    return coalesceCoord;
  }

} //namespace MueLu

#endif //if defined(HAVE_MUELU_ZOLTAN) && defined(HAVE_MPI)

#endif // MUELU_ZOLTANINTERFACE_DEF_HPP
