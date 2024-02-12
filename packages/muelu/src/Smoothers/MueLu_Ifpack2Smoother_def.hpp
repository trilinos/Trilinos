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
#ifndef MUELU_IFPACK2SMOOTHER_DEF_HPP
#define MUELU_IFPACK2SMOOTHER_DEF_HPP

#include "MueLu_ConfigDefs.hpp"

#if defined(HAVE_MUELU_IFPACK2)

#include <Teuchos_ParameterList.hpp>

#include <Tpetra_RowMatrix.hpp>

#include <Ifpack2_Chebyshev.hpp>
#include <Ifpack2_Hiptmair.hpp>
#include <Ifpack2_RILUK.hpp>
#include <Ifpack2_Relaxation.hpp>
#include <Ifpack2_ILUT.hpp>
#include <Ifpack2_BlockRelaxation.hpp>
#include <Ifpack2_Factory.hpp>
#include <Ifpack2_Parameters.hpp>

#include <Xpetra_BlockedCrsMatrix.hpp>
#include <Xpetra_CrsMatrix.hpp>
#include <Xpetra_CrsMatrixWrap.hpp>
#include <Xpetra_TpetraBlockCrsMatrix.hpp>
#include <Xpetra_Matrix.hpp>
#include <Xpetra_MatrixMatrix.hpp>
#include <Xpetra_MultiVectorFactory.hpp>
#include <Xpetra_TpetraMultiVector.hpp>

#include <Tpetra_BlockCrsMatrix_Helpers.hpp>

#include "MueLu_Ifpack2Smoother_decl.hpp"
#include "MueLu_Level.hpp"
#include "MueLu_Utilities.hpp"
#include "MueLu_Monitor.hpp"
#include "MueLu_Aggregates.hpp"

#ifdef HAVE_MUELU_INTREPID2
#include "MueLu_IntrepidPCoarsenFactory_decl.hpp"
#include "MueLu_IntrepidPCoarsenFactory_def.hpp"
#include "Intrepid2_Basis.hpp"
#include "Kokkos_DynRankView.hpp"
#endif

namespace MueLu {

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
Ifpack2Smoother<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Ifpack2Smoother(const std::string& type, const Teuchos::ParameterList& paramList, const LO& overlap)
  : type_(type)
  , overlap_(overlap) {
  typedef Tpetra::RowMatrix<SC, LO, GO, NO> tRowMatrix;
  bool isSupported = Ifpack2::Factory::isSupported<tRowMatrix>(type_) || (type_ == "LINESMOOTHING_TRIDI_RELAXATION" ||
                                                                          type_ == "LINESMOOTHING_TRIDI RELAXATION" ||
                                                                          type_ == "LINESMOOTHING_TRIDIRELAXATION" ||
                                                                          type_ == "LINESMOOTHING_TRIDIAGONAL_RELAXATION" ||
                                                                          type_ == "LINESMOOTHING_TRIDIAGONAL RELAXATION" ||
                                                                          type_ == "LINESMOOTHING_TRIDIAGONALRELAXATION" ||
                                                                          type_ == "LINESMOOTHING_BANDED_RELAXATION" ||
                                                                          type_ == "LINESMOOTHING_BANDED RELAXATION" ||
                                                                          type_ == "LINESMOOTHING_BANDEDRELAXATION" ||
                                                                          type_ == "LINESMOOTHING_BLOCK_RELAXATION" ||
                                                                          type_ == "LINESMOOTHING_BLOCK RELAXATION" ||
                                                                          type_ == "LINESMOOTHING_BLOCKRELAXATION" ||
                                                                          type_ == "TOPOLOGICAL" ||
                                                                          type_ == "AGGREGATE");
  this->declareConstructionOutcome(!isSupported, "Ifpack2 does not provide the smoother '" + type_ + "'.");
  if (isSupported)
    SetParameterList(paramList);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void Ifpack2Smoother<Scalar, LocalOrdinal, GlobalOrdinal, Node>::SetParameterList(const Teuchos::ParameterList& paramList) {
  Factory::SetParameterList(paramList);

  if (SmootherPrototype::IsSetup()) {
    // It might be invalid to change parameters after the setup, but it depends entirely on Ifpack implementation.
    // TODO: I don't know if Ifpack returns an error code or exception or ignore parameters modification in this case...
    SetPrecParameters();
  }
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void Ifpack2Smoother<Scalar, LocalOrdinal, GlobalOrdinal, Node>::SetPrecParameters(const Teuchos::ParameterList& list) const {
  std::string prefix       = this->ShortClassName() + ": SetPrecParameters";
  RCP<TimeMonitor> tM      = rcp(new TimeMonitor(*this, prefix, Timings0));
  ParameterList& paramList = const_cast<ParameterList&>(this->GetParameterList());

  paramList.setParameters(list);

  RCP<ParameterList> precList = this->RemoveFactoriesFromList(this->GetParameterList());

  // Do we want an Ifpack2 apply timer?
  precList->set("timer for apply", this->IsPrint(Timings));

  if (!precList.is_null() && precList->isParameter("partitioner: type") && precList->get<std::string>("partitioner: type") == "linear" &&
      !precList->isParameter("partitioner: local parts")) {
    LO matrixBlockSize = 1;
    int lclSize        = A_->getRangeMap()->getLocalNumElements();
    RCP<Matrix> matA   = rcp_dynamic_cast<Matrix>(A_);
    if (!matA.is_null()) {
      lclSize         = matA->getLocalNumRows();
      matrixBlockSize = matA->GetFixedBlockSize();
    }
    precList->set("partitioner: local parts", lclSize / matrixBlockSize);
  }

  prec_->setParameters(*precList);

  paramList.setParameters(*precList);  // what about that??
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void Ifpack2Smoother<Scalar, LocalOrdinal, GlobalOrdinal, Node>::DeclareInput(Level& currentLevel) const {
  this->Input(currentLevel, "A");

  if (type_ == "LINESMOOTHING_TRIDI_RELAXATION" ||
      type_ == "LINESMOOTHING_TRIDI RELAXATION" ||
      type_ == "LINESMOOTHING_TRIDIRELAXATION" ||
      type_ == "LINESMOOTHING_TRIDIAGONAL_RELAXATION" ||
      type_ == "LINESMOOTHING_TRIDIAGONAL RELAXATION" ||
      type_ == "LINESMOOTHING_TRIDIAGONALRELAXATION" ||
      type_ == "LINESMOOTHING_BANDED_RELAXATION" ||
      type_ == "LINESMOOTHING_BANDED RELAXATION" ||
      type_ == "LINESMOOTHING_BANDEDRELAXATION" ||
      type_ == "LINESMOOTHING_BLOCK_RELAXATION" ||
      type_ == "LINESMOOTHING_BLOCK RELAXATION" ||
      type_ == "LINESMOOTHING_BLOCKRELAXATION") {
    this->Input(currentLevel, "CoarseNumZLayers");           // necessary for fallback criterion
    this->Input(currentLevel, "LineDetection_VertLineIds");  // necessary to feed block smoother
  } else if (type_ == "BLOCK RELAXATION" ||
             type_ == "BLOCK_RELAXATION" ||
             type_ == "BLOCKRELAXATION" ||
             // Banded
             type_ == "BANDED_RELAXATION" ||
             type_ == "BANDED RELAXATION" ||
             type_ == "BANDEDRELAXATION" ||
             // Tridiagonal
             type_ == "TRIDI_RELAXATION" ||
             type_ == "TRIDI RELAXATION" ||
             type_ == "TRIDIRELAXATION" ||
             type_ == "TRIDIAGONAL_RELAXATION" ||
             type_ == "TRIDIAGONAL RELAXATION" ||
             type_ == "TRIDIAGONALRELAXATION") {
    // We need to check for the "partitioner type" = "line"
    ParameterList precList = this->GetParameterList();
    if (precList.isParameter("partitioner: type") &&
        precList.get<std::string>("partitioner: type") == "line") {
      this->Input(currentLevel, "Coordinates");
    }
  } else if (type_ == "TOPOLOGICAL") {
    // for the topological smoother, we require an element to node map:
    this->Input(currentLevel, "pcoarsen: element to node map");
  } else if (type_ == "AGGREGATE") {
    // Aggregate smoothing needs aggregates
    this->Input(currentLevel, "Aggregates");
  } else if (type_ == "HIPTMAIR") {
    // Hiptmair needs D0 and NodeMatrix
    this->Input(currentLevel, "NodeMatrix");
    this->Input(currentLevel, "D0");
  }
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void Ifpack2Smoother<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Setup(Level& currentLevel) {
  FactoryMonitor m(*this, "Setup Smoother", currentLevel);
  A_                       = Factory::Get<RCP<Operator>>(currentLevel, "A");
  ParameterList& paramList = const_cast<ParameterList&>(this->GetParameterList());

  // If the user asked us to convert the matrix into BlockCrsMatrix form, we do that now

  if (paramList.isParameter("smoother: use blockcrsmatrix storage") && paramList.get<bool>("smoother: use blockcrsmatrix storage")) {
    int blocksize    = 1;
    RCP<Matrix> matA = rcp_dynamic_cast<Matrix>(A_);
    if (!matA.is_null())
      blocksize = matA->GetFixedBlockSize();
    if (blocksize) {
      // NOTE: Don't think you can move this out of the if block.  You can't. The test MueLu_MeshTyingBlocked_SimpleSmoother_2dof_medium_MPI_1 will fail

      RCP<CrsMatrixWrap> AcrsWrap = rcp_dynamic_cast<CrsMatrixWrap>(A_);
      if (AcrsWrap.is_null())
        throw std::runtime_error("Ifpack2Smoother: Cannot convert matrix A to CrsMatrixWrap object.");
      RCP<CrsMatrix> Acrs = AcrsWrap->getCrsMatrix();
      if (Acrs.is_null())
        throw std::runtime_error("Ifpack2Smoother: Cannot extract CrsMatrix from matrix A.");
      RCP<TpetraCrsMatrix> At = rcp_dynamic_cast<TpetraCrsMatrix>(Acrs);
      if (At.is_null()) {
        if (!Xpetra::Helpers<Scalar, LO, GO, Node>::isTpetraBlockCrs(matA))
          throw std::runtime_error("Ifpack2Smoother: Cannot extract CrsMatrix or BlockCrsMatrix from matrix A.");
        this->GetOStream(Statistics0) << "Ifpack2Smoother: Using (native) BlockCrsMatrix storage with blocksize " << blocksize << std::endl;
        paramList.remove("smoother: use blockcrsmatrix storage");
      } else {
        RCP<Tpetra::BlockCrsMatrix<Scalar, LO, GO, Node>> blockCrs = Tpetra::convertToBlockCrsMatrix(*At->getTpetra_CrsMatrix(), blocksize);
        RCP<CrsMatrix> blockCrs_as_crs                             = rcp(new TpetraBlockCrsMatrix(blockCrs));
        RCP<CrsMatrixWrap> blockWrap                               = rcp(new CrsMatrixWrap(blockCrs_as_crs));
        A_                                                         = blockWrap;
        this->GetOStream(Statistics0) << "Ifpack2Smoother: Using BlockCrsMatrix storage with blocksize " << blocksize << std::endl;

        paramList.remove("smoother: use blockcrsmatrix storage");
      }
    }
  }

  if (type_ == "SCHWARZ")
    SetupSchwarz(currentLevel);

  else if (type_ == "LINESMOOTHING_TRIDI_RELAXATION" ||
           type_ == "LINESMOOTHING_TRIDI RELAXATION" ||
           type_ == "LINESMOOTHING_TRIDIRELAXATION" ||
           type_ == "LINESMOOTHING_TRIDIAGONAL_RELAXATION" ||
           type_ == "LINESMOOTHING_TRIDIAGONAL RELAXATION" ||
           type_ == "LINESMOOTHING_TRIDIAGONALRELAXATION" ||
           type_ == "LINESMOOTHING_BANDED_RELAXATION" ||
           type_ == "LINESMOOTHING_BANDED RELAXATION" ||
           type_ == "LINESMOOTHING_BANDEDRELAXATION" ||
           type_ == "LINESMOOTHING_BLOCK_RELAXATION" ||
           type_ == "LINESMOOTHING_BLOCK RELAXATION" ||
           type_ == "LINESMOOTHING_BLOCKRELAXATION")
    SetupLineSmoothing(currentLevel);

  else if (type_ == "BLOCK_RELAXATION" ||
           type_ == "BLOCK RELAXATION" ||
           type_ == "BLOCKRELAXATION" ||
           // Banded
           type_ == "BANDED_RELAXATION" ||
           type_ == "BANDED RELAXATION" ||
           type_ == "BANDEDRELAXATION" ||
           // Tridiagonal
           type_ == "TRIDI_RELAXATION" ||
           type_ == "TRIDI RELAXATION" ||
           type_ == "TRIDIRELAXATION" ||
           type_ == "TRIDIAGONAL_RELAXATION" ||
           type_ == "TRIDIAGONAL RELAXATION" ||
           type_ == "TRIDIAGONALRELAXATION")
    SetupBlockRelaxation(currentLevel);

  else if (type_ == "CHEBYSHEV")
    SetupChebyshev(currentLevel);

  else if (type_ == "TOPOLOGICAL") {
#ifdef HAVE_MUELU_INTREPID2
    SetupTopological(currentLevel);
#else
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::invalid_argument, "'TOPOLOGICAL' smoother choice requires Intrepid2");
#endif
  } else if (type_ == "AGGREGATE")
    SetupAggregate(currentLevel);

  else if (type_ == "HIPTMAIR")
    SetupHiptmair(currentLevel);

  else {
    SetupGeneric(currentLevel);
  }

  SmootherPrototype::IsSetup(true);

  this->GetOStream(Statistics1) << description() << std::endl;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void Ifpack2Smoother<Scalar, LocalOrdinal, GlobalOrdinal, Node>::SetupSchwarz(Level& /* currentLevel */) {
  typedef Tpetra::RowMatrix<SC, LO, GO, NO> tRowMatrix;

  bool reusePreconditioner = false;
  if (this->IsSetup() == true) {
    // Reuse the constructed preconditioner
    this->GetOStream(Runtime1) << "MueLu::Ifpack2Smoother::SetupSchwarz(): Setup() has already been called, assuming reuse" << std::endl;

    bool isTRowMatrix = true;
    RCP<const tRowMatrix> tA;
    try {
      tA = Utilities::Op2NonConstTpetraRow(A_);
    } catch (Exceptions::BadCast&) {
      isTRowMatrix = false;
    }

    RCP<Ifpack2::Details::CanChangeMatrix<tRowMatrix>> prec = rcp_dynamic_cast<Ifpack2::Details::CanChangeMatrix<tRowMatrix>>(prec_);
    if (!prec.is_null() && isTRowMatrix) {
      prec->setMatrix(tA);
      reusePreconditioner = true;
    } else {
      this->GetOStream(Warnings0) << "MueLu::Ifpack2Smoother::SetupSchwarz(): reuse of this type is not available "
                                     "(either failed cast to CanChangeMatrix, or to Tpetra Row Matrix), reverting to full construction"
                                  << std::endl;
    }
  }

  if (!reusePreconditioner) {
    ParameterList& paramList = const_cast<ParameterList&>(this->GetParameterList());

    bool isBlockedMatrix = false;
    RCP<Matrix> merged2Mat;

    std::string sublistName = "subdomain solver parameters";
    if (paramList.isSublist(sublistName)) {
      // If we are doing "user" partitioning, we assume that what the user
      // really wants to do is make tiny little subdomains with one row
      // assigned to each subdomain. The rows used for these little
      // subdomains correspond to those in the 2nd block row. Then,
      // if we overlap these mini-subdomains, we will do something that
      // looks like Vanka (grabbing all velocities associated with each
      // each pressure unknown). In addition, we put all Dirichlet points
      // as a little mini-domain.
      ParameterList& subList = paramList.sublist(sublistName);

      std::string partName = "partitioner: type";
      // Pretty sure no one has been using this. Unfortunately, old if
      // statement (which checked for equality with "user") prevented
      // anyone from employing other types of Ifpack2 user partition
      // options. Leaving this and switching if to "vanka user" just
      // in case some day someone might want to use this.
      if (subList.isParameter(partName) && subList.get<std::string>(partName) == "vanka user") {
        isBlockedMatrix = true;

        RCP<BlockedCrsMatrix> bA = rcp_dynamic_cast<BlockedCrsMatrix>(A_);
        TEUCHOS_TEST_FOR_EXCEPTION(bA.is_null(), Exceptions::BadCast,
                                   "Matrix A must be of type BlockedCrsMatrix.");

        size_t numVels = bA->getMatrix(0, 0)->getLocalNumRows();
        size_t numPres = bA->getMatrix(1, 0)->getLocalNumRows();
        size_t numRows = rcp_dynamic_cast<Matrix>(A_, true)->getLocalNumRows();

        ArrayRCP<LocalOrdinal> blockSeeds(numRows, Teuchos::OrdinalTraits<LocalOrdinal>::invalid());

        size_t numBlocks = 0;
        for (size_t rowOfB = numVels; rowOfB < numVels + numPres; ++rowOfB)
          blockSeeds[rowOfB] = numBlocks++;

        RCP<BlockedCrsMatrix> bA2 = rcp_dynamic_cast<BlockedCrsMatrix>(A_);
        TEUCHOS_TEST_FOR_EXCEPTION(bA2.is_null(), Exceptions::BadCast,
                                   "Matrix A must be of type BlockedCrsMatrix.");

        merged2Mat = bA2->Merge();

        // Add Dirichlet rows to the list of seeds
        ArrayRCP<const bool> boundaryNodes = Utilities::DetectDirichletRows(*merged2Mat, 0.0);
        bool haveBoundary                  = false;
        for (LO i = 0; i < boundaryNodes.size(); i++)
          if (boundaryNodes[i]) {
            // FIXME:
            // 1. would not this [] overlap with some in the previos blockSeed loop?
            // 2. do we need to distinguish between pressure and velocity Dirichlet b.c.
            blockSeeds[i] = numBlocks;
            haveBoundary  = true;
          }
        if (haveBoundary)
          numBlocks++;

        subList.set("partitioner: type", "user");
        subList.set("partitioner: map", blockSeeds);
        subList.set("partitioner: local parts", as<int>(numBlocks));

      } else {
        RCP<BlockedCrsMatrix> bA = rcp_dynamic_cast<BlockedCrsMatrix>(A_);
        if (!bA.is_null()) {
          isBlockedMatrix = true;
          merged2Mat      = bA->Merge();
        }
      }
    }

    RCP<const tRowMatrix> tA;
    if (isBlockedMatrix == true)
      tA = Utilities::Op2NonConstTpetraRow(merged2Mat);
    else
      tA = Utilities::Op2NonConstTpetraRow(A_);

    prec_ = Ifpack2::Factory::create(type_, tA, overlap_);
    SetPrecParameters();

    prec_->initialize();
  }

  prec_->compute();
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void Ifpack2Smoother<Scalar, LocalOrdinal, GlobalOrdinal, Node>::SetupAggregate(Level& currentLevel) {
  ParameterList& paramList = const_cast<ParameterList&>(this->GetParameterList());

  if (this->IsSetup() == true) {
    this->GetOStream(Warnings0) << "MueLu::Ifpack2Smoother::SetupAggregate(): Setup() has already been called" << std::endl;
    this->GetOStream(Warnings0) << "MueLu::Ifpack2Smoother::SetupAggregate(): reuse of this type is not available, reverting to full construction" << std::endl;
  }

  this->GetOStream(Statistics0) << "Ifpack2Smoother: Using Aggregate Smoothing" << std::endl;

  RCP<Aggregates> aggregates = Factory::Get<RCP<Aggregates>>(currentLevel, "Aggregates");

  RCP<const LOMultiVector> vertex2AggId = aggregates->GetVertex2AggId();
  ArrayRCP<LO> aggregate_ids            = rcp_const_cast<LOMultiVector>(vertex2AggId)->getDataNonConst(0);
  ArrayRCP<LO> dof_ids;

  // We need to unamalgamate, if the FixedBlockSize > 1
  LO blocksize     = 1;
  RCP<Matrix> matA = rcp_dynamic_cast<Matrix>(A_);
  if (!matA.is_null())
    blocksize = matA->GetFixedBlockSize();
  if (blocksize > 1) {
    dof_ids.resize(aggregate_ids.size() * blocksize);
    for (LO i = 0; i < (LO)aggregate_ids.size(); i++) {
      for (LO j = 0; j < (LO)blocksize; j++)
        dof_ids[i * blocksize + j] = aggregate_ids[i];
    }
  } else {
    dof_ids = aggregate_ids;
  }

  paramList.set("partitioner: map", dof_ids);
  paramList.set("partitioner: type", "user");
  paramList.set("partitioner: overlap", 0);
  paramList.set("partitioner: local parts", (int)aggregates->GetNumAggregates());

  RCP<const Tpetra::RowMatrix<SC, LO, GO, NO>> tA = Utilities::Op2NonConstTpetraRow(A_);

  type_ = "BLOCKRELAXATION";
  prec_ = Ifpack2::Factory::create(type_, tA, overlap_);
  SetPrecParameters();
  prec_->initialize();
  prec_->compute();
}

#ifdef HAVE_MUELU_INTREPID2
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void Ifpack2Smoother<Scalar, LocalOrdinal, GlobalOrdinal, Node>::SetupTopological(Level& currentLevel) {
  /*

   basic notion:

   Look for user input indicating topo dimension, something like "topological domain type: {node|edge|face}"
   Call something like what you can find in Poisson example line 1180 to set seeds for a smoother

   */
  if (this->IsSetup() == true) {
    this->GetOStream(Warnings0) << "MueLu::Ifpack2Smoother::SetupTopological(): Setup() has already been called" << std::endl;
    this->GetOStream(Warnings0) << "MueLu::Ifpack2Smoother::SetupTopological(): reuse of this type is not available, reverting to full construction" << std::endl;
  }

  ParameterList& paramList = const_cast<ParameterList&>(this->GetParameterList());

  typedef typename Node::device_type::execution_space ES;

  typedef Kokkos::DynRankView<LocalOrdinal, typename Node::device_type> FCO;  //

  LocalOrdinal lo_invalid = Teuchos::OrdinalTraits<LO>::invalid();

  using namespace std;

  const Teuchos::RCP<FCO> elemToNode = Factory::Get<Teuchos::RCP<FCO>>(currentLevel, "pcoarsen: element to node map");

  string basisString = paramList.get<string>("pcoarsen: hi basis");
  int degree;
  // NOTE: To make sure Stokhos works we only instantiate these guys with double.  There's a lot
  // of stuff in the guts of Intrepid2 that doesn't play well with Stokhos as of yet.  Here, we only
  // care about the assignment of basis ordinals to topological entities, so this code is actually
  // independent of the Scalar type--hard-coding double here won't hurt us.
  auto basis = MueLuIntrepid::BasisFactory<double, ES>(basisString, degree);

  string topologyTypeString = paramList.get<string>("smoother: neighborhood type");
  int dimension;
  if (topologyTypeString == "node")
    dimension = 0;
  else if (topologyTypeString == "edge")
    dimension = 1;
  else if (topologyTypeString == "face")
    dimension = 2;
  else if (topologyTypeString == "cell")
    dimension = basis->getBaseCellTopology().getDimension();
  else
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::invalid_argument, "Unrecognized smoother neighborhood type.  Supported types are node, edge, face.");
  RCP<Matrix> matA = rcp_dynamic_cast<Matrix>(A_, true);
  vector<vector<LocalOrdinal>> seeds;
  MueLuIntrepid::FindGeometricSeedOrdinals(basis, *elemToNode, seeds, *matA->getRowMap(), *matA->getColMap());

  // Ifpack2 wants the seeds in an array of the same length as the number of local elements,
  // with local partition #s marked for the ones that are seeds, and invalid for the rest
  int myNodeCount = matA->getRowMap()->getLocalNumElements();
  ArrayRCP<LocalOrdinal> nodeSeeds(myNodeCount, lo_invalid);
  int localPartitionNumber = 0;
  for (LocalOrdinal seed : seeds[dimension]) {
    nodeSeeds[seed] = localPartitionNumber++;
  }

  paramList.remove("smoother: neighborhood type");
  paramList.remove("pcoarsen: hi basis");

  paramList.set("partitioner: map", nodeSeeds);
  paramList.set("partitioner: type", "user");
  paramList.set("partitioner: overlap", 1);
  paramList.set("partitioner: local parts", int(seeds[dimension].size()));

  RCP<const Tpetra::RowMatrix<SC, LO, GO, NO>> tA = Utilities::Op2NonConstTpetraRow(A_);

  type_ = "BLOCKRELAXATION";
  prec_ = Ifpack2::Factory::create(type_, tA, overlap_);
  SetPrecParameters();
  prec_->initialize();
  prec_->compute();
}
#endif

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void Ifpack2Smoother<Scalar, LocalOrdinal, GlobalOrdinal, Node>::SetupLineSmoothing(Level& currentLevel) {
  if (this->IsSetup() == true) {
    this->GetOStream(Warnings0) << "MueLu::Ifpack2Smoother::SetupLineSmoothing(): Setup() has already been called" << std::endl;
    this->GetOStream(Warnings0) << "MueLu::Ifpack2Smoother::SetupLineSmoothing(): reuse of this type is not available, reverting to full construction" << std::endl;
  }

  ParameterList& myparamList = const_cast<ParameterList&>(this->GetParameterList());

  LO CoarseNumZLayers = Factory::Get<LO>(currentLevel, "CoarseNumZLayers");
  if (CoarseNumZLayers > 0) {
    Teuchos::ArrayRCP<LO> TVertLineIdSmoo = Factory::Get<Teuchos::ArrayRCP<LO>>(currentLevel, "LineDetection_VertLineIds");

    // determine number of local parts
    LO maxPart = 0;
    for (size_t k = 0; k < Teuchos::as<size_t>(TVertLineIdSmoo.size()); k++) {
      if (maxPart < TVertLineIdSmoo[k]) maxPart = TVertLineIdSmoo[k];
    }
    RCP<Matrix> matA    = rcp_dynamic_cast<Matrix>(A_, true);
    size_t numLocalRows = matA->getLocalNumRows();

    TEUCHOS_TEST_FOR_EXCEPTION(numLocalRows % TVertLineIdSmoo.size() != 0, Exceptions::RuntimeError,
                               "MueLu::Ifpack2Smoother::Setup(): the number of local nodes is incompatible with the TVertLineIdsSmoo.");

    // actualDofsPerNode is the actual number of matrix rows per mesh element.
    // It is encoded in either the MueLu Level, or in the Xpetra matrix block size.
    // This value is needed by Ifpack2 to do decoupled block relaxation.
    int actualDofsPerNode = numLocalRows / TVertLineIdSmoo.size();
    LO matrixBlockSize    = matA->GetFixedBlockSize();
    if (matrixBlockSize > 1 && actualDofsPerNode > 1) {
      TEUCHOS_TEST_FOR_EXCEPTION(actualDofsPerNode != matrixBlockSize, Exceptions::RuntimeError,
                                 "MueLu::Ifpack2Smoother::Setup(): A is a block matrix but its block size and DOFs/node from partitioner disagree");
    } else if (matrixBlockSize > 1) {
      actualDofsPerNode = matrixBlockSize;
    }
    myparamList.set("partitioner: PDE equations", actualDofsPerNode);

    if (numLocalRows == Teuchos::as<size_t>(TVertLineIdSmoo.size())) {
      myparamList.set("partitioner: type", "user");
      myparamList.set("partitioner: map", TVertLineIdSmoo);
      myparamList.set("partitioner: local parts", maxPart + 1);
    } else {
      // we assume a constant number of DOFs per node
      size_t numDofsPerNode = numLocalRows / TVertLineIdSmoo.size();

      // Create a new Teuchos::ArrayRCP<LO> of size numLocalRows and fill it with the corresponding information
      Teuchos::ArrayRCP<LO> partitionerMap(numLocalRows, Teuchos::OrdinalTraits<LocalOrdinal>::invalid());
      for (size_t blockRow = 0; blockRow < Teuchos::as<size_t>(TVertLineIdSmoo.size()); ++blockRow)
        for (size_t dof = 0; dof < numDofsPerNode; dof++)
          partitionerMap[blockRow * numDofsPerNode + dof] = TVertLineIdSmoo[blockRow];
      myparamList.set("partitioner: type", "user");
      myparamList.set("partitioner: map", partitionerMap);
      myparamList.set("partitioner: local parts", maxPart + 1);
    }

    if (type_ == "LINESMOOTHING_BANDED_RELAXATION" ||
        type_ == "LINESMOOTHING_BANDED RELAXATION" ||
        type_ == "LINESMOOTHING_BANDEDRELAXATION")
      type_ = "BANDEDRELAXATION";
    else if (type_ == "LINESMOOTHING_TRIDI_RELAXATION" ||
             type_ == "LINESMOOTHING_TRIDI RELAXATION" ||
             type_ == "LINESMOOTHING_TRIDIRELAXATION" ||
             type_ == "LINESMOOTHING_TRIDIAGONAL_RELAXATION" ||
             type_ == "LINESMOOTHING_TRIDIAGONAL RELAXATION" ||
             type_ == "LINESMOOTHING_TRIDIAGONALRELAXATION")
      type_ = "TRIDIAGONALRELAXATION";
    else
      type_ = "BLOCKRELAXATION";
  } else {
    // line detection failed -> fallback to point-wise relaxation
    this->GetOStream(Runtime0) << "Line detection failed: fall back to point-wise relaxation" << std::endl;
    myparamList.remove("partitioner: type", false);
    myparamList.remove("partitioner: map", false);
    myparamList.remove("partitioner: local parts", false);
    type_ = "RELAXATION";
  }

  RCP<const Tpetra::RowMatrix<SC, LO, GO, NO>> tA = Utilities::Op2NonConstTpetraRow(A_);

  prec_ = Ifpack2::Factory::create(type_, tA, overlap_);
  SetPrecParameters();
  prec_->initialize();
  prec_->compute();
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void Ifpack2Smoother<Scalar, LocalOrdinal, GlobalOrdinal, Node>::SetupBlockRelaxation(Level& currentLevel) {
  typedef Tpetra::RowMatrix<SC, LO, GO, NO> tRowMatrix;

  RCP<BlockedCrsMatrix> bA = rcp_dynamic_cast<BlockedCrsMatrix>(A_);
  if (!bA.is_null())
    A_ = bA->Merge();

  RCP<const tRowMatrix> tA = Utilities::Op2NonConstTpetraRow(A_);

  bool reusePreconditioner = false;
  if (this->IsSetup() == true) {
    // Reuse the constructed preconditioner
    this->GetOStream(Runtime1) << "MueLu::Ifpack2Smoother::SetupBlockRelaxation(): Setup() has already been called, assuming reuse" << std::endl;

    RCP<Ifpack2::Details::CanChangeMatrix<tRowMatrix>> prec = rcp_dynamic_cast<Ifpack2::Details::CanChangeMatrix<tRowMatrix>>(prec_);
    if (!prec.is_null()) {
      prec->setMatrix(tA);
      reusePreconditioner = true;
    } else {
      this->GetOStream(Warnings0) << "MueLu::Ifpack2Smoother::SetupBlockRelaxation(): reuse of this type is not available (failed cast to CanChangeMatrix), "
                                     "reverting to full construction"
                                  << std::endl;
    }
  }

  if (!reusePreconditioner) {
    ParameterList& myparamList = const_cast<ParameterList&>(this->GetParameterList());
    myparamList.print();
    if (myparamList.isParameter("partitioner: type") &&
        myparamList.get<std::string>("partitioner: type") == "line") {
      Teuchos::RCP<Xpetra::MultiVector<typename Teuchos::ScalarTraits<Scalar>::magnitudeType, LO, GO, NO>> xCoordinates =
          Factory::Get<Teuchos::RCP<Xpetra::MultiVector<typename Teuchos::ScalarTraits<Scalar>::magnitudeType, LO, GO, NO>>>(currentLevel, "Coordinates");
      Teuchos::RCP<Tpetra::MultiVector<typename Teuchos::ScalarTraits<Scalar>::magnitudeType, LO, GO, NO>> coordinates = Teuchos::rcpFromRef(Xpetra::toTpetra<typename Teuchos::ScalarTraits<Scalar>::magnitudeType, LO, GO, NO>(*xCoordinates));

      RCP<Matrix> matA = rcp_dynamic_cast<Matrix>(A_);
      size_t lclSize   = A_->getRangeMap()->getLocalNumElements();
      if (!matA.is_null())
        lclSize = matA->getLocalNumRows();
      size_t numDofsPerNode = lclSize / xCoordinates->getMap()->getLocalNumElements();
      myparamList.set("partitioner: coordinates", coordinates);
      myparamList.set("partitioner: PDE equations", (int)numDofsPerNode);
    }

    prec_ = Ifpack2::Factory::create(type_, tA, overlap_);
    SetPrecParameters();
    prec_->initialize();
  }

  prec_->compute();
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void Ifpack2Smoother<Scalar, LocalOrdinal, GlobalOrdinal, Node>::SetupChebyshev(Level& currentLevel) {
  typedef Tpetra::RowMatrix<SC, LO, GO, NO> tRowMatrix;
  using STS                = Teuchos::ScalarTraits<SC>;
  RCP<BlockedCrsMatrix> bA = rcp_dynamic_cast<BlockedCrsMatrix>(A_);
  if (!bA.is_null())
    A_ = bA->Merge();

  RCP<const tRowMatrix> tA = Utilities::Op2NonConstTpetraRow(A_);

  bool reusePreconditioner = false;

  if (this->IsSetup() == true) {
    // Reuse the constructed preconditioner
    this->GetOStream(Runtime1) << "MueLu::Ifpack2Smoother::SetupChebyshev(): Setup() has already been called, assuming reuse" << std::endl;

    RCP<Ifpack2::Details::CanChangeMatrix<tRowMatrix>> prec = rcp_dynamic_cast<Ifpack2::Details::CanChangeMatrix<tRowMatrix>>(prec_);
    if (!prec.is_null()) {
      prec->setMatrix(tA);
      reusePreconditioner = true;
    } else {
      this->GetOStream(Warnings0) << "MueLu::Ifpack2Smoother::SetupChebyshev(): reuse of this type is not available (failed cast to CanChangeMatrix), "
                                     "reverting to full construction"
                                  << std::endl;
    }
  }

  // Take care of the eigenvalues
  ParameterList& paramList = const_cast<ParameterList&>(this->GetParameterList());
  SC negone                = -STS::one();
  SC lambdaMax             = SetupChebyshevEigenvalues(currentLevel, "A", "", paramList);

  if (!reusePreconditioner) {
    prec_ = Ifpack2::Factory::create(type_, tA, overlap_);
    SetPrecParameters();
    {
      SubFactoryMonitor(*this, "Preconditioner init", currentLevel);
      prec_->initialize();
    }
  } else
    SetPrecParameters();

  {
    SubFactoryMonitor(*this, "Preconditioner compute", currentLevel);
    prec_->compute();
  }

  if (lambdaMax == negone) {
    typedef Tpetra::RowMatrix<SC, LO, GO, NO> MatrixType;

    Teuchos::RCP<Ifpack2::Chebyshev<MatrixType>> chebyPrec = rcp_dynamic_cast<Ifpack2::Chebyshev<MatrixType>>(prec_);
    if (chebyPrec != Teuchos::null) {
      lambdaMax        = chebyPrec->getLambdaMaxForApply();
      RCP<Matrix> matA = rcp_dynamic_cast<Matrix>(A_);
      if (!matA.is_null())
        matA->SetMaxEigenvalueEstimate(lambdaMax);
      this->GetOStream(Statistics1) << "chebyshev: max eigenvalue (calculated by Ifpack2)"
                                    << " = " << lambdaMax << std::endl;
    }
    TEUCHOS_TEST_FOR_EXCEPTION(lambdaMax == negone, Exceptions::RuntimeError, "MueLu::Ifpack2Smoother::Setup(): no maximum eigenvalue estimate");
  }
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void Ifpack2Smoother<Scalar, LocalOrdinal, GlobalOrdinal, Node>::SetupHiptmair(Level& currentLevel) {
  typedef Tpetra::RowMatrix<SC, LO, GO, NO> tRowMatrix;
  RCP<BlockedCrsMatrix> bA = rcp_dynamic_cast<BlockedCrsMatrix>(A_);
  if (!bA.is_null())
    A_ = bA->Merge();

  RCP<const tRowMatrix> tA = Utilities::Op2NonConstTpetraRow(A_);

  bool reusePreconditioner = false;
  if (this->IsSetup() == true) {
    // Reuse the constructed preconditioner
    this->GetOStream(Runtime1) << "MueLu::Ifpack2Smoother::SetupHiptmair(): Setup() has already been called, assuming reuse" << std::endl;

    RCP<Ifpack2::Details::CanChangeMatrix<tRowMatrix>> prec = rcp_dynamic_cast<Ifpack2::Details::CanChangeMatrix<tRowMatrix>>(prec_);
    if (!prec.is_null()) {
      prec->setMatrix(tA);
      reusePreconditioner = true;
    } else {
      this->GetOStream(Warnings0) << "MueLu::Ifpack2Smoother::SetupHiptmair(): reuse of this type is not available (failed cast to CanChangeMatrix), "
                                     "reverting to full construction"
                                  << std::endl;
    }
  }

  // If we're doing Chebyshev subsmoothing, we'll need to get the eigenvalues
  SC negone                = -Teuchos::ScalarTraits<Scalar>::one();
  ParameterList& paramList = const_cast<ParameterList&>(this->GetParameterList());
  std::string smoother1    = paramList.get("hiptmair: smoother type 1", "CHEBYSHEV");
  std::string smoother2    = paramList.get("hiptmair: smoother type 2", "CHEBYSHEV");
  SC lambdaMax11           = negone;

  if (smoother1 == "CHEBYSHEV") {
    ParameterList& list1 = paramList.sublist("hiptmair: smoother list 1");
    lambdaMax11          = SetupChebyshevEigenvalues(currentLevel, "A", "EdgeMatrix ", list1);
  }
  if (smoother2 == "CHEBYSHEV") {
    ParameterList& list2 = paramList.sublist("hiptmair: smoother list 2");
    SetupChebyshevEigenvalues(currentLevel, "A", "EdgeMatrix ", list2);
  }

  // FIXME: Should really add some checks to make sure the eigenvalue calcs worked like in
  // the regular SetupChebyshev

  // Grab the auxillary matrices and stick them on the list
  RCP<Operator> NodeMatrix = currentLevel.Get<RCP<Operator>>("NodeMatrix");
  RCP<Operator> D0         = currentLevel.Get<RCP<Operator>>("D0");

  RCP<tRowMatrix> tNodeMatrix = Utilities::Op2NonConstTpetraRow(NodeMatrix);
  RCP<tRowMatrix> tD0         = Utilities::Op2NonConstTpetraRow(D0);

  Teuchos::ParameterList newlist;
  newlist.set("P", tD0);
  newlist.set("PtAP", tNodeMatrix);
  if (!reusePreconditioner) {
    prec_ = Ifpack2::Factory::create(type_, tA, overlap_);
    SetPrecParameters(newlist);
    prec_->initialize();
  }

  prec_->compute();

  // Post-processing the (1,1) eigenvalue, if we have one
  if (smoother1 == "CHEBYSHEV" && lambdaMax11 == negone) {
    using Teuchos::rcp_dynamic_cast;
    typedef Tpetra::RowMatrix<SC, LO, GO, NO> MatrixType;
    auto hiptmairPrec = rcp_dynamic_cast<Ifpack2::Hiptmair<MatrixType>>(prec_);
    if (hiptmairPrec != Teuchos::null) {
      auto chebyPrec = rcp_dynamic_cast<Ifpack2::Chebyshev<MatrixType>>(hiptmairPrec->getPrec1());
      if (chebyPrec != Teuchos::null) {
        lambdaMax11      = chebyPrec->getLambdaMaxForApply();
        RCP<Matrix> matA = rcp_dynamic_cast<Matrix>(A_);
        if (!matA.is_null())
          matA->SetMaxEigenvalueEstimate(lambdaMax11);
        this->GetOStream(Statistics1) << "chebyshev: max eigenvalue (calculated by Ifpack2)"
                                      << " = " << lambdaMax11 << std::endl;
      }
      TEUCHOS_TEST_FOR_EXCEPTION(lambdaMax11 == negone, Exceptions::RuntimeError, "MueLu::Ifpack2Smoother::Setup(): no maximum eigenvalue estimate");
    }
  }
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
Scalar Ifpack2Smoother<Scalar, LocalOrdinal, GlobalOrdinal, Node>::SetupChebyshevEigenvalues(Level& currentLevel, const std::string& matrixName, const std::string& label, ParameterList& paramList) const {
  // Helper: This gets used for smoothers that want to set up Chebyhev
  typedef Teuchos::ScalarTraits<SC> STS;
  SC negone                  = -STS::one();
  RCP<const Matrix> currentA = currentLevel.Get<RCP<Matrix>>(matrixName);
  SC lambdaMax               = negone;

  std::string maxEigString   = "chebyshev: max eigenvalue";
  std::string eigRatioString = "chebyshev: ratio eigenvalue";

  // Get/calculate the maximum eigenvalue
  if (paramList.isParameter(maxEigString)) {
    if (paramList.isType<double>(maxEigString))
      lambdaMax = paramList.get<double>(maxEigString);
    else
      lambdaMax = paramList.get<SC>(maxEigString);
    this->GetOStream(Statistics1) << label << maxEigString << " (cached with smoother parameter list) = " << lambdaMax << std::endl;
    RCP<Matrix> matA = rcp_dynamic_cast<Matrix>(A_);
    if (!matA.is_null())
      matA->SetMaxEigenvalueEstimate(lambdaMax);

  } else {
    lambdaMax = currentA->GetMaxEigenvalueEstimate();
    if (lambdaMax != negone) {
      this->GetOStream(Statistics1) << label << maxEigString << " (cached with matrix) = " << lambdaMax << std::endl;
      paramList.set(maxEigString, lambdaMax);
    }
  }

  // Calculate the eigenvalue ratio
  const SC defaultEigRatio = 20;

  SC ratio = defaultEigRatio;
  if (paramList.isParameter(eigRatioString)) {
    if (paramList.isType<double>(eigRatioString))
      ratio = paramList.get<double>(eigRatioString);
    else
      ratio = paramList.get<SC>(eigRatioString);
  }
  if (currentLevel.GetLevelID()) {
    // Update ratio to be
    //   ratio = max(number of fine DOFs / number of coarse DOFs, defaultValue)
    //
    // NOTE: We don't need to request previous level matrix as we know for sure it was constructed
    RCP<const Matrix> fineA = currentLevel.GetPreviousLevel()->Get<RCP<Matrix>>(matrixName);
    size_t nRowsFine        = fineA->getGlobalNumRows();
    size_t nRowsCoarse      = currentA->getGlobalNumRows();

    SC levelRatio = as<SC>(as<float>(nRowsFine) / nRowsCoarse);
    if (STS::magnitude(levelRatio) > STS::magnitude(ratio))
      ratio = levelRatio;
  }

  this->GetOStream(Statistics1) << label << eigRatioString << " (computed) = " << ratio << std::endl;
  paramList.set(eigRatioString, ratio);

  if (paramList.isParameter("chebyshev: use rowsumabs diagonal scaling")) {
    this->GetOStream(Runtime1) << "chebyshev: using rowsumabs diagonal scaling" << std::endl;
    bool doScale = false;
    doScale      = paramList.get<bool>("chebyshev: use rowsumabs diagonal scaling");
    paramList.remove("chebyshev: use rowsumabs diagonal scaling");
    double chebyReplaceTol = Teuchos::ScalarTraits<Scalar>::eps() * 100;
    std::string paramName  = "chebyshev: rowsumabs diagonal replacement tolerance";
    if (paramList.isParameter(paramName)) {
      chebyReplaceTol = paramList.get<double>(paramName);
      paramList.remove(paramName);
    }
    double chebyReplaceVal = Teuchos::ScalarTraits<double>::zero();
    paramName              = "chebyshev: rowsumabs diagonal replacement value";
    if (paramList.isParameter(paramName)) {
      chebyReplaceVal = paramList.get<double>(paramName);
      paramList.remove(paramName);
    }
    bool chebyReplaceSingleEntryRowWithZero = false;
    paramName                               = "chebyshev: rowsumabs replace single entry row with zero";
    if (paramList.isParameter(paramName)) {
      chebyReplaceSingleEntryRowWithZero = paramList.get<bool>(paramName);
      paramList.remove(paramName);
    }
    bool useAverageAbsDiagVal = false;
    paramName                 = "chebyshev: rowsumabs use automatic diagonal tolerance";
    if (paramList.isParameter(paramName)) {
      useAverageAbsDiagVal = paramList.get<bool>(paramName);
      paramList.remove(paramName);
    }
    if (doScale) {
      const bool doReciprocal                                                       = true;
      RCP<Vector> lumpedDiagonal                                                    = Utilities::GetLumpedMatrixDiagonal(*currentA, doReciprocal, chebyReplaceTol, chebyReplaceVal, chebyReplaceSingleEntryRowWithZero, useAverageAbsDiagVal);
      const Xpetra::TpetraVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& tmpVec = dynamic_cast<const Xpetra::TpetraVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>&>(*lumpedDiagonal);
      paramList.set("chebyshev: operator inv diagonal", tmpVec.getTpetra_Vector());
    }
  }

  return lambdaMax;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void Ifpack2Smoother<Scalar, LocalOrdinal, GlobalOrdinal, Node>::SetupGeneric(Level& /* currentLevel */) {
  typedef Tpetra::RowMatrix<SC, LO, GO, NO> tRowMatrix;
  RCP<BlockedCrsMatrix> bA = rcp_dynamic_cast<BlockedCrsMatrix>(A_);
  if (!bA.is_null())
    A_ = bA->Merge();

  RCP<const tRowMatrix> tA = Utilities::Op2NonConstTpetraRow(A_);

  bool reusePreconditioner = false;
  if (this->IsSetup() == true) {
    // Reuse the constructed preconditioner
    this->GetOStream(Runtime1) << "MueLu::Ifpack2Smoother::SetupGeneric(): Setup() has already been called, assuming reuse" << std::endl;

    RCP<Ifpack2::Details::CanChangeMatrix<tRowMatrix>> prec = rcp_dynamic_cast<Ifpack2::Details::CanChangeMatrix<tRowMatrix>>(prec_);
    if (!prec.is_null()) {
      prec->setMatrix(tA);
      reusePreconditioner = true;
    } else {
      this->GetOStream(Warnings0) << "MueLu::Ifpack2Smoother::SetupGeneric(): reuse of this type is not available (failed cast to CanChangeMatrix), "
                                     "reverting to full construction"
                                  << std::endl;
    }
  }

  if (!reusePreconditioner) {
    prec_ = Ifpack2::Factory::create(type_, tA, overlap_);
    SetPrecParameters();
    prec_->initialize();
  }

  prec_->compute();
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void Ifpack2Smoother<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Apply(MultiVector& X, const MultiVector& B, bool InitialGuessIsZero) const {
  TEUCHOS_TEST_FOR_EXCEPTION(SmootherPrototype::IsSetup() == false, Exceptions::RuntimeError, "MueLu::Ifpack2Smoother::Apply(): Setup() has not been called");

  // Forward the InitialGuessIsZero option to Ifpack2
  // TODO:  It might be nice to switch back the internal
  //        "zero starting solution" option of the ifpack2 object prec_ to his
  //        initial value at the end but there is no way right now to get
  //        the current value of the "zero starting solution" in ifpack2.
  //        It's not really an issue, as prec_  can only be used by this method.
  Teuchos::ParameterList paramList;
  bool supportInitialGuess            = false;
  const Teuchos::ParameterList params = this->GetParameterList();

  if (prec_->supportsZeroStartingSolution()) {
    prec_->setZeroStartingSolution(InitialGuessIsZero);
    supportInitialGuess = true;
  } else if (type_ == "SCHWARZ") {
    paramList.set("schwarz: zero starting solution", InitialGuessIsZero);
    // Because additive Schwarz has "delta" semantics, it's sufficient to
    // toggle only the zero initial guess flag, and not pass in already
    // set parameters.  If we call SetPrecParameters, the subdomain solver
    // will be destroyed.
    prec_->setParameters(paramList);
    supportInitialGuess = true;
  }

  // TODO JJH 30Apr2014  Calling SetPrecParameters(paramList) when the smoother
  // is Ifpack2::AdditiveSchwarz::setParameterList() will destroy the subdomain
  //(aka inner) solver.  This behavior is documented but a departure from what
  // it previously did, and what other Ifpack2 solvers currently do.  So I have
  // moved SetPrecParameters(paramList) into the if-else block above.

  // Apply
  if (InitialGuessIsZero || supportInitialGuess) {
    Tpetra::MultiVector<SC, LO, GO, NO>& tpX       = Utilities::MV2NonConstTpetraMV(X);
    const Tpetra::MultiVector<SC, LO, GO, NO>& tpB = Utilities::MV2TpetraMV(B);
    prec_->apply(tpB, tpX);
  } else {
    typedef Teuchos::ScalarTraits<Scalar> TST;

    RCP<MultiVector> Residual;
    {
      std::string prefix  = this->ShortClassName() + ": Apply: ";
      RCP<TimeMonitor> tM = rcp(new TimeMonitor(*this, prefix + "residual calculation", Timings0));
      Residual            = Utilities::Residual(*A_, X, B);
    }

    RCP<MultiVector> Correction = MultiVectorFactory::Build(A_->getDomainMap(), X.getNumVectors());

    Tpetra::MultiVector<SC, LO, GO, NO>& tpX       = Utilities::MV2NonConstTpetraMV(*Correction);
    const Tpetra::MultiVector<SC, LO, GO, NO>& tpB = Utilities::MV2TpetraMV(*Residual);

    prec_->apply(tpB, tpX);

    X.update(TST::one(), *Correction, TST::one());
  }
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<MueLu::SmootherPrototype<Scalar, LocalOrdinal, GlobalOrdinal, Node>> Ifpack2Smoother<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Copy() const {
  RCP<Ifpack2Smoother> smoother = rcp(new Ifpack2Smoother(*this));
  smoother->SetParameterList(this->GetParameterList());
  return smoother;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
std::string Ifpack2Smoother<Scalar, LocalOrdinal, GlobalOrdinal, Node>::description() const {
  std::ostringstream out;
  if (SmootherPrototype::IsSetup()) {
    out << prec_->description();
  } else {
    out << SmootherPrototype::description();
    out << "{type = " << type_ << "}";
  }
  return out.str();
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void Ifpack2Smoother<Scalar, LocalOrdinal, GlobalOrdinal, Node>::print(Teuchos::FancyOStream& out, const VerbLevel verbLevel) const {
  MUELU_DESCRIBE;

  if (verbLevel & Parameters0)
    out0 << "Prec. type: " << type_ << std::endl;

  if (verbLevel & Parameters1) {
    out0 << "Parameter list: " << std::endl;
    Teuchos::OSTab tab2(out);
    out << this->GetParameterList();
    out0 << "Overlap: " << overlap_ << std::endl;
  }

  if (verbLevel & External)
    if (prec_ != Teuchos::null) {
      Teuchos::OSTab tab2(out);
      out << *prec_ << std::endl;
    }

  if (verbLevel & Debug) {
    out0 << "IsSetup: " << Teuchos::toString(SmootherPrototype::IsSetup()) << std::endl
         << "-" << std::endl
         << "RCP<prec_>: " << prec_ << std::endl;
  }
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
size_t Ifpack2Smoother<Scalar, LocalOrdinal, GlobalOrdinal, Node>::getNodeSmootherComplexity() const {
  typedef Tpetra::RowMatrix<SC, LO, GO, NO> MatrixType;
  // NOTE: Only works for a subset of Ifpack2's smoothers
  RCP<Ifpack2::Relaxation<MatrixType>> pr = rcp_dynamic_cast<Ifpack2::Relaxation<MatrixType>>(prec_);
  if (!pr.is_null()) return pr->getNodeSmootherComplexity();

  RCP<Ifpack2::Chebyshev<MatrixType>> pc = rcp_dynamic_cast<Ifpack2::Chebyshev<MatrixType>>(prec_);
  if (!pc.is_null()) return pc->getNodeSmootherComplexity();

  RCP<Ifpack2::BlockRelaxation<MatrixType>> pb = rcp_dynamic_cast<Ifpack2::BlockRelaxation<MatrixType>>(prec_);
  if (!pb.is_null()) return pb->getNodeSmootherComplexity();

  RCP<Ifpack2::ILUT<MatrixType>> pi = rcp_dynamic_cast<Ifpack2::ILUT<MatrixType>>(prec_);
  if (!pi.is_null()) return pi->getNodeSmootherComplexity();

  RCP<Ifpack2::RILUK<MatrixType>> pk = rcp_dynamic_cast<Ifpack2::RILUK<MatrixType>>(prec_);
  if (!pk.is_null()) return pk->getNodeSmootherComplexity();

  return Teuchos::OrdinalTraits<size_t>::invalid();
}

}  // namespace MueLu

#endif  // HAVE_MUELU_IFPACK2
#endif  // MUELU_IFPACK2SMOOTHER_DEF_HPP
