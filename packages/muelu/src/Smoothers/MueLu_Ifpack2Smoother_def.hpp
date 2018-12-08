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

#if defined(HAVE_MUELU_TPETRA) && defined(HAVE_MUELU_IFPACK2)

#include <Teuchos_ParameterList.hpp>

#include <Tpetra_RowMatrix.hpp>

#include <Ifpack2_Chebyshev.hpp>
#include <Ifpack2_Relaxation.hpp>
#include <Ifpack2_ILUT.hpp>
#include <Ifpack2_BlockRelaxation.hpp>
#include <Ifpack2_Factory.hpp>
#include <Ifpack2_Parameters.hpp>

#include <Xpetra_BlockedCrsMatrix.hpp>
#include <Xpetra_CrsMatrix.hpp>
#include <Xpetra_CrsMatrixWrap.hpp>
#include <Xpetra_Matrix.hpp>
#include <Xpetra_MultiVectorFactory.hpp>
#include <Xpetra_TpetraMultiVector.hpp>

#include "MueLu_Ifpack2Smoother_decl.hpp"
#include "MueLu_Level.hpp"
#include "MueLu_FactoryManagerBase.hpp"
#include "MueLu_Utilities.hpp"
#include "MueLu_Monitor.hpp"

#ifdef HAVE_MUELU_INTREPID2
#include "MueLu_IntrepidPCoarsenFactory_decl.hpp"
#include "MueLu_IntrepidPCoarsenFactory_def.hpp"
#include "Intrepid2_Basis.hpp"
#include "Kokkos_DynRankView.hpp"
#endif

// #define IFPACK2_HAS_PROPER_REUSE

namespace MueLu {

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node>
  Ifpack2Smoother<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Ifpack2Smoother(const std::string& type, const Teuchos::ParameterList& paramList, const LO& overlap)
    : type_(type), overlap_(overlap)
  {
    SetParameterList(paramList);
  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node>
  void Ifpack2Smoother<Scalar, LocalOrdinal, GlobalOrdinal, Node>::SetParameterList(const Teuchos::ParameterList& paramList) {
    Factory::SetParameterList(paramList);

    if (SmootherPrototype::IsSetup()) {
      // It might be invalid to change parameters after the setup, but it depends entirely on Ifpack implementation.
      // TODO: I don't know if Ifpack returns an error code or exception or ignore parameters modification in this case...
      SetPrecParameters();
    }
  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node>
  void Ifpack2Smoother<Scalar, LocalOrdinal, GlobalOrdinal, Node>::SetPrecParameters(const Teuchos::ParameterList& list) const {
    ParameterList& paramList = const_cast<ParameterList&>(this->GetParameterList());
    paramList.setParameters(list);

    RCP<ParameterList> precList = this->RemoveFactoriesFromList(this->GetParameterList());

    prec_->setParameters(*precList);

    paramList.setParameters(*precList); // what about that??
  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node>
  void Ifpack2Smoother<Scalar, LocalOrdinal, GlobalOrdinal, Node>::DeclareInput(Level& currentLevel) const {
    this->Input(currentLevel, "A");

    if (type_ == "LINESMOOTHING_TRIDI_RELAXATION"        ||
        type_ == "LINESMOOTHING_TRIDI RELAXATION"        ||
        type_ == "LINESMOOTHING_TRIDIRELAXATION"         ||
        type_ == "LINESMOOTHING_TRIDIAGONAL_RELAXATION"  ||
        type_ == "LINESMOOTHING_TRIDIAGONAL RELAXATION"  ||
        type_ == "LINESMOOTHING_TRIDIAGONALRELAXATION"   ||
        type_ == "LINESMOOTHING_BANDED_RELAXATION"       ||
        type_ == "LINESMOOTHING_BANDED RELAXATION"       ||
        type_ == "LINESMOOTHING_BANDEDRELAXATION"        ||
        type_ == "LINESMOOTHING_BLOCK_RELAXATION"        ||
        type_ == "LINESMOOTHING_BLOCK RELAXATION"        ||
        type_ == "LINESMOOTHING_BLOCKRELAXATION") {
      this->Input(currentLevel, "CoarseNumZLayers");            // necessary for fallback criterion
      this->Input(currentLevel, "LineDetection_VertLineIds");   // necessary to feed block smoother
    }
    else if (type_ == "BLOCK RELAXATION" ||
             type_ == "BLOCK_RELAXATION" ||
             type_ == "BLOCKRELAXATION")
    {
      //We need to check for the "partitioner type" = "line"
      ParameterList precList = this->GetParameterList();
      if(precList.isParameter("partitioner: type") &&
         precList.get<std::string>("partitioner: type") == "line") {
        this->Input(currentLevel, "Coordinates");
      }
    }
    else if (type_ == "TOPOLOGICAL")
    {
      // for the topological smoother, we require an element to node map:
      this->Input(currentLevel, "pcoarsen: element to node map");
    }
  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node>
  void Ifpack2Smoother<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Setup(Level& currentLevel) {
    FactoryMonitor m(*this, "Setup Smoother", currentLevel);

    A_ = Factory::Get< RCP<Matrix> >(currentLevel, "A");

    if      (type_ == "SCHWARZ")
      SetupSchwarz(currentLevel);

    else if (type_ == "LINESMOOTHING_TRIDI_RELAXATION"       ||
             type_ == "LINESMOOTHING_TRIDI RELAXATION"       ||
             type_ == "LINESMOOTHING_TRIDIRELAXATION"        ||
             type_ == "LINESMOOTHING_TRIDIAGONAL_RELAXATION" ||
             type_ == "LINESMOOTHING_TRIDIAGONAL RELAXATION" ||
             type_ == "LINESMOOTHING_TRIDIAGONALRELAXATION"  ||
             type_ == "LINESMOOTHING_BANDED_RELAXATION"      ||
             type_ == "LINESMOOTHING_BANDED RELAXATION"      ||
             type_ == "LINESMOOTHING_BANDEDRELAXATION"       ||
             type_ == "LINESMOOTHING_BLOCK_RELAXATION"       ||
             type_ == "LINESMOOTHING_BLOCK RELAXATION"       ||
             type_ == "LINESMOOTHING_BLOCKRELAXATION")
      SetupLineSmoothing(currentLevel);

    else if (type_ == "BLOCK_RELAXATION" ||
             type_ == "BLOCK RELAXATION" ||
             type_ == "BLOCKRELAXATION")
      SetupBlockRelaxation(currentLevel);

    else if (type_ == "CHEBYSHEV")
      SetupChebyshev(currentLevel);

    else if (type_ == "TOPOLOGICAL")
    {
#ifdef HAVE_MUELU_INTREPID2
      SetupTopological(currentLevel);
#else
      TEUCHOS_TEST_FOR_EXCEPTION(true, std::invalid_argument, "'TOPOLOGICAL' smoother choice requires Intrepid2");
#endif
    }
    else
    {
      SetupGeneric(currentLevel);
    }

    SmootherPrototype::IsSetup(true);

    this->GetOStream(Statistics1) << description() << std::endl;
  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node>
  void Ifpack2Smoother<Scalar, LocalOrdinal, GlobalOrdinal, Node>::SetupSchwarz(Level& currentLevel) {
    typedef Tpetra::RowMatrix<SC,LO,GO,NO> tRowMatrix;

    bool reusePreconditioner = false;
    if (this->IsSetup() == true) {
      // Reuse the constructed preconditioner
      this->GetOStream(Runtime1) << "MueLu::Ifpack2Smoother::SetupSchwarz(): Setup() has already been called, assuming reuse" << std::endl;

      bool isTRowMatrix = true;
      RCP<const tRowMatrix> tA;
      try {
        tA = Utilities::Op2NonConstTpetraRow(A_);
      } catch (Exceptions::BadCast) {
        isTRowMatrix = false;
      }

      RCP<Ifpack2::Details::CanChangeMatrix<tRowMatrix> > prec = rcp_dynamic_cast<Ifpack2::Details::CanChangeMatrix<tRowMatrix> >(prec_);
      if (!prec.is_null() && isTRowMatrix) {
#ifdef IFPACK2_HAS_PROPER_REUSE
        prec->resetMatrix(tA);
        reusePreconditioner = true;
#else
        this->GetOStream(Errors) << "Ifpack2 does not have proper reuse yet." << std::endl;
#endif

      } else {
        this->GetOStream(Warnings0) << "MueLu::Ifpack2Smoother::SetupSchwarz(): reuse of this type is not available "
            "(either failed cast to CanChangeMatrix, or to Tpetra Row Matrix), reverting to full construction" << std::endl;
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
        if (subList.isParameter(partName) && subList.get<std::string>(partName) == "user") {
          isBlockedMatrix = true;

          RCP<BlockedCrsMatrix> bA = rcp_dynamic_cast<BlockedCrsMatrix>(A_);
          TEUCHOS_TEST_FOR_EXCEPTION(bA.is_null(), Exceptions::BadCast,
                                     "Matrix A must be of type BlockedCrsMatrix.");

          size_t numVels = bA->getMatrix(0,0)->getNodeNumRows();
          size_t numPres = bA->getMatrix(1,0)->getNodeNumRows();
          size_t numRows = A_->getNodeNumRows();

          ArrayRCP<LocalOrdinal> blockSeeds(numRows, Teuchos::OrdinalTraits<LocalOrdinal>::invalid());

          size_t numBlocks = 0;
          for (size_t rowOfB = numVels; rowOfB < numVels+numPres; ++rowOfB)
            blockSeeds[rowOfB] = numBlocks++;

          RCP<BlockedCrsMatrix> bA2 = rcp_dynamic_cast<BlockedCrsMatrix>(A_);
          TEUCHOS_TEST_FOR_EXCEPTION(bA2.is_null(), Exceptions::BadCast,
                                     "Matrix A must be of type BlockedCrsMatrix.");

          merged2Mat = bA2->Merge();

          // Add Dirichlet rows to the list of seeds
          ArrayRCP<const bool> boundaryNodes = Utilities::DetectDirichletRows(*merged2Mat, 0.0);
          bool haveBoundary = false;
          for (LO i = 0; i < boundaryNodes.size(); i++)
            if (boundaryNodes[i]) {
              // FIXME:
              // 1. would not this [] overlap with some in the previos blockSeed loop?
              // 2. do we need to distinguish between pressure and velocity Dirichlet b.c.
              blockSeeds[i] = numBlocks;
              haveBoundary = true;
            }
          if (haveBoundary)
            numBlocks++;

          subList.set("partitioner: map",         blockSeeds);
          subList.set("partitioner: local parts", as<int>(numBlocks));

        } else {
          RCP<BlockedCrsMatrix> bA = rcp_dynamic_cast<BlockedCrsMatrix>(A_);
          if (!bA.is_null()) {
            isBlockedMatrix = true;
            merged2Mat = bA->Merge();
          }
        }
      }

      RCP<const tRowMatrix> tA;
      if (isBlockedMatrix == true) tA = Utilities::Op2NonConstTpetraRow(merged2Mat);
      else                         tA = Utilities::Op2NonConstTpetraRow(A_);

      prec_ = Ifpack2::Factory::create(type_, tA, overlap_);
      SetPrecParameters();

      prec_->initialize();
    }

    prec_->compute();
  }

#ifdef HAVE_MUELU_INTREPID2
  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node>
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
    
    typedef Kokkos::DynRankView<LocalOrdinal,typename Node::device_type> FCO; //
    
    LocalOrdinal  lo_invalid = Teuchos::OrdinalTraits<LO>::invalid();
    
    using namespace std;
    
    const Teuchos::RCP<FCO> elemToNode = Factory::Get<Teuchos::RCP<FCO> >(currentLevel,"pcoarsen: element to node map");
    
    string basisString = paramList.get<string>("pcoarsen: hi basis");
    int degree;
    // NOTE: To make sure Stokhos works we only instantiate these guys with double.  There's a lot
    // of stuff in the guts of Intrepid2 that doesn't play well with Stokhos as of yet.  Here, we only
    // care about the assignment of basis ordinals to topological entities, so this code is actually
    // independent of the Scalar type--hard-coding double here won't hurt us.
    auto basis = MueLuIntrepid::BasisFactory<double,ES>(basisString, degree);
    
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
      TEUCHOS_TEST_FOR_EXCEPTION(true,std::invalid_argument,"Unrecognized smoother neighborhood type.  Supported types are node, edge, face.");
    vector<vector<LocalOrdinal>> seeds;
    MueLuIntrepid::FindGeometricSeedOrdinals(basis, *elemToNode, seeds, *A_->getRowMap(), *A_->getColMap());
    
    // Ifpack2 wants the seeds in an array of the same length as the number of local elements,
    // with local partition #s marked for the ones that are seeds, and invalid for the rest
    int myNodeCount = A_->getRowMap()->getNodeNumElements();
    ArrayRCP<LocalOrdinal> nodeSeeds(myNodeCount,lo_invalid);
    int localPartitionNumber = 0;
    for (LocalOrdinal seed : seeds[dimension])
    {
      nodeSeeds[seed] = localPartitionNumber++;
    }
    
    paramList.remove("smoother: neighborhood type");
    paramList.remove("pcoarsen: hi basis");
    
    paramList.set("partitioner: map", nodeSeeds);
    paramList.set("partitioner: type", "user");
    paramList.set("partitioner: overlap", 1);
    paramList.set("partitioner: local parts", int(seeds[dimension].size()));

    RCP<const Tpetra::RowMatrix<SC, LO, GO, NO> > tA = Utilities::Op2NonConstTpetraRow(A_);
    
    type_ = "BLOCKRELAXATION";
    prec_ = Ifpack2::Factory::create(type_, tA, overlap_);
    SetPrecParameters();
    prec_->initialize();
    prec_->compute();
  }
#endif

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node>
  void Ifpack2Smoother<Scalar, LocalOrdinal, GlobalOrdinal, Node>::SetupLineSmoothing(Level& currentLevel) {
    if (this->IsSetup() == true) {
      this->GetOStream(Warnings0) << "MueLu::Ifpack2Smoother::SetupLineSmoothing(): Setup() has already been called" << std::endl;
      this->GetOStream(Warnings0) << "MueLu::Ifpack2Smoother::SetupLineSmoothing(): reuse of this type is not available, reverting to full construction" << std::endl;
    }

    ParameterList& myparamList = const_cast<ParameterList&>(this->GetParameterList());

    LO CoarseNumZLayers = Factory::Get<LO>(currentLevel,"CoarseNumZLayers");
    if (CoarseNumZLayers > 0) {
      Teuchos::ArrayRCP<LO> TVertLineIdSmoo = Factory::Get< Teuchos::ArrayRCP<LO> >(currentLevel, "LineDetection_VertLineIds");

      // determine number of local parts
      LO maxPart = 0;
      for(size_t k = 0; k < Teuchos::as<size_t>(TVertLineIdSmoo.size()); k++) {
        if(maxPart < TVertLineIdSmoo[k]) maxPart = TVertLineIdSmoo[k];
      }

      size_t numLocalRows = A_->getNodeNumRows();
      TEUCHOS_TEST_FOR_EXCEPTION(numLocalRows % TVertLineIdSmoo.size() != 0, Exceptions::RuntimeError,
        "MueLu::Ifpack2Smoother::Setup(): the number of local nodes is incompatible with the TVertLineIdsSmoo.");

      if (numLocalRows == Teuchos::as<size_t>(TVertLineIdSmoo.size())) {
        myparamList.set("partitioner: type","user");
        myparamList.set("partitioner: map",TVertLineIdSmoo);
        myparamList.set("partitioner: local parts",maxPart+1);
      } else {
        // we assume a constant number of DOFs per node
        size_t numDofsPerNode = numLocalRows / TVertLineIdSmoo.size();

        // Create a new Teuchos::ArrayRCP<LO> of size numLocalRows and fill it with the corresponding information
        Teuchos::ArrayRCP<LO> partitionerMap(numLocalRows, Teuchos::OrdinalTraits<LocalOrdinal>::invalid());
        for (size_t blockRow = 0; blockRow < Teuchos::as<size_t>(TVertLineIdSmoo.size()); ++blockRow)
          for (size_t dof = 0; dof < numDofsPerNode; dof++)
            partitionerMap[blockRow * numDofsPerNode + dof] = TVertLineIdSmoo[blockRow];
        myparamList.set("partitioner: type","user");
        myparamList.set("partitioner: map",partitionerMap);
        myparamList.set("partitioner: local parts",maxPart + 1);
      }

      if (type_ == "LINESMOOTHING_BANDED_RELAXATION" ||
          type_ == "LINESMOOTHING_BANDED RELAXATION" ||
          type_ == "LINESMOOTHING_BANDEDRELAXATION")
        type_ = "BANDEDRELAXATION";
      else if (type_ == "LINESMOOTHING_TRIDI_RELAXATION"       ||
               type_ == "LINESMOOTHING_TRIDI RELAXATION"       ||
               type_ == "LINESMOOTHING_TRIDIRELAXATION"        ||
               type_ == "LINESMOOTHING_TRIDIAGONAL_RELAXATION" ||
               type_ == "LINESMOOTHING_TRIDIAGONAL RELAXATION" ||
               type_ == "LINESMOOTHING_TRIDIAGONALRELAXATION")
        type_ = "TRIDIAGONALRELAXATION";
      else
        type_ = "BLOCKRELAXATION";
    } else {
      // line detection failed -> fallback to point-wise relaxation
      this->GetOStream(Runtime0) << "Line detection failed: fall back to point-wise relaxation" << std::endl;
      myparamList.remove("partitioner: type",false);
      myparamList.remove("partitioner: map", false);
      myparamList.remove("partitioner: local parts",false);
      type_ = "RELAXATION";
    }

    RCP<const Tpetra::RowMatrix<SC, LO, GO, NO> > tA = Utilities::Op2NonConstTpetraRow(A_);

    prec_ = Ifpack2::Factory::create(type_, tA, overlap_);
    SetPrecParameters();
    prec_->initialize();
    prec_->compute();
  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node>
  void Ifpack2Smoother<Scalar, LocalOrdinal, GlobalOrdinal, Node>::SetupBlockRelaxation(Level& currentLevel) {
    typedef Tpetra::RowMatrix<SC,LO,GO,NO> tRowMatrix;

    RCP<BlockedCrsMatrix> bA = rcp_dynamic_cast<BlockedCrsMatrix>(A_);
    if (!bA.is_null())
      A_ = bA->Merge();

    RCP<const tRowMatrix> tA = Utilities::Op2NonConstTpetraRow(A_);

    bool reusePreconditioner = false;
    if (this->IsSetup() == true) {
      // Reuse the constructed preconditioner
      this->GetOStream(Runtime1) << "MueLu::Ifpack2Smoother::SetupGeneric(): Setup() has already been called, assuming reuse" << std::endl;

      RCP<Ifpack2::Details::CanChangeMatrix<tRowMatrix> > prec = rcp_dynamic_cast<Ifpack2::Details::CanChangeMatrix<tRowMatrix> >(prec_);
      if (!prec.is_null()) {
#ifdef IFPACK2_HAS_PROPER_REUSE
        prec->resetMatrix(tA);
        reusePreconditioner = true;
#else
        this->GetOStream(Errors) << "Ifpack2 does not have proper reuse yet." << std::endl;
#endif

      } else {
        this->GetOStream(Warnings0) << "MueLu::Ifpack2Smoother::SetupSchwarz(): reuse of this type is not available (failed cast to CanChangeMatrix), "
            "reverting to full construction" << std::endl;
      }
    }

    if (!reusePreconditioner) {
      ParameterList& myparamList = const_cast<ParameterList&>(this->GetParameterList());
      myparamList.print();
      if(myparamList.isParameter("partitioner: type") &&
         myparamList.get<std::string>("partitioner: type") == "line") {
        Teuchos::RCP<Xpetra::MultiVector<typename Teuchos::ScalarTraits<Scalar>::magnitudeType,LO,GO,NO> > xCoordinates =
          Factory::Get<Teuchos::RCP<Xpetra::MultiVector<typename Teuchos::ScalarTraits<Scalar>::magnitudeType,LO,GO,NO> > >(currentLevel, "Coordinates");
        Teuchos::RCP<Tpetra::MultiVector<typename Teuchos::ScalarTraits<Scalar>::magnitudeType,LO,GO,NO> > coordinates = Teuchos::rcpFromRef(Xpetra::toTpetra<typename Teuchos::ScalarTraits<Scalar>::magnitudeType,LO,GO,NO>(*xCoordinates));
        myparamList.set("partitioner: coordinates", coordinates);
      }

      prec_ = Ifpack2::Factory::create(type_, tA, overlap_);
      SetPrecParameters();
      prec_->initialize();
    }

    prec_->compute();
  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node>
  void Ifpack2Smoother<Scalar, LocalOrdinal, GlobalOrdinal, Node>::SetupChebyshev(Level& currentLevel) {
    if (this->IsSetup() == true) {
      this->GetOStream(Warnings0) << "MueLu::Ifpack2Smoother::SetupChebyshev(): SetupChebyshev() has already been called" << std::endl;
      this->GetOStream(Warnings0) << "MueLu::Ifpack2Smoother::SetupChebyshev(): reuse of this type is not available, reverting to full construction" << std::endl;
    }

    typedef Teuchos::ScalarTraits<SC> STS;
    SC negone = -STS::one();

    SC lambdaMax = negone;
    {
      std::string maxEigString   = "chebyshev: max eigenvalue";
      std::string eigRatioString = "chebyshev: ratio eigenvalue";

      ParameterList& paramList = const_cast<ParameterList&>(this->GetParameterList());

      // Get/calculate the maximum eigenvalue
      if (paramList.isParameter(maxEigString)) {
        if (paramList.isType<double>(maxEigString))
          lambdaMax = paramList.get<double>(maxEigString);
        else
          lambdaMax = paramList.get<SC>(maxEigString);
        this->GetOStream(Statistics1) << maxEigString << " (cached with smoother parameter list) = " << lambdaMax << std::endl;

      } else {
        lambdaMax = A_->GetMaxEigenvalueEstimate();
        if (lambdaMax != negone) {
          this->GetOStream(Statistics1) << maxEigString << " (cached with matrix) = " << lambdaMax << std::endl;
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
        RCP<const Matrix> fineA = currentLevel.GetPreviousLevel()->Get<RCP<Matrix> >("A");
        size_t nRowsFine   = fineA->getGlobalNumRows();
        size_t nRowsCoarse = A_->getGlobalNumRows();

        SC levelRatio = as<SC>(as<float>(nRowsFine)/nRowsCoarse);
        if (STS::magnitude(levelRatio) > STS::magnitude(ratio))
          ratio = levelRatio;
      }

      this->GetOStream(Statistics1) << eigRatioString << " (computed) = " << ratio << std::endl;
      paramList.set(eigRatioString, ratio);
    }

    RCP<const Tpetra::RowMatrix<SC, LO, GO, NO> > tA = Utilities::Op2NonConstTpetraRow(A_);

    prec_ = Ifpack2::Factory::create(type_, tA, overlap_);
    SetPrecParameters();
    {
      SubFactoryMonitor(*this, "Preconditioner init", currentLevel);
      prec_->initialize();
    }
    {
      SubFactoryMonitor(*this, "Preconditioner compute", currentLevel);
      prec_->compute();
    }

    if (lambdaMax == negone) {
      typedef Tpetra::RowMatrix<SC, LO, GO, NO> MatrixType;

      Teuchos::RCP<Ifpack2::Chebyshev<MatrixType> > chebyPrec = rcp_dynamic_cast<Ifpack2::Chebyshev<MatrixType> >(prec_);
      if (chebyPrec != Teuchos::null) {
        lambdaMax = chebyPrec->getLambdaMaxForApply();
        A_->SetMaxEigenvalueEstimate(lambdaMax);
        this->GetOStream(Statistics1) << "chebyshev: max eigenvalue (calculated by Ifpack2)" << " = " << lambdaMax << std::endl;
      }
      TEUCHOS_TEST_FOR_EXCEPTION(lambdaMax == negone, Exceptions::RuntimeError, "MueLu::Ifpack2Smoother::Setup(): no maximum eigenvalue estimate");
    }
  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node>
  void Ifpack2Smoother<Scalar, LocalOrdinal, GlobalOrdinal, Node>::SetupGeneric(Level& currentLevel) {
    typedef Tpetra::RowMatrix<SC,LO,GO,NO> tRowMatrix;

    RCP<BlockedCrsMatrix> bA = rcp_dynamic_cast<BlockedCrsMatrix>(A_);
    if (!bA.is_null())
      A_ = bA->Merge();

    RCP<const tRowMatrix> tA = Utilities::Op2NonConstTpetraRow(A_);

    bool reusePreconditioner = false;
    if (this->IsSetup() == true) {
      // Reuse the constructed preconditioner
      this->GetOStream(Runtime1) << "MueLu::Ifpack2Smoother::SetupGeneric(): Setup() has already been called, assuming reuse" << std::endl;

      RCP<Ifpack2::Details::CanChangeMatrix<tRowMatrix> > prec = rcp_dynamic_cast<Ifpack2::Details::CanChangeMatrix<tRowMatrix> >(prec_);
      if (!prec.is_null()) {
#ifdef IFPACK2_HAS_PROPER_REUSE
        prec->resetMatrix(tA);
        reusePreconditioner = true;
#else
        this->GetOStream(Errors) << "Ifpack2 does not have proper reuse yet." << std::endl;
#endif

      } else {
        this->GetOStream(Warnings0) << "MueLu::Ifpack2Smoother::SetupSchwarz(): reuse of this type is not available (failed cast to CanChangeMatrix), "
            "reverting to full construction" << std::endl;
      }
    }

    if (!reusePreconditioner) {
      prec_ = Ifpack2::Factory::create(type_, tA, overlap_);
      SetPrecParameters();
      prec_->initialize();
    }

    prec_->compute();
  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node>
  void Ifpack2Smoother<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Apply(MultiVector& X, const MultiVector& B, bool InitialGuessIsZero) const {
    TEUCHOS_TEST_FOR_EXCEPTION(SmootherPrototype::IsSetup() == false, Exceptions::RuntimeError, "MueLu::Ifpack2Smoother::Apply(): Setup() has not been called");

    // Forward the InitialGuessIsZero option to Ifpack2
    // TODO:  It might be nice to switch back the internal
    //        "zero starting solution" option of the ifpack2 object prec_ to his
    //        initial value at the end but there is no way right now to get
    //        the current value of the "zero starting solution" in ifpack2.
    //        It's not really an issue, as prec_  can only be used by this method.
    // TODO: When https://software.sandia.gov/bugzilla/show_bug.cgi?id=5283#c2 is done
    // we should remove the if/else/elseif and just test if this
    // option is supported by current ifpack2 preconditioner
    Teuchos::ParameterList paramList;
    bool supportInitialGuess = false;
    if (type_ == "CHEBYSHEV") {
      paramList.set("chebyshev: zero starting solution", InitialGuessIsZero);
      SetPrecParameters(paramList);
      supportInitialGuess = true;

    } else if (type_ == "RELAXATION") {
      paramList.set("relaxation: zero starting solution", InitialGuessIsZero);
      SetPrecParameters(paramList);
      supportInitialGuess = true;

    } else if (type_ == "KRYLOV") {
      paramList.set("krylov: zero starting solution", InitialGuessIsZero);
      SetPrecParameters(paramList);
      supportInitialGuess = true;

    } else if (type_ == "SCHWARZ") {
      paramList.set("schwarz: zero starting solution", InitialGuessIsZero);
      //Because additive Schwarz has "delta" semantics, it's sufficient to
      //toggle only the zero initial guess flag, and not pass in already
      //set parameters.  If we call SetPrecParameters, the subdomain solver
      //will be destroyed.
      prec_->setParameters(paramList);
      supportInitialGuess = true;
    }

    //TODO JJH 30Apr2014  Calling SetPrecParameters(paramList) when the smoother
    //is Ifpack2::AdditiveSchwarz::setParameterList() will destroy the subdomain
    //(aka inner) solver.  This behavior is documented but a departure from what
    //it previously did, and what other Ifpack2 solvers currently do.  So I have
    //moved SetPrecParameters(paramList) into the if-else block above.

    // Apply
    if (InitialGuessIsZero || supportInitialGuess) {
      Tpetra::MultiVector<SC,LO,GO,NO>&       tpX = Utilities::MV2NonConstTpetraMV(X);
      const Tpetra::MultiVector<SC,LO,GO,NO>& tpB = Utilities::MV2TpetraMV(B);
      prec_->apply(tpB, tpX);
    } else {
      typedef Teuchos::ScalarTraits<Scalar> TST;
      RCP<MultiVector> Residual   = Utilities::Residual(*A_, X, B);
      RCP<MultiVector> Correction = MultiVectorFactory::Build(A_->getDomainMap(), X.getNumVectors());

      Tpetra::MultiVector<SC,LO,GO,NO>&       tpX = Utilities::MV2NonConstTpetraMV(*Correction);
      const Tpetra::MultiVector<SC,LO,GO,NO>& tpB = Utilities::MV2TpetraMV(*Residual);

      prec_->apply(tpB, tpX);

      X.update(TST::one(), *Correction, TST::one());
    }
  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node>
  RCP<MueLu::SmootherPrototype<Scalar, LocalOrdinal, GlobalOrdinal, Node> > Ifpack2Smoother<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Copy() const {
    RCP<Ifpack2Smoother> smoother = rcp(new Ifpack2Smoother(*this) );
    smoother->SetParameterList(this->GetParameterList());
    return smoother;
  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node>
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

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node>
  void Ifpack2Smoother<Scalar, LocalOrdinal, GlobalOrdinal, Node>::print(Teuchos::FancyOStream &out, const VerbLevel verbLevel) const {
    MUELU_DESCRIBE;

    if (verbLevel & Parameters0)
      out0 << "Prec. type: " << type_ << std::endl;

    if (verbLevel & Parameters1) {
      out0 << "Parameter list: " << std::endl;
      Teuchos::OSTab tab2(out);
      out << this->GetParameterList();
      out0 << "Overlap: "        << overlap_ << std::endl;
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

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node>
  size_t Ifpack2Smoother<Scalar, LocalOrdinal, GlobalOrdinal, Node>::getNodeSmootherComplexity() const {
    typedef Tpetra::RowMatrix<SC,LO,GO,NO> MatrixType;
    // NOTE: Only works for a subset of Ifpack2's smoothers
    RCP<Ifpack2::Relaxation<MatrixType> > pr     = rcp_dynamic_cast<Ifpack2::Relaxation<MatrixType> >(prec_);
    if(!pr.is_null()) return pr->getNodeSmootherComplexity();

    RCP<Ifpack2::Chebyshev<MatrixType> > pc       = rcp_dynamic_cast<Ifpack2::Chebyshev<MatrixType> >(prec_);
    if(!pc.is_null()) return pc->getNodeSmootherComplexity();

    RCP<Ifpack2::BlockRelaxation<MatrixType> > pb = rcp_dynamic_cast<Ifpack2::BlockRelaxation<MatrixType> >(prec_);
    if(!pb.is_null()) return pb->getNodeSmootherComplexity();

    RCP<Ifpack2::ILUT<MatrixType> > pi            = rcp_dynamic_cast<Ifpack2::ILUT<MatrixType> >(prec_);
    if(!pi.is_null()) return pi->getNodeSmootherComplexity();

    RCP<Ifpack2::RILUK<MatrixType> > pk            = rcp_dynamic_cast<Ifpack2::RILUK<MatrixType> >(prec_);
    if(!pk.is_null()) return pk->getNodeSmootherComplexity();


    return Teuchos::OrdinalTraits<size_t>::invalid();
  }


} // namespace MueLu

#endif // HAVE_MUELU_TPETRA && HAVE_MUELU_IFPACK2
#endif // MUELU_IFPACK2SMOOTHER_DEF_HPP
