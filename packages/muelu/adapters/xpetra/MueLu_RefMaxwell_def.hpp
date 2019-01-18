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
#ifndef MUELU_REFMAXWELL_DEF_HPP
#define MUELU_REFMAXWELL_DEF_HPP

#include <sstream>

#include "MueLu_ConfigDefs.hpp"

#include "Xpetra_Map.hpp"
#include "Xpetra_MatrixMatrix.hpp"
#include "Xpetra_TripleMatrixMultiply.hpp"
#include "Xpetra_CrsMatrixUtils.hpp"
#include "Xpetra_MatrixUtils.hpp"

#include "MueLu_RefMaxwell_decl.hpp"

#include "MueLu_AmalgamationFactory.hpp"
#include "MueLu_RAPFactory.hpp"
#include "MueLu_ThresholdAFilterFactory.hpp"
#include "MueLu_TransPFactory.hpp"
#include "MueLu_SmootherFactory.hpp"

#include "MueLu_CoalesceDropFactory.hpp"
#include "MueLu_CoarseMapFactory.hpp"
#include "MueLu_CoordinatesTransferFactory.hpp"
#include "MueLu_UncoupledAggregationFactory.hpp"
#include "MueLu_TentativePFactory.hpp"
#include "MueLu_Utilities.hpp"

#ifdef HAVE_MUELU_KOKKOS_REFACTOR
#include "MueLu_CoalesceDropFactory_kokkos.hpp"
#include "MueLu_CoarseMapFactory_kokkos.hpp"
#include "MueLu_CoordinatesTransferFactory_kokkos.hpp"
#include "MueLu_UncoupledAggregationFactory_kokkos.hpp"
#include "MueLu_TentativePFactory_kokkos.hpp"
#include "MueLu_Utilities_kokkos.hpp"
#include <Kokkos_Core.hpp>
#include <KokkosSparse_CrsMatrix.hpp>
#endif

#include "MueLu_ZoltanInterface.hpp"
#include "MueLu_Zoltan2Interface.hpp"
#include "MueLu_RepartitionHeuristicFactory.hpp"
#include "MueLu_RepartitionFactory.hpp"
#include "MueLu_RebalanceAcFactory.hpp"
#include "MueLu_RebalanceTransferFactory.hpp"

#include "MueLu_VerbosityLevel.hpp"

#include <MueLu_CreateXpetraPreconditioner.hpp>

#ifdef HAVE_MUELU_CUDA
#include "cuda_profiler_api.h"
#endif


namespace MueLu {

  template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  Teuchos::RCP<const Xpetra::Map<LocalOrdinal,GlobalOrdinal,Node> > RefMaxwell<Scalar,LocalOrdinal,GlobalOrdinal,Node>::getDomainMap() const {
    return SM_Matrix_->getDomainMap();
  }


  template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  Teuchos::RCP<const Xpetra::Map<LocalOrdinal,GlobalOrdinal,Node> > RefMaxwell<Scalar,LocalOrdinal,GlobalOrdinal,Node>::getRangeMap() const {
    return SM_Matrix_->getRangeMap();
  }


  template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void RefMaxwell<Scalar,LocalOrdinal,GlobalOrdinal,Node>::setParameters(Teuchos::ParameterList& list) {

    parameterList_    = list;
    disable_addon_    = list.get("refmaxwell: disable addon",true);
    mode_             = list.get("refmaxwell: mode","additive");
    use_as_preconditioner_ = list.get<bool>("refmaxwell: use as preconditioner");
    dump_matrices_    = list.get("refmaxwell: dump matrices",false);

    if(list.isSublist("refmaxwell: 11list"))
      precList11_     =  list.sublist("refmaxwell: 11list");

    if(list.isSublist("refmaxwell: 22list"))
      precList22_     =  list.sublist("refmaxwell: 22list");

    if(list.isSublist("smoother: params")) {
      smootherList_ = list.sublist("smoother: params");
    }

#if !defined(HAVE_MUELU_KOKKOS_REFACTOR)
    useKokkos_ = false;
#elif defined(HAVE_MUELU_KOKKOS_REFACTOR_USE_BY_DEFAULT)
    useKokkos_ = list.get("use kokkos refactor",true);
#else
    useKokkos_ = list.get("use kokkos refactor",false);
#endif

  }


  template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void RefMaxwell<Scalar,LocalOrdinal,GlobalOrdinal,Node>::compute() {


    typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType realType;

#ifdef HAVE_MUELU_CUDA
    if (parameterList_.get<bool>("refmaxwell: cuda profile setup", false)) cudaProfilerStart();
#endif

    Teuchos::TimeMonitor tmCompute(*Teuchos::TimeMonitor::getNewTimer("MueLu RefMaxwell: compute"));

    std::map<std::string, MsgType> verbMap;
    verbMap["none"]    = None;
    verbMap["low"]     = Low;
    verbMap["medium"]  = Medium;
    verbMap["high"]    = High;
    verbMap["extreme"] = Extreme;
    verbMap["test"]    = Test;

    VerbLevel oldVerbLevel = VerboseObject::GetDefaultVerbLevel();
    std::string verbosityLevel = parameterList_.get<std::string>("verbosity", "medium");

    TEUCHOS_TEST_FOR_EXCEPTION(verbMap.count(verbosityLevel) == 0, Exceptions::RuntimeError,
                               "Invalid verbosity level: \"" << verbosityLevel << "\"");
    VerboseObject::SetDefaultVerbLevel(verbMap[verbosityLevel]);


    bool defaultFilter = false;

    // Remove zero entries from D0 if necessary.
    // In the construction of the prolongator we use the graph of the
    // matrix, so zero entries mess it up.
    if (parameterList_.get<bool>("refmaxwell: filter D0", true) && D0_Matrix_->getNodeMaxNumRowEntries()>2) {
      Level fineLevel;
      fineLevel.SetFactoryManager(null);
      fineLevel.SetLevelID(0);
      fineLevel.Set("A",D0_Matrix_);
      fineLevel.setlib(D0_Matrix_->getDomainMap()->lib());
      // We expect D0 to have entries +-1, so any threshold value will do.
      RCP<ThresholdAFilterFactory> ThreshFact = rcp(new ThresholdAFilterFactory("A",1.0e-8,/*keepDiagonal=*/false,/*expectedNNZperRow=*/2));
      fineLevel.Request("A",ThreshFact.get());
      ThreshFact->Build(fineLevel);
      D0_Matrix_ = fineLevel.Get< RCP<Matrix> >("A",ThreshFact.get());

      // If D0 has too many zeros, maybe SM and M1 do as well.
      defaultFilter = true;
    }

    if (parameterList_.get<bool>("refmaxwell: filter SM", defaultFilter)) {
      RCP<Vector> diag = VectorFactory::Build(SM_Matrix_->getRowMap());
      // find a reasonable absolute value threshold
      SM_Matrix_->getLocalDiagCopy(*diag);
      magnitudeType threshold = 1.0e-8 * diag->normInf();

      Level fineLevel;
      fineLevel.SetFactoryManager(null);
      fineLevel.SetLevelID(0);
      fineLevel.Set("A",SM_Matrix_);
      fineLevel.setlib(SM_Matrix_->getDomainMap()->lib());
      RCP<ThresholdAFilterFactory> ThreshFact = rcp(new ThresholdAFilterFactory("A",threshold,/*keepDiagonal=*/true));
      fineLevel.Request("A",ThreshFact.get());
      ThreshFact->Build(fineLevel);
      SM_Matrix_ = fineLevel.Get< RCP<Matrix> >("A",ThreshFact.get());
    }

    if (parameterList_.get<bool>("refmaxwell: filter M1", defaultFilter)) {
      RCP<Vector> diag = VectorFactory::Build(M1_Matrix_->getRowMap());
      // find a reasonable absolute value threshold
      M1_Matrix_->getLocalDiagCopy(*diag);
      magnitudeType threshold = 1.0e-8 * diag->normInf();

      Level fineLevel;
      fineLevel.SetFactoryManager(null);
      fineLevel.SetLevelID(0);
      fineLevel.Set("A",M1_Matrix_);
      fineLevel.setlib(M1_Matrix_->getDomainMap()->lib());
      RCP<ThresholdAFilterFactory> ThreshFact = rcp(new ThresholdAFilterFactory("A",threshold,/*keepDiagonal=*/true));
      fineLevel.Request("A",ThreshFact.get());
      ThreshFact->Build(fineLevel);
      M1_Matrix_ = fineLevel.Get< RCP<Matrix> >("A",ThreshFact.get());
    }

    // clean rows associated with boundary conditions
    // Find rows with only 1 or 2 nonzero entries, record them in BCrows_.
    // BCrows_[i] is true, iff i is a boundary row
    // BCcols_[i] is true, iff i is a boundary column
#ifdef HAVE_MUELU_KOKKOS_REFACTOR
    if (useKokkos_) {
      BCrowsKokkos_ = Utilities_kokkos::DetectDirichletRows(*SM_Matrix_,Teuchos::ScalarTraits<magnitudeType>::eps(),/*count_twos_as_dirichlet=*/true);
      BCcolsKokkos_ = Utilities_kokkos::DetectDirichletCols(*D0_Matrix_,BCrowsKokkos_);
    } else {
      BCrows_ = Utilities::DetectDirichletRows(*SM_Matrix_,Teuchos::ScalarTraits<magnitudeType>::eps(),/*count_twos_as_dirichlet=*/true);
      BCcols_ = Utilities::DetectDirichletCols(*D0_Matrix_,BCrows_);
    }
#else
    BCrows_ = Utilities::DetectDirichletRows(*SM_Matrix_,Teuchos::ScalarTraits<magnitudeType>::eps(),/*count_twos_as_dirichlet=*/true);
    BCcols_ = Utilities::DetectDirichletCols(*D0_Matrix_,BCrows_);
#endif

    // build nullspace if necessary
    if(Nullspace_ != null) {
      // no need to do anything - nullspace is built
    }
    else if(Nullspace_ == null && Coords_ != null) {
      // normalize coordinates
      typedef typename RealValuedMultiVector::scalar_type realScalarType;
      typedef typename Teuchos::ScalarTraits<realScalarType>::magnitudeType realMagnitudeType;
      Array<realMagnitudeType> norms(Coords_->getNumVectors());
      Coords_->norm2(norms);
      for (size_t i=0;i<Coords_->getNumVectors();i++)
        norms[i] = ((realMagnitudeType)1.0)/norms[i];
      Coords_->scale(norms());
      Nullspace_ = MultiVectorFactory::Build(SM_Matrix_->getRowMap(),Coords_->getNumVectors());

      // Cast coordinates to Scalar so they can be multiplied against D0
#ifdef HAVE_MUELU_KOKKOS_REFACTOR
      RCP<MultiVector> CoordsSC;
      if (useKokkos_)
        CoordsSC = Utilities_kokkos::RealValuedToScalarMultiVector(Coords_);
      else
        CoordsSC = Utilities::RealValuedToScalarMultiVector(Coords_);
#else
      RCP<MultiVector> CoordsSC = Utilities::RealValuedToScalarMultiVector(Coords_);
#endif
      D0_Matrix_->apply(*CoordsSC,*Nullspace_);
    }
    else {
      GetOStream(Errors) << "MueLu::RefMaxwell::compute(): either the nullspace or the nodal coordinates must be provided." << std::endl;
    }

    // Nuke the BC edges in nullspace
#ifdef HAVE_MUELU_KOKKOS_REFACTOR
    if (useKokkos_)
      Utilities_kokkos::ZeroDirichletRows(Nullspace_,BCrowsKokkos_);
    else
      Utilities::ZeroDirichletRows(Nullspace_,BCrows_);
#else
    Utilities::ZeroDirichletRows(Nullspace_,BCrows_);
#endif

    if (dump_matrices_)
      Xpetra::IO<SC, LO, GlobalOrdinal, Node>::Write(std::string("D0_clean.mat"), *D0_Matrix_);

    // build special prolongator for (1,1)-block
    if(P11_==null) {
      // Form A_nodal = D0* M1 D0  (aka TMT_agg)
      Level fineLevel, coarseLevel;
      fineLevel.SetFactoryManager(null);
      coarseLevel.SetFactoryManager(null);
      coarseLevel.SetPreviousLevel(rcpFromRef(fineLevel));
      fineLevel.SetLevelID(0);
      coarseLevel.SetLevelID(1);
      fineLevel.Set("A",M1_Matrix_);
      coarseLevel.Set("P",D0_Matrix_);
      coarseLevel.setlib(M1_Matrix_->getDomainMap()->lib());
      fineLevel.setlib(M1_Matrix_->getDomainMap()->lib());
      coarseLevel.setObjectLabel("RefMaxwell (1,1) A_nodal");
      fineLevel.setObjectLabel("RefMaxwell (1,1) A_nodal");

      RCP<RAPFactory> rapFact = rcp(new RAPFactory());
      ParameterList rapList = *(rapFact->GetValidParameterList());
      rapList.set("transpose: use implicit", parameterList_.get<bool>("transpose: use implicit", false));
      rapList.set("rap: fix zero diagonals", parameterList_.get<bool>("rap: fix zero diagonals", true));
      rapList.set("rap: triple product", parameterList_.get<bool>("rap: triple product", false));
      rapFact->SetParameterList(rapList);

      RCP<TransPFactory> transPFactory;
      if (!parameterList_.get<bool>("transpose: use implicit", false)) {
        transPFactory = rcp(new TransPFactory());
        rapFact->SetFactory("R", transPFactory);
      }

      coarseLevel.Request("A", rapFact.get());

      A_nodal_Matrix_ = coarseLevel.Get< RCP<Matrix> >("A", rapFact.get());

      // build special prolongator
      GetOStream(Runtime0) << "RefMaxwell::compute(): building special prolongator" << std::endl;
      buildProlongator();

#ifdef HAVE_MUELU_KOKKOS_REFACTOR
      if (useKokkos_)
        R11_ = Utilities_kokkos::Transpose(*P11_);
      else
        R11_ = Utilities::Transpose(*P11_);
#else
      R11_ = Utilities::Transpose(*P11_);
#endif
    }

    bool doRebalancing = parameterList_.get<bool>("refmaxwell: subsolves on subcommunicators", false);
    int numProcsAH, numProcsA22;
    {
      // build coarse grid operator for (1,1)-block
      formCoarseMatrix();

#ifdef HAVE_MPI
      int numProcs = SM_Matrix_->getDomainMap()->getComm()->getSize();
      if (doRebalancing && numProcs > 1) {
        GlobalOrdinal globalNumRowsAH = AH_->getRowMap()->getGlobalNumElements();
        GlobalOrdinal globalNumRowsA22 = D0_Matrix_->getDomainMap()->getGlobalNumElements();
        double ratio = parameterList_.get<double>("refmaxwell: ratio AH / A22 subcommunicators", 1.0);
        numProcsAH = numProcs * globalNumRowsAH / (globalNumRowsAH + ratio*globalNumRowsA22);
        numProcsA22 = numProcs * ratio * globalNumRowsA22 / (globalNumRowsAH + ratio*globalNumRowsA22);
        if (numProcsAH + numProcsA22 < numProcs)
          ++numProcsAH;
        if (numProcsAH + numProcsA22 < numProcs)
          ++numProcsA22;
        numProcsAH = std::max(numProcsAH, 1);
        numProcsA22 = std::max(numProcsA22, 1);
      } else
        doRebalancing = false;

      if (doRebalancing) { // rebalance AH
        Teuchos::TimeMonitor tm(*Teuchos::TimeMonitor::getNewTimer("MueLu RefMaxwell: Rebalance AH"));

        Level fineLevel, coarseLevel;
        fineLevel.SetFactoryManager(null);
        coarseLevel.SetFactoryManager(null);
        coarseLevel.SetPreviousLevel(rcpFromRef(fineLevel));
        fineLevel.SetLevelID(0);
        coarseLevel.SetLevelID(1);
        coarseLevel.Set("A",AH_);
        coarseLevel.Set("P",P11_);
        coarseLevel.Set("R",R11_);
        coarseLevel.Set("Coordinates",CoordsH_);
        coarseLevel.Set("number of partitions", numProcsAH);

        coarseLevel.setlib(AH_->getDomainMap()->lib());
        fineLevel.setlib(AH_->getDomainMap()->lib());
        coarseLevel.setObjectLabel("RefMaxwell (1,1)");
        fineLevel.setObjectLabel("RefMaxwell (1,1)");

        // auto repartheurFactory = rcp(new RepartitionHeuristicFactory());
        // ParameterList repartheurParams;
        // repartheurParams.set("repartition: start level",0);
        // repartheurParams.set("repartition: min rows per proc", precList11_.get<int>("repartition: min rows per proc", 1024));
        // repartheurParams.set("repartition: target rows per proc", precList11_.get<int>("repartition: target rows per proc", 0));
        // repartheurParams.set("repartition: max imbalance", precList11_.get<double>("repartition: max imbalance", 1.1));
        // repartheurFactory->SetParameterList(repartheurParams);

        std::string partName = precList11_.get<std::string>("repartition: partitioner", "zoltan2");
        RCP<Factory> partitioner;
        if (partName == "zoltan") {
#ifdef HAVE_MUELU_ZOLTAN
          partitioner = rcp(new ZoltanInterface());
          // NOTE: ZoltanInteface ("zoltan") does not support external parameters through ParameterList
          // partitioner->SetFactory("number of partitions", repartheurFactory);
#else
          throw Exceptions::RuntimeError("Zoltan interface is not available");
#endif
        } else if (partName == "zoltan2") {
#ifdef HAVE_MUELU_ZOLTAN2
          partitioner = rcp(new Zoltan2Interface());
          ParameterList partParams;
          RCP<const ParameterList> partpartParams = rcp(new ParameterList(precList11_.sublist("repartition: params", false)));
          partParams.set("ParameterList", partpartParams);
          partitioner->SetParameterList(partParams);
          // partitioner->SetFactory("number of partitions", repartheurFactory);
#else
          throw Exceptions::RuntimeError("Zoltan2 interface is not available");
#endif
        }

        auto repartFactory = rcp(new RepartitionFactory());
        ParameterList repartParams;
        repartParams.set("repartition: print partition distribution", precList11_.get<bool>("repartition: print partition distribution", false));
        repartParams.set("repartition: remap parts", precList11_.get<bool>("repartition: remap parts", true));
        repartFactory->SetParameterList(repartParams);
        // repartFactory->SetFactory("number of partitions", repartheurFactory);
        repartFactory->SetFactory("Partition", partitioner);

        auto newP = rcp(new RebalanceTransferFactory());
        ParameterList newPparams;
        newPparams.set("type", "Interpolation");
        newPparams.set("repartition: rebalance P and R", precList11_.get<bool>("repartition: rebalance P and R", false));
        newPparams.set("repartition: use subcommunicators", true);
        newPparams.set("repartition: rebalance Nullspace", false);
        newP->SetFactory("Coordinates", NoFactory::getRCP());
        newP->SetParameterList(newPparams);
        newP->SetFactory("Importer", repartFactory);

        // Rebalanced R
        auto newR = rcp(new RebalanceTransferFactory());
        ParameterList newRparams;
        newRparams.set("type", "Restriction");
        newRparams.set("repartition: rebalance P and R", precList11_.get<bool>("repartition: rebalance P and R", false));
        newRparams.set("repartition: use subcommunicators", true);
        newR->SetParameterList(newRparams);
        newR->SetFactory("Importer", repartFactory);

        auto newA = rcp(new RebalanceAcFactory());
        ParameterList rebAcParams;
        rebAcParams.set("repartition: use subcommunicators", true);
        newA->SetParameterList(rebAcParams);
        newA->SetFactory("Importer", repartFactory);

        coarseLevel.Request("R", newR.get());
        coarseLevel.Request("P", newP.get());
        coarseLevel.Request("Importer", repartFactory.get());
        coarseLevel.Request("A", newA.get());
        coarseLevel.Request("Coordinates", newP.get());
        repartFactory->Build(coarseLevel);

        if (!precList11_.get<bool>("repartition: rebalance P and R", false))
          ImporterH_ = coarseLevel.Get< RCP<const Import> >("Importer", repartFactory.get());
        P11_ = coarseLevel.Get< RCP<Matrix> >("P", newP.get());
        R11_ = coarseLevel.Get< RCP<Matrix> >("R", newR.get());
        AH_ = coarseLevel.Get< RCP<Matrix> >("A", newA.get());
        if (!AH_.is_null())
          AH_->setObjectLabel("RefMaxwell (1,1)");
        CoordsH_ = coarseLevel.Get< RCP<RealValuedMultiVector> >("Coordinates", newP.get());
      }
#endif // HAVE_MPI

#ifdef HAVE_MUELU_KOKKOS_REFACTOR
      // This should be taken out again as soon as
      // CoalesceDropFactory_kokkos supports BlockSize > 1 and
      // drop tol != 0.0
      if (useKokkos_ && precList11_.isParameter("aggregation: drop tol") && precList11_.get<double>("aggregation: drop tol") != 0.0) {
        GetOStream(Warnings0) << "RefMaxwell::compute(): Setting \"aggregation: drop tol\". to 0.0, since CoalesceDropFactory_kokkos does not "
                              << "support BlockSize > 1 and drop tol != 0.0" << std::endl;
        precList11_.set("aggregation: drop tol", 0.0);
      }
#endif
      if (!AH_.is_null()) {
        int oldRank = SetProcRankVerbose(AH_->getDomainMap()->getComm()->getRank());
        HierarchyH_ = MueLu::CreateXpetraPreconditioner(AH_, precList11_, CoordsH_);
        SetProcRankVerbose(oldRank);
      }
      VerboseObject::SetDefaultVerbLevel(verbMap[verbosityLevel]);

    }

    {
      GetOStream(Runtime0) << "RefMaxwell::compute(): nuking BC edges of D0" << std::endl;

      D0_Matrix_->resumeFill();
      // Scalar replaceWith = Teuchos::ScalarTraits<SC>::eps();
      Scalar replaceWith = Teuchos::ScalarTraits<SC>::zero();
#ifdef HAVE_MUELU_KOKKOS_REFACTOR
      if (useKokkos_) {
        Utilities_kokkos::ZeroDirichletRows(D0_Matrix_,BCrowsKokkos_,replaceWith);
        Utilities_kokkos::ZeroDirichletCols(D0_Matrix_,BCcolsKokkos_,replaceWith);
      } else {
        Utilities::ZeroDirichletRows(D0_Matrix_,BCrows_,replaceWith);
        Utilities::ZeroDirichletCols(D0_Matrix_,BCcols_,replaceWith);
      }
#else
      Utilities::ZeroDirichletRows(D0_Matrix_,BCrows_,replaceWith);
      Utilities::ZeroDirichletCols(D0_Matrix_,BCcols_,replaceWith);
#endif
      D0_Matrix_->fillComplete(D0_Matrix_->getDomainMap(),D0_Matrix_->getRangeMap());
    }

    {
      GetOStream(Runtime0) << "RefMaxwell::compute(): building MG for (2,2)-block" << std::endl;

      { // build fine grid operator for (2,2)-block, D0* SM D0  (aka TMT)
        Teuchos::TimeMonitor tm(*Teuchos::TimeMonitor::getNewTimer("MueLu RefMaxwell: Build A22"));

        Level fineLevel, coarseLevel;
        fineLevel.SetFactoryManager(null);
        coarseLevel.SetFactoryManager(null);
        coarseLevel.SetPreviousLevel(rcpFromRef(fineLevel));
        fineLevel.SetLevelID(0);
        coarseLevel.SetLevelID(1);
        fineLevel.Set("A",SM_Matrix_);
        coarseLevel.Set("P",D0_Matrix_);
        coarseLevel.Set("Coordinates",Coords_);
        coarseLevel.Set("number of partitions", numProcsA22);

        coarseLevel.setlib(SM_Matrix_->getDomainMap()->lib());
        fineLevel.setlib(SM_Matrix_->getDomainMap()->lib());
        coarseLevel.setObjectLabel("RefMaxwell (2,2)");
        fineLevel.setObjectLabel("RefMaxwell (2,2)");

        RCP<RAPFactory> rapFact = rcp(new RAPFactory());
        ParameterList rapList = *(rapFact->GetValidParameterList());
        rapList.set("transpose: use implicit", false);
        rapList.set("rap: fix zero diagonals", parameterList_.get<bool>("rap: fix zero diagonals", true));
        rapList.set("rap: triple product", parameterList_.get<bool>("rap: triple product", false));
        rapFact->SetParameterList(rapList);

        RCP<TransPFactory> transPFactory;
        transPFactory = rcp(new TransPFactory());
        rapFact->SetFactory("R", transPFactory);

#ifdef HAVE_MPI
        if (doRebalancing) {

          // auto repartheurFactory = rcp(new RepartitionHeuristicFactory());
          // ParameterList repartheurParams;
          // repartheurParams.set("repartition: start level",0);
          // repartheurParams.set("repartition: min rows per proc", precList22_.get<int>("repartition: min rows per proc", 1024));
          // repartheurParams.set("repartition: target rows per proc", precList22_.get<int>("repartition: target rows per proc", 0));
          // repartheurParams.set("repartition: max imbalance", precList22_.get<double>("repartition: max imbalance", 1.1));
          // repartheurFactory->SetParameterList(repartheurParams);
          // repartheurFactory->SetFactory("A", rapFact);

          std::string partName = precList22_.get<std::string>("repartition: partitioner", "zoltan2");
          RCP<Factory> partitioner;
          if (partName == "zoltan") {
#ifdef HAVE_MUELU_ZOLTAN
            partitioner = rcp(new ZoltanInterface());
            partitioner->SetFactory("A", rapFact);
            // partitioner->SetFactory("number of partitions", repartheurFactory);
            // NOTE: ZoltanInteface ("zoltan") does not support external parameters through ParameterList
#else
            throw Exceptions::RuntimeError("Zoltan interface is not available");
#endif
          } else if (partName == "zoltan2") {
#ifdef HAVE_MUELU_ZOLTAN2
            partitioner = rcp(new Zoltan2Interface());
            ParameterList partParams;
            RCP<const ParameterList> partpartParams = rcp(new ParameterList(precList22_.sublist("repartition: params", false)));
            partParams.set("ParameterList", partpartParams);
            partitioner->SetParameterList(partParams);
            partitioner->SetFactory("A", rapFact);
            // partitioner->SetFactory("number of partitions", repartheurFactory);
#else
            throw Exceptions::RuntimeError("Zoltan2 interface is not available");
#endif
          }

          auto repartFactory = rcp(new RepartitionFactory());
          ParameterList repartParams;
          repartParams.set("repartition: print partition distribution", precList22_.get<bool>("repartition: print partition distribution", false));
          repartParams.set("repartition: remap parts", precList22_.get<bool>("repartition: remap parts", true));
          repartParams.set("repartition: remap accept partition", AH_.is_null());
          repartFactory->SetParameterList(repartParams);
          repartFactory->SetFactory("A", rapFact);
          // repartFactory->SetFactory("number of partitions", repartheurFactory);
          repartFactory->SetFactory("Partition", partitioner);

          auto newP = rcp(new RebalanceTransferFactory());
          ParameterList newPparams;
          newPparams.set("type", "Interpolation");
          newPparams.set("repartition: rebalance P and R", precList22_.get<bool>("repartition: rebalance P and R", false));
          newPparams.set("repartition: use subcommunicators", true);
          newPparams.set("repartition: rebalance Nullspace", false);
          newP->SetFactory("Coordinates", NoFactory::getRCP());
          newP->SetParameterList(newPparams);
          newP->SetFactory("Importer", repartFactory);

          // Rebalanced R
          auto newR = rcp(new RebalanceTransferFactory());
          ParameterList newRparams;
          newRparams.set("type", "Restriction");
          newRparams.set("repartition: rebalance P and R", precList22_.get<bool>("repartition: rebalance P and R", false));
          newRparams.set("repartition: use subcommunicators", true);
          newR->SetParameterList(newRparams);
          newR->SetFactory("Importer", repartFactory);
          newR->SetFactory("R", transPFactory);

          auto newA = rcp(new RebalanceAcFactory());
          ParameterList rebAcParams;
          rebAcParams.set("repartition: use subcommunicators", true);
          newA->SetParameterList(rebAcParams);
          newA->SetFactory("A", rapFact);
          newA->SetFactory("Importer", repartFactory);

          coarseLevel.Request("R", newR.get());
          coarseLevel.Request("P", newP.get());
          coarseLevel.Request("Importer", repartFactory.get());
          coarseLevel.Request("A", newA.get());
          coarseLevel.Request("Coordinates", newP.get());
          rapFact->Build(fineLevel,coarseLevel);
          repartFactory->Build(coarseLevel);

          if (!precList22_.get<bool>("repartition: rebalance P and R", false))
            Importer22_ = coarseLevel.Get< RCP<const Import> >("Importer", repartFactory.get());
          D0_Matrix_ = coarseLevel.Get< RCP<Matrix> >("P", newP.get());
          D0_T_Matrix_ = coarseLevel.Get< RCP<Matrix> >("R", newR.get());
          A22_ = coarseLevel.Get< RCP<Matrix> >("A", newA.get());
          if (!A22_.is_null())
            A22_->setObjectLabel("RefMaxwell (2,2)");
          Coords_ = coarseLevel.Get< RCP<RealValuedMultiVector> >("Coordinates", newP.get());
        } else
#endif // HAVE_MPI
        {
          coarseLevel.Request("A", rapFact.get());
          coarseLevel.Request("R", transPFactory.get());

          A22_ = coarseLevel.Get< RCP<Matrix> >("A", rapFact.get());
          A22_->setObjectLabel("RefMaxwell (2,2)");
          D0_T_Matrix_ = coarseLevel.Get< RCP<Matrix> >("R", transPFactory.get());
        }
      }

      if (!A22_.is_null()) {
        int oldRank = SetProcRankVerbose(A22_->getDomainMap()->getComm()->getRank());
        Hierarchy22_ = MueLu::CreateXpetraPreconditioner(A22_, precList22_, Coords_);
        SetProcRankVerbose(oldRank);
      }
      VerboseObject::SetDefaultVerbLevel(verbMap[verbosityLevel]);

    }

    {
      if (parameterList_.isType<std::string>("smoother: type") &&
          parameterList_.get<std::string>("smoother: type") == "hiptmair" &&
          SM_Matrix_->getDomainMap()->lib() == Xpetra::UseTpetra &&
          A22_->getDomainMap()->lib() == Xpetra::UseTpetra &&
          D0_Matrix_->getDomainMap()->lib() == Xpetra::UseTpetra) {
#if defined(HAVE_MUELU_IFPACK2) && (!defined(HAVE_MUELU_EPETRA) || (defined(HAVE_MUELU_INST_DOUBLE_INT_INT)))
        ParameterList hiptmairPreList, hiptmairPostList, smootherPreList, smootherPostList;

        if (smootherList_.isSublist("smoother: pre params"))
          smootherPreList = smootherList_.sublist("smoother: pre params");
        else if (smootherList_.isSublist("smoother: params"))
          smootherPreList = smootherList_.sublist("smoother: params");
        hiptmairPreList.set("hiptmair: smoother type 1",
                            smootherPreList.get<std::string>("hiptmair: smoother type 1", "CHEBYSHEV"));
        hiptmairPreList.set("hiptmair: smoother type 2",
                            smootherPreList.get<std::string>("hiptmair: smoother type 2", "CHEBYSHEV"));
        if(smootherPreList.isSublist("hiptmair: smoother list 1"))
          hiptmairPreList.set("hiptmair: smoother list 1", smootherPreList.sublist("hiptmair: smoother list 1"));
        if(smootherPreList.isSublist("hiptmair: smoother list 2"))
          hiptmairPreList.set("hiptmair: smoother list 2", smootherPreList.sublist("hiptmair: smoother list 2"));
        hiptmairPreList.set("hiptmair: pre or post",
                            smootherPreList.get<std::string>("hiptmair: pre or post", "pre"));
        hiptmairPreList.set("hiptmair: zero starting solution",
                            smootherPreList.get<bool>("hiptmair: zero starting solution", true));

        if (smootherList_.isSublist("smoother: post params"))
          smootherPostList = smootherList_.sublist("smoother: post params");
        else if (smootherList_.isSublist("smoother: params"))
          smootherPostList = smootherList_.sublist("smoother: params");
        hiptmairPostList.set("hiptmair: smoother type 1",
                             smootherPostList.get<std::string>("hiptmair: smoother type 1", "CHEBYSHEV"));
        hiptmairPostList.set("hiptmair: smoother type 2",
                             smootherPostList.get<std::string>("hiptmair: smoother type 2", "CHEBYSHEV"));
        if(smootherPostList.isSublist("hiptmair: smoother list 1"))
          hiptmairPostList.set("hiptmair: smoother list 1", smootherPostList.sublist("hiptmair: smoother list 1"));
        if(smootherPostList.isSublist("hiptmair: smoother list 2"))
          hiptmairPostList.set("hiptmair: smoother list 2", smootherPostList.sublist("hiptmair: smoother list 2"));
        hiptmairPostList.set("hiptmair: pre or post",
                             smootherPostList.get<std::string>("hiptmair: pre or post", "post"));
        hiptmairPostList.set("hiptmair: zero starting solution",
                             smootherPostList.get<bool>("hiptmair: zero starting solution", false));

        typedef Tpetra::RowMatrix<SC, LO, GO, NO> TROW;
        RCP<const TROW > EdgeMatrix = Utilities::Op2NonConstTpetraRow(SM_Matrix_);
        RCP<const TROW > NodeMatrix = Utilities::Op2NonConstTpetraRow(A22_);
        RCP<const TROW > PMatrix = Utilities::Op2NonConstTpetraRow(D0_Matrix_);

        hiptmairPreSmoother_  = rcp( new Ifpack2::Hiptmair<TROW>(EdgeMatrix,NodeMatrix,PMatrix) );
        hiptmairPreSmoother_ -> setParameters(hiptmairPreList);
        hiptmairPreSmoother_ -> initialize();
        hiptmairPreSmoother_ -> compute();
        hiptmairPostSmoother_ = rcp( new Ifpack2::Hiptmair<TROW>(EdgeMatrix,NodeMatrix,PMatrix) );
        hiptmairPostSmoother_ -> setParameters(hiptmairPostList);
        hiptmairPostSmoother_ -> initialize();
        hiptmairPostSmoother_ -> compute();
        useHiptmairSmoothing_ = true;
#else
        throw(Xpetra::Exceptions::RuntimeError("MueLu must be compiled with Ifpack2 for Hiptmair smoothing."));
#endif  // defined(HAVE_MUELU_IFPACK2) && (!defined(HAVE_MUELU_EPETRA) || defined(HAVE_MUELU_INST_DOUBLE_INT_INT))
      } else {
        if (parameterList_.isType<std::string>("smoother: pre type") && parameterList_.isType<std::string>("smoother: post type")) {
          std::string preSmootherType = parameterList_.get<std::string>("smoother: pre type");
          std::string postSmootherType = parameterList_.get<std::string>("smoother: post type");

          ParameterList preSmootherList, postSmootherList;
          if (parameterList_.isSublist("smoother: pre params"))
            preSmootherList = parameterList_.sublist("smoother: pre params");
          if (parameterList_.isSublist("smoother: post params"))
            postSmootherList = parameterList_.sublist("smoother: post params");

          Level level;
          RCP<MueLu::FactoryManagerBase> factoryHandler = rcp(new FactoryManager());
          level.SetFactoryManager(factoryHandler);
          level.SetLevelID(0);
          level.setObjectLabel("RefMaxwell (1,1)");
          level.Set("A",SM_Matrix_);
          level.setlib(SM_Matrix_->getDomainMap()->lib());

          RCP<SmootherPrototype> preSmootherPrototype = rcp(new TrilinosSmoother(preSmootherType, preSmootherList));
          RCP<SmootherFactory> preSmootherFact = rcp(new SmootherFactory(preSmootherPrototype));

          RCP<SmootherPrototype> postSmootherPrototype = rcp(new TrilinosSmoother(postSmootherType, postSmootherList));
          RCP<SmootherFactory> postSmootherFact = rcp(new SmootherFactory(postSmootherPrototype));

          level.Request("PreSmoother",preSmootherFact.get());
          preSmootherFact->Build(level);
          PreSmoother_ = level.Get<RCP<SmootherBase> >("PreSmoother",preSmootherFact.get());

          level.Request("PostSmoother",postSmootherFact.get());
          postSmootherFact->Build(level);
          PostSmoother_ = level.Get<RCP<SmootherBase> >("PostSmoother",postSmootherFact.get());
        } else {
          std::string smootherType = parameterList_.get<std::string>("smoother: type", "CHEBYSHEV");
          Level level;
          RCP<MueLu::FactoryManagerBase> factoryHandler = rcp(new FactoryManager());
          level.SetFactoryManager(factoryHandler);
          level.SetLevelID(0);
          level.setObjectLabel("RefMaxwell (1,1)");
          level.Set("A",SM_Matrix_);
          level.setlib(SM_Matrix_->getDomainMap()->lib());
          RCP<SmootherPrototype> smootherPrototype = rcp(new TrilinosSmoother(smootherType, smootherList_));
          RCP<SmootherFactory> SmootherFact = rcp(new SmootherFactory(smootherPrototype));
          level.Request("PreSmoother",SmootherFact.get());
          SmootherFact->Build(level);
          PreSmoother_ = level.Get<RCP<SmootherBase> >("PreSmoother",SmootherFact.get());
          PostSmoother_ = PreSmoother_;
        }
        useHiptmairSmoothing_ = false;
      }
    }

    // Allocate temporary MultiVectors for solve
    P11res_    = MultiVectorFactory::Build(R11_->getRangeMap(), 1);
    if (!ImporterH_.is_null()) {
      P11resTmp_ = MultiVectorFactory::Build(ImporterH_->getTargetMap(), 1);
      P11xTmp_   = MultiVectorFactory::Build(ImporterH_->getSourceMap(), 1);
      P11x_      = MultiVectorFactory::Build(ImporterH_->getTargetMap(), 1);
    } else
      P11x_      = MultiVectorFactory::Build(P11_->getDomainMap(), 1);
    D0res_     = MultiVectorFactory::Build(D0_T_Matrix_->getRangeMap(), 1);
    if (!Importer22_.is_null()) {
      D0resTmp_ = MultiVectorFactory::Build(Importer22_->getTargetMap(), 1);
      D0xTmp_   = MultiVectorFactory::Build(Importer22_->getSourceMap(), 1);
      D0x_      = MultiVectorFactory::Build(Importer22_->getTargetMap(), 1);
    } else
      D0x_      = MultiVectorFactory::Build(D0_Matrix_->getDomainMap(), 1);
    residual_  = MultiVectorFactory::Build(SM_Matrix_->getDomainMap(), 1);

#ifdef HAVE_MUELU_CUDA
    if (parameterList_.get<bool>("refmaxwell: cuda profile setup", false)) cudaProfilerStop();
#endif

    describe(GetOStream(Runtime0));

    if (dump_matrices_) {
      GetOStream(Runtime0) << "RefMaxwell::compute(): dumping data" << std::endl;
      Xpetra::IO<SC, LO, GlobalOrdinal, Node>::Write(std::string("SM.mat"), *SM_Matrix_);
      Xpetra::IO<SC, LO, GlobalOrdinal, Node>::Write(std::string("M1.mat"), *M1_Matrix_);
      Xpetra::IO<SC, LO, GlobalOrdinal, Node>::Write(std::string("M0inv.mat"), *M0inv_Matrix_);
#ifndef HAVE_MUELU_KOKKOS_REFACTOR
      std::ofstream outBCrows("BCrows.mat");
      std::copy(BCrows_.begin(), BCrows_.end(), std::ostream_iterator<LO>(outBCrows, "\n"));
      std::ofstream outBCcols("BCcols.mat");
      std::copy(BCcols_.begin(), BCcols_.end(), std::ostream_iterator<LO>(outBCcols, "\n"));
#else
      if (useKokkos_) {
        std::ofstream outBCrows("BCrows.mat");
        auto BCrows = Kokkos::create_mirror_view (BCrowsKokkos_);
        Kokkos::deep_copy(BCrows , BCrowsKokkos_);
        for (size_t i = 0; i < BCrows.size(); i++)
          outBCrows << BCrows[i] << "\n";

        std::ofstream outBCcols("BCcols.mat");
        auto BCcols = Kokkos::create_mirror_view (BCcolsKokkos_);
        Kokkos::deep_copy(BCcols , BCcolsKokkos_);
        for (size_t i = 0; i < BCcols.size(); i++)
          outBCcols << BCcols[i] << "\n";
      } else {
        std::ofstream outBCrows("BCrows.mat");
        std::copy(BCrows_.begin(), BCrows_.end(), std::ostream_iterator<LO>(outBCrows, "\n"));
        std::ofstream outBCcols("BCcols.mat");
        std::copy(BCcols_.begin(), BCcols_.end(), std::ostream_iterator<LO>(outBCcols, "\n"));
      }
#endif
      Xpetra::IO<SC, LO, GlobalOrdinal, Node>::Write(std::string("nullspace.mat"), *Nullspace_);
      if (Coords_ != null)
        Xpetra::IO<realType, LO, GlobalOrdinal, Node>::Write(std::string("coords.mat"), *Coords_);
      Xpetra::IO<SC, LO, GlobalOrdinal, Node>::Write(std::string("D0_nuked.mat"), *D0_Matrix_);
      Xpetra::IO<SC, LO, GlobalOrdinal, Node>::Write(std::string("A_nodal.mat"), *A_nodal_Matrix_);
      Xpetra::IO<SC, LO, GO, NO>::Write(std::string("P11.mat"), *P11_);
      if (!AH_.is_null())
        Xpetra::IO<SC, LO, GlobalOrdinal, Node>::Write(std::string("AH.mat"), *AH_);
      if (!A22_.is_null())
        Xpetra::IO<SC, LO, GlobalOrdinal, Node>::Write(std::string("A22.mat"), *A22_);
    }

    VerboseObject::SetDefaultVerbLevel(oldVerbLevel);
  }


  template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void RefMaxwell<Scalar,LocalOrdinal,GlobalOrdinal,Node>::buildProlongator() {
    // The P11 matrix maps node based aggregrates { A_j } to edges { e_i }.
    //
    // The old implementation used
    // P11(i, j*dim+k) = sum_{nodes n_l in e_i intersected with A_j}  0.5 * phi_k(e_i) * P(n_l, A_j)
    // yet the paper gives
    // P11(i, j*dim+k) = sum_{nodes n_l in e_i intersected with A_j}  0.5 * phi_k(e_i)
    // where phi_k is the k-th nullspace vector.
    //
    // The graph of D0 contains the incidence from nodes to edges.
    // The nodal prolongator P maps aggregates to nodes.

    const SC SC_ZERO = Teuchos::ScalarTraits<SC>::zero();
    const SC SC_ONE = Teuchos::ScalarTraits<SC>::one();
    const Scalar half = SC_ONE / (SC_ONE + SC_ONE);
    size_t dim = Nullspace_->getNumVectors();
    size_t numLocalRows = SM_Matrix_->getNodeNumRows();

    // build prolongator: algorithm 1 in the reference paper
    // First, build nodal unsmoothed prolongator using the matrix A_nodal
    RCP<Matrix> P_nodal;
    bool read_P_from_file = parameterList_.get("refmaxwell: read_P_from_file",false);
    if (read_P_from_file) {
      // This permits to read in an ML prolongator, so that we get the same hierarchy.
      // (ML and MueLu typically produce different aggregates.)
      std::string P_filename = parameterList_.get("refmaxwell: P_filename",std::string("P.mat"));
      P_nodal = Xpetra::IO<SC, LO, GO, NO>::Read(P_filename, A_nodal_Matrix_->getDomainMap());
    } else {
      Level fineLevel, coarseLevel;
      fineLevel.SetFactoryManager(null);
      coarseLevel.SetFactoryManager(null);
      coarseLevel.SetPreviousLevel(rcpFromRef(fineLevel));
      fineLevel.SetLevelID(0);
      coarseLevel.SetLevelID(1);
      fineLevel.Set("A",A_nodal_Matrix_);
      fineLevel.Set("Coordinates",Coords_);
      fineLevel.Set("DofsPerNode",1);
      coarseLevel.setlib(A_nodal_Matrix_->getDomainMap()->lib());
      fineLevel.setlib(A_nodal_Matrix_->getDomainMap()->lib());
      coarseLevel.setObjectLabel("RefMaxwell (1,1)");
      fineLevel.setObjectLabel("RefMaxwell (1,1)");

      LocalOrdinal NSdim = 1;
      RCP<MultiVector> nullSpace = MultiVectorFactory::Build(A_nodal_Matrix_->getRowMap(),NSdim);
      nullSpace->putScalar(SC_ONE);
      fineLevel.Set("Nullspace",nullSpace);

      RCP<AmalgamationFactory> amalgFact = rcp(new AmalgamationFactory());
#ifdef HAVE_MUELU_KOKKOS_REFACTOR
      RCP<Factory> dropFact, UncoupledAggFact, coarseMapFact, TentativePFact, Tfact;
      if (useKokkos_) {
        dropFact = rcp(new CoalesceDropFactory_kokkos());
        UncoupledAggFact = rcp(new UncoupledAggregationFactory_kokkos());
        coarseMapFact = rcp(new CoarseMapFactory_kokkos());
        TentativePFact = rcp(new TentativePFactory_kokkos());
        Tfact = rcp(new CoordinatesTransferFactory_kokkos());
      } else {
        dropFact = rcp(new CoalesceDropFactory());
        UncoupledAggFact = rcp(new UncoupledAggregationFactory());
        coarseMapFact = rcp(new CoarseMapFactory());
        TentativePFact = rcp(new TentativePFactory());
        Tfact = rcp(new CoordinatesTransferFactory());
      }
#else
      RCP<CoalesceDropFactory> dropFact = rcp(new CoalesceDropFactory());
      RCP<UncoupledAggregationFactory> UncoupledAggFact = rcp(new UncoupledAggregationFactory());
      RCP<CoarseMapFactory> coarseMapFact = rcp(new CoarseMapFactory());
      RCP<TentativePFactory> TentativePFact = rcp(new TentativePFactory());
      RCP<CoordinatesTransferFactory> Tfact = rcp(new CoordinatesTransferFactory());
#endif
      dropFact->SetFactory("UnAmalgamationInfo", amalgFact);
      double dropTol = parameterList_.get("aggregation: drop tol",0.0);
      dropFact->SetParameter("aggregation: drop tol",Teuchos::ParameterEntry(dropTol));

      UncoupledAggFact->SetFactory("Graph", dropFact);

      coarseMapFact->SetFactory("Aggregates", UncoupledAggFact);

      TentativePFact->SetFactory("Aggregates", UncoupledAggFact);
      TentativePFact->SetFactory("UnAmalgamationInfo", amalgFact);
      TentativePFact->SetFactory("CoarseMap", coarseMapFact);

      Tfact->SetFactory("Aggregates", UncoupledAggFact);
      Tfact->SetFactory("CoarseMap", coarseMapFact);

      coarseLevel.Request("P",TentativePFact.get());
      coarseLevel.Request("Coordinates",Tfact.get());

      coarseLevel.Get("P",P_nodal,TentativePFact.get());
      coarseLevel.Get("Coordinates",CoordsH_,Tfact.get());
    }
    if (dump_matrices_)
      Xpetra::IO<SC, LO, GO, NO>::Write(std::string("P_nodal.mat"), *P_nodal);

    RCP<CrsMatrix> D0Crs = rcp_dynamic_cast<CrsMatrixWrap>(D0_Matrix_)->getCrsMatrix();

    // Import off-rank rows of P_nodal into P_nodal_imported
    RCP<CrsMatrix> P_nodal_imported;
    int numProcs = P_nodal->getDomainMap()->getComm()->getSize();
    if (numProcs > 1) {
      RCP<CrsMatrixWrap> P_nodal_temp;
      RCP<const Map> targetMap = D0Crs->getColMap();
      P_nodal_temp = rcp(new CrsMatrixWrap(targetMap, 0));
      RCP<const Import> importer = D0Crs->getCrsGraph()->getImporter();
      P_nodal_temp->doImport(*P_nodal, *importer, Xpetra::INSERT);
      P_nodal_temp->fillComplete(rcp_dynamic_cast<CrsMatrixWrap>(P_nodal)->getCrsMatrix()->getDomainMap(),
                                 rcp_dynamic_cast<CrsMatrixWrap>(P_nodal)->getCrsMatrix()->getRangeMap());
      P_nodal_imported = P_nodal_temp->getCrsMatrix();
      if (dump_matrices_)
        Xpetra::IO<SC, LO, GO, NO>::Write(std::string("P_nodal_imported.mat"), *P_nodal_temp);
    } else
      P_nodal_imported = rcp_dynamic_cast<CrsMatrixWrap>(P_nodal)->getCrsMatrix();

#ifdef HAVE_MUELU_KOKKOS_REFACTOR
    if (useKokkos_) {

      using ATS        = Kokkos::ArithTraits<SC>;
      using range_type = Kokkos::RangePolicy<LO, typename NO::execution_space>;

      typedef typename Matrix::local_matrix_type KCRS;
      typedef typename KCRS::device_type device_t;
      typedef typename KCRS::StaticCrsGraphType graph_t;
      typedef typename graph_t::row_map_type::non_const_type lno_view_t;
      typedef typename graph_t::entries_type::non_const_type lno_nnz_view_t;
      typedef typename KCRS::values_type::non_const_type scalar_view_t;

      // Get data out of P_nodal_imported and D0.
      auto localP = P_nodal_imported->getLocalMatrix();
      auto localD0 = D0_Matrix_->getLocalMatrix();

      // Which algorithm should we use for the construction of the special prolongator?
      // Option "mat-mat":
      //   Multiply D0 * P_nodal, take graph, blow up the domain space and compute the entries.
      std::string defaultAlgo = "mat-mat";
      std::string algo = parameterList_.get("refmaxwell: prolongator compute algorithm",defaultAlgo);

      if (algo == "mat-mat") {
        RCP<Matrix> D0_P_nodal = MatrixFactory::Build(SM_Matrix_->getRowMap(),0);
        Xpetra::MatrixMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>::Multiply(*D0_Matrix_,false,*P_nodal,false,*D0_P_nodal,true,true);

        // Get data out of D0*P.
        auto localD0P = D0_P_nodal->getLocalMatrix();

        // Create the matrix object
        RCP<Map> blockColMap    = Xpetra::MapFactory<LO,GO,NO>::Build(P_nodal_imported->getColMap(), dim);
        RCP<Map> blockDomainMap = Xpetra::MapFactory<LO,GO,NO>::Build(P_nodal->getDomainMap(), dim);

        lno_view_t P11rowptr("P11_rowptr", numLocalRows+1);
        lno_nnz_view_t P11colind("P11_colind",dim*localD0P.graph.entries.size());
        scalar_view_t P11vals("P11_vals",dim*localD0P.graph.entries.size());

        // adjust rowpointer
        Kokkos::parallel_for("MueLu:RefMaxwell::buildProlongator_adjustRowptr", range_type(0,numLocalRows+1),
                             KOKKOS_LAMBDA(const size_t i) {
                               P11rowptr(i) = dim*localD0P.graph.row_map(i);
                             });

        // adjust column indices
        Kokkos::parallel_for("MueLu:RefMaxwell::buildProlongator_adjustColind", range_type(0,localD0P.graph.entries.size()),
                             KOKKOS_LAMBDA(const size_t jj) {
                               for (size_t k = 0; k < dim; k++) {
                                 P11colind(dim*jj+k) = dim*localD0P.graph.entries(jj)+k;
                                 P11vals(dim*jj+k) = SC_ZERO;
                               }
                             });

        auto localNullspace = Nullspace_->template getLocalView<device_t>();

        // enter values
        if (D0_Matrix_->getNodeMaxNumRowEntries()>2) {
          // The matrix D0 has too many entries per row.
          // Therefore we need to check whether its entries are actually non-zero.
          // This is the case for the matrices built by MiniEM.
          GetOStream(Warnings0) << "RefMaxwell::buildProlongator(): D0 matrix has more than 2 entries per row. Taking inefficient code path." << std::endl;

          magnitudeType tol = Teuchos::ScalarTraits<magnitudeType>::eps();

          Kokkos::parallel_for("MueLu:RefMaxwell::buildProlongator_enterValues_D0wZeros", range_type(0,numLocalRows),
                               KOKKOS_LAMBDA(const size_t i) {
                                 for (size_t ll = localD0.graph.row_map(i); ll < localD0.graph.row_map(i+1); ll++) {
                                   LO l = localD0.graph.entries(ll);
                                   SC p = localD0.values(ll);
                                   if (ATS::magnitude(p) < tol)
                                     continue;
                                   for (size_t jj = localP.graph.row_map(l); jj < localP.graph.row_map(l+1); jj++) {
                                     LO j = localP.graph.entries(jj);
                                     SC v = localP.values(jj);
                                     for (size_t k = 0; k < dim; k++) {
                                       LO jNew = dim*j+k;
                                       SC n = localNullspace(i,k);
                                       size_t m;
                                       for (m = P11rowptr(i); m < P11rowptr(i+1); m++)
                                         if (P11colind(m) == jNew)
                                           break;
#if defined(HAVE_MUELU_DEBUG) && !defined(HAVE_MUELU_CUDA)
                                       TEUCHOS_ASSERT_EQUALITY(P11colind(m),jNew);
#endif
                                       P11vals(m) += half * v * n;
                                     }
                                   }
                                 }
                               });

        } else {
          Kokkos::parallel_for("MueLu:RefMaxwell::buildProlongator_enterValues", range_type(0,numLocalRows),
                               KOKKOS_LAMBDA(const size_t i) {
                                 for (size_t ll = localD0.graph.row_map(i); ll < localD0.graph.row_map(i+1); ll++) {
                                   LO l = localD0.graph.entries(ll);
                                   for (size_t jj = localP.graph.row_map(l); jj < localP.graph.row_map(l+1); jj++) {
                                     LO j = localP.graph.entries(jj);
                                     SC v = localP.values(jj);
                                     for (size_t k = 0; k < dim; k++) {
                                       LO jNew = dim*j+k;
                                       SC n = localNullspace(i,k);
                                       size_t m;
                                       for (m = P11rowptr(i); m < P11rowptr(i+1); m++)
                                         if (P11colind(m) == jNew)
                                           break;
#if defined(HAVE_MUELU_DEBUG) && !defined(HAVE_MUELU_CUDA)
                                       TEUCHOS_ASSERT_EQUALITY(P11colind(m),jNew);
#endif
                                       P11vals(m) += half * v * n;
                                     }
                                   }
                                 }
                               });
        }

        P11_ = rcp(new CrsMatrixWrap(SM_Matrix_->getRowMap(), blockColMap, 0, Xpetra::StaticProfile));
        RCP<CrsMatrix> P11Crs = rcp_dynamic_cast<CrsMatrixWrap>(P11_)->getCrsMatrix();
        P11Crs->setAllValues(P11rowptr, P11colind, P11vals);
        P11Crs->expertStaticFillComplete(blockDomainMap, SM_Matrix_->getRangeMap());

      } else
        TEUCHOS_TEST_FOR_EXCEPTION(false,std::invalid_argument,algo << " is not a valid option for \"refmaxwell: prolongator compute algorithm\"");
    } else {
#endif // ifdef(HAVE_MUELU_KOKKOS_REFACTOR)

    // get nullspace vectors
    ArrayRCP<ArrayRCP<const SC> > nullspaceRCP(dim);
    ArrayRCP<ArrayView<const SC> > nullspace(dim);
    for(size_t i=0; i<dim; i++) {
      nullspaceRCP[i] = Nullspace_->getData(i);
      nullspace[i] = nullspaceRCP[i]();
    }

    // Get data out of P_nodal_imported and D0.
    ArrayRCP<const size_t>      Prowptr_RCP, D0rowptr_RCP;
    ArrayRCP<const LO>          Pcolind_RCP, D0colind_RCP;
    ArrayRCP<const SC>          Pvals_RCP, D0vals_RCP;
    ArrayRCP<size_t>            P11rowptr_RCP;
    ArrayRCP<LO>                P11colind_RCP;
    ArrayRCP<SC>                P11vals_RCP;

    P_nodal_imported->getAllValues(Prowptr_RCP, Pcolind_RCP, Pvals_RCP);
    rcp_dynamic_cast<CrsMatrixWrap>(D0_Matrix_)->getCrsMatrix()->getAllValues(D0rowptr_RCP, D0colind_RCP, D0vals_RCP);

    // For efficiency
    // Refers to an issue where Teuchos::ArrayRCP::operator[] may be
    // slower than Teuchos::ArrayView::operator[].
    ArrayView<const size_t>     Prowptr, D0rowptr;
    ArrayView<const LO>         Pcolind, D0colind;
    ArrayView<const SC>         Pvals, D0vals;
    Prowptr  = Prowptr_RCP();   Pcolind  = Pcolind_RCP();   Pvals = Pvals_RCP();
    D0rowptr = D0rowptr_RCP();  D0colind = D0colind_RCP();  D0vals = D0vals_RCP();

    // Which algorithm should we use for the construction of the special prolongator?
    // Option "mat-mat":
    //   Multiply D0 * P_nodal, take graph, blow up the domain space and compute the entries.
    // Option "gustavson":
    //   Loop over D0, P and nullspace and allocate directly. (Gustavson-like)
    //   More efficient, but only available for serial node.
    std::string defaultAlgo = "gustavson";
    std::string algo = parameterList_.get("refmaxwell: prolongator compute algorithm",defaultAlgo);

    if (algo == "mat-mat") {
      RCP<Matrix> D0_P_nodal = MatrixFactory::Build(SM_Matrix_->getRowMap(),0);
      Xpetra::MatrixMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>::Multiply(*D0_Matrix_,false,*P_nodal,false,*D0_P_nodal,true,true);

      // Get data out of D0*P.
      ArrayRCP<const size_t>      D0Prowptr_RCP;
      ArrayRCP<const LO>          D0Pcolind_RCP;
      ArrayRCP<const SC>          D0Pvals_RCP;
      rcp_dynamic_cast<CrsMatrixWrap>(D0_P_nodal)->getCrsMatrix()->getAllValues(D0Prowptr_RCP, D0Pcolind_RCP, D0Pvals_RCP);

      // For efficiency
      // Refers to an issue where Teuchos::ArrayRCP::operator[] may be
      // slower than Teuchos::ArrayView::operator[].
      ArrayView<const size_t>     D0Prowptr;
      ArrayView<const LO>         D0Pcolind;
      D0Prowptr = D0Prowptr_RCP(); D0Pcolind = D0Pcolind_RCP();

      // Create the matrix object
      RCP<Map> blockColMap    = Xpetra::MapFactory<LO,GO,NO>::Build(P_nodal_imported->getColMap(), dim);
      RCP<Map> blockDomainMap = Xpetra::MapFactory<LO,GO,NO>::Build(P_nodal->getDomainMap(), dim);
      P11_ = rcp(new CrsMatrixWrap(SM_Matrix_->getRowMap(), blockColMap, 0, Xpetra::StaticProfile));
      RCP<CrsMatrix> P11Crs = rcp_dynamic_cast<CrsMatrixWrap>(P11_)->getCrsMatrix();
      P11Crs->allocateAllValues(dim*D0Pcolind.size(), P11rowptr_RCP, P11colind_RCP, P11vals_RCP);

      ArrayView<size_t> P11rowptr = P11rowptr_RCP();
      ArrayView<LO>     P11colind = P11colind_RCP();
      ArrayView<SC>     P11vals   = P11vals_RCP();

      // adjust rowpointer
      for (size_t i = 0; i < numLocalRows+1; i++) {
        P11rowptr[i] = dim*D0Prowptr[i];
      }

      // adjust column indices
      size_t nnz = 0;
      for (size_t jj = 0; jj < (size_t) D0Pcolind.size(); jj++)
        for (size_t k = 0; k < dim; k++) {
          P11colind[nnz] = dim*D0Pcolind[jj]+k;
          P11vals[nnz] = SC_ZERO;
          nnz++;
        }

      // enter values
      if (D0_Matrix_->getNodeMaxNumRowEntries()>2) {
        // The matrix D0 has too many entries per row.
        // Therefore we need to check whether its entries are actually non-zero.
        // This is the case for the matrices built by MiniEM.
        GetOStream(Warnings0) << "RefMaxwell::buildProlongator(): D0 matrix has more than 2 entries per row. Taking inefficient code path." << std::endl;

        magnitudeType tol = Teuchos::ScalarTraits<magnitudeType>::eps();
        for (size_t i = 0; i < numLocalRows; i++) {
          for (size_t ll = D0rowptr[i]; ll < D0rowptr[i+1]; ll++) {
            LO l = D0colind[ll];
            SC p = D0vals[ll];
            if (Teuchos::ScalarTraits<Scalar>::magnitude(p) < tol)
              continue;
            for (size_t jj = Prowptr[l]; jj < Prowptr[l+1]; jj++) {
              LO j = Pcolind[jj];
              SC v = Pvals[jj];
              for (size_t k = 0; k < dim; k++) {
                LO jNew = dim*j+k;
                SC n = nullspace[k][i];
                size_t m;
                for (m = P11rowptr[i]; m < P11rowptr[i+1]; m++)
                  if (P11colind[m] == jNew)
                    break;
#ifdef HAVE_MUELU_DEBUG
                TEUCHOS_ASSERT_EQUALITY(P11colind[m],jNew);
#endif
                  P11vals[m] += half * v * n;
              }
            }
          }
        }
      } else {
        // enter values
        for (size_t i = 0; i < numLocalRows; i++) {
          for (size_t ll = D0rowptr[i]; ll < D0rowptr[i+1]; ll++) {
            LO l = D0colind[ll];
            for (size_t jj = Prowptr[l]; jj < Prowptr[l+1]; jj++) {
              LO j = Pcolind[jj];
              SC v = Pvals[jj];
              for (size_t k = 0; k < dim; k++) {
                LO jNew = dim*j+k;
                SC n = nullspace[k][i];
                size_t m;
                for (m = P11rowptr[i]; m < P11rowptr[i+1]; m++)
                  if (P11colind[m] == jNew)
                    break;
#ifdef HAVE_MUELU_DEBUG
                TEUCHOS_ASSERT_EQUALITY(P11colind[m],jNew);
#endif
                  P11vals[m] += half * v * n;
              }
            }
          }
        }
      }

      P11Crs->setAllValues(P11rowptr_RCP, P11colind_RCP, P11vals_RCP);
      P11Crs->expertStaticFillComplete(blockDomainMap, SM_Matrix_->getRangeMap());

    } else if (algo == "gustavson") {

      LO maxP11col = dim * P_nodal_imported->getColMap()->getMaxLocalIndex();
      const size_t ST_INVALID = Teuchos::OrdinalTraits<LO>::invalid();
      Array<size_t> P11_status(dim*maxP11col, ST_INVALID);
      // This is ad-hoc and should maybe be replaced with some better heuristics.
      size_t nnz_alloc = dim*D0vals_RCP.size();

      // Create the matrix object
      RCP<Map> blockColMap    = Xpetra::MapFactory<LO,GO,NO>::Build(P_nodal_imported->getColMap(), dim);
      RCP<Map> blockDomainMap = Xpetra::MapFactory<LO,GO,NO>::Build(P_nodal->getDomainMap(), dim);
      P11_ = rcp(new CrsMatrixWrap(SM_Matrix_->getRowMap(), blockColMap, 0, Xpetra::StaticProfile));
      RCP<CrsMatrix> P11Crs = rcp_dynamic_cast<CrsMatrixWrap>(P11_)->getCrsMatrix();
      P11Crs->allocateAllValues(nnz_alloc, P11rowptr_RCP, P11colind_RCP, P11vals_RCP);

      ArrayView<size_t> P11rowptr = P11rowptr_RCP();
      ArrayView<LO>     P11colind = P11colind_RCP();
      ArrayView<SC>     P11vals   = P11vals_RCP();

      size_t nnz;
      if (D0_Matrix_->getNodeMaxNumRowEntries()>2) {
        // The matrix D0 has too many entries per row.
        // Therefore we need to check whether its entries are actually non-zero.
        // This is the case for the matrices built by MiniEM.
        GetOStream(Warnings0) << "RefMaxwell::buildProlongator(): D0 matrix has more than 2 entries per row. Taking inefficient code path." << std::endl;

        magnitudeType tol = Teuchos::ScalarTraits<magnitudeType>::eps();
        nnz = 0;
        size_t nnz_old = 0;
        for (size_t i = 0; i < numLocalRows; i++) {
          P11rowptr[i] = nnz;
          for (size_t ll = D0rowptr[i]; ll < D0rowptr[i+1]; ll++) {
            LO l = D0colind[ll];
            SC p = D0vals[ll];
            if (Teuchos::ScalarTraits<Scalar>::magnitude(p) < tol)
              continue;
            for (size_t jj = Prowptr[l]; jj < Prowptr[l+1]; jj++) {
              LO j = Pcolind[jj];
              SC v = Pvals[jj];
              for (size_t k = 0; k < dim; k++) {
                LO jNew = dim*j+k;
                SC n = nullspace[k][i];
                // do we already have an entry for (i, jNew)?
                if (P11_status[jNew] == ST_INVALID || P11_status[jNew] < nnz_old) {
                  P11_status[jNew] = nnz;
                  P11colind[nnz] = jNew;
                  P11vals[nnz] = half * v * n;
                  // or should it be
                  // P11vals[nnz] = half * n;
                  nnz++;
                } else {
                  P11vals[P11_status[jNew]] += half * v * n;
                  // or should it be
                  // P11vals[P11_status[jNew]] += half * n;
                }
              }
            }
          }
          nnz_old = nnz;
        }
        P11rowptr[numLocalRows] = nnz;
      } else {
        nnz = 0;
        size_t nnz_old = 0;
        for (size_t i = 0; i < numLocalRows; i++) {
          P11rowptr[i] = nnz;
          for (size_t ll = D0rowptr[i]; ll < D0rowptr[i+1]; ll++) {
            LO l = D0colind[ll];
            for (size_t jj = Prowptr[l]; jj < Prowptr[l+1]; jj++) {
              LO j = Pcolind[jj];
              SC v = Pvals[jj];
              for (size_t k = 0; k < dim; k++) {
                LO jNew = dim*j+k;
                SC n = nullspace[k][i];
                // do we already have an entry for (i, jNew)?
                if (P11_status[jNew] == ST_INVALID || P11_status[jNew] < nnz_old) {
                  P11_status[jNew] = nnz;
                  P11colind[nnz] = jNew;
                  P11vals[nnz] = half * v * n;
                  // or should it be
                  // P11vals[nnz] = half * n;
                  nnz++;
                } else {
                  P11vals[P11_status[jNew]] += half * v * n;
                  // or should it be
                  // P11vals[P11_status[jNew]] += half * n;
                }
              }
            }
          }
          nnz_old = nnz;
        }
        P11rowptr[numLocalRows] = nnz;
      }

      if (blockDomainMap->lib() == Xpetra::UseTpetra) {
        // Downward resize
        // - Cannot resize for Epetra, as it checks for same pointers
        // - Need to resize for Tpetra, as it checks ().size() == P11rowptr[numLocalRows]
        P11vals_RCP.resize(nnz);
        P11colind_RCP.resize(nnz);
      }

      P11Crs->setAllValues(P11rowptr_RCP, P11colind_RCP, P11vals_RCP);
      P11Crs->expertStaticFillComplete(blockDomainMap, SM_Matrix_->getRangeMap());
    } else
      TEUCHOS_TEST_FOR_EXCEPTION(false,std::invalid_argument,algo << " is not a valid option for \"refmaxwell: prolongator compute algorithm\"");
#ifdef HAVE_MUELU_KOKKOS_REFACTOR
    }
#endif
  }

  template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void RefMaxwell<Scalar,LocalOrdinal,GlobalOrdinal,Node>::formCoarseMatrix() {
    Teuchos::TimeMonitor tm(*Teuchos::TimeMonitor::getNewTimer("MueLu RefMaxwell: Build coarse (1,1) matrix"));
    
    // coarse matrix for P11* (M1 + D1* M2 D1) P11
    RCP<Matrix> Matrix1 = MatrixFactory::Build(P11_->getDomainMap(),0);
    if (parameterList_.get<bool>("rap: triple product", false) == false) {
      RCP<Matrix> C = MatrixFactory::Build(SM_Matrix_->getRowMap(),0);
      // construct (M1 + D1* M2 D1) P11
      Xpetra::MatrixMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>::Multiply(*SM_Matrix_,false,*P11_,false,*C,true,true);
      // construct P11* (M1 + D1* M2 D1) P11
      Xpetra::MatrixMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>::Multiply(*R11_,false,*C,false,*Matrix1,true,true);
    } else {
      Xpetra::TripleMatrixMultiply<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
        MultiplyRAP(*P11_, true, *SM_Matrix_, false, *P11_, false, *Matrix1, true, true);
    }
    if (parameterList_.get<bool>("rap: fix zero diagonals", true)) {
      const double threshold = parameterList_.get<double>("rap: fix zero diagonals threshold", Teuchos::ScalarTraits<double>::eps());
      Xpetra::MatrixUtils<SC,LO,GO,NO>::CheckRepairMainDiagonal(Matrix1, true, GetOStream(Warnings1), threshold);
    }

    if(disable_addon_==true) {
      // if add-on is not chosen
      AH_=Matrix1;
    }
    else {
      Teuchos::TimeMonitor tmAddon(*Teuchos::TimeMonitor::getNewTimer("MueLu RefMaxwell: Build coarse addon matrix"));
      // catch a failure
      TEUCHOS_TEST_FOR_EXCEPTION(M0inv_Matrix_==Teuchos::null,std::invalid_argument,
                                 "MueLu::RefMaxwell::formCoarseMatrix(): Inverse of "
                                 "lumped mass matrix required for add-on (i.e. M0inv_Matrix is null)");

      // coarse matrix for add-on, i.e P11* (M1 D0 M0inv D0* M1) P11
      RCP<Matrix> Zaux = MatrixFactory::Build(M1_Matrix_->getRowMap(),0);
      RCP<Matrix> Z = MatrixFactory::Build(D0_Matrix_->getDomainMap(),0);

      // construct Zaux = M1 P11
      Xpetra::MatrixMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>::Multiply(*M1_Matrix_,false,*P11_,false,*Zaux,true,true);
      // construct Z = D0* M1 P11 = D0* Zaux
      Xpetra::MatrixMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>::Multiply(*D0_Matrix_,true,*Zaux,false,*Z,true,true);

      // construct Z* M0inv Z
      RCP<Matrix> Matrix2 = MatrixFactory::Build(Z->getDomainMap(),0);
      if (M0inv_Matrix_->getGlobalMaxNumRowEntries()<=1) {
        // We assume that if M0inv has at most one entry per row then
        // these are all diagonal entries.
#ifdef HAVE_MUELU_KOKKOS_REFACTOR
        RCP<Matrix> ZT;
        if (useKokkos_)
          ZT = Utilities_kokkos::Transpose(*Z);
        else
          ZT = Utilities::Transpose(*Z);
#else
        RCP<Matrix> ZT = Utilities::Transpose(*Z);
#endif
        RCP<Vector> diag = VectorFactory::Build(M0inv_Matrix_->getRowMap());
        M0inv_Matrix_->getLocalDiagCopy(*diag);
        if (Z->getRowMap()->isSameAs(*(diag->getMap())))
          Z->leftScale(*diag);
        else {
          RCP<Import> importer = ImportFactory::Build(diag->getMap(),Z->getRowMap());
          RCP<Vector> diag2 = VectorFactory::Build(Z->getRowMap());
          diag2->doImport(*diag,*importer,Xpetra::INSERT);
          Z->leftScale(*diag2);
        }
        Xpetra::MatrixMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>::Multiply(*ZT,false,*Z,false,*Matrix2,true,true);
      } else if (parameterList_.get<bool>("rap: triple product", false) == false) {
        RCP<Matrix> C2 = MatrixFactory::Build(M0inv_Matrix_->getRowMap(),0);
        // construct C2 = M0inv Z
        Xpetra::MatrixMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>::Multiply(*M0inv_Matrix_,false,*Z,false,*C2,true,true);
        // construct Matrix2 = Z* M0inv Z = Z* C2
        Xpetra::MatrixMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>::Multiply(*Z,true,*C2,false,*Matrix2,true,true);
      } else {
        // construct Matrix2 = Z* M0inv Z
        Xpetra::TripleMatrixMultiply<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
          MultiplyRAP(*Z, true, *M0inv_Matrix_, false, *Z, false, *Matrix2, true, true);
      }
      // add matrices together
      Xpetra::MatrixMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>::TwoMatrixAdd(*Matrix1,false,(Scalar)1.0,*Matrix2,false,(Scalar)1.0,AH_,GetOStream(Runtime0));
      AH_->fillComplete();
    }

    // set fixed block size for vector nodal matrix
    size_t dim = Nullspace_->getNumVectors();
    AH_->SetFixedBlockSize(dim);
    AH_->setObjectLabel("RefMaxwell (1,1)");

  }


  template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void RefMaxwell<Scalar,LocalOrdinal,GlobalOrdinal,Node>::resetMatrix(RCP<Matrix> SM_Matrix_new) {
    SM_Matrix_ = SM_Matrix_new;
  }


  template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void RefMaxwell<Scalar,LocalOrdinal,GlobalOrdinal,Node>::solveH() const {

    // iterate on coarse (1, 1) block
    if (!ImporterH_.is_null()) {
      Teuchos::TimeMonitor tmH(*Teuchos::TimeMonitor::getNewTimer("MueLu RefMaxwell: import coarse (1,1)"));
      P11resTmp_->doImport(*P11res_, *ImporterH_, Xpetra::INSERT);
      P11res_.swap(P11resTmp_);
    }
    if (!AH_.is_null()) {
      Teuchos::TimeMonitor tmH(*Teuchos::TimeMonitor::getNewTimer("MueLu RefMaxwell: solve coarse (1,1)"));

      RCP<const Map> origXMap = P11x_->getMap();
      RCP<const Map> origRhsMap = P11res_->getMap();

      // Replace maps with maps with a subcommunicator
      P11res_->replaceMap(AH_->getRangeMap());
      P11x_  ->replaceMap(AH_->getDomainMap());
      HierarchyH_->Iterate(*P11res_, *P11x_, 1, true);
      P11x_  ->replaceMap(origXMap);
      P11res_->replaceMap(origRhsMap);
    }
    if (!ImporterH_.is_null()) {
      Teuchos::TimeMonitor tmH(*Teuchos::TimeMonitor::getNewTimer("MueLu RefMaxwell: export coarse (1,1)"));
      P11xTmp_->doExport(*P11x_, *ImporterH_, Xpetra::INSERT);
      P11x_.swap(P11xTmp_);
    }
  }


  template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void RefMaxwell<Scalar,LocalOrdinal,GlobalOrdinal,Node>::solve22() const {
    // iterate on (2, 2) block
    if (!Importer22_.is_null()) {
      Teuchos::TimeMonitor tm22(*Teuchos::TimeMonitor::getNewTimer("MueLu RefMaxwell: import (2,2)"));
      D0resTmp_->doImport(*D0res_, *Importer22_, Xpetra::INSERT);
      D0res_.swap(D0resTmp_);
    }
    if (!A22_.is_null()) {
      Teuchos::TimeMonitor tm22(*Teuchos::TimeMonitor::getNewTimer("MueLu RefMaxwell: solve (2,2)"));

      RCP<const Map> origXMap = D0x_->getMap();
      RCP<const Map> origRhsMap = D0res_->getMap();

      // Replace maps with maps with a subcommunicator
      D0res_->replaceMap(A22_->getRangeMap());
      D0x_  ->replaceMap(A22_->getDomainMap());
      Hierarchy22_->Iterate(*D0res_, *D0x_, 1, true);
      D0x_  ->replaceMap(origXMap);
      D0res_->replaceMap(origRhsMap);
    }
    if (!Importer22_.is_null()) {
      Teuchos::TimeMonitor tm22(*Teuchos::TimeMonitor::getNewTimer("MueLu RefMaxwell: export (2,2)"));
      D0xTmp_->doExport(*D0x_, *Importer22_, Xpetra::INSERT);
      D0x_.swap(D0xTmp_);
    }
  }


  template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void RefMaxwell<Scalar,LocalOrdinal,GlobalOrdinal,Node>::applyInverseAdditive(const MultiVector& RHS, MultiVector& X) const {

    Scalar one = Teuchos::ScalarTraits<Scalar>::one(), negone = -one, zero = Teuchos::ScalarTraits<Scalar>::zero();

    { // compute residuals

      Teuchos::TimeMonitor tmRes(*Teuchos::TimeMonitor::getNewTimer("MueLu RefMaxwell: residual calculation"));
      SM_Matrix_->apply(X, *residual_, Teuchos::NO_TRANS, one, zero);
      residual_->update(one, RHS, negone);
      R11_->apply(*residual_,*P11res_,Teuchos::NO_TRANS);
      D0_T_Matrix_->apply(*residual_,*D0res_,Teuchos::NO_TRANS);
    }

    // block diagonal preconditioner on 2x2 (V-cycle for diagonal blocks)

    if (!ImporterH_.is_null()) {
      Teuchos::TimeMonitor tmH(*Teuchos::TimeMonitor::getNewTimer("MueLu RefMaxwell: import coarse (1,1)"));
      P11resTmp_->doImport(*P11res_, *ImporterH_, Xpetra::INSERT);
      P11res_.swap(P11resTmp_);
    }
    if (!Importer22_.is_null()) {
      Teuchos::TimeMonitor tm22(*Teuchos::TimeMonitor::getNewTimer("MueLu RefMaxwell: import (2,2)"));
      D0resTmp_->doImport(*D0res_, *Importer22_, Xpetra::INSERT);
      D0res_.swap(D0resTmp_);
    }

    // iterate on coarse (1, 1) block
    if (!AH_.is_null()) {
      Teuchos::TimeMonitor tmH(*Teuchos::TimeMonitor::getNewTimer("MueLu RefMaxwell: solve coarse (1,1)"));

      RCP<const Map> origXMap = P11x_->getMap();
      RCP<const Map> origRhsMap = P11res_->getMap();

      // Replace maps with maps with a subcommunicator
      P11res_->replaceMap(AH_->getRangeMap());
      P11x_  ->replaceMap(AH_->getDomainMap());
      HierarchyH_->Iterate(*P11res_, *P11x_, 1, true);
      P11x_  ->replaceMap(origXMap);
      P11res_->replaceMap(origRhsMap);
    }

    // iterate on (2, 2) block
    if (!A22_.is_null()) {
      Teuchos::TimeMonitor tm22(*Teuchos::TimeMonitor::getNewTimer("MueLu RefMaxwell: solve (2,2)"));

      RCP<const Map> origXMap = D0x_->getMap();
      RCP<const Map> origRhsMap = D0res_->getMap();

      // Replace maps with maps with a subcommunicator
      D0res_->replaceMap(A22_->getRangeMap());
      D0x_  ->replaceMap(A22_->getDomainMap());
      Hierarchy22_->Iterate(*D0res_, *D0x_, 1, true);
      D0x_  ->replaceMap(origXMap);
      D0res_->replaceMap(origRhsMap);
    }

    if (!Importer22_.is_null()) {
      Teuchos::TimeMonitor tm22(*Teuchos::TimeMonitor::getNewTimer("MueLu RefMaxwell: export (2,2)"));
      D0xTmp_->doExport(*D0x_, *Importer22_, Xpetra::INSERT);
      D0x_.swap(D0xTmp_);
    }

    if (!ImporterH_.is_null()) {
      Teuchos::TimeMonitor tmH(*Teuchos::TimeMonitor::getNewTimer("MueLu RefMaxwell: export coarse (1,1)"));
      P11xTmp_->doExport(*P11x_, *ImporterH_, Xpetra::INSERT);
      P11x_.swap(P11xTmp_);
    }

    { // update current solution
      Teuchos::TimeMonitor tmUp(*Teuchos::TimeMonitor::getNewTimer("MueLu RefMaxwell: update"));
      P11_->apply(*P11x_,*residual_,Teuchos::NO_TRANS);
      D0_Matrix_->apply(*D0x_,*residual_,Teuchos::NO_TRANS,one,one);
      X.update(one, *residual_, one);

      if (!ImporterH_.is_null()) {
        P11res_.swap(P11resTmp_);
        P11x_.swap(P11xTmp_);
      }
      if (!Importer22_.is_null()) {
        D0res_.swap(D0resTmp_);
        D0x_.swap(D0xTmp_);
      }

    }
  }


  template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void RefMaxwell<Scalar,LocalOrdinal,GlobalOrdinal,Node>::applyInverse121(const MultiVector& RHS, MultiVector& X) const {

    // precondition (1,1)-block
    Scalar one = Teuchos::ScalarTraits<Scalar>::one(), negone = -one, zero = Teuchos::ScalarTraits<Scalar>::zero();
    SM_Matrix_->apply(X, *residual_, Teuchos::NO_TRANS, one, zero);
    residual_->update(one, RHS, negone);
    R11_->apply(*residual_,*P11res_,Teuchos::NO_TRANS);
    solveH();
    P11_->apply(*P11x_,*residual_,Teuchos::NO_TRANS);
    X.update(one, *residual_, one);
    if (!ImporterH_.is_null()) {
      P11res_.swap(P11resTmp_);
      P11x_.swap(P11xTmp_);
    }

    // precondition (2,2)-block
    SM_Matrix_->apply(X, *residual_, Teuchos::NO_TRANS, one, zero);
    residual_->update(one, RHS, negone);
    D0_T_Matrix_->apply(*residual_,*D0res_,Teuchos::NO_TRANS);
    solve22();
    D0_Matrix_->apply(*D0x_,*residual_,Teuchos::NO_TRANS);
    X.update(one, *residual_, one);
    if (!Importer22_.is_null()) {
      D0res_.swap(D0resTmp_);
      D0x_.swap(D0xTmp_);
    }

    // precondition (1,1)-block
    SM_Matrix_->apply(X, *residual_, Teuchos::NO_TRANS, one, zero);
    residual_->update(one, RHS, negone);
    R11_->apply(*residual_,*P11res_,Teuchos::NO_TRANS);
    solveH();
    P11_->apply(*P11x_,*residual_,Teuchos::NO_TRANS);
    X.update(one, *residual_, one);
    if (!ImporterH_.is_null()) {
      P11res_.swap(P11resTmp_);
      P11x_.swap(P11xTmp_);
    }

  }


  template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void RefMaxwell<Scalar,LocalOrdinal,GlobalOrdinal,Node>::applyInverse212(const MultiVector& RHS, MultiVector& X) const {

    // precondition (2,2)-block
    Scalar one = Teuchos::ScalarTraits<Scalar>::one(), negone = -one, zero = Teuchos::ScalarTraits<Scalar>::zero();
    SM_Matrix_->apply(X, *residual_, Teuchos::NO_TRANS, one, zero);
    residual_->update(one, RHS, negone);
    D0_T_Matrix_->apply(*residual_,*D0res_,Teuchos::NO_TRANS);
    solve22();
    D0_Matrix_->apply(*D0x_,*residual_,Teuchos::NO_TRANS);
    X.update(one, *residual_, one);
    if (!Importer22_.is_null()) {
      D0res_.swap(D0resTmp_);
      D0x_.swap(D0xTmp_);
    }

    // precondition (1,1)-block
    SM_Matrix_->apply(X, *residual_, Teuchos::NO_TRANS, one, zero);
    residual_->update(one, RHS, negone);
    R11_->apply(*residual_,*P11res_,Teuchos::NO_TRANS);
    solveH();
    P11_->apply(*P11x_,*residual_,Teuchos::NO_TRANS);
    X.update(one, *residual_, one);
    if (!ImporterH_.is_null()) {
      P11res_.swap(P11resTmp_);
      P11x_.swap(P11xTmp_);
    }

    // precondition (2,2)-block
    SM_Matrix_->apply(X, *residual_, Teuchos::NO_TRANS, one, zero);
    residual_->update(one, RHS, negone);
    D0_T_Matrix_->apply(*residual_,*D0res_,Teuchos::NO_TRANS);
    solve22();
    D0_Matrix_->apply(*D0x_,*residual_,Teuchos::NO_TRANS);
    X.update(one, *residual_, one);
    if (!Importer22_.is_null()) {
      D0res_.swap(D0resTmp_);
      D0x_.swap(D0xTmp_);
    }

  }

  template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void RefMaxwell<Scalar,LocalOrdinal,GlobalOrdinal,Node>::applyInverse11only(const MultiVector& RHS, MultiVector& X) const {

    // compute residuals
    Scalar one = Teuchos::ScalarTraits<Scalar>::one(), negone = -one, zero = Teuchos::ScalarTraits<Scalar>::zero();
    SM_Matrix_->apply(X, *residual_, Teuchos::NO_TRANS, one, zero);
    residual_->update(one, RHS, negone);
    R11_->apply(*residual_,*P11res_,Teuchos::NO_TRANS);

    solveH();

    // update current solution
    P11_->apply(*P11x_,*residual_,Teuchos::NO_TRANS);
    X.update(one, *residual_, one);

    if (!ImporterH_.is_null()) {
      P11res_.swap(P11resTmp_);
      P11x_.swap(P11xTmp_);
    }

  }


  template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void RefMaxwell<Scalar,LocalOrdinal,GlobalOrdinal,Node>::apply (const MultiVector& RHS, MultiVector& X,
                                                                  Teuchos::ETransp mode,
                                                                  Scalar alpha,
                                                                  Scalar beta) const {

    Teuchos::TimeMonitor tm(*Teuchos::TimeMonitor::getNewTimer("MueLu RefMaxwell: solve"));

    // make sure that we have enough temporary memory
    if (X.getNumVectors() != P11res_->getNumVectors()) {
      P11res_    = MultiVectorFactory::Build(R11_->getRangeMap(), X.getNumVectors());
      if (!ImporterH_.is_null()) {
        P11resTmp_ = MultiVectorFactory::Build(ImporterH_->getTargetMap(), X.getNumVectors());
        P11xTmp_   = MultiVectorFactory::Build(ImporterH_->getSourceMap(), X.getNumVectors());
        P11x_      = MultiVectorFactory::Build(ImporterH_->getTargetMap(), X.getNumVectors());
      } else
        P11x_      = MultiVectorFactory::Build(P11_->getDomainMap(), X.getNumVectors());
      if (!Importer22_.is_null()) {
        D0resTmp_ = MultiVectorFactory::Build(Importer22_->getTargetMap(), X.getNumVectors());
        D0xTmp_   = MultiVectorFactory::Build(Importer22_->getSourceMap(), X.getNumVectors());
        D0x_      = MultiVectorFactory::Build(Importer22_->getTargetMap(), X.getNumVectors());
      } else
        D0x_      = MultiVectorFactory::Build(D0_Matrix_->getDomainMap(), X.getNumVectors());
      residual_  = MultiVectorFactory::Build(SM_Matrix_->getDomainMap(), X.getNumVectors());

    }

    { // apply pre-smoothing

      Teuchos::TimeMonitor tmSm(*Teuchos::TimeMonitor::getNewTimer("MueLu RefMaxwell: smoothing"));

#if defined(HAVE_MUELU_IFPACK2) && (!defined(HAVE_MUELU_EPETRA) || defined(HAVE_MUELU_INST_DOUBLE_INT_INT))
      if (useHiptmairSmoothing_) {
        Tpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> tX = Utilities::MV2NonConstTpetraMV(X);
        Tpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> tRHS = Utilities::MV2TpetraMV(RHS);
        hiptmairPreSmoother_->apply(tRHS, tX);
      }
      else
#endif
        PreSmoother_->Apply(X, RHS, use_as_preconditioner_);
    }

    // do solve for the 2x2 block system
    if(mode_=="additive")
      applyInverseAdditive(RHS,X);
    else if(mode_=="121")
      applyInverse121(RHS,X);
    else if(mode_=="212")
      applyInverse212(RHS,X);
    else if(mode_=="11only")
      applyInverse11only(RHS,X);
    else if(mode_=="none") {
      // do nothing
    }
    else
      applyInverseAdditive(RHS,X);

    { // apply post-smoothing

      Teuchos::TimeMonitor tmSm(*Teuchos::TimeMonitor::getNewTimer("MueLu RefMaxwell: smoothing"));

#if defined(HAVE_MUELU_IFPACK2) && (!defined(HAVE_MUELU_EPETRA) || defined(HAVE_MUELU_INST_DOUBLE_INT_INT))
      if (useHiptmairSmoothing_)
        {
          Tpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> tX = Utilities::MV2NonConstTpetraMV(X);
          Tpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> tRHS = Utilities::MV2TpetraMV(RHS);
          hiptmairPostSmoother_->apply(tRHS, tX);
        }
      else
#endif
        PostSmoother_->Apply(X, RHS, false);
    }

  }


  template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  bool RefMaxwell<Scalar,LocalOrdinal,GlobalOrdinal,Node>::hasTransposeApply() const {
    return false;
  }


  template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void RefMaxwell<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
  initialize(const Teuchos::RCP<Matrix> & D0_Matrix,
             const Teuchos::RCP<Matrix> & M0inv_Matrix,
             const Teuchos::RCP<Matrix> & M1_Matrix,
             const Teuchos::RCP<MultiVector>  & Nullspace,
             const Teuchos::RCP<RealValuedMultiVector>  & Coords,
             Teuchos::ParameterList& List)
  {
    // some pre-conditions
    TEUCHOS_ASSERT(D0_Matrix!=Teuchos::null);
    TEUCHOS_ASSERT(M1_Matrix!=Teuchos::null);

    HierarchyH_ = Teuchos::null;
    Hierarchy22_ = Teuchos::null;
    PreSmoother_ = Teuchos::null;
    PostSmoother_ = Teuchos::null;
    disable_addon_ = false;
    mode_ = "additive";

    // set parameters
    setParameters(List);

    D0_Matrix_ = D0_Matrix;
    M0inv_Matrix_ = M0inv_Matrix;
    M1_Matrix_ = M1_Matrix;
    Coords_ = Coords;
    Nullspace_ = Nullspace;
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void RefMaxwell<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
  describe(Teuchos::FancyOStream& out, const Teuchos::EVerbosityLevel verbLevel) const {

    std::ostringstream oss;

    RCP<const Teuchos::Comm<int> > comm = SM_Matrix_->getDomainMap()->getComm();

    int root;
    if (!A22_.is_null())
      root = comm->getRank();
    else
      root = -1;

#ifdef HAVE_MPI
    int actualRoot;
    reduceAll(*comm, Teuchos::REDUCE_MAX, root, Teuchos::ptr(&actualRoot));
    root = actualRoot;
#endif


    oss << "\n--------------------------------------------------------------------------------\n" <<
      "---                            RefMaxwell Summary                            ---\n"
      "--------------------------------------------------------------------------------" << std::endl;
    oss << std::endl;

    GlobalOrdinal numRows;
    GlobalOrdinal nnz;

    SM_Matrix_->getRowMap()->getComm()->barrier();

    numRows = SM_Matrix_->getGlobalNumRows();
    nnz = SM_Matrix_->getGlobalNumEntries();

    Xpetra::global_size_t tt = numRows;
    int rowspacer = 3; while (tt != 0) { tt /= 10; rowspacer++; }
    tt = nnz;
    int nnzspacer = 2; while (tt != 0) { tt /= 10; nnzspacer++; }

    oss  << "block " << std::setw(rowspacer) << " rows " << std::setw(nnzspacer) << " nnz " << std::setw(9) << " nnz/row" << std::endl;
    oss << "(1, 1)" << std::setw(rowspacer) << numRows << std::setw(nnzspacer) << nnz << std::setw(9) << as<double>(nnz) / numRows << std::endl;

    if (!A22_.is_null()) {
      // ToDo: make sure that this is printed correctly
      numRows = A22_->getGlobalNumRows();
      nnz = A22_->getGlobalNumEntries();

      oss << "(2, 2)" << std::setw(rowspacer) << numRows << std::setw(nnzspacer) << nnz << std::setw(9) << as<double>(nnz) / numRows << std::endl;
    }

    oss << std::endl;

    std::string outstr = oss.str();

#ifdef HAVE_MPI
    RCP<const Teuchos::MpiComm<int> > mpiComm = rcp_dynamic_cast<const Teuchos::MpiComm<int> >(comm);
    MPI_Comm rawComm = (*mpiComm->getRawMpiComm())();

    int strLength = outstr.size();
    MPI_Bcast(&strLength, 1, MPI_INT, root, rawComm);
    if (comm->getRank() != root)
      outstr.resize(strLength);
    MPI_Bcast(&outstr[0], strLength, MPI_CHAR, root, rawComm);
#endif

    out << outstr;

    if (IsPrint(Statistics2)) {
      // Print the grid of processors
      std::ostringstream oss2;

      oss2 << "Sub-solver distribution over ranks" << std::endl;
      oss2 << "( (1,1) block only is indicated by '1', (2,2) block only by '2', and both blocks by 'B' and none by '.')" << std::endl;

      int numProcs = comm->getSize();
#ifdef HAVE_MPI
      RCP<const Teuchos::MpiComm<int> > tmpic = rcp_dynamic_cast<const Teuchos::MpiComm<int> >(comm);
      TEUCHOS_TEST_FOR_EXCEPTION(tmpic == Teuchos::null, Exceptions::RuntimeError, "Cannot cast base Teuchos::Comm to Teuchos::MpiComm object.");
      RCP<const Teuchos::OpaqueWrapper<MPI_Comm> > rawMpiComm = tmpic->getRawMpiComm();
#endif

      char status = 0;
      if (!AH_.is_null())
        status += 1;
      if (!A22_.is_null())
        status += 2;
      std::vector<char> states(numProcs, 0);
#ifdef HAVE_MPI
      MPI_Gather(&status, 1, MPI_CHAR, &states[0], 1, MPI_CHAR, 0, *rawMpiComm);
#else
      states.push_back(status);
#endif

      int rowWidth = std::min(Teuchos::as<int>(ceil(sqrt(numProcs))), 100);
      for (int proc = 0; proc < numProcs; proc += rowWidth) {
        for (int j = 0; j < rowWidth; j++)
          if (proc + j < numProcs)
            if (states[proc+j] == 0)
              oss2 << ".";
            else if (states[proc+j] == 1)
              oss2 << "1";
            else if (states[proc+j] == 2)
              oss2 << "2";
            else
              oss2 << "B";
          else
            oss2 << " ";

        oss2 << "      " << proc << ":" << std::min(proc + rowWidth, numProcs) - 1 << std::endl;
      }
      oss2 << std::endl;
      GetOStream(Statistics2) << oss2.str();
    }


  }

} // namespace

#define MUELU_REFMAXWELL_SHORT
#endif //ifdef MUELU_REFMAXWELL_DEF_HPP
