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
#include "MueLu_SaPFactory.hpp"
#include "MueLu_AggregationExportFactory.hpp"
#include "MueLu_Utilities.hpp"
#include "MueLu_Maxwell_Utils.hpp"

#include "MueLu_AmalgamationFactory_kokkos.hpp"
#include "MueLu_CoalesceDropFactory_kokkos.hpp"
#include "MueLu_CoarseMapFactory_kokkos.hpp"
#include "MueLu_CoordinatesTransferFactory_kokkos.hpp"
#include "MueLu_UncoupledAggregationFactory_kokkos.hpp"
#include "MueLu_TentativePFactory_kokkos.hpp"
#include "MueLu_SaPFactory_kokkos.hpp"
#include "MueLu_Utilities_kokkos.hpp"
#include <Kokkos_Core.hpp>
#include <KokkosSparse_CrsMatrix.hpp>

#include "MueLu_ZoltanInterface.hpp"
#include "MueLu_Zoltan2Interface.hpp"
#include "MueLu_RepartitionHeuristicFactory.hpp"
#include "MueLu_RepartitionFactory.hpp"
#include "MueLu_RebalanceAcFactory.hpp"
#include "MueLu_RebalanceTransferFactory.hpp"

#include "MueLu_VerbosityLevel.hpp"

#include <MueLu_CreateXpetraPreconditioner.hpp>
#include <MueLu_ML2MueLuParameterTranslator.hpp>

#ifdef HAVE_MUELU_CUDA
#include "cuda_profiler_api.h"
#endif

// Stratimikos
#if defined(HAVE_MUELU_STRATIMIKOS) && defined(HAVE_MUELU_THYRA)
#include <Xpetra_ThyraLinearOp.hpp>
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

    if (list.isType<std::string>("parameterlist: syntax") && list.get<std::string>("parameterlist: syntax") == "ml") {
      Teuchos::ParameterList newList = *Teuchos::getParametersFromXmlString(MueLu::ML2MueLuParameterTranslator::translate(list,"refmaxwell"));
      if(list.isSublist("refmaxwell: 11list") && list.sublist("refmaxwell: 11list").isSublist("edge matrix free: coarse"))
        newList.sublist("refmaxwell: 11list") = *Teuchos::getParametersFromXmlString(MueLu::ML2MueLuParameterTranslator::translate(list.sublist("refmaxwell: 11list").sublist("edge matrix free: coarse"),"SA"));
      if(list.isSublist("refmaxwell: 22list"))
        newList.sublist("refmaxwell: 22list") = *Teuchos::getParametersFromXmlString(MueLu::ML2MueLuParameterTranslator::translate(list.sublist("refmaxwell: 22list"),"SA"));
      list = newList;
    }

    parameterList_             = list;
    std::string verbosityLevel = parameterList_.get<std::string>("verbosity", MasterList::getDefault<std::string>("verbosity"));
    VerboseObject::SetDefaultVerbLevel(toVerbLevel(verbosityLevel));
    std::string outputFilename = parameterList_.get<std::string>("output filename", MasterList::getDefault<std::string>("output filename"));
    if (outputFilename != "")
      VerboseObject::SetMueLuOFileStream(outputFilename);
    if (parameterList_.isType<Teuchos::RCP<Teuchos::FancyOStream> >("output stream"))
      VerboseObject::SetMueLuOStream(parameterList_.get<Teuchos::RCP<Teuchos::FancyOStream> >("output stream"));

    if (parameterList_.get("print initial parameters",MasterList::getDefault<bool>("print initial parameters")))
      GetOStream(static_cast<MsgType>(Runtime1), 0) << parameterList_ << std::endl;
    disable_addon_             = list.get("refmaxwell: disable addon",         MasterList::getDefault<bool>("refmaxwell: disable addon"));
    mode_                      = list.get("refmaxwell: mode",                  MasterList::getDefault<std::string>("refmaxwell: mode"));
    use_as_preconditioner_     = list.get("refmaxwell: use as preconditioner", MasterList::getDefault<bool>("refmaxwell: use as preconditioner"));
    dump_matrices_             = list.get("refmaxwell: dump matrices",         MasterList::getDefault<bool>("refmaxwell: dump matrices"));
    enable_reuse_              = list.get("refmaxwell: enable reuse",          MasterList::getDefault<bool>("refmaxwell: enable reuse"));
    implicitTranspose_         = list.get("transpose: use implicit",           MasterList::getDefault<bool>("transpose: use implicit"));
    fuseProlongationAndUpdate_ = list.get("fuse prolongation and update",      MasterList::getDefault<bool>("fuse prolongation and update"));
    skipFirstLevel_            = list.get("refmaxwell: skip first (1,1) level", MasterList::getDefault<bool>("refmaxwell: skip first (1,1) level"));
    syncTimers_                = list.get("sync timers",                       false);
    numItersH_                 = list.get("refmaxwell: num iters H",           1);
    numIters22_                = list.get("refmaxwell: num iters 22",          1);
    applyBCsToAnodal_          = list.get("refmaxwell: apply BCs to Anodal",   false);
    applyBCsToH_               = list.get("refmaxwell: apply BCs to H",        true);
    applyBCsTo22_              = list.get("refmaxwell: apply BCs to 22",       true);

    if(list.isSublist("refmaxwell: 11list"))
      precList11_     =  list.sublist("refmaxwell: 11list");
    if(!precList11_.isType<std::string>("Preconditioner Type") &&
       !precList11_.isType<std::string>("smoother: type") &&
       !precList11_.isType<std::string>("smoother: pre type") &&
       !precList11_.isType<std::string>("smoother: post type")) {
      precList11_.set("smoother: type", "CHEBYSHEV");
      precList11_.sublist("smoother: params").set("chebyshev: degree",2);
      precList11_.sublist("smoother: params").set("chebyshev: ratio eigenvalue",5.4);
      precList11_.sublist("smoother: params").set("chebyshev: eigenvalue max iterations",30);
    }

    if(list.isSublist("refmaxwell: 22list"))
      precList22_     =  list.sublist("refmaxwell: 22list");
    if(!precList22_.isType<std::string>("Preconditioner Type") &&
       !precList22_.isType<std::string>("smoother: type") &&
       !precList22_.isType<std::string>("smoother: pre type") &&
       !precList22_.isType<std::string>("smoother: post type")) {
      precList22_.set("smoother: type", "CHEBYSHEV");
      precList22_.sublist("smoother: params").set("chebyshev: degree",2);
      precList22_.sublist("smoother: params").set("chebyshev: ratio eigenvalue",7.0);
      precList22_.sublist("smoother: params").set("chebyshev: eigenvalue max iterations",30);
    }

    if(!list.isType<std::string>("smoother: type") && !list.isType<std::string>("smoother: pre type") && !list.isType<std::string>("smoother: post type")) {
      list.set("smoother: type", "CHEBYSHEV");
      list.sublist("smoother: params").set("chebyshev: degree",2);
      list.sublist("smoother: params").set("chebyshev: ratio eigenvalue",20.0);
      list.sublist("smoother: params").set("chebyshev: eigenvalue max iterations",30);
    }

    if(list.isSublist("smoother: params")) {
      smootherList_ = list.sublist("smoother: params");
    }

    if (enable_reuse_ &&
        !precList11_.isType<std::string>("Preconditioner Type") &&
        !precList11_.isParameter("reuse: type"))
      precList11_.set("reuse: type", "full");
    if (enable_reuse_ &&
        !precList22_.isType<std::string>("Preconditioner Type") &&
        !precList22_.isParameter("reuse: type"))
      precList22_.set("reuse: type", "full");

# ifdef HAVE_MUELU_SERIAL
    if (typeid(Node).name() == typeid(Kokkos::Compat::KokkosSerialWrapperNode).name())
      useKokkos_ = false;
# endif
# ifdef HAVE_MUELU_OPENMP
    if (typeid(Node).name() == typeid(Kokkos::Compat::KokkosOpenMPWrapperNode).name())
      useKokkos_ = true;
# endif
# ifdef HAVE_MUELU_CUDA
    if (typeid(Node).name() == typeid(Kokkos::Compat::KokkosCudaWrapperNode).name())
      useKokkos_ = true;
# endif
# ifdef HAVE_MUELU_HIP
    if (typeid(Node).name() == typeid(Kokkos::Compat::KokkosHIPWrapperNode).name())
      useKokkos_ = true;
# endif
    useKokkos_ = list.get("use kokkos refactor",useKokkos_);
  }



  template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void RefMaxwell<Scalar,LocalOrdinal,GlobalOrdinal,Node>::compute(bool reuse) {

#ifdef HAVE_MUELU_CUDA
    if (parameterList_.get<bool>("refmaxwell: cuda profile setup", false)) cudaProfilerStart();
#endif

    std::string timerLabel;
    if (reuse)
      timerLabel = "MueLu RefMaxwell: compute (reuse)";
    else
      timerLabel = "MueLu RefMaxwell: compute";
    RCP<Teuchos::TimeMonitor> tmCompute = getTimer(timerLabel);

    ////////////////////////////////////////////////////////////////////////////////
    // Remove explicit zeros from matrices
    Maxwell_Utils<SC,LO,GO,NO>::removeExplicitZeros(parameterList_,D0_Matrix_,SM_Matrix_,M1_Matrix_,Ms_Matrix_);

    if (IsPrint(Statistics2)) {
      RCP<ParameterList> params = rcp(new ParameterList());;
      params->set("printLoadBalancingInfo", true);
      params->set("printCommInfo",          true);
      GetOStream(Statistics2) << PerfUtils::PrintMatrixInfo(*SM_Matrix_, "SM_Matrix", params);
    }

    ////////////////////////////////////////////////////////////////////////////////
    // Detect Dirichlet boundary conditions
    if (!reuse) {
      magnitudeType rowSumTol = parameterList_.get("refmaxwell: row sum drop tol (1,1)",-1.0);
      Maxwell_Utils<SC,LO,GO,NO>::detectBoundaryConditionsSM(SM_Matrix_,D0_Matrix_,rowSumTol,
                                                             useKokkos_,BCrowsKokkos_,BCcolsKokkos_,BCdomainKokkos_,
                                                             BCedges_,BCnodes_,BCrows_,BCcols_,BCdomain_,
                                                             allEdgesBoundary_,allNodesBoundary_);
      if (IsPrint(Statistics2)) {
        GetOStream(Statistics2) << "MueLu::RefMaxwell::compute(): Detected " << BCedges_ << " BC rows and " << BCnodes_ << " BC columns." << std::endl;
      }
    }

    if (allEdgesBoundary_) {
      // All edges have been detected as boundary edges.
      // Do not attempt to construct sub-hierarchies, but just set up a single level preconditioner.
      GetOStream(Warnings0) << "All edges are detected as boundary edges!" << std::endl;
      mode_ = "none";
      setFineLevelSmoother();
      return;
    }

    ////////////////////////////////////////////////////////////////////////////////
    // build nullspace if necessary

    if(Nullspace_ != null) {
      // no need to do anything - nullspace is built
      TEUCHOS_ASSERT(Nullspace_->getMap()->isCompatible(*(SM_Matrix_->getRowMap())));
    }
    else if(Nullspace_ == null && Coords_ != null) {
      RCP<MultiVector> CoordsSC;
      if (useKokkos_)
        CoordsSC = Utilities_kokkos::RealValuedToScalarMultiVector(Coords_);
      else
        CoordsSC = Utilities::RealValuedToScalarMultiVector(Coords_);
      Nullspace_ = MultiVectorFactory::Build(SM_Matrix_->getRowMap(),Coords_->getNumVectors());
      D0_Matrix_->apply(*CoordsSC,*Nullspace_);

      bool normalize = parameterList_.get<bool>("refmaxwell: normalize nullspace", MasterList::getDefault<bool>("refmaxwell: normalize nullspace"));
      
      coordinateType minLen, maxLen, meanLen;
      if (IsPrint(Statistics2) || normalize){
        // compute edge lengths
        ArrayRCP<ArrayRCP<const Scalar> > localNullspace(Nullspace_->getNumVectors());
        for (size_t i = 0; i < Nullspace_->getNumVectors(); i++)
          localNullspace[i] = Nullspace_->getData(i);
        coordinateType localMinLen = Teuchos::ScalarTraits<coordinateType>::rmax();
        coordinateType localMeanLen = Teuchos::ScalarTraits<coordinateType>::zero();
        coordinateType localMaxLen = Teuchos::ScalarTraits<coordinateType>::zero();
        for (size_t j=0; j < Nullspace_->getMap()->getLocalNumElements(); j++) {
          Scalar lenSC = Teuchos::ScalarTraits<Scalar>::zero();
          for (size_t i=0; i < Nullspace_->getNumVectors(); i++)
           lenSC += localNullspace[i][j]*localNullspace[i][j];
          coordinateType len = Teuchos::as<coordinateType>(Teuchos::ScalarTraits<Scalar>::real(Teuchos::ScalarTraits<Scalar>::squareroot(lenSC)));
          localMinLen = std::min(localMinLen, len);
          localMaxLen = std::max(localMaxLen, len);
          localMeanLen += len;
        }
        
        RCP<const Teuchos::Comm<int> > comm = Nullspace_->getMap()->getComm();
        MueLu_minAll(comm, localMinLen,  minLen);
        MueLu_sumAll(comm, localMeanLen, meanLen);
        MueLu_maxAll(comm, localMaxLen,  maxLen);
        meanLen /= Nullspace_->getMap()->getGlobalNumElements();
      }
      
      if (IsPrint(Statistics2)) {
        GetOStream(Statistics0) << "Edge length (min/mean/max): " << minLen << " / " << meanLen << " / " << maxLen << std::endl;
      }
      
      if (normalize) {
        // normalize the nullspace
        GetOStream(Runtime0) << "RefMaxwell::compute(): normalizing nullspace" << std::endl;

        const Scalar one = Teuchos::ScalarTraits<Scalar>::one();

        Array<Scalar> normsSC(Coords_->getNumVectors(), one / Teuchos::as<Scalar>(meanLen));
        Nullspace_->scale(normsSC());
      }
    }
    else {
      GetOStream(Errors) << "MueLu::RefMaxwell::compute(): either the nullspace or the nodal coordinates must be provided." << std::endl;
    }

    if (!reuse && skipFirstLevel_) {
      // Nuke the BC edges in nullspace
      if (useKokkos_)
        Utilities_kokkos::ZeroDirichletRows(Nullspace_,BCrowsKokkos_);
      else
        Utilities::ZeroDirichletRows(Nullspace_,BCrows_);
      dump(*Nullspace_, "nullspace.m");
    }

    ////////////////////////////////////////////////////////////////////////////////
    // build special prolongator for (1,1)-block

    if(P11_.is_null()) {
      if (skipFirstLevel_) {
        // Form A_nodal = D0* Ms D0  (aka TMT_agg)
        std::string label("D0*Ms*D0");
        A_nodal_Matrix_ = Maxwell_Utils<SC,LO,GO,NO>::PtAPWrapper(Ms_Matrix_,D0_Matrix_,parameterList_,label);

        if (applyBCsToAnodal_) {
          // Apply boundary conditions to A_nodal
          if (useKokkos_)
            Utilities_kokkos::ApplyOAZToMatrixRows(A_nodal_Matrix_,BCdomainKokkos_);
          else
            Utilities::ApplyOAZToMatrixRows(A_nodal_Matrix_,BCdomain_);
        }
        dump(*A_nodal_Matrix_, "A_nodal.m");
      }

      // build special prolongator
      GetOStream(Runtime0) << "RefMaxwell::compute(): building special prolongator" << std::endl;
      buildProlongator();
    }

#ifdef HAVE_MPI
    bool doRebalancing = false;
    doRebalancing = parameterList_.get<bool>("refmaxwell: subsolves on subcommunicators", MasterList::getDefault<bool>("refmaxwell: subsolves on subcommunicators"));
    int rebalanceStriding = parameterList_.get<int>("refmaxwell: subsolves striding", -1);
    int numProcsAH, numProcsA22;
    int numProcs = SM_Matrix_->getDomainMap()->getComm()->getSize();
    if (numProcs == 1)
      doRebalancing = false;
#endif
    {
      // build coarse grid operator for (1,1)-block
      formCoarseMatrix();

      if (!reuse) {
#ifdef HAVE_MPI
        if (doRebalancing) {

          {
            // decide on number of ranks for coarse (1, 1) problem

            Level level;
            level.SetFactoryManager(null);
            level.SetLevelID(0);
            level.Set("A",AH_);

            auto repartheurFactory = rcp(new RepartitionHeuristicFactory());
            ParameterList repartheurParams;
            repartheurParams.set("repartition: start level",            0);
            // Setting min == target on purpose.
            int defaultTargetRows = 10000;
            repartheurParams.set("repartition: min rows per proc",      precList11_.get<int>("repartition: target rows per proc", defaultTargetRows));
            repartheurParams.set("repartition: target rows per proc",   precList11_.get<int>("repartition: target rows per proc", defaultTargetRows));
            repartheurParams.set("repartition: min rows per thread",    precList11_.get<int>("repartition: target rows per thread", defaultTargetRows));
            repartheurParams.set("repartition: target rows per thread", precList11_.get<int>("repartition: target rows per thread", defaultTargetRows));
            repartheurParams.set("repartition: max imbalance",          precList11_.get<double>("repartition: max imbalance", 1.1));
            repartheurFactory->SetParameterList(repartheurParams);

            level.Request("number of partitions", repartheurFactory.get());
            repartheurFactory->Build(level);
            numProcsAH = level.Get<int>("number of partitions", repartheurFactory.get());
            numProcsAH = std::min(numProcsAH,numProcs);
          }

          {
            // decide on number of ranks for (2, 2) problem

            Level level;
            level.SetFactoryManager(null);
            level.SetLevelID(0);

            level.Set("Map",D0_Matrix_->getDomainMap());

            auto repartheurFactory = rcp(new RepartitionHeuristicFactory());
            ParameterList repartheurParams;
            repartheurParams.set("repartition: start level",            0);
            repartheurParams.set("repartition: use map",                true);
            // Setting min == target on purpose.
            int defaultTargetRows = 10000;
            repartheurParams.set("repartition: min rows per proc",      precList22_.get<int>("repartition: target rows per proc", defaultTargetRows));
            repartheurParams.set("repartition: target rows per proc",   precList22_.get<int>("repartition: target rows per proc", defaultTargetRows));
            repartheurParams.set("repartition: min rows per thread",    precList22_.get<int>("repartition: target rows per thread", defaultTargetRows));
            repartheurParams.set("repartition: target rows per thread", precList22_.get<int>("repartition: target rows per thread", defaultTargetRows));
            // repartheurParams.set("repartition: max imbalance",        precList22_.get<double>("repartition: max imbalance", 1.1));
            repartheurFactory->SetParameterList(repartheurParams);

            level.Request("number of partitions", repartheurFactory.get());
            repartheurFactory->Build(level);
            numProcsA22 = level.Get<int>("number of partitions", repartheurFactory.get());
            numProcsA22 = std::min(numProcsA22,numProcs);
          }

          if (rebalanceStriding >= 1) {
            TEUCHOS_ASSERT(rebalanceStriding*numProcsAH<=numProcs);
            TEUCHOS_ASSERT(rebalanceStriding*numProcsA22<=numProcs);
            if (rebalanceStriding*(numProcsAH+numProcsA22)>numProcs) {
              GetOStream(Warnings0) << "RefMaxwell::compute(): Disabling striding = " << rebalanceStriding << ", since AH needs " << numProcsAH
                                    << " procs and A22 needs " << numProcsA22 << " procs."<< std::endl;
              rebalanceStriding = -1;
            }
          }

          // }

          if ((numProcsAH < 0) || (numProcsA22 < 0) || (numProcsAH + numProcsA22 > numProcs)) {
            GetOStream(Warnings0) << "RefMaxwell::compute(): Disabling rebalancing of subsolves, since partition heuristic resulted "
                                  << "in undesirable number of partitions: " << numProcsAH << ", " << numProcsA22 << std::endl;
            doRebalancing = false;
          }

          // check again, as we could have changed the value above
          if (doRebalancing) {
            // rebalance AH
            RCP<Teuchos::TimeMonitor> tm = getTimer("MueLu RefMaxwell: Rebalance AH");

            Level fineLevel, coarseLevel;
            fineLevel.SetFactoryManager(null);
            coarseLevel.SetFactoryManager(null);
            coarseLevel.SetPreviousLevel(rcpFromRef(fineLevel));
            fineLevel.SetLevelID(0);
            coarseLevel.SetLevelID(1);
            coarseLevel.Set("A",AH_);
            coarseLevel.Set("P",P11_);
            coarseLevel.Set("Coordinates",CoordsH_);
            if (!NullspaceH_.is_null())
              coarseLevel.Set("Nullspace",NullspaceH_);
            coarseLevel.Set("number of partitions", numProcsAH);
            coarseLevel.Set("repartition: heuristic target rows per process", 1000);

            coarseLevel.setlib(AH_->getDomainMap()->lib());
            fineLevel.setlib(AH_->getDomainMap()->lib());
            coarseLevel.setObjectLabel("RefMaxwell coarse (1,1)");
            fineLevel.setObjectLabel("RefMaxwell coarse (1,1)");

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
            if (rebalanceStriding >= 1) {
              bool acceptPart = (SM_Matrix_->getDomainMap()->getComm()->getRank() % rebalanceStriding) == 0;
              if (SM_Matrix_->getDomainMap()->getComm()->getRank() >= numProcsAH*rebalanceStriding)
                acceptPart = false;
              repartParams.set("repartition: remap accept partition", acceptPart);
            }
            repartFactory->SetParameterList(repartParams);
            // repartFactory->SetFactory("number of partitions", repartheurFactory);
            repartFactory->SetFactory("Partition", partitioner);

            auto newP = rcp(new RebalanceTransferFactory());
            ParameterList newPparams;
            newPparams.set("type", "Interpolation");
            newPparams.set("repartition: rebalance P and R", precList11_.get<bool>("repartition: rebalance P and R", false));
            newPparams.set("repartition: use subcommunicators", true);
            newPparams.set("repartition: rebalance Nullspace", !NullspaceH_.is_null());
            newP->SetFactory("Coordinates", NoFactory::getRCP());
            if (!NullspaceH_.is_null())
              newP->SetFactory("Nullspace", NoFactory::getRCP());
            newP->SetParameterList(newPparams);
            newP->SetFactory("Importer", repartFactory);

            auto newA = rcp(new RebalanceAcFactory());
            ParameterList rebAcParams;
            rebAcParams.set("repartition: use subcommunicators", true);
            newA->SetParameterList(rebAcParams);
            newA->SetFactory("Importer", repartFactory);

            coarseLevel.Request("P", newP.get());
            coarseLevel.Request("Importer", repartFactory.get());
            coarseLevel.Request("A", newA.get());
            coarseLevel.Request("Coordinates", newP.get());
            if (!NullspaceH_.is_null())
              coarseLevel.Request("Nullspace", newP.get());
            repartFactory->Build(coarseLevel);

            if (!precList11_.get<bool>("repartition: rebalance P and R", false))
              ImporterH_ = coarseLevel.Get< RCP<const Import> >("Importer", repartFactory.get());
            P11_ = coarseLevel.Get< RCP<Matrix> >("P", newP.get());
            AH_ = coarseLevel.Get< RCP<Matrix> >("A", newA.get());
            CoordsH_ = coarseLevel.Get< RCP<RealValuedMultiVector> >("Coordinates", newP.get());
            if (!NullspaceH_.is_null())
              NullspaceH_ = coarseLevel.Get< RCP<MultiVector> >("Nullspace", newP.get());

            AH_AP_reuse_data_ = Teuchos::null;
            AH_RAP_reuse_data_ = Teuchos::null;

            if (!disable_addon_ && enable_reuse_) {
              // Rebalance the addon for next setup
              RCP<const Import> ImporterH = coarseLevel.Get< RCP<const Import> >("Importer", repartFactory.get());
              RCP<const Map> targetMap = ImporterH->getTargetMap();
              ParameterList XpetraList;
              XpetraList.set("Restrict Communicator",true);
              Addon_Matrix_ = MatrixFactory::Build(Addon_Matrix_, *ImporterH, *ImporterH, targetMap, targetMap, rcp(&XpetraList,false));
            }
          }
        }
#endif // HAVE_MPI

        // This should be taken out again as soon as
        // CoalesceDropFactory_kokkos supports BlockSize > 1 and
        // drop tol != 0.0
        if (useKokkos_ && precList11_.isParameter("aggregation: drop tol") && precList11_.get<double>("aggregation: drop tol") != 0.0) {
          GetOStream(Warnings0) << "RefMaxwell::compute(): Setting \"aggregation: drop tol\". to 0.0, since CoalesceDropFactory_kokkos does not "
                                << "support BlockSize > 1 and drop tol != 0.0" << std::endl;
          precList11_.set("aggregation: drop tol", 0.0);
        }
        dump(*P11_, "P11.m");

        if (!implicitTranspose_) {
          if (useKokkos_)
            R11_ = Utilities_kokkos::Transpose(*P11_);
          else
            R11_ = Utilities::Transpose(*P11_);
          dump(*R11_, "R11.m");
        }
      }

      VerbLevel verbosityLevel = VerboseObject::GetDefaultVerbLevel();
      if (!AH_.is_null()) {
        // set fixed block size for vector nodal matrix
        size_t dim = Nullspace_->getNumVectors();
        AH_->SetFixedBlockSize(dim);
        AH_->setObjectLabel("RefMaxwell coarse (1,1)");
        dump(*AH_, "AH.m");
        int oldRank = SetProcRankVerbose(AH_->getDomainMap()->getComm()->getRank());
        if (IsPrint(Statistics2)) {
          RCP<ParameterList> params = rcp(new ParameterList());;
          params->set("printLoadBalancingInfo", true);
          params->set("printCommInfo",          true);
          GetOStream(Statistics2) << PerfUtils::PrintMatrixInfo(*AH_, "AH", params);
        }
#if defined(HAVE_MUELU_STRATIMIKOS) && defined(HAVE_MUELU_THYRA)
        if (precList11_.isType<std::string>("Preconditioner Type")) {
          // build a Stratimikos preconditioner
          if (precList11_.get<std::string>("Preconditioner Type") == "MueLu") {
            ParameterList& userParamList = precList11_.sublist("Preconditioner Types").sublist("MueLu").sublist("user data");
            userParamList.set<RCP<RealValuedMultiVector> >("Coordinates", CoordsH_);
          }
          thyraPrecOpH_ = rcp(new XpetraThyraLinearOp<Scalar,LocalOrdinal,GlobalOrdinal,Node>(AH_, rcp(&precList11_, false)));
        } else
#endif
          {
          // build a MueLu hierarchy

          if (!reuse) {
            dumpCoords(*CoordsH_, "coordsH.m");
            if (!NullspaceH_.is_null())
              dump(*NullspaceH_, "NullspaceH.m");
            ParameterList& userParamList = precList11_.sublist("user data");
            userParamList.set<RCP<RealValuedMultiVector> >("Coordinates", CoordsH_);
            if (!NullspaceH_.is_null())
              userParamList.set<RCP<MultiVector> >("Nullspace", NullspaceH_);
            HierarchyH_ = MueLu::CreateXpetraPreconditioner(AH_, precList11_);
          } else {
            RCP<MueLu::Level> level0 = HierarchyH_->GetLevel(0);
            level0->Set("A", AH_);
            HierarchyH_->SetupRe();
          }
        }
        SetProcRankVerbose(oldRank);
      }
      VerboseObject::SetDefaultVerbLevel(verbosityLevel);

    }

    if(!reuse && applyBCsTo22_) {
      GetOStream(Runtime0) << "RefMaxwell::compute(): nuking BC nodes of D0" << std::endl;

      D0_Matrix_->resumeFill();
      Scalar replaceWith;
      if (D0_Matrix_->getRowMap()->lib() == Xpetra::UseEpetra)
        replaceWith= Teuchos::ScalarTraits<SC>::eps();
      else
        replaceWith = Teuchos::ScalarTraits<SC>::zero();
      if (useKokkos_) {
        Utilities_kokkos::ZeroDirichletCols(D0_Matrix_,BCcolsKokkos_,replaceWith);
      } else {
        Utilities::ZeroDirichletCols(D0_Matrix_,BCcols_,replaceWith);
      }
      D0_Matrix_->fillComplete(D0_Matrix_->getDomainMap(),D0_Matrix_->getRangeMap());
    }

    ////////////////////////////////////////////////////////////////////////////////
    // Build A22 = D0* SM D0 and hierarchy for A22
    if (!allNodesBoundary_) {
      GetOStream(Runtime0) << "RefMaxwell::compute(): building MG for (2,2)-block" << std::endl;

      if (!reuse) { // build fine grid operator for (2,2)-block, D0* SM D0  (aka TMT)
        RCP<Teuchos::TimeMonitor> tm = getTimer("MueLu RefMaxwell: Build A22");

        Level fineLevel, coarseLevel;
        fineLevel.SetFactoryManager(null);
        coarseLevel.SetFactoryManager(null);
        coarseLevel.SetPreviousLevel(rcpFromRef(fineLevel));
        fineLevel.SetLevelID(0);
        coarseLevel.SetLevelID(1);
        fineLevel.Set("A",SM_Matrix_);
        coarseLevel.Set("P",D0_Matrix_);
        coarseLevel.Set("Coordinates",Coords_);

        coarseLevel.setlib(SM_Matrix_->getDomainMap()->lib());
        fineLevel.setlib(SM_Matrix_->getDomainMap()->lib());
        coarseLevel.setObjectLabel("RefMaxwell (2,2)");
        fineLevel.setObjectLabel("RefMaxwell (2,2)");

        RCP<RAPFactory> rapFact = rcp(new RAPFactory());
        ParameterList rapList = *(rapFact->GetValidParameterList());
        rapList.set("transpose: use implicit", true);
        rapList.set("rap: fix zero diagonals", parameterList_.get<bool>("rap: fix zero diagonals", true));
        rapList.set("rap: fix zero diagonals threshold", parameterList_.get<double>("rap: fix zero diagonals threshold", Teuchos::ScalarTraits<double>::eps()));
        rapList.set("rap: triple product", parameterList_.get<bool>("rap: triple product", false));
        rapFact->SetParameterList(rapList);

        if (!A22_AP_reuse_data_.is_null()) {
          coarseLevel.AddKeepFlag("AP reuse data", rapFact.get());
          coarseLevel.Set<Teuchos::RCP<Teuchos::ParameterList> >("AP reuse data", A22_AP_reuse_data_, rapFact.get());
        }
        if (!A22_RAP_reuse_data_.is_null()) {
          coarseLevel.AddKeepFlag("RAP reuse data", rapFact.get());
          coarseLevel.Set<Teuchos::RCP<Teuchos::ParameterList> >("RAP reuse data", A22_RAP_reuse_data_, rapFact.get());
        }

#ifdef HAVE_MPI
        if (doRebalancing) {

          coarseLevel.Set("number of partitions", numProcsA22);
          coarseLevel.Set("repartition: heuristic target rows per process", 1000);

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
          if (rebalanceStriding >= 1) {
            bool acceptPart = ((SM_Matrix_->getDomainMap()->getComm()->getSize()-1-SM_Matrix_->getDomainMap()->getComm()->getRank()) % rebalanceStriding) == 0;
            if (SM_Matrix_->getDomainMap()->getComm()->getSize()-1-SM_Matrix_->getDomainMap()->getComm()->getRank() >= numProcsA22*rebalanceStriding)
              acceptPart = false;
            if (acceptPart)
              TEUCHOS_ASSERT(AH_.is_null());
            repartParams.set("repartition: remap accept partition", acceptPart);
          } else
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

          auto newA = rcp(new RebalanceAcFactory());
          ParameterList rebAcParams;
          rebAcParams.set("repartition: use subcommunicators", true);
          newA->SetParameterList(rebAcParams);
          newA->SetFactory("A", rapFact);
          newA->SetFactory("Importer", repartFactory);

          coarseLevel.Request("P", newP.get());
          coarseLevel.Request("Importer", repartFactory.get());
          coarseLevel.Request("A", newA.get());
          coarseLevel.Request("Coordinates", newP.get());
          rapFact->Build(fineLevel,coarseLevel);
          repartFactory->Build(coarseLevel);

          if (!precList22_.get<bool>("repartition: rebalance P and R", false))
            Importer22_ = coarseLevel.Get< RCP<const Import> >("Importer", repartFactory.get());
          D0_Matrix_ = coarseLevel.Get< RCP<Matrix> >("P", newP.get());
          A22_ = coarseLevel.Get< RCP<Matrix> >("A", newA.get());
          Coords_ = coarseLevel.Get< RCP<RealValuedMultiVector> >("Coordinates", newP.get());

        } else
#endif // HAVE_MPI
        {
          coarseLevel.Request("A", rapFact.get());
          if (enable_reuse_) {
            coarseLevel.Request("AP reuse data", rapFact.get());
            coarseLevel.Request("RAP reuse data", rapFact.get());
          }

          A22_ = coarseLevel.Get< RCP<Matrix> >("A", rapFact.get());

          if (enable_reuse_) {
            if (coarseLevel.IsAvailable("AP reuse data", rapFact.get()))
              A22_AP_reuse_data_ = coarseLevel.Get< RCP<ParameterList> >("AP reuse data", rapFact.get());
            if (coarseLevel.IsAvailable("RAP reuse data", rapFact.get()))
              A22_RAP_reuse_data_ = coarseLevel.Get< RCP<ParameterList> >("RAP reuse data", rapFact.get());
          }
        }
      } else {
        RCP<Teuchos::TimeMonitor> tm = getTimer("MueLu RefMaxwell: Build A22");
        if (Importer22_.is_null()) {
          RCP<Matrix> temp;
          temp = Xpetra::MatrixMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>::Multiply(*SM_Matrix_,false,*D0_Matrix_,false,temp,GetOStream(Runtime0),true,true);
          if (!implicitTranspose_)
            A22_ = Xpetra::MatrixMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>::Multiply(*D0_T_Matrix_,false,*temp,false,A22_,GetOStream(Runtime0),true,true);
          else
            A22_ = Xpetra::MatrixMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>::Multiply(*D0_Matrix_,true,*temp,false,A22_,GetOStream(Runtime0),true,true);
        } else {
          // we replaced domain map and importer on D0, reverse that
          RCP<const Import> D0importer = rcp_dynamic_cast<CrsMatrixWrap>(D0_Matrix_)->getCrsMatrix()->getCrsGraph()->getImporter();
          rcp_dynamic_cast<CrsMatrixWrap>(D0_Matrix_)->getCrsMatrix()->replaceDomainMapAndImporter(D0origDomainMap_, D0origImporter_);

          RCP<Matrix> temp, temp2;
          temp = Xpetra::MatrixMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>::Multiply(*SM_Matrix_,false,*D0_Matrix_,false,temp,GetOStream(Runtime0),true,true);
          if (!implicitTranspose_)
            temp2 = Xpetra::MatrixMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>::Multiply(*D0_T_Matrix_,false,*temp,false,temp2,GetOStream(Runtime0),true,true);
          else
            temp2 = Xpetra::MatrixMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>::Multiply(*D0_Matrix_,true,*temp,false,temp2,GetOStream(Runtime0),true,true);

          // and back again
          rcp_dynamic_cast<CrsMatrixWrap>(D0_Matrix_)->getCrsMatrix()->replaceDomainMapAndImporter(Importer22_->getTargetMap(), D0importer);

          ParameterList XpetraList;
          XpetraList.set("Restrict Communicator",true);
          XpetraList.set("Timer Label","MueLu::RebalanceA22");
          RCP<const Map> targetMap = Importer22_->getTargetMap();
          A22_ = MatrixFactory::Build(temp2, *Importer22_, *Importer22_, targetMap, targetMap, rcp(&XpetraList,false));
        }
      }

      if (!implicitTranspose_ && !reuse) {
        if (useKokkos_)
          D0_T_Matrix_ = Utilities_kokkos::Transpose(*D0_Matrix_);
        else
          D0_T_Matrix_ = Utilities::Transpose(*D0_Matrix_);
      }

      VerbLevel verbosityLevel = VerboseObject::GetDefaultVerbLevel();
      if (!A22_.is_null()) {
        dump(*A22_, "A22.m");
        A22_->setObjectLabel("RefMaxwell (2,2)");
        int oldRank = SetProcRankVerbose(A22_->getDomainMap()->getComm()->getRank());
        if (IsPrint(Statistics2)) {
          RCP<ParameterList> params = rcp(new ParameterList());;
          params->set("printLoadBalancingInfo", true);
          params->set("printCommInfo",          true);
          GetOStream(Statistics2) << PerfUtils::PrintMatrixInfo(*A22_, "A22", params);
        }
#if defined(HAVE_MUELU_STRATIMIKOS) && defined(HAVE_MUELU_THYRA)
        if (precList22_.isType<std::string>("Preconditioner Type")) {
          // build a Stratimikos preconditioner
          if (precList22_.get<std::string>("Preconditioner Type") == "MueLu") {
            ParameterList& userParamList = precList22_.sublist("Preconditioner Types").sublist("MueLu").sublist("user data");
            userParamList.set<RCP<RealValuedMultiVector> >("Coordinates", Coords_);
          }
          thyraPrecOp22_ = rcp(new XpetraThyraLinearOp<Scalar,LocalOrdinal,GlobalOrdinal,Node>(A22_, rcp(&precList22_, false)));
        } else
#endif
          {
          // build a MueLu hierarchy
          if (!reuse) {
            ParameterList& userParamList = precList22_.sublist("user data");
            userParamList.set<RCP<RealValuedMultiVector> >("Coordinates", Coords_);
            // If we detected no boundary conditions, the (2,2) problem is singular.
            // Therefore, if we want to use a direct coarse solver, we need to fix up the nullspace.
            std::string coarseType = "";
            if (precList22_.isParameter("coarse: type")) {
              coarseType = precList22_.get<std::string>("coarse: type");
              // Transform string to "Abcde" notation
              std::transform(coarseType.begin(),   coarseType.end(),   coarseType.begin(), ::tolower);
              std::transform(coarseType.begin(), ++coarseType.begin(), coarseType.begin(), ::toupper);
            }
            if (BCedges_ == 0 &&
                (coarseType == "" ||
                 coarseType == "Klu" ||
                 coarseType == "Klu2") &&
                (!precList22_.isSublist("coarse: params") ||
                 !precList22_.sublist("coarse: params").isParameter("fix nullspace")))
              precList22_.sublist("coarse: params").set("fix nullspace",true);
            Hierarchy22_ = MueLu::CreateXpetraPreconditioner(A22_, precList22_);
          } else {
            RCP<MueLu::Level> level0 = Hierarchy22_->GetLevel(0);
            level0->Set("A", A22_);
            Hierarchy22_->SetupRe();
          }
        }
        SetProcRankVerbose(oldRank);
      }
      VerboseObject::SetDefaultVerbLevel(verbosityLevel);

    }

    if(!reuse && !allNodesBoundary_ && applyBCsTo22_) {
      GetOStream(Runtime0) << "RefMaxwell::compute(): nuking BC edges of D0" << std::endl;

      D0_Matrix_->resumeFill();
      Scalar replaceWith;
      if (D0_Matrix_->getRowMap()->lib() == Xpetra::UseEpetra)
        replaceWith= Teuchos::ScalarTraits<SC>::eps();
      else
        replaceWith = Teuchos::ScalarTraits<SC>::zero();
      if (useKokkos_) {
        Utilities_kokkos::ZeroDirichletRows(D0_Matrix_,BCrowsKokkos_,replaceWith);
      } else {
        Utilities::ZeroDirichletRows(D0_Matrix_,BCrows_,replaceWith);
      }
      D0_Matrix_->fillComplete(D0_Matrix_->getDomainMap(),D0_Matrix_->getRangeMap());
      dump(*D0_Matrix_, "D0_nuked.m");
    }

    setFineLevelSmoother();

    if (!reuse) {
      if (!ImporterH_.is_null()) {
        RCP<const Import> ImporterP11 = ImportFactory::Build(ImporterH_->getTargetMap(),P11_->getColMap());
        rcp_dynamic_cast<CrsMatrixWrap>(P11_)->getCrsMatrix()->replaceDomainMapAndImporter(ImporterH_->getTargetMap(), ImporterP11);
      }

      if (!Importer22_.is_null()) {
        if (enable_reuse_) {
          D0origDomainMap_ = D0_Matrix_->getDomainMap();
          D0origImporter_ = rcp_dynamic_cast<CrsMatrixWrap>(D0_Matrix_)->getCrsMatrix()->getCrsGraph()->getImporter();
        }
        RCP<const Import> ImporterD0 = ImportFactory::Build(Importer22_->getTargetMap(),D0_Matrix_->getColMap());
        rcp_dynamic_cast<CrsMatrixWrap>(D0_Matrix_)->getCrsMatrix()->replaceDomainMapAndImporter(Importer22_->getTargetMap(), ImporterD0);
      }

#ifdef HAVE_MUELU_TPETRA
      if ((!D0_T_Matrix_.is_null()) &&
          (!R11_.is_null()) &&
          (!rcp_dynamic_cast<CrsMatrixWrap>(D0_T_Matrix_)->getCrsMatrix()->getCrsGraph()->getImporter().is_null()) &&
          (!rcp_dynamic_cast<CrsMatrixWrap>(R11_)->getCrsMatrix()->getCrsGraph()->getImporter().is_null()) &&
          (D0_T_Matrix_->getColMap()->lib() == Xpetra::UseTpetra) &&
          (R11_->getColMap()->lib() == Xpetra::UseTpetra))
        D0_T_R11_colMapsMatch_ = D0_T_Matrix_->getColMap()->isSameAs(*R11_->getColMap());
      else
#endif
        D0_T_R11_colMapsMatch_ = false;
      if (D0_T_R11_colMapsMatch_)
        GetOStream(Runtime0) << "RefMaxwell::compute(): D0_T and R11 have matching colMaps" << std::endl;

      // Allocate temporary MultiVectors for solve
      allocateMemory(1);

      if (parameterList_.isSublist("matvec params"))
        {
          RCP<ParameterList> matvecParams = rcpFromRef(parameterList_.sublist("matvec params"));
          Maxwell_Utils<SC,LO,GO,NO>::setMatvecParams(*D0_Matrix_, matvecParams);
          Maxwell_Utils<SC,LO,GO,NO>::setMatvecParams(*P11_, matvecParams);
          if (!D0_T_Matrix_.is_null()) Maxwell_Utils<SC,LO,GO,NO>::setMatvecParams(*D0_T_Matrix_, matvecParams);
          if (!R11_.is_null())         Maxwell_Utils<SC,LO,GO,NO>::setMatvecParams(*R11_, matvecParams);
          if (!ImporterH_.is_null())   ImporterH_->setDistributorParameters(matvecParams);
          if (!Importer22_.is_null())  Importer22_->setDistributorParameters(matvecParams);
        }
      if (!ImporterH_.is_null() && parameterList_.isSublist("refmaxwell: ImporterH params")){
        RCP<ParameterList> importerParams = rcpFromRef(parameterList_.sublist("refmaxwell: ImporterH params"));
        ImporterH_->setDistributorParameters(importerParams);
      }
      if (!Importer22_.is_null() && parameterList_.isSublist("refmaxwell: Importer22 params")){
        RCP<ParameterList> importerParams = rcpFromRef(parameterList_.sublist("refmaxwell: Importer22 params"));
        Importer22_->setDistributorParameters(importerParams);
      }
    }

    describe(GetOStream(Runtime0));

#ifdef HAVE_MUELU_CUDA
    if (parameterList_.get<bool>("refmaxwell: cuda profile setup", false)) cudaProfilerStop();
#endif
  }


  template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void RefMaxwell<Scalar,LocalOrdinal,GlobalOrdinal,Node>::setFineLevelSmoother() {

    Level level;
    RCP<MueLu::FactoryManagerBase> factoryHandler = rcp(new FactoryManager());
    level.SetFactoryManager(factoryHandler);
    level.SetLevelID(0);
    level.setObjectLabel("RefMaxwell (1,1)");
    level.Set("A",SM_Matrix_);
    level.setlib(SM_Matrix_->getDomainMap()->lib());
    // For Hiptmair
    level.Set("NodeMatrix", A22_);
    level.Set("D0", D0_Matrix_);

    if (parameterList_.isType<std::string>("smoother: pre type") && parameterList_.isType<std::string>("smoother: post type")) {
      std::string preSmootherType = parameterList_.get<std::string>("smoother: pre type");
      std::string postSmootherType = parameterList_.get<std::string>("smoother: post type");

      ParameterList preSmootherList, postSmootherList;
      if (parameterList_.isSublist("smoother: pre params"))
        preSmootherList = parameterList_.sublist("smoother: pre params");
      if (parameterList_.isSublist("smoother: post params"))
        postSmootherList = parameterList_.sublist("smoother: post params");

      RCP<SmootherPrototype> preSmootherPrototype = rcp(new TrilinosSmoother(preSmootherType, preSmootherList));
      RCP<SmootherPrototype> postSmootherPrototype = rcp(new TrilinosSmoother(postSmootherType, postSmootherList));
      RCP<SmootherFactory> smootherFact = rcp(new SmootherFactory(preSmootherPrototype, postSmootherPrototype));

      level.Request("PreSmoother",smootherFact.get());
      level.Request("PostSmoother",smootherFact.get());
      if (enable_reuse_) {
        ParameterList smootherFactoryParams;
        smootherFactoryParams.set("keep smoother data", true);
        smootherFact->SetParameterList(smootherFactoryParams);
        level.Request("PreSmoother data", smootherFact.get());
        level.Request("PostSmoother data", smootherFact.get());
        if (!PreSmootherData_.is_null())
          level.Set("PreSmoother data", PreSmootherData_, smootherFact.get());
        if (!PostSmootherData_.is_null())
          level.Set("PostSmoother data", PostSmootherData_, smootherFact.get());
      }
      smootherFact->Build(level);
      PreSmoother_ = level.Get<RCP<SmootherBase> >("PreSmoother",smootherFact.get());
      PostSmoother_ = level.Get<RCP<SmootherBase> >("PostSmoother",smootherFact.get());
      if (enable_reuse_) {
        PreSmootherData_ = level.Get<RCP<SmootherPrototype> >("PreSmoother data",smootherFact.get());
        PostSmootherData_ = level.Get<RCP<SmootherPrototype> >("PostSmoother data",smootherFact.get());
      }
    } else {
      std::string smootherType = parameterList_.get<std::string>("smoother: type", "CHEBYSHEV");

      RCP<SmootherPrototype> smootherPrototype = rcp(new TrilinosSmoother(smootherType, smootherList_));
      RCP<SmootherFactory> smootherFact = rcp(new SmootherFactory(smootherPrototype));
      level.Request("PreSmoother",smootherFact.get());
      if (enable_reuse_) {
        ParameterList smootherFactoryParams;
        smootherFactoryParams.set("keep smoother data", true);
        smootherFact->SetParameterList(smootherFactoryParams);
        level.Request("PreSmoother data", smootherFact.get());
        if (!PreSmootherData_.is_null())
          level.Set("PreSmoother data", PreSmootherData_, smootherFact.get());
      }
      smootherFact->Build(level);
      PreSmoother_ = level.Get<RCP<SmootherBase> >("PreSmoother",smootherFact.get());
      PostSmoother_ = PreSmoother_;
      if (enable_reuse_)
        PreSmootherData_ = level.Get<RCP<SmootherPrototype> >("PreSmoother data",smootherFact.get());
    }

  }


  template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void RefMaxwell<Scalar,LocalOrdinal,GlobalOrdinal,Node>::allocateMemory(int numVectors) const {

    RCP<Teuchos::TimeMonitor> tmAlloc = getTimer("MueLu RefMaxwell: Allocate MVs");

    if (!R11_.is_null())
      P11res_    = MultiVectorFactory::Build(R11_->getRangeMap(), numVectors);
    else
      P11res_    = MultiVectorFactory::Build(P11_->getDomainMap(), numVectors);
    P11res_->setObjectLabel("P11res");
    if (D0_T_R11_colMapsMatch_) {
      D0TR11Tmp_ = MultiVectorFactory::Build(R11_->getColMap(), numVectors);
      D0TR11Tmp_->setObjectLabel("D0TR11Tmp");
    }
    if (!ImporterH_.is_null()) {
      P11resTmp_ = MultiVectorFactory::Build(ImporterH_->getTargetMap(), numVectors);
      P11resTmp_->setObjectLabel("P11resTmp");
      P11x_      = MultiVectorFactory::Build(ImporterH_->getTargetMap(), numVectors);
    } else
      P11x_      = MultiVectorFactory::Build(P11_->getDomainMap(), numVectors);
    P11x_->setObjectLabel("P11x");
    if (!D0_T_Matrix_.is_null())
      D0res_     = MultiVectorFactory::Build(D0_T_Matrix_->getRangeMap(), numVectors);
    else
      D0res_     = MultiVectorFactory::Build(D0_Matrix_->getDomainMap(), numVectors);
    D0res_->setObjectLabel("D0res");
    if (!Importer22_.is_null()) {
      D0resTmp_ = MultiVectorFactory::Build(Importer22_->getTargetMap(), numVectors);
      D0resTmp_->setObjectLabel("D0resTmp");
      D0x_      = MultiVectorFactory::Build(Importer22_->getTargetMap(), numVectors);
    } else
      D0x_      = MultiVectorFactory::Build(D0_Matrix_->getDomainMap(), numVectors);
    D0x_->setObjectLabel("D0x");
    if (!AH_.is_null()) {
      if (!ImporterH_.is_null() && !implicitTranspose_)
        P11resSubComm_ = MultiVectorFactory::Build(P11resTmp_, Teuchos::View);
      else
        P11resSubComm_ = MultiVectorFactory::Build(P11res_, Teuchos::View);
      P11resSubComm_->replaceMap(AH_->getRangeMap());
      P11resSubComm_->setObjectLabel("P11resSubComm");

      P11xSubComm_ = MultiVectorFactory::Build(P11x_, Teuchos::View);
      P11xSubComm_->replaceMap(AH_->getDomainMap());
      P11xSubComm_->setObjectLabel("P11xSubComm");
    }
    if (!A22_.is_null()) {
      if (!Importer22_.is_null() && !implicitTranspose_)
        D0resSubComm_ = MultiVectorFactory::Build(D0resTmp_, Teuchos::View);
      else
        D0resSubComm_ = MultiVectorFactory::Build(D0res_, Teuchos::View);
      D0resSubComm_->replaceMap(A22_->getRangeMap());
      D0resSubComm_->setObjectLabel("D0resSubComm");

      D0xSubComm_ = MultiVectorFactory::Build(D0x_, Teuchos::View);
      D0xSubComm_->replaceMap(A22_->getDomainMap());
      D0xSubComm_->setObjectLabel("D0xSubComm");
    }
    residual_  = MultiVectorFactory::Build(SM_Matrix_->getDomainMap(), numVectors);
    residual_->setObjectLabel("residual");
  }


  template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void RefMaxwell<Scalar,LocalOrdinal,GlobalOrdinal,Node>::dump(const Matrix& A, std::string name) const {
    if (dump_matrices_) {
      GetOStream(Runtime0) << "Dumping to " << name << std::endl;
      Xpetra::IO<SC, LO, GO, NO>::Write(name, A);
    }
  }


  template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void RefMaxwell<Scalar,LocalOrdinal,GlobalOrdinal,Node>::dump(const MultiVector& X, std::string name) const {
    if (dump_matrices_) {
      GetOStream(Runtime0) << "Dumping to " << name << std::endl;
      Xpetra::IO<SC, LO, GO, NO>::Write(name, X);
    }
  }


  template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void RefMaxwell<Scalar,LocalOrdinal,GlobalOrdinal,Node>::dumpCoords(const RealValuedMultiVector& X, std::string name) const {
    if (dump_matrices_) {
      GetOStream(Runtime0) << "Dumping to " << name << std::endl;
      Xpetra::IO<coordinateType, LO, GO, NO>::Write(name, X);
    }
  }


  template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void RefMaxwell<Scalar,LocalOrdinal,GlobalOrdinal,Node>::dump(const Teuchos::ArrayRCP<bool>& v, std::string name) const {
    if (dump_matrices_) {
      GetOStream(Runtime0) << "Dumping to " << name << std::endl;
      std::ofstream out(name);
      for (size_t i = 0; i < Teuchos::as<size_t>(v.size()); i++)
        out << v[i] << "\n";
    }
  }

  template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void RefMaxwell<Scalar,LocalOrdinal,GlobalOrdinal,Node>::dump(const Kokkos::View<bool*, typename Node::device_type>& v, std::string name) const {
    if (dump_matrices_) {
      GetOStream(Runtime0) << "Dumping to " << name << std::endl;
      std::ofstream out(name);
      auto vH = Kokkos::create_mirror_view (v);
          Kokkos::deep_copy(vH , v);
          for (size_t i = 0; i < vH.size(); i++)
            out << vH[i] << "\n";
    }
  }

  template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  Teuchos::RCP<Teuchos::TimeMonitor> RefMaxwell<Scalar,LocalOrdinal,GlobalOrdinal,Node>::getTimer(std::string name, RCP<const Teuchos::Comm<int> > comm) const {
    if (IsPrint(Timings)) {
      if (!syncTimers_)
        return Teuchos::rcp(new Teuchos::TimeMonitor(*Teuchos::TimeMonitor::getNewTimer(name)));
      else {
        if (comm.is_null())
          return Teuchos::rcp(new Teuchos::SyncTimeMonitor(*Teuchos::TimeMonitor::getNewTimer(name), SM_Matrix_->getRowMap()->getComm().ptr()));
        else
          return Teuchos::rcp(new Teuchos::SyncTimeMonitor(*Teuchos::TimeMonitor::getNewTimer(name), comm.ptr()));
      }
    } else
      return Teuchos::null;
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
    size_t numLocalRows = SM_Matrix_->getLocalNumRows();

    RCP<Matrix> P_nodal;
    RCP<CrsMatrix> P_nodal_imported;
    RCP<MultiVector> Nullspace_nodal;
    if (skipFirstLevel_) {
      // build prolongator: algorithm 1 in the reference paper
      // First, build nodal unsmoothed prolongator using the matrix A_nodal
      bool read_P_from_file = parameterList_.get("refmaxwell: read_P_from_file",false);
      if (read_P_from_file) {
        // This permits to read in an ML prolongator, so that we get the same hierarchy.
        // (ML and MueLu typically produce different aggregates.)
        std::string P_filename = parameterList_.get("refmaxwell: P_filename",std::string("P.m"));
        std::string domainmap_filename = parameterList_.get("refmaxwell: P_domainmap_filename",std::string("domainmap_P.m"));
        std::string colmap_filename = parameterList_.get("refmaxwell: P_colmap_filename",std::string("colmap_P.m"));
        std::string coords_filename = parameterList_.get("refmaxwell: CoordsH",std::string("coordsH.m"));
        RCP<const Map> colmap = Xpetra::IO<SC, LO, GO, NO>::ReadMap(colmap_filename, A_nodal_Matrix_->getDomainMap()->lib(),A_nodal_Matrix_->getDomainMap()->getComm());
        RCP<const Map> domainmap = Xpetra::IO<SC, LO, GO, NO>::ReadMap(domainmap_filename, A_nodal_Matrix_->getDomainMap()->lib(),A_nodal_Matrix_->getDomainMap()->getComm());
        P_nodal = Xpetra::IO<SC, LO, GO, NO>::Read(P_filename, A_nodal_Matrix_->getDomainMap(), colmap, domainmap, A_nodal_Matrix_->getDomainMap());
        CoordsH_ = Xpetra::IO<typename RealValuedMultiVector::scalar_type, LO, GO, NO>::ReadMultiVector(coords_filename, domainmap);
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
        coarseLevel.setObjectLabel("RefMaxwell (1,1) A_nodal");
        fineLevel.setObjectLabel("RefMaxwell (1,1) A_nodal");

        LocalOrdinal NSdim = 1;
        RCP<MultiVector> nullSpace = MultiVectorFactory::Build(A_nodal_Matrix_->getRowMap(),NSdim);
        nullSpace->putScalar(SC_ONE);
        fineLevel.Set("Nullspace",nullSpace);

        RCP<Factory> amalgFact, dropFact, UncoupledAggFact, coarseMapFact, TentativePFact, Tfact, SaPFact;
        if (useKokkos_) {
          amalgFact = rcp(new AmalgamationFactory_kokkos());
          dropFact = rcp(new CoalesceDropFactory_kokkos());
          UncoupledAggFact = rcp(new UncoupledAggregationFactory_kokkos());
          coarseMapFact = rcp(new CoarseMapFactory_kokkos());
          TentativePFact = rcp(new TentativePFactory_kokkos());
          if (parameterList_.get("multigrid algorithm","unsmoothed") == "sa")
            SaPFact = rcp(new SaPFactory_kokkos());
          Tfact = rcp(new CoordinatesTransferFactory_kokkos());
        } else
          {
            amalgFact = rcp(new AmalgamationFactory());
            dropFact = rcp(new CoalesceDropFactory());
            UncoupledAggFact = rcp(new UncoupledAggregationFactory());
            coarseMapFact = rcp(new CoarseMapFactory());
            TentativePFact = rcp(new TentativePFactory());
            if (parameterList_.get("multigrid algorithm","unsmoothed") == "sa")
              SaPFact = rcp(new SaPFactory());
            Tfact = rcp(new CoordinatesTransferFactory());
          }
        dropFact->SetFactory("UnAmalgamationInfo", amalgFact);
        double dropTol = parameterList_.get("aggregation: drop tol",0.0);
        std::string dropScheme = parameterList_.get("aggregation: drop scheme","classical");
        std::string distLaplAlgo = parameterList_.get("aggregation: distance laplacian algo","default");
        dropFact->SetParameter("aggregation: drop tol",Teuchos::ParameterEntry(dropTol));
        dropFact->SetParameter("aggregation: drop scheme",Teuchos::ParameterEntry(dropScheme));
        if (!useKokkos_)
          dropFact->SetParameter("aggregation: distance laplacian algo",Teuchos::ParameterEntry(distLaplAlgo));

        UncoupledAggFact->SetFactory("Graph", dropFact);
        int minAggSize = parameterList_.get("aggregation: min agg size",2);
        UncoupledAggFact->SetParameter("aggregation: min agg size",Teuchos::ParameterEntry(minAggSize));
        int maxAggSize = parameterList_.get("aggregation: max agg size",-1);
        UncoupledAggFact->SetParameter("aggregation: max agg size",Teuchos::ParameterEntry(maxAggSize));

        coarseMapFact->SetFactory("Aggregates", UncoupledAggFact);

        TentativePFact->SetFactory("Aggregates", UncoupledAggFact);
        TentativePFact->SetFactory("UnAmalgamationInfo", amalgFact);
        TentativePFact->SetFactory("CoarseMap", coarseMapFact);

        Tfact->SetFactory("Aggregates", UncoupledAggFact);
        Tfact->SetFactory("CoarseMap", coarseMapFact);

        if (parameterList_.get("multigrid algorithm","unsmoothed") == "sa") {
          SaPFact->SetFactory("P", TentativePFact);
          coarseLevel.Request("P", SaPFact.get());
        } else
          coarseLevel.Request("P",TentativePFact.get());
        coarseLevel.Request("Nullspace",TentativePFact.get());
        coarseLevel.Request("Coordinates",Tfact.get());

        RCP<AggregationExportFactory> aggExport;
        if (parameterList_.get("aggregation: export visualization data",false)) {
          aggExport = rcp(new AggregationExportFactory());
          ParameterList aggExportParams;
          aggExportParams.set("aggregation: output filename", "aggs.vtk");
          aggExportParams.set("aggregation: output file: agg style", "Jacks");
          aggExport->SetParameterList(aggExportParams);

          aggExport->SetFactory("Aggregates", UncoupledAggFact);
          aggExport->SetFactory("UnAmalgamationInfo", amalgFact);
          fineLevel.Request("Aggregates",UncoupledAggFact.get());
          fineLevel.Request("UnAmalgamationInfo",amalgFact.get());
        }

        if (parameterList_.get("multigrid algorithm","unsmoothed") == "sa")
          coarseLevel.Get("P",P_nodal,SaPFact.get());
        else
          coarseLevel.Get("P",P_nodal,TentativePFact.get());
        coarseLevel.Get("Nullspace",Nullspace_nodal,TentativePFact.get());
        coarseLevel.Get("Coordinates",CoordsH_,Tfact.get());


        if (parameterList_.get("aggregation: export visualization data",false))
          aggExport->Build(fineLevel, coarseLevel);
      }
      dump(*P_nodal, "P_nodal.m");
      dump(*Nullspace_nodal, "Nullspace_nodal.m");


      RCP<CrsMatrix> D0Crs = rcp_dynamic_cast<CrsMatrixWrap>(D0_Matrix_)->getCrsMatrix();

      // Import off-rank rows of P_nodal into P_nodal_imported
      int numProcs = P_nodal->getDomainMap()->getComm()->getSize();
      if (numProcs > 1) {
        RCP<CrsMatrixWrap> P_nodal_temp;
        RCP<const Map> targetMap = D0Crs->getColMap();
        P_nodal_temp = rcp(new CrsMatrixWrap(targetMap));
        RCP<const Import> importer = D0Crs->getCrsGraph()->getImporter();
        P_nodal_temp->doImport(*P_nodal, *importer, Xpetra::INSERT);
        P_nodal_temp->fillComplete(rcp_dynamic_cast<CrsMatrixWrap>(P_nodal)->getCrsMatrix()->getDomainMap(),
                                   rcp_dynamic_cast<CrsMatrixWrap>(P_nodal)->getCrsMatrix()->getRangeMap());
        P_nodal_imported = P_nodal_temp->getCrsMatrix();
        dump(*P_nodal_temp, "P_nodal_imported.m");
      } else
        P_nodal_imported = rcp_dynamic_cast<CrsMatrixWrap>(P_nodal)->getCrsMatrix();

    }

    if (useKokkos_) {

      using ATS        = Kokkos::ArithTraits<SC>;
      using impl_Scalar = typename ATS::val_type;
      using impl_ATS = Kokkos::ArithTraits<impl_Scalar>;
      using range_type = Kokkos::RangePolicy<LO, typename NO::execution_space>;

      typedef typename Matrix::local_matrix_type KCRS;
      typedef typename KCRS::StaticCrsGraphType graph_t;
      typedef typename graph_t::row_map_type::non_const_type lno_view_t;
      typedef typename graph_t::entries_type::non_const_type lno_nnz_view_t;
      typedef typename KCRS::values_type::non_const_type scalar_view_t;

      const impl_Scalar impl_SC_ZERO = impl_ATS::zero();
      const impl_Scalar impl_SC_ONE = impl_ATS::one();
      const impl_Scalar impl_half = impl_SC_ONE / (impl_SC_ONE + impl_SC_ONE);


      // Which algorithm should we use for the construction of the special prolongator?
      // Option "mat-mat":
      //   Multiply D0 * P_nodal, take graph, blow up the domain space and compute the entries.
      std::string defaultAlgo = "mat-mat";
      std::string algo = parameterList_.get("refmaxwell: prolongator compute algorithm",defaultAlgo);

      if (skipFirstLevel_) {
	// Get data out of P_nodal_imported and D0.

        if (algo == "mat-mat") {
          RCP<Matrix> D0_P_nodal = MatrixFactory::Build(SM_Matrix_->getRowMap());
          Xpetra::MatrixMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>::Multiply(*D0_Matrix_,false,*P_nodal,false,*D0_P_nodal,true,true);

#ifdef HAVE_MUELU_DEBUG
          TEUCHOS_ASSERT(D0_P_nodal->getColMap()->isSameAs(*P_nodal_imported->getColMap()));
#endif

          // Get data out of D0*P.
          auto localD0P = D0_P_nodal->getLocalMatrixDevice();

          // Create the matrix object
          RCP<Map> blockColMap    = Xpetra::MapFactory<LO,GO,NO>::Build(P_nodal_imported->getColMap(), dim);
          RCP<Map> blockDomainMap = Xpetra::MapFactory<LO,GO,NO>::Build(P_nodal->getDomainMap(), dim);

          size_t nnzEstimate = dim*localD0P.graph.entries.size();
          lno_view_t P11rowptr("P11_rowptr", numLocalRows+1);
          lno_nnz_view_t P11colind("P11_colind",nnzEstimate);
          scalar_view_t P11vals("P11_vals",nnzEstimate);

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
                                   P11vals(dim*jj+k) = impl_SC_ZERO;
                                 }
                               });

          auto localNullspace = Nullspace_->getDeviceLocalView(Xpetra::Access::ReadOnly);

          // enter values
          if (D0_Matrix_->getLocalMaxNumRowEntries()>2) {
            // The matrix D0 has too many entries per row.
            // Therefore we need to check whether its entries are actually non-zero.
            // This is the case for the matrices built by MiniEM.
            GetOStream(Warnings0) << "RefMaxwell::buildProlongator(): D0 matrix has more than 2 entries per row. Taking inefficient code path." << std::endl;

            magnitudeType tol = Teuchos::ScalarTraits<magnitudeType>::eps();

	    auto localD0 = D0_Matrix_->getLocalMatrixDevice();
	    auto localP = P_nodal_imported->getLocalMatrixDevice();
            Kokkos::parallel_for("MueLu:RefMaxwell::buildProlongator_enterValues_D0wZeros", range_type(0,numLocalRows),
                                 KOKKOS_LAMBDA(const size_t i) {
                                   for (size_t ll = localD0.graph.row_map(i); ll < localD0.graph.row_map(i+1); ll++) {
                                     LO l = localD0.graph.entries(ll);
                                     impl_Scalar p = localD0.values(ll);
                                     if (impl_ATS::magnitude(p) < tol)
                                       continue;
                                     for (size_t jj = localP.graph.row_map(l); jj < localP.graph.row_map(l+1); jj++) {
                                       LO j = localP.graph.entries(jj);
                                       impl_Scalar v = localP.values(jj);
                                       for (size_t k = 0; k < dim; k++) {
                                         LO jNew = dim*j+k;
                                         impl_Scalar n = localNullspace(i,k);
                                         size_t m;
                                         for (m = P11rowptr(i); m < P11rowptr(i+1); m++)
                                           if (P11colind(m) == jNew)
                                             break;
#if defined(HAVE_MUELU_DEBUG) && !defined(HAVE_MUELU_CUDA) && !defined(HAVE_MUELU_HIP)
                                         TEUCHOS_ASSERT_EQUALITY(P11colind(m),jNew);
#endif
                                         P11vals(m) += impl_half * v * n;
                                       }
                                     }
                                   }
                                 });

          } else {
	    auto localD0 = D0_Matrix_->getLocalMatrixDevice();
	    auto localP = P_nodal_imported->getLocalMatrixDevice();
            Kokkos::parallel_for("MueLu:RefMaxwell::buildProlongator_enterValues", range_type(0,numLocalRows),
                                 KOKKOS_LAMBDA(const size_t i) {
                                   for (size_t ll = localD0.graph.row_map(i); ll < localD0.graph.row_map(i+1); ll++) {
                                     LO l = localD0.graph.entries(ll);
                                     for (size_t jj = localP.graph.row_map(l); jj < localP.graph.row_map(l+1); jj++) {
                                       LO j = localP.graph.entries(jj);
                                       impl_Scalar v = localP.values(jj);
                                       for (size_t k = 0; k < dim; k++) {
                                         LO jNew = dim*j+k;
                                         impl_Scalar n = localNullspace(i,k);
                                         size_t m;
                                         for (m = P11rowptr(i); m < P11rowptr(i+1); m++)
                                           if (P11colind(m) == jNew)
                                             break;
#if defined(HAVE_MUELU_DEBUG) && !defined(HAVE_MUELU_CUDA) && !defined(HAVE_MUELU_HIP)
                                         TEUCHOS_ASSERT_EQUALITY(P11colind(m),jNew);
#endif
                                         P11vals(m) += impl_half * v * n;
                                       }
                                     }
                                   }
                                 });
          }

          P11_ = rcp(new CrsMatrixWrap(SM_Matrix_->getRowMap(), blockColMap, 0));
          RCP<CrsMatrix> P11Crs = rcp_dynamic_cast<CrsMatrixWrap>(P11_)->getCrsMatrix();
          P11Crs->setAllValues(P11rowptr, P11colind, P11vals);
          P11Crs->expertStaticFillComplete(blockDomainMap, SM_Matrix_->getRangeMap());

        } else
          TEUCHOS_TEST_FOR_EXCEPTION(false,std::invalid_argument,algo << " is not a valid option for \"refmaxwell: prolongator compute algorithm\"");

        NullspaceH_ = MultiVectorFactory::Build(P11_->getDomainMap(), dim);

        auto localNullspace_nodal = Nullspace_nodal->getDeviceLocalView(Xpetra::Access::ReadOnly);
        auto localNullspaceH = NullspaceH_->getDeviceLocalView(Xpetra::Access::ReadWrite);
        Kokkos::parallel_for("MueLu:RefMaxwell::buildProlongator_nullspace", range_type(0,Nullspace_nodal->getLocalLength()),
                             KOKKOS_LAMBDA(const size_t i) {
                               impl_Scalar val = localNullspace_nodal(i,0);
                               for (size_t j = 0; j < dim; j++)
                                 localNullspaceH(dim*i+j, j) = val;
                             });

      } else { // !skipFirstLevel_
	// Get data out of P_nodal_imported and D0.
	auto localD0 = D0_Matrix_->getLocalMatrixDevice();

        CoordsH_ = Coords_;

        if (algo == "mat-mat") {

          // Create the matrix object
          RCP<Map> blockColMap    = Xpetra::MapFactory<LO,GO,NO>::Build(D0_Matrix_->getColMap(), dim);
          RCP<Map> blockDomainMap = Xpetra::MapFactory<LO,GO,NO>::Build(D0_Matrix_->getDomainMap(), dim);

          size_t nnzEstimate = dim*localD0.graph.entries.size();
          lno_view_t P11rowptr("P11_rowptr", numLocalRows+1);
          lno_nnz_view_t P11colind("P11_colind",nnzEstimate);
          scalar_view_t P11vals("P11_vals",nnzEstimate);

          // adjust rowpointer
          Kokkos::parallel_for("MueLu:RefMaxwell::buildProlongator_adjustRowptr", range_type(0,numLocalRows+1),
                               KOKKOS_LAMBDA(const size_t i) {
                                 P11rowptr(i) = dim*localD0.graph.row_map(i);
                               });

          // adjust column indices
          Kokkos::parallel_for("MueLu:RefMaxwell::buildProlongator_adjustColind", range_type(0,localD0.graph.entries.size()),
                               KOKKOS_LAMBDA(const size_t jj) {
                                 for (size_t k = 0; k < dim; k++) {
                                   P11colind(dim*jj+k) = dim*localD0.graph.entries(jj)+k;
                                   P11vals(dim*jj+k) = impl_SC_ZERO;
                                 }
                               });

          auto localNullspace = Nullspace_->getDeviceLocalView(Xpetra::Access::ReadOnly);

          // enter values
          if (D0_Matrix_->getLocalMaxNumRowEntries()>2) {
            // The matrix D0 has too many entries per row.
            // Therefore we need to check whether its entries are actually non-zero.
            // This is the case for the matrices built by MiniEM.
            GetOStream(Warnings0) << "RefMaxwell::buildProlongator(): D0 matrix has more than 2 entries per row. Taking inefficient code path." << std::endl;

            magnitudeType tol = Teuchos::ScalarTraits<magnitudeType>::eps();

            Kokkos::parallel_for("MueLu:RefMaxwell::buildProlongator_enterValues_D0wZeros", range_type(0,numLocalRows),
                                 KOKKOS_LAMBDA(const size_t i) {
                                   for (size_t jj = localD0.graph.row_map(i); jj < localD0.graph.row_map(i+1); jj++) {
                                     LO j = localD0.graph.entries(jj);
                                     impl_Scalar p = localD0.values(jj);
                                     if (impl_ATS::magnitude(p) < tol)
                                       continue;
                                     for (size_t k = 0; k < dim; k++) {
                                       LO jNew = dim*j+k;
                                       impl_Scalar n = localNullspace(i,k);
                                       size_t m;
                                       for (m = P11rowptr(i); m < P11rowptr(i+1); m++)
                                         if (P11colind(m) == jNew)
                                           break;
#if defined(HAVE_MUELU_DEBUG) && !defined(HAVE_MUELU_CUDA) && !defined(HAVE_MUELU_HIP)
                                       TEUCHOS_ASSERT_EQUALITY(P11colind(m),jNew);
#endif
                                       P11vals(m) += impl_half * n;
                                     }
                                   }
                                 });

          } else {
            Kokkos::parallel_for("MueLu:RefMaxwell::buildProlongator_enterValues", range_type(0,numLocalRows),
                                 KOKKOS_LAMBDA(const size_t i) {
                                   for (size_t jj = localD0.graph.row_map(i); jj < localD0.graph.row_map(i+1); jj++) {
                                     LO j = localD0.graph.entries(jj);
                                     for (size_t k = 0; k < dim; k++) {
                                       LO jNew = dim*j+k;
                                       impl_Scalar n = localNullspace(i,k);
                                       size_t m;
                                       for (m = P11rowptr(i); m < P11rowptr(i+1); m++)
                                         if (P11colind(m) == jNew)
                                           break;
#if defined(HAVE_MUELU_DEBUG) && !defined(HAVE_MUELU_CUDA) && !defined(HAVE_MUELU_HIP)
                                       TEUCHOS_ASSERT_EQUALITY(P11colind(m),jNew);
#endif
                                       P11vals(m) += impl_half * n;
                                     }
                                   }
                                 });
          }

          P11_ = rcp(new CrsMatrixWrap(SM_Matrix_->getRowMap(), blockColMap, 0));
          RCP<CrsMatrix> P11Crs = rcp_dynamic_cast<CrsMatrixWrap>(P11_)->getCrsMatrix();
          P11Crs->setAllValues(P11rowptr, P11colind, P11vals);
          P11Crs->expertStaticFillComplete(blockDomainMap, SM_Matrix_->getRangeMap());
        } else
          TEUCHOS_TEST_FOR_EXCEPTION(false,std::invalid_argument,algo << " is not a valid option for \"refmaxwell: prolongator compute algorithm\"");

      }
    } else
      {
        // get nullspace vectors
        ArrayRCP<ArrayRCP<const SC> > nullspaceRCP(dim);
        ArrayRCP<ArrayView<const SC> > nullspace(dim);
        for(size_t i=0; i<dim; i++) {
          nullspaceRCP[i] = Nullspace_->getData(i);
          nullspace[i] = nullspaceRCP[i]();
        }

        // Get data out of P_nodal_imported and D0.
        ArrayRCP<size_t>            P11rowptr_RCP;
        ArrayRCP<LO>                P11colind_RCP;
        ArrayRCP<SC>                P11vals_RCP;


        // Which algorithm should we use for the construction of the special prolongator?
        // Option "mat-mat":
        //   Multiply D0 * P_nodal, take graph, blow up the domain space and compute the entries.
        // Option "gustavson":
        //   Loop over D0, P and nullspace and allocate directly. (Gustavson-like)
        //   More efficient, but only available for serial node.
        std::string defaultAlgo = "mat-mat";
        std::string algo = parameterList_.get("refmaxwell: prolongator compute algorithm",defaultAlgo);

        if (skipFirstLevel_) {

          if (algo == "mat-mat") {
            RCP<Matrix> D0_P_nodal = MatrixFactory::Build(SM_Matrix_->getRowMap());
            Xpetra::MatrixMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>::Multiply(*D0_Matrix_,false,*P_nodal,false,*D0_P_nodal,true,true);


	    ArrayRCP<const size_t>      D0rowptr_RCP;
	    ArrayRCP<const LO>          D0colind_RCP;
	    ArrayRCP<const SC>          D0vals_RCP;
	    rcp_dynamic_cast<CrsMatrixWrap>(D0_Matrix_)->getCrsMatrix()->getAllValues(D0rowptr_RCP, D0colind_RCP, D0vals_RCP);
	    // For efficiency
	    // Refers to an issue where Teuchos::ArrayRCP::operator[] may be
	    // slower than Teuchos::ArrayView::operator[].
	    ArrayView<const size_t>     D0rowptr;
	    ArrayView<const LO>         D0colind;
	    ArrayView<const SC>         D0vals;
	    D0rowptr = D0rowptr_RCP();  D0colind = D0colind_RCP();  D0vals = D0vals_RCP();

	    // Get data out of P_nodal_imported and D0.
	    ArrayRCP<const size_t>      Prowptr_RCP;
	    ArrayRCP<const LO>          Pcolind_RCP;
	    ArrayRCP<const SC>          Pvals_RCP;
	    P_nodal_imported->getAllValues(Prowptr_RCP, Pcolind_RCP, Pvals_RCP);
	    ArrayView<const size_t>     Prowptr;
	    ArrayView<const LO>         Pcolind;
	    ArrayView<const SC>         Pvals;
	    Prowptr  = Prowptr_RCP();   Pcolind  = Pcolind_RCP();   Pvals = Pvals_RCP();

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
            P11_ = rcp(new CrsMatrixWrap(SM_Matrix_->getRowMap(), blockColMap, 0));
            RCP<CrsMatrix> P11Crs = rcp_dynamic_cast<CrsMatrixWrap>(P11_)->getCrsMatrix();
            size_t nnzEstimate = dim*D0Prowptr[numLocalRows];
            P11Crs->allocateAllValues(nnzEstimate, P11rowptr_RCP, P11colind_RCP, P11vals_RCP);

            ArrayView<size_t> P11rowptr = P11rowptr_RCP();
            ArrayView<LO>     P11colind = P11colind_RCP();
            ArrayView<SC>     P11vals   = P11vals_RCP();

            // adjust rowpointer
            for (size_t i = 0; i < numLocalRows+1; i++) {
              P11rowptr[i] = dim*D0Prowptr[i];
            }

            // adjust column indices
            for (size_t jj = 0; jj < (size_t) D0Prowptr[numLocalRows]; jj++)
              for (size_t k = 0; k < dim; k++) {
                P11colind[dim*jj+k] = dim*D0Pcolind[jj]+k;
                P11vals[dim*jj+k] = SC_ZERO;
              }

            RCP<const Map> P_nodal_imported_colmap = P_nodal_imported->getColMap();
            RCP<const Map> D0_P_nodal_colmap = D0_P_nodal->getColMap();
            // enter values
            if (D0_Matrix_->getLocalMaxNumRowEntries()>2) {
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
                    j = D0_P_nodal_colmap->getLocalElement(P_nodal_imported_colmap->getGlobalElement(j));
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
                    j = D0_P_nodal_colmap->getLocalElement(P_nodal_imported_colmap->getGlobalElement(j));
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
	    ArrayRCP<const size_t>      D0rowptr_RCP;
	    ArrayRCP<const LO>          D0colind_RCP;
	    ArrayRCP<const SC>          D0vals_RCP;
	    rcp_dynamic_cast<CrsMatrixWrap>(D0_Matrix_)->getCrsMatrix()->getAllValues(D0rowptr_RCP, D0colind_RCP, D0vals_RCP);
	    // For efficiency
	    // Refers to an issue where Teuchos::ArrayRCP::operator[] may be
	    // slower than Teuchos::ArrayView::operator[].
	    ArrayView<const size_t>     D0rowptr;
	    ArrayView<const LO>         D0colind;
	    ArrayView<const SC>         D0vals;
	    D0rowptr = D0rowptr_RCP();  D0colind = D0colind_RCP();  D0vals = D0vals_RCP();

	    // Get data out of P_nodal_imported and D0.
	    ArrayRCP<const size_t>      Prowptr_RCP;
	    ArrayRCP<const LO>          Pcolind_RCP;
	    ArrayRCP<const SC>          Pvals_RCP;
	    P_nodal_imported->getAllValues(Prowptr_RCP, Pcolind_RCP, Pvals_RCP);
	    ArrayView<const size_t>     Prowptr;
	    ArrayView<const LO>         Pcolind;
	    ArrayView<const SC>         Pvals;
	    Prowptr  = Prowptr_RCP();   Pcolind  = Pcolind_RCP();   Pvals = Pvals_RCP();

            LO maxP11col = dim * P_nodal_imported->getColMap()->getMaxLocalIndex();
            const size_t ST_INVALID = Teuchos::OrdinalTraits<LO>::invalid();
            Array<size_t> P11_status(dim*maxP11col, ST_INVALID);
            // This is ad-hoc and should maybe be replaced with some better heuristics.
            size_t nnz_alloc = dim*D0vals_RCP.size();

            // Create the matrix object
            RCP<Map> blockColMap    = Xpetra::MapFactory<LO,GO,NO>::Build(P_nodal_imported->getColMap(), dim);
            RCP<Map> blockDomainMap = Xpetra::MapFactory<LO,GO,NO>::Build(P_nodal->getDomainMap(), dim);
            P11_ = rcp(new CrsMatrixWrap(SM_Matrix_->getRowMap(), blockColMap, 0));
            RCP<CrsMatrix> P11Crs = rcp_dynamic_cast<CrsMatrixWrap>(P11_)->getCrsMatrix();
            P11Crs->allocateAllValues(nnz_alloc, P11rowptr_RCP, P11colind_RCP, P11vals_RCP);

            ArrayView<size_t> P11rowptr = P11rowptr_RCP();
            ArrayView<LO>     P11colind = P11colind_RCP();
            ArrayView<SC>     P11vals   = P11vals_RCP();

            size_t nnz;
            if (D0_Matrix_->getLocalMaxNumRowEntries()>2) {
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

          NullspaceH_ = MultiVectorFactory::Build(P11_->getDomainMap(), dim);

          ArrayRCP<const Scalar> ns_rcp = Nullspace_nodal->getData(0);
          ArrayView<const Scalar> ns_view = ns_rcp();
          for (size_t i = 0; i < Nullspace_nodal->getLocalLength(); i++) {
            Scalar val = ns_view[i];
            for (size_t j = 0; j < dim; j++)
              NullspaceH_->replaceLocalValue(dim*i+j, j, val);
          }


        } else { // !skipFirstLevel_
	  ArrayRCP<const size_t>      D0rowptr_RCP;
	  ArrayRCP<const LO>          D0colind_RCP;
	  ArrayRCP<const SC>          D0vals_RCP;
	  rcp_dynamic_cast<CrsMatrixWrap>(D0_Matrix_)->getCrsMatrix()->getAllValues(D0rowptr_RCP, D0colind_RCP, D0vals_RCP);
	  // For efficiency
	  // Refers to an issue where Teuchos::ArrayRCP::operator[] may be
	  // slower than Teuchos::ArrayView::operator[].
	  ArrayView<const size_t>     D0rowptr;
	  ArrayView<const LO>         D0colind;
	  ArrayView<const SC>         D0vals;
	  D0rowptr = D0rowptr_RCP();  D0colind = D0colind_RCP();  D0vals = D0vals_RCP();

          CoordsH_ = Coords_;
          if (algo == "mat-mat") {

            // Create the matrix object
            RCP<Map> blockColMap    = Xpetra::MapFactory<LO,GO,NO>::Build(D0_Matrix_->getColMap(), dim);
            RCP<Map> blockDomainMap = Xpetra::MapFactory<LO,GO,NO>::Build(D0_Matrix_->getDomainMap(), dim);
            P11_ = rcp(new CrsMatrixWrap(SM_Matrix_->getRowMap(), blockColMap, 0));
            RCP<CrsMatrix> P11Crs = rcp_dynamic_cast<CrsMatrixWrap>(P11_)->getCrsMatrix();
            size_t nnzEstimate = dim*D0rowptr[numLocalRows];
            P11Crs->allocateAllValues(nnzEstimate, P11rowptr_RCP, P11colind_RCP, P11vals_RCP);

            ArrayView<size_t> P11rowptr = P11rowptr_RCP();
            ArrayView<LO>     P11colind = P11colind_RCP();
            ArrayView<SC>     P11vals   = P11vals_RCP();

            // adjust rowpointer
            for (size_t i = 0; i < numLocalRows+1; i++) {
              P11rowptr[i] = dim*D0rowptr[i];
            }

            // adjust column indices
            for (size_t jj = 0; jj < (size_t) D0rowptr[numLocalRows]; jj++)
              for (size_t k = 0; k < dim; k++) {
                P11colind[dim*jj+k] = dim*D0colind[jj]+k;
                P11vals[dim*jj+k] = SC_ZERO;
              }

            // enter values
            if (D0_Matrix_->getLocalMaxNumRowEntries()>2) {
              // The matrix D0 has too many entries per row.
              // Therefore we need to check whether its entries are actually non-zero.
              // This is the case for the matrices built by MiniEM.
              GetOStream(Warnings0) << "RefMaxwell::buildProlongator(): D0 matrix has more than 2 entries per row. Taking inefficient code path." << std::endl;

              magnitudeType tol = Teuchos::ScalarTraits<magnitudeType>::eps();
              for (size_t i = 0; i < numLocalRows; i++) {
                for (size_t jj = D0rowptr[i]; jj < D0rowptr[i+1]; jj++) {
                  LO j = D0colind[jj];
                  SC p = D0vals[jj];
                  if (Teuchos::ScalarTraits<Scalar>::magnitude(p) < tol)
                    continue;
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
                    P11vals[m] += half * n;
                  }
                }
              }
            } else {
              // enter values
              for (size_t i = 0; i < numLocalRows; i++) {
                for (size_t jj = D0rowptr[i]; jj < D0rowptr[i+1]; jj++) {
                  LO j = D0colind[jj];

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
                    P11vals[m] += half * n;
                  }
                }
              }
            }

            P11Crs->setAllValues(P11rowptr_RCP, P11colind_RCP, P11vals_RCP);
            P11Crs->expertStaticFillComplete(blockDomainMap, SM_Matrix_->getRangeMap());

          } else
            TEUCHOS_TEST_FOR_EXCEPTION(false,std::invalid_argument,algo << " is not a valid option for \"refmaxwell: prolongator compute algorithm\"");
        }
      }
  }

  template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void RefMaxwell<Scalar,LocalOrdinal,GlobalOrdinal,Node>::formCoarseMatrix() {
    RCP<Teuchos::TimeMonitor> tm = getTimer("MueLu RefMaxwell: Build coarse (1,1) matrix");

    const Scalar one = Teuchos::ScalarTraits<Scalar>::one();

    // coarse matrix for P11* (M1 + D1* M2 D1) P11
    RCP<Matrix> temp;
    temp = Xpetra::MatrixMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>::Multiply(*SM_Matrix_,false,*P11_,false,temp,GetOStream(Runtime0),true,true);
    if (ImporterH_.is_null())
      AH_ = Xpetra::MatrixMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>::Multiply(*P11_,true,*temp,false,AH_,GetOStream(Runtime0),true,true);
    else {

      RCP<Matrix> temp2;
      temp2 = Xpetra::MatrixMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>::Multiply(*P11_,true,*temp,false,temp2,GetOStream(Runtime0),true,true);

      RCP<const Map> map = ImporterH_->getTargetMap()->removeEmptyProcesses();
      temp2->removeEmptyProcessesInPlace(map);
      if (!temp2.is_null() && temp2->getRowMap().is_null())
        temp2 = Teuchos::null;
      AH_ = temp2;
    }

    if (!disable_addon_) {

      RCP<Matrix> addon;

      if (!AH_.is_null() && Addon_Matrix_.is_null()) {
        // construct addon
        RCP<Teuchos::TimeMonitor> tmAddon = getTimer("MueLu RefMaxwell: Build coarse addon matrix");
        // catch a failure
        TEUCHOS_TEST_FOR_EXCEPTION(M0inv_Matrix_==Teuchos::null,std::invalid_argument,
                                   "MueLu::RefMaxwell::formCoarseMatrix(): Inverse of "
                                   "lumped mass matrix required for add-on (i.e. M0inv_Matrix is null)");

        // coarse matrix for add-on, i.e P11* (M1 D0 M0inv D0* M1) P11
        RCP<Matrix> Zaux, Z;

        // construct Zaux = M1 P11
        Zaux = Xpetra::MatrixMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>::Multiply(*M1_Matrix_,false,*P11_,false,Zaux,GetOStream(Runtime0),true,true);
        // construct Z = D0* M1 P11 = D0* Zaux
        Z = Xpetra::MatrixMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>::Multiply(*D0_Matrix_,true,*Zaux,false,Z,GetOStream(Runtime0),true,true);

        // construct Z* M0inv Z
        if (M0inv_Matrix_->getGlobalMaxNumRowEntries()<=1) {
          // We assume that if M0inv has at most one entry per row then
          // these are all diagonal entries.
          RCP<Vector> diag = VectorFactory::Build(M0inv_Matrix_->getRowMap());
          M0inv_Matrix_->getLocalDiagCopy(*diag);
	  {
	    ArrayRCP<Scalar> diagVals = diag->getDataNonConst(0);
	    for (size_t j=0; j < diag->getMap()->getLocalNumElements(); j++) {
	      diagVals[j] = Teuchos::ScalarTraits<Scalar>::squareroot(diagVals[j]);
	    }
	  }
          if (Z->getRowMap()->isSameAs(*(diag->getMap())))
            Z->leftScale(*diag);
          else {
            RCP<Import> importer = ImportFactory::Build(diag->getMap(),Z->getRowMap());
            RCP<Vector> diag2 = VectorFactory::Build(Z->getRowMap());
            diag2->doImport(*diag,*importer,Xpetra::INSERT);
            Z->leftScale(*diag2);
          }
          addon = Xpetra::MatrixMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>::Multiply(*Z,true,*Z,false,addon,GetOStream(Runtime0),true,true);
        } else if (parameterList_.get<bool>("rap: triple product", false) == false) {
          RCP<Matrix> C2;
          // construct C2 = M0inv Z
          C2 = Xpetra::MatrixMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>::Multiply(*M0inv_Matrix_,false,*Z,false,C2,GetOStream(Runtime0),true,true);
          // construct Matrix2 = Z* M0inv Z = Z* C2
          addon = Xpetra::MatrixMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>::Multiply(*Z,true,*C2,false,addon,GetOStream(Runtime0),true,true);
        } else {
          addon = MatrixFactory::Build(Z->getDomainMap());
          // construct Matrix2 = Z* M0inv Z
          Xpetra::TripleMatrixMultiply<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
            MultiplyRAP(*Z, true, *M0inv_Matrix_, false, *Z, false, *addon, true, true);
        }
        // Should we keep the addon for next setup?
        if (enable_reuse_)
          Addon_Matrix_ = addon;
      } else
        addon = Addon_Matrix_;

      if (!AH_.is_null()) {
        // add matrices together
        RCP<Matrix> newAH;
        Xpetra::MatrixMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>::TwoMatrixAdd(*AH_,false,one,*addon,false,one,newAH,GetOStream(Runtime0));
        newAH->fillComplete();
        AH_ = newAH;
      }
    }

    if (!AH_.is_null() && !skipFirstLevel_) {
      ArrayRCP<bool> AHBCrows;
      AHBCrows.resize(AH_->getRowMap()->getLocalNumElements());
      size_t dim = Nullspace_->getNumVectors();
      if (useKokkos_)
        for (size_t i = 0; i < BCdomainKokkos_.size(); i++)
          for (size_t k = 0; k < dim; k++)
            AHBCrows[i*dim+k] = BCdomainKokkos_(i);
      else
        for (size_t i = 0; i < static_cast<size_t>(BCdomain_.size()); i++)
          for (size_t k = 0; k < dim; k++)
            AHBCrows[i*dim+k] = BCdomain_[i];
      magnitudeType rowSumTol = parameterList_.get("refmaxwell: row sum drop tol (1,1)",-1.0);
      if (rowSumTol > 0.)
        Utilities::ApplyRowSumCriterion(*AH_, rowSumTol, AHBCrows);
      if (applyBCsToH_)
        Utilities::ApplyOAZToMatrixRows(AH_, AHBCrows);
    }

    if (!AH_.is_null()) {
      // If we already applied BCs to A_nodal, we likely do not need
      // to fix up AH.
      // If we did not apply BCs to A_nodal, we now need to correct
      // the zero diagonals of AH, since we did nuke the nullspace.

      bool fixZeroDiagonal = !applyBCsToAnodal_;
      if (precList11_.isParameter("rap: fix zero diagonals"))
        fixZeroDiagonal = precList11_.get<bool>("rap: fix zero diagonals");

      if (fixZeroDiagonal) {
        magnitudeType threshold = 1e-16;
        Scalar replacement = 1.0;
        if (precList11_.isType<magnitudeType>("rap: fix zero diagonals threshold"))
          threshold = precList11_.get<magnitudeType>("rap: fix zero diagonals threshold");
        else if (precList11_.isType<double>("rap: fix zero diagonals threshold"))
          threshold = Teuchos::as<magnitudeType>(precList11_.get<double>("rap: fix zero diagonals threshold"));
        if (precList11_.isType<double>("rap: fix zero diagonals replacement"))
          replacement = Teuchos::as<Scalar>(precList11_.get<double>("rap: fix zero diagonals replacement"));
        Xpetra::MatrixUtils<SC,LO,GO,NO>::CheckRepairMainDiagonal(AH_, true, GetOStream(Warnings1), threshold, replacement);
      }

      // Set block size
      size_t dim = Nullspace_->getNumVectors();
      AH_->SetFixedBlockSize(dim);
    }
  }


  template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void RefMaxwell<Scalar,LocalOrdinal,GlobalOrdinal,Node>::resetMatrix(RCP<Matrix> SM_Matrix_new, bool ComputePrec) {
    bool reuse = !SM_Matrix_.is_null();
    SM_Matrix_ = SM_Matrix_new;
    dump(*SM_Matrix_, "SM.m");
    if (ComputePrec)
      compute(reuse);
  }


  template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void RefMaxwell<Scalar,LocalOrdinal,GlobalOrdinal,Node>::applyInverseAdditive(const MultiVector& RHS, MultiVector& X) const {

    Scalar one = Teuchos::ScalarTraits<Scalar>::one();

    { // compute residual

      RCP<Teuchos::TimeMonitor> tmRes = getTimer("MueLu RefMaxwell: residual calculation");
      Utilities::Residual(*SM_Matrix_, X, RHS, *residual_);
    }

    { // restrict residual to sub-hierarchies

      if (implicitTranspose_) {
        {
          RCP<Teuchos::TimeMonitor> tmRes = getTimer("MueLu RefMaxwell: restriction coarse (1,1) (implicit)");
          P11_->apply(*residual_,*P11res_,Teuchos::TRANS);
        }
        if (!allNodesBoundary_) {
          RCP<Teuchos::TimeMonitor> tmD0 = getTimer("MueLu RefMaxwell: restriction (2,2) (implicit)");
          D0_Matrix_->apply(*residual_,*D0res_,Teuchos::TRANS);
        }
      } else {
#ifdef MUELU_HAVE_TPETRA
        if (D0_T_R11_colMapsMatch_) {
          // Column maps of D0_T and R11 match, and we're running Tpetra
          {
            RCP<Teuchos::TimeMonitor> tmD0 = getTimer("MueLu RefMaxwell: restrictions import");
            D0TR11Tmp_->doImport(*residual_, *rcp_dynamic_cast<CrsMatrixWrap>(R11_)->getCrsMatrix()->getCrsGraph()->getImporter(), Xpetra::INSERT);
          }
          if (!allNodesBoundary_) {
            RCP<Teuchos::TimeMonitor> tmD0 = getTimer("MueLu RefMaxwell: restriction (2,2) (explicit)");
            rcp_dynamic_cast<TpetraCrsMatrix>(rcp_dynamic_cast<CrsMatrixWrap>(D0_T_Matrix_)->getCrsMatrix())->getTpetra_CrsMatrix()->localApply(toTpetra(*D0TR11Tmp_),toTpetra(*D0res_),Teuchos::NO_TRANS);
          }
          {
            RCP<Teuchos::TimeMonitor> tmP11 = getTimer("MueLu RefMaxwell: restriction coarse (1,1) (explicit)");
            rcp_dynamic_cast<TpetraCrsMatrix>(rcp_dynamic_cast<CrsMatrixWrap>(R11_)->getCrsMatrix())->getTpetra_CrsMatrix()->localApply(toTpetra(*D0TR11Tmp_),toTpetra(*P11res_),Teuchos::NO_TRANS);
          }
        } else
#endif
        {
          {
            RCP<Teuchos::TimeMonitor> tmP11 = getTimer("MueLu RefMaxwell: restriction coarse (1,1) (explicit)");
            R11_->apply(*residual_,*P11res_,Teuchos::NO_TRANS);
          }
          if (!allNodesBoundary_) {
            RCP<Teuchos::TimeMonitor> tmD0 = getTimer("MueLu RefMaxwell: restriction (2,2) (explicit)");
            D0_T_Matrix_->apply(*residual_,*D0res_,Teuchos::NO_TRANS);
          }
        }
      }
    }

    {
      RCP<Teuchos::TimeMonitor> tmSubSolves = getTimer("MueLu RefMaxwell: subsolves");

      // block diagonal preconditioner on 2x2 (V-cycle for diagonal blocks)

      if (!ImporterH_.is_null() && !implicitTranspose_) {
        RCP<Teuchos::TimeMonitor> tmH = getTimer("MueLu RefMaxwell: import coarse (1,1)");
        P11resTmp_->beginImport(*P11res_, *ImporterH_, Xpetra::INSERT);
      }
      if (!allNodesBoundary_ && !Importer22_.is_null() && !implicitTranspose_) {
        RCP<Teuchos::TimeMonitor> tm22 = getTimer("MueLu RefMaxwell: import (2,2)");
        D0resTmp_->beginImport(*D0res_, *Importer22_, Xpetra::INSERT);
      }

      // iterate on coarse (1, 1) block
      if (!AH_.is_null()) {
        if (!ImporterH_.is_null() && !implicitTranspose_)
          P11resTmp_->endImport(*P11res_, *ImporterH_, Xpetra::INSERT);

        RCP<Teuchos::TimeMonitor> tmH = getTimer("MueLu RefMaxwell: solve coarse (1,1)", AH_->getRowMap()->getComm());

#if defined(HAVE_MUELU_STRATIMIKOS) && defined(HAVE_MUELU_THYRA)
        if (!thyraPrecOpH_.is_null()) {
          Scalar zero = Teuchos::ScalarTraits<Scalar>::zero();
          thyraPrecOpH_->apply(*P11resSubComm_, *P11xSubComm_, Teuchos::NO_TRANS, one, zero);
        }
        else
#endif
           HierarchyH_->Iterate(*P11resSubComm_, *P11xSubComm_, numItersH_, true);
      }

      // iterate on (2, 2) block
      if (!A22_.is_null()) {
        if (!allNodesBoundary_ && !Importer22_.is_null() && !implicitTranspose_)
          D0resTmp_->endImport(*D0res_, *Importer22_, Xpetra::INSERT);

        RCP<Teuchos::TimeMonitor> tm22 = getTimer("MueLu RefMaxwell: solve (2,2)", A22_->getRowMap()->getComm());

#if defined(HAVE_MUELU_STRATIMIKOS) && defined(HAVE_MUELU_THYRA)
        if (!thyraPrecOp22_.is_null()) {
          Scalar zero = Teuchos::ScalarTraits<Scalar>::zero();
          thyraPrecOp22_->apply(*D0resSubComm_, *D0xSubComm_, Teuchos::NO_TRANS, one, zero);
        }
        else
#endif
          Hierarchy22_->Iterate(*D0resSubComm_, *D0xSubComm_, numIters22_, true);
      }

      if (AH_.is_null() && !ImporterH_.is_null() && !implicitTranspose_)
        P11resTmp_->endImport(*P11res_, *ImporterH_, Xpetra::INSERT);
      if (A22_.is_null() && !allNodesBoundary_ && !Importer22_.is_null() && !implicitTranspose_)
        D0resTmp_->endImport(*D0res_, *Importer22_, Xpetra::INSERT);
    }

    if (fuseProlongationAndUpdate_) {
      { // prolongate (1,1) block
        RCP<Teuchos::TimeMonitor> tmP11 = getTimer("MueLu RefMaxwell: prolongation coarse (1,1) (fused)");
        P11_->apply(*P11x_,X,Teuchos::NO_TRANS,one,one);
      }

      if (!allNodesBoundary_) { // prolongate (2,2) block
        RCP<Teuchos::TimeMonitor> tmD0 = getTimer("MueLu RefMaxwell: prolongation (2,2) (fused)");
        D0_Matrix_->apply(*D0x_,X,Teuchos::NO_TRANS,one,one);
      }
    } else {
      { // prolongate (1,1) block
        RCP<Teuchos::TimeMonitor> tmP11 = getTimer("MueLu RefMaxwell: prolongation coarse (1,1) (unfused)");
        P11_->apply(*P11x_,*residual_,Teuchos::NO_TRANS);
      }

      if (!allNodesBoundary_) { // prolongate (2,2) block
        RCP<Teuchos::TimeMonitor> tmD0 = getTimer("MueLu RefMaxwell: prolongation (2,2) (unfused)");
        D0_Matrix_->apply(*D0x_,*residual_,Teuchos::NO_TRANS,one,one);
      }

      { // update current solution
        RCP<Teuchos::TimeMonitor> tmUpdate = getTimer("MueLu RefMaxwell: update");
        X.update(one, *residual_, one);
      }
    }

  }


  template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void RefMaxwell<Scalar,LocalOrdinal,GlobalOrdinal,Node>::solveH(const MultiVector& RHS, MultiVector& X) const {

    Scalar one = Teuchos::ScalarTraits<Scalar>::one();

    { // compute residual
      RCP<Teuchos::TimeMonitor> tmRes = getTimer("MueLu RefMaxwell: residual calculation");
      Utilities::Residual(*SM_Matrix_, X, RHS,*residual_);
      if (implicitTranspose_)
        P11_->apply(*residual_,*P11res_,Teuchos::TRANS);
      else
        R11_->apply(*residual_,*P11res_,Teuchos::NO_TRANS);
    }

    { // solve coarse (1,1) block
      if (!ImporterH_.is_null() && !implicitTranspose_) {
        RCP<Teuchos::TimeMonitor> tmH = getTimer("MueLu RefMaxwell: import coarse (1,1)");
        P11resTmp_->doImport(*P11res_, *ImporterH_, Xpetra::INSERT);
      }
      if (!AH_.is_null()) {
        RCP<Teuchos::TimeMonitor> tmH = getTimer("MueLu RefMaxwell: solve coarse (1,1)", AH_->getRowMap()->getComm());
        HierarchyH_->Iterate(*P11resSubComm_, *P11xSubComm_, numItersH_, true);
      }
    }

    { // update current solution
      RCP<Teuchos::TimeMonitor> tmUp = getTimer("MueLu RefMaxwell: update");
      P11_->apply(*P11x_,*residual_,Teuchos::NO_TRANS);
      X.update(one, *residual_, one);
    }

  }


  template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void RefMaxwell<Scalar,LocalOrdinal,GlobalOrdinal,Node>::solve22(const MultiVector& RHS, MultiVector& X) const {

    if (allNodesBoundary_)
      return;

    Scalar one = Teuchos::ScalarTraits<Scalar>::one();

    { // compute residual
      RCP<Teuchos::TimeMonitor> tmRes = getTimer("MueLu RefMaxwell: residual calculation");
      Utilities::Residual(*SM_Matrix_, X, RHS, *residual_);
      if (implicitTranspose_)
        D0_Matrix_->apply(*residual_,*D0res_,Teuchos::TRANS);
      else
        D0_T_Matrix_->apply(*residual_,*D0res_,Teuchos::NO_TRANS);
    }

    { // solve (2,2) block
      if (!Importer22_.is_null() && !implicitTranspose_) {
        RCP<Teuchos::TimeMonitor> tm22 = getTimer("MueLu RefMaxwell: import (2,2)");
        D0resTmp_->doImport(*D0res_, *Importer22_, Xpetra::INSERT);
      }
      if (!A22_.is_null()) {
        RCP<Teuchos::TimeMonitor> tm22 = getTimer("MueLu RefMaxwell: solve (2,2)", A22_->getRowMap()->getComm());
        Hierarchy22_->Iterate(*D0resSubComm_, *D0xSubComm_, numIters22_, true);
      }
    }

    { //update current solution
      RCP<Teuchos::TimeMonitor> tmUp = getTimer("MueLu RefMaxwell: update");
      D0_Matrix_->apply(*D0x_,*residual_,Teuchos::NO_TRANS);
      X.update(one, *residual_, one);
    }

  }


  template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void RefMaxwell<Scalar,LocalOrdinal,GlobalOrdinal,Node>::apply (const MultiVector& RHS, MultiVector& X,
                                                                  Teuchos::ETransp /* mode */,
                                                                  Scalar /* alpha */,
                                                                  Scalar /* beta */) const {

    RCP<Teuchos::TimeMonitor> tm = getTimer("MueLu RefMaxwell: solve");

    // make sure that we have enough temporary memory
    if (!allEdgesBoundary_ && X.getNumVectors() != P11res_->getNumVectors())
      allocateMemory(X.getNumVectors());

    { // apply pre-smoothing

      RCP<Teuchos::TimeMonitor> tmSm = getTimer("MueLu RefMaxwell: smoothing");

      PreSmoother_->Apply(X, RHS, use_as_preconditioner_);
    }

    // do solve for the 2x2 block system
    if(mode_=="additive")
      applyInverseAdditive(RHS,X);
    else if(mode_=="121") {
      solveH(RHS,X);
      solve22(RHS,X);
      solveH(RHS,X);
    } else if(mode_=="212") {
      solve22(RHS,X);
      solveH(RHS,X);
      solve22(RHS,X);
    } else if(mode_=="1")
      solveH(RHS,X);
    else if(mode_=="2")
      solve22(RHS,X);
    else if(mode_=="7") {
      solveH(RHS,X);
      { // apply pre-smoothing

        RCP<Teuchos::TimeMonitor> tmSm = getTimer("MueLu RefMaxwell: smoothing");

        PreSmoother_->Apply(X, RHS, false);
      }
      solve22(RHS,X);
      { // apply post-smoothing

        RCP<Teuchos::TimeMonitor> tmSm = getTimer("MueLu RefMaxwell: smoothing");

        PostSmoother_->Apply(X, RHS, false);
      }
      solveH(RHS,X);
    } else if(mode_=="none") {
      // do nothing
    }
    else
      applyInverseAdditive(RHS,X);

    { // apply post-smoothing

      RCP<Teuchos::TimeMonitor> tmSm = getTimer("MueLu RefMaxwell: smoothing");

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
             const Teuchos::RCP<Matrix> & Ms_Matrix,
             const Teuchos::RCP<Matrix> & M0inv_Matrix,
             const Teuchos::RCP<Matrix> & M1_Matrix,
             const Teuchos::RCP<MultiVector>  & Nullspace,
             const Teuchos::RCP<RealValuedMultiVector>  & Coords,
             Teuchos::ParameterList& List)
  {
    // some pre-conditions
    TEUCHOS_ASSERT(D0_Matrix!=Teuchos::null);
    TEUCHOS_ASSERT(Ms_Matrix!=Teuchos::null);
    TEUCHOS_ASSERT(M1_Matrix!=Teuchos::null);

#ifdef HAVE_MUELU_DEBUG
    TEUCHOS_ASSERT(M1_Matrix->getDomainMap()->isSameAs(*M1_Matrix->getRangeMap()));
    TEUCHOS_ASSERT(M1_Matrix->getDomainMap()->isSameAs(*M1_Matrix->getRowMap()));
    TEUCHOS_ASSERT(M1_Matrix->getDomainMap()->isSameAs(*D0_Matrix->getRangeMap()));

    TEUCHOS_ASSERT(Ms_Matrix->getDomainMap()->isSameAs(*Ms_Matrix->getRangeMap()));
    TEUCHOS_ASSERT(Ms_Matrix->getDomainMap()->isSameAs(*Ms_Matrix->getRowMap()));
    TEUCHOS_ASSERT(Ms_Matrix->getDomainMap()->isSameAs(*D0_Matrix->getRangeMap()));

    TEUCHOS_ASSERT(D0_Matrix->getRangeMap()->isSameAs(*D0_Matrix->getRowMap()));

    if (!M0inv_Matrix) {
      TEUCHOS_ASSERT(M0inv_Matrix->getDomainMap()->isSameAs(*M0inv_Matrix->getRowMap()));
      TEUCHOS_ASSERT(M0inv_Matrix->getDomainMap()->isSameAs(*M0inv_Matrix->getRangeMap()));
      TEUCHOS_ASSERT(M0inv_Matrix->getDomainMap()->isSameAs(*D0_Matrix->getDomainMap()));
    }
#endif

    HierarchyH_    = Teuchos::null;
    Hierarchy22_   = Teuchos::null;
    PreSmoother_   = Teuchos::null;
    PostSmoother_  = Teuchos::null;
    disable_addon_ = false;
    mode_          = "additive";

    // set parameters
    setParameters(List);

    if (D0_Matrix->getRowMap()->lib() == Xpetra::UseTpetra) {
      // We will remove boundary conditions from D0, and potentially change maps, so we copy the input.
      // Fortunately, D0 is quite sparse.
      // We cannot use the Tpetra copy constructor, since it does not copy the graph.

      RCP<Matrix> D0copy = MatrixFactory::Build(D0_Matrix->getRowMap(), D0_Matrix->getColMap(), 0);
      RCP<CrsMatrix> D0copyCrs = rcp_dynamic_cast<CrsMatrixWrap>(D0copy,true)->getCrsMatrix();
      ArrayRCP<const size_t> D0rowptr_RCP;
      ArrayRCP<const LO>     D0colind_RCP;
      ArrayRCP<const SC>     D0vals_RCP;
      rcp_dynamic_cast<CrsMatrixWrap>(D0_Matrix,true)->getCrsMatrix()->getAllValues(D0rowptr_RCP,
                                                                                    D0colind_RCP,
                                                                                    D0vals_RCP);

      ArrayRCP<size_t> D0copyrowptr_RCP;
      ArrayRCP<LO>     D0copycolind_RCP;
      ArrayRCP<SC>     D0copyvals_RCP;
      D0copyCrs->allocateAllValues(D0vals_RCP.size(),D0copyrowptr_RCP,D0copycolind_RCP,D0copyvals_RCP);
      D0copyrowptr_RCP.deepCopy(D0rowptr_RCP());
      D0copycolind_RCP.deepCopy(D0colind_RCP());
      D0copyvals_RCP.deepCopy(D0vals_RCP());
      D0copyCrs->setAllValues(D0copyrowptr_RCP,
                              D0copycolind_RCP,
                              D0copyvals_RCP);
      D0copyCrs->expertStaticFillComplete(D0_Matrix->getDomainMap(), D0_Matrix->getRangeMap(),
                                          rcp_dynamic_cast<CrsMatrixWrap>(D0_Matrix,true)->getCrsMatrix()->getCrsGraph()->getImporter(),
                                          rcp_dynamic_cast<CrsMatrixWrap>(D0_Matrix,true)->getCrsMatrix()->getCrsGraph()->getExporter());
      D0_Matrix_ = D0copy;
    } else
      D0_Matrix_ = MatrixFactory::BuildCopy(D0_Matrix);


    M0inv_Matrix_ = M0inv_Matrix;
    Ms_Matrix_    = Ms_Matrix;
    M1_Matrix_    = M1_Matrix;
    Coords_       = Coords;
    Nullspace_    = Nullspace;

    dump(*D0_Matrix_, "D0_clean.m");
    dump(*Ms_Matrix_, "Ms.m");
    dump(*M1_Matrix_, "M1.m");
    if (!M0inv_Matrix_.is_null()) dump(*M0inv_Matrix_, "M0inv.m");
    if (!Coords_.is_null())       dumpCoords(*Coords_, "coords.m");

  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void RefMaxwell<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
  describe(Teuchos::FancyOStream& out, const Teuchos::EVerbosityLevel /* verbLevel */) const {

    std::ostringstream oss;

    RCP<const Teuchos::Comm<int> > comm = SM_Matrix_->getDomainMap()->getComm();

#ifdef HAVE_MPI
    int root;
    if (!AH_.is_null())
      root = comm->getRank();
    else
      root = -1;

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

    {
      if (PreSmoother_ != null && PreSmoother_ == PostSmoother_)
        oss << "Smoother both : " << PreSmoother_->description() << std::endl;
      else {
        oss << "Smoother pre  : "
            << (PreSmoother_ != null ?  PreSmoother_->description() : "no smoother") << std::endl;
        oss << "Smoother post : "
            << (PostSmoother_ != null ?  PostSmoother_->description() : "no smoother") << std::endl;
      }
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

    if (!HierarchyH_.is_null())
      HierarchyH_->describe(out, GetVerbLevel());

    if (!Hierarchy22_.is_null())
      Hierarchy22_->describe(out, GetVerbLevel());

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
