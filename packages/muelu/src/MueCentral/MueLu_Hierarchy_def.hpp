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
#ifndef MUELU_HIERARCHY_DEF_HPP
#define MUELU_HIERARCHY_DEF_HPP

#include <sstream>
#include "MueLu_BoostGraphviz.hpp"

#include <Xpetra_MultiVectorFactory.hpp>
#include <Xpetra_Matrix.hpp>

#include "MueLu_Hierarchy_decl.hpp"
#include "MueLu_Level.hpp"
#include "MueLu_HierarchyHelpers.hpp"
#include "MueLu_SmootherFactoryBase.hpp"
#include "MueLu_SmootherFactory.hpp"
#include "MueLu_SmootherBase.hpp"
#include "MueLu_FactoryManager.hpp"
#include "MueLu_PFactory.hpp"
#include "MueLu_Utilities.hpp"
#include "MueLu_Monitor.hpp"

namespace MueLu {

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  Hierarchy<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::Hierarchy()
    : maxCoarseSize_(50), implicitTranspose_(false), isPreconditioner_(true), isDumpingEnabled_(false), dumpLevel_(-1),
      totalNnz_(0), operatorComplexity_(0.)
  {
    AddLevel(rcp( new Level() ));
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  Hierarchy<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::Hierarchy(const RCP<Matrix> & A)
    : maxCoarseSize_(50), implicitTranspose_(false), isPreconditioner_(true), isDumpingEnabled_(false), dumpLevel_(-1),
      totalNnz_(0), operatorComplexity_(0.)
  {
    RCP<Level> Finest = rcp( new Level() );
    AddLevel(Finest);

    Finest->Set("A", A);
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  Hierarchy<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::~Hierarchy() { }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void Hierarchy<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::SetMaxCoarseSize(Xpetra::global_size_t const &maxCoarseSize) { maxCoarseSize_ = maxCoarseSize; }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  Xpetra::global_size_t Hierarchy<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::GetMaxCoarseSize() const { return maxCoarseSize_; }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  int Hierarchy<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::LastLevelID() const { return Levels_.size() - 1; }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void Hierarchy<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::AddLevel(const RCP<Level> & level) {
    int levelID = LastLevelID() + 1; // ID of the inserted level

    if (level->GetLevelID() != -1 && (level->GetLevelID() != levelID))
      GetOStream(Warnings1, 0) << "Warning: Hierarchy::AddLevel(): Level with ID=" << level->GetLevelID() << " have been added at the end of the hierarchy" << std::endl
                               << "         but its ID have been redefined because last level ID of the hierarchy was " << LastLevelID() << "." << std::endl;

    Levels_.push_back(level);
    level->SetLevelID(levelID);

    if (levelID == 0)
      level->SetPreviousLevel(Teuchos::null);
    else
      level->SetPreviousLevel(Levels_[LastLevelID() - 1]);
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void Hierarchy<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::AddNewLevel() {
    RCP<Level> newLevel = Levels_[LastLevelID()]->Build(); // new coarse level, using copy constructor
    this->AddLevel(newLevel);                              // add to hierarchy
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  RCP<Level> & Hierarchy<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::GetLevel(const int levelID) {
    TEUCHOS_TEST_FOR_EXCEPTION(levelID < 0 || levelID > LastLevelID(), Exceptions::RuntimeError, "MueLu::Hierarchy::GetLevel(): invalid input parameter value: LevelID = " << levelID);
    return Levels_[levelID];
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  LocalOrdinal Hierarchy<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::GetNumLevels() const { return Levels_.size(); }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void Hierarchy<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::SetImplicitTranspose(const bool &implicit) { implicitTranspose_ = implicit; }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  bool Hierarchy<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::GetImplicitTranspose() const { return implicitTranspose_; }

  // Coherence checks todo in Setup() (using an helper function):
  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void Hierarchy<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::CheckLevel(Level& level, int levelID) {
    TEUCHOS_TEST_FOR_EXCEPTION(level.GetLevelID() != levelID, Exceptions::RuntimeError, "MueLu::Hierarchy::CheckLevel(): wrong level ID");
    TEUCHOS_TEST_FOR_EXCEPTION(levelID != 0 && level.GetPreviousLevel() != Levels_[levelID-1], Exceptions::RuntimeError, "MueLu::Hierarchy::Setup(): wrong level parent");
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  bool Hierarchy<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::Setup(int coarseLevelID, const Teuchos::Ptr<const FactoryManagerBase> fineLevelManager, const Teuchos::Ptr<const FactoryManagerBase> coarseLevelManager,
               const Teuchos::Ptr<const FactoryManagerBase> nextLevelManager) {

    // Use PrintMonitor/TimerMonitor instead of just a FactoryMonitor to print "Level 0" instead of Hierarchy(0)
    // Print is done after the requests for next coarse level
    TimeMonitor m1(*this, this->ShortClassName() + ": " + "Setup");
    TimeMonitor m2(*this, this->ShortClassName() + ": " + "Setup" + " (level=" + Teuchos::Utils::toString(coarseLevelID) + ")");

    TEUCHOS_TEST_FOR_EXCEPTION(coarseLevelManager == Teuchos::null, Exceptions::RuntimeError, "MueLu::Hierarchy::Setup(): argument coarseLevelManager cannot be null"); //So, it should not be passed as a pointer but as a reference

    //TODO
    typedef MueLu::TopRAPFactory<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps> TopRAPFactory;
    typedef MueLu::TopSmootherFactory<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps> TopSmootherFactory;

    //
    // Init
    //

    RCP<const FactoryManagerBase> rcpfineLevelManager  = rcpFromPtr(fineLevelManager);
    RCP<const FactoryManagerBase> rcpcoarseLevelManager= rcpFromPtr(coarseLevelManager);
    RCP<const FactoryManagerBase> rcpnextLevelManager  = rcpFromPtr(nextLevelManager);

    //    int coarseLevelID = LastLevelID() - 1; // Level built by this function
    TEUCHOS_TEST_FOR_EXCEPTION(LastLevelID() < coarseLevelID, Exceptions::RuntimeError, "MueLu::Hierarchy:Setup(): level " << coarseLevelID << " (specified by coarseLevelID argument) must be build before calling this function.");
    CheckLevel(*Levels_[coarseLevelID], coarseLevelID);

    bool isLastLevel = false;
    bool isFinestLevel = false;
    if(fineLevelManager == Teuchos::null) isFinestLevel = true;
    if(nextLevelManager == Teuchos::null) isLastLevel = true;

    // Attach FactoryManager to coarse and fine level
    SetFactoryManager SFMCoarse(Levels_[coarseLevelID], rcpcoarseLevelManager);
    RCP<SetFactoryManager> SFMFine, SFMNext;
    if (!isFinestLevel)
      SFMFine = rcp(new SetFactoryManager(Levels_[coarseLevelID-1], rcpfineLevelManager));

    //
    // Requests for finest level
    //

    if (isFinestLevel) {
      Levels_[coarseLevelID]->Request(TopSmootherFactory(rcpcoarseLevelManager, "Smoother")); // TODO: skip this line if we know that it is the lastLevel
      Levels_[coarseLevelID]->Request(TopSmootherFactory(rcpcoarseLevelManager, "CoarseSolver"));
    }

    if (isDumpingEnabled_ && dumpLevel_ == 0 && coarseLevelID == 1)
      DumpCurrentGraph();

    //
    // Requests for next coarse level
    //

    int nextLevelID = coarseLevelID + 1;

    if (!isLastLevel) {
      if (nextLevelID > LastLevelID()) { AddNewLevel(); }
      CheckLevel(*Levels_[nextLevelID], nextLevelID);
      SFMNext = rcp(new SetFactoryManager(Levels_[coarseLevelID+1], rcpnextLevelManager)); // Attach FactoryManager

      GetOStream(Debug, 0) << "Debug: Level: " << nextLevelID << " + R/S/C" << std::endl;
      Levels_[nextLevelID]->Request(TopRAPFactory(rcpcoarseLevelManager, rcpnextLevelManager));
      Levels_[nextLevelID]->Request(TopSmootherFactory(rcpnextLevelManager, "Smoother")); // TODO: skip this line if we know that it is the lastLevel
      Levels_[nextLevelID]->Request(TopSmootherFactory(rcpnextLevelManager, "CoarseSolver"));
    }

    PrintMonitor m0(*this, "Level " +  Teuchos::Utils::toString(coarseLevelID));

    //
    // Build coarse level
    //

    Level & level = *Levels_[coarseLevelID];

    // Build coarse level hierarchy
    if (!isFinestLevel) {
      TopRAPFactory coarseRAPFactory(rcpfineLevelManager, rcpcoarseLevelManager);
      coarseRAPFactory.Build(*level.GetPreviousLevel(), level);
      GetOStream(Debug, 0) << "Debug: Level: " << coarseLevelID << " - R" << std::endl;
      level.Release(coarseRAPFactory);
    }

    if (isDumpingEnabled_ && dumpLevel_ > 0 && coarseLevelID == dumpLevel_)
      DumpCurrentGraph();

    // Test if we reach the end of the hierarchy
    {
      RCP<Matrix> Ac;
      if (level.IsAvailable("A")) {
        Ac = level.Get<RCP<Matrix> >("A");
      } else {
        //TODO: happen when Ac factory = do nothing (ie: SetSmoothers)
      }

      if (Ac != Teuchos::null && Ac->getRowMap()->getGlobalNumElements() <= maxCoarseSize_) { // or if (coarseLevel == lastLevel
        if (isLastLevel == false) {
          GetOStream(Debug, 0) << "Debug: Level: " << nextLevelID << " - R/S/C" << std::endl;
          Levels_[nextLevelID]->Release(TopRAPFactory(rcpcoarseLevelManager, rcpnextLevelManager));
          Levels_[nextLevelID]->Release(TopSmootherFactory(rcpnextLevelManager, "Smoother"));
          Levels_[nextLevelID]->Release(TopSmootherFactory(rcpnextLevelManager, "CoarseSolver"));
          Levels_.pop_back(); // remove next level
        }

        isLastLevel = true;
      }
    }

    // Build coarse level smoother
    TopSmootherFactory smootherFact      (rcpcoarseLevelManager, "Smoother");
    TopSmootherFactory coarsestSolverFact(rcpcoarseLevelManager, "CoarseSolver");

    if (!isLastLevel) {
      smootherFact.Build(level);
    } else {
      RCP<Teuchos::FancyOStream> fos = Utils::MakeFancy(std::cout);
      level.print(*fos);
      coarsestSolverFact.Build(level); //TODO: PRE?POST
    }

    GetOStream(Debug, 0) << "Debug: Level: " << coarseLevelID << " - S/C" << std::endl;
    level.Release(smootherFact);
    level.Release(coarsestSolverFact);

    return isLastLevel;
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void Hierarchy<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::Setup(const FactoryManagerBase & manager, const int &startLevel, const int &numDesiredLevels) {
    PrintMonitor m0(*this, "Setup (" + this->MueLu::BaseClass::description() + ")"); // Use MueLu::BaseClass::description()to avoid printing "{numLevels = 1}" (numLevels is increasing...)

    RCP<const FactoryManagerBase> rcpManager = rcpFromRef(manager);

    // 2011/12 JG: Requests on the fine level are now posted at the beginning of the subroutine: Setup(fineLevelManager, coarseLevelManager, nextLevelManager)

    TEUCHOS_TEST_FOR_EXCEPTION(numDesiredLevels < 2, Exceptions::RuntimeError, "MueLu::Hierarchy::Setup(): numDesiredLevels < 2"); //FIXME: it has to work for one level method (numDesiredLevels==1)!!

    //TODO: check Levels_[startLevel] exists.

    // Check for fine level matrix A
    TEUCHOS_TEST_FOR_EXCEPTION(!Levels_[startLevel]->IsAvailable("A"), Exceptions::RuntimeError, "MueLu::Hierarchy::Setup(): no fine level matrix A! Set fine level matrix A using Level.Set()");

    // Check coarse levels
    // TODO: check if Ac available. If yes, issue a warning (bcse level already built...)

    //
    const int lastLevel = startLevel + numDesiredLevels - 1;
    int iLevel = 0;  // counter for the current number of multigrid levels after Setup phase
    GetOStream(Runtime0, 0) << "Loop: startLevel=" << startLevel << ", lastLevel=" << lastLevel << " (stop if numLevels = " << numDesiredLevels << " or Ac.size() = " << maxCoarseSize_ << ")" << std::endl;

    // set multigrid levels
    Teuchos::Ptr<const FactoryManagerBase> ptrmanager = Teuchos::ptrInArg(manager);
    bool bIsLastLevel = Setup(startLevel, Teuchos::null, ptrmanager, ptrmanager); // setup finest level (=level0)
    if(bIsLastLevel == false) {
      for(iLevel=startLevel + 1; iLevel < lastLevel; iLevel++) {                  // setup intermediate levels
        bIsLastLevel = Setup(iLevel, ptrmanager, ptrmanager, ptrmanager);
        if(bIsLastLevel == true) break;
      }
      if(bIsLastLevel == false) Setup(lastLevel, ptrmanager, ptrmanager, Teuchos::null); // setup coarsest level
    }

    // Levels_.resize(iLevel + 1);  // resize array of multigrid levels. add 1 to iLevel for the finest level (=level0) //TODO: still useful? Crop done in Setup() subroutine
    TEUCHOS_TEST_FOR_EXCEPTION(Levels_.size() != iLevel + 1, Exceptions::RuntimeError, "MueLu::Hierarchy::Setup(): number of level"); // some check like this should be done at the beginning of the routine

    // TODO: this is not exception safe: manager will still hold default factories if you exit this function with an exception
    manager.Clean();

  } // Setup()

  // ---------------------------------------- Summarize Hierarchy information -------------------------------

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  Teuchos::ParameterList Hierarchy<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::Summarize(MsgType verbLevel) {
     // TODO smoother information from each level
  
    Teuchos::ParameterList status;
    status.set("number of levels", Teuchos::as<int>(Levels_.size()));
    totalNnz_ = 0;
    for (int i=0; i<Levels_.size(); ++i) {
      TEUCHOS_TEST_FOR_EXCEPTION(!(Levels_[i]->IsAvailable("A")) , Exceptions::RuntimeError, "Operator complexity cannot be calculated because A is unavailable on level " << i);
      totalNnz_ += Levels_[i]->Get<RCP<Matrix> >("A")->getGlobalNumEntries();
    }
    operatorComplexity_ = Teuchos::as<double>(totalNnz_) / Levels_[0]->Get< RCP<Matrix> >("A")->getGlobalNumEntries();
    status.set("complexity", operatorComplexity_);

    GetOStream(verbLevel, 0) << "Number of levels    = " << status.get<int>("number of levels") << std::endl;
    GetOStream(verbLevel, 0) << "Operator complexity = " << std::setprecision(2) << std::setiosflags(std::ios::fixed) << status.get<double>("complexity") << std::endl;

    return status;
  } //Summarize()

  // ---------------------------------------- Iterate -------------------------------------------------------

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void Hierarchy<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::Iterate(MultiVector const &B, LO nIts, MultiVector &X, //TODO: move parameter nIts and default value = 1
                                                                                  const bool &InitialGuessIsZero, const CycleType &Cycle, const LO &startLevel)
  {

    RCP<Monitor> h;
    if (startLevel == 0)                                                               // -> Timing and msg only if startLevel == 0
      h = rcp(new Monitor(*this, "Iterate", (nIts == 1) ? None : Runtime0, Timings0)); // -> Do not issue msg if part of an iterative method (but always start timer)

    //Teuchos::Array<Magnitude> norms(1);
    bool zeroGuess=InitialGuessIsZero;

    RCP<Level> Fine = Levels_[startLevel];

    // Print residual information before iterating
    if (startLevel == 0 && IsPrint(Statistics1) && !isPreconditioner_) {

      Teuchos::Array<typename Teuchos::ScalarTraits<Scalar>::magnitudeType> rn;
      rn = Utils::ResidualNorm(*(Fine->Get< RCP<Matrix> >("A")), X, B);
      GetOStream(Statistics1, 0) << "iter:    "
                                 << std::setiosflags(std::ios::left)
                                 << std::setprecision(3) << 0 /* iter 0 */
                                 << "           residual = "
                                 << std::setprecision(10) << rn
                                 << std::endl;
    }

    for (LO i=1; i<=nIts; i++) {

      //X.norm2(norms);
      if (Fine->Get< RCP<Matrix> >("A")->getDomainMap()->isCompatible(*(X.getMap())) == false) {
        std::ostringstream buf;
        buf << startLevel;
        std::string msg = "Level " + buf.str() + ": level A's domain map is not compatible with X";
        throw(Exceptions::Incompatible(msg));
      }

      if (Fine->Get< RCP<Matrix> >("A")->getRangeMap()->isCompatible(*(B.getMap())) == false) {
        std::ostringstream buf;
        buf << startLevel;
        std::string msg = "Level " + buf.str() + ": level A's range map is not compatible with B";
        throw(Exceptions::Incompatible(msg));
      }

      //If on the coarse level, do either smoothing (if defined) or a direct solve.
      if (startLevel == ((LO)Levels_.size())-1) //FIXME is this right?
        {
          bool emptySolve = true;
          if (Fine->IsAvailable("PreSmoother")) { // important to use IsAvailable before Get here. It avoids building default smoother
            RCP<SmootherBase> preSmoo = Fine->Get< RCP<SmootherBase> >("PreSmoother");
            preSmoo->Apply(X, B, false);
            emptySolve=false;
          }
          if (Fine->IsAvailable("PostSmoother")) { // important to use IsAvailable before Get here. It avoids building default smoother
            RCP<SmootherBase> postSmoo = Fine->Get< RCP<SmootherBase> >("PostSmoother");
            postSmoo->Apply(X, B, false);
            emptySolve=false;
          }
          if (emptySolve==true)
            GetOStream(Warnings0, 0) << "Warning: No coarse grid solver" << std::endl;
        } else {
        //on an intermediate level
        RCP<Level> Coarse = Levels_[startLevel+1];

        //TODO: add IsAvailable test to avoid building default smoother
        if (Fine->IsAvailable("PreSmoother")) {
          RCP<SmootherBase> preSmoo = Fine->Get< RCP<SmootherBase> >("PreSmoother");
          preSmoo->Apply(X, B, zeroGuess);
        } else {
          GetOStream(Warnings0, 0) << "Warning: Level " <<  startLevel << ": No PreSmoother!" << std::endl;
        }

        RCP<MultiVector> residual = Utils::Residual(*(Fine->Get< RCP<Matrix> >("A")), X, B);

        RCP<Matrix> P = Coarse->Get< RCP<Matrix> >("P");
        RCP<Matrix> R;
        RCP<MultiVector> coarseRhs, coarseX;
        if (implicitTranspose_) {
          coarseRhs = MultiVectorFactory::Build(P->getDomainMap(), X.getNumVectors());
          coarseX   = MultiVectorFactory::Build(P->getDomainMap(), X.getNumVectors());
          P->apply(*residual, *coarseRhs, Teuchos::TRANS, 1.0, 0.0);
        } else {
          R = Coarse->Get< RCP<Matrix> >("R");
          coarseRhs = MultiVectorFactory::Build(R->getRangeMap(), X.getNumVectors());
          coarseX   = MultiVectorFactory::Build(R->getRangeMap(), X.getNumVectors());
          R->apply(*residual, *coarseRhs, Teuchos::NO_TRANS, 1.0, 0.0);
        }
        coarseX->putScalar(0.);

        Iterate(*coarseRhs, 1, *coarseX, true, Cycle, startLevel+1);
        // ^^ zero initial guess
        if (Cycle>1)
          Iterate(*coarseRhs, 1, *coarseX, false, Cycle, startLevel+1);
        // ^^ nonzero initial guess

        // update X+=P * coarseX
        //P->apply(*coarseX, X, Teuchos::NO_TRANS, 1.0, 1.0);  //Xpetra throws an error if linAlgebra==Epetra
        RCP<MultiVector> correction = MultiVectorFactory::Build(P->getRangeMap(), X.getNumVectors());
        P->apply(*coarseX, *correction, Teuchos::NO_TRANS, 1.0, 0.0);
        X.update(1.0, *correction, 1.0);

        //X.norm2(norms);
        //TODO: add IsAvailable test to avoid building default smoother
        if (Fine->IsAvailable("PostSmoother")) {
          RCP<SmootherBase> postSmoo = Fine->Get< RCP<SmootherBase> >("PostSmoother");
          postSmoo->Apply(X, B, false);
        } else {
          GetOStream(Warnings0, 0) << "Warning: Level " <<  startLevel << ": No PostSmoother!" << std::endl;
        }
      }
      zeroGuess=false;

      if (startLevel == 0 && IsPrint(Statistics1) && !isPreconditioner_) {
        Teuchos::Array<typename Teuchos::ScalarTraits<Scalar>::magnitudeType> rn;
        rn = Utils::ResidualNorm(*(Fine->Get< RCP<Matrix> >("A")), X, B);
        GetOStream(Statistics1, 0) << "iter:    "
                                   << std::setiosflags(std::ios::left)
                                   << std::setprecision(3) << i
                                   << "           residual = "
                                   << std::setprecision(10) << rn
                                   << std::endl;
      }

    } //for (LO i=0; i<nIts; i++)

  } //Iterate()

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void Hierarchy<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::Write(const LO &start, const LO &end)
  {

    LO startLevel = start;
    LO endLevel = end;

    if (startLevel==-1) startLevel = 0;
    if (endLevel==-1)   endLevel = Levels_.size();

    TEUCHOS_TEST_FOR_EXCEPTION(startLevel > endLevel, Exceptions::RuntimeError, "MueLu::Hierarchy::Write : startLevel must be <= endLevel");

    TEUCHOS_TEST_FOR_EXCEPTION(startLevel < 0 || endLevel > Levels_.size(), Exceptions::RuntimeError, "MueLu::Hierarchy::Write bad start or end level");

    for (LO i=startLevel; i<endLevel; ++i) {

      std::ostringstream buf; buf << i;
      std::string fileName = "A_" + buf.str() + ".m";
      Utils::Write( fileName,*(Levels_[i]-> template Get< RCP< Matrix> >("A")) );
      if (i>startLevel) {
        fileName = "P_" + buf.str() + ".m";
        Utils::Write( fileName,*(Levels_[i]-> template Get< RCP< Matrix> >("P")) );
        if (!implicitTranspose_) {
          fileName = "R_" + buf.str() + ".m";
          Utils::Write( fileName,*(Levels_[i]-> template Get< RCP< Matrix> >("R")) );
        }
      }
    }

  } //Write()

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void Hierarchy<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::Keep(const std::string & ename, const FactoryBase* factory) {
    for (Array<RCP<Level> >::iterator it = Levels_.begin(); it != Levels_.end(); ++it)
      (*it)->Keep(ename, factory);
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void Hierarchy<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::Delete(const std::string& ename, const FactoryBase* factory) {
    for (Array<RCP<Level> >::iterator it = Levels_.begin(); it != Levels_.end(); ++it)
      (*it)->Delete(ename, factory);
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void Hierarchy<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::AddKeepFlag(const std::string & ename, const FactoryBase* factory, KeepType keep) {
    for (Array<RCP<Level> >::iterator it = Levels_.begin(); it != Levels_.end(); ++it)
      (*it)->AddKeepFlag(ename, factory, keep);
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void Hierarchy<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::RemoveKeepFlag(const std::string & ename, const FactoryBase* factory, KeepType keep) {
    for (Array<RCP<Level> >::iterator it = Levels_.begin(); it != Levels_.end(); ++it)
      (*it)->RemoveKeepFlag(ename, factory, keep);
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  std::string Hierarchy<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::description() const {
    std::ostringstream out;
    out << BaseClass::description();
    out << "{numLevels = " << GetNumLevels() << "}";
    return out.str();
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void Hierarchy<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::print(Teuchos::FancyOStream &out, const VerbLevel verbLevel) const {
    MUELU_DESCRIBE;

    if (verbLevel & Parameters0) {
      out0 << "Number of Levels: " << GetNumLevels() << endl;
    }

    if (verbLevel & Parameters1) {
      out0 << "Max Coarse Size: "    << maxCoarseSize_ << endl;
      out0 << "Implicit Transpose: " << implicitTranspose_ << endl;
    }

    if (verbLevel & Statistics1) {
      Teuchos::OSTab tab2(out);
      for(int i = 0; i < GetNumLevels(); i++) {
        Levels_[i]->print(out0, verbLevel);
        out0 << std::endl;
      }
    }
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  Hierarchy<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::Hierarchy(const Hierarchy &h) { }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void Hierarchy<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::IsPreconditioner(const bool flag) {
    isPreconditioner_ = flag;
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void Hierarchy<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::DumpCurrentGraph() const {
#if defined(HAVE_MUELU_BOOST) && defined(BOOST_VERSION) && (BOOST_VERSION >= 104400)
    BoostGraph      graph;

    BoostProperties dp;
    dp.property("label", boost::get(boost::vertex_name,  graph));
    dp.property("id",    boost::get(boost::vertex_index, graph));
    dp.property("label", boost::get(boost::edge_name,    graph));
    dp.property("color", boost::get(boost::edge_color,   graph));

    // create local maps
    std::map<const FactoryBase*, BoostVertex>                                     vindices;
    typedef std::map<std::pair<BoostVertex,BoostVertex>, std::string> emap; emap  edges;

    for (int i = dumpLevel_; i <= dumpLevel_+1 && i < GetNumLevels(); i++) {
      edges.clear();
      Levels_[i]->UpdateGraph(vindices, edges, dp, graph);

      for (emap::const_iterator eit = edges.begin(); eit != edges.end(); eit++) {
        std::pair<BoostEdge, bool> boost_edge = boost::add_edge(eit->first.first, eit->first.second, graph);
        boost::put("label", dp, boost_edge.first, eit->second);
        if (i == dumpLevel_)
          boost::put("color", dp, boost_edge.first, std::string("red"));
        else
          boost::put("color", dp, boost_edge.first, std::string("blue"));
      }
    }

    // add legend
    std::ostringstream legend;
    legend << "< <TABLE BORDER=\"0\" CELLBORDER=\"1\" CELLSPACING=\"0\" CELLPADDING=\"4\"> \
               <TR><TD COLSPAN=\"2\">Legend</TD></TR> \
               <TR><TD><FONT color=\"red\">Level " << dumpLevel_ << "</FONT></TD><TD><FONT color=\"blue\">Level " << dumpLevel_+1 << "</FONT></TD></TR> \
               </TABLE> >";
    BoostVertex boost_vertex = boost::add_vertex(graph);
    boost::put("label", dp, boost_vertex, legend.str());

    std::ofstream out(dumpFile_.c_str());
    boost::write_graphviz_dp(out, graph, dp, std::string("id"));
#endif
  }

} //namespace MueLu

// TODO: We need a Set/Get function to change the CycleType (for when Iterate() calls are embedded in a Belos Preconditionner for instance).

#endif // MUELU_HIERARCHY_DEF_HPP
