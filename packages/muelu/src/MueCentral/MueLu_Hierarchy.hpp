#ifndef MUELU_HIERARCHY_HPP
#define MUELU_HIERARCHY_HPP

#include "MueLu_ConfigDefs.hpp"
#include "MueLu_BaseClass.hpp"
#include "MueLu_FactoryManagerBase.hpp"
#include "MueLu_Types.hpp"
#include "MueLu_Level.hpp"
#include "MueLu_Monitor.hpp"

#include "MueLu_FactoryManager.hpp"

#include "MueLu_Exceptions.hpp"
#include "MueLu_SmootherFactory.hpp"

#include "MueLu_HierarchyHelpers.hpp"

namespace MueLu {
/*!
    @class Hierarchy
    @brief Provides methods to build a multigrid hierarchy and apply multigrid cycles.

    Allows users to manually populate operators at different levels within 
    a multigrid method and push them into the hierarchy via SetLevel() 
    and/or to supply factories for automatically generating prolongators, 
    restrictors, and coarse level discretizations.  Additionally, this class contains 
    an apply method that supports V and W cycles.
 */
template <class Scalar = double, class LocalOrdinal = int, class GlobalOrdinal = LocalOrdinal, class Node = Kokkos::DefaultNode::DefaultNodeType, class LocalMatOps = typename Kokkos::DefaultKernels<void,LocalOrdinal,Node>::SparseOps> //TODO: or BlockSparseOp ?
class Hierarchy : public BaseClass {

#include "MueLu_UseShortNames.hpp"

public:

  //! @name Constructors/Destructors
  //@{

  //! Default constructor.
  Hierarchy()
  : maxCoarseSize_(50), implicitTranspose_(false)
  {
    AddLevel(rcp( new Level() ));
  }

  //! Constructor
  Hierarchy(const RCP<Operator> & A)
  :  maxCoarseSize_(50), implicitTranspose_(false)
  {
    RCP<Level> Finest = rcp( new Level() );
    AddLevel(Finest);

    Finest->Set("A", A);
  }

  //! Destructor.
  virtual ~Hierarchy() {}

  //@}

  //! @name Set/Get Methods.
  //@{

  //!
  void SetMaxCoarseSize(Xpetra::global_size_t const &maxCoarseSize) { maxCoarseSize_ = maxCoarseSize; }

  //!
  Xpetra::global_size_t GetMaxCoarseSize() const { return maxCoarseSize_; }

private:
  int LastLevelID() const { return Levels_.size() - 1; }

public:

  //! Add a level at the end of the hierarchy
  void AddLevel(const RCP<Level> & level) {
    int levelID = LastLevelID() + 1; // ID of the inserted level

    if (level->GetLevelID() != -1 && (level->GetLevelID() != levelID))
      GetOStream(Warnings1, 0) << "Warning: Hierarchy::AddLevel(): Level with ID=" << level->GetLevelID() << " have been added at the end of the hierarchy" << endl
			       << "         but its ID have been redefined because last level ID of the hierarchy was " << LastLevelID() << "." << std::endl;

    Levels_.push_back(level);
    level->SetLevelID(levelID);

    if (levelID == 0)
      level->SetPreviousLevel(Teuchos::null);
    else
      level->SetPreviousLevel(Levels_[LastLevelID() - 1]);
  }

  //! Add a new level at the end of the hierarchy
  void AddNewLevel() {
    RCP<Level> newLevel = Levels_[LastLevelID()]->Build(); // new coarse level, using copy constructor
    this->AddLevel(newLevel);                              // add to hierarchy
  }


  //! Retrieve a certain level from hierarchy.
  RCP<Level> & GetLevel(const int levelID = 0) {
    TEST_FOR_EXCEPTION(levelID < 0 || levelID > LastLevelID(), Exceptions::RuntimeError, "MueLu::Hierarchy::GetLevel(): invalid input parameter value: LevelID = " << levelID);
    return Levels_[levelID];
  }

  LO GetNumLevels() const { return Levels_.size(); }

  //! Indicate that Iterate should use tranpose of prolongator for restriction operations.
  void SetImplicitTranspose(const bool &implicit) { implicitTranspose_ = implicit; }

  //! If true is returned, iterate will use tranpose of prolongator for restriction operations.
  bool GetImplicitTranspose() const { return implicitTranspose_; }

  //@}

  //! @name Populate Methods.
  //@{

  /*!
      @brief Constructs components of the hierarchy.

      Invoke a set of factories to populate (construct prolongation,
      restriction, coarse level discretizations, and smoothers in this
      order) a multigrid Hierarchy starting with information on 'startLevel'
      and continuing for at most 'numDesiredLevels'.
   */
  Teuchos::ParameterList FullPopulate(const FactoryBase & PFact,
      const FactoryBase & RFact,
      const TwoLevelFactoryBase & AcFact,
      const SmootherFactory & SmooFact,
      const int &startLevel = 0, const int &numDesiredLevels = 10) {

    // Note: It's OK to use rcpFromRef here, because data will only be kept by the FactoryManager
    //       and the FactoryManager is deleted at the end of this function.

    FactoryManager manager(rcpFromRef(PFact), rcpFromRef(RFact), rcpFromRef(AcFact));
    manager.SetFactory("Smoother", rcpFromRef(SmooFact));

    return Setup(manager, startLevel, numDesiredLevels);
  }

  /*! @brief Populate hierarchy with A's, R's, and P's.

    Invoke a set of factories to populate (construct prolongation,
    restriction, and coarse level discretizations in this
    order) a multigrid Hierarchy starting with information on 'startLevel' 
    and continuing for at most 'numDesiredLevels'. 

    @return  List containing starting and ending level numbers, operator complexity, \#nonzeros in the fine
    matrix, and the sum of nonzeros all matrices (including the fine).
   */
  Teuchos::ParameterList FillHierarchy(const PFactory & PFact, const RFactory & RFact,
      const TwoLevelFactoryBase & AcFact,
      const int startLevel = 0, const int numDesiredLevels = 10) {

    FactoryManager manager(rcpFromRef(PFact), rcpFromRef(RFact), rcpFromRef(AcFact));

    return Setup(manager, startLevel, numDesiredLevels);

  }
  // FillHierarchy

  Teuchos::ParameterList Setup(const FactoryManager & manager, const int &startLevel = 0, const int &numDesiredLevels = 10) {
    RCP<const FactoryManager> rcpManager = rcpFromRef(manager);

    TopRAPFactory<SC,LO,GO,NO>      rapFactory           (rcpManager); //TODO: remove SC,LO,GO,NO
    TopSmootherFactory<SC,LO,GO,NO> smootherFactory      (rcpManager, "Smoother");
    TopSmootherFactory<SC,LO,GO,NO> coarsestSolverFactory(rcpManager, "CoarseSolver");

    Monitor h(*this, "Setup");
    Xpetra::global_size_t sumCoarseNnz = 0;

    TEST_FOR_EXCEPTION(numDesiredLevels < 2, Exceptions::RuntimeError, "MueLu::Hierarchy::Setup(): numDesiredLevels < 2");

    //TODO: check Levels_[startLevel] exists.

    // Check for fine level matrix A
    TEST_FOR_EXCEPTION(!Levels_[startLevel]->IsAvailable("A"), Exceptions::RuntimeError, "MueLu::Hierarchy::Setup(): no fine level matrix A! Set fine level matrix A using Level.Set()");

    // Check coarse levels
    // TODO: check if Ac available. If yes, issue a warning (bcse level already built...)

    Levels_[startLevel]->Request(smootherFactory);
    // Levels_[startLevel]->Request(coarsestSolverFactory);

    //
    const int lastLevel = startLevel + numDesiredLevels - 1;
    int iLevel;
    GetOStream(Runtime0, 0) << "Loop: startLevel=" << startLevel << ", lastLevel=" << lastLevel << " (numLevels=" << numDesiredLevels << ")" << std::endl;

    for (iLevel = startLevel; iLevel <= lastLevel; iLevel++) {
      SubMonitor m(*this, "Level " + Teuchos::toString(iLevel));

      // Prepare next level (to keep useful data of current level!)
      int nextCoarseLevelID = iLevel + 1;
      if (nextCoarseLevelID <= lastLevel) {
        if (nextCoarseLevelID > LastLevelID()) { AddNewLevel(); }

        std::cout << "Level " << nextCoarseLevelID << ": Request RAP" << std::endl; std::cout.flush();
        Levels_[nextCoarseLevelID]->Request(rapFactory);
        Levels_[nextCoarseLevelID]->Request(smootherFactory); //TODO: skip if lastLevel
        // Levels_[nextCoarseLevelID]->Request(coarsestSolverFactory);
      }

      Level & level = *Levels_[iLevel];
      TEST_FOR_EXCEPTION(level.GetLevelID() != iLevel, Exceptions::RuntimeError, "MueLu::Hierarchy::Setup(): wrong level ID");

      //
      if (iLevel != startLevel) {
        TEST_FOR_EXCEPTION(level.GetPreviousLevel() != Levels_[iLevel-1], Exceptions::RuntimeError, "MueLu::Hierarchy::Setup(): wrong level parent");

        rapFactory.Build(*level.GetPreviousLevel(), level);
        std::cout << "Level " << iLevel << ": Release RAP" << std::endl; std::cout.flush();
        level.Release(rapFactory);
      }
      //

      RCP<Operator> Ac;
      if (level.IsAvailable("A")) {
        Ac = level.Get<RCP<Operator> >("A");
        sumCoarseNnz += Ac->getGlobalNumEntries();
      } else {
        //TODO: happen when Ac factory = do nothing (ie: SetSmoothers)
      }

      if (iLevel == lastLevel || (Ac != Teuchos::null && Ac->getRowMap()->getGlobalNumElements() <= maxCoarseSize_)) {
        if (nextCoarseLevelID <= lastLevel) {
          std::cout << "Level " << nextCoarseLevelID << ": Release RAP" << std::endl; std::cout.flush();
          // Levels_[nextCoarseLevelID]->Release(rapFactory);
          // Levels_[nextCoarseLevelID]->Release(smootherFactory);
          // Levels_[nextCoarseLevelID]->Release(coarsestSolverFactory);
        }

        std::cout << "BUILD COARSE" << std::endl; std::cout.flush();
        coarsestSolverFactory.Build(level); //TODO: PRE?POST
        // level.Release(smootherFactory);
        level.Release(coarsestSolverFactory);

        break;

      } else {

        std::cout << "BUILD SMOO" << std::endl; std::cout.flush();
        smootherFactory.Build(level);
        level.Release(smootherFactory);
        // level.Release(coarsestSolverFactory);

      }

    } // for loop

    // Crop. TODO: add a warning
    Levels_.resize(iLevel + 1);

    // TODO: not exception safe: manager will still hold default factories if you exit this function with an exception
    manager.Clean();

    // Gather statistics
    Xpetra::global_size_t fineNnz = Levels_[startLevel]->Get< RCP<Operator> >("A")->getGlobalNumEntries();
    Xpetra::global_size_t totalNnz = fineNnz + sumCoarseNnz;

    Teuchos::ParameterList status;
    status.set("fine nnz", fineNnz);
    status.set("total nnz", totalNnz);
    status.set("start level", startLevel);
    status.set("end level", iLevel); //TODO: check if it's actually correct (exit by break vs. end of loop)
    status.set("operator complexity", ((SC)totalNnz) / fineNnz);

    return status;

  } // Setup()

  /*! @brief Set solve method for coarsest level.

    @param smooFact  fully constructed SmootherFactory 
    @param pop       whether to use pre,post, or both pre and post smoothing 

    Note: Whether the SmootherFactory builds both a pre- and post-smoother can be also be
    controlled by SmootherFactory::SetSmootherPrototypes. This approach is a bit cumbersome,
    however.
   */
  //TODO: remove PRE/POST

  void SetCoarsestSolver(SmootherFactoryBase const &smooFact, PreOrPost const &pop = BOTH) {
    Level & level = *Levels_[LastLevelID()];
    RCP<const FactoryManagerBase> manager = rcp(new FactoryManager());
    SetFactoryManager SFM(level, manager);

    level.Request(smooFact);
    smooFact.BuildSmoother(level, pop);

    if (level.IsAvailable("PreSmoother", &smooFact)) {
      RCP<SmootherBase> Pre  = level.Get<RCP<SmootherBase> >("PreSmoother", &smooFact);
      level.Set("PreSmoother", Pre);
    }

    if (level.IsAvailable("PostSmoother", &smooFact)) {
      RCP<SmootherBase> Post = level.Get<RCP<SmootherBase> >("PostSmoother", &smooFact);
      level.Set("PostSmoother", Post);
    }

    level.Release(smooFact);
  }

  /*! @brief Construct smoothers on all levels but the coarsest.

      Invoke a set of factories to construct smoothers within 
      a multigrid Hierarchy starting with information on 'startLevel' 
      and continuing for at most 'numDesiredLevels'. 

      Note: last level smoother will not be set here. Use SetCoarsestSolver()
      to define a smoother for the last level. Otherwise, a direct solve is
      assumed
   */
  void SetSmoothers(SmootherFactory const & smooFact, LO const & startLevel = 0, LO numDesiredLevels = -1) {
    Monitor h(*this, "SetSmoothers");

    if (numDesiredLevels == -1)
      numDesiredLevels = GetNumLevels() - startLevel - 1;
    LO lastLevel = startLevel + numDesiredLevels - 1;

    //checks
    if (startLevel >= GetNumLevels())
      throw(Exceptions::RuntimeError("startLevel >= actual number of levels"));

    if (startLevel == GetNumLevels() - 1)
      throw(Exceptions::RuntimeError("startLevel == coarse level. Use SetCoarseSolver()"));

    if (lastLevel >= GetNumLevels() - 1) {
      lastLevel = GetNumLevels() - 2;
      GetOStream(Warnings0, 0) << "Warning: coarsest level solver will not be changed!" << std::endl;
    }

    FactoryManager manager;
    manager.SetFactory("Smoother", rcpFromRef(smooFact));
    manager.SetFactory("P", Teuchos::null);
    manager.SetFactory("R", Teuchos::null);
    manager.SetFactory("A", Teuchos::null);

    lastLevel++; // hack: nothing will be done on the last level in Setup() because coarse solver of manager == Teuchos::null. TODO: print() of Setup() will be confusing
    numDesiredLevels = lastLevel - startLevel + 1;

    // std::cout << "startLevel=" << startLevel << ", nummDesiredLevels=" << numDesiredLevels << std::endl;
    Setup(manager, startLevel, numDesiredLevels);

  } //SetSmoothers()

  // #define GimmeNorm(someVec,someLabel) { (someVec).norm2(norms); GetOStream(Statistics1, 0) << someLabel << " = " << norms << std::endl; }

  /*!
      @brief Apply the multigrid preconditioner.

      In theory, more general cycle types than just V- and W-cycles are possible.  However,
      the enumerated type CycleType would have to be extended.

      @param B right-hand side of linear problem
      @param nIts number of multigrid cycles to perform
      @param X initial and final (approximate) solution of linear problem
      @param InitialGuessIsZero Indicates whether the initial guess is zero
      @param Cycle Supports VCYCLE and WCYCLE types.
   */
  void Iterate(MultiVector const &B, LO nIts, MultiVector &X, //TODO: move parameter nIts and default value = 1
      const bool &InitialGuessIsZero = false, const CycleType &Cycle = VCYCLE, const LO &startLevel = 0)
  {
    //Teuchos::Array<Magnitude> norms(1);
    bool zeroGuess=InitialGuessIsZero;

    for (LO i=0; i<nIts; i++) {

      RCP<Level> Fine = Levels_[startLevel];

      if (startLevel == 0 && IsPrint(Statistics1)) {
        GetOStream(Statistics1, 0) << "iter:    "
				   << std::setiosflags(std::ios::left)
				   << std::setprecision(3) << i
				   << "           residual = "
				   << std::setprecision(10) << Utils::ResidualNorm(*(Fine->Get< RCP<Operator> >("A")),X,B)
				   << std::endl;
      }

      //X.norm2(norms);
      if (Fine->Get< RCP<Operator> >("A")->getDomainMap()->isCompatible(*(X.getMap())) == false) {
        std::ostringstream buf;
        buf << startLevel;
        std::string msg = "Level " + buf.str() + ": level A's domain map is not compatible with X";
        throw(Exceptions::Incompatible(msg));
      }

      if (Fine->Get< RCP<Operator> >("A")->getRangeMap()->isCompatible(*(B.getMap())) == false) {
        std::ostringstream buf;
        buf << startLevel;
        std::string msg = "Level " + buf.str() + ": level A's range map is not compatible with B";
        throw(Exceptions::Incompatible(msg));
      }

      //If on the coarse level, do either smoothing (if defined) or a direct solve.
      if (startLevel == ((LO)Levels_.size())-1) //FIXME is this right?
      {
        bool emptySolve = true;
        if (Fine->IsAvailable("PreSmoother")) { // important do use IsAvailable before Get here. Avoid building default smoother
          RCP<SmootherBase> preSmoo = Fine->Get< RCP<SmootherBase> >("PreSmoother");
          preSmoo->Apply(X, B, false);
          emptySolve=false;
        }
        if (Fine->IsAvailable("PostSmoother")) { // important do use IsAvailable before Get here. Avoid building default smoother
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
          GetOStream(Errors, 0) << "Error: Level " <<  startLevel << ": No Smoother!" << std::endl;
        }

        RCP<MultiVector> residual = Utils::Residual(*(Fine->Get< RCP<Operator> >("A")),X,B);

        RCP<Operator> P = Coarse->Get< RCP<Operator> >("P");
        RCP<Operator> R;
        RCP<MultiVector> coarseRhs, coarseX;
        if (implicitTranspose_) {
          coarseRhs = MultiVectorFactory::Build(P->getDomainMap(),X.getNumVectors());
          coarseX   = MultiVectorFactory::Build(P->getDomainMap(),X.getNumVectors());
          P->apply(*residual,*coarseRhs,Teuchos::TRANS,1.0,0.0);
        } else {
          R = Coarse->Get< RCP<Operator> >("R");
          coarseRhs = MultiVectorFactory::Build(R->getRangeMap(),X.getNumVectors());
          coarseX   = MultiVectorFactory::Build(R->getRangeMap(),X.getNumVectors());
          R->apply(*residual,*coarseRhs,Teuchos::NO_TRANS,1.0,0.0);
        }
        coarseX->putScalar(0.);

        Iterate(*coarseRhs,1,*coarseX,true,Cycle,startLevel+1);
        // ^^ zero initial guess
        if (Cycle>1)
          Iterate(*coarseRhs,1,*coarseX,false,Cycle,startLevel+1);
        // ^^ nonzero initial guess

        // update X+=P * coarseX
        //P->apply(*coarseX,X,Teuchos::NO_TRANS,1.0,1.0);  //Xpetra throws an error if linAlgebra==0
        RCP<MultiVector> correction = MultiVectorFactory::Build(P->getRangeMap(),X.getNumVectors());
        P->apply(*coarseX,*correction,Teuchos::NO_TRANS,1.0,0.0);
        X.update(1.0,*correction,1.0);

        //X.norm2(norms);
        //TODO: add IsAvailable test to avoid building default smoother
        if (Fine->IsAvailable("PostSmoother")) {
          RCP<SmootherBase> postSmoo = Fine->Get< RCP<SmootherBase> >("PostSmoother");
          postSmoo->Apply(X, B, false);
        } else {
          GetOStream(Errors,0) << "Error: Level " <<  startLevel << ": No Smoother!" << std::endl;
        }
      }
      zeroGuess=false;
    } //for (LO i=0; i<nIts; i++)

  } //Iterate()

  //@}

private:
  //! Copy constructor is not implemented.
  Hierarchy(const Hierarchy &h) { }

  //! vector of Level objects
  Array<RCP<Level> > Levels_;

  Xpetra::global_size_t maxCoarseSize_;
  bool implicitTranspose_;

}; //class Hierarchy

} //namespace MueLu

#define MUELU_HIERARCHY_SHORT

#endif //ifndef MUELU_HIERARCHY_HPP

// TODO: We need a Set/Get function to change the CycleType (for when Iterate() calls are embedded in a Belos Preconditionner for instance).
