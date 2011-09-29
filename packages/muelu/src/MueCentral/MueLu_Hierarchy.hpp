#ifndef MUELU_HIERARCHY_HPP
#define MUELU_HIERARCHY_HPP

#include "MueLu_ConfigDefs.hpp"
#include "MueLu_BaseClass.hpp"
#include "MueLu_DefaultFactoryHandler.hpp"
#include "MueLu_Types.hpp"
#include "MueLu_Level.hpp"
#include "MueLu_Monitor.hpp"

#include "MueLu_Exceptions.hpp"
#include "MueLu_SmootherFactory.hpp"

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

    template<class AA, class BB, class CC, class DD, class EE>
    inline friend std::ostream& operator<<(std::ostream& os, Hierarchy<AA,BB,CC,DD,EE> &hierarchy);

  private:

    //! vector of Level objects
    std::vector<RCP<Level> > Levels_;

    Xpetra::global_size_t maxCoarseSize_;
    bool implicitTranspose_;

  public:

    //! @name Constructors/Destructors
    //@{

    //! Default constructor.
    Hierarchy() 
      : maxCoarseSize_(50), implicitTranspose_(false)
    { }

    //!
    Hierarchy(const RCP<Operator> & A) :  maxCoarseSize_(50), implicitTranspose_(false) {
      RCP<Level> Finest = rcp( new Level() );
      Finest->Set< RCP<Operator> >("A", A);
      SetLevel(Finest);
    }

    //! Copy constructor.
    Hierarchy(Hierarchy const &inHierarchy) {
      std::cerr << "Not implemented yet." << std::endl;
    }

    //! Destructor.
    virtual ~Hierarchy() {}

    //@}

    //! @name Set/Get Methods.
    //@{

    //! @name Set methods.
    //@{
    void SetMaxCoarseSize(Xpetra::global_size_t const &maxCoarseSize) {
      maxCoarseSize_ = maxCoarseSize;
    }

    //@}

    //! @name Get methods.
    //@{
  
    Xpetra::global_size_t GetMaxCoarseSize() const {
      return maxCoarseSize_;
    }

    //! Assign a level to hierarchy.
    // TODO from JG: behavior should be explain or changed. 
    //               Right now, we change the LevelID of the input level and push it at the end of hierarchy.
    //               Certainly better to take the LevelID of the input level into account
    void SetLevel(RCP<Level> const& level) {
      Levels_.push_back(level);
      level->SetLevelID(Levels_.size());

      if (Levels_.size() < 2)
        level->SetPreviousLevel(Teuchos::null);
      else
        level->SetPreviousLevel(Levels_[Levels_.size()-2]); // new level = size-1, previous = size-2
    }

    //! Retrieve a certain level from hierarchy.
    RCP<Level>& GetLevel(int const levelID) /* const */ {
      TEST_FOR_EXCEPTION(levelID < 1, Exceptions::RuntimeError, "MueLu::Hierarchy::GetLevel(): invalid input parameter value");
      return Levels_[levelID-1];
    }

    LO GetNumberOfLevels() {
      return Levels_.size();
    }

    //! Indicate that Iterate should use tranpose of prolongator for restriction operations.
    void SetImplicitTranspose(bool const &implicit) {
      implicitTranspose_ = implicit;
    }

    //! If true is returned, iterate will use tranpose of prolongator for restriction operations.
    bool GetImplicitTranspose() {
      return implicitTranspose_;
    }

    //@}

    //! @name Populate Methods.
    //@{

    /*!
      @brief Constructs components of the hierarchy.

      Invoke a set of factories to populate (construct prolongation,
      restriction, coarse level discretizations, and smoothers in this
      order) a multigrid Hierarchy starting with information on 'startLevel'
      and continuing for at most 'numDesiredLevels'.

      Note: Empty factories are simply skipped.
    */
    Teuchos::ParameterList FullPopulate(const PFactory & PFact,
                                        const RFactory & RFact,
                                        const TwoLevelFactoryBase & AcFact,
                                        const SmootherFactory & SmooFact,
                                        int const &startLevel=0, int const &numDesiredLevels=10 )
    {
      Teuchos::ParameterList status;
      status = FillHierarchy(PFact,RFact,AcFact,startLevel,numDesiredLevels);

      SetSmoothers(SmooFact,startLevel,numDesiredLevels-1);
      return status;

    } //FullPopulate()

    /*! @brief Populate hierarchy with A's, R's, and P's.

    Invoke a set of factories to populate (construct prolongation,
    restriction, and coarse level discretizations in this
    order) a multigrid Hierarchy starting with information on 'startLevel' 
    and continuing for at most 'numDesiredLevels'. 

    @return  List containing starting and ending level numbers, operator complexity, \#nonzeros in the fine
    matrix, and the sum of nonzeros all matrices (including the fine).
    */
    Teuchos::ParameterList FillHierarchy(PFactory const &PFact, RFactory const &RFact,
                                         TwoLevelFactoryBase const &AcFact,
                                         int const startLevel=0, int const numDesiredLevels=10 ) //TODO: startLevel should be 1!! Because a) it's the way it is in MueMat; b) according to SetLevel(), LevelID of first level=1, not 0
    {
      Monitor h(*this, "FillHierarchy");
      
      // Check for fine level matrix A
      TEST_FOR_EXCEPTION(!Levels_[startLevel]->IsAvailable("A"), Exceptions::RuntimeError, "MueLu::Hierarchy::FillHierarchy(): no fine level matrix A! Set fine level matrix A using Level.Set()");

      Xpetra::global_size_t fineNnz = -1;
      {
        RCP<Operator> A = Levels_[startLevel]->Get< RCP<Operator> >("A");
        fineNnz = A->getGlobalNumEntries();
      }
      Xpetra::global_size_t totalNnz = fineNnz;

      RCP<DefaultFactoryHandler> factoryHandler = rcp(new DefaultFactoryHandler());
      {
        Level & fineLevel = *Levels_[startLevel];
        fineLevel.SetDefaultFactoryHandler(factoryHandler);
      }

      // Setup level structure
      int i = startLevel;
      while (i < startLevel + numDesiredLevels - 1)
        {
    	  // 3.1) setup levels
          Level & fineLevel = *Levels_[i];

          if ((i+1) >= (int) Levels_.size() || Levels_[i+1] == Teuchos::null) {
            RCP<Level> coarseLevel = fineLevel.Build(); // new coarse level, using copy constructor
            this->SetLevel(coarseLevel);                // add to hierarchy
          }

          Level & coarseLevel = *Levels_[i+1];
          coarseLevel.SetDefaultFactoryHandler(factoryHandler);

          // Warning: shift of 1 between i and LevelID. Weird...
          TEST_FOR_EXCEPTION(fineLevel.GetLevelID()   != i+1, Exceptions::RuntimeError, "MueLu::Hierarchy::FillHierarchy(): FineLevel have a wrong level ID");
          TEST_FOR_EXCEPTION(coarseLevel.GetLevelID() != i+2, Exceptions::RuntimeError, "MueLu::Hierarchy::FillHierarchy(): CoarseLevel have a wrong level ID");
          TEST_FOR_EXCEPTION(coarseLevel.GetPreviousLevel() != Levels_[i], Exceptions::RuntimeError, "MueLu::Hierarchy::FillHierarchy(): coarseLevel parent is not fineLevel");

          // 3.2) declare input for levels (recursively)
          GetOStream(Debug, 0) << "declareInput for P's, R's and RAP" << std::endl;
          coarseLevel.Request(PFact);  // TAW: corresponds to SetNeeds
          coarseLevel.Request(RFact);  // TAW: corresponds to SetNeeds
          coarseLevel.Request(AcFact);  // TAW: corresponds to SetNeeds

          ++i;
        } //while

      // Build levels
      i = startLevel;
      while (i < startLevel + numDesiredLevels - 1)
        {
          SubMonitor m(*this, "Level " + Teuchos::toString(i));

          Level & fineLevel   = *Levels_[i];
          Level & coarseLevel = *Levels_[i+1];

          {
            RCP<Operator> A = fineLevel.Get< RCP<Operator> >("A");
            if (A->getRowMap()->getGlobalNumElements() <= maxCoarseSize_) {
              Levels_.resize(i+1);
              break;
            }
          }

          // Request P, R and A at the beginning, because they might be generated by first Build() called, if Pfact==Rfact==AFact...
          coarseLevel.Request("P", &PFact);
          coarseLevel.Request("R", &RFact);
          coarseLevel.Request("A", &AcFact);

          PFact.Build(fineLevel, coarseLevel);
          RCP<Operator> P = coarseLevel.Get< RCP<Operator> >("P", &PFact);
          coarseLevel.Release("P", &PFact);
          coarseLevel.Set("P", P);

          RFact.Build(fineLevel, coarseLevel);
          RCP<Operator> R = coarseLevel.Get< RCP<Operator> >("R", &RFact);
          coarseLevel.Release("R", &RFact);
          coarseLevel.Set("R", R);

          //          TEST_FOR_EXCEPTION((int)Levels_.size() > i, Exceptions::RuntimeError, "WhaT?");

          //TODO: group with the break at the beginning of the loop by changing i to i-1?
          if ((int)Levels_.size() <= i) { //TODO is this the right way to cast?
            Levels_.resize(i+1); //keep only entries 0..i
            break;
          }

          RCP<Operator> Ac = coarseLevel.Get< RCP<Operator> >("A", &AcFact);
          coarseLevel.Release("A", &AcFact);
          coarseLevel.Set("A", Ac);
          
          totalNnz += Ac->getGlobalNumEntries();

          ++i;
        } // while

      for(size_t j = startLevel; j < Levels_.size(); j++) {
        Levels_[j]->SetDefaultFactoryHandler(Teuchos::null);
      }

      ////////////////////////////////////////////////////////////
      //i = startLevel;
      //while (i < startLevel + Levels_.size() - 1)
      //  {
      //    Level & fineLevel = *Levels_[i];
      //    Level & coarseLevel = *Levels_[i+1];
      //    std::cout << "GenericPRFactory.Build: (start)" << std::endl;
      //    std::cout << "FineLevel:" << std::endl;
      //    fineLevel.print(std::cout);
      //    std::cout << "CoarseLevel:" << std::endl;
      //    coarseLevel.print(std::cout);
      //    ++i;
      //  }
      ////////////////////////////////////////////////////////////

      // Gather statistics
      Teuchos::ParameterList status;
      status.set("fine nnz",fineNnz);
      status.set("total nnz",totalNnz);
      status.set("start level",startLevel);
      status.set("end level",i);
      status.set("operator complexity",((SC)totalNnz) / fineNnz);
      return status;
    } //FillHierarchy

    /*! @brief Set solve method for coarsest level.

    @param smooFact  fully constructed SmootherFactory 
    @param pop       whether to use pre,post, or both pre and post smoothing 

    Note: Whether the SmootherFactory builds both a pre- and post-smoother can be also be
    controlled by SmootherFactory::SetSmootherPrototypes. This approach is a bit cumbersome,
    however.
    */
    void SetCoarsestSolver(SmootherFactoryBase const &smooFact, PreOrPost const &pop = BOTH) {
      
      LO clevel = GetNumberOfLevels()-1;
      
      RCP<DefaultFactoryHandler> factoryHandler = rcp(new DefaultFactoryHandler());
      Levels_[clevel]->SetDefaultFactoryHandler(factoryHandler);

      smooFact.BuildSmoother(*Levels_[clevel], pop);

      Levels_[clevel]->SetDefaultFactoryHandler(Teuchos::null);
      //TODO: setdefaultfactoryhandler
    }

    /*! @brief Construct smoothers on all levels but the coarsest.
      TODO need to check whether using Tpetra or Epetra

      Invoke a set of factories to construct smoothers within 
      a multigrid Hierarchy starting with information on 'startLevel' 
      and continuing for at most 'numDesiredLevels'. 

      Note: last level smoother will not be set here. Use SetCoarsestSolver()
      to define a smoother for the last level. Otherwise, a direct solve is
      assumed
    */
    void SetSmoothers(SmootherFactory const &smooFact, LO const &startLevel=0, LO numDesiredLevels=-1)
    {
      Monitor h(*this, "SetSmoothers");

      if (numDesiredLevels == -1)
        numDesiredLevels = GetNumberOfLevels()-startLevel-1;
      LO lastLevel = startLevel + numDesiredLevels - 1;

      //checks
      if (startLevel >= GetNumberOfLevels())
        throw(Exceptions::RuntimeError("startLevel >= actual number of levels"));

      if (startLevel == GetNumberOfLevels() - 1)
        throw(Exceptions::RuntimeError("startLevel == coarse level. Use SetCoarseSolver()"));

      if (lastLevel >= GetNumberOfLevels() - 1) {
        lastLevel = GetNumberOfLevels() - 2;
        GetOStream(Warnings0, 0) << "Warning: coarsest level will have a direct solve!" << std::endl;
      }

      for(int j = startLevel; j<=lastLevel; j++) {
        Levels_[j]->SetDefaultFactoryHandler(Teuchos::null);
      }

      {
        RCP<DefaultFactoryHandler> factoryHandler = rcp(new DefaultFactoryHandler());
        for(int j = startLevel; j<=lastLevel; j++) {
          Levels_[j]->SetDefaultFactoryHandler(factoryHandler);
        }
      }

      for (int i=startLevel; i<=lastLevel; i++) {
        SubMonitor m(*this, "Level " + Teuchos::toString(i));
        smooFact.Build(*Levels_[i]);
      }

      for(int j = startLevel; j<=lastLevel; j++) {
        Levels_[j]->SetDefaultFactoryHandler(Teuchos::null);
      }

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
    void Iterate(MultiVector const &B, LO nIts, MultiVector &X, //TODO: move parameter nIts and default value=1
                 bool const &InitialGuessIsZero=false, CycleType const &Cycle=VCYCLE, LO const &startLevel=0)
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
              GetOStream(Warnings0, 0) << "WARNING:  no coarse grid solver" << std::endl;
          } else {
          //on an intermediate level
          RCP<Level> Coarse = Levels_[startLevel+1];

          //TODO: add IsAvailable test to avoid building default smoother
          RCP<SmootherBase> preSmoo = Fine->Get< RCP<SmootherBase> >("PreSmoother");
          preSmoo->Apply(X, B, zeroGuess);

          RCP<MultiVector> residual = Utils::Residual(*(Fine->Get< RCP<Operator> >("A")),X,B);

          RCP<Operator> P = Coarse->Get< RCP<Operator> >("P");
          RCP<Operator> R;
          RCP<MultiVector> coarseRhs, coarseX;
          if (implicitTranspose_) {
            coarseRhs = MultiVectorFactory::Build(P->getDomainMap(),X.getNumVectors());
            coarseX = MultiVectorFactory::Build(P->getDomainMap(),X.getNumVectors());
            P->apply(*residual,*coarseRhs,Teuchos::TRANS,1.0,0.0);
          } else {
            R = Coarse->Get< RCP<Operator> >("R");
            coarseRhs = MultiVectorFactory::Build(R->getRangeMap(),X.getNumVectors());
            coarseX = MultiVectorFactory::Build(R->getRangeMap(),X.getNumVectors());
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
          RCP<SmootherBase> postSmoo = Fine->Get< RCP<SmootherBase> >("PostSmoother");
          postSmoo->Apply(X, B, false);
        }
        zeroGuess=false;
      } //for (LO i=0; i<nIts; i++)

    } //Iterate()

    //@}
    
  }; //class Hierarchy

  //TODO
  template<class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  std::ostream& operator<<(std::ostream& os, Hierarchy<Scalar,LocalOrdinal,GlobalOrdinal,Node, LocalMatOps> &hierarchy) {
    os << "Printing Hierarchy object" << std::endl;
    typename std::vector< RCP<Level> >::const_iterator i;
    for (i=hierarchy.Levels_.begin(); i != hierarchy.Levels_.end(); ++i)
      os << *(*i) << std::endl;
    return os;
  }

} //namespace MueLu

#define MUELU_HIERARCHY_SHORT

#endif //ifndef MUELU_HIERARCHY_HPP
// TODO: We need a Set/Get function to change the CycleType (for when Iterate() calls are embedded in a Belos Preconditionner for instance).

//TODO: GetNumbertOfLevels() -> GetNumLevels()
