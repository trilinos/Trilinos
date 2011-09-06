#ifndef MUELU_HIERARCHY_HPP
#define MUELU_HIERARCHY_HPP

#include "MueLu_ConfigDefs.hpp"
#include "MueLu_BaseClass.hpp"
#include "MueLu_DefaultFactoryHandler.hpp"
#include "MueLu_Types.hpp"
#include "MueLu_Level.hpp"

#include "MueLu_Exceptions.hpp"
#include "MueLu_GenericPRFactory.hpp"
#include "MueLu_SmootherFactory.hpp"

// used as default:
#include "MueLu_SaPFactory.hpp"
#include "MueLu_TransPFactory.hpp"
#include "MueLu_RAPFactory.hpp"
#include "MueLu_IfpackSmoother.hpp"
#include "MueLu_Ifpack2Smoother.hpp"

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
    //! print residual history during iteration
    bool printResidualHistory_;
    bool implicitTranspose_;

    RCP<DefaultFactoryHandlerBase> defaultFactoryHandler_;

  protected:
    RCP<Teuchos::FancyOStream> out_;

  public:

    //! @name Constructors/Destructors
    //@{

    //! Default constructor.
    Hierarchy() : printResidualHistory_(false), implicitTranspose_(false), defaultFactoryHandler_(rcp(new DefaultFactoryHandler())), out_(this->getOStream()) {}

    //! Copy constructor.
    Hierarchy(Hierarchy const &inHierarchy) {
      std::cerr << "Not implemented yet." << std::endl;
    }

    //! Destructor.
    virtual ~Hierarchy() {}

    //@}

    //! @name Set/Get Methods.
    //@{

    //! Assign a level to hierarchy.
    // TODO from JG: behavior should be explain or changed. 
    //               Right now, we change the LevelID of the input level and push it at the end of hierarchy.
    //               Certainly better to take the LevelID of the input level into account
    void SetLevel(RCP<Level> const& level) {
      Levels_.push_back(level);
      level->SetLevelID(Levels_.size());

      level->SetDefaultFactoryHandler(defaultFactoryHandler_);

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

    //! Indicate whether to print residual history to a FancyOStream.
    void PrintResidualHistory(bool value) {
      printResidualHistory_ = value;
    }

    //! If this returns true, the residual history will print to a FancyOStream.
    bool PrintResidualHistory() {
      return printResidualHistory_;
    }

    //TODO: allow users to change default factory handler.
    // void SetDefaultFactory(const std::string&, RCP<Factory>&) {
    //   defaultFactoryHandler_->SetDefaultFactory(...)
    // }

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
    Teuchos::ParameterList FullPopulate(RCP<PRFactory> PRFact=Teuchos::null,
                                        RCP<TwoLevelFactoryBase> AcFact=Teuchos::null,
                                        RCP<SmootherFactory> SmooFact=Teuchos::null,
                                        int const &startLevel=0, int const &numDesiredLevels=10 )
    {
      if (PRFact == Teuchos::null) {
        RCP<SaPFactory> SaPFact = rcp( new SaPFactory() );
        PRFact = rcp( new GenericPRFactory(SaPFact));
      }
      if (AcFact == Teuchos::null) AcFact = rcp( new RAPFactory());

      Teuchos::ParameterList status;
      status = FillHierarchy(*PRFact,*AcFact,startLevel,numDesiredLevels /*,status*/);
      if (SmooFact != Teuchos::null) {
        SetSmoothers(*SmooFact,startLevel,numDesiredLevels-1);
      }
      return status;

    } //FullPopulate()


    //TODO JG to JHU: do we need this ??
    /*! @brief Populate hierarchy with A's, R's, and P's using default factories.

    The prolongator factory defaults to Smoothed Aggregation, and the coarse matrix factory defaults
    to the Galerkin product.

    @return  List containing starting and ending level numbers, operator complexity, \#nonzeros in the fine
    matrix, and the sum of nonzeros all matrices (including the fine).
    */
    Teuchos::ParameterList FillHierarchy() {
      RCP<SaPFactory> PFact = rcp(new SaPFactory());
      RCP<TransPFactory> RFact = rcp(new TransPFactory());
      RCP<GenericPRFactory>  PRFact = rcp(new GenericPRFactory(PFact,RFact));
      RCP<RAPFactory> AcFact = rcp(new RAPFactory());
      Teuchos::ParameterList status;
      status = FillHierarchy(*PRFact,*AcFact);
      return status;
    } //FillHierarchy()

    /*! @brief Populate hierarchy with A's, R's, and P's.

    Populate hierarchy with provided PRFactory.  The coarse matrix factory defaults
    to the Galerkin product.

    @return  List containing starting and ending level numbers, operator complexity, \#nonzeros in the fine
    matrix, and the sum of nonzeros all matrices (including the fine).
    */
    Teuchos::ParameterList FillHierarchy(PRFactory const &PRFact) {
      RAPFactory AcFact;
      Teuchos::ParameterList status;
      status = FillHierarchy(PRFact,AcFact);
      return status;
    }

    //TODO should there be another version of FillHierarchy:
    //TODO   Teuchos::ParameterList FillHierarchy(TwoLevelFactoryBase const &AcFact)
    //TODO where the PRFactory is automatically generated?
     
    /*! @brief Populate hierarchy with A's, R's, and P's.

    Invoke a set of factories to populate (construct prolongation,
    restriction, and coarse level discretizations in this
    order) a multigrid Hierarchy starting with information on 'startLevel' 
    and continuing for at most 'numDesiredLevels'. 

    @return  List containing starting and ending level numbers, operator complexity, \#nonzeros in the fine
    matrix, and the sum of nonzeros all matrices (including the fine).
    */
    Teuchos::ParameterList FillHierarchy(PRFactory const &PRFact,
                                         TwoLevelFactoryBase const &AcFact,
                                         int startLevel=0, int numDesiredLevels=10 ) //TODO: startLevel should be 1!! Because a) it's the way it is in MueMat; b) according to SetLevel(), LevelID of first level=1, not 0
    {

      RCP<Operator> A = Levels_[startLevel]->Get< RCP<Operator> >("A");

      // Set default, very important to do that! (Otherwise, factory use default factories instead of user defined factories - ex: RAPFactory will request a new "P" from default factory)
      defaultFactoryHandler_->SetDefaultFactory("P", rcpFromRef(PRFact)); // TODO: remove rcpFromRef
      defaultFactoryHandler_->SetDefaultFactory("R", rcpFromRef(PRFact));
      defaultFactoryHandler_->SetDefaultFactory("A", rcpFromRef(AcFact));

      Xpetra::global_size_t fineNnz = A->getGlobalNumEntries();
      Xpetra::global_size_t totalNnz = fineNnz;

      bool goodBuild=true;
      int i = startLevel;
      while (i < startLevel + numDesiredLevels - 1)
        {
          Level & fineLevel = *Levels_[i];

          if ((i+1) >= (int) Levels_.size() || Levels_[i+1] == Teuchos::null) {
            RCP<Level> coarseLevel = fineLevel.Build(*out_); // new coarse level, using copy constructor
            this->SetLevel(coarseLevel);                     // add to hierarchy
          }
          
          Level & coarseLevel = *Levels_[i+1];

          // Warning: shift of 1 between i and LevelID. Weird...
          TEST_FOR_EXCEPTION(fineLevel.GetLevelID()   != i+1, Exceptions::RuntimeError, "MueLu::Hierarchy::FillHierarchy(): FineLevel have a wrong level ID");
          TEST_FOR_EXCEPTION(coarseLevel.GetLevelID() != i+2, Exceptions::RuntimeError, "MueLu::Hierarchy::FillHierarchy(): CoarseLevel have a wrong level ID");
          TEST_FOR_EXCEPTION(coarseLevel.GetPreviousLevel() != Levels_[i], Exceptions::RuntimeError, "MueLu::Hierarchy::FillHierarchy(): coarseLevel parent is not fineLevel");

          *out_ << "starting build of P's and R's"  << std::endl;
          PRFact.DeclareInput(fineLevel, coarseLevel);  // TAW: corresponds to SetNeeds
          goodBuild = PRFact.Build(fineLevel, coarseLevel);
          if ((int)Levels_.size() <= i) goodBuild=false; //TODO is this the right way to cast?
          if (!goodBuild) {
            Levels_.resize(i+1); //keep only entries 0..i
            break;
          }
          *out_ << "starting build of RAP"  << std::endl;
          AcFact.DeclareInput(fineLevel, coarseLevel); // TAW: corresponds to SetNeeds
          if ( !AcFact.Build(fineLevel, coarseLevel) ) {
            Levels_.resize(i+1); //keep only entries 0..i
            break;
          }
          //RCP<Operator> A = coarseLevel.Get< RCP<Operator> >("A");
          totalNnz += coarseLevel.Get< RCP<Operator> >("A")->getGlobalNumEntries();

          ++i;
        } //while

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
      smooFact.BuildSmoother(*Levels_[clevel], pop);
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
        MueLu_cout(Teuchos::VERB_HIGH)
          << "Warning: Coarsest Level will have a direct solve!" << std::endl;
      }

      for (int i=startLevel; i<=lastLevel; i++) {
        smooFact.Build(*Levels_[i]);
      }

    } //SetSmoothers()

    // JG to JHU: Do we need this ??
    //      /*! @brief Construct smoothers on all levels but the coarsest.
    //        Defaults to using IfpackSmoother factory that generates Gauss-Seidel smoothers.
    //      */
    //      void SetSmoothers()
    //      {
    // #ifdef HAVE_MUELU_IFPACK
    //        Teuchos::ParameterList  ifpackList;
    // //FIXME #ifdef we're using tpetra
    // //FIXME    throw(Exceptions::NotImplemented("No default smoother is defined"));
    // //FIXME #else we're using epetra
    //        ifpackList.set("relaxation: type", "Gauss-Seidel");
    //        ifpackList.set("relaxation: sweeps", (int) 1);
    //        ifpackList.set("relaxation: damping factor", (double) 1.0);
    //        ifpackList.set("relaxation: zero starting solution", false);
    //        RCP<IfpackSmoother>  smoother = rcp( new IfpackSmoother("point relaxation stand-alone",ifpackList) );
    //        SmootherFactory smooFact(smoother);
    // //FIXME #endif
    //        SetSmoothers(smooFact);
    // #endif
    //      }

    //         typedef typename Teuchos::ScalarTraits<SC>::magnitudeType Magnitude;
    //FIXME delete this macro
#define GimmeNorm(someVec,someLabel) {(someVec).norm2(norms);   \
      *out_ << someLabel << " = " << norms<< std::endl;}

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
      //Teuchos::OSTab tab(*out_);
      bool zeroGuess=InitialGuessIsZero;

      for (LO i=0; i<nIts; i++) {

        RCP<Level> Fine = Levels_[startLevel];

        if (printResidualHistory_ && startLevel==0) {
          *out_ << "iter:    "
                << std::setiosflags(std::ios::left)
                << std::setprecision(3) << i;
          *out_ << "           residual = "
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
            if (Fine->IsAvailable("PreSmoother")) {
              RCP<SmootherBase> preSmoo = Fine->Get< RCP<SmootherBase> >("PreSmoother");
              preSmoo->Apply(X, B, false);
              emptySolve=false;
            }
            if (Fine->IsAvailable("PostSmoother")) {
              RCP<SmootherBase> postSmoo = Fine->Get< RCP<SmootherBase> >("PostSmoother");
              postSmoo->Apply(X, B, false); 
              emptySolve=false;
            }
            if (emptySolve==true)
              *out_ << "WARNING:  no coarse grid solve!!" << std::endl;
          } else {
          //on an intermediate level
          RCP<Level> Coarse = Levels_[startLevel+1];

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
          RCP<SmootherBase> postSmoo = Fine->Get< RCP<SmootherBase> >("PostSmoother");
          postSmoo->Apply(X, B, false);
        }
        zeroGuess=false;
      } //for (LO i=0; i<nIts; i++)

    } //Iterate()

    //@}
    
  }; //class Hierarchy

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
