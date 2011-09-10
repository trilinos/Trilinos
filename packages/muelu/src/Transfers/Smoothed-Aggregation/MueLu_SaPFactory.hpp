#ifndef MUELU_SAPFACTORY_HPP
#define MUELU_SAPFACTORY_HPP

#include <iostream>

#include <Xpetra_Map.hpp>
#include <Xpetra_CrsMatrix.hpp>
#include <Xpetra_CrsOperator.hpp>
#include <Xpetra_Vector.hpp>
#include <Xpetra_VectorFactory.hpp>

#include "MueLu_ConfigDefs.hpp"
#include "MueLu_PFactory.hpp"
#include "MueLu_Utilities.hpp"
#include "MueLu_MatrixFactory.hpp"
#include "MueLu_TentativePFactory.hpp"
#include "MueLu_UCAggregationFactory.hpp"
#include "MueLu_Exceptions.hpp"

namespace MueLu {

  /*!
    @class SaPFactory class.
    @brief Factory for building Smoothed Aggregation prolongators.

    Right now this factory assumes a 1D problem.  Aggregation is hard-coded to divide
    the # fine dofs by 3.
  */

  template <class Scalar = double, class LocalOrdinal = int, class GlobalOrdinal = LocalOrdinal, class Node = Kokkos::DefaultNode::DefaultNodeType, class LocalMatOps = typename Kokkos::DefaultKernels<void,LocalOrdinal,Node>::SparseOps> //TODO: or BlockSparseOp ?
  class SaPFactory : public PFactory {
#include "MueLu_UseShortNames.hpp"

    template<class AA, class BB, class CC, class DD, class EE>
    inline friend std::ostream& operator<<(std::ostream& os, SaPFactory<AA,BB,CC,DD, EE> &factory);

  public:

    //! @name Constructors/Destructors.
    //@{
  
    /*! @brief Constructor.
      User can supply a factory for generating the tentative prolongator.
    */
    SaPFactory(RCP<PFactory> InitialPFact = Teuchos::null, RCP<SingleLevelFactoryBase> AFact = Teuchos::null)
      : initialPFact_(InitialPFact), AFact_(AFact),
        dampingFactor_(4./3), diagonalView_("current") {
      if(initialPFact_ == Teuchos::null)
      {
        // todo: if no initial PFactory, set it to TentativePFactory
        //std::string msg = "SaPFactory: no initial P factory defined";
        //throw(Exceptions::RuntimeError(msg));
        initialPFact_ = rcp(new TentativePFactory());
      }
    }
  
    //! Destructor.
    virtual ~SaPFactory() {}
  
    //@}

    //! @name Set methods.
    //@{

    //! Set prolongator smoother damping factor.
    void SetDampingFactor(Scalar dampingFactor) {
      dampingFactor_ = dampingFactor;
    }

    //! Change view of diagonal.
    void SetDiagonalView(std::string const& diagView) {
      diagonalView_ = diagView;
    }
    //@}

    //! @name Get methods.
    //@{

    //! Returns prolongator smoother damping factor.
    Scalar GetDampingFactor() {
      return dampingFactor_;
    }

    //! Returns current view of diagonal.
    std::string GetDiagonalView() {
      return diagonalView_;
    }

    //@}

    //! Input
    //@{

    void DeclareInput(Level &fineLevel, Level &coarseLevel) const {
      if(!coarseLevel.IsRequested("P",initialPFact_))
        initialPFact_->DeclareInput(fineLevel,coarseLevel);
      if(!fineLevel.IsRequested("A",AFact_))
        AFact_->DeclareInput(fineLevel);
      coarseLevel.Request("P",initialPFact_);
      fineLevel.Request("A",AFact_);
    };

    //@}

    //! @name Build methods.
    //@{
  
    /*!
      @brief Build method.

      Builds smoothed aggregation prolongator and returns it in <tt>coarseLevel</tt>.
      //FIXME what does the return code mean (unclear in MueMat)?
      //FIXME how should nullspace be stored?
      */
    bool Build(Level& fineLevel, Level &coarseLevel) const {
      return BuildP(fineLevel,coarseLevel);
    }

    bool BuildP(Level &fineLevel, Level &coarseLevel) const {
      std::ostringstream buf; buf << coarseLevel.GetLevelID();
      RCP<Teuchos::Time> timer = rcp(new Teuchos::Time("SaPFactory::BuildP_"+buf.str()));
      timer->start(true);

      // Level Get
      RCP<Operator> A     = fineLevel.  Get< RCP<Operator> >("A", AFact_);
      RCP<Operator> Ptent = coarseLevel.Get< RCP<Operator> >("P", initialPFact_);

      fineLevel.Release("A", AFact_);
      coarseLevel.Release("P", initialPFact_);

      //Build final prolongator
      RCP<Operator> finalP; // output

      //FIXME Xpetra::Operator should calculate/stash max eigenvalue
      //FIXME SC lambdaMax = A->GetDinvALambda();

      if (dampingFactor_ != 0) {

        RCP<Teuchos::Time> sapTimer;
        //sapTimer = rcp(new Teuchos::Time("SaPFactory:I * Ptent"));
        //sapTimer->start(true);
        //Teuchos::ParameterList matrixList;
        //RCP<Operator> I = MueLu::Gallery::CreateCrsMatrix<SC,LO,GO, Map,CrsOperator>("Identity",fineLevel.Get< RCP<Operator> >("A")->getRowMap(),matrixList);
        //RCP<Operator> newPtent = Utils::TwoMatrixMultiply(I,false,Ptent,false);
        //Ptent = newPtent; //I tried a checkout of the original Ptent, and it seems to be gone now (which is good)
        //sapTimer->stop();
        //MemUtils::ReportTimeAndMemory(*sapTimer, *(A->getRowMap()->getComm()));

        sapTimer = rcp(new Teuchos::Time("SaPFactory:APtent_"+buf.str()));
        sapTimer->start(true);

        //JJH -- If I switch doFillComplete to false, the resulting matrix seems weird when printed with describe.
        //JJH -- The final prolongator is wrong, to boot.  So right now, I fillComplete AP, but avoid fillComplete
        //JJH -- in the scaling.  Long story short, we're doing 2 fillCompletes, where ideally we'd do just one.
        bool doFillComplete=true;
        bool optimizeStorage=false;
        RCP<Operator> AP = Utils::TwoMatrixMultiply(A,false,Ptent,false,doFillComplete,optimizeStorage);
        sapTimer->stop();
        MemUtils::ReportTimeAndMemory(*sapTimer, *(A->getRowMap()->getComm()));

        sapTimer = rcp(new Teuchos::Time("SaPFactory:Dinv_APtent_"+buf.str()));
        sapTimer->start(true);
        doFillComplete=false;
        optimizeStorage=false;
        Teuchos::ArrayRCP<SC> diag = Utils::GetMatrixDiagonal(A);
        Utils::MyOldScaleMatrix(AP,diag,true,doFillComplete,optimizeStorage); //scale matrix with reciprocal of diag
        sapTimer->stop();
        MemUtils::ReportTimeAndMemory(*sapTimer, *(A->getRowMap()->getComm()));

        sapTimer = rcp(new Teuchos::Time("SaPFactory:eigen_estimate_"+buf.str()));
        sapTimer->start(true);
        Scalar lambdaMax = Utils::PowerMethod(*A, true, (LO) 10,(Scalar)1e-4);
        sapTimer->stop();
        MemUtils::ReportTimeAndMemory(*sapTimer, *(A->getRowMap()->getComm()));
        RCP<const Teuchos::Comm<int> > comm = A->getRowMap()->getComm();
        if (comm->getRank() == 0)
          std::cout << "damping factor = " << dampingFactor_/lambdaMax << " ("
                    << dampingFactor_ << " / " << lambdaMax << ")" << std::endl;

        sapTimer = rcp(new Teuchos::Time("SaPFactory:Pt_plus_DinvAPtent_"+buf.str()));
        sapTimer->start(true);

        bool doTranspose=false; 
        if (AP->isFillComplete())
          Utils::TwoMatrixAdd(Ptent,doTranspose,1.0,AP,doTranspose,-dampingFactor_/lambdaMax,finalP);
        else {
          Utils::TwoMatrixAdd(Ptent,doTranspose,1.0,AP,-dampingFactor_/lambdaMax);
          finalP = AP;
        }
        sapTimer->stop();
        MemUtils::ReportTimeAndMemory(*sapTimer, *(A->getRowMap()->getComm()));

        sapTimer = rcp(new Teuchos::Time("SaPFactory:finalP_fillComplete_"+buf.str()));
        sapTimer->start(true);
        finalP->fillComplete( Ptent->getDomainMap(), Ptent->getRangeMap() );
        sapTimer->stop();
        MemUtils::ReportTimeAndMemory(*sapTimer, *(A->getRowMap()->getComm()));
      }
      else {
        finalP = Ptent;
      }

      timer->stop();
      MemUtils::ReportTimeAndMemory(*timer, *(finalP->getRowMap()->getComm()));

      // Level Set
      coarseLevel.Set("P", finalP, this);

      return true;
    } //Build()
    //@}

    /*
    //TODO
    function [this] = SaPFactory(CoalesceFact,AggFact, diagonalView) //copy ctor
    function SetDiagonalView(this, diagonalView)
    */


  private:

    //! Input factories
    RCP<PFactory> initialPFact_; //! Ptentative Factory
    RCP<SingleLevelFactoryBase> AFact_;     //! A Factory
    
    //! Factory parameters
    Scalar dampingFactor_;
    std::string diagonalView_;

    // TODO: test ReUse

  }; //class SaPFactory

  //! Friend print function.
  template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  std::ostream& operator<<(std::ostream& os, SaPFactory<Scalar,LocalOrdinal,GlobalOrdinal,Node, LocalMatOps> &factory) {
    os << "Printing an SaPFactory object" << std::endl;
    return os;
  }

} //namespace MueLu

#define MUELU_SAPFACTORY_SHORT

#endif //ifndef MUELU_SAPFACTORY_HPP
