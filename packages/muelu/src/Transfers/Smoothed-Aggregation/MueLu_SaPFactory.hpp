#ifndef MUELU_SAPFACTORY_HPP
#define MUELU_SAPFACTORY_HPP

#include <Cthulhu_Map.hpp>
#include <Cthulhu_CrsMatrix.hpp>
#include <Cthulhu_CrsOperator.hpp>
#include <Cthulhu_Vector.hpp>
#include <Cthulhu_VectorFactory.hpp>

#ifdef HAVE_MUELU_EPETRAEXT
#include "EpetraExt_MatrixMatrix.h"
#endif

#include "MueLu_PFactory.hpp"
#include "MueLu_Utilities.hpp"
#include "MueLu_MatrixFactory.hpp"
#include "MueLu_TentativePFactory.hpp"
#include "MueLu_UCAggregationFactory.hpp"
#include "MueLu_Exceptions.hpp"

#include <iostream>

namespace MueLu {

/*!
  @class SaPFactory class.
  @brief Factory for building Smoothed Aggregation prolongators.

  Right now this factory assumes a 1D problem.  Aggregation is hard-coded to divide
  the # fine dofs by 3.
*/

template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
class SaPFactory : public PFactory<Scalar,LocalOrdinal,GlobalOrdinal,Node, LocalMatOps> {
#include "MueLu_UseShortNames.hpp"

  template<class AA, class BB, class CC, class DD, class EE>
  inline friend std::ostream& operator<<(std::ostream& os, SaPFactory<AA,BB,CC,DD, EE> &factory);

  private:
/*
     TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO 
     RCP<MueLu::AggregationFactory<LO,GO,NO,LMO> > AggFact_;
     CoalesceFact_
     TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO 
*/
     Teuchos::RCP<PFactory> initialPFact_;
     std::string diagonalView_;
     bool doQR_;
     Scalar dampingFactor_;
     bool useAFiltered_;
     bool reUseP_;
     bool reUsePtent_;

  public:
    //! @name Constructors/Destructors.
    //@{

    /*! @brief Default constructor.

    */
    SaPFactory() : diagonalView_("current"),
                   //AggFact_(Teuchos::null),
                   doQR_(false), dampingFactor_(4./3), useAFiltered_(false), reUseP_(false),
                   reUsePtent_(false)
                   //, PFactory::reUseGraph_(false), PFactory::reUseAggregates_(false)
    {
      initialPFact_ = Teuchos::rcp(new TentativePFactory()),
      PFactory::reUseGraph_=false;
      PFactory::reUseAggregates_=false;
      //Teuchos::OSTab tab(this->out_);
      //MueLu_cout(Teuchos::VERB_HIGH) << "SaPFactory: Instantiating a new factory" << std::endl;
    }

    /*! @brief Constructor.

        User can supply a factory for generating the tentative prolongator.
    */
    SaPFactory(RCP<PFactory> InitialPFact) : initialPFact_(InitialPFact), diagonalView_("current"),
                   //AggFact_(Teuchos::null),
                   doQR_(false), dampingFactor_(4./3), useAFiltered_(false), reUseP_(false),
                   reUsePtent_(false)
                   //, PFactory::reUseGraph_(false), PFactory::reUseAggregates_(false)
    {
      PFactory::reUseGraph_=false;
      PFactory::reUseAggregates_=false;
    }

    //! Destructor.
    virtual ~SaPFactory() {}
    //@}

    //! @name Build methods.
    //@{

    /*!
      @brief Build method.

      Builds smoothed aggregation prolongator and returns it in <tt>coarseLevel</tt>.
      //FIXME what does the return code mean (unclear in MueMat)?
      //FIXME how should nullspace be stored?
    */
    bool BuildP(Level &fineLevel, Level &coarseLevel) const {
      Teuchos::OSTab tab(this->out_);

      Teuchos::RCP<Teuchos::Time> timer = rcp(new Teuchos::Time("SaPFactory::BuildP"));
      timer->start(true);

      RCP<Operator> finalP;

      if (reUseP_) {
        if (coarseLevel.GetP() == Teuchos::null)
          throw(std::runtime_error("SaPFactory: you have asked to reuse P, but it doesn't exist"));
        if (coarseLevel.IsSaved("Nullspace") == false)
          throw(std::runtime_error("SaPFactory: you have asked to reuse cnull, but it doesn't exist"));
        return true;
      }

      //TODO get or generate fine grid nullspace here
      RCP<MultiVector> fineNullspace;
      if (fineLevel.IsSaved("Nullspace"))
        fineLevel.CheckOut("Nullspace",fineNullspace);
      else {
        //TODO add this functionality
        //throw(Exceptions::NotImplemented("SaPFactory.Build():  nullspace generation not implemented yet"));
        std::cout << "nullspace generation not implemented yet" << std::endl;
      }

      coarseLevel.Request("Ptent");
      coarseLevel.Request("Nullspace");
      initialPFact_->BuildP(fineLevel,coarseLevel);
      RCP<Operator> Ptent;
      coarseLevel.CheckOut("Ptent",Ptent);
      RCP<MultiVector> coarseNullspace;
      if (reUsePtent_) {
        try {
          coarseLevel.CheckOut("Ptent",Ptent); //FIXME throws an error, replace with recomputation
          coarseLevel.CheckOut("Nullspace",coarseNullspace); //FIXME throws an error, replace with recomputation
        }
        catch(...) {
          throw(Exceptions::NotImplemented("SaPFactory.Build(): regeneration of Ptent/nullspace not implemented yet"));
        }
      }

      if (coarseLevel.IsRequested("Ptent"))
        coarseLevel.Save("Ptent",Ptent);
      

      //Build final prolongator

      SC lambdaMax = 2.0; //FIXME hard-coded for right now for 1D constant-coefficient Poisson
      //FIXME Cthulhu::Operator should calculate/stash max eigenvalue
      //FIXME SC lambdaMax = Op->GetDinvALambda();

      if (dampingFactor_ != 0) {
        Teuchos::RCP< Operator > Op = fineLevel.GetA();
        RCP<Operator> D = Utils::BuildMatrixInverseDiagonal(Op);
        RCP<Operator> AP = Utils::TwoMatrixMultiply(Op,Ptent);
        RCP<Operator> DAP = Utils::TwoMatrixMultiply(D,AP);
        finalP = Utils::TwoMatrixAdd(Ptent,DAP,1.0,-dampingFactor_/lambdaMax);
      }
      else {
        finalP = Ptent;
      }

      coarseLevel.SetP(finalP);
      //coarseLevel.Save("Nullspace",coarseNullspace);

      //Utils::MatrixPrint(finalP);

      timer->stop();
      Utils::ReportTimeAndMemory(*timer, *(finalP->getRowMap()->getComm()));

      return true;
    } //Build()
    //@}

    //! @name Set methods.
    //@{

    //! Set prolongator smoother damping factor.
    void SetDampingFactor(Scalar dampingFactor) {
      dampingFactor_ = dampingFactor;
    }

    void TentativeWithQR(bool value) {
      doQR_ = value;
    }

    //! Change view of diagonal.
    void SetDiagonalView(std::string const& diagView) {
      diagonalView_ = diagView;
    }

    void SetUseAFiltered(bool value) {
      throw(Exceptions::NotImplemented("SetUseAFiltered not fully implemented"));
      useAFiltered_ = value;
      //FIXME add to needs?
    }

    void ReUseP(bool value) {
      reUseP_ = value;
    }

    void ReUsePtent(bool value) {
      reUsePtent_ = value;
    }

    //@}

    //! @name Get methods.
    //@{

    //! Returns prolongator smoother damping factor.
    Scalar GetDampingFactor() {
      return dampingFactor_;
    }

    bool TentativeWithQR() {
      return doQR_;
    }

    bool ReUseP() {
      return reUseP_;
    }

    bool ReUsePtent() {
      return reUsePtent_;
    }

    //! Returns current view of diagonal.
    std::string GetDiagonalView() {
      return diagonalView_;
    }

    //@}
/*
//TODO
function [this] = SaPFactory(CoalesceFact,AggFact, diagonalView) //copy ctor
function SetDiagonalView(this, diagonalView)
*/

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
