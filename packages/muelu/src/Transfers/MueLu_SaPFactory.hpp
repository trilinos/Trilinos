#ifndef MUELU_SAPFACTORY_HPP
#define MUELU_SAPFACTORY_HPP

#include <Cthulhu_Map.hpp>
#include <Cthulhu_EpetraMap.hpp>
#include <Cthulhu_CrsMatrix.hpp>
#include <Cthulhu_EpetraCrsMatrix.hpp>
#include <Cthulhu_CrsOperator.hpp>
#include <Cthulhu.hpp>
#include <Cthulhu_Vector.hpp>
#include <Cthulhu_VectorFactory.hpp>

#include "EpetraExt_MatrixMatrix.h"

#include "MueLu_OperatorFactory.hpp"
#include "MueLu_Utilities.hpp"
#include "MueLu_MatrixFactory.hpp"
#include "MueLu_TentativePFactory.hpp"

#include <iostream>

namespace MueLu {

/*!
  @class SaPFactory class.
  @brief Factory for building Smoothed Aggregation prolongators.

  Right now this factory assumes a 1D problem.  Aggregation is hard-coded to divide
  the # fine dofs by 3.
*/

template<class ScalarType, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
class SaPFactory : public OperatorFactory<ScalarType,LocalOrdinal,GlobalOrdinal,Node, LocalMatOps> {
#include "MueLu_UseShortNames.hpp"

  template<class AA, class BB, class CC, class DD, class EE>
  inline friend std::ostream& operator<<(std::ostream& os, SaPFactory<AA,BB,CC,DD, EE> &factory);

  private:
/*
     AggFact_
     CoalesceFact_
     diagonalView_ = 'current' % diagonal view label (default == current view)
*/
     GO maxCoarseSize_;
     ScalarType dampingFactor_;
     bool doQR_;
     bool reUseP_;
     bool reUsePtent_;
     bool reUseGraph_;
     bool reUseAggregates_;

  public:
    //! @name Constructors/Destructors.
    //@{

    //! Constructor.
    SaPFactory() {
      Teuchos::OSTab tab(this->out_);
      MueLu_cout(Teuchos::VERB_HIGH) << "SaPFactory: Instantiating a new factory" << std::endl;
    }

    //! Destructor.
    virtual ~SaPFactory() {}
    //@}

    //! @name Build methods.
    //@{

    /*!
      @brief Build method.

      Builds smoothed aggregation prolongator and returns it in <tt>coarseLevel</tt>.
    */
    bool Build(Level &fineLevel, Level &coarseLevel) {
      Teuchos::OSTab tab(this->out_);
      MueLu_cout(Teuchos::VERB_HIGH) << "SaPFactory: Building a prolongator" << std::endl;

      RCP<Operator> Ptent = TentativePFactory::MakeTentative(fineLevel);

      //Build the smoother prolongator
      Teuchos::RCP< Operator > Op = fineLevel.GetA();
      RCP<Operator> D = MueLu::Utils<SC, LO, GO, NO, LMO>::BuildMatrixInverseDiagonal(Op);
      RCP<Operator> AP = MueLu::Utils<SC, LO, GO, NO, LMO>::TwoMatrixMultiply(Op,Ptent);
      RCP<Operator> DAP = MueLu::Utils<SC, LO, GO, NO, LMO>::TwoMatrixMultiply(D,AP);
      RCP<Operator> smP = MueLu::Utils<SC, LO, GO, NO, LMO>::TwoMatrixAdd(Ptent,DAP,1.0,-2.0/3.0);

      coarseLevel.SetP(smP);

      //MatrixPrint(smP);

      return true;
    }
    //@}

    //! @name Set methods.
    //@{
    void SetMaxCoarseSize(GO maxCoarseSize) {
      maxCoarseSize_ = maxCoarseSize;
    }

    void SetDampingFactor(ScalarType dampingFactor) {
      dampingFactor_ = dampingFactor;
    }

    void TentativeWithQR(bool value) {
      doQR_ = value;
    }

    void ReUseP(bool value) {
      reUseP_ = value;
    }

    void ReUsePtent(bool value) {
      reUsePtent_ = value;
    }

    void ReUseAggregates(bool value) {
      reUseAggregates_ = value;
    }

    void ReUseGraph(bool value) {
      reUseGraph_ = value;
    }
    //@}

    //! @name Get methods.
    //@{

    GO GetMaxCoarseSize() {
      return maxCoarseSize_;
    }

    ScalarType GetDampingFactor() {
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

    bool ReUseAggregates() {
      return reUseAggregates_;
    }

    bool ReUseGraph() {
      return reUseGraph_;
    }

    //@}
/*
//TODO
function [this] = SaPFactory(CoalesceFact,AggFact, diagonalView) //copy ctor
function SetDiagonalView(this, diagonalView)
function [z] = GetNeeds(this)
function  AggInfo = BuildAggregates(FineLevel, ...
function  [P] = MakeNoQRTentative(AggInfo,Amat,nullspace,OutputLevel)
*/

}; //class SaPFactory

//! Friend print function.
template<class ScalarType, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
std::ostream& operator<<(std::ostream& os, SaPFactory<ScalarType,LocalOrdinal,GlobalOrdinal,Node, LocalMatOps> &factory) {
  os << "Printing an SaPFactory object" << std::endl;
  return os;
}

} //namespace MueLu

#define MUELU_SAPFACTORY_SHORT

#endif //ifndef MUELU_SAPFACTORY_HPP
