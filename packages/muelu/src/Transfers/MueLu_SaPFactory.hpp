#ifndef MUELU_SAPFACTORY_HPP
#define MUELU_SAPFACTORY_HPP

#include <iostream>
#include "MueLu_OperatorFactory.hpp"

namespace MueLu {

/*!
  @class SaPFactory class.
  @brief Factory for building Smoothed Aggregation prolongators.
*/

template<class Scalar, class LO, class GO, class Node>
class SaPFactory : public OperatorFactory<Scalar,LO,GO,Node> {
  template<class AA, class BB, class CC, class DD>
  inline friend std::ostream& operator<<(std::ostream& os, SaPFactory<AA,BB,CC,DD> &factory);

  private:
/*
     CoalesceFact_
     AggFact_
     diagonalView_ = 'current' % diagonal view label (default == current view)
*/
     GO maxCoarseSize_;
     Scalar dampingFactor_;
     bool doQR_;
     bool reUseP_;
     bool reUsePtent_;
     bool reUseGraph_;
     bool reUseAggregates_;

  public:
    //@{ Constructors/Destructors.
    //! Constructor.
    SaPFactory() {
      Teuchos::OSTab tab(this->out_);
      *(this->out_) << "SaPFactory: Instantiating a new factory" << std::endl;
    }

    //! Destructor.
    virtual ~SaPFactory() {}
    //@}

    //@{ Build methods.

    //! Build method.
    bool Build(Level<Scalar,LO,GO,Node> &fineLevel, Level<Scalar,LO,GO,Node> &coarseLevel) {
      Teuchos::OSTab tab(this->out_); *(this->out_) << "SaPFactory: Building a prolongator" << std::endl; return true;
    }
    //@}

    //@{
    // @name Set methods.
    void SetMaxCoarseSize(GO maxCoarseSize) {
      maxCoarseSize_ = maxCoarseSize;
    }

    void SetDampingFactor(Scalar dampingFactor) {
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

    //@{
    // @name Get methods.

    GO GetMaxCoarseSize() {
      return maxCoarseSize_;
    }

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
function  [P,CoarseNull] = MakeTentative(AggInfo,Amat,nullspace,QROnOrOff,OutputLevel)
function  [P] = MakeNoQRTentative(AggInfo,Amat,nullspace,OutputLevel)
*/

}; //class SaPFactory

//! Friend print function.
template<class Scalar, class LO, class GO, class Node>
std::ostream& operator<<(std::ostream& os, SaPFactory<Scalar,LO,GO,Node> &factory) {
  os << "Printing an SaPFactory object" << std::endl;
  return os;
}

} //namespace MueLu

#endif //ifndef MUELU_SAPFACTORY_HPP
