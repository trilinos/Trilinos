#ifndef MUELU_SALEVEL_HPP
#define MUELU_SALEVEL_HPP

#include "Teuchos_RefCountPtr.hpp"
#include "MueLu_Level.hpp"

namespace MueLu {
/*! @class SaLevel
    @brief Data associated with one level of a smoothed aggregation method

     In addition to standard AMG operators, SaLevel also provides
     a near null space associated with the discretization operator.
*/
template<class SC,class LO, class GO, class NO>
class SaLevel : public Level<SC,LO,GO,NO> {

  typedef Tpetra::CrsMatrix<SC,LO,GO,NO> Operator;
  typedef Tpetra::Vector<SC,LO,GO,NO> Vector;

  template<class AA, class BB, class CC, class DD>
  friend inline std::ostream& operator<<(std::ostream &os, SaLevel<AA,BB,CC,DD> const & level);

  private:
    Teuchos::RCP<Vector>   nullSpace_;
    //aggInfo_;
    Teuchos::RCP<Operator> Ptent_;
    Teuchos::RCP<Operator> Rtent_;

  public:

  //! Constructor
  SaLevel() {
    //Teuchos::OSTab tab(this->out_);
    Teuchos::OSTab newTab = this->getOSTab();
    *(this->out_) << "SaLevel: Instantiating an SaLevel" << std::endl;
  }

  //! Destructor
  virtual ~SaLevel() {}

  //! Set nullspace.
  void SetNullSpace(Teuchos::RCP<Vector> const &nullspace) {
    nullSpace_ = nullspace;
  }

  //! Get nullspace.
  Teuchos::RCP<Vector> GetNullSpace() {
    return nullSpace_;
  }

  //! SetPtent
  void SetPtent(Teuchos::RCP<Operator> Ptent) {
    Ptent_ = Ptent;
  }

  //! SetRtent
  void SetRtent(Teuchos::RCP<Operator> Rtent) {
    Rtent_ = Rtent;
  }

  //! GetPtent
  Teuchos::RCP<Operator> GetPtent() {
    return Ptent_;
  }

  //! GetRtent
  Teuchos::RCP<Operator> GetRtent() {
    return Rtent_;
  }

/* //TODO
     function SetAggregates(this, AggInfo)
     function [z] = GetAggregates(this)
*/

}; //class SaLevel

template<class SC,class LO, class GO, class NO>
std::ostream& operator<<(std::ostream &os, SaLevel<SC,LO,GO,NO> const & level)
{
  os << level.nullSpace_;
  os << level.Ptent_;
  os << level.Rtent_;
  return os;
}

} //namespace MueLu


#endif //ifndef MUELU_SALEVEL_HPP
