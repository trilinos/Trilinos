#ifndef MUELU_SALEVEL_HPP
#define MUELU_SALEVEL_HPP

#include "Teuchos_RefCountPtr.hpp"
#include "MueLu_Level.hpp"

namespace MueLu {
  using Teuchos::RCP;

  /*! @class SaLevel
    @brief Data associated with one level of a smoothed aggregation method

    In addition to standard AMG operators, SaLevel also provides
    a near null space associated with the discretization operator.
  */
  template<class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  class SaLevel : public Level<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps> {

#include "MueLu_UseShortNames.hpp"

    //JG to JJH: use Teuchos::Describable instead ?
    template<class AA, class BB, class CC, class DD, class EE>
    friend inline std::ostream& operator<<(std::ostream &os, SaLevel<AA,BB,CC,DD,EE> const & level);

  private:
    RCP<Vector>   nullSpace_;
    //aggInfo_;
    RCP<Operator> Ptent_;
    RCP<Operator> Rtent_;

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
    void SetNullSpace(RCP<Vector> const &nullspace) {
      nullSpace_ = nullspace;
    }

    //! Get nullspace.
    RCP<Vector> GetNullSpace() {
      return nullSpace_;
    }

    //! SetPtent
    void SetPtent(RCP<Operator> Ptent) {
      Ptent_ = Ptent;
    }

    //! SetRtent
    void SetRtent(RCP<Operator> Rtent) {
      Rtent_ = Rtent;
    }

    //! GetPtent
    RCP<Operator> GetPtent() {
      return Ptent_;
    }

    //! GetRtent
    RCP<Operator> GetRtent() {
      return Rtent_;
    }

    /* //TODO
       function SetAggregates(this, AggInfo)
       function [z] = GetAggregates(this)
    */

  }; //class SaLevel

  template<class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  std::ostream& operator<<(std::ostream &os, SaLevel<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps> const & level)
  {
    os << level.nullSpace_;
    os << level.Ptent_;
    os << level.Rtent_;
    return os;
  }

} //namespace MueLu

#define MUELU_SALEVEL_SHORT

#endif //ifndef MUELU_SALEVEL_HPP

//JG: TODO: add inline keyword ?
