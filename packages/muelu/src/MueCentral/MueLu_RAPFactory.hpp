#ifndef MUELU_RAPFACTORY_HPP
#define MUELU_RAPFACTORY_HPP

#include <iostream>
#include "MueLu_TwoLevelFactoryBase.hpp"
#include "MueLu_Exceptions.hpp"
#include "MueLu_Utilities.hpp"

namespace MueLu {
/*!
  @class RAPFactory class.
  @brief Factory for building coarse matrices.
*/
template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
class RAPFactory : public TwoLevelFactoryBase<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps> {

#include "MueLu_UseShortNames.hpp"


  //JG to JJH: use Teuchos::Describable instead ?
  template<class AA, class BB, class CC, class DD, class EE>
  inline friend std::ostream& operator<<(std::ostream& os, RAPFactory<AA,BB,CC,DD,EE> &factory);

  public:
    //@{ Constructors/Destructors.
    RAPFactory() {}

    virtual ~RAPFactory() {}
    //@}

    //@{ Build methods.
    bool Build(Level &fineLevel, Level &coarseLevel) const {  //FIXME make fineLevel const!!

      Teuchos::RCP<Teuchos::Time> timer = rcp(new Teuchos::Time("RAP::Build"));
      timer->start(true);

      Teuchos::OSTab tab(this->getOStream());
      //MueLu_cout(Teuchos::VERB_LOW) << "call the Epetra matrix-matrix multiply here" << std::endl;
      RCP<Operator> P = coarseLevel.GetP();
      RCP<Operator> A = fineLevel.GetA();
      RCP<Operator> AP = Utils::TwoMatrixMultiply(A,P);
      RCP<Operator> R = coarseLevel.GetR();
      RCP<Operator> RAP = Utils::TwoMatrixMultiply(R,AP);
      coarseLevel.SetA(RAP);

      timer->stop();
      MemUtils::ReportTimeAndMemory(*timer, *(RAP->getRowMap()->getComm()));

      return true;
    }
    //@}


}; //class RAPFactory

//! Friend print method.
  template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
std::ostream& operator<<(std::ostream& os, RAPFactory<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps> &factory) {
  os << "Printing RAPFactory object" << std::endl;
  return os;
}

} //namespace MueLu

#define MUELU_RAPFACTORY_SHORT

#endif //ifndef MUELU_RAPFACTORY_HPP
