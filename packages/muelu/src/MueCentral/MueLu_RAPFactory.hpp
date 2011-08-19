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

  private:
    bool implicitTranspose_;

  public:
    //@{ Constructors/Destructors.
    RAPFactory() : implicitTranspose_(false) {}

    virtual ~RAPFactory() {}
    //@}

    //@{ Build methods.
    bool Build(Level &fineLevel, Level &coarseLevel) const {  //FIXME make fineLevel const!!

      Teuchos::RCP<Teuchos::Time> timer = rcp(new Teuchos::Time("RAP::Build"));
      timer->start(true);

      Teuchos::OSTab tab(this->getOStream());
      //MueLu_cout(Teuchos::VERB_LOW) << "call the Epetra matrix-matrix multiply here" << std::endl;
      RCP<Operator> P = coarseLevel.template Get< Teuchos::RCP<Operator> >("P");
      RCP<Operator> A = fineLevel.template Get< Teuchos::RCP<Operator> >("A");
Teuchos::RCP<Teuchos::Time> apTimer = rcp(new Teuchos::Time("RAP::A_times_P"));
apTimer->start(true);
      RCP<Operator> AP = Utils::TwoMatrixMultiply(A,false,P,false);
apTimer->stop();
MemUtils::ReportTimeAndMemory(*apTimer, *(P->getRowMap()->getComm()));
      //std::string filename="AP.dat";
      //Utils::Write(filename,AP);

      if (implicitTranspose_) {
        //RCP<Operator> RA = Utils::TwoMatrixMultiply(P,true,A,false);
        //filename = "PtA.dat";
        //Utils::Write(filename,AP);
        RCP<Operator> RAP = Utils::TwoMatrixMultiply(P,true,AP,false);
        coarseLevel.Set("A",RAP);
      } else {
        RCP<Operator> R = coarseLevel.template Get< Teuchos::RCP<Operator> >("R");
Teuchos::RCP<Teuchos::Time> rapTimer = rcp(new Teuchos::Time("RAP::R_times_AP"));
rapTimer->start(true);
        RCP<Operator> RAP = Utils::TwoMatrixMultiply(R,false,AP,false);
rapTimer->stop();
MemUtils::ReportTimeAndMemory(*rapTimer, *(P->getRowMap()->getComm()));
        coarseLevel.Set("A",RAP);
      }

      timer->stop();
      MemUtils::ReportTimeAndMemory(*timer, *(P->getRowMap()->getComm()));

      return true;
    }
    //@}

    void SetImplicitTranspose(bool const &implicit) {
      implicitTranspose_ = implicit;
    }

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
