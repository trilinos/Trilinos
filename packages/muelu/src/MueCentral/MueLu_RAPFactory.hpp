#ifndef MUELU_RAPFACTORY_HPP
#define MUELU_RAPFACTORY_HPP

#include <iostream>

#include "MueLu_ConfigDefs.hpp"
#include "MueLu_TwoLevelFactoryBase.hpp"
#include "MueLu_Level.hpp"
#include "MueLu_Exceptions.hpp"
#include "MueLu_Utilities.hpp"

namespace MueLu {
/*!
  @class RAPFactory class.
  @brief Factory for building coarse matrices.
*/
template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
class RAPFactory : public TwoLevelFactoryBase {

#include "MueLu_UseShortNames.hpp"


  //JG to JJH: use Teuchos::Describable instead ?
  template<class AA, class BB, class CC, class DD, class EE>
  inline friend std::ostream& operator<<(std::ostream& os, RAPFactory<AA,BB,CC,DD,EE> &factory);

  public:
    //@{ Constructors/Destructors.
  /*RAPFactory(RCP<FactoryBase> PFact = Teuchos::null, RCP<FactoryBase> RFact = Teuchos::null)
    : PFact_(PFact), RFact_(RFact),
      implicitTranspose_(false) {}*/

    RAPFactory(RCP<FactoryBase> PRFact = Teuchos::null, RCP<FactoryBase> AFact = Teuchos::null)
    : PRFact_(PRFact), AFact_(AFact), implicitTranspose_(false) {}

    virtual ~RAPFactory() {}
    //@}

    //! Input
    //@{

    void DeclareInput(Level &fineLevel, Level &coarseLevel) const {
      fineLevel.Request("A", AFact_.get());     // AFact per default Teuchos::null -> default factory for this
      coarseLevel.Request("P",PRFact_.get());   // transfer operators (from PRFactory, not from PFactory and RFactory!)
      coarseLevel.Request("R",PRFact_.get());
    }

    //@}

    //@{ Build methods.
    bool Build(Level &fineLevel, Level &coarseLevel) const {  //FIXME make fineLevel const!!

      std::ostringstream buf; buf << coarseLevel.GetLevelID();
      RCP<Teuchos::Time> timer = rcp(new Teuchos::Time("RAP::Build_"+buf.str()));
      timer->start(true);

      Teuchos::OSTab tab(this->getOStream());
      //MueLu_cout(Teuchos::VERB_LOW) << "call the Epetra matrix-matrix multiply here" << std::endl;
      RCP<Operator> P = coarseLevel.Get< RCP<Operator> >("P", PRFact_.get());

      RCP<Operator> A = fineLevel.Get< RCP<Operator> >("A",AFact_.get());

RCP<Teuchos::Time> apTimer = rcp(new Teuchos::Time("RAP::A_times_P_"+buf.str()));
apTimer->start(true);
      RCP<Operator> AP = Utils::TwoMatrixMultiply(A,false,P,false);
apTimer->stop();
MemUtils::ReportTimeAndMemory(*apTimer, *(P->getRowMap()->getComm()));
      //std::string filename="AP.dat";
      //Utils::Write(filename,AP);

      RCP<Operator> RAP;
      if (implicitTranspose_) {
        //RCP<Operator> RA = Utils::TwoMatrixMultiply(P,true,A,false);
        //filename = "PtA.dat";
        //Utils::Write(filename,AP);
        RAP = Utils::TwoMatrixMultiply(P,true,AP,false);
      } else {
        RCP<Operator> R = coarseLevel.Get< RCP<Operator> >("R", PRFact_.get());
RCP<Teuchos::Time> rapTimer = rcp(new Teuchos::Time("RAP::R_times_AP_"+buf.str()));
rapTimer->start(true);
        RAP = Utils::TwoMatrixMultiply(R,false,AP,false);
rapTimer->stop();
MemUtils::ReportTimeAndMemory(*rapTimer, *(P->getRowMap()->getComm()));

      }
      
      coarseLevel.Set("A", RAP, AFact_.get());

      timer->stop();
      MemUtils::ReportTimeAndMemory(*timer, *(P->getRowMap()->getComm()));

      fineLevel.Release("A",AFact_.get());
      coarseLevel.Release("P",PRFact_.get());
      coarseLevel.Release("R",PRFact_.get());

      return true;
    }
    //@}

    void SetImplicitTranspose(bool const &implicit) {
      implicitTranspose_ = implicit;
    }

    void SetPRFact(RCP<FactoryBase> PRFact) {
      PRFact_ = PRFact;
    }

private:
  //! P Factory
  //RCP<FactoryBase> PFact_;

  //! R Factory
  //RCP<FactoryBase> RFact_;

  //! PR Factory
  RCP<FactoryBase> PRFact_;
  
  //! A Factory
  RCP<FactoryBase> AFact_;

  bool implicitTranspose_;

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
