#ifndef MUELU_RAPFACTORY_HPP
#define MUELU_RAPFACTORY_HPP

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

  public:
    //@{ Constructors/Destructors.
  /*RAPFactory(RCP<FactoryBase> PFact = Teuchos::null, RCP<FactoryBase> RFact = Teuchos::null)
    : PFact_(PFact), RFact_(RFact),
      implicitTranspose_(false) {}*/

    RAPFactory(RCP<FactoryBase> PFact = Teuchos::null, RCP<FactoryBase> RFact = Teuchos::null, RCP<FactoryBase> AFact = Teuchos::null)
    : PFact_(PFact), RFact_(RFact), AFact_(AFact), implicitTranspose_(false) {}

    virtual ~RAPFactory() {}
    //@}

    //! Input
    //@{

    void DeclareInput(Level &fineLevel, Level &coarseLevel) const {
      fineLevel.DeclareInput("A", AFact_.get());   // AFact per default Teuchos::null -> default factory for this
      coarseLevel.DeclareInput("P",PFact_.get());  // transfer operators (from PRFactory, not from PFactory and RFactory!)
      coarseLevel.DeclareInput("R",RFact_.get()); //TODO: must be request according to (implicitTranspose flag!!!!!
    }

    //@}

    //@{ Build methods.
    void Build(Level &fineLevel, Level &coarseLevel) const {  //FIXME make fineLevel const!!
      std::ostringstream buf; buf << coarseLevel.GetLevelID();
      RCP<Teuchos::Time> timer = rcp(new Teuchos::Time("RAP::Build_"+buf.str()));
      timer->start(true);

      Teuchos::OSTab tab(this->getOStream());
      RCP<Operator> P = coarseLevel.Get< RCP<Operator> >("P", PFact_.get());
      RCP<Operator> A = fineLevel.Get< RCP<Operator> >("A",AFact_.get());

RCP<Teuchos::Time> apTimer = rcp(new Teuchos::Time("RAP::A_times_P_"+buf.str()));
apTimer->start(true);
      RCP<Operator> AP = Utils::TwoMatrixMultiply(A,false,P,false);
apTimer->stop();
MemUtils::ReportTimeAndMemory(*apTimer, *(P->getRowMap()->getComm()));
      //std::string filename="AP.dat";
      //Utils::Write(filename,AP);

      Monitor m(*this, "Computing Ac = RAP");

      RCP<Operator> RAP;
      if (implicitTranspose_) {
        //RCP<Operator> RA = Utils::TwoMatrixMultiply(P,true,A,false);
        //filename = "PtA.dat";
        //Utils::Write(filename,AP);

        // avoid implicitTranspose for Epetra, since EpetraExt matrix-matrix multiplication
        // with implicit transpose flags has bugs. This will hopefully be fixed, soon. (see bug #5363)
        //if(RAP->getRangeMap()->lib() == Xpetra::UseEpetra)
          GetOStream(Warnings0, 0) << "The implicitTranspose_ flag within RAPFactory for Epetra in parallel produces wrong results" << std::endl;

        RAP = Utils::TwoMatrixMultiply(P,true,AP,false);
      } else {
        RCP<Operator> R = coarseLevel.Get< RCP<Operator> >("R", RFact_.get());
RCP<Teuchos::Time> rapTimer = rcp(new Teuchos::Time("RAP::R_times_AP_"+buf.str()));
rapTimer->start(true);
        RAP = Utils::TwoMatrixMultiply(R,false,AP,false);
rapTimer->stop();
MemUtils::ReportTimeAndMemory(*rapTimer, *(P->getRowMap()->getComm()));

      }
      
      coarseLevel.Set("A", RAP, this);

      timer->stop();
      MemUtils::ReportTimeAndMemory(*timer, *(P->getRowMap()->getComm()));

      GetOStream(Statistics0, 0) << "Ac: # global rows = " << RAP->getGlobalNumRows() << ", estim. global nnz = " << RAP->getGlobalNumEntries() << std::endl;
    }
    //@}

    void SetImplicitTranspose(bool const &implicit) {
      implicitTranspose_ = implicit;
    }

private:
  //! P Factory
  RCP<FactoryBase> PFact_;

  //! R Factory
  RCP<FactoryBase> RFact_;
  
  //! A Factory
  RCP<FactoryBase> AFact_;

  bool implicitTranspose_;

}; //class RAPFactory

} //namespace MueLu

#define MUELU_RAPFACTORY_SHORT

#endif //ifndef MUELU_RAPFACTORY_HPP
