#ifndef MUELU_TENTATIVEPFACTORY_HPP
#define MUELU_TENTATIVEPFACTORY_HPP

#include <iostream>
#include "MueLu_OperatorFactory.hpp"
#include "MueLu_Level.hpp"

namespace MueLu {

  /*!
    @class TentativePFactory class.
    @brief Factory for building tentative prolongator.
  */

template<class ScalarType, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
class TentativePFactory : public OperatorFactory<ScalarType,LocalOrdinal,GlobalOrdinal,Node, LocalMatOps> {
#include "MueLu_UseShortNames.hpp"

  private:

#ifdef MUELU_NOT_IMPLEMENTED
    CoalesceFact_;
    AggFact_;
#endif //MUELU_NOT_IMPLEMENTED
    //! use QR decomposition for improving nullspace information per default
    bool QR_;
    bool reUseGraph_;
    bool reUseAggregates_;

  public:
    //! @name Constructors/Destructors.
    //@{

    //! Constructor.
    TentativePFactory() : QR_(false),reUseGraph_(false),reUseAggregates_(false) {
      Teuchos::OSTab tab(this->out_);
      MueLu_cout(Teuchos::VERB_HIGH) << "TentativePFactory: Instantiating a new factory" << std::endl;
    }

    //! Destructor.
    virtual ~TentativePFactory() {}
    //@}

    //! @name Set methods.
    //@{
#ifdef MUELU_NOT_IMPLEMENTED
    void TentativeWithQR() {}
#endif //MUELU_NOT_IMPLEMENTED
    bool ReUseAggregates() {
      return reUseAggregates_;
    }

    void ReUseAggregates(bool const ToF) {
      reUseAggregates_ = ToF;
    }

    bool ReUseGraph() {
      return reUseGraph_;
    }

    void ReUseGraph(bool ToF) {
      reUseGraph_ = ToF;
    }
    //@}

    //! @name Build methods.
    //@{
    bool Build(Level & fineLevel, Level & coarseLevel) {
      Teuchos::OSTab tab(this->out_); MueLu_cout(Teuchos::VERB_HIGH) << "TentativePFactory: Building a restriction operator" << std::endl; return true;
    }
    //@}

    //! @name Static methods.
    //@{
    //TODO this signature does not match MueMat
    static RCP<Operator> MakeTentative(Level const &currentLevel)
    {
      Teuchos::RCP< Operator > Op = currentLevel.GetA();
      GO nFineDofs = Op->getGlobalNumRows();
      GO nCoarseDofs = nFineDofs/3;
      Teuchos::RCP< Operator > Ptent = Teuchos::rcp( new CrsOperator(Op->getRowMap(), 2) );
      std::vector<SC> Values(1);
      Values[0] = 1.0/sqrt(3.0);
      std::vector<LO> Indices(1);
      Teuchos::ArrayView<SC> av(&Values[0],1);
      Teuchos::ArrayView<GO> iv(&Indices[0],1);
      for (int i=0; i<nFineDofs; i++) {
        Indices[0] = i / 3;
        Ptent->insertGlobalValues(i,iv,av);
      }
      //TODO replace this with a factory or something else
      RCP<Map> domainMap = rcp( new Cthulhu::EpetraMap(nCoarseDofs,Op->getRowMap()->getIndexBase(),Op->getRowMap()->getComm()) );
      Ptent->fillComplete(domainMap, Op->getRowMap());
      //MatrixPrint(Op);
      return Ptent;
    } //MakeTentative()

#ifdef MUELU_NOT_IMPLEMENTED //TODO
    static BuildAggregates() {};
    static MakeNoQRTentative() {};
#endif //MUELU_NOT_IMPLEMENTED
    //@}

}; //class TentativePFactory

} //namespace MueLu

#define MUELU_TENTATIVEPFACTORY_SHORT

#endif //ifndef MUELU_TENTATIVEPFACTORY_HPP
