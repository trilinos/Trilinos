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
class TentativePFactory : public PFactory<ScalarType,LocalOrdinal,GlobalOrdinal,Node, LocalMatOps> {
#include "MueLu_UseShortNames.hpp"

  private:

    /*
    TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO 
    coalesceFact_;
    aggregationFact_;
    TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO 
    */
    //! use QR decomposition for improving nullspace information per default
    bool QR_;

  public:
    //! @name Constructors/Destructors.
    //@{

    //! Constructor.
    TentativePFactory(/*TODO*/ /*CoalesceFact, AggregationFact*/ /*TODO*/) : QR_(false) {
      Teuchos::OSTab tab(this->out_);
      MueLu_cout(Teuchos::VERB_HIGH) << "TentativePFactory: Instantiating a new factory" << std::endl;
      /*
        if (CoalesceFact != Teuchos::null) coalesceFact_ = CoalesceFact;
        else                               coalesceFact_ = rcp(new CoalesceDropFactory());
        if (AggregationFact != Teuchos::null) aggregationFact_ = AggregationFact;
        else                                  aggregationFact_ = rcp(new AggregationFactory());
      */
    }

    //! Destructor.
    virtual ~TentativePFactory() {}
    //@}

    //! @name Set methods.
    //@{

    void TentativeWithQR(bool value) {
      QR_ = value;
    }

    //@}

    //! @name Build methods.
    //@{
    bool BuildP(Level & fineLevel, Level & coarseLevel) {
      Teuchos::OSTab tab(this->out_); MueLu_cout(Teuchos::VERB_HIGH) << "TentativePFactory: Building a restriction operator" << std::endl; return true;
    }
    //@}

    //! @name Static methods.
    //@{
    /*! @brief Make tenative prolongator.
        TODO this signature does not match MueMat
    */
    static RCP<Operator> MakeTentative(Level const &currentLevel)
    {
      Teuchos::RCP< Operator > Op = currentLevel.GetA();
      GO nFineDofs = Op->getGlobalNumRows();
      GO nCoarseDofs = nFineDofs/3;
      if (nCoarseDofs*3 != nFineDofs)
        throw(Exceptions::NotImplemented("MakeTentative: currently #fine DOFS must be a multiple of 3"));
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
#ifdef HAVE_CTHULHU_EPETRA
      RCP<Map> domainMap = rcp( new Cthulhu::EpetraMap(nCoarseDofs,Op->getRowMap()->getIndexBase(),Op->getRowMap()->getComm()) );

      Ptent->fillComplete(domainMap, Op->getRowMap());
#else
#warning
#endif

      //MatrixPrint(Op);
      return Ptent;
    } //MakeTentative()

    static void BuildAggregates() {
      throw(Exceptions::NotImplemented());
    }
    static void MakeNoQRTentative() {
      throw(Exceptions::NotImplemented());
    }
    //@}

}; //class TentativePFactory

} //namespace MueLu

#define MUELU_TENTATIVEPFACTORY_SHORT

#endif //ifndef MUELU_TENTATIVEPFACTORY_HPP
