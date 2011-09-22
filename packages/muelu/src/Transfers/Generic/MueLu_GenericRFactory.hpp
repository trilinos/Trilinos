/*
 * MueLu_GenericRFactory.hpp
 *
 *  Created on: 20.09.2011
 *      Author: tobias
 */

#ifndef MUELU_GENERICRFACTORY_HPP_
#define MUELU_GENERICRFACTORY_HPP_

#include <iostream>

#include "Xpetra_CrsOperator.hpp"

#include "MueLu_ConfigDefs.hpp"
#include "MueLu_Level.hpp"
#include "MueLu_RFactory.hpp"
#include "MueLu_Exceptions.hpp"
#include "MueLu_Utilities.hpp"

namespace MueLu {

  /*!
    @class GenericRFactory class.
    @brief Factory for building restriction operators using a prolongator factory
  */

  template < class Scalar = double, class LocalOrdinal = int, class GlobalOrdinal = LocalOrdinal, class Node = Kokkos::DefaultNode::DefaultNodeType, class LocalMatOps = typename Kokkos::DefaultKernels<void,LocalOrdinal,Node>::SparseOps >
  class GenericRFactory : public RFactory {

#include "MueLu_UseShortNames.hpp"

    template<class AA, class BB, class CC, class DD, class EE>
    inline friend std::ostream& operator<<(std::ostream& os, GenericRFactory<AA,BB,CC,DD,EE> &factory);

  public:
    //! @name Constructors/Destructors.
    //@{

    //! Constructor.
    GenericRFactory(RCP<PFactory> PFact = Teuchos::null)
      : PFact_(PFact)
    {
        TEST_FOR_EXCEPTION(PFact_==Teuchos::null, MueLu::Exceptions::RuntimeError, "GenericRFactory: no valid PFactory provided.");
    }

    //! Destructor.
    virtual ~GenericRFactory() {}
    //@}

    //! Input
    //@{

    void DeclareInput(Level &fineLevel, Level &coarseLevel) const {
      //coarseLevel.Request("P",PFact_.get());
        // todo: declare input
        bool rmode = PFact_->isRestrictionModeSet();
        PFact_->setRestrictionMode(true);             // set restriction mode
        //PFact_->DeclareInput(fineLevel, coarseLevel); // call DeclareInput explicitely
        coarseLevel.Request("R", PFact_.get());       // we expect the prolongation operator factory to produce "R" as output
                                                      // call declareInput is called within Request call
        PFact_->setRestrictionMode(rmode);            // reset restriciton mode flag
    }

    //@}

    //! @name Build methods.
    //@{

    bool Build(Level & fineLevel, Level & coarseLevel) const {
      return BuildR(fineLevel,coarseLevel);
    }

    bool BuildR(Level & fineLevel, Level & coarseLevel) const {

      std::ostringstream buf; buf << coarseLevel.GetLevelID();
      RCP<Teuchos::Time> timer = rcp(new Teuchos::Time("GenericRFactory::BuildR_"+buf.str()));
      timer->start(true);

      // BuildR
      bool rmode = PFact_->isRestrictionModeSet();
      PFact_->setRestrictionMode(true);     // switch prolongator factory to restriction mode

      //PFact_->Build(fineLevel, coarseLevel);  // call PFactory::Build explicitely

      RCP<Operator> R = coarseLevel.Get<RCP<Operator> >("R",PFact_.get());
      coarseLevel.Release("R",PFact_.get());

      PFact_->setRestrictionMode(rmode);    // reset restriction mode flag

      coarseLevel.Set("R", R, this);

      timer->stop();
      MemUtils::ReportTimeAndMemory(*timer, *(R->getRowMap()->getComm()));

      return true;
    } //BuildR

    //@}

  private:

    //! P Factory
    RCP<PFactory> PFact_;

  }; //class TransPFactory

  //! Friend print function.
  template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  std::ostream& operator<<(std::ostream& os, GenericRFactory<Scalar,LocalOrdinal,GlobalOrdinal,Node, LocalMatOps> &factory) {
    os << "Printing a GenericRFactory object" << std::endl;
    return os;
  }

} //namespace MueLu

#define MUELU_GENERICRFACTORY_SHORT



#endif /* MUELU_GENERICRFACTORY_HPP_ */
