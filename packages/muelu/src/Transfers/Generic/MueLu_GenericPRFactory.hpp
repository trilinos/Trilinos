#ifndef MUELU_GENERICPRFACTORY_HPP
#define MUELU_GENERICPRFACTORY_HPP

#include <iostream>

#include "MueLu_ConfigDefs.hpp"
#include "MueLu_PRFactory.hpp"
#include "MueLu_PFactory.hpp"
#include "MueLu_RFactory.hpp"
#include "MueLu_TransPFactory.hpp"
#include "MueLu_SaPFactory.hpp"
#include "MueLu_Exceptions.hpp"

namespace MueLu {

/*!
  @class GenericPRFactory
  @brief Concrete implementation of PRFactory.

  This class combines a prolongation operator object (derived from
  PFactory) and a restriction operator object (derived from RFactory) and
  creates a PRFactory, that handles both the prolongator and the
  restrictor.
*/

template <class Scalar = double, class LocalOrdinal = int, class GlobalOrdinal = LocalOrdinal, class Node = Kokkos::DefaultNode::DefaultNodeType, class LocalMatOps = typename Kokkos::DefaultKernels<void,LocalOrdinal,Node>::SparseOps> //TODO: or BlockSparseOp ?
class GenericPRFactory : public PRFactory {
#include "MueLu_UseShortNames.hpp"

  private:

    RCP<PFactory> PFact_;
    RCP<RFactory> RFact_;

  public:
    //! @name Constructors/Destructors.
    //@{

    /* @brief Constructor.

       \param PFact prolongator factory
       \param RFact restriction factory

       Defaults to plain aggregation.
    */
    GenericPRFactory(RCP<PFactory> PFact = Teuchos::null, RCP<RFactory> RFact = Teuchos::null)
    {
      //Teuchos::OSTab tab(this->out_);
      //MueLu_cout(Teuchos::VERB_HIGH) << "GenericPRFactory: Instantiating a new factory" << std::endl;

      if (PFact != Teuchos::null)
        PFact_ = PFact;
      else
        PFact_ = rcp(new SaPFactory());

      if (RFact != Teuchos::null)
        RFact_ = RFact;
      else
        RFact_ = rcp(new TransPFactory());

      PRFactory::reUseAggregates_ = PFact_->ReUseAggregates();
      PRFactory::reUseGraph_      = PFact_->ReUseGraph();

    }

    //! Destructor.
    virtual ~GenericPRFactory() {}
    //@}

    //! Input
    //@{

    void DeclareInput(Level &fineLevel, Level &coarseLevel) const { }

    //@}

    //! @name Build methods.
    //@{

    /*!
      @brief Build method.
    */
    bool Build(Level &fineLevel, Level &coarseLevel) const {

      //FIXME what if map is a block map .... I'm pretty sure maxCoarseSize_ will always be point DOFs
      if (fineLevel.Get< RCP<Operator> >("A")->getRowMap()->getComm()->getRank() == 0)
        std::cout << "warning: if map is blocked, this comparison to maxCoarseSize_ will be wrong!" << std::endl;
      if (fineLevel.Get< RCP<Operator> >("A")->getRowMap()->getGlobalNumElements() <= PRFactory::maxCoarseSize_)
        return false;
      
      //FIXME cache output level here

      RCP<Operator> P = coarseLevel.Get< RCP<Operator> >("P", PFact_);
      coarseLevel.Set< RCP<Operator> >("P", P, this);
      //      PFact_->BuildP(fineLevel,coarseLevel);
      //      if (RFact_ != Teuchos::null) {
        RCP<Operator> R = coarseLevel.Get< RCP<Operator> >("R", RFact_);
        coarseLevel.Set< RCP<Operator> >("R", R, this);
        //        RFact_->BuildR(fineLevel,coarseLevel);
        //}
      //FIXME restore output level here
      return true;
    }
    //@}

    //! @name Set methods.
    //@{

    void ReUseAggregates(bool const &value) {
      PRFactory::reUseAggregates_ = value;
      PFact_->ReUseAggregates(value);
    }

    void ReUseGraph(bool const &value) {
      PRFactory::reUseGraph_ = value;
      PFact_->ReUseGraph(value);
    }
    //@}

    //! @name Get methods.
    //@{

    bool ReUseAggregates() const {
      return PRFactory::reUseAggregates_;
    }

    bool ReUseGraph() const {
      return PRFactory::reUseGraph_;
    }

    //@}

}; //class GenericPRFactory

} //namespace MueLu

#define MUELU_GENERICPRFACTORY_SHORT

#endif //ifndef MUELU_GENERICPRFACTORY_HPP
