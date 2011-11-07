#ifndef MUELU_RAPFACTORY_DECL_HPP
#define MUELU_RAPFACTORY_DECL_HPP

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
  template <class Scalar = double, class LocalOrdinal = int, class GlobalOrdinal = LocalOrdinal, class Node = Kokkos::DefaultNode::DefaultNodeType, class LocalMatOps = typename Kokkos::DefaultKernels<void, LocalOrdinal, Node>::SparseOps>
  class RAPFactory : public TwoLevelFactoryBase {

#include "MueLu_UseShortNames.hpp"

  public:
    //@{ Constructors/Destructors.
  /*RAPFactory(RCP<FactoryBase> PFact = Teuchos::null, RCP<FactoryBase> RFact = Teuchos::null)
    : PFact_(PFact), RFact_(RFact),
      implicitTranspose_(false);*/

    RAPFactory(RCP<FactoryBase> PFact = Teuchos::null, RCP<FactoryBase> RFact = Teuchos::null, RCP<FactoryBase> AFact = Teuchos::null);

    virtual ~RAPFactory();
    //@}

    //! Input
    //@{

    void DeclareInput(Level &fineLevel, Level &coarseLevel) const;

    //@}

    //@{ Build methods.
    void Build(Level &fineLevel, Level &coarseLevel) const;
    //@}

    void SetImplicitTranspose(bool const &implicit);

    //@{ Handling of user-defined transfer factories

    //! add transfer factory in the end of list of transfer factories in RAPFactory
    //! Transfer factories are derived from TwoLevelFactoryBase and project some data from the fine level to the next coarser level
    void AddTransferFactory(const RCP<FactoryBase>& factory);

    // TODO add a function to remove a specific transfer factory?

    //! returns number of transfer factories
    size_t NumTransferFactories() const;

    //@}
private:
  //! P Factory
  RCP<FactoryBase> PFact_;

  //! R Factory
  RCP<FactoryBase> RFact_;
  
  //! A Factory
  RCP<FactoryBase> AFact_;

  //! list of user-defined transfer Factories
  std::vector<RCP<FactoryBase> > TransferFacts_;

  bool implicitTranspose_;

}; //class RAPFactory

} //namespace MueLu

#define MUELU_RAPFACTORY_SHORT
#endif // MUELU_RAPFACTORY_DECL_HPP
