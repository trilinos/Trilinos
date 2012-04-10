#ifndef MUELU_RAPFACTORY_DECL_HPP
#define MUELU_RAPFACTORY_DECL_HPP

#include <Xpetra_Operator_fwd.hpp>

#include "MueLu_ConfigDefs.hpp"
#include "MueLu_TwoLevelFactoryBase.hpp"
#include "MueLu_RAPFactory_fwd.hpp"

#include "MueLu_Level_fwd.hpp"
#include "MueLu_FactoryBase_fwd.hpp"
#include "MueLu_Utilities_fwd.hpp"

namespace MueLu {
  /*!
    @class RAPFactory
    @brief Factory for building coarse matrices.
  */
  template <class Scalar = double, class LocalOrdinal = int, class GlobalOrdinal = LocalOrdinal, class Node = Kokkos::DefaultNode::DefaultNodeType, class LocalMatOps = typename Kokkos::DefaultKernels<void, LocalOrdinal, Node>::SparseOps>
  class RAPFactory : public TwoLevelFactoryBase {
#undef MUELU_RAPFACTORY_SHORT
#include "MueLu_UseShortNames.hpp"

  public:
    //! @name Constructors/Destructors.
    //@{

    RAPFactory(RCP<const FactoryBase> PFact = Teuchos::null, RCP<const FactoryBase> RFact = Teuchos::null, RCP<const FactoryBase> AFact = Teuchos::null);

    virtual ~RAPFactory();
    //@}

    //! @name Input
    //@{

    void DeclareInput(Level &fineLevel, Level &coarseLevel) const;

    //@}

    //! @name Build methods.
    //@{
    void Build(Level &fineLevel, Level &coarseLevel) const;
    //@}

    //! @name Handling of user-defined transfer factories
    //@{

    //! Indicate that the restriction operator action should be implicitly defined by the transpose of the prolongator.
    void SetImplicitTranspose(bool const &implicit);

    /*! @brief Add transfer factory in the end of list of transfer factories in RAPFactory.

    Transfer factories are derived from TwoLevelFactoryBase and project some data from the fine level to
    the next coarser level.
    */
    void AddTransferFactory(const RCP<const FactoryBase>& factory);

    // TODO add a function to remove a specific transfer factory?

    //! Returns number of transfer factories.
    size_t NumTransferFactories() const;

    //! Returns factory that generates P.
    RCP<const FactoryBase> GetPFactory() const;

    //! Returns factory that generates R.
    RCP<const FactoryBase> GetRFactory() const;

    //@}
  private:
    //! @name internal Build methods.
    //@{
    RCP<Operator> BuildRAPExplicit(const RCP<Operator>& R, const RCP<Operator>& A, const RCP<Operator>& P, int levelID) const; //TODO: remove RCP for input args.
    RCP<Operator> BuildRAPImplicit(const RCP<Operator>& A, const RCP<Operator>& P, int levelID) const;
    RCP<Operator> BuildRAPBlock   (const RCP<Operator>& R, const RCP<Operator>& A, const RCP<Operator>& P, int levelID) const;
    //@}

    //! P Factory
    RCP<const FactoryBase> PFact_;

    //! R Factory
    RCP<const FactoryBase> RFact_;
  
    //! A Factory
    RCP<const FactoryBase> AFact_;

    //! list of user-defined transfer Factories
    std::vector<RCP<const FactoryBase> > TransferFacts_;

    //! If true, the action of the restriction operator action is implicitly defined by the transpose of the prolongator.
    bool implicitTranspose_;

  }; //class RAPFactory

} //namespace MueLu

#define MUELU_RAPFACTORY_SHORT
#endif // MUELU_RAPFACTORY_DECL_HPP
