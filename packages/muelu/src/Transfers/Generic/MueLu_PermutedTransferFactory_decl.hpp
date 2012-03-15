#ifndef MUELU_PERMUTEDTRANSFER_FACTORY_DECL_HPP
#define MUELU_PERMUTEDTRANSFER_FACTORY_DECL_HPP

#include <Xpetra_Operator_fwd.hpp>
#include "Xpetra_MultiVector_fwd.hpp"
#include "Xpetra_MultiVectorFactory_fwd.hpp"

#include "MueLu_ConfigDefs.hpp"
#include "MueLu_TwoLevelFactoryBase.hpp"
#include "MueLu_PermutedTransferFactory_fwd.hpp"
#include "MueLu_Utilities_fwd.hpp"
#include "MueLu_Types.hpp"

namespace MueLu {

  /*!
    @class PermutedTransferFactory class.
    @brief Applies permutation to grid transfer operators.
    @ingroup MueLuTransferClasses
  */

  template <class Scalar = double, class LocalOrdinal = int, class GlobalOrdinal = LocalOrdinal, class Node = Kokkos::DefaultNode::DefaultNodeType, class LocalMatOps = typename Kokkos::DefaultKernels<void,LocalOrdinal,Node>::SparseOps>
  class PermutedTransferFactory : public TwoLevelFactoryBase {
#undef MUELU_PERMUTEDTRANSFER_FACTORY_SHORT
#include "MueLu_UseShortNames.hpp"

  public:
    //! @name Constructors/Destructors.
    //@{

    //! Constructor.
    PermutedTransferFactory(RCP<FactoryBase> repartitionFact=Teuchos::null,
                             RCP<FactoryBase> initialAFact=Teuchos::null,
                             RCP<FactoryBase> initialTransferFact=Teuchos::null,
                             TransferType PorR = MueLu::INTERPOLATION,
                             RCP<FactoryBase> nullspaceFact=Teuchos::null,
                             RCP<FactoryBase> coordinateFact=Teuchos::null);

    //! Destructor.
    virtual ~PermutedTransferFactory();

    //@}

    //! @name Input
    //@{

    /*! @brief Specifies the data that this class needs, and the factories that generate that data.

        If the Build method of this class requires some data, but the generating factory is not specified in DeclareInput, then this class
        will fall back to the settings in FactoryManager.
    */
    void DeclareInput(Level &fineLevel, Level &coarseLevel) const;

    //@}

    //! @name Build methods.
    //@{

    //! Build an object with this factory.
    void Build(Level &fineLevel, Level &coarseLevel) const;

    //@}

  private:
    //! Factory that builds the permutation matrix.
    RCP<FactoryBase> repartitionFact_;
    //! Factory that builds the A matrix.
    RCP<FactoryBase> initialAFact_;
    //! Factory that builds the unpermuted grid transfer operator.
    RCP<FactoryBase> initialTransferFact_;
    //! Indicate that the transfer factory is for interpolation or restriction.
    TransferType     PorR_;
    //! Factory that builds the unpermuted nullspace.
    RCP<FactoryBase> nullspaceFact_;
    //! Factory that builds the unpermuted coordinates.
    RCP<FactoryBase> coordinateFact_;

  }; // class PermutedTransferFactory

} // namespace MueLu

#define MUELU_PERMUTEDTRANSFER_FACTORY_SHORT
#endif // MUELU_PERMUTEDTRANSFER_FACTORY_DECL_HPP
