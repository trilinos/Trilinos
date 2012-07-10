#ifndef MUELU_MULTIVECTORTRANSFER_FACTORY_DECL_HPP
#define MUELU_MULTIVECTORTRANSFER_FACTORY_DECL_HPP

#include "MueLu_ConfigDefs.hpp"
#include "MueLu_TwoLevelFactoryBase.hpp"
#include "Xpetra_MultiVector_fwd.hpp"
#include "Xpetra_MultiVectorFactory_fwd.hpp"
#include "Xpetra_Operator.hpp"
#include "MueLu_MultiVectorTransferFactory_fwd.hpp"

namespace MueLu {

  /*!
    @class MultiVectorTransferFactory class.
    @brief Class for restricting a MultiVector from a finer to a coarser level.

    This is to be used in conjunction with Muelu::RAPFactory::AddTransferFactory().
  */

  template <class Scalar = double, class LocalOrdinal = int, class GlobalOrdinal = LocalOrdinal, class Node = Kokkos::DefaultNode::DefaultNodeType, class LocalMatOps = typename Kokkos::DefaultKernels<void,LocalOrdinal,Node>::SparseOps>
  class MultiVectorTransferFactory : public TwoLevelFactoryBase {
#undef MUELU_MULTIVECTORTRANSFERFACTORY_SHORT
#include "MueLu_UseShortNames.hpp"

  public:
    //! @name Constructors/Destructors.
    //@{

    /*! @brief Constructor.

       @param vectorName The name of the quantity to be restricted.
       @param restrictionName The name of the restriction Operator.

       The operator associated with <tt>projectionName</tt> will be applied to the MultiVector associated with
       <tt>vectorName</tt>.
    */
    MultiVectorTransferFactory(std::string const & vectorName, std::string const & restrictionName, RCP<const FactoryBase> const &restrictionFact=Teuchos::null);

    //! Destructor.
    virtual ~MultiVectorTransferFactory();

    //@}

    //! @name Input
    //@{

    /*! @brief Specifies the data that this class needs, and the factories that generate that data.

        If the Build method of this class requires some data, but the generating factory is not specified in DeclareInput, then this class
        will fall back to the settings in FactoryManager.
    */
    void DeclareInput(Level &finelevel, Level &coarseLevel) const;

    //@}

    //! @name Build methods.
    //@{

    //! Build an object with this factory.
    void Build(Level & fineLevel, Level &coarseLevel) const;

    //@}

  private:
    //! name of MultiVector to be transfered.
    std::string vectorName_;
    //! name of Operator that will do the transfering.
    std::string restrictionName_;
    //! factory that generates the transfer Operator.
    RCP<const FactoryBase> restrictionFact_;

    static ArrayRCP<SC> expandCoordinates(ArrayRCP<SC> coord, LocalOrdinal blksize);

  }; // class MultiVectorTransferFactory

} // namespace MueLu

#define MUELU_MULTIVECTORTRANSFERFACTORY_SHORT
#endif // MUELU_MULTIVECTORTRANSFER_FACTORY_DECL_HPP
