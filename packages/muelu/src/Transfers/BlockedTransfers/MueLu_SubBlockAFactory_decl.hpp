/*
 * MueLu_SubBlockAFactory_decl.hpp
 *
 *  Created on: 02.01.2012
 *      Author: tobias
 */

#ifndef MUELU_SUBBLOCKAFACTORY_DECL_HPP_
#define MUELU_SUBBLOCKAFACTORY_DECL_HPP_

#include "Xpetra_Map_fwd.hpp"
#include "Xpetra_StridedMap_fwd.hpp"

#include "MueLu_ConfigDefs.hpp"
#include "MueLu_SingleLevelFactoryBase.hpp"
#include "MueLu_SubBlockAFactory_fwd.hpp"


namespace MueLu {

  /*!
    @class SubBlockAFactory class.
    @brief Factory for building a thresholded operator.

    This is a very simple class to access a single matrix block in a blocked operator A.

    Example
    \code
    Teuchos::RCP<Xpetra::BlockedCrsOperator<Scalar,LO,GO,Node> > bOp = Teuchos::rcp(new Xpetra::BlockedCrsOperator<Scalar,LO,GO>(mapExtractor,mapExtractor,10));
    // ... let bOp be a 2x2 blocked operator ...
    bOp->fillComplete();

    // define factory for accessing block (0,0) in blocked operator A (assuming that the blocked operator is stored in Level class with NoFactory as generating factory)
    RCP<SubBlockAFactory> A11Fact = Teuchos::rcp(new SubBlockAFactory(MueLu::NoFactory::getRCP(), 0, 0));

    // define factory for accessing block (1,1) in blocked operator A
    RCP<SubBlockAFactory> A22Fact = Teuchos::rcp(new SubBlockAFactory(MueLu::NoFactory::getRCP(), 1, 1));

    RCP<Operator> A11 = level.Get<RCP<Operator> >("A", A11Fact); // extract (0,0) block from blocked operator A
    RCP<Operator> A22 = level.Get<RCP<Operator> >("A", A22Fact); // extract (1,1) block from blocked operator A
    \endcode
  */

  template <class Scalar = double, class LocalOrdinal = int, class GlobalOrdinal = LocalOrdinal, class Node = Kokkos::DefaultNode::DefaultNodeType, class LocalMatOps = typename Kokkos::DefaultKernels<void,LocalOrdinal,Node>::SparseOps>
  class SubBlockAFactory : public SingleLevelFactoryBase {
#undef MUELU_SUBBLOCKAFACTORY_SHORT
    #include "MueLu_UseShortNames.hpp"

  public:
    //! @name Constructors/Destructors.
    //@{

    //! Constructor.
    SubBlockAFactory(Teuchos::RCP<const FactoryBase> Afact, size_t row, size_t col, LocalOrdinal blksize = 1);

    //! Destructor.
    virtual ~SubBlockAFactory();
    //@}

    //! Input
    //@{

    void DeclareInput(Level &currentLevel) const;

    //@}

    //@{
    //! @name Build methods.

    //! Build an object with this factory.
    void Build(Level & currentLevel) const;

    //@}

    // ----------------------------------------------------------------------------------
    // "TEMPORARY" VIEW MECHANISM
    // TODO: the view mechanism should be implemented as in MueMat.
    void         SetFixedBlockSize(LocalOrdinal blksize);
    LocalOrdinal GetFixedBlockSize() const;
    // ----------------------------------------------------------------------------------
  private:
    Teuchos::RCP<const FactoryBase> Afact_;   ///< generating factory of input variable
    const size_t                    row_;     ///< row id
    const size_t                    col_;     ///< column id

    // ----------------------------------------------------------------------------------
    // "TEMPORARY" VIEW MECHANISM
    // TODO: the view mechanism should be implemented as in MueMat.
    LocalOrdinal blksize_;
    // RCP<GOVector> variableBlockSizeInfo_; TODO: should be moved from CoalesceDropFactory to here.
    // ----------------------------------------------------------------------------------
  }; // class SubBlockAFactory

} // namespace MueLu

#define MUELU_SUBBLOCKAFACTORY_SHORT
#endif /* MUELU_SUBBLOCKAFACTORY_DECL_HPP_ */
