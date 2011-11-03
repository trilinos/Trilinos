#ifndef MUELU_NULLSPACEFACTORY_DECL_HPP
#define MUELU_NULLSPACEFACTORY_DECL_HPP

#include "MueLu_ConfigDefs.hpp"

#include "Xpetra_Operator.hpp"
#include "Xpetra_MultiVectorFactory.hpp"

#include "MueLu_SingleLevelFactoryBase.hpp"
#include "MueLu_Level.hpp"

namespace MueLu {

  template <class Scalar = double, class LocalOrdinal = int, class GlobalOrdinal = LocalOrdinal, class Node = Kokkos::DefaultNode::DefaultNodeType, class LocalMatOps = typename Kokkos::DefaultKernels<void,LocalOrdinal,Node>::SparseOps> //TODO: or BlockSparseOp ?
  class NullspaceFactory : public SingleLevelFactoryBase {

#include "MueLu_UseShortNames.hpp"

  public:

    //! @name Constructors/Destructors.
    //@{

    //! Constructor
    NullspaceFactory(RCP<const FactoryBase> AFact = Teuchos::null)
    ;

    //! Destructor
    virtual ~NullspaceFactory() ;

    //@}

    //! Input
    //@{

    void DeclareInput(Level &currentLevel) const ;

    //@}

    void Build(Level &currentLevel) const ; // Build

  private:

    //! A Factory
    RCP<const FactoryBase> AFact_;

  }; //class NullspaceFactory

} //namespace MueLu

#define MUELU_NULLSPACEFACTORY_SHORT
#endif // MUELU_NULLSPACEFACTORY_DECL_HPP
