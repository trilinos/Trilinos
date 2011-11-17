#ifndef MUELU_COALESCEDROPFACTORY_DECL_HPP
#define MUELU_COALESCEDROPFACTORY_DECL_HPP

#include <Xpetra_MapFactory.hpp>
#include <Xpetra_Operator.hpp>
#include <Xpetra_CrsGraphFactory.hpp>

#include "MueLu_ConfigDefs.hpp"
#include "MueLu_SingleLevelFactoryBase.hpp"
#include "MueLu_CoalesceDropFactory_fwd.hpp"

#include "MueLu_Level_fwd.hpp"
#include "MueLu_Graph_fwd.hpp"
#include "MueLu_PreDropFunctionBaseClass_fwd.hpp"

namespace MueLu {

  /*!
    @class CoalesceDropFactory
    @brief Factory for creating a graph base on a given matrix.

    Factory for creating graphs from matrices with entries selectively dropped.
  
    - TODO This factory is very incomplete.
    - TODO The Build method simply builds the matrix graph with no dropping.
  */

  template <class Scalar = double, class LocalOrdinal = int, class GlobalOrdinal = LocalOrdinal, class Node = Kokkos::DefaultNode::DefaultNodeType, class LocalMatOps = typename Kokkos::DefaultKernels<void,LocalOrdinal,Node>::SparseOps> //TODO: or BlockSparseOp ?
  class CoalesceDropFactory : public SingleLevelFactoryBase {
#undef MUELU_COALESCEDROPFACTORY_SHORT
#include "MueLu_UseShortNames.hpp"

  public:

    //! @name Constructors/Destructors.
    //@{

    //! Constructor
    CoalesceDropFactory(RCP<const FactoryBase> AFact = Teuchos::null, RCP<const FactoryBase> nullspaceFact = Teuchos::null);

    //! Destructor
    virtual ~CoalesceDropFactory() { }

    //@}

    //! Input
    //@{

    void DeclareInput(Level &currentLevel) const;

    /// set fixed block size
    void SetFixedBlockSize(LocalOrdinal blksize);

    /// set predrop function
    void SetPreDropFunction(const RCP<MueLu::PreDropFunctionBaseClass<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps> > &predrop);

    // todo: method that takes a block map...

    //@}

    void Build(Level &currentLevel) const; // Build

    void Amalgamate(const RCP<Operator>& A, const LocalOrdinal blocksize, RCP<Graph>& graph) const; // Amalgamate

  private:
    //! A Factory
    RCP<const FactoryBase> AFact_;

    //! nullspace factory
    //! The nullspace dimension is necessary for setting the block size in amalgamation routine
    RCP<const FactoryBase> nullspaceFact_;

    /// blocksize for fixed blocksize setup
    LocalOrdinal blksize_;

    /// are we doing fixed or variable blocks
    bool fixedBlkSize_;

    /// pre-drop function
    RCP<PreDropFunctionBaseClass> predrop_;

  }; //class CoalesceDropFactory

} //namespace MueLu

#define MUELU_COALESCEDROPFACTORY_SHORT
#endif // MUELU_COALESCEDROPFACTORY_DECL_HPP
