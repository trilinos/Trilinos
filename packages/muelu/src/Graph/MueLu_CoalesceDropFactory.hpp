#ifndef MUELU_COALESCEDROPFACTORY_HPP
#define MUELU_COALESCEDROPFACTORY_HPP

#include "Xpetra_Operator.hpp"

#include "MueLu_ConfigDefs.hpp"
#include "MueLu_SingleLevelFactoryBase.hpp"
#include "MueLu_Level.hpp"
#include "MueLu_Graph.hpp"

namespace MueLu {

static const std::string color_esc = "\x1b[";
static const std::string color_std = "39;49;00m";
static const std::string color_purple = "35m";

  /*!
    @class CoalesceDropFactory
    @brief Factory for creating a graph base on a given matrix.

    Factory for creating graphs from matrices with entries selectively dropped.
  
    - TODO This factory is very incomplete.
    - TODO The Build method simply builds the matrix graph with no dropping.
  */

  template <class Scalar = double, class LocalOrdinal = int, class GlobalOrdinal = LocalOrdinal, class Node = Kokkos::DefaultNode::DefaultNodeType, class LocalMatOps = typename Kokkos::DefaultKernels<void,LocalOrdinal,Node>::SparseOps> //TODO: or BlockSparseOp ?
  class CoalesceDropFactory : public SingleLevelFactoryBase {

#include "MueLu_UseShortNames.hpp"

  public:

    //! @name Constructors/Destructors.
    //@{

    //! Constructor
    CoalesceDropFactory(RCP<const FactoryBase> AFact = Teuchos::null)
      : AFact_(AFact), fixedBlkSize_(true)
    { }

    //! Destructor
    virtual ~CoalesceDropFactory() {}
    //@}

    //! Input
    //@{

    void DeclareInput(Level &currentLevel) const {
      currentLevel.DeclareInput("A", AFact_.get());
    }

    void SetFixedBlockSize(GlobalOrdinal blksize) {
    	blksize_ = blksize;
    	fixedBlkSize_ = true;
    	GetOStream(Low, 0) << color_esc << color_purple << "CoalesceDropFactory::SetFixedBlockSize()" << color_esc << color_std << endl;
    }

    // todo: method that takes a block map...

    //@}

    void Build(Level &currentLevel) const {
      RCP<Operator> A = currentLevel.Get< RCP<Operator> >("A", AFact_.get());

      RCP<Graph> graph = rcp(new Graph(A->getCrsGraph(), "Graph of A"));

      currentLevel.Set("Graph", graph, this);

    } // Build

  private:
    //! A Factory
    RCP<const FactoryBase> AFact_;

    /// blocksize for fixed blocksize setup
    GO blksize_;

    /// are we doing fixed or variabled blocks
    bool fixedBlkSize_;

  }; //class CoalesceDropFactory

} //namespace MueLu

#define MUELU_COALESCEDROPFACTORY_SHORT

#endif //ifndef MUELU_COALESCEDROPFACTORY_HPP
