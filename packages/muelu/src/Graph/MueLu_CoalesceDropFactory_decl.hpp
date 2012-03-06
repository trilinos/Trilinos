#ifndef MUELU_COALESCEDROPFACTORY_DECL_HPP
#define MUELU_COALESCEDROPFACTORY_DECL_HPP

#include <Xpetra_Map_fwd.hpp>
#include <Xpetra_MapFactory_fwd.hpp>
#include <Xpetra_Vector_fwd.hpp>
#include <Xpetra_MultiVector_fwd.hpp>
#include <Xpetra_CrsGraph_fwd.hpp>
#include <Xpetra_CrsGraphFactory.hpp> //TODO
#include <Xpetra_Operator_fwd.hpp>

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

    /// set variable block size
    void SetVariableBlockSize();

    /// set predrop function
    void SetPreDropFunction(const RCP<MueLu::PreDropFunctionBaseClass<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps> > &predrop);

    // todo: method that takes a block map...

    //@}

    void Build(Level &currentLevel) const; // Build

    void Amalgamate(const RCP<Operator>& A, const LocalOrdinal blocksize, RCP<Graph>& graph) const; // Amalgamate

  private:

    //! translate global (row/column) id to global amalgamation block id
    // @param gid (GlobalOrdinal): input global id (row gid or column gid)
    // @param A: input operator (just used to check the maps for validity)
    // @param globalgid2globalamalblockid_vector: Xpetra vector which holds block amalgamation gids for all column gids (vector lives on overlapping column map of A! needed for variable block size)
    // @param blockSize (LocalOrdinal): block size (needed for constant block size)
    GlobalOrdinal GlobalId2GlobalAmalBlockId(GlobalOrdinal gid, const RCP<Operator>& A, const RCP<GOVector>& globalgid2globalamalblockid_vector, LocalOrdinal blockSize) const;

    //! setup amalgamation data
    // This routine fills the private members lobalamalblockid2myrowid_ and lobalamalblockid2globalrowid_
    // These are then used in Graph and Aggregates for unamalgamation
    // @param A: input operator (just used to check the maps for validity)
    // @param globalgid2globalamalblockid_vector: Xpetra vector which holds block amalgamation gids for all column gids (vector lives on overlapping column map of A! needed for variable block size)
    // @param blockSize (LocalOrdinal): block size (needed for constant block size)
    // returns amalgamated map
    const Teuchos::RCP<Map> SetupAmalgamationData(const RCP<Operator>& A, const RCP<GOVector>& globalgid2globalamalblockid_vector, LocalOrdinal blockSize) const;

    //! A Factory
    RCP<const FactoryBase> AFact_;

    //! nullspace factory
    //! The nullspace dimension is necessary for setting the block size in amalgamation routine
    RCP<const FactoryBase> nullspaceFact_;

    /// are we doing fixed or variable blocks
    mutable bool fixedBlkSize_;

    /// blocksize vector for variable blocksize setup
    mutable RCP<GOVector> blkSizeInfo_; // lives on overlapping column map of A

    /// pre-drop function
    RCP<PreDropFunctionBaseClass> predrop_;

    /// amalgamation information
    mutable RCP<std::map<GlobalOrdinal,std::vector<LocalOrdinal > > > globalamalblockid2myrowid_;
    mutable RCP<std::map<GlobalOrdinal,std::vector<GlobalOrdinal> > > globalamalblockid2globalrowid_;


  }; //class CoalesceDropFactory

} //namespace MueLu

#define MUELU_COALESCEDROPFACTORY_SHORT
#endif // MUELU_COALESCEDROPFACTORY_DECL_HPP
