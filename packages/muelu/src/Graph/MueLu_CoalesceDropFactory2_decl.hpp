/*
 * MueLu_CoalesceDropFactory2_decl.hpp
 *
 *  Created on: 10.07.2012
 *      Author: tobias
 */

#ifndef MUELU_COALESCEDROPFACTORY2_DECL_HPP_
#define MUELU_COALESCEDROPFACTORY2_DECL_HPP_

#include <Xpetra_Map_fwd.hpp>
#include <Xpetra_MapFactory_fwd.hpp>
#include <Xpetra_Vector_fwd.hpp>
#include <Xpetra_MultiVector_fwd.hpp>
#include <Xpetra_CrsGraph_fwd.hpp>
#include <Xpetra_CrsGraphFactory.hpp> //TODO
#include <Xpetra_Operator_fwd.hpp>
#include <Xpetra_StridedMap_fwd.hpp>

#include "MueLu_ConfigDefs.hpp"
#include "MueLu_SingleLevelFactoryBase.hpp"
#include "MueLu_CoalesceDropFactory2_fwd.hpp"

#include "MueLu_Level_fwd.hpp"
#include "MueLu_Graph_fwd.hpp"
#include "MueLu_AmalgamationInfo_fwd.hpp"
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
  class CoalesceDropFactory2 : public SingleLevelFactoryBase {
#undef MUELU_COALESCEDROPFACTORY2_SHORT
#include "MueLu_UseShortNames.hpp"

  public:

    //! @name Constructors/Destructors.
    //@{

    //! Constructor
    CoalesceDropFactory2(RCP<const FactoryBase> AFact = Teuchos::null);

    //! Destructor
    virtual ~CoalesceDropFactory2() { }

    //@}

    //! Input
    //@{

    void DeclareInput(Level &currentLevel) const;

    //@}

    void Build(Level &currentLevel) const; // Build

    //void Amalgamate(const RCP<Operator>& A, RCP<Graph>& graph) const; // Amalgamate

  private:

    //! translate global (row/column) id to global amalgamation block id
    // @param gid (GlobalOrdinal): input global id (row gid or column gid)
    // @param A: input operator (just used to check the maps for validity)
    // @param blockSize (LocalOrdinal): block size (needed for constant block size)
    // @param offset (GlobalOrdinal): global offset for dofs (stored in strided map, default = 0)
    GlobalOrdinal DOFGid2NodeId(GlobalOrdinal gid, const RCP<Operator>& A, LocalOrdinal blockSize, const GlobalOrdinal offset = 0) const;

    //! setup amalgamation data
    // This routine fills the private members lobalamalblockid2myrowid_ and lobalamalblockid2globalrowid_
    // These are then used in Graph and Aggregates for unamalgamation
    // @param A: input operator (just used to check the maps for validity)
    // returns amalgamated map
    //const Teuchos::RCP<Map> SetupAmalgamationData(const RCP<Operator>& A) const;

    //! A Factory
    RCP<const FactoryBase> AFact_;

    /// amalgamation information
    mutable RCP<std::map<GlobalOrdinal,std::vector<GlobalOrdinal> > > nodegid2dofgids_;

  }; //class CoalesceDropFactory2

} //namespace MueLu

#define MUELU_COALESCEDROPFACTORY2_SHORT


#endif /* MUELU_COALESCEDROPFACTORY2_DECL_HPP_ */
