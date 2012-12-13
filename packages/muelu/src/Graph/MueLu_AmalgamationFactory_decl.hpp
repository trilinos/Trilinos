// @HEADER
//
// ***********************************************************************
//
//        MueLu: A package for multigrid based preconditioning
//                  Copyright 2012 Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact
//                    Jeremie Gaidamour (jngaida@sandia.gov)
//                    Jonathan Hu       (jhu@sandia.gov)
//                    Ray Tuminaro      (rstumin@sandia.gov)
//
// ***********************************************************************
//
// @HEADER
#ifndef MUELU_AMALGAMATIONFACTORY_DECL_HPP
#define MUELU_AMALGAMATIONFACTORY_DECL_HPP

#include <Xpetra_Matrix_fwd.hpp>
#include <Xpetra_Map_fwd.hpp>

#include "MueLu_ConfigDefs.hpp"
#include "MueLu_SingleLevelFactoryBase.hpp"
#include "MueLu_Aggregates_fwd.hpp"
#include "MueLu_AmalgamationInfo_fwd.hpp"

#include "MueLu_Level_fwd.hpp"
#include "MueLu_Exceptions.hpp"

namespace MueLu {

  /*!
    @class AmalgamationFactory
    @brief AmalgamationFactory for subblocks of strided map based amalgamation data

    Class generates unamalgamation information using matrix A with strided maps.
    It stores the output information within an AmalgamationInfo object as "UnAmalgamationInfo".
    This object contains

    \li \c nodegid2dofgids_ a map of all node ids of which the current proc has corresponding DOF gids (used by \c TentativePFactory).
    \li \c gNodeIds vector of all node ids on the current proc (may be less than nodegid2dofgids_.size()). These nodes are stored on the current proc.

  */

  template <class Scalar = double, class LocalOrdinal = int, class GlobalOrdinal = LocalOrdinal, class Node = Kokkos::DefaultNode::DefaultNodeType, class LocalMatOps = typename Kokkos::DefaultKernels<void,LocalOrdinal,Node>::SparseOps>
  class AmalgamationFactory : public SingleLevelFactoryBase {
#undef MUELU_AMALGAMATIONFACTORY_SHORT
#include "MueLu_UseShortNames.hpp"

  public:

    //! @name Constructors/Destructors.
    //@{

    //! Constructor
    AmalgamationFactory() { }

    //! Destructor
    virtual ~AmalgamationFactory() { }

    RCP<const ParameterList> GetValidParameterList(const ParameterList& paramList = ParameterList()) const;

    //@}

    //! Input
    //@{

    void DeclareInput(Level &currentLevel) const;

    //@}

    void Build(Level &currentLevel) const; // Build

    //! translate global (row/column) id to global amalgamation block id
    // @param gid (GlobalOrdinal): input global id (row gid or column gid)
    // @param A: input operator (just used to check the maps for validity)
    // @param blockSize (LocalOrdinal): block size (needed for constant block size)
    // @param offset (GlobalOrdinal): global offset for dofs (stored in strided map, default = 0)
    static const GlobalOrdinal DOFGid2NodeId(GlobalOrdinal gid, const RCP<Matrix>& A, LocalOrdinal blockSize, const GlobalOrdinal offset = 0);

    /*! @brief ComputeUnamalgamatedAggregateSizes
     * computes the size of the aggregates (in DOFs)
     */
    static void ComputeUnamalgamatedAggregateSizes(const Aggregates& aggregates, const AmalgamationInfo& amalgInfo, Teuchos::ArrayRCP<LocalOrdinal> & aggSizes);

    /*! @brief UnamalgamateAggregates
     * puts all dofs for aggregate \c i in aggToRowMap[\c i]
     */
    static void UnamalgamateAggregates(const Aggregates& aggregates, const AmalgamationInfo& amalgInfo, const Teuchos::ArrayRCP<LocalOrdinal> & aggSizes, Teuchos::ArrayRCP<Teuchos::ArrayRCP<GlobalOrdinal> > & aggToRowMap);

    /*! @brief ComputeUnamalgamatedImportDofMap
     * build overlapping dof row map from aggregates needed for overlapping null space
     */
    static Teuchos::RCP< Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> > ComputeUnamalgamatedImportDofMap(const Aggregates& aggregates, const AmalgamationInfo& amalgInfo, const Teuchos::ArrayRCP<LocalOrdinal> & aggSizes);


  private:

    // amalgamation information
    mutable RCP<std::map<GlobalOrdinal,std::vector<GlobalOrdinal> > > nodegid2dofgids_;

  }; //class SubBlockUnAmalgamationFactory

} //namespace MueLu

#define MUELU_AMALGAMATIONFACTORY_SHORT
#endif // MUELU_AMALGAMATIONFACTORY_DECL_HPP
