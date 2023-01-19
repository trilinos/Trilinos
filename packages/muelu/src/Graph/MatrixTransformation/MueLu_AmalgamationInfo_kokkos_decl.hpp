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
//                    Jonathan Hu       (jhu@sandia.gov)
//                    Andrey Prokopenko (aprokop@sandia.gov)
//                    Ray Tuminaro      (rstumin@sandia.gov)
//
// ***********************************************************************
//
// @HEADER
/*
 * MueLu_AmalgamationInfo_decl.hpp
 *
 *  Created on: Mar 28, 2012
 *      Author: wiesner
 */

#ifndef MUELU_AMALGAMATIONINFO_KOKKOS_DECL_HPP_
#define MUELU_AMALGAMATIONINFO_KOKKOS_DECL_HPP_

#include <Xpetra_ConfigDefs.hpp>   // global_size_t
#include <Xpetra_Map_fwd.hpp>
#include <Xpetra_MapFactory_fwd.hpp>
#include <Xpetra_Vector_fwd.hpp>

#include "MueLu_ConfigDefs.hpp"

#include "MueLu_BaseClass.hpp"

#include "MueLu_AmalgamationInfo_kokkos_fwd.hpp"
#include "MueLu_Aggregates_kokkos_fwd.hpp"

namespace MueLu {

/*!
  @class AmalgamationInfo_kokkos
  @brief minimal container class for storing amalgamation information

  Helps create a mapping between global node id on current processor to global DOFs ids on
  current processor.  That mapping is used for unamalgamation.
*/

  template<class LocalOrdinal = DefaultLocalOrdinal,
           class GlobalOrdinal = DefaultGlobalOrdinal,
           class Node = DefaultNode>
  class AmalgamationInfo_kokkos
    : public BaseClass {
#undef MUELU_AMALGAMATIONINFO_KOKKOS_SHORT
#include "MueLu_UseShortNamesOrdinal.hpp"

  public:

    AmalgamationInfo_kokkos(RCP<Array<LO> > rowTranslation,
                            RCP<Array<LO> > colTranslation,
                            RCP<const Map> nodeRowMap,
                            RCP<const Map> nodeColMap,
                            RCP< const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> > const &columnMap,
                            LO fullblocksize, GO offset, LO blockid, LO nStridedOffset, LO stridedblocksize) :
      rowTranslation_(rowTranslation),
      colTranslation_(colTranslation),
      nodeRowMap_(nodeRowMap),
      nodeColMap_(nodeColMap),
      columnMap_(columnMap),
      fullblocksize_(fullblocksize),
      offset_(offset),
      blockid_(blockid),
      nStridedOffset_(nStridedOffset),
      stridedblocksize_(stridedblocksize),
      indexBase_(columnMap->getIndexBase())
    {}

    virtual ~AmalgamationInfo_kokkos() {}

    /// Return a simple one-line description of this object.
    std::string description() const { return "AmalgamationInfo"; }

    //! Print the object with some verbosity level to an FancyOStream object.
    //using MueLu::Describable::describe; // overloading, not hiding
    //void describe(Teuchos::FancyOStream &out, const VerbLevel verbLevel = Default) const;;
    void print(Teuchos::FancyOStream &out, const VerbLevel verbLevel = Default) const;

    RCP<const Map> getNodeRowMap() const { return nodeRowMap_; } //! < returns the node row map for the graph
    RCP<const Map> getNodeColMap() const { return nodeColMap_; } //! < returns the node column map for the graph

    /* @brief Translation arrays
     *
     * Returns translation arrays providing local node ids given local dof ids built from either
     * the non-overlapping (unique) row map or the overlapping (non-unique) column map.
     * The getColTranslation routine, e.g., is used for the MergeRows routine in CoalesceDropFactory.
     */
    //@{
    RCP<Array<LO> > getRowTranslation() const { return rowTranslation_; }
    RCP<Array<LO> > getColTranslation() const { return colTranslation_; }
    //@}

    /*! @brief UnamalgamateAggregates

       Puts all dofs for aggregate \c i in aggToRowMap[\c i].  Also calculate aggregate sizes.
    */
    void UnamalgamateAggregates(const Aggregates_kokkos& aggregates, Teuchos::ArrayRCP<LocalOrdinal>& aggStart, Teuchos::ArrayRCP<GlobalOrdinal>& aggToRowMap) const;
    void UnamalgamateAggregatesLO(const Aggregates_kokkos& aggregates, Teuchos::ArrayRCP<LocalOrdinal>& aggStart, Teuchos::ArrayRCP<LO>& aggToRowMap) const;

    /*! @brief ComputeUnamalgamatedImportDofMap
     * build overlapping dof row map from aggregates needed for overlapping null space
     */
    Teuchos::RCP< Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> > ComputeUnamalgamatedImportDofMap(const Aggregates_kokkos& aggregates) const;

    /*! @brief ComputeGlobalDOF
     * return global dof id associated with global node id gNodeID and dof index k
     *
     * @param (GO): global node id
     * @param (LO): local dof index within node
     * @return (GO): global dof id
     */
    GO ComputeGlobalDOF(GO const &gNodeID, LO const &k=0) const;

    /*! Access routines */

    /// returns offset of global dof ids
    GO GlobalOffset() { return offset_; }

    /// returns striding information
    void GetStridingInformation(LO& fullBlockSize, LO& blockID, LO& stridingOffset, LO& stridedBlockSize, GO& indexBase) {
      fullBlockSize = fullblocksize_;
      blockID = blockid_;
      stridingOffset = nStridedOffset_;
      stridedBlockSize = stridedblocksize_;
      indexBase = indexBase_;
    }

  private:

    //! @name amalgamation information variables
    //@{

    // arrays containing local node ids given local dof ids
    RCP<Array<LO> > rowTranslation_;
    RCP<Array<LO> > colTranslation_;

    // node row and column map of graph (built from row and column map of A)
    RCP<const Map> nodeRowMap_;
    RCP<const Map> nodeColMap_;

    //! @brief DOF map (really column map of A)
    // keep an RCP on the column map to make sure that the map is still valid when it is used
    Teuchos::RCP< const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> > columnMap_;

    //@}

    //! @name Strided map information.
    //@{
    LO fullblocksize_;
    GO offset_;
    LO blockid_;
    LO nStridedOffset_;
    LO stridedblocksize_;
    GO indexBase_;
    //@}

  };

} // namespace MueLu

#define MUELU_AMALGAMATIONINFO_KOKKOS_SHORT
#endif /* MUELU_AMALGAMATIONINFO_KOKKOS_DECL_HPP_ */
