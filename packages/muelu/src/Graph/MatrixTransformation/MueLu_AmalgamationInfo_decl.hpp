// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/*
 * MueLu_AmalgamationInfo_decl.hpp
 *
 *  Created on: Mar 28, 2012
 *      Author: wiesner
 */

#ifndef MUELU_AMALGAMATIONINFO_DECL_HPP_
#define MUELU_AMALGAMATIONINFO_DECL_HPP_

#include <Xpetra_ConfigDefs.hpp>  // global_size_t
#include <Xpetra_Map_fwd.hpp>
#include <Xpetra_Vector_fwd.hpp>
#include <Xpetra_MapFactory_fwd.hpp>

#include "MueLu_ConfigDefs.hpp"

#include "MueLu_BaseClass.hpp"

#include "MueLu_AmalgamationInfo_fwd.hpp"
#include "MueLu_Aggregates_fwd.hpp"

namespace MueLu {

/*!
  @class AmalgamationInfo
  @brief minimal container class for storing amalgamation information

  Helps create a mapping between local node id on current processor to local DOFs ids on
  current processor.  That mapping is used for unamalgamation.
*/

template <class LocalOrdinal  = DefaultLocalOrdinal,
          class GlobalOrdinal = DefaultGlobalOrdinal,
          class Node          = DefaultNode>
class AmalgamationInfo
  : public BaseClass {
#undef MUELU_AMALGAMATIONINFO_SHORT
#include "MueLu_UseShortNamesOrdinal.hpp"

 public:
  /// Constructor
  AmalgamationInfo(RCP<Array<LO> > rowTranslation,
                   RCP<Array<LO> > colTranslation,
                   RCP<const Map> nodeRowMap,
                   RCP<const Map> nodeColMap,
                   RCP<const Map> const& columnMap,
                   LO fullblocksize, GO offset, LO blockid, LO nStridedOffset, LO stridedblocksize)
    : rowTranslation_(rowTranslation)
    , colTranslation_(colTranslation)
    , nodeRowMap_(nodeRowMap)
    , nodeColMap_(nodeColMap)
    , columnMap_(columnMap)
    , fullblocksize_(fullblocksize)
    , offset_(offset)
    , blockid_(blockid)
    , nStridedOffset_(nStridedOffset)
    , stridedblocksize_(stridedblocksize)
    , indexBase_(columnMap->getIndexBase()) {}

  /// Destructor
  virtual ~AmalgamationInfo() {}

  /// Return a simple one-line description of this object.
  std::string description() const { return "AmalgamationInfo"; }

  //! Print the object with some verbosity level to an FancyOStream object.
  // using MueLu::Describable::describe; // overloading, not hiding
  // void describe(Teuchos::FancyOStream &out, const VerbLevel verbLevel = Default) const;;
  void print(Teuchos::FancyOStream& out, const VerbLevel verbLevel = Default) const;

  RCP<const Map> getNodeRowMap() const { return nodeRowMap_; }  //! < returns the node row map for the graph
  RCP<const Map> getNodeColMap() const { return nodeColMap_; }  //! < returns the node column map for the graph

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
  void UnamalgamateAggregates(const Aggregates& aggregates, Teuchos::ArrayRCP<LocalOrdinal>& aggStart, Teuchos::ArrayRCP<GlobalOrdinal>& aggToRowMap) const;
  void UnamalgamateAggregatesLO(const Aggregates& aggregates, Teuchos::ArrayRCP<LocalOrdinal>& aggStart, Teuchos::ArrayRCP<LO>& aggToRowMap) const;

  /*! @brief ComputeUnamalgamatedImportDofMap
   * build overlapping dof row map from aggregates needed for overlapping null space
   */
  Teuchos::RCP<Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> > ComputeUnamalgamatedImportDofMap(const Aggregates& aggregates) const;

 private:
  void UnamalgamateAggregates(const Teuchos::RCP<const Map>& nodeMap,
                              const RCP<LOVector>& procWinnerVec,
                              const RCP<LOMultiVector>& vertex2AggIdVec,
                              const GO numAggregates,
                              Teuchos::ArrayRCP<LocalOrdinal>& aggStart,
                              Teuchos::ArrayRCP<GlobalOrdinal>& aggToRowMap) const;

  void UnamalgamateAggregatesLO(const Teuchos::RCP<const Map>& nodeMap,
                                const RCP<LOVector>& procWinnerVec,
                                const RCP<LOMultiVector>& vertex2AggIdVec,
                                const GO numAggregates,
                                Teuchos::ArrayRCP<LocalOrdinal>& aggStart,
                                Teuchos::ArrayRCP<LO>& aggToRowMap) const;

  Teuchos::RCP<Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> > ComputeUnamalgamatedImportDofMap(const Teuchos::RCP<const Map>& nodeMap) const;

 public:
  /*! @brief ComputeGlobalDOF
   *
   * Return global dof id associated with global node id gNodeID and dof index k
   *
   * \note We assume that \c indexBase_ is valid for both the node and the dof map.
   *
   * @param (GO): global node id
   * @param (LO): local dof index within node
   * @return (GO): global dof id
   */
  GO ComputeGlobalDOF(GO const& gNodeID, LO const& k = 0) const;

  /*! @brief ComputeLocalDOF
   * return locbal dof id associated with local node id lNodeID and dof index k
   *
   * @param (LO): local node id
   * @param (LO): local dof index within node
   * @return (LO): local dof id
   */
  LO ComputeLocalDOF(LocalOrdinal const& lNodeID, LocalOrdinal const& k) const;

  LO ComputeLocalNode(LocalOrdinal const& ldofID) const;

  /*! Access routines */

  /// returns offset of global dof ids
  GO GlobalOffset() { return offset_; }

  /// returns striding information
  void GetStridingInformation(LO& fullBlockSize, LO& blockID, LO& stridingOffset, LO& stridedBlockSize, GO& indexBase) {
    fullBlockSize    = fullblocksize_;
    blockID          = blockid_;
    stridingOffset   = nStridedOffset_;
    stridedBlockSize = stridedblocksize_;
    indexBase        = indexBase_;
  }

 private:
  //! @name amalgamation information variables
  //@{

  //! Arrays containing local node ids given local dof ids
  RCP<Array<LO> > rowTranslation_;
  RCP<Array<LO> > colTranslation_;

  //! node row and column map of graph (built from row and column map of A)
  RCP<const Map> nodeRowMap_;
  RCP<const Map> nodeColMap_;

  /*! @brief DOF map (really column map of A)

  We keep a RCP on the column map to make sure that the map is still valid when it is used.
  */
  RCP<const Map> columnMap_;

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

}  // namespace MueLu

#define MUELU_AMALGAMATIONINFO_SHORT
#endif /* MUELU_AMALGAMATIONINFO_DECL_HPP_ */
