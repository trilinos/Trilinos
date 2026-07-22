// @HEADER
// *****************************************************************************
//             Xpetra: A linear algebra interface package
//
// Copyright 2012 NTESS and the Xpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PACKAGES_XPETRA_SUP_MAP_UTILS_HPP_
#define PACKAGES_XPETRA_SUP_MAP_UTILS_HPP_

#include "Xpetra_ConfigDefs.hpp"

#include "Xpetra_Exceptions.hpp"
#include "Xpetra_Map.hpp"
#include "Xpetra_MapFactory.hpp"

namespace Xpetra {

#ifndef DOXYGEN_SHOULD_SKIP_THIS
// forward declaration of BlockedMap, needed to prevent circular inclusions
template <class LO, class GO, class N>
class BlockedMap;
#endif

/*!
  @class MapUtils
  @brief Xpetra utility class for common map-related routines

  The routines should be independent from Epetra/Tpetra and be purely implemented in Xpetra.

*/
template <class LocalOrdinal,
          class GlobalOrdinal,
          class Node = Tpetra::KokkosClassic::DefaultNode::DefaultNodeType>
class MapUtils {
#undef XPETRA_MAPUTILS_SHORT
#include "Xpetra_UseShortNamesOrdinal.hpp"

 public:
  /*! @brief Helper function to concatenate several maps

    @param  subMaps    vector of maps which are concatenated
    @return            concatenated map

    The routine builds a global map by concatenating all provided maps in the ordering defined by the vector.
    The GIDs are just appended in the same ordering as in the subMaps. No reordering or sorting is performed.
    This routine is supposed to generate the full map in an Xpetra::MapExtractor for a block operator. Note, it
    should not be used for strided maps since the GIDs are not reordered.

    Example: subMap[0] = { 0, 1, 3, 4 };
             subMap[1] = { 2, 5 };
             concatenated map = { 0, 1, 3, 4, 2 ,5 };
    */
  static Teuchos::RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> > concatenateMaps(const std::vector<Teuchos::RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> > >& subMaps) {
    // ToDo Resolve header issues to allow for using this routing in Xpetra::BlockedMap.

    // merge submaps to global map
    std::vector<GlobalOrdinal> gids;
    for (size_t tt = 0; tt < subMaps.size(); ++tt) {
      Teuchos::RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> > subMap = subMaps[tt];
      Teuchos::ArrayView<const GlobalOrdinal> subMapGids                         = subMap->getLocalElementList();
      gids.insert(gids.end(), subMapGids.begin(), subMapGids.end());
    }

    const GlobalOrdinal INVALID = Teuchos::OrdinalTraits<Xpetra::global_size_t>::invalid();
    // std::sort(gids.begin(), gids.end());
    // gids.erase(std::unique(gids.begin(), gids.end()), gids.end());
    Teuchos::ArrayView<GlobalOrdinal> gidsView(&gids[0], gids.size());
    Teuchos::RCP<Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> > fullMap = Xpetra::MapFactory<LocalOrdinal, GlobalOrdinal, Node>::Build(subMaps[0]->lib(), INVALID, gidsView, subMaps[0]->getIndexBase(), subMaps[0]->getComm());
    return fullMap;
  }

  /*! @brief Helper function to shrink the GIDs and generate a standard map whith GIDs starting at 0
   *
    @param  input       Input map (may be overlapping) containing all GIDs. Think of it as a column map.
    @param  nonOvlInput Non-overlapping version of "input" map. Think of it is the corresponding domain map associated with the column map "input"
    @return             New map with unique continuous global ids starting with GID 0

    This helper routine may be useful for the transformation of MapExtractors in Xpetra-style GID ordering to the Thyra-style ordering.

    Example: input = { 10, 15, 26, 37, 48 }; on proc 0
             input = { 37, 48, 59, 60, 70 }; on proc 1
             nonOvlInput = { 10, 15, 26, 37 }; on proc 0
             nonOvlInput = { 48, 59, 60, 70 }: on proc 1
             result = { 0, 1, 2, 3, 4 }; on proc 0
             result = { 3, 4, 5, 6, 7 }; on proc 1
    */
  static Teuchos::RCP<Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> >
  shrinkMapGIDs(const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node>& input,
                const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node>& nonOvlInput) {
    TEUCHOS_TEST_FOR_EXCEPTION(nonOvlInput.getLocalNumElements() > input.getLocalNumElements(),
                               Xpetra::Exceptions::Incompatible,
                               "Xpetra::MatrixUtils::shrinkMapGIDs: the non-overlapping map must not have more local ids than the overlapping map.")

    TEUCHOS_TEST_FOR_EXCEPTION(nonOvlInput.getMaxAllGlobalIndex() != input.getMaxAllGlobalIndex(),
                               Xpetra::Exceptions::Incompatible,
                               "Xpetra::MatrixUtils::shrinkMapGIDs: the maximum GIDs of the overlapping and non-overlapping maps must be the same.")

    RCP<const Teuchos::Comm<int> > comm = input.getComm();

    // we expect input to be the potentially overlapping map associated with nonOvlInput as the non-overlapping
    // map with the same GIDs over all processors (e.g. column map and domain map). We use the nonOvlInput map
    // to determine which GIDs are owned by which processor.

    // calculate offset for new global Ids
    std::vector<int> myGIDs(comm->getSize(), 0);
    std::vector<int> numGIDs(comm->getSize(), 0);
    myGIDs[comm->getRank()] = nonOvlInput.getLocalNumElements();
    Teuchos::reduceAll(*comm, Teuchos::REDUCE_MAX, comm->getSize(), &myGIDs[0], &numGIDs[0]);
    size_t gidOffset = 0;
    for (int p = 0; p < comm->getRank(); p++) gidOffset += numGIDs[p];

    // we use nonOvlInput to assign the globally unique shrinked GIDs and communicate them to input.
    std::map<const GlobalOrdinal, GlobalOrdinal> origGID2newGID;
    for (size_t i = 0; i < nonOvlInput.getLocalNumElements(); i++) {
      origGID2newGID[nonOvlInput.getGlobalElement(i)] = Teuchos::as<GlobalOrdinal>(i) + Teuchos::as<GlobalOrdinal>(gidOffset);
    }
    // build an overlapping version of mySpecialMap
    Teuchos::Array<GlobalOrdinal> ovlUnknownStatusGids;
    Teuchos::Array<GlobalOrdinal> ovlFoundStatusGids;
    // loop over global column map of A and find all GIDs where it is not sure, whether they are special or not
    for (size_t i = 0; i < input.getLocalNumElements(); i++) {
      GlobalOrdinal gcid = input.getGlobalElement(i);
      if (nonOvlInput.isNodeGlobalElement(gcid) == false) {
        ovlUnknownStatusGids.push_back(gcid);
      }
    }

    // Communicate the number of DOFs on each processor
    std::vector<int> myUnknownDofGIDs(comm->getSize(), 0);
    std::vector<int> numUnknownDofGIDs(comm->getSize(), 0);
    myUnknownDofGIDs[comm->getRank()] = ovlUnknownStatusGids.size();
    Teuchos::reduceAll(*comm, Teuchos::REDUCE_MAX, comm->getSize(), &myUnknownDofGIDs[0], &numUnknownDofGIDs[0]);

    // create array containing all DOF GIDs
    size_t cntUnknownDofGIDs = 0;
    for (int p = 0; p < comm->getSize(); p++) cntUnknownDofGIDs += numUnknownDofGIDs[p];
    std::vector<GlobalOrdinal> lUnknownDofGIDs(cntUnknownDofGIDs, 0);  // local version to be filled
    std::vector<GlobalOrdinal> gUnknownDofGIDs(cntUnknownDofGIDs, 0);  // global version after communication
    // calculate the offset and fill chunk of memory with local data on each processor
    size_t cntUnknownOffset = 0;
    for (int p = 0; p < comm->getRank(); p++) cntUnknownOffset += numUnknownDofGIDs[p];
    for (size_t k = 0; k < Teuchos::as<size_t>(ovlUnknownStatusGids.size()); k++) {
      lUnknownDofGIDs[k + cntUnknownOffset] = ovlUnknownStatusGids[k];
    }
    if (cntUnknownDofGIDs > 0)  // only perform communication if there are unknown DOF GIDs
      Teuchos::reduceAll(*comm, Teuchos::REDUCE_MAX, Teuchos::as<int>(cntUnknownDofGIDs), &lUnknownDofGIDs[0], &gUnknownDofGIDs[0]);
    std::vector<GlobalOrdinal> lTranslatedDofGIDs(cntUnknownDofGIDs, 0);  // local version to be filled
    std::vector<GlobalOrdinal> gTranslatedDofGIDs(cntUnknownDofGIDs, 0);  // global version after communication
    // loop through all GIDs with unknown status
    for (size_t k = 0; k < gUnknownDofGIDs.size(); k++) {
      GlobalOrdinal curgid = gUnknownDofGIDs[k];
      if (nonOvlInput.isNodeGlobalElement(curgid)) {
        lTranslatedDofGIDs[k] = origGID2newGID[curgid];  // curgid is in special map (on this processor)
      }
    }
    if (cntUnknownDofGIDs > 0)  // only perform communication if there are unknown DOF GIDs
      Teuchos::reduceAll(*comm, Teuchos::REDUCE_MAX, Teuchos::as<int>(cntUnknownDofGIDs), &lTranslatedDofGIDs[0], &gTranslatedDofGIDs[0]);

    for (size_t k = 0; k < Teuchos::as<size_t>(ovlUnknownStatusGids.size()); k++) {
      origGID2newGID[ovlUnknownStatusGids[k]] = gTranslatedDofGIDs[k + cntUnknownOffset];
    }
    Teuchos::Array<GlobalOrdinal> ovlDomainMapArray;
    for (size_t i = 0; i < input.getLocalNumElements(); i++) {
      GlobalOrdinal gcid = input.getGlobalElement(i);
      ovlDomainMapArray.push_back(origGID2newGID[gcid]);
    }
    RCP<Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> > ovlDomainMap =
        Xpetra::MapFactory<LocalOrdinal, GlobalOrdinal, Node>::Build(nonOvlInput.lib(), Teuchos::OrdinalTraits<GlobalOrdinal>::invalid(), ovlDomainMapArray(), 0, comm);
    return ovlDomainMap;
  }

  /*! @brief replace set of global ids by new global ids
   *
    @param  input                Overlapping input map.
    @param  nonOvlInput          Non-overlapping map containing GIDs corresponding to "input". Think of it is the corresponding domain map associated with the column map "input"
    @param  nonOvlReferenceInput Non-overlapping reference map containing new GIDs.
    @return                      Overlapping map compatible to "input" using the GIDs as defined by "nonOvlReferenceInput"

    Example: input = { 0, 1, 2, 3 }; on proc 0
             input = { 2, 3, 4, 5 }; on proc 1
             nonOvlInput = { 0, 1, 2 }; on proc 0
             nonOvlInput = { 3, 4, 5 }: on proc 1
             nonOvlReferenceInput = { 33, 44, 55 }; on proc 0
             nonOvlReferenceInput = { 101, 102, 103 }; on proc 1
             result = { 33, 44, 55, 101 }; on proc 0
             result = { 55, 101, 102, 103}; on proc 1
    */
  static Teuchos::RCP<Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> > transformThyra2XpetraGIDs(
      const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node>& input,
      const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node>& nonOvlInput,
      const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node>& nonOvlReferenceInput) {
    // TEUCHOS_TEST_FOR_EXCEPTION(nonOvlInput.getLocalNumElements() > input.getLocalNumElements(), Xpetra::Exceptions::Incompatible, "Xpetra::MatrixUtils::transformThyra2XpetraGIDs: the non-overlapping map must not have more local ids than the overlapping map.");
    TEUCHOS_TEST_FOR_EXCEPTION(nonOvlInput.getLocalNumElements() != nonOvlReferenceInput.getLocalNumElements(), Xpetra::Exceptions::Incompatible, "Xpetra::MatrixUtils::transformThyra2XpetraGIDs: the number of local Xpetra reference GIDs and local Thyra GIDs of the non-overlapping maps must be the same!");
    // TEUCHOS_TEST_FOR_EXCEPTION(nonOvlInput.getMaxAllGlobalIndex() != input.getMaxAllGlobalIndex(), Xpetra::Exceptions::Incompatible, "Xpetra::MatrixUtils::transformThyra2XpetraGIDs: the maximum GIDs of the overlapping and non-overlapping maps must be the same. nonOvlInput.getMaxAllGlobalIndex() = " << nonOvlInput.getMaxAllGlobalIndex() << " ovlInput.getMaxAllGlobalIndex() = " << input.getMaxAllGlobalIndex());

    RCP<const Teuchos::Comm<int> > comm = input.getComm();

    // fill translation map as far as possible
    std::map<const GlobalOrdinal, GlobalOrdinal> thyra2xpetraGID;
    for (size_t i = 0; i < nonOvlInput.getLocalNumElements(); i++) {
      thyra2xpetraGID[nonOvlInput.getGlobalElement(i)] =
          nonOvlReferenceInput.getGlobalElement(i);
    }

    // find all GIDs of the overlapping Thyra map which are not owned by this proc
    Teuchos::Array<GlobalOrdinal> ovlUnknownStatusGids;
    // loop over global column map of A and find all GIDs where it is not sure, whether they are special or not
    for (size_t i = 0; i < input.getLocalNumElements(); i++) {
      GlobalOrdinal gcid = input.getGlobalElement(i);
      if (nonOvlInput.isNodeGlobalElement(gcid) == false) {
        ovlUnknownStatusGids.push_back(gcid);
      }
    }

    // Communicate the number of DOFs on each processor
    std::vector<int> myUnknownDofGIDs(comm->getSize(), 0);
    std::vector<int> numUnknownDofGIDs(comm->getSize(), 0);
    myUnknownDofGIDs[comm->getRank()] = ovlUnknownStatusGids.size();
    Teuchos::reduceAll(*comm, Teuchos::REDUCE_MAX, comm->getSize(), &myUnknownDofGIDs[0], &numUnknownDofGIDs[0]);

    // create array containing all DOF GIDs
    size_t cntUnknownDofGIDs = 0;
    for (int p = 0; p < comm->getSize(); p++) cntUnknownDofGIDs += numUnknownDofGIDs[p];
    std::vector<GlobalOrdinal> lUnknownDofGIDs(cntUnknownDofGIDs, 0);  // local version to be filled
    std::vector<GlobalOrdinal> gUnknownDofGIDs(cntUnknownDofGIDs, 0);  // global version after communication
    // calculate the offset and fill chunk of memory with local data on each processor
    size_t cntUnknownOffset = 0;
    for (int p = 0; p < comm->getRank(); p++) cntUnknownOffset += numUnknownDofGIDs[p];
    for (size_t k = 0; k < Teuchos::as<size_t>(ovlUnknownStatusGids.size()); k++) {
      lUnknownDofGIDs[k + cntUnknownOffset] = ovlUnknownStatusGids[k];
    }
    if (cntUnknownDofGIDs > 0)  // only perform communication if there are unknown DOF GIDs
      Teuchos::reduceAll(*comm, Teuchos::REDUCE_MAX, Teuchos::as<int>(cntUnknownDofGIDs), &lUnknownDofGIDs[0], &gUnknownDofGIDs[0]);
    std::vector<GlobalOrdinal> lTranslatedDofGIDs(cntUnknownDofGIDs, 0);  // local version to be filled
    std::vector<GlobalOrdinal> gTranslatedDofGIDs(cntUnknownDofGIDs, 0);  // global version after communication
    // loop through all GIDs with unknown status
    for (size_t k = 0; k < gUnknownDofGIDs.size(); k++) {
      GlobalOrdinal curgid = gUnknownDofGIDs[k];
      if (nonOvlInput.isNodeGlobalElement(curgid)) {
        lTranslatedDofGIDs[k] = thyra2xpetraGID[curgid];
      }
    }
    if (cntUnknownDofGIDs > 0)  // only perform communication if there are unknown DOF GIDs
      Teuchos::reduceAll(*comm, Teuchos::REDUCE_MAX, Teuchos::as<int>(cntUnknownDofGIDs), &lTranslatedDofGIDs[0], &gTranslatedDofGIDs[0]);

    for (size_t k = 0; k < Teuchos::as<size_t>(ovlUnknownStatusGids.size()); k++) {
      thyra2xpetraGID[ovlUnknownStatusGids[k]] = gTranslatedDofGIDs[k + cntUnknownOffset];
    }
    Teuchos::Array<GlobalOrdinal> ovlDomainMapArray;
    for (size_t i = 0; i < input.getLocalNumElements(); i++) {
      GlobalOrdinal gcid = input.getGlobalElement(i);
      ovlDomainMapArray.push_back(thyra2xpetraGID[gcid]);
    }
    RCP<Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> > ovlDomainMap =
        Xpetra::MapFactory<LocalOrdinal, GlobalOrdinal, Node>::Build(nonOvlInput.lib(), Teuchos::OrdinalTraits<GlobalOrdinal>::invalid(), ovlDomainMapArray(), 0, comm);

    TEUCHOS_TEST_FOR_EXCEPTION(input.getLocalNumElements() != ovlDomainMap->getLocalNumElements(), Xpetra::Exceptions::Incompatible, "Xpetra::MatrixUtils::transformThyra2XpetraGIDs: the number of local Thyra reference GIDs (overlapping) and local Xpetra GIDs (overlapping) must be the same!");
    // TEUCHOS_TEST_FOR_EXCEPTION(nonOvlReferenceInput.getMaxAllGlobalIndex() != ovlDomainMap->getMaxAllGlobalIndex(), Xpetra::Exceptions::Incompatible, "Xpetra::MatrixUtils::transformThyra2XpetraGIDs: the maximum GIDs of the overlapping and non-overlapping Xpetra maps must be the same.");

    return ovlDomainMap;
  }

  /*! @brief replace set of global ids by new global ids
   *
    @param  input map (either Map or BlockedMap) containing Thyra GIDs
    @param  offset GID offset for resulting Xpetra GIDs
    @return Map (or BlockedMap) containing Xpetra GIDs
  */
  static Teuchos::RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> > transformThyra2XpetraGIDs(
      const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node>& input, GlobalOrdinal offset) {
    const GO INVALID                    = Teuchos::OrdinalTraits<Xpetra::global_size_t>::invalid();
    RCP<const Teuchos::Comm<int> > comm = input.getComm();

    RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> > rcpInput = Teuchos::rcpFromRef(input);

    // check whether input map is a BlockedMap or a standard Map
    RCP<const Xpetra::BlockedMap<LocalOrdinal, GlobalOrdinal, Node> > rcpBlockedInput = Teuchos::rcp_dynamic_cast<const Xpetra::BlockedMap<LocalOrdinal, GlobalOrdinal, Node> >(rcpInput);
    if (rcpBlockedInput.is_null() == true) {
      // create a new map with Xpetra GIDs (may start not from GID = 0)
      std::vector<GlobalOrdinal> gids;
      for (LocalOrdinal l = 0; l < Teuchos::as<LocalOrdinal>(rcpInput->getLocalNumElements()); ++l) {
        GlobalOrdinal gid = rcpInput->getGlobalElement(l) + offset;
        gids.push_back(gid);
      }
      Teuchos::ArrayView<GO> gidsView(&gids[0], gids.size());
      RCP<Map> fullMap = MapFactory::Build(rcpInput->lib(), INVALID, gidsView, rcpInput->getIndexBase(), comm);
      return fullMap;
    }

    // SPECIAL CASE: input is a blocked map
    // we have to recursively call this routine to get a BlockedMap containing (unique) Xpetra style GIDs

    size_t numMaps = rcpBlockedInput->getNumMaps();

    // first calucale GID offsets in submaps
    // we need that for generating Xpetra GIDs
    std::vector<GlobalOrdinal> gidOffsets(numMaps, 0);
    for (size_t i = 1; i < numMaps; ++i) {
      gidOffsets[i] = rcpBlockedInput->getMap(i - 1, true)->getMaxAllGlobalIndex() + gidOffsets[i - 1] + 1;
    }

    std::vector<RCP<const Map> > mapsXpetra(rcpBlockedInput->getNumMaps(), Teuchos::null);
    std::vector<RCP<const Map> > mapsThyra(rcpBlockedInput->getNumMaps(), Teuchos::null);
    for (size_t b = 0; b < rcpBlockedInput->getNumMaps(); ++b) {
      // extract sub map with Thyra style gids
      // this can be an underlying Map or BlockedMap object
      RCP<const Map> subMapThyra  = rcpBlockedInput->getMap(b, true);
      RCP<const Map> subMapXpetra = MapUtils::transformThyra2XpetraGIDs(*subMapThyra, gidOffsets[b] + offset);  // recursive call
      mapsXpetra[b]               = subMapXpetra;                                                               // map can be of type Map or BlockedMap
      mapsThyra[b]                = subMapThyra;                                                                // map can be of type Map or BlockedMap
    }

    Teuchos::RCP<Map> resultMap = Teuchos::rcp(new Xpetra::BlockedMap<LocalOrdinal, GlobalOrdinal, Node>(mapsXpetra, mapsThyra));
    return resultMap;
  }
};

}  // end namespace Xpetra

#define XPETRA_MAPUTILS_SHORT

#endif  // PACKAGES_XPETRA_SUP_MAP_UTILS_HPP_
