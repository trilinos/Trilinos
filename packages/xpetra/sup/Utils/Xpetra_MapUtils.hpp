// @HEADER
//
// ***********************************************************************
//
//             Xpetra: A linear algebra interface package
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
//                    Tobias Wiesner    (tawiesn@sandia.gov)
//
// ***********************************************************************
//
// @HEADER
#ifndef PACKAGES_XPETRA_SUP_MAP_UTILS_HPP_
#define PACKAGES_XPETRA_SUP_MAP_UTILS_HPP_

#include "Xpetra_ConfigDefs.hpp"

#include "Xpetra_Map.hpp"

namespace Xpetra {

/*!
  @class MapUtils
  @brief Xpetra utility class for common map-related routines

  The routines should be independent from Epetra/Tpetra and be purely implemented in Xpetra.

*/
template <class LocalOrdinal,
          class GlobalOrdinal,
          class Node>
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
  static Teuchos::RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> > concatenateMaps(std::vector<Teuchos::RCP<const Xpetra::Map<LocalOrdinal,GlobalOrdinal,Node> > > & subMaps) {

    // merge submaps to global map
    std::vector<GlobalOrdinal> gids;
    for(size_t tt = 0; tt<subMaps.size(); ++tt) {
      Teuchos::RCP<const Xpetra::Map<LocalOrdinal,GlobalOrdinal,Node> > subMap = subMaps[tt];
      for(LocalOrdinal l = 0; l < Teuchos::as<LocalOrdinal>(subMap->getNodeNumElements()); ++l) {
        GlobalOrdinal gid = subMap->getGlobalElement(l);
        gids.push_back(gid);
      }
    }

    const GlobalOrdinal INVALID = Teuchos::OrdinalTraits<Xpetra::global_size_t>::invalid();
    //std::sort(gids.begin(), gids.end());
    //gids.erase(std::unique(gids.begin(), gids.end()), gids.end());
    Teuchos::ArrayView<GlobalOrdinal> gidsView(&gids[0], gids.size());
    Teuchos::RCP<Xpetra::Map<LocalOrdinal,GlobalOrdinal,Node> > fullMap = Xpetra::MapFactory<LocalOrdinal,GlobalOrdinal,Node>::Build(subMaps[0]->lib(), INVALID, gidsView, subMaps[0]->getIndexBase(), subMaps[0]->getComm());
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
  static Teuchos::RCP<Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> > shrinkMapGIDs(
      const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node>& input,
      const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node>& nonOvlInput ) {
    TEUCHOS_TEST_FOR_EXCEPTION(nonOvlInput.getNodeNumElements() > input.getNodeNumElements(), Xpetra::Exceptions::Incompatible, "Xpetra::MatrixUtils::shrinkMapGIDs: the non-overlapping map must not have more local ids than the overlapping map.")
    TEUCHOS_TEST_FOR_EXCEPTION(nonOvlInput.getMaxAllGlobalIndex() != input.getMaxAllGlobalIndex(), Xpetra::Exceptions::Incompatible, "Xpetra::MatrixUtils::shrinkMapGIDs: the maximum GIDs of the overlapping and non-overlapping maps must be the same.")

    RCP< const Teuchos::Comm<int> > comm = input.getComm();

    // we expect input to be the potentially overlapping map associated with nonOvlInput as the non-overlapping
    // map with the same GIDs over all processors (e.g. column map and domain map). We use the nonOvlInput map
    // to determine which GIDs are owned by which processor.

    // calculate offset for new global Ids
    std::vector<int> myGIDs(comm->getSize(),0);
    std::vector<int> numGIDs(comm->getSize(),0);
    myGIDs[comm->getRank()] = nonOvlInput.getNodeNumElements();
    Teuchos::reduceAll(*comm,Teuchos::REDUCE_MAX,comm->getSize(),&myGIDs[0],&numGIDs[0]);
    size_t gidOffset = 0;
    for(int p = 0; p < comm->getRank(); p++) gidOffset += numGIDs[p];

    // we use nonOvlInput to assing the globally unique shrinked GIDs and communicate them to input.
    std::map<const GlobalOrdinal, GlobalOrdinal> origGID2newGID;
    for(size_t i = 0; i < nonOvlInput.getNodeNumElements(); i++) {
      origGID2newGID[nonOvlInput.getGlobalElement(i)] = Teuchos::as<GlobalOrdinal>(i) + Teuchos::as<GlobalOrdinal>(gidOffset);
    }
    // build an overlapping version of mySpecialMap
    Teuchos::Array<GlobalOrdinal> ovlUnknownStatusGids;
    Teuchos::Array<GlobalOrdinal> ovlFoundStatusGids;
    // loop over global column map of A and find all GIDs where it is not sure, whether they are special or not
    for(size_t i = 0; i<input.getNodeNumElements(); i++) {
      GlobalOrdinal gcid = input.getGlobalElement(i);
      if( nonOvlInput.isNodeGlobalElement(gcid) == false) {
        ovlUnknownStatusGids.push_back(gcid);
      }
    }

    // Communicate the number of DOFs on each processor
    std::vector<int> myUnknownDofGIDs(comm->getSize(),0);
    std::vector<int> numUnknownDofGIDs(comm->getSize(),0);
    myUnknownDofGIDs[comm->getRank()] = ovlUnknownStatusGids.size();
    Teuchos::reduceAll(*comm,Teuchos::REDUCE_MAX,comm->getSize(),&myUnknownDofGIDs[0],&numUnknownDofGIDs[0]);

    // create array containing all DOF GIDs
    size_t cntUnknownDofGIDs = 0;
    for(int p = 0; p < comm->getSize(); p++) cntUnknownDofGIDs += numUnknownDofGIDs[p];
    std::vector<GlobalOrdinal> lUnknownDofGIDs(cntUnknownDofGIDs,0); // local version to be filled
    std::vector<GlobalOrdinal> gUnknownDofGIDs(cntUnknownDofGIDs,0); // global version after communication
    // calculate the offset and fill chunk of memory with local data on each processor
    size_t cntUnknownOffset = 0;
    for(int p = 0; p < comm->getRank(); p++) cntUnknownOffset += numUnknownDofGIDs[p];
    for(size_t k=0; k < Teuchos::as<size_t>(ovlUnknownStatusGids.size()); k++) {
      lUnknownDofGIDs[k+cntUnknownOffset] = ovlUnknownStatusGids[k];
    }
    if(cntUnknownDofGIDs > 0) // only perform communication if there are unknown DOF GIDs
      Teuchos::reduceAll(*comm,Teuchos::REDUCE_MAX,Teuchos::as<int>(cntUnknownDofGIDs),&lUnknownDofGIDs[0],&gUnknownDofGIDs[0]);
    std::vector<GlobalOrdinal> lTranslatedDofGIDs(cntUnknownDofGIDs,0); // local version to be filled
    std::vector<GlobalOrdinal> gTranslatedDofGIDs(cntUnknownDofGIDs,0); // global version after communication
    // loop through all GIDs with unknown status
    for(size_t k=0; k < gUnknownDofGIDs.size(); k++) {
      GlobalOrdinal curgid = gUnknownDofGIDs[k];
      if(nonOvlInput.isNodeGlobalElement(curgid)) {
        lTranslatedDofGIDs[k] = origGID2newGID[curgid]; // curgid is in special map (on this processor)
      }
    }
    if(cntUnknownDofGIDs > 0) // only perform communication if there are unknown DOF GIDs
      Teuchos::reduceAll(*comm,Teuchos::REDUCE_MAX,Teuchos::as<int>(cntUnknownDofGIDs),&lTranslatedDofGIDs[0],&gTranslatedDofGIDs[0]);

    for(size_t k=0; k < Teuchos::as<size_t>(ovlUnknownStatusGids.size()); k++) {
      origGID2newGID[ovlUnknownStatusGids[k]] = gTranslatedDofGIDs[k+cntUnknownOffset];
    }
    Teuchos::Array<GlobalOrdinal> ovlDomainMapArray;
    for(size_t i = 0; i<input.getNodeNumElements(); i++) {
      GlobalOrdinal gcid = input.getGlobalElement(i);
      ovlDomainMapArray.push_back(origGID2newGID[gcid]);
    }
    RCP<Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> > ovlDomainMap =
        Xpetra::MapFactory<LocalOrdinal, GlobalOrdinal, Node>::Build
        (nonOvlInput.lib(),Teuchos::OrdinalTraits<GlobalOrdinal>::invalid(),ovlDomainMapArray(),0,comm);
    return ovlDomainMap;
  }



};

} // end namespace Xpetra

#define XPETRA_MAPUTILS_SHORT

#endif // PACKAGES_XPETRA_SUP_MAP_UTILS_HPP_
