// @HEADER
// ***********************************************************************
// 
//          Tpetra: Templated Linear Algebra Services Package
//                 Copyright (2004) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
// @HEADER

// TODO: make sure that Ordinal values in constructors aren't invalid()

#ifndef TPETRA_MAP_HPP
#define TPETRA_MAP_HPP

#include <Teuchos_TypeNameTraits.hpp>
#include <Teuchos_OrdinalTraits.hpp>
#include <Teuchos_CommHelpers.hpp>
#include <Teuchos_as.hpp>
#include "Tpetra_MapDecl.hpp"
#include "Tpetra_MapData.hpp"
#include "Tpetra_Directory.hpp"

//
// compute isDistributed. it will be global.
// min/max GID are always computed (from indexBase and numGlobal), and are therefore globally coherent as long as deps are.
// indexBase and numGlobal must always be verified.
// isContiguous is true for the "easy" constructors, assume false for the expert constructor
// 
// so we explicitly check    : isCont, numGlobal, indexBase
// then we explicitly compute: MinAllGID, MaxAllGID
// Data explicitly checks    : isDistributed

namespace Tpetra {

  template<class LocalOrdinal, class GlobalOrdinal>
  Map<LocalOrdinal,GlobalOrdinal>::Map(GlobalOrdinal numGlobalEntries, Teuchos_Ordinal indexBase, const Teuchos::RCP<const Teuchos::Comm<int> > &comm, bool local) 
    : MapData_()
  {
    // distribute the entries across the nodes so that they are 
    // - non-overlapping
    // - contiguous
    // - as evenly distributed as possible
    const GlobalOrdinal GONE = Teuchos::OrdinalTraits<GlobalOrdinal>::one();
    const GlobalOrdinal GZERO = Teuchos::OrdinalTraits<GlobalOrdinal>::zero();

    std::string errPrefix;
    errPrefix = Teuchos::typeName(*this) + "::constructor(numGlobal,indexBase,comm,local): ";

    if (local == false) 
    {
      int numImages = comm->getSize();
      int myImageID = comm->getRank();

      // check that numGlobalEntries,indexBase is equivalent across images
      std::vector<GlobalOrdinal> root_entries(2);
      root_entries[0] = numGlobalEntries;
      root_entries[1] = indexBase;
      Teuchos::broadcast<int,GlobalOrdinal>(*comm,0,2,&root_entries[0]);   // broadcast 2 ordinals from node 0
      int localChecks[2], globalChecks[2];
      localChecks[0] = -1;  // fail or pass
      localChecks[1] = 0;    // fail reason
      if (numGlobalEntries != root_entries[0]) {
        localChecks[0] = myImageID;
        localChecks[1] = 1;
      }
      else if (indexBase != root_entries[1]) {
        localChecks[0] = myImageID;
        localChecks[1] = 2;
      }
      // REDUCE_MAX will give us the image ID of the highest rank proc that DID NOT pass, as well as the reason
      // these will be -1 and 0 if all procs passed
      Teuchos::reduceAll<int,int>(*comm,Teuchos::REDUCE_MAX,2,localChecks,globalChecks);
      if (globalChecks[0] != -1) {
        if (globalChecks[1] == 1) {
          TEST_FOR_EXCEPTION(true,std::invalid_argument,
              errPrefix << "numGlobal must be the same on all nodes (examine node " << globalChecks[0] << ").");
        }
        else if (globalChecks[1] == 2) {
          TEST_FOR_EXCEPTION(true,std::invalid_argument,
              errPrefix << "indexBase must be the same on all nodes (examine node " << globalChecks[0] << ").");
        }
        else {
          // logic error on my part
          TEST_FOR_EXCEPTION(true,std::logic_error,
              errPrefix << "logic error. Please contact the Tpetra team.");
        }
      }
      // numgGlobalEntries is coherent, but is it valid? this comparison looks funny, but it avoids compiler warnings on unsigned types.
      TEST_FOR_EXCEPTION(numGlobalEntries < GONE && numGlobalEntries != GZERO, std::invalid_argument,
          errPrefix << "numGlobal == " << numGlobalEntries << ". Must be >= 0.");

      /* compute numLocalEntries
         We can write numGlobalEntries as follows:
         numGlobalEntries == numImages * B + remainder
         Each image is allocated entries as follows:
         [ B+1    iff myImageID <  remainder
         numLocalEntries = [
         [ B      iff myImageID >= remainder
         In the case that remainder == 0, then all images fall into the 
         latter case: numLocalEntries == B == numGlobalEntries / numImages
         It can then be shown that 
         numImages
         \Sum      numLocalEntries_i  == numGlobalEntries
         i=0
         This strategy is simple, requires no communication, and is optimal vis-a-vis
         uniform distribution of entries.
         This strategy is valid for any value of numGlobalEntries and numImages, 
         including the following border cases:
         - numImages == 1         -> remainder == 0 && numGlobalEntries == numLocalEntries
         - numEntries < numImages -> remainder == numGlobalEntries && numLocalEntries \in [0,1]
       */
      LocalOrdinal numLocalEntries = Teuchos::as<LocalOrdinal>(numGlobalEntries / numImages);
      LocalOrdinal remainder = Teuchos::as<LocalOrdinal>(numGlobalEntries % numImages);
#ifdef HAVE_TEUCHOS_DEBUG
      // the above code assumes truncation. is that safe?
      SHARED_TEST_FOR_EXCEPTION(numLocalEntries * numImages + remainder != numGlobalEntries,
          std::logic_error, "Tpetra::Map::constructor(numGlobal,indexBase,platform): GlobalOrdinal does not implement division with truncation."
          << " Please contact Tpetra team.",*comm);
#endif
      GlobalOrdinal start_index;
      if (myImageID < remainder) {
        ++numLocalEntries;
        /* the myImageID images before were each allocated 
           numGlobalEntries/numImages+1
           ergo, my offset is:
           myImageID * (numGlobalEntries/numImages+1)
           == myImageID * numLocalEntries
         */
        start_index = myImageID * numLocalEntries;
      }
      else {
        /* a quantity (equal to remainder) of the images before me
           were each allocated 
           numGlobalEntries/numImages+1
           entries. a quantity (equal to myImageID-remainder) of the remaining 
           images before me were each allocated 
           numGlobalEntries/numImages
           entries. ergo, my offset is:
           remainder*(numGlobalEntries/numImages+1) + (myImageID-remainder)*numGlobalEntries/numImages
           == remainder*numLocalEntries + remainder + myImageID*numLocalEntries - remainder*numLocalEntries
           == myImageID*numLocalEntries + remainder
         */
        start_index = myImageID*numLocalEntries + remainder;
      }

      // create empty maps between local and global entries: let the MapData constructor fill them
      Teuchos::ArrayRCP<GlobalOrdinal> lgMap = Teuchos::null;
      std::map<GlobalOrdinal,LocalOrdinal> glMap;

      // compute the min/max global IDs
      GlobalOrdinal minAllGID = indexBase;
      GlobalOrdinal maxAllGID = indexBase + numGlobalEntries - GONE;
      GlobalOrdinal minMyGID  = start_index + indexBase;
      GlobalOrdinal maxMyGID  = minMyGID + numLocalEntries - GONE;

      // create a MapData structure 
      MapData_ = Teuchos::rcp( new MapData<LocalOrdinal,GlobalOrdinal>(indexBase,numGlobalEntries,numLocalEntries,
            minAllGID,maxAllGID,minMyGID,maxMyGID,
            lgMap,glMap,true,  // we are contiguous
            comm) );

      // initialize the directory
      directorySetup();
    }
    else 
    {
      Teuchos::ArrayRCP<GlobalOrdinal> lgMap = Teuchos::null;
      std::map<GlobalOrdinal,LocalOrdinal>  glMap;

      // compute the min/max global IDs
      GlobalOrdinal minAllGID = indexBase;
      GlobalOrdinal maxAllGID = indexBase + numGlobalEntries - GONE;
      GlobalOrdinal minMyGID  = minAllGID;
      GlobalOrdinal maxMyGID  = maxAllGID;

      // create a MapData structure 
      MapData_ = Teuchos::rcp( new MapData<LocalOrdinal,GlobalOrdinal>(indexBase,numGlobalEntries,numGlobalEntries,
            minAllGID,maxAllGID,minMyGID,maxMyGID,
            lgMap,glMap,true,  // we are contiguous
            comm,true) );
    }
  }

  template<class LocalOrdinal, class GlobalOrdinal>
  Map<LocalOrdinal,GlobalOrdinal>::Map(GlobalOrdinal numGlobalEntries, LocalOrdinal numLocalEntries, Teuchos_Ordinal indexBase, 
                                       const Teuchos::RCP<const Teuchos::Comm<int> > &comm)
    : MapData_()
  {
    // Distribute the entries across the nodes so that they are 
    // - non-overlapping
    // - contiguous
    // This differs from Map(Ord,Ord,Plat) in that the user has specified the number of entries 
    // per node, so that they are not (necessarily) evenly distributed

    const LocalOrdinal LONE  = Teuchos::OrdinalTraits<LocalOrdinal>::one();
    const LocalOrdinal LZERO = Teuchos::OrdinalTraits<LocalOrdinal>::zero();
    const GlobalOrdinal GONE  = Teuchos::OrdinalTraits<GlobalOrdinal>::one();
    const GlobalOrdinal GZERO = Teuchos::OrdinalTraits<GlobalOrdinal>::zero();
    const GlobalOrdinal GINVALID = Teuchos::OrdinalTraits<GlobalOrdinal>::invalid();

    std::string errPrefix;
    errPrefix = Teuchos::typeName(*this) + "::constructor(numGlobal,numLocal,indexBase,platform): ";

    // get a internodal communicator from the Platform
    int myImageID = comm->getRank();
    // for communicating failures 
    int localChecks[2], globalChecks[2];

    /* compute the global size 
       we are computing the number of global entries because exactly ONE of the following is true:
       - the user didn't specify it, and we need it
       - the user did specify it, but we need to 
         + validate it against the sum of the local sizes, and
         + ensure that it is the same on all nodes
     */
    GlobalOrdinal global_sum;
    Teuchos::reduceAll<int,GlobalOrdinal>(*comm,Teuchos::REDUCE_SUM,Teuchos::as<GlobalOrdinal>(numLocalEntries),&global_sum);
    /* there are three errors we should be detecting:
       - numGlobalEntries != invalid() and it is incorrect/invalid
       - numLocalEntries invalid (<0)
    */
    localChecks[0] = -1;
    localChecks[1] = 0;
    if (numLocalEntries < LONE && numLocalEntries != LZERO) {
      // invalid
      localChecks[0] = myImageID;
      localChecks[1] = 1;
    }
    else if (numGlobalEntries < GONE && numGlobalEntries != GZERO && numGlobalEntries != GINVALID) {
      // invalid
      localChecks[0] = myImageID;
      localChecks[1] = 2;
    }
    else if (numGlobalEntries != GINVALID && numGlobalEntries != global_sum) {
      // incorrect
      localChecks[0] = myImageID;
      localChecks[1] = 3;
    }
    // now check that indexBase is equivalent across images
    GlobalOrdinal root_indexbase = indexBase;
    Teuchos::broadcast<int,GlobalOrdinal>(*comm,0,&root_indexbase);   // broadcast one ordinal from node 0
    if (indexBase != root_indexbase) {
      localChecks[0] = myImageID;
      localChecks[1] = 4;
    }
    // REDUCE_MAX will give us the image ID of the highest rank proc that DID NOT pass
    // this will be -1 if all procs passed
    Teuchos::reduceAll<int,int>(*comm,Teuchos::REDUCE_MAX,2,localChecks,globalChecks);
    if (globalChecks[0] != -1) {
      if (globalChecks[1] == 1) {
        TEST_FOR_EXCEPTION(true,std::invalid_argument,
            errPrefix << "numLocal is not valid on at least one node (possibly node " 
            << globalChecks[0] << ").");
      }
      else if (globalChecks[1] == 2) {
        TEST_FOR_EXCEPTION(true,std::invalid_argument,
            errPrefix << "numGlobal is not valid on at least one node (possibly node " 
            << globalChecks[0] << ").");
      }
      else if (globalChecks[1] == 3) {
        TEST_FOR_EXCEPTION(true,std::invalid_argument,
            errPrefix << "numGlobal doesn't match sum of numLocal (== " 
            << global_sum << ") on at least one node (possibly node " 
            << globalChecks[0] << ").");
      }
      else if (globalChecks[1] == 4) {
        TEST_FOR_EXCEPTION(true,std::invalid_argument,
            errPrefix << "indexBase is not the same on all nodes (examine node " 
            << globalChecks[0] << ").");
      }
      else {
        // logic error on my part
        TEST_FOR_EXCEPTION(true,std::logic_error,
            errPrefix << "logic error. Please contact the Tpetra team.");
      }
    }
    // set numGlobalEntries
    if (numGlobalEntries == GINVALID) {
      numGlobalEntries = global_sum;
    }

    // compute my local offset
    GlobalOrdinal start_index;
    Teuchos::scan<int,GlobalOrdinal>(*comm,Teuchos::REDUCE_SUM,numLocalEntries,&start_index);
    start_index -= numLocalEntries;

    // create empty maps between local and global entries: let the MapData constructor fill them
    Teuchos::ArrayRCP<GlobalOrdinal> lgMap = Teuchos::null;
    std::map<GlobalOrdinal,LocalOrdinal> glMap;

    // compute the min/max global IDs
    GlobalOrdinal minAllGID = indexBase;
    GlobalOrdinal maxAllGID = indexBase + numGlobalEntries - GONE;
    GlobalOrdinal minMyGID = start_index + indexBase;
    GlobalOrdinal maxMyGID = minMyGID + numLocalEntries - GONE;

    // create a MapData structure 
    MapData_ = Teuchos::rcp( new MapData<LocalOrdinal,GlobalOrdinal>(indexBase,numGlobalEntries,numLocalEntries,
                                                      minAllGID,maxAllGID,minMyGID,maxMyGID,
                                                      lgMap,glMap,true,  // we are contiguous
                                                      comm) );

    // initialize the directory
    directorySetup();
  }

  template<class LocalOrdinal, class GlobalOrdinal>
  Map<LocalOrdinal,GlobalOrdinal>::Map (GlobalOrdinal numGlobalEntries, const Teuchos::ArrayView<const GlobalOrdinal> &entryList, 
                     Teuchos_Ordinal indexBase, const Teuchos::RCP<const Teuchos::Comm<int> > &comm)
    : MapData_()
  {
    // Distribute the entries across the nodes in an arbitrary user-specified manner
    // They are not necessarily contiguous or evenly distributed
    const LocalOrdinal LZERO = Teuchos::OrdinalTraits<LocalOrdinal>::zero();
    const GlobalOrdinal GONE  = Teuchos::OrdinalTraits<GlobalOrdinal>::one();
    const GlobalOrdinal GZERO = Teuchos::OrdinalTraits<GlobalOrdinal>::zero();
    const GlobalOrdinal GINVALID = Teuchos::OrdinalTraits<GlobalOrdinal>::invalid();

    LocalOrdinal numLocalEntries = Teuchos::as<LocalOrdinal>(entryList.size());

    std::string errPrefix;
    errPrefix = Teuchos::typeName(*this) + "::constructor(numGlobal,entryList,indexBase,platform): ";

    int myImageID = comm->getRank();
    // for communicating failures 
    int localChecks[2], globalChecks[2];

    /* compute the global size 
       we are computing the number of global entries because exactly ONE of the following is true:
       - the user didn't specify it, and we need it
       - the user did specify it, but we need to 
         + validate it against the sum of the local sizes, and
         + ensure that it is the same on all nodes
     */
    GlobalOrdinal global_sum;
    Teuchos::reduceAll<int,GlobalOrdinal>(*comm,Teuchos::REDUCE_SUM,numLocalEntries,&global_sum);
    localChecks[0] = -1;
    localChecks[1] = 0;
    if (numGlobalEntries < GONE && numGlobalEntries != GZERO && numGlobalEntries != GINVALID) {
      // invalid
      localChecks[0] = myImageID;
      localChecks[1] = 1;
    }
    else if (numGlobalEntries != GINVALID && numGlobalEntries != global_sum) {
      // incorrect
      localChecks[0] = myImageID;
      localChecks[1] = 2;
    }
    // now check that indexBase is equivalent across images
    GlobalOrdinal root_indexbase = indexBase;
    Teuchos::broadcast<int,GlobalOrdinal>(*comm,0,&root_indexbase);   // broadcast one ordinal from node 0
    if (indexBase != root_indexbase) {
      localChecks[0] = myImageID;
      localChecks[1] = 3;
    }
    // REDUCE_MAX will give us the image ID of the highest rank proc that DID NOT pass
    // this will be -1 if all procs passed
    Teuchos::reduceAll<int,int>(*comm,Teuchos::REDUCE_MAX,2,localChecks,globalChecks);
    if (globalChecks[0] != -1) {
      if (globalChecks[1] == 1) {
        TEST_FOR_EXCEPTION(true,std::invalid_argument,
            errPrefix << "numGlobal is not valid on at least one node (possibly node "
            << globalChecks[0] << ").");
      }
      else if (globalChecks[1] == 2) {
        TEST_FOR_EXCEPTION(true,std::invalid_argument,
            errPrefix << "numGlobal doesn't match sum of numLocal (" 
            << global_sum << ") on at least one node (possibly node "
            << globalChecks[0] << ").");
      }
      else if (globalChecks[1] == 3) {
        TEST_FOR_EXCEPTION(true,std::invalid_argument,
            errPrefix << "indexBase is not the same on all nodes (possibly node "
            << globalChecks[0] << ").");
      }
      else {
        // logic error on my part
        TEST_FOR_EXCEPTION(true,std::logic_error,
            errPrefix << "logic error. Please contact the Tpetra team.");
      }
    }
    // set numGlobalEntries
    if (numGlobalEntries == GINVALID) {
      numGlobalEntries = global_sum;
    }

    // setup lgmap and glmap, and min/maxMyGIDs
    std::map<GlobalOrdinal,LocalOrdinal> glMap;
    // assume for now that there are numLocalEntries (there may be less, if some
    // GIDs are duplicated in entryList)
    Teuchos::ArrayRCP<GlobalOrdinal> lgMap;
    GlobalOrdinal minMyGID = indexBase;
    GlobalOrdinal maxMyGID = indexBase;
    // create the GID to LID map; do not assume GID in entryList are distinct.
    // in the case that a GID is duplicated, keep the previous LID
    // this is necessary so that LIDs are in [0,numLocal)
    // FINISH: make sure this is legal
    {
      Teuchos_Ordinal numLIDs = 0;
      if (numLocalEntries > LZERO) {
        lgMap = Teuchos::arcp<GlobalOrdinal>(numLocalEntries);
        for(Teuchos_Ordinal i=0; i < numLocalEntries; i++) {
          if (glMap.find(entryList[i]) == glMap.end()) {
            lgMap[numLIDs] = entryList[i];   // lgMap: LID to GID
            glMap[entryList[i]] = numLIDs;   // glMap: GID to LID
            numLIDs++;
          }
        }
        minMyGID = *std::min_element(entryList.begin(), entryList.end());
        maxMyGID = *std::max_element(entryList.begin(), entryList.end());
        // shrink lgMap appropriately
        if (numLocalEntries != Teuchos::as<Teuchos_Ordinal>(numLIDs)) {
          numLocalEntries = Teuchos::as<Teuchos_Ordinal>(numLIDs);
          lgMap = lgMap.persistingView(0,numLocalEntries);
        }
      }
    } 

    // set min/maxAllGIDs
    GlobalOrdinal minAllGID, maxAllGID;
    {
      GlobalOrdinal minmaxAllGIDlcl[2], minmaxAllGIDgbl[2];
      minmaxAllGIDlcl[0] = -minMyGID;  // negative allows us to do a single
      minmaxAllGIDlcl[1] =  maxMyGID;  // reduction below
      Teuchos::reduceAll<int,GlobalOrdinal>(*comm,Teuchos::REDUCE_MAX,2,minmaxAllGIDlcl,minmaxAllGIDgbl);
      minAllGID = -minmaxAllGIDgbl[0];
      maxAllGID =  minmaxAllGIDgbl[1];
    }
    TEST_FOR_EXCEPTION(minAllGID < indexBase, std::invalid_argument,
        errPrefix << "minimum GID == " << minAllGID << " is less than indexBase == " << indexBase);

    // call ESData constructor
    MapData_ = Teuchos::rcp(new MapData<LocalOrdinal,GlobalOrdinal>(indexBase, numGlobalEntries, 
          numLocalEntries, minAllGID, maxAllGID, 
          minMyGID, maxMyGID, lgMap, glMap, 
          false, comm));

    // initialize directory
    directorySetup();
  }

  template<class LocalOrdinal, class GlobalOrdinal>
  Map<LocalOrdinal,GlobalOrdinal>::Map (const Map<LocalOrdinal,GlobalOrdinal> &map) 
    : MapData_(map.MapData_)
  {}

  template<class LocalOrdinal, class GlobalOrdinal>
  Map<LocalOrdinal,GlobalOrdinal>::~Map () 
  {}

  template<class LocalOrdinal, class GlobalOrdinal>
  GlobalOrdinal Map<LocalOrdinal,GlobalOrdinal>::getNumGlobalEntries() const {
    return MapData_->numGlobalEntries_;
  }

  template<class LocalOrdinal, class GlobalOrdinal>
  LocalOrdinal Map<LocalOrdinal,GlobalOrdinal>::getNumMyEntries() const {
    return MapData_->numMyEntries_;
  }

  template<class LocalOrdinal, class GlobalOrdinal>
  Teuchos_Ordinal Map<LocalOrdinal,GlobalOrdinal>::getIndexBase() const {
    return MapData_->indexBase_;
  }

  template<class LocalOrdinal, class GlobalOrdinal>
  LocalOrdinal Map<LocalOrdinal,GlobalOrdinal>::getMinLocalIndex() const {
    return 0;
  }

  template<class LocalOrdinal, class GlobalOrdinal>
  LocalOrdinal Map<LocalOrdinal,GlobalOrdinal>::getMaxLocalIndex() const {
    return MapData_->numMyEntries_-1;
  }

  template<class LocalOrdinal, class GlobalOrdinal>
  GlobalOrdinal Map<LocalOrdinal,GlobalOrdinal>::getMinGlobalIndex() const {
    return MapData_->minMyGID_;
  }

  template<class LocalOrdinal, class GlobalOrdinal>
  GlobalOrdinal Map<LocalOrdinal,GlobalOrdinal>::getMaxGlobalIndex() const {
    return MapData_->maxMyGID_;
  }

  template<class LocalOrdinal, class GlobalOrdinal>
  GlobalOrdinal Map<LocalOrdinal,GlobalOrdinal>::getMinAllGlobalIndex() const {
    return MapData_->minAllGID_;
  }

  template<class LocalOrdinal, class GlobalOrdinal>
  GlobalOrdinal Map<LocalOrdinal,GlobalOrdinal>::getMaxAllGlobalIndex() const {
    return MapData_->maxAllGID_;
  }

  template<class LocalOrdinal, class GlobalOrdinal>
  LocalOrdinal Map<LocalOrdinal,GlobalOrdinal>::getLocalIndex(GlobalOrdinal globalIndex) const {
    if (MapData_->contiguous_) {
      if (globalIndex < getMinGlobalIndex() || globalIndex > getMaxGlobalIndex()) {
        return Teuchos::OrdinalTraits<LocalOrdinal>::invalid();
      }
      return Teuchos::as<LocalOrdinal>(globalIndex - getMinGlobalIndex());
    }
    else {
      typename std::map<GlobalOrdinal,LocalOrdinal>::const_iterator i;
      i = MapData_->glMap_.find(globalIndex);
      if (i == MapData_->glMap_.end()) {
        return Teuchos::OrdinalTraits<GlobalOrdinal>::invalid();
      }
      return i->second;
    }
  }

  template<class LocalOrdinal, class GlobalOrdinal>
  GlobalOrdinal Map<LocalOrdinal,GlobalOrdinal>::getGlobalIndex(LocalOrdinal localIndex) const {
    if (localIndex < getMinLocalIndex() || localIndex > getMaxLocalIndex()) {
      return Teuchos::OrdinalTraits<GlobalOrdinal>::invalid();
    }
    if (MapData_->contiguous_) {
      return getMinGlobalIndex() + localIndex;
    }
    else {
      return MapData_->lgMap_[localIndex];
    }
  }

  template<class LocalOrdinal, class GlobalOrdinal>
  bool Map<LocalOrdinal,GlobalOrdinal>::isMyLocalIndex(LocalOrdinal localIndex) const {
    if (getMinLocalIndex() <= localIndex && localIndex <= getMaxLocalIndex()) {
      return true;
    }
    return false;
  }

  template<class LocalOrdinal, class GlobalOrdinal>
  bool Map<LocalOrdinal,GlobalOrdinal>::isMyGlobalIndex(GlobalOrdinal globalIndex) const {
    if (MapData_->contiguous_) {
      return (getMinGlobalIndex() <= globalIndex) && (globalIndex <= getMaxGlobalIndex());
    }
    else {
      typename std::map<GlobalOrdinal,LocalOrdinal>::iterator i;
      i = MapData_->glMap_.find(globalIndex);
      return (i != MapData_->glMap_.end());
    }
  }

  template<class LocalOrdinal, class GlobalOrdinal>
  bool Map<LocalOrdinal,GlobalOrdinal>::isContiguous() const {
    return MapData_->contiguous_;
  }

  template<class LocalOrdinal, class GlobalOrdinal>
  bool Map<LocalOrdinal,GlobalOrdinal>::isCompatible (const Map< LocalOrdinal,GlobalOrdinal> &map) const {
    // check to make sure distribution is the same
    char iscompat_lcl;
    if (getNumGlobalEntries() != map.getNumGlobalEntries() ||
        getNumMyEntries() != map.getNumMyEntries()) 
    {
      // NOT compat on this node
      iscompat_lcl = 0;
    }
    else {
      // compat on this node
      iscompat_lcl = 1;
    }
    char iscompat_gbl;
    Teuchos::reduceAll<int,char>(*MapData_->comm_,Teuchos::REDUCE_MIN,iscompat_lcl,&iscompat_gbl);
    return(iscompat_gbl == 1);
  }

  template<class LocalOrdinal, class GlobalOrdinal>
  bool Map<LocalOrdinal,GlobalOrdinal>::isSameAs (const Map<LocalOrdinal,GlobalOrdinal> &map) const {
    if (MapData_ == map.MapData_) {
      // we should assume that this is globally coherent
      // if they share the same underlying MapData, then they must be equivalent maps
      return true;
    }

    // check all other globally coherent properties
    // if they do not share each of these properties, then they cannot be 
    // equivalent maps
    if ( (getMinGlobalIndex()   != map.getMinGlobalIndex())   ||
         (getMaxGlobalIndex()   != map.getMaxGlobalIndex())   ||
         (getNumGlobalEntries() != map.getNumGlobalEntries()) ||
         (isDistributed()       != map.isDistributed())       || 
         (getIndexBase()        != map.getIndexBase())        ){
      return false;
    }

    // If we get this far, we need to check local properties and the 
    // communicate same-ness across all nodes
    // we prefer local work over communication, ergo, we will perform all
    // comparisons and conduct a single communication
    char isSame_lcl = 1;

    // check number of entries owned by this node
    if (getNumMyEntries() != map.getNumMyEntries()) {
      isSame_lcl = 0;
    }

    // check the identity of the entries owned by this node
    // only do this if we haven't already determined not-same-ness
    if (isSame_lcl == 1) {
      // if they are contiguous, we can check the ranges easily
      // if they are not contiguous, we must check the individual LID -> GID mappings
      // the latter approach is valid in either case, but the former is faster
      if (MapData_->contiguous_ == true && map.MapData_->contiguous_ == true) {
        if ( (getMinGlobalIndex() != map.getMinGlobalIndex()) ||
             (getMaxGlobalIndex() != map.getMaxGlobalIndex()) ){
          isSame_lcl = 0;
        }
      }
      else {
        /* Note: std::equal requires that the latter range is as large as the former.
         * We know the ranges have equal length, because they have the same number of 
         * local entries. 
         * However, we only know that lgMap_ has been filled for the Map that is not
         * contiguous (one is potentially contiguous.) Calling getMyGlobalEntries()
         * will create it. */
        Teuchos::ArrayView<const GlobalOrdinal> ge1, ge2;
        ge1 =     getMyGlobalEntries();
        ge2 = map.getMyGlobalEntries();
        if (!std::equal(ge1.begin(),ge1.end(),ge2.begin())) {
          isSame_lcl = 0;
        }
      }
    }

    // now, determine if we detected not-same-ness on any node
    char isSame_gbl;
    Teuchos::reduceAll<int,char>(*MapData_->comm_,Teuchos::REDUCE_MIN,isSame_lcl,&isSame_gbl);
    return(isSame_gbl == 1);
  }

  template<class LocalOrdinal, class GlobalOrdinal>
  bool Map<LocalOrdinal,GlobalOrdinal>::operator== (const Map< LocalOrdinal,GlobalOrdinal > &map) const {
    return isSameAs(map);
  }

  template<class LocalOrdinal, class GlobalOrdinal>
  bool Map<LocalOrdinal,GlobalOrdinal>::operator!= (const Map< LocalOrdinal,GlobalOrdinal > &map) const {
    return !isSameAs(map);
  }

  template<class LocalOrdinal, class GlobalOrdinal>
  Map<LocalOrdinal,GlobalOrdinal>& Map<LocalOrdinal,GlobalOrdinal>::operator = (const Map<LocalOrdinal,GlobalOrdinal> & Source) {
    MapData_ = Source.MapData_;
    return *this;
  }

  template<class LocalOrdinal, class GlobalOrdinal>
  Teuchos::ArrayView<const GlobalOrdinal>
  Map<LocalOrdinal,GlobalOrdinal>::getMyGlobalEntries() const {
    const LocalOrdinal LZERO = Teuchos::OrdinalTraits<LocalOrdinal>::zero();
    // check to see if lgMap is empty
    // if so (and we have local entries), then fill it.
    if (MapData_->lgMap_ == Teuchos::null && MapData_->numMyEntries_ > LZERO) {
      // this would have been set up for a non-contiguous map
#ifdef HAVE_TEUCHOS_DEBUG
      TEST_FOR_EXCEPTION(MapData_->contiguous_ != true, std::logic_error,
          "Tpetra::Map::getMyGlobalEntries: logic error. Please notify the Tpetra team.");
#endif
      Teuchos::ArrayRCP<GlobalOrdinal> lgMap = Teuchos::arcp<GlobalOrdinal>(MapData_->numMyEntries_);
      for (GlobalOrdinal gid=MapData_->minMyGID_; gid <= MapData_->maxMyGID_; ++gid) {
        *(lgMap++) = gid;
      }
      MapData_->lgMap_ = lgMap;
    }
    return MapData_->lgMap_();
  }

  template<class LocalOrdinal, class GlobalOrdinal>
  bool Map<LocalOrdinal,GlobalOrdinal>::isDistributed() const {
    return MapData_->distributed_;
  }

  template<class LocalOrdinal, class GlobalOrdinal>
  void Map<LocalOrdinal,GlobalOrdinal>::print(std::ostream& os) const {
    const LocalOrdinal nME = getNumMyEntries();
    
    using std::endl;

    getMyGlobalEntries(); // neglect the output, we call this here only to ensure that sure list is generated
    int myImageID = MapData_->comm_->getRank();
    int numImages = MapData_->comm_->getSize();

    for(int imageCtr = 0; imageCtr < numImages; ++imageCtr) {
      if(myImageID == imageCtr) {
        if(myImageID == 0) { // this is the root node (only output this info once)
          os << endl 
             << "Number of Global Entries = " << getNumGlobalEntries()  << endl
             << "Maximum of all GIDs      = " << getMaxAllGlobalIndex() << endl
             << "Minimum of all GIDs      = " << getMinAllGlobalIndex() << endl
             << "Index Base               = " << getIndexBase()         << endl;
        }
        os << endl;
        os << "Number of Local Elements   = " << nME           << endl
           << "Maximum of my GIDs         = " << getMaxGlobalIndex() << endl
           << "Minimum of my GIDs         = " << getMinGlobalIndex() << endl;
        os << endl;
        os << std::setw(16) << "ImageID"
           << std::setw(16) << "Local Index"
           << std::setw(16) << "Global Index"
           << endl;
        Teuchos::ArrayView<const GlobalOrdinal> myEntries = getMyGlobalEntries();
        for(Teuchos_Ordinal i=0; i < nME; i++) {
          os << std::setw(16) << myImageID 
             << std::setw(16) << i
             << std::setw(16) << myEntries[i]
             << endl;
        }
        os << std::flush;
      }
      // Do a few global ops to give I/O a chance to complete
      MapData_->comm_->barrier();
      MapData_->comm_->barrier();
      MapData_->comm_->barrier();
    }
  }

  template<class LocalOrdinal, class GlobalOrdinal>
  void Map<LocalOrdinal,GlobalOrdinal>::directorySetup() {
    if (getNumGlobalEntries() != Teuchos::OrdinalTraits<GlobalOrdinal>::zero()) {
      if (MapData_->directory_ == Teuchos::null) {
        MapData_->directory_ = Teuchos::rcp( new Directory<LocalOrdinal,GlobalOrdinal>(*this) );
      }
    }
  }

  template<class LocalOrdinal, class GlobalOrdinal>
  bool Map<LocalOrdinal,GlobalOrdinal>::getRemoteIndexList(
      const Teuchos::ArrayView<const GlobalOrdinal> & GIDList, 
      const Teuchos::ArrayView<int> & imageIDList, 
      const Teuchos::ArrayView<LocalOrdinal> & LIDList) const 
  {
    if (GIDList.size() == 0) return false;
    TEST_FOR_EXCEPTION(getNumGlobalEntries() == 0, std::runtime_error,
        Teuchos::typeName(*this) << "::getRemoteIndexList(): getRemoteIndexList() cannot be called, zero entries on node.");
    return MapData_->directory_->getDirectoryEntries(GIDList, imageIDList, LIDList);
  }

  template<class LocalOrdinal, class GlobalOrdinal>
  bool Map<LocalOrdinal,GlobalOrdinal>::getRemoteIndexList(
      const Teuchos::ArrayView<const GlobalOrdinal> & GIDList, 
      const Teuchos::ArrayView<int> & imageIDList) const 
  {
    if (GIDList.size() == 0) return false;
    TEST_FOR_EXCEPTION(getNumGlobalEntries() == 0, std::runtime_error,
        Teuchos::typeName(*this) << "::getRemoteIndexList(): getRemoteIndexList() cannot be called, zero entries on node.");
    return MapData_->directory_->getDirectoryEntries(GIDList, imageIDList);
  }

  template<class LocalOrdinal, class GlobalOrdinal>
  Teuchos::RCP<const Teuchos::Comm<int> >
  Map<LocalOrdinal,GlobalOrdinal>::getComm() const {
    return MapData_->comm_;
  }

} // Tpetra namespace

#endif // TPETRA_MAP_HPP

