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

#ifndef TPETRA_MAP_HPP
#define TPETRA_MAP_HPP

#include <Teuchos_OrdinalTraits.hpp>
#include <Teuchos_CommHelpers.hpp>
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

  template<typename Ordinal>
  Map<Ordinal>::Map(Ordinal numGlobalEntries, Ordinal indexBase, const Platform<Ordinal> &platform, bool local) 
    : Teuchos::Object(local == false ? "Tpetra::Map" : "Teptra::Map(local)")
    , MapData_()
  {
    // distribute the entries across the nodes so that they are 
    // - non-overlapping
    // - contiguous
    // - as evenly distributed as possible
    const Ordinal ONE = Teuchos::OrdinalTraits<Ordinal>::one();
    const Ordinal ZERO = Teuchos::OrdinalTraits<Ordinal>::zero();

    std::string errPrefix;
    errPrefix = "Tpetra::Map<" + Teuchos::OrdinalTraits<Ordinal>::name() 
              + ">::constructor(numGlobal,indexBase,platform,local): ";

    // get a internodal communicator from the Platform
    Teuchos::RCP< Teuchos::Comm<Ordinal> > comm = platform.createComm();
    if (local == false) 
    {
      Ordinal numImages = comm->getSize();
      Ordinal myImageID = comm->getRank();

      // check that numGlobalEntries,indexBase is equivalent across images
      std::vector<Ordinal> root_entries(2);
      root_entries[0] = numGlobalEntries;
      root_entries[1] = indexBase;
      Teuchos::broadcast(*comm,(Ordinal)0,(Ordinal)2,&root_entries[0]);   // broadcast 2 ordinals from node 0
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
      Teuchos::reduceAll<Ordinal,int>(*comm,Teuchos::REDUCE_MAX,2,localChecks,globalChecks);
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
      TEST_FOR_EXCEPTION(numGlobalEntries < ONE && numGlobalEntries != ZERO, std::invalid_argument,
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
      Ordinal numLocalEntries = numGlobalEntries / numImages;    
      Ordinal remainder = numGlobalEntries % numImages;
#ifdef TEUCHOS_DEBUG
      // the above code assumes truncation. is that safe?
      SHARED_TEST_FOR_EXCEPTION(numLocalEntries * numImages + remainder != numGlobalEntries,
          std::logic_error, "Tpetra::Map::constructor(numGlobal,indexBase,platform): Ordinal does not implement division with truncation."
          << " Please contact Tpetra team.",*comm);
#endif
      Ordinal start_index;
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
      Teuchos::ArrayRCP<Ordinal> lgMap = Teuchos::null;
      std::map<Ordinal,Ordinal> glMap;

      // compute the min/max global IDs
      Ordinal minAllGID = indexBase;
      Ordinal maxAllGID = indexBase + numGlobalEntries - ONE;
      Ordinal minMyGID  = start_index + indexBase;
      Ordinal maxMyGID  = minMyGID + numLocalEntries - ONE;

      Teuchos::RCP< Platform<Ordinal> > platform_clone = platform.clone();

      // create a MapData structure 
      MapData_ = Teuchos::rcp( new MapData<Ordinal>(indexBase,numGlobalEntries,numLocalEntries,
            minAllGID,maxAllGID,minMyGID,maxMyGID,
            lgMap,glMap,true,  // we are contiguous
            platform_clone, comm) );

      // initialize the directory
      directorySetup();
    }
    else 
    {
      Teuchos::ArrayRCP<Ordinal> lgMap = Teuchos::null;
      std::map<Ordinal,Ordinal>  glMap;

      // compute the min/max global IDs
      Ordinal minAllGID = indexBase;
      Ordinal maxAllGID = indexBase + numGlobalEntries - ONE;
      Ordinal minMyGID  = minAllGID;
      Ordinal maxMyGID  = maxAllGID;

      Teuchos::RCP< Platform<Ordinal> > platform_clone = platform.clone();

      // create a MapData structure 
      MapData_ = Teuchos::rcp( new MapData<Ordinal>(indexBase,numGlobalEntries,numGlobalEntries,
            minAllGID,maxAllGID,minMyGID,maxMyGID,
            lgMap,glMap,true,  // we are contiguous
            platform_clone, comm,true) );
    }
  }

  template<typename Ordinal>
  Map<Ordinal>::Map (Ordinal numGlobalEntries, Ordinal numLocalEntries, Ordinal indexBase, 
                         const Platform<Ordinal> &platform)
    : Teuchos::Object("Tpetra::Map")
    , MapData_()
  {
    // Distribute the entries across the nodes so that they are 
    // - non-overlapping
    // - contiguous
    // This differs from Map(Ord,Ord,Plat) in that the user has specified the number of entries 
    // per node, so that they are not (necessarily) evenly distributed

    const Ordinal ONE = Teuchos::OrdinalTraits<Ordinal>::one();
    const Ordinal ZERO = Teuchos::OrdinalTraits<Ordinal>::zero();
    const Ordinal INVALID = Teuchos::OrdinalTraits<Ordinal>::invalid();

    std::string errPrefix;
    errPrefix = "Tpetra::Map<" + Teuchos::OrdinalTraits<Ordinal>::name() 
              + ">::constructor(numGlobal,numLocal,indexBase,platform): ";

    // get a internodal communicator from the Platform
    Teuchos::RCP< Teuchos::Comm<Ordinal> > comm = platform.createComm();
    Ordinal myImageID = comm->getRank();
    // for communicating failures 
    int localChecks[2], globalChecks[2];

    /* compute the global size 
       we are computing the number of global entries because exactly ONE of the following is true:
       - the user didn't specify it, and we need it
       - the user did specify it, but we need to 
         + validate it against the sum of the local sizes, and
         + ensure that it is the same on all nodes
     */
    Ordinal global_sum;
    Teuchos::reduceAll(*comm,Teuchos::REDUCE_SUM,numLocalEntries,&global_sum);
    /* there are three errors we should be detecting:
       - numGlobalEntries != invalid() and it is incorrect/invalid
       - numLocalEntries invalid (<0)
    */
    localChecks[0] = -1;
    if (numLocalEntries < ONE && numLocalEntries != ZERO) {
      // invalid
      localChecks[0] = myImageID;
      localChecks[1] = 1;
    }
    else if (numGlobalEntries < ONE && numGlobalEntries != ZERO && numGlobalEntries != INVALID) {
      // invalid
      localChecks[0] = myImageID;
      localChecks[1] = 2;
    }
    else if (numGlobalEntries != INVALID && numGlobalEntries != global_sum) {
      // incorrect
      localChecks[0] = myImageID;
      localChecks[1] = 3;
    }
    // now check that indexBase is equivalent across images
    Ordinal root_indexbase = indexBase;
    Teuchos::broadcast(*comm,0,&root_indexbase);   // broadcast one ordinal from node 0
    if (indexBase != root_indexbase) {
      localChecks[0] = myImageID;
      localChecks[1] = 4;
    }
    // REDUCE_MAX will give us the image ID of the highest rank proc that DID NOT pass
    // this will be -1 if all procs passed
    Teuchos::reduceAll<Ordinal,int>(*comm,Teuchos::REDUCE_MAX,2,localChecks,globalChecks);
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
    if (numGlobalEntries == INVALID) {
      numGlobalEntries = global_sum;
    }

    // compute my local offset
    Ordinal start_index;
    Teuchos::scan(*comm,Teuchos::REDUCE_SUM,numLocalEntries,&start_index);
    start_index -= numLocalEntries;

    // create empty maps between local and global entries: let the MapData constructor fill them
    Teuchos::ArrayRCP<Ordinal> lgMap = Teuchos::null;
    std::map<Ordinal,Ordinal> glMap;

    // compute the min/max global IDs
    Ordinal minAllGID = indexBase;
    Ordinal maxAllGID = indexBase + numGlobalEntries - ONE;
    Ordinal minMyGID = start_index + indexBase;
    Ordinal maxMyGID = minMyGID + numLocalEntries - ONE;

    Teuchos::RCP< Platform<Ordinal> > platform_clone = platform.clone();

    // create a MapData structure 
    MapData_ = Teuchos::rcp( new MapData<Ordinal>(indexBase,numGlobalEntries,numLocalEntries,
                                                      minAllGID,maxAllGID,minMyGID,maxMyGID,
                                                      lgMap,glMap,true,  // we are contiguous
                                                      platform_clone, comm) );

    // initialize the directory
    directorySetup();
  }

  template<typename Ordinal>
  Map<Ordinal>::Map (Ordinal numGlobalEntries, const Teuchos::ArrayView<const Ordinal> &entryList, 
                     Ordinal indexBase, const Platform<Ordinal> &platform) 
    : Teuchos::Object("Tpetra::Map")
    , MapData_()
  {
    // Distribute the entries across the nodes in an arbitrary user-specified manner
    // They are not necessarily contiguous or evenly distributed
    const Ordinal ONE = Teuchos::OrdinalTraits<Ordinal>::one();
    const Ordinal ZERO = Teuchos::OrdinalTraits<Ordinal>::zero();
    const Ordinal INVALID = Teuchos::OrdinalTraits<Ordinal>::invalid();

    Ordinal numLocalEntries = entryList.size();

    std::string errPrefix;
    errPrefix = "Tpetra::Map<" + Teuchos::OrdinalTraits<Ordinal>::name() 
              + ">::constructor(numGlobal,entryList,indexBase,platform): ";

    // platform & comm setup
    Teuchos::RCP< Teuchos::Comm<Ordinal> > comm = platform.createComm();
    Ordinal myImageID = comm->getRank();
    // for communicating failures 
    int localChecks[2], globalChecks[2];

    /* compute the global size 
       we are computing the number of global entries because exactly ONE of the following is true:
       - the user didn't specify it, and we need it
       - the user did specify it, but we need to 
         + validate it against the sum of the local sizes, and
         + ensure that it is the same on all nodes
     */
    Ordinal global_sum;
    Teuchos::reduceAll(*comm,Teuchos::REDUCE_SUM,numLocalEntries,&global_sum);
    localChecks[0] = -1;
    if (numGlobalEntries < ONE && numGlobalEntries != ZERO && numGlobalEntries != INVALID) {
      // invalid
      localChecks[0] = myImageID;
      localChecks[1] = 1;
    }
    else if (numGlobalEntries != INVALID && numGlobalEntries != global_sum) {
      // incorrect
      localChecks[0] = myImageID;
      localChecks[1] = 2;
    }
    // now check that indexBase is equivalent across images
    Ordinal root_indexbase = indexBase;
    Teuchos::broadcast(*comm,0,&root_indexbase);   // broadcast one ordinal from node 0
    if (indexBase != root_indexbase) {
      localChecks[0] = myImageID;
      localChecks[1] = 3;
    }
    // REDUCE_MAX will give us the image ID of the highest rank proc that DID NOT pass
    // this will be -1 if all procs passed
    Teuchos::reduceAll<Ordinal,int>(*comm,Teuchos::REDUCE_MAX,2,localChecks,globalChecks);
    if (globalChecks[0] != INVALID) {
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
    if (numGlobalEntries == INVALID) {
      numGlobalEntries = global_sum;
    }

    // setup lgmap and glmap, and min/maxMyGIDs
    std::map<Ordinal,Ordinal> glMap;
    // assume for now that there are numLocalEntries (there may be less, if some
    // GIDs are duplicated in entryList)
    Teuchos::ArrayRCP<Ordinal> lgMap;
    Ordinal minMyGID = indexBase;
    Ordinal maxMyGID = indexBase;
    // create the GID to LID map; do not assume GID in entryList are distinct.
    // in the case that a GID is duplicated, keep the previous LID
    // this is necessary so that LIDs are in [0,numLocal)
    // FINISH: make sure this is legal
    {
      Ordinal numLIDs = ZERO;
      if (numLocalEntries > ZERO) {
        lgMap = Teuchos::arcp<Ordinal>(numLocalEntries);
        for(Ordinal i = ZERO; i < numLocalEntries; i++) {
          if (glMap.find(entryList[i]) == glMap.end()) {
            lgMap[numLIDs] = entryList[i];   // lgMap: LID to GID
            glMap[entryList[i]] = numLIDs;   // glMap: GID to LID
            numLIDs++;
          }
        }
        minMyGID = *std::min_element(entryList.begin(), entryList.end());
        maxMyGID = *std::max_element(entryList.begin(), entryList.end());
        // shrink lgMap appropriately
        if (numLocalEntries != numLIDs) {
          numLocalEntries = numLIDs;
          lgMap = lgMap.persistingView(0,numLocalEntries);
        }
      }
    } 

    // set min/maxAllGIDs
    Ordinal minAllGID, maxAllGID;
    {
      Ordinal minmaxAllGIDlcl[2], minmaxAllGIDgbl[2];
      minmaxAllGIDlcl[0] = -minMyGID;  // negative allows us to do a single
      minmaxAllGIDlcl[1] =  maxMyGID;  // reduction below
      Teuchos::reduceAll<Ordinal>(*comm,Teuchos::REDUCE_MAX,2,minmaxAllGIDlcl,minmaxAllGIDgbl);
      minAllGID = -minmaxAllGIDgbl[0];
      maxAllGID =  minmaxAllGIDgbl[1];
    }
    TEST_FOR_EXCEPTION(minAllGID < indexBase, std::invalid_argument,
        errPrefix << "minimum GID == " << minAllGID << " is less than indexBase == " << indexBase);

    Teuchos::RCP< Platform<Ordinal> > platform_clone = platform.clone();

    // call ESData constructor
    MapData_ = Teuchos::rcp(new MapData<Ordinal>(indexBase, numGlobalEntries, 
          numLocalEntries, minAllGID, maxAllGID, 
          minMyGID, maxMyGID, lgMap, glMap, 
          false, platform_clone, comm));

    // initialize directory
    directorySetup();
  }

  template<typename Ordinal>
  Map<Ordinal>::Map (const Map<Ordinal> &map) 
    : Teuchos::Object(map.label())
    , MapData_(map.MapData_)
  {}

  template<typename Ordinal>
  Map<Ordinal>::~Map () 
  { }

  template<typename Ordinal>
  Ordinal Map<Ordinal>::getNumGlobalEntries() const {
    return MapData_->numGlobalEntries_;
  }

  template<typename Ordinal>
  Ordinal Map<Ordinal>::getNumMyEntries() const {
    return MapData_->numMyEntries_;
  }

  template<typename Ordinal>
  Ordinal Map<Ordinal>::getIndexBase() const {
    return MapData_->indexBase_;
  }

  template<typename Ordinal>
  Ordinal Map<Ordinal>::getMinLocalIndex() const {
    return Teuchos::OrdinalTraits<Ordinal>::zero();
  }

  template<typename Ordinal>
  Ordinal Map<Ordinal>::getMaxLocalIndex() const {
    return MapData_->numMyEntries_-1;
  }

  template<typename Ordinal>
  Ordinal Map<Ordinal>::getMinGlobalIndex() const {
    return MapData_->minMyGID_;
  }

  template<typename Ordinal>
  Ordinal Map<Ordinal>::getMaxGlobalIndex() const {
    return MapData_->maxMyGID_;
  }

  template<typename Ordinal>
  Ordinal Map<Ordinal>::getMinAllGlobalIndex() const {
    return MapData_->minAllGID_;
  }

  template<typename Ordinal>
  Ordinal Map<Ordinal>::getMaxAllGlobalIndex() const {
    return MapData_->maxAllGID_;
  }

  template<typename Ordinal>
  Ordinal Map<Ordinal>::getLocalIndex(Ordinal globalIndex) const {
    if (MapData_->contiguous_) {
      TEST_FOR_EXCEPTION(
        globalIndex < getMinGlobalIndex() || globalIndex > getMaxGlobalIndex(), 
        std::invalid_argument,
        "Tpetra::Map<" << Teuchos::OrdinalTraits<Ordinal>::name() 
                       << ">::getLocalIndex(gid): gid does not belong to the map."
      );
      return globalIndex - getMinGlobalIndex();
    }
    else {
      typename std::map<Ordinal,Ordinal>::const_iterator i;
      i = MapData_->glMap_.find(globalIndex);
      TEST_FOR_EXCEPTION(
        i == MapData_->glMap_.end(), std::invalid_argument,
        "Tpetra::Map<" << Teuchos::OrdinalTraits<Ordinal>::name() 
                       << ">::getLocalIndex(gid): gid does not belong to the map."
      );
      return i->second;
    }
  }

  template<typename Ordinal>
  Ordinal Map<Ordinal>::getGlobalIndex(Ordinal localIndex) const {
    TEST_FOR_EXCEPTION(
        localIndex < 0 || localIndex > MapData_->numMyEntries_-1,
        std::invalid_argument,
        "Tpetra::Map<" << Teuchos::OrdinalTraits<Ordinal>::name() 
           << ">::getGlobalIndex(lid): lid not valid.");
    if (MapData_->contiguous_) {
      return localIndex + getMinGlobalIndex();
    }
    else {
      return MapData_->lgMap_[localIndex];
    }
  }

  template<typename Ordinal>
  bool Map<Ordinal>::isMyLocalIndex(Ordinal localIndex) const {
    // we have local indices in [0,numLocalEntries)
    if (localIndex >= 0 && localIndex < MapData_->numMyEntries_) {
      return true;
    }
    return false;
  }

  template<typename Ordinal>
  bool Map<Ordinal>::isMyGlobalIndex(Ordinal globalIndex) const {
    if (MapData_->contiguous_) {
      return (MapData_->minMyGID_ <= globalIndex) && (globalIndex <= MapData_->maxMyGID_);
    }
    else {
      typename std::map<Ordinal,Ordinal>::iterator i;
      i = MapData_->glMap_.find(globalIndex);
      return (i != MapData_->glMap_.end());
    }
  }

  template<typename Ordinal>
  bool Map<Ordinal>::isContiguous() const {
    return MapData_->contiguous_;
  }

  template<typename Ordinal>
  bool Map<Ordinal>::isCompatible (const Map< Ordinal> &map) const {
    // this is wrong: the non-compat procs will return false, then the compat 
    // procs will stall on this communication. this is demonstrated by tpetra/test/Map/Map_test.exe
    /*
    if (getNumGlobalEntries() != map.getNumGlobalEntries() ||
        getNumMyEntries() != map.getNumMyEntries()) 
    { return false; }
    */

    // check to make sure distribution is the same
    char iscompat_lcl;
    if (getNumGlobalEntries() != map.getNumGlobalEntries() ||
        getNumMyEntries() != map.getNumMyEntries()) 
    {
      // NOT compat on this node
      iscompat_lcl = Teuchos::ScalarTraits<char>::zero();
    }
    else {
      // compat on this node
      iscompat_lcl = Teuchos::ScalarTraits<char>::one();
    }
    char iscompat_gbl;
    Teuchos::reduceAll(*MapData_->comm_,Teuchos::REDUCE_MIN,iscompat_lcl,&iscompat_gbl);
    return(iscompat_gbl == Teuchos::ScalarTraits<char>::one());
  }

  template<typename Ordinal>
  bool Map<Ordinal>::isSameAs (const Map<Ordinal> &map) const {
    if (MapData_ == map.MapData_) {
      // we should assume that this is globally coherent
      // if they share the same underlying MapData, then they must be equivalent maps
      return true;
    }

    // check all other globally coherent properties
    // if they do not share each of these properties, then they cannot be 
    // equivalent maps
    if ( (getMinGlobalIndex()   != map.getMinGlobalIndex()) ||
         (getMaxGlobalIndex()   != map.getMaxGlobalIndex()) ||
         (getNumGlobalEntries() != map.getNumGlobalEntries()) ||
         (isDistributed()       != map.isDistributed()) ) {
      return false;
    }

    // If we get this far, we need to check local properties and the 
    // communicate same-ness across all nodes
    // we prefer local work over communication, ergo, we will perform all
    // comparisons and conduct a single communication
    char isSame_lcl = Teuchos::ScalarTraits<char>::one();

    // check number of entries owned by this node
    if (getNumMyEntries() != map.getNumMyEntries()) {
      isSame_lcl = Teuchos::ScalarTraits<char>::zero();
    }

    // check the identity of the entries owned by this node
    // only do this if we haven't already determined not-same-ness
    if (isSame_lcl == Teuchos::ScalarTraits<char>::one()) {
      // if they are contiguous, we can check the ranges easily
      // if they are not contiguous, we must check the individual LID -> GID mappings
      // the latter approach is valid in either case, but the former is faster
      if (MapData_->contiguous_ == true && map.MapData_->contiguous_ == true) {
        if (MapData_->minMyGID_ != map.MapData_->minMyGID_ ||
            MapData_->maxMyGID_ != map.MapData_->maxMyGID_) {
          isSame_lcl = Teuchos::ScalarTraits<char>::zero();
        }
      }
      else {
        if (!std::equal(    MapData_->lgMap_.begin(),MapData_->lgMap_.end(),
                        map.MapData_->lgMap_.begin())) {
          isSame_lcl = Teuchos::ScalarTraits<char>::zero();
        }
      }
    }

    // now, determine if we detected not-same-ness on any node
    char isSame_gbl;
    Teuchos::reduceAll(*MapData_->comm_,Teuchos::REDUCE_MIN,isSame_lcl,&isSame_gbl);
    return(isSame_gbl == Teuchos::ScalarTraits<char>::one());
  }

  template<typename Ordinal>
  bool Map<Ordinal>::operator== (const Map< Ordinal > &map) const {
    return isSameAs(map);
  }

  template<typename Ordinal>
  bool Map<Ordinal>::operator!= (const Map< Ordinal > &map) const {
    return !isSameAs(map);
  }

  template<typename Ordinal>
  Map<Ordinal>& Map<Ordinal>::operator = (const Map<Ordinal> & Source) {
    MapData_ = Source.MapData_;
    return *this;
  }

  template<typename Ordinal>
  Teuchos::ArrayView<const Ordinal>
  Map<Ordinal>::getMyGlobalEntries() const {
    const Ordinal ZERO = Teuchos::OrdinalTraits<Ordinal>::zero();
    // check to see if lgMap is empty
    // if so (and we have local entries), then fill it.
    if (MapData_->lgMap_ == Teuchos::null && MapData_->numMyEntries_ > ZERO) {
      // this would have been set up for a non-contiguous map
#ifdef TEUCHOS_DEBUG
      TEST_FOR_EXCEPTION(MapData_->contiguous_ != true, std::logic_error,
          "Tpetra::Map::getMyGlobalEntries: logic error. Please notify the Tpetra team.");
#endif
      Teuchos::ArrayRCP<Ordinal> lgMap = Teuchos::arcp<Ordinal>(MapData_->numMyEntries_);
      MapData_->lgMap_ = lgMap;
      for (Ordinal lid=MapData_->minMyGID_; lid <= MapData_->maxMyGID_; ++lid) {
        *(lgMap++) = lid;
      }
    }
    return MapData_->lgMap_();
  }

  template<typename Ordinal>
  bool Map<Ordinal>::isDistributed() const {
    return MapData_->distributed_;
  }

  template<typename Ordinal>
  void Map<Ordinal>::print(std::ostream& os) const {
    const Ordinal zero = Teuchos::OrdinalTraits<Ordinal>::zero();
    const Ordinal nME = getNumMyEntries();
    
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
        Teuchos::ArrayView<const Ordinal> myEntries = getMyGlobalEntries();
        for(Ordinal i = zero; i < nME; i++) {
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

  template<typename Ordinal>
  void Map<Ordinal>::directorySetup() {
    if (getNumGlobalEntries() != Teuchos::OrdinalTraits<Ordinal>::zero()) {
      if (MapData_->directory_ == Teuchos::null) {
        MapData_->directory_ = Teuchos::rcp( new Directory<Ordinal>(*this) );
      }
    }
  }

  template<typename Ordinal>
  bool Map<Ordinal>::getRemoteIndexList(
      const Teuchos::ArrayView<const Ordinal> & GIDList, 
      const Teuchos::ArrayView<Ordinal> & imageIDList, 
      const Teuchos::ArrayView<Ordinal> & LIDList) const 
  {
    if (GIDList.size() == 0) return false;
    TEST_FOR_EXCEPTION(getNumGlobalEntries() == 0, std::runtime_error,
        "Tpetra::Map<" + Teuchos::OrdinalTraits<Ordinal>::name() 
        + ">::getRemoteIndexList(): getRemoteIndexList() cannot be called, zero entries on node.");
    return MapData_->directory_->getDirectoryEntries(GIDList, imageIDList, LIDList);
  }

  template<typename Ordinal>
  bool Map<Ordinal>::getRemoteIndexList(
      const Teuchos::ArrayView<const Ordinal> & GIDList, 
      const Teuchos::ArrayView<Ordinal> & imageIDList) const 
  {
    if (GIDList.size() == 0) return false;
    TEST_FOR_EXCEPTION(getNumGlobalEntries() == 0, std::runtime_error,
        "Tpetra::Map<" + Teuchos::OrdinalTraits<Ordinal>::name() 
        + ">::getRemoteIndexList(): getRemoteIndexList() cannot be called, zero entries on node.");
    return MapData_->directory_->getDirectoryEntries(GIDList, imageIDList);
  }


  template<typename Ordinal>
  Teuchos::RCP<const Platform<Ordinal> >
  Map<Ordinal>::getPlatform() const {
    return MapData_->platform_;
  }

  template<typename Ordinal>
  Teuchos::RCP<const Teuchos::Comm<Ordinal> >
  Map<Ordinal>::getComm() const {
    return MapData_->comm_;
  }


} // Tpetra namespace

#endif // TPETRA_MAP_HPP

