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

//
// compute isDistributed. it will be global.
// min/max GID are always computed (from indexBase and numGlobal), and are therefore globally coherent as long as deps are.
// indexBase and numGlobal must always be verified.
// isContiguous is true for the "easy" constructors, assume false for the expert constructor
// 
// so we explicitly check    : isCont, numGlobal, indexBase
// then we explicitly compute: MinAllGID, MaxAllGID
// Data explicitly checks    : isDistributed
//
// FINISH:
// might consider allowing user to specify isDistributed, for supporting local maps.
// the question is, why not just get rid of local map? think about this, and ask mike.
// the point is, we would like to be able to create a local map that does not require any communication

namespace Tpetra {

  template<typename OrdinalType>
  Map<OrdinalType>::Map(OrdinalType numGlobalEntries, OrdinalType indexBase, const Platform<OrdinalType> &platform) 
    : Teuchos::Object("Tpetra::Map")
    , MapData_()
  {
    // distribute the entries across the nodes so that they are 
    // - non-overlapping
    // - contiguous
    // - as evenly distributed as possible
    const OrdinalType one = Teuchos::OrdinalTraits<OrdinalType>::one();
    const OrdinalType zero = Teuchos::OrdinalTraits<OrdinalType>::zero();
    const OrdinalType negOne = zero - one;

    std::string errPrefix;
    errPrefix = "Tpetra::Map<" + Teuchos::OrdinalTraits<OrdinalType>::name() 
              + ">::constructor(numGlobal,indexBase,platform): ";

    // get a internodal communicator from the Platform
    Teuchos::RCP< Teuchos::Comm<OrdinalType> > comm = platform.createComm();
    OrdinalType numImages = comm->getSize();
    OrdinalType myImageID = comm->getRank();

    // check that numGlobalEntries,indexBase is equivalent across images
    std::vector<OrdinalType> root_entries(2);
    root_entries[0] = numGlobalEntries;
    root_entries[1] = indexBase;
    Teuchos::broadcast(*comm,(OrdinalType)0,(OrdinalType)2,&root_entries[0]);   // broadcast 2 ordinals from node 0
    std::vector<OrdinalType> localChecks(2);
    localChecks[0] = negOne;  // fail or pass
    localChecks[1] = zero;    // fail reason
    if (numGlobalEntries != root_entries[0]) {
      localChecks[0] = myImageID;
      localChecks[1] = one;
    }
    else if (indexBase != root_entries[1]) {
      localChecks[0] = myImageID;
      localChecks[1] = one + one;
    }
    // if checkPassed == negOne, then it passed on this proc
    // otherwise, checkPassed == myImageID signifies it did not pass
    // REDUCE_MAX will give us the image ID of the highest rank proc that DID NOT pass, as well as the reason
    // this will be negOne and zero if all procs passed
    OrdinalType globalChecks[2];
    Teuchos::reduceAll(*comm,Teuchos::REDUCE_MAX,(OrdinalType)2,&localChecks[0],&globalChecks[0]);
    if (globalChecks[0] != negOne) {
      if (globalChecks[1] == one) {
        TEST_FOR_EXCEPTION(true,std::invalid_argument,
            errPrefix << "numGlobal must be the same on all nodes (examine node " << globalChecks[0] << ").");
      }
      else if (globalChecks[1] == one+one) {
        TEST_FOR_EXCEPTION(true,std::invalid_argument,
            errPrefix << "indexBase must be the same on all nodes (examine node " << globalChecks[0] << ").");
      }
      else {
        // logic error on my part
        TEST_FOR_EXCEPTION(true,std::invalid_argument,
            errPrefix << "logic error. Please contact the Tpetra team.");
      }
    }
    // numgGlobalEntries is coherent, but is it valid?
    TEST_FOR_EXCEPTION(numGlobalEntries < zero, std::invalid_argument,
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
    OrdinalType numLocalEntries = numGlobalEntries / numImages;    // the above code assumes truncation. is that safe? FINISH
    OrdinalType remainder = numGlobalEntries % numImages;
    OrdinalType start_index;
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
    std::vector<OrdinalType> lgMap;
    std::map<OrdinalType,OrdinalType> glMap;

    // compute the min/max global IDs
    OrdinalType minAllGID = indexBase;
    OrdinalType maxAllGID = indexBase + numGlobalEntries - one;
    OrdinalType minMyGID  = start_index + indexBase;
    OrdinalType maxMyGID  = minMyGID + numLocalEntries - one;

    Teuchos::RCP< Platform<OrdinalType> > platform_clone = platform.clone();

    // create a MapData structure 
    MapData_ = Teuchos::rcp( new MapData<OrdinalType>(indexBase,numGlobalEntries,numLocalEntries,
                                                      minAllGID,maxAllGID,minMyGID,maxMyGID,
                                                      lgMap,glMap,true,  // we are contiguous
                                                      platform_clone, comm) );

    // initialize the directory
    directorySetup();
  }

  template<typename OrdinalType>
  Map<OrdinalType>::Map (OrdinalType numGlobalEntries, OrdinalType numLocalEntries, OrdinalType indexBase, 
                         const Platform<OrdinalType> &platform)
    : Teuchos::Object("Tpetra::Map")
    , MapData_()
  {
    // Distribute the entries across the nodes so that they are 
    // - non-overlapping
    // - contiguous
    // This differs from Map(Ord,Ord,Plat) in that the user has specified the number of entries 
    // per node, so that they are not (necessarily) evenly distributed

    const OrdinalType one = Teuchos::OrdinalTraits<OrdinalType>::one();
    const OrdinalType zero = Teuchos::OrdinalTraits<OrdinalType>::zero();
    const OrdinalType negOne = zero - one;

    std::string errPrefix;
    errPrefix = "Tpetra::Map<" + Teuchos::OrdinalTraits<OrdinalType>::name() 
              + ">::constructor(numGlobal,numLocal,indexBase,platform): ";

    // get a internodal communicator from the Platform
    Teuchos::RCP< Teuchos::Comm<OrdinalType> > comm = platform.createComm();
    OrdinalType myImageID = comm->getRank();
    // for communicating failures 
    std::vector<OrdinalType> localChecks(2), globalChecks(2);

    /* compute the global size 
       we are computing the number of global entries because exactly ONE of the following is true:
       - the user didn't specify it, and we need it
       - the user did specify it, but we need to 
         + validate it against the sum of the local sizes, and
         + ensure that it is the same on all nodes
     */
    OrdinalType global_sum;
    Teuchos::reduceAll(*comm,Teuchos::REDUCE_SUM,numLocalEntries,&global_sum);
    /* there are three errors we should be detecting:
       - numGlobalEntries != -1 and it is incorrect/invalid
       - numLocalEntries invalid (<0)
    */
    localChecks[0] = negOne;
    if (numLocalEntries < zero) {
      localChecks[0] = myImageID;
      localChecks[1] = one;
    }
    else if (numGlobalEntries < negOne) {
      localChecks[0] = myImageID;
      localChecks[1] = one+one;
    }
    else if (numGlobalEntries > negOne && numGlobalEntries != global_sum) {
      localChecks[0] = myImageID;
      localChecks[1] = one+one+one;
    }
    // now check that indexBase is equivalent across images
    OrdinalType root_indexbase = indexBase;
    Teuchos::broadcast(*comm,0,&root_indexbase);   // broadcast one ordinal from node 0
    if (indexBase != root_indexbase) {
      localChecks[0] = myImageID;
      localChecks[1] = one+one+one+one;
    }
    // if checkPassed == negOne, then it passed on this proc
    // otherwise, checkPassed == myImageID signifies it did not pass
    // REDUCE_MAX will give us the image ID of the highest rank proc that DID NOT pass
    // this will be negOne if all procs passed
    Teuchos::reduceAll(*comm,Teuchos::REDUCE_MAX,(OrdinalType)2,&localChecks[0],&globalChecks[0]);
    if (globalChecks[0] != negOne) {
      if (globalChecks[1] == one) {
        TEST_FOR_EXCEPTION(true,std::invalid_argument,
            errPrefix << "numLocal is invalid (< 0) on at least one node (possibly node " 
            << globalChecks[0] << ").");
      }
      else if (globalChecks[1] == one+one) {
        TEST_FOR_EXCEPTION(true,std::invalid_argument,
            errPrefix << "numGlobal is invalid (< -1) on at least one node (possibly node " 
            << globalChecks[0] << ").");
      }
      else if (globalChecks[1] == one+one+one) {
        TEST_FOR_EXCEPTION(true,std::invalid_argument,
            errPrefix << "numGlobal doesn't match sum of numLocal (== " 
            << global_sum << ") on at least one node (possibly node " 
            << globalChecks[0] << ").");
      }
      else if (globalChecks[1] == one+one+one+one) {
        TEST_FOR_EXCEPTION(true,std::invalid_argument,
            errPrefix << "indexBase is not the same on all nodes (examine node " 
            << globalChecks[0] << ").");
      }
      else {
        // logic error on my part
        TEST_FOR_EXCEPTION(true,std::invalid_argument,
            errPrefix << "logic error. Please contact the Tpetra team.");
      }
    }
    // set numGlobalEntries
    if (numGlobalEntries == negOne) {
      numGlobalEntries = global_sum;
    }

    // compute my local offset
    OrdinalType start_index;
    Teuchos::scan(*comm,Teuchos::REDUCE_SUM,numLocalEntries,&start_index);
    start_index -= numLocalEntries;

    // create empty maps between local and global entries: let the MapData constructor fill them
    std::vector<OrdinalType> lgMap;
    std::map<OrdinalType,OrdinalType> glMap;

    // compute the min/max global IDs
    OrdinalType minAllGID = indexBase;
    OrdinalType maxAllGID = indexBase + numGlobalEntries - one;
    OrdinalType minMyGID = start_index + indexBase;
    OrdinalType maxMyGID = minMyGID + numLocalEntries - one;

    Teuchos::RCP< Platform<OrdinalType> > platform_clone = platform.clone();

    // create a MapData structure 
    MapData_ = Teuchos::rcp( new MapData<OrdinalType>(indexBase,numGlobalEntries,numLocalEntries,
                                                      minAllGID,maxAllGID,minMyGID,maxMyGID,
                                                      lgMap,glMap,true,  // we are contiguous
                                                      platform_clone, comm) );

    // initialize the directory
    directorySetup();
  }

  template<typename OrdinalType>
  Map<OrdinalType>::Map (OrdinalType numGlobalEntries, OrdinalType numLocalEntries, const std::vector<OrdinalType> &entryList, OrdinalType indexBase, 
                         const Platform<OrdinalType> &platform) 
    : Teuchos::Object("Tpetra::Map")
    , MapData_()
  {
    // Distribute the entries across the nodes in an arbitrary user-specified manner
    // They are not necessarily contiguous or evenly distributed
    const OrdinalType one = Teuchos::OrdinalTraits<OrdinalType>::one();
    const OrdinalType zero = Teuchos::OrdinalTraits<OrdinalType>::zero();
    const OrdinalType negOne = zero - one;

    std::string errPrefix;
    errPrefix = "Tpetra::Map<" + Teuchos::OrdinalTraits<OrdinalType>::name() 
              + ">::constructor(numGlobal,numLocalEntries,entryList,indexBase,platform): ";

    // platform & comm setup
    Teuchos::RCP< Teuchos::Comm<OrdinalType> > comm = platform.createComm();
    OrdinalType myImageID = comm->getRank();
    // for communicating failures 
    std::vector<OrdinalType> localChecks(2), globalChecks(2);

    /* compute the global size 
       we are computing the number of global entries because exactly ONE of the following is true:
       - the user didn't specify it, and we need it
       - the user did specify it, but we need to 
         + validate it against the sum of the local sizes, and
         + ensure that it is the same on all nodes
     */
    OrdinalType global_sum;
    Teuchos::reduceAll(*comm,Teuchos::REDUCE_SUM,numLocalEntries,&global_sum);
    localChecks[0] = negOne;
    if (numLocalEntries < zero || numLocalEntries > (OrdinalType)entryList.size()) {
      localChecks[0] = myImageID;
      localChecks[1] = one;
    }
    else if (numGlobalEntries < negOne) {
      localChecks[0] = myImageID;
      localChecks[1] = one + one;
    }
    else if (numGlobalEntries > negOne && numGlobalEntries != global_sum) {
      localChecks[0] = myImageID;
      localChecks[1] = one + one + one;
    }
    // now check that indexBase is equivalent across images
    OrdinalType root_indexbase = indexBase;
    Teuchos::broadcast(*comm,0,&root_indexbase);   // broadcast one ordinal from node 0
    if (indexBase != root_indexbase) {
      localChecks[0] = myImageID;
      localChecks[1] = one + one + one + one;
    }
    // if checkPassed == negOne, then it passed on this proc
    // otherwise, checkPassed == myImageID signifies it did not pass
    // REDUCE_MAX will give us the image ID of the highest rank proc that DID NOT pass
    // this will be negOne if all procs passed
    Teuchos::reduceAll(*comm,Teuchos::REDUCE_MAX,(OrdinalType)2,&localChecks[0],&globalChecks[0]);
    if (globalChecks[0] != negOne) {
      if (globalChecks[1] == one) {
        TEST_FOR_EXCEPTION(true,std::invalid_argument,
            errPrefix << "numLocal is invalid (< 0 or > entryList.size()) on at least one node (possibly node " 
            << globalChecks[0] << ").");
      }
      else if (globalChecks[1] == one+one) {
        TEST_FOR_EXCEPTION(true,std::invalid_argument,
            errPrefix << "numGlobal is invalid (< -1) on at least one node (possibly node "
            << globalChecks[0] << ").");
      }
      else if (globalChecks[1] == one+one+one) {
        TEST_FOR_EXCEPTION(true,std::invalid_argument,
            errPrefix << "numGlobal doesn't match sum of numLocal (" 
            << global_sum << ") on at least one node (possibly node "
            << globalChecks[0] << ").");
      }
      else if (globalChecks[1] == one+one+one+one) {
        TEST_FOR_EXCEPTION(true,std::invalid_argument,
            errPrefix << "indexBase is not the same on all nodes (possibly node "
            << globalChecks[0] << ").");
      }
      else {
        // logic error on my part
        TEST_FOR_EXCEPTION(true,std::invalid_argument,
            errPrefix << "logic error. Please contact the Tpetra team.");
      }
    }
    // set numGlobalEntries
    if (numGlobalEntries == negOne) {
      numGlobalEntries = global_sum;
    }

    // setup lgmap and glmap, and min/maxMyGIDs
    std::map<OrdinalType,OrdinalType> glMap;
    std::vector<OrdinalType> lgMap;
    OrdinalType minMyGID = indexBase;
    OrdinalType maxMyGID = indexBase;
    // create the GID to LID map; do not assume GID in entryList are distinct.
    // in the case that a GID is duplicated, keep the previous LID
    // this is necessary so that LIDs are in [0,numLocal)
    // FINISH: make sure this is legal
    {
      OrdinalType numLIDs = zero;
      if (numLocalEntries > zero) {
        for(OrdinalType i = zero; i < numLocalEntries; i++) {
          if (glMap.find(entryList[i]) == glMap.end()) {
            lgMap.push_back(entryList[i]);   // lgMap: LID to GID
            glMap[entryList[i]] = numLIDs;   // glMap: GID to LID
            numLIDs++;
          }
        }
        minMyGID = *min_element(entryList.begin(), entryList.end());
        maxMyGID = *max_element(entryList.begin(), entryList.end());
      }
      numLocalEntries = numLIDs;
    } 

    // set min/maxAllGIDs
    OrdinalType minAllGID, maxAllGID;
    {
      OrdinalType minmaxAllGIDlcl[2], minmaxAllGIDgbl[2];
      minmaxAllGIDlcl[0] = -minMyGID;  // negative allows us to do a single
      minmaxAllGIDlcl[1] =  maxMyGID;  // reduction below
      Teuchos::reduceAll<OrdinalType>(*comm,Teuchos::REDUCE_MAX,2,minmaxAllGIDlcl,minmaxAllGIDgbl);
      minAllGID = -minmaxAllGIDgbl[0];
      maxAllGID =  minmaxAllGIDgbl[1];
    }
    TEST_FOR_EXCEPTION(minAllGID < indexBase, std::invalid_argument,
        errPrefix << "minimum GID == " << minAllGID << " is less than indexBase == " << indexBase);

    Teuchos::RCP< Platform<OrdinalType> > platform_clone = platform.clone();

    // call ESData constructor
    MapData_ = Teuchos::rcp(new MapData<OrdinalType>(indexBase, numGlobalEntries, 
          numLocalEntries, minAllGID, maxAllGID, 
          minMyGID, maxMyGID, lgMap, glMap, 
          false, platform_clone, comm));

    // initialize directory
    directorySetup();
  }

  template<typename OrdinalType>
  Map<OrdinalType>::Map (const Map<OrdinalType> &map) 
    : Teuchos::Object(map.label())
    , MapData_(map.MapData_)
  {}

  template<typename OrdinalType>
  Map<OrdinalType>::~Map () 
  { }

  template<typename OrdinalType>
  OrdinalType Map<OrdinalType>::getNumGlobalEntries() const {
    return MapData_->numGlobalEntries_;
  }

  template<typename OrdinalType>
  OrdinalType Map<OrdinalType>::getNumMyEntries() const {
    return MapData_->numMyEntries_;
  }

  template<typename OrdinalType>
  OrdinalType Map<OrdinalType>::getIndexBase() const {
    return MapData_->indexBase_;
  }

  template<typename OrdinalType>
  OrdinalType Map<OrdinalType>::getMinLocalIndex() const {
    return Teuchos::OrdinalTraits<OrdinalType>::zero();
  }

  template<typename OrdinalType>
  OrdinalType Map<OrdinalType>::getMaxLocalIndex() const {
    return MapData_->numMyEntries_-1;
  }

  template<typename OrdinalType>
  OrdinalType Map<OrdinalType>::getMinGlobalIndex() const {
    return MapData_->minMyGID_;
  }

  template<typename OrdinalType>
  OrdinalType Map<OrdinalType>::getMaxGlobalIndex() const {
    return MapData_->maxMyGID_;
  }

  template<typename OrdinalType>
  OrdinalType Map<OrdinalType>::getMinAllGlobalIndex() const {
    return MapData_->minAllGID_;
  }

  template<typename OrdinalType>
  OrdinalType Map<OrdinalType>::getMaxAllGlobalIndex() const {
    return MapData_->maxAllGID_;
  }

  template<typename OrdinalType>
  OrdinalType Map<OrdinalType>::getLocalIndex(OrdinalType globalIndex) const {
    if (MapData_->contiguous_) {
      TEST_FOR_EXCEPTION(
        globalIndex < getMinGlobalIndex() || globalIndex > getMaxGlobalIndex(), 
        std::invalid_argument,
        "Tpetra::Map<" << Teuchos::OrdinalTraits<OrdinalType>::name() 
                       << ">::getLocalIndex(gid): gid does not belong to this map."
      );
      return globalIndex - getMinGlobalIndex();
    }
    else {
      typename std::map<OrdinalType,OrdinalType>::const_iterator i;
      i = MapData_->glMap_.find(globalIndex);
      TEST_FOR_EXCEPTION(
        i == MapData_->glMap_.end(), std::invalid_argument,
        "Tpetra::Map<" << Teuchos::OrdinalTraits<OrdinalType>::name() 
                       << ">::getLocalIndex(gid): gid does not belong to this map."
      );
      return i->second;
    }
  }

  template<typename OrdinalType>
  OrdinalType Map<OrdinalType>::getGlobalIndex(OrdinalType localIndex) const {
    TEST_FOR_EXCEPTION(
        localIndex < 0 || localIndex > MapData_->numMyEntries_-1,
        std::invalid_argument,
        "Tpetra::Map<" << Teuchos::OrdinalTraits<OrdinalType>::name() 
           << ">::getGlobalIndex(lid): lid invalid.");
    if (MapData_->contiguous_) {
      return localIndex + getMinGlobalIndex();
    }
    else {
      return MapData_->lgMap_[localIndex];
    }
  }

  template<typename OrdinalType>
  bool Map<OrdinalType>::isMyLocalIndex(OrdinalType localIndex) const {
    // we have local indices in [0,numLocalEntries)
    if (localIndex >= 0 && localIndex < MapData_->numMyEntries_) {
      return true;
    }
    return false;
  }

  template<typename OrdinalType>
  bool Map<OrdinalType>::isMyGlobalIndex(OrdinalType globalIndex) const {
    typename std::map<OrdinalType,OrdinalType>::iterator i;
    i = MapData_->glMap_.find(globalIndex);
    return (i != MapData_->glMap_.end());
  }

  template<typename OrdinalType>
  bool Map<OrdinalType>::isContiguous() const {
    return MapData_->contiguous_;
  }

  template<typename OrdinalType>
  bool Map<OrdinalType>::isCompatible (const Map< OrdinalType> &map) const {
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

  template<typename OrdinalType>
  bool Map<OrdinalType>::isSameAs (const Map<OrdinalType> &map) const {
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

    // check number of entries own by this node
    if (getNumMyEntries() != map.getNumMyEntries()) {
      isSame_lcl = Teuchos::ScalarTraits<char>::zero();
    }

    // check the identity of the entries owned by this node
    // only do this if we haven't already determined not-same-ness
    if (isSame_lcl == Teuchos::ScalarTraits<char>::one()) {
      // if they are contiguous, we can check the ranges easily
      // if they are not contiguous, we must check the individual LID -> GID mappings
      // the latter approach is valid in either case, but the former is faster
      if (MapData_->contiguous_ == map.MapData_->contiguous_) {
        if (MapData_->minMyGID_ != map.MapData_->minMyGID_ ||
            MapData_->maxMyGID_ != map.MapData_->maxMyGID_) {
          isSame_lcl = Teuchos::ScalarTraits<char>::zero();
        }
      }
      else {
        if (MapData_->lgMap_ != map.MapData_->lgMap_) {
          isSame_lcl = Teuchos::ScalarTraits<char>::zero();
        }
      }
    }

    // now, determine if we detected not-same-ness on any node
    char isSame_gbl;
    Teuchos::reduceAll(*MapData_->comm_,Teuchos::REDUCE_MIN,isSame_lcl,&isSame_gbl);
    return(isSame_gbl == Teuchos::ScalarTraits<char>::one());
  }

  template<typename OrdinalType>
  bool Map<OrdinalType>::operator== (const Map< OrdinalType > &map) const {
    return isSameAs(map);
  }

  template<typename OrdinalType>
  bool Map<OrdinalType>::operator!= (const Map< OrdinalType > &map) const {
    return !isSameAs(map);
  }

  template<typename OrdinalType>
  Map<OrdinalType>& Map<OrdinalType>::operator = (const Map<OrdinalType> & Source) {
    MapData_ = Source.MapData_;
    return *this;
  }

  template<typename OrdinalType>
  const std::vector<OrdinalType> & 
  Map<OrdinalType>::getMyGlobalEntries() const {
    return MapData_->lgMap_;
  }

  template<typename OrdinalType>
  bool Map<OrdinalType>::isDistributed() const {
    return MapData_->distributed_;
  }

  template<typename OrdinalType>
  void Map<OrdinalType>::print(ostream& os) const {
    const OrdinalType zero = Teuchos::OrdinalTraits<OrdinalType>::zero();
    const OrdinalType nME = getNumMyEntries();
    
    using std::endl;

    getMyGlobalEntries(); // throw away output, we call this to make sure list is generated
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
        std::vector<OrdinalType> myEntries = getMyGlobalEntries();
        for(OrdinalType i = zero; i < nME; i++) {
          os << std::setw(16) << myImageID 
             << std::setw(16) << i
             << std::setw(16) << MapData_->lgMap_[i]
             << endl;
        }
        os << flush;
      }
      // Do a few global ops to give I/O a chance to complete
      MapData_->comm_->barrier();
      MapData_->comm_->barrier();
      MapData_->comm_->barrier();
    }
  }

  template<typename OrdinalType>
  void Map<OrdinalType>::directorySetup() {
    if (getNumGlobalEntries() != Teuchos::OrdinalTraits<OrdinalType>::zero()) {
      if (MapData_->directory_ == Teuchos::null) {
        MapData_->directory_ = Teuchos::rcp( new Directory<OrdinalType>(*this) );
      }
    }
  }

  template<typename OrdinalType>
  void Map<OrdinalType>::getRemoteIndexList(
      const std::vector<OrdinalType> & GIDList, 
            std::vector<OrdinalType> & imageIDList, 
            std::vector<OrdinalType> & LIDList) const 
  {
    if (GIDList.size() == 0) return;
    TEST_FOR_EXCEPTION(getNumGlobalEntries() == 0, std::runtime_error,
        "Tpetra::Map<" + Teuchos::OrdinalTraits<OrdinalType>::name() 
        + ">::getRemoteIndexList(): getRemoteIndexList() cannot be called.");
    MapData_->directory_->getDirectoryEntries(GIDList, imageIDList, LIDList);
  }

  template<typename OrdinalType>
  void Map<OrdinalType>::getRemoteIndexList(
      std::vector<OrdinalType> const& GIDList, 
      std::vector<OrdinalType>& imageIDList) const 
  {
    if (GIDList.size() == 0) return;
    TEST_FOR_EXCEPTION(getNumGlobalEntries() == 0, std::runtime_error,
        "Tpetra::Map<" + Teuchos::OrdinalTraits<OrdinalType>::name() 
        + ">::getRemoteIndexList(): getRemoteIndexList() cannot be called.");
    MapData_->directory_->getDirectoryEntries(GIDList, imageIDList);
  }


  template<typename OrdinalType>
  Teuchos::RCP<const Platform<OrdinalType> >
  Map<OrdinalType>::getPlatform() const {
    return MapData_->platform_;
  }

  template<typename OrdinalType>
  Teuchos::RCP< Teuchos::Comm<OrdinalType> >
  Map<OrdinalType>::getComm() const {
    return MapData_->comm_;
  }


} // Tpetra namespace

#endif // TPETRA_MAP_HPP

