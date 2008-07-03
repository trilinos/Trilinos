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
#include "Tpetra_MapDecl.hpp"

// FINISH: get a list of all globally coherent data members

namespace Tpetra {

  template<typename OrdinalType>
  Map<OrdinalType>::Map (OrdinalType numGlobalEntries, OrdinalType indexBase, const Platform &platform) 
    : Teuchos::Object("Tpetra::Map")
    , MapData_()
  {
    // distribute the entries across the nodes so that they are 
    // - non-overlapping
    // - contiguous
    // - as evenly distributed as possible
    const OrdinalType one = Teuchos::OrdinalTraits<OrdinalType>::one();
    const OrdinalType zero = Teuchos::OrdinalTraits<OrdinalType>::zero();

    // initial tests
    TEST_FOR_EXCEPTION(numGlobalEntries < 0, std::invalid_argument,
          "Tpetra::Map<" << Teuchos::OrdinalTraits<OrindalType>::name() << ">::constructor(numGlobal,base,platform): numGlobal == "
          << numGlobalEntries << ". Must be >= 0.");

    // get a internodal communicator from the Platform
    Teuchos::RCP< Teuchos::Comm<OrdinalType> > comm = platform.createComm();
    OrdinalType numImages = comm->getSize();
    OrdinalType myImageID = comm->getRank();

    // check that numGlobalEntries is equivalent across images
    OrdinalType root_numGlobalEntries = numGlobalEntries;
    Teuchos::broadcast(*comm,0,&root_numGlobalEntries);
    bool checkPassed = (numGlobalEntries == root_numGlobalEntries);
    bool globalCheckPassed;
    Teuchos::reduceAll(*comm,Teuchos::REDUCE_AND,checkPassed,&globalCheckPassed);
    TEST_FOR_EXCEPTION(globalCheckPassed == false,std::invalid_argument,
          "Tpetra::Map<" << Teuchos::OrdinalTraits<OrindalType>::name() << ">::constructor(numGlobal,base,platform): numGlobal == "
          << numGlobalEntries << ". Must be the same on all nodes.");

    /* compute numMyEntries
       We can write numGlobalEntries as follows:
          numGlobalEntries == numImages * B + remainder
       Each image is allocated entries as follows:
                         [ B+1    iff myImageID <  remainder
          numMyEntries = [
                         [ B      iff myImageID >= remainder
       In the case that remainder == 0, then all images fall into the 
       latter case: numMyEntries == B == numGlobalEntries / numImages
       It can then be shown that 
          numImages
            \Sum      numMyEntries_i  == numGlobalEntries
             i=0
       This strategy is simple, requires no communication, and is optimal vis-a-vis
       uniform distribution of entries.
       This strategy is valid for any value of numGlobalEntries and numImages, 
       including the following border cases:
         - numImages == 1         -> remainder == 0 && numGlobalEntries == numMyEntries
         - numEntries < numImages -> remainder == numGlobalEntries && numMyEntries \in [0,1]
     */
    OrdinalType numMyEntries = numGlobalEntries / numImages;    // the above code assumes truncation: FINISH
    OrdinalType remainder = numGlobalEntries % numImages;
    OrdinalType start_index;
    if (myImageID < remainder) {
      ++numMyEntries;
      /* the myImageID images before were each allocated 
            numGlobalEntries/numImages+1
         ergo, my offset is:
            myImageID * (numGlobalEntries/numImages+1)
            == myImageID * numMyEntries
       */
      start_index = myImageID * numMyEntries;
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
           == remainder*numMyEntries + remainder + myImageID*numMyEntries - remainder*numMyEntries
           == myImageID*numMyEntries + remainder
       */
      start_index = myImageID*numMyEntries + remainder;
    }

    // create empty maps between local and global entries: let the MapData constructor fill them
    std::map<OrdinalType,OrdinalType> lgMap, glMap;

    // compute the min/max global IDs
    OrdinalType minAllGID = indexBase;
    OrdinalType maxAllGID = indexBase + numGlobalEntries - one;
    OrdinalType minMyGID  = start_index + indexBase;
    OrdinalType maxMyGID  = minMyGID + numMyEntries - one;

    Teuchos::RCP< Platform<OrdinalType> > platform_clone = platform.clone();

    // create a MapData structure 
    MapData_ = Teuchos::rcp( new MapData<OrdinalType>(indexBase,numGlobalEntries,numMyEntries,
                                                      minAllGID,maxAllGID,minMyGID,maxMyGID,
                                                      lgMap,glMap,true,  // we are contiguous
                                                      platform_clone, comm) );

    // initialize the directory
    directorySetup();
  }

  template<typename OrdinalType>
  Map<OrdinalType>::Map (OrdinalType numGlobalEntries, OrdinalType numMyEntries, OrdinalType indexBase, 
                         const Platform &platform) 
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

    // initial tests
    TEST_FOR_EXCEPTION(numGlobalEntries < negOne, std::invalid_argument,
          "Tpetra::Map<" << Teuchos::OrdinalTraits<OrindalType>::name() << ">::constructor(numGlobal,numLocal,base,platform): numGlobal == "
          << numGlobalEntries << ". Must be >= 0 (to specify) or -1 (to compute).");
    TEST_FOR_EXCEPTION(numLocalEntries < zero, std::invalid_argument,
          "Tpetra::Map<" << Teuchos::OrdinalTraits<OrindalType>::name() << ">::constructor(numGlobal,numLocal,base,platform): numLocal == "
          << numLocalEntries << ". Must be >= 0.");

    // get a internodal communicator from the Platform
    Teuchos::RCP< Teuchos::Comm<OrdinalType> > comm = platform.createComm();

    /* compute the global size 
       we are computing the number of global entries because exactly ONE of the following is true:
       - the user didn't specify it, and we need it
       - the user did specify it, but we need to 
         + validate it against the sum of the local sizes, and
         + ensure that it is the same on all nodes
     */
    int global_sum;
    Teuchos::reduceAll(*comm,Teuchos::REDUCE_SUM,numMyEntries,&global_sum);
    // set it or check it
    if (numGlobalEntries == -1) {
      numGlobalEntries = global_sum;
    }
    else {
      // if the global entries don't agree across all nodes, then at least one node will exit here.
      TEST_FOR_EXCEPTION(numGlobalEntries != global_sum, std::invalid_argument,
          "Tpetra::Map<" << Teuchos::OrdinalTraits<OrindalType>::name() << ">::constructor(numGlobal,numLocal,base,platform): numGlobal ("
          << numGlobalEntries << ") did not equal the sum of numLocal across all nodes (" << global_sum << ").");
    }

    // compute my local offset
    OrdinalType start_index;
    Teuchos::scan(*comm,Teuchos::REDUCE_SUM,numMyEntries,&start_index);
    start_index -= numMyEntries;

    // create empty maps between local and global entries: let the MapData constructor fill them
    std::map<OrdinalType,OrdinalType> lgMap, glMap;

    // compute the min/max global IDs
    OrdinalType minAllGID = indexBase;
    OrdinalType maxAllGID = indexBase + numGlobalElements - one;
    OrdinalType minMyGID = start_index + indexBase;
    OrdinalType maxMyGID = minMyGID + numMyEntries - one;

    Teuchos::RCP< Platform<OrdinalType> > platform_clone = platform.clone();

    // create a MapData structure 
    MapData_ = Teuchos::rcp( new MapData<OrdinalType>(indexBase,numGlobalEntries,numMyEntries,
                                                      minAllGID,maxAllGID,minMyGID,maxMyGID,
                                                      lgMap,glMap,true,  // we are contiguous
                                                      platform_clone, comm) );

    // initialize the directory
    directorySetup();
  }

  template<typename OrdinalType>
  Map<OrdinalType>::Map (OrdinalType numGlobalEntries, const std::vector< OrdinalType > &entryList, OrdinalType indexBase, 
                         const Platform &platform) 
    : Teuchos::Object("Tpetra::Map")
    , MapSpaceData_()
  {
    // Distribute the entries across the nodes in an arbitrary user-specified manner
    // They are not necessarily contiguous or evenly distributed
    const OrdinalType one = Teuchos::OrdinalTraits<OrdinalType>::one();
    const OrdinalType zero = Teuchos::OrdinalTraits<OrdinalType>::zero();
    const OrdinalType negOne = zero - one;

    // initial tests
    TEST_FOR_EXCEPTION(numGlobalEntries < negOne, std::invalid_argument,
          "Tpetra::Map<" << Teuchos::OrdinalTraits<OrindalType>::name() << ">::constructor(numGlobal,localList,base,platform): numGlobal == "
          << numGlobalEntries << ". Must be >= 0 (to specify) or -1 (to compute).");

    // platform & comm setup
    Teuchos::RCP< Teuchos::Comm<OrdinalType> > comm = platform.createOrdinalComm();

    // get the number of local entries: the size of the entry list
    OrdinalType numMyEntries = (OrdinalType)entryList.size();

    /* compute the global size 
       we are computing the number of global entries because exactly ONE of the following is true:
       - the user didn't specify it, and we need it
       - the user did specify it, but we need to 
         + validate it against the sum of the local sizes, and
         + ensure that it is the same on all nodes
     */
    int global_sum;
    Teuchos::reduceAll(*comm,Teuchos::REDUCE_SUM,numMyEntries,&global_sum);
    // set it or check it
    if (numGlobalEntries == -1) {
      numGlobalEntries = global_sum;
    }
    else {
      // if the global entries don't agree across all nodes, then at least one node will exit here.
      TEST_FOR_EXCEPTION(numGlobalEntries != global_sum, std::invalid_argument,
          "Tpetra::Map<" << Teuchos::OrdinalTraits<OrindalType>::name() << ">::constructor(numGlobal,localList,base,platform): numGlobal ("
          << numGlobalEntries << ") did not equal the sum of numLocal across all nodes (" << global_sum << ").");
    }

    // FINISH: this effort should also check the validity of the local entry list, IF we require that the list is "valid"
    // setup lgmap and glmap, and min/maxMyGIDs
    map<OrdinalType, OrdinalType> lgMap;
    map<OrdinalType, OrdinalType> glMap;
    // FINISH: make sure this doesn't break anything elsewhere
    // set minMyGID and maxMyGID to empty, i.e., so that there are no numbers satisfying 
    //         minMyGID <= i <= maxMyGID
    // this can only be done if maxMyGID < minMyGID
    OrdinalType minMyGID = indexBase;
    OrdinalType maxMyGID = indexBase-1;
    // FINISH: modify this code to check for contiguity
    if (numMyElements > zero) {
      for(OrdinalType i = zero; i < numMyElements; i++) {
        lgMap[i] = elementList[i];   // lgmap: LID to GID
        glMap[elementList[i]] = i;   // glmap: GID to LID
      }
      minMyGID = *min_element(elementList.begin(), elementList.end());
      maxMyGID = *max_element(elementList.begin(), elementList.end());
    }

    // set min/maxAllGIDs
    OrdinalType minAllGID;
    OrdinalType maxAllGID;
    Teuchos::reduceAll(*comm,Teuchos::REDUCE_MIN,minMyGid,&minAllGID);
    Teuchos::reduceAll(*comm,Teuchos::REDUCE_MAX,maxMyGid,&maxAllGID);
    TEST_FOR_EXCEPTION(minAllGID < indexBase, std::invalid_argument,
        "Tpetra::Map<" << Teuchos::OrdinalTraits<OrindalType>::name() << ">::constructor(numGlobal,entryList,base,platform): minimum global element = " 
        << minAllGID << " is less than the index base = " << indexBase);

    Teuchos::RCP< Platform<OrdinalType, OrdinalType> > platform_clone = platform.clone();

    // call ESData constructor
    ElementSpaceData_ = Teuchos::rcp(new ElementSpaceData<OrdinalType>(indexBase, numGlobalElements, 
          numMyElements, minAllGID, maxAllGID, 
          minMyGID, maxMyGID, lgMap, glMap, 
          false, platform_clone, comm));

    // initialize directory
    directorySetup();
  }

  template<typename OrdinalType>
  Map<OrdinalType>::Map (const Map<OrdinalType> &map) 
    : Teuchos::Object(map.labe())
    , MapData_(map.MapData_)
  {}

  template<typename OrdinalType>
  Map<OrdinalType>::~Map () 
  { 
    // FINISH: what about the directory? i don't think it needs any special treatment
  }

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
    return MapData_->minGID_;
  }

  template<typename OrdinalType>
  OrdinalType Map<OrdinalType>::getMaxGlobalIndex() const {
    return MapData_->maxGID_;
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
      std::map<OrdinalType,OrdinalType>::iterator i;
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
    if (MapData_->contiguous_) {
      TEST_FOR_EXCEPTION(
        localIndex < 0 || localIndex > MapData_->numMyEntries_-1,
        std::invalid_argument,
        "Tpetra::Map<" << Teuchos::OrdinalTraits<OrdinalType>::name() 
                       << ">::getGlobalIndex(lid): lid invalid."
      );
      return globalIndex - getMinGlobalIndex();
    }
    else {
      std::map<OrdinalType,OrdinalType>::iterator i;
      i = MapData_->lgMap_.find(globalIndex);
      TEST_FOR_EXCEPTION(
        i == MapData_->lgMap_.end(), std::invalid_argument,
        "Tpetra::Map<" << Teuchos::OrdinalTraits<OrdinalType>::name() 
                       << ">::getGlobalIndex(lid): lid invalid."
      );
      return i->second;
    }
  }

  template<typename OrdinalType>
  bool Map<OrdinalType>::isMyLocalIndex(OrdinalType localIndex) const {
    // we have local indices in [0,numMyEntries)
    if (localIndex >= 0 && localIndex < MapData_->numMyEntries_) {
      return true;
    }
    return false;
  }

  template<typename OrdinalType>
  bool Map<OrdinalType>::isMyGlobalIndex(OrdinalType globalIndex) const {
    std::map<OrdinalType,OrdinalType>::iterator i;
    i = MapData_->glMap_.find(globalIndex);
    return (i != MapData_->glMap_.end());
  }

  template<typename OrdinalType>
  bool Map<OrdinalType>::isContiguous() const {
    return MapData_->contiguous_;
  }

  template<typename OrdinalType>
  bool Map<OrdinalType>::isCompatible (const Map< OrdinalType> &map) const {
    // FINISH: figure out what "compatible" means
    // i think it means that, on every node, the maps have the same number of global and local entries
  }

  template<typename OrdinalType>
  bool Map<OrdinalType>::isSameAs (const Map<OrdinalType> &map) const {
    // FINISH: make sure the latter is written so that two maps can be "same", even if one 
    // is marked contiguous and the other is not
    return (MapData_ == map.MapData_);
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
  void Map<OrdinalType>::print(ostream& os) const {
    // FINISH
  }

} // Tpetra namespace

#endif // TPETRA_MAP_HPP

