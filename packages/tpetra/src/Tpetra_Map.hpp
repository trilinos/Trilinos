// @HEADER
// ***********************************************************************
// 
//          Tpetra: Templated Linear Algebra Services Package
//                 Copyright (2008) Sandia Corporation
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

  template <class LocalOrdinal, class GlobalOrdinal>
  Map<LocalOrdinal,GlobalOrdinal>::Map(
                        global_size_t numGlobalElements_in, 
                        GlobalOrdinal indexBase_in, 
                        const Teuchos::RCP<const Teuchos::Comm<int> > &comm_in, 
                        LocalGlobal lOrG) : 
    comm_(comm_in) {
    // distribute the elements across the nodes so that they are 
    // - non-overlapping
    // - contiguous
    // - as evenly distributed as possible
    using Teuchos::as;
    const global_size_t GST0 = Teuchos::OrdinalTraits<global_size_t>::zero();
    const global_size_t GST1 = Teuchos::OrdinalTraits<global_size_t>::one();
    const GlobalOrdinal G1 = Teuchos::OrdinalTraits<GlobalOrdinal>::one();

    std::string errPrefix;
    errPrefix = Teuchos::typeName(*this) + "::constructor(numGlobal,indexBase,comm,lOrG): ";

    if (lOrG == GloballyDistributed) {
      const int numImages = comm_->getSize();
      const int myImageID = comm_->getRank();

      // check that numGlobalElements,indexBase is equivalent across images
      global_size_t rootNGE = numGlobalElements_in;
      GlobalOrdinal rootIB  = indexBase_in;
      Teuchos::broadcast<int,global_size_t>(*comm_,0,&rootNGE);
      Teuchos::broadcast<int,GlobalOrdinal>(*comm_,0,&rootIB);
      int localChecks[2], globalChecks[2];
      localChecks[0] = -1;   // fail or pass
      localChecks[1] = 0;    // fail reason
      if (numGlobalElements_in != rootNGE) {
        localChecks[0] = myImageID;
        localChecks[1] = 1;
      }
      else if (indexBase_in != rootIB) {
        localChecks[0] = myImageID;
        localChecks[1] = 2;
      }
      // REDUCE_MAX will give us the image ID of the highest rank proc that DID NOT pass, as well as the reason
      // these will be -1 and 0 if all procs passed
      Teuchos::reduceAll<int,int>(*comm_,Teuchos::REDUCE_MAX,2,localChecks,globalChecks);
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
          // logic error on our part
          TEST_FOR_EXCEPTION(true,std::logic_error,
              errPrefix << "logic error. Please contact the Tpetra team.");
        }
      }
      // numGlobalElements is coherent, but is it valid? this comparison looks funny, but it avoids compiler warnings on unsigned types.
      TEST_FOR_EXCEPTION(numGlobalElements_in < GST1 && numGlobalElements_in != GST0, std::invalid_argument,
          errPrefix << "numGlobalElements (== " << rootNGE << ") must be >= 0.");

      indexBase_ = rootIB;
      numGlobalElements_ = rootNGE;

      /* compute numLocalElements
         We can write numGlobalElements as follows:
         numGlobalElements == numImages * B + remainder
         Each image is allocated elements as follows:
         [ B+1    iff myImageID <  remainder
         numLocalElements = [
         [ B      iff myImageID >= remainder
         In the case that remainder == 0, then all images fall into the 
         latter case: numLocalElements == B == numGlobalElements / numImages
         It can then be shown that 
         numImages
         \Sum      numLocalElements_i  == numGlobalElements
         i=0
         This strategy is simple, requires no communication, and is optimal vis-a-vis
         uniform distribution of elements.
         This strategy is valid for any value of numGlobalElements and numImages, 
         including the following border cases:
         - numImages == 1         -> remainder == 0 && numGlobalElements == numLocalElements
         - numelements < numImages -> remainder == numGlobalElements && numLocalElements \in [0,1]
       */
      numLocalElements_ = as<size_t>(numGlobalElements_ / as<global_size_t>(numImages));
      int remainder = as<int>(numGlobalElements_ % as<global_size_t>(numImages));
#ifdef HAVE_TEUCHOS_DEBUG
      // the above code assumes truncation. is that safe?
      SHARED_TEST_FOR_EXCEPTION(numLocalElements_ * numImages + remainder != numGlobalElements_,
          std::logic_error, "Tpetra::Map::constructor(numGlobal,indexBase,platform): GlobalOrdinal does not implement division with truncation."
          << " Please contact Tpetra team.",*comm_);
#endif
      GlobalOrdinal start_index;
      if (myImageID < remainder) {
        ++numLocalElements_;
        /* the myImageID images before were each allocated 
           numGlobalElements/numImages+1
           ergo, my offset is:
           myImageID * (numGlobalElements/numImages+1)
           == myImageID * numLocalElements
         */
        start_index = as<GlobalOrdinal>(myImageID) * as<GlobalOrdinal>(numLocalElements_);
      }
      else {
        /* a quantity (equal to remainder) of the images before me
           were each allocated 
           numGlobalElements/numImages+1
           elements. a quantity (equal to myImageID-remainder) of the remaining 
           images before me were each allocated 
           numGlobalElements/numImages
           elements. ergo, my offset is:
           remainder*(numGlobalElements/numImages+1) + (myImageID-remainder)*numGlobalElements/numImages
           == remainder*numLocalElements + remainder + myImageID*numLocalElements - remainder*numLocalElements
           == myImageID*numLocalElements + remainder
         */
        start_index = as<GlobalOrdinal>(myImageID)*as<GlobalOrdinal>(numLocalElements_) + as<GlobalOrdinal>(remainder);
      }

      // compute the min/max global IDs
      minMyGID_  = start_index + indexBase_;
      maxMyGID_  = minMyGID_ + numLocalElements_ - G1;
      minAllGID_ = indexBase_;
      maxAllGID_ = indexBase_ + numGlobalElements_ - G1;
      contiguous_ = true;
      distributed_ = (numImages > 1 ? true : false);
      setupDirectory();
    }
    else {  // lOrG == LocallyReplicated
      // compute the min/max global IDs
      indexBase_ = indexBase_;
      numGlobalElements_ = numGlobalElements_in;
      numLocalElements_  = as<size_t>(numGlobalElements_in);
      minAllGID_ = indexBase_;
      maxAllGID_ = indexBase_ + numGlobalElements_ - G1;
      minMyGID_  = minAllGID_;
      maxMyGID_  = maxAllGID_;
      contiguous_ = true;
      distributed_ = false;
    }
  }

  template <class LocalOrdinal, class GlobalOrdinal>
  Map<LocalOrdinal,GlobalOrdinal>::Map(global_size_t numGlobalElements_in, size_t numLocalElements_in, GlobalOrdinal indexBase_in, 
                                       const Teuchos::RCP<const Teuchos::Comm<int> > &comm_in) :
    comm_(comm_in) {
    // Distribute the elements across the nodes so that they are 
    // - non-overlapping
    // - contiguous
    // This differs from Map(Ord,Ord,Plat) in that the user has specified the number of elements 
    // per node, so that they are not (necessarily) evenly distributed

    const size_t  L0 = Teuchos::OrdinalTraits<size_t>::zero();
    const size_t  L1 = Teuchos::OrdinalTraits<size_t>::one();
    const global_size_t GST0 = Teuchos::OrdinalTraits<global_size_t>::zero();
    const global_size_t GST1 = Teuchos::OrdinalTraits<global_size_t>::one();
    const global_size_t GSTI = Teuchos::OrdinalTraits<global_size_t>::invalid();
    const GlobalOrdinal G1 = Teuchos::OrdinalTraits<GlobalOrdinal>::one();

    std::string errPrefix;
    errPrefix = Teuchos::typeName(*this) + "::constructor(numGlobal,numLocal,indexBase,platform): ";

    // get a internodal communicator from the Platform
    const int myImageID = comm_->getRank();

    { // begin scoping block
      // for communicating failures 
      int localChecks[2], globalChecks[2];
      /* compute the global size 
         we are computing the number of global elements because exactly ONE of the following is true:
         - the user didn't specify it, and we need it
         - the user did specify it, but we need to 
           + validate it against the sum of the local sizes, and
           + ensure that it is the same on all nodes
       */
      global_size_t global_sum;
      Teuchos::reduceAll<int,global_size_t>(*comm_,Teuchos::REDUCE_SUM,Teuchos::as<global_size_t>(numLocalElements_in),&global_sum);
      /* there are three errors we should be detecting:
         - numGlobalElements != invalid() and it is incorrect/invalid
         - numLocalElements invalid (<0)
      */
      localChecks[0] = -1;
      localChecks[1] = 0;
      if (numLocalElements_in < L1 && numLocalElements_in != L0) {
        // invalid
        localChecks[0] = myImageID;
        localChecks[1] = 1;
      }
      else if (numGlobalElements_in < GST1 && numGlobalElements_in != GST0 && numGlobalElements_in != GSTI) {
        // invalid
        localChecks[0] = myImageID;
        localChecks[1] = 2;
      }
      else if (numGlobalElements_in != GSTI && numGlobalElements_in != global_sum) {
        // incorrect
        localChecks[0] = myImageID;
        localChecks[1] = 3;
      }
      // now check that indexBase is equivalent across images
      GlobalOrdinal rootIB = indexBase_in;
      Teuchos::broadcast<int,GlobalOrdinal>(*comm_,0,&rootIB);   // broadcast one ordinal from node 0
      if (indexBase_in != rootIB) {
        localChecks[0] = myImageID;
        localChecks[1] = 4;
      }
      // REDUCE_MAX will give us the image ID of the highest rank proc that DID NOT pass
      // this will be -1 if all procs passed
      Teuchos::reduceAll<int,int>(*comm_,Teuchos::REDUCE_MAX,2,localChecks,globalChecks);
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
      // set numGlobalElements
      if (numGlobalElements_in == GSTI) {
        numGlobalElements_ = global_sum;
      }
      else {
        numGlobalElements_ = numGlobalElements_in;
      }
      numLocalElements_ = numLocalElements_in;
      indexBase_ = indexBase_in;
    } // end of scoping block

    // compute my local offset
    GlobalOrdinal start_index;
    Teuchos::scan<int,GlobalOrdinal>(*comm_,Teuchos::REDUCE_SUM,numLocalElements_,&start_index);
    start_index -= numLocalElements_;

    minAllGID_ = indexBase_;
    maxAllGID_ = indexBase_ + numGlobalElements_ - G1;
    minMyGID_ = start_index + indexBase_;
    maxMyGID_ = minMyGID_ + numLocalElements_ - G1;
    contiguous_ = true;
    distributed_ = checkIsDist();
    setupDirectory();
  }

  template <class LocalOrdinal, class GlobalOrdinal>
  Map<LocalOrdinal,GlobalOrdinal>::Map (global_size_t numGlobalElements_in, const Teuchos::ArrayView<const GlobalOrdinal> &entryList, 
                                        GlobalOrdinal indexBase_in, const Teuchos::RCP<const Teuchos::Comm<int> > &comm_in) :
    comm_(comm_in) {
    using Teuchos::as;
    // Distribute the elements across the nodes in an arbitrary user-specified manner
    // They are not necessarily contiguous or evenly distributed
    const size_t  L0 = Teuchos::OrdinalTraits<size_t>::zero();
    const global_size_t GST0 = Teuchos::OrdinalTraits<global_size_t>::zero();
    const global_size_t GST1 = Teuchos::OrdinalTraits<global_size_t>::one();
    const global_size_t GSTI = Teuchos::OrdinalTraits<global_size_t>::invalid();

    LocalOrdinal numLocalElements_in = Teuchos::as<LocalOrdinal>(entryList.size());

    std::string errPrefix;
    errPrefix = Teuchos::typeName(*this) + "::constructor(numGlobal,entryList,indexBase,platform): ";

    const int myImageID = comm_->getRank();
    { // begin scoping block
      // for communicating failures 
      int localChecks[2], globalChecks[2];

      /* compute the global size 
         we are computing the number of global elements because exactly ONE of the following is true:
         - the user didn't specify it, and we need it
         - the user did specify it, but we need to 
           + validate it against the sum of the local sizes, and
           + ensure that it is the same on all nodes
       */
      global_size_t global_sum;
      Teuchos::reduceAll<int,global_size_t>(*comm_,Teuchos::REDUCE_SUM,as<global_size_t>(numLocalElements_in),&global_sum);
      localChecks[0] = -1;
      localChecks[1] = 0;
      if (numGlobalElements_in < GST1 && numGlobalElements_in != GST0 && numGlobalElements_in != GSTI) {
        // invalid
        localChecks[0] = myImageID;
        localChecks[1] = 1;
      }
      else if (numGlobalElements_in != GSTI && numGlobalElements_in != global_sum) {
        // incorrect
        localChecks[0] = myImageID;
        localChecks[1] = 2;
      }
      // now check that indexBase is equivalent across images
      GlobalOrdinal rootIB = indexBase_in;
      Teuchos::broadcast<int,GlobalOrdinal>(*comm_,0,&rootIB);   // broadcast one ordinal from node 0
      if (indexBase_in != rootIB) {
        localChecks[0] = myImageID;
        localChecks[1] = 3;
      }
      // REDUCE_MAX will give us the image ID of the highest rank proc that DID NOT pass
      // this will be -1 if all procs passed
      Teuchos::reduceAll<int,int>(*comm_,Teuchos::REDUCE_MAX,2,localChecks,globalChecks);
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

      // these are all validated/computed now
      if (numGlobalElements_in == GSTI) {
        numGlobalElements_ = global_sum;
      }
      else {
        numGlobalElements_ = numGlobalElements_in;
      }
      numLocalElements_ = numLocalElements_in;
      indexBase_ = indexBase_in;
    } // end scoping block

    // assume for now that there are numLocalElements (there may be less, if some
    // GIDs are duplicated in entryList)
    minMyGID_ = indexBase_;
    maxMyGID_ = indexBase_;
    // create the GID to LID map; do not assume GID in entryList are distinct.
    // in the case that a GID is duplicated, keep the previous LID
    // this is necessary so that LIDs are in [0,numLocal)
    size_t numUniqueGIDs = 0;
    if (numLocalElements_ > L0) {
      lgMap_ = Teuchos::arcp<GlobalOrdinal>(numLocalElements_);
      for (size_t i=0; i < numLocalElements_; i++) {
        if (glMap_.find(entryList[i]) == glMap_.end()) {
          lgMap_[numUniqueGIDs] = entryList[i];   // lgMap_:  LID to GID
          glMap_[entryList[i]] = numUniqueGIDs;   // glMap_: GID to LID
          numUniqueGIDs++;
        }
      }
      // shrink lgMap appropriately
      if (numLocalElements_ != numUniqueGIDs) {
        numLocalElements_ = numUniqueGIDs;
        lgMap_ = lgMap_.persistingView(0,numLocalElements_);
      }
      minMyGID_ = *std::min_element(lgMap_.begin(), lgMap_.end());
      maxMyGID_ = *std::max_element(lgMap_.begin(), lgMap_.end());
    }

    // set min/maxAllGIDs
    {
      GlobalOrdinal minmaxAllGIDlcl[2], minmaxAllGIDgbl[2];
      minmaxAllGIDlcl[0] = -minMyGID_;  // negative allows us to do a single
      minmaxAllGIDlcl[1] =  maxMyGID_;  // reduction below
      Teuchos::reduceAll<int,GlobalOrdinal>(*comm_,Teuchos::REDUCE_MAX,2,minmaxAllGIDlcl,minmaxAllGIDgbl);
      minAllGID_ = -minmaxAllGIDgbl[0];
      maxAllGID_ =  minmaxAllGIDgbl[1];
    }
    contiguous_  = false;
    distributed_ = checkIsDist();
    TEST_FOR_EXCEPTION(minAllGID_ < indexBase_, std::invalid_argument,
        errPrefix << "minimum GID (== " << minAllGID_ << ") is less than indexBase (== " << indexBase_ << ")");
    setupDirectory();
  }

  template <class LocalOrdinal, class GlobalOrdinal>
  Map<LocalOrdinal,GlobalOrdinal>::~Map () 
  {}

  template <class LocalOrdinal, class GlobalOrdinal>
  global_size_t Map<LocalOrdinal,GlobalOrdinal>::getGlobalNumElements() const {
    return numGlobalElements_;
  }

  template <class LocalOrdinal, class GlobalOrdinal>
  size_t Map<LocalOrdinal,GlobalOrdinal>::getLocalNumElements() const {
    return numLocalElements_;
  }

  template <class LocalOrdinal, class GlobalOrdinal>
  GlobalOrdinal Map<LocalOrdinal,GlobalOrdinal>::getIndexBase() const {
    return indexBase_;
  }

  template <class LocalOrdinal, class GlobalOrdinal>
  LocalOrdinal Map<LocalOrdinal,GlobalOrdinal>::getMinLocalIndex() const {
    return Teuchos::OrdinalTraits<LocalOrdinal>::zero();
  }

  template <class LocalOrdinal, class GlobalOrdinal>
  LocalOrdinal Map<LocalOrdinal,GlobalOrdinal>::getMaxLocalIndex() const {
    return Teuchos::as<LocalOrdinal>(numLocalElements_-1);
  }

  template <class LocalOrdinal, class GlobalOrdinal>
  GlobalOrdinal Map<LocalOrdinal,GlobalOrdinal>::getMinGlobalIndex() const {
    return minMyGID_;
  }

  template <class LocalOrdinal, class GlobalOrdinal>
  GlobalOrdinal Map<LocalOrdinal,GlobalOrdinal>::getMaxGlobalIndex() const {
    return maxMyGID_;
  }

  template <class LocalOrdinal, class GlobalOrdinal>
  GlobalOrdinal Map<LocalOrdinal,GlobalOrdinal>::getMinAllGlobalIndex() const {
    return minAllGID_;
  }

  template <class LocalOrdinal, class GlobalOrdinal>
  GlobalOrdinal Map<LocalOrdinal,GlobalOrdinal>::getMaxAllGlobalIndex() const {
    return maxAllGID_;
  }

  template <class LocalOrdinal, class GlobalOrdinal>
  LocalOrdinal Map<LocalOrdinal,GlobalOrdinal>::getLocalIndex(GlobalOrdinal globalIndex) const {
    if (contiguous_) {
      if (globalIndex < getMinGlobalIndex() || globalIndex > getMaxGlobalIndex()) {
        return Teuchos::OrdinalTraits<LocalOrdinal>::invalid();
      }
      return Teuchos::as<LocalOrdinal>(globalIndex - getMinGlobalIndex());
    }
    else {
      typename std::map<GlobalOrdinal,LocalOrdinal>::const_iterator i;
      i = glMap_.find(globalIndex);
      if (i == glMap_.end()) {
        return Teuchos::OrdinalTraits<LocalOrdinal>::invalid();
      }
      return i->second;
    }
  }

  template <class LocalOrdinal, class GlobalOrdinal>
  GlobalOrdinal Map<LocalOrdinal,GlobalOrdinal>::getGlobalIndex(LocalOrdinal localIndex) const {
    if (localIndex < getMinLocalIndex() || localIndex > getMaxLocalIndex()) {
      return Teuchos::OrdinalTraits<GlobalOrdinal>::invalid();
    }
    if (contiguous_) {
      return getMinGlobalIndex() + localIndex;
    }
    else {
      return lgMap_[localIndex];
    }
  }

  template <class LocalOrdinal, class GlobalOrdinal>
  bool Map<LocalOrdinal,GlobalOrdinal>::isMyLocalIndex(LocalOrdinal localIndex) const {
    if (localIndex < getMinLocalIndex() || localIndex > getMaxLocalIndex()) {
      return false;
    }
    return true;
  }

  template <class LocalOrdinal, class GlobalOrdinal>
  bool Map<LocalOrdinal,GlobalOrdinal>::isMyGlobalIndex(GlobalOrdinal globalIndex) const {
    if (contiguous_) {
      return (getMinGlobalIndex() <= globalIndex) && (globalIndex <= getMaxGlobalIndex());
    }
    else {
      typename std::map<GlobalOrdinal,LocalOrdinal>::const_iterator i;
      i = glMap_.find(globalIndex);
      return (i != glMap_.end());
    }
  }

  template <class LocalOrdinal, class GlobalOrdinal>
  bool Map<LocalOrdinal,GlobalOrdinal>::isContiguous() const {
    return contiguous_;
  }

  template <class LocalOrdinal, class GlobalOrdinal>
  bool Map<LocalOrdinal,GlobalOrdinal>::isCompatible (const Map< LocalOrdinal,GlobalOrdinal> &map) const {
    // check to make sure distribution is the same
    char iscompat_lcl;
    if (getGlobalNumElements() != map.getGlobalNumElements() ||
         getLocalNumElements() != map.getLocalNumElements()) {
      // NOT compat on this node
      iscompat_lcl = 0;
    }
    else {
      // compat on this node
      iscompat_lcl = 1;
    }
    char iscompat_gbl;
    Teuchos::reduceAll<int,char>(*comm_,Teuchos::REDUCE_MIN,iscompat_lcl,&iscompat_gbl);
    return (iscompat_gbl == 1);
  }

  template <class LocalOrdinal, class GlobalOrdinal>
  bool Map<LocalOrdinal,GlobalOrdinal>::isSameAs (const Map<LocalOrdinal,GlobalOrdinal> &map) const {
    if (this == &map) {
      // we should assume that this is globally coherent
      // if they share the same underlying MapData, then they must be equivalent maps
      return true;
    }

    // check all other globally coherent properties
    // if they do not share each of these properties, then they cannot be 
    // equivalent maps
    if ( (getMinGlobalIndex()   != map.getMinGlobalIndex())   ||
         (getMaxGlobalIndex()   != map.getMaxGlobalIndex())   ||
         (getGlobalNumElements() != map.getGlobalNumElements()) ||
         (isDistributed()       != map.isDistributed())       || 
         (getIndexBase()        != map.getIndexBase())         )  {
      return false;
    }

    // If we get this far, we need to check local properties and the 
    // communicate same-ness across all nodes
    // we prefer local work over communication, ergo, we will perform all
    // comparisons and conduct a single communication
    char isSame_lcl = 1;

    // check number of entries owned by this node
    if (getLocalNumElements() != map.getLocalNumElements()) {
      isSame_lcl = 0;
    }

    // check the identity of the entries owned by this node
    // only do this if we haven't already determined not-same-ness
    if (isSame_lcl == 1) {
      // if they are contiguous, we can check the ranges easily
      // if they are not contiguous, we must check the individual LID -> GID mappings
      // the latter approach is valid in either case, but the former is faster
      if (contiguous_ == true && map.contiguous_ == true) {
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
         * contiguous (one is potentially contiguous.) Calling getElementList()
         * will create it. */
        Teuchos::ArrayView<const GlobalOrdinal> ge1, ge2;
        ge1 =     getElementList();
        ge2 = map.getElementList();
        if (!std::equal(ge1.begin(),ge1.end(),ge2.begin())) {
          isSame_lcl = 0;
        }
      }
    }

    // now, determine if we detected not-same-ness on any node
    char isSame_gbl;
    Teuchos::reduceAll<int,char>(*comm_,Teuchos::REDUCE_MIN,isSame_lcl,&isSame_gbl);
    return (isSame_gbl == 1);
  }

  template <class LocalOrdinal, class GlobalOrdinal>
  Teuchos::ArrayView<const GlobalOrdinal>
  Map<LocalOrdinal,GlobalOrdinal>::getElementList() const {
    // if so (and we have local entries), then fill it.
    if (lgMap_ == Teuchos::null && numLocalElements_ > 0) {
#ifdef HAVE_TEUCHOS_DEBUG
      // this would have been set up for a non-contiguous map
      TEST_FOR_EXCEPTION(contiguous_ != true, std::logic_error,
          "Tpetra::Map::getElementList: logic error. Please notify the Tpetra team.");
#endif
      lgMap_ = Teuchos::arcp<GlobalOrdinal>(numLocalElements_);
      Teuchos::ArrayRCP<GlobalOrdinal> lgptr = lgMap_;
      for (GlobalOrdinal gid=minMyGID_; gid <= maxMyGID_; ++gid) {
        *(lgptr++) = gid;
      }
    }
    return lgMap_();
  }

  template <class LocalOrdinal, class GlobalOrdinal>
  bool Map<LocalOrdinal,GlobalOrdinal>::isDistributed() const {
    return distributed_;
  }

  template <class LocalOrdinal, class GlobalOrdinal>
  std::string Map<LocalOrdinal,GlobalOrdinal>::description() const {
    std::ostringstream oss;
    oss << Teuchos::Describable::description();
    oss << "{getGlobalNumElements() = " << getGlobalNumElements()
        << ", getLocalNumElements() = " << getLocalNumElements()
        << ", isContiguous() = " << isContiguous()
        << ", isDistributed() = " << isDistributed()
        << "}";
    return oss.str();
  }

  template <class LocalOrdinal, class GlobalOrdinal>
  void Map<LocalOrdinal,GlobalOrdinal>::describe( Teuchos::FancyOStream &out, const Teuchos::EVerbosityLevel verbLevel) const {
    using std::endl;
    using std::setw;
    using Teuchos::VERB_DEFAULT;
    using Teuchos::VERB_NONE;
    using Teuchos::VERB_LOW;
    using Teuchos::VERB_MEDIUM;
    using Teuchos::VERB_HIGH;
    using Teuchos::VERB_EXTREME;

    const size_t nME = getLocalNumElements();
    Teuchos::ArrayView<const GlobalOrdinal> myEntries = getElementList();
    int myImageID = comm_->getRank();
    int numImages = comm_->getSize();

    Teuchos::EVerbosityLevel vl = verbLevel;
    if (vl == VERB_DEFAULT) vl = VERB_LOW;

    size_t width = 1;
    for (size_t dec=10; dec<getGlobalNumElements(); dec *= 10) {
      ++width;
    }
    width = std::max<size_t>(width,14);

    Teuchos::OSTab tab(out);

    if (vl == VERB_NONE) {
      // do nothing
    }
    else if (vl == VERB_LOW) {
      out << this->description() << endl;
    }
    else {  // MEDIUM, HIGH or EXTREME
      for (int imageCtr = 0; imageCtr < numImages; ++imageCtr) {
        if (myImageID == imageCtr) {
          if (myImageID == 0) { // this is the root node (only output this info once)
            out << endl 
                << "Number of Global Entries = " << getGlobalNumElements()  << endl
                << "Maximum of all GIDs      = " << getMaxAllGlobalIndex() << endl
                << "Minimum of all GIDs      = " << getMinAllGlobalIndex() << endl
                << "Index Base               = " << getIndexBase()         << endl;
          }
          out << endl;
          if (vl == VERB_HIGH || vl == VERB_EXTREME) {
            out << "Number of Local Elements   = " << nME           << endl
                << "Maximum of my GIDs         = " << getMaxGlobalIndex() << endl
                << "Minimum of my GIDs         = " << getMinGlobalIndex() << endl;
            out << endl;
          }
          if (vl == VERB_EXTREME) {
            out << std::setw(width) << "ImageID"
                << std::setw(width) << "Local Index"
                << std::setw(width) << "Global Index"
                << endl;
            for (size_t i=0; i < nME; i++) {
              out << std::setw(width) << myImageID 
                  << std::setw(width) << i
                  << std::setw(width) << myEntries[i]
                  << endl;
            }
            out << std::flush;
          }
        }
        // Do a few global ops to give I/O a chance to complete
        comm_->barrier();
        comm_->barrier();
        comm_->barrier();
      }
    }
  }

  template <class LocalOrdinal, class GlobalOrdinal>
  void Map<LocalOrdinal,GlobalOrdinal>::setupDirectory() {
    if (getGlobalNumElements() != Teuchos::OrdinalTraits<global_size_t>::zero()) {
      if (directory_ == Teuchos::null) {
        directory_ = Teuchos::rcp( new Directory<LocalOrdinal,GlobalOrdinal>(Teuchos::rcp(this,false)) );
      }
    }
  }

  template <class LocalOrdinal, class GlobalOrdinal>
  LookupStatus Map<LocalOrdinal,GlobalOrdinal>::getRemoteIndexList(
                    const Teuchos::ArrayView<const GlobalOrdinal> & GIDList, 
                    const Teuchos::ArrayView<int> & imageIDList, 
                    const Teuchos::ArrayView<LocalOrdinal> & LIDList) const {
    if (GIDList.size() == 0) return AllIDsPresent;
    TEST_FOR_EXCEPTION(getGlobalNumElements() == 0, std::runtime_error,
        Teuchos::typeName(*this) << "::getRemoteIndexList(): getRemoteIndexList() cannot be called, zero entries in Map.");
    return directory_->getDirectoryEntries(GIDList, imageIDList, LIDList);
  }

  template <class LocalOrdinal, class GlobalOrdinal>
  LookupStatus Map<LocalOrdinal,GlobalOrdinal>::getRemoteIndexList(
                    const Teuchos::ArrayView<const GlobalOrdinal> & GIDList, 
                    const Teuchos::ArrayView<int> & imageIDList) const {
    if (GIDList.size() == 0) return AllIDsPresent;
    TEST_FOR_EXCEPTION(getGlobalNumElements() == 0, std::runtime_error,
        Teuchos::typeName(*this) << "::getRemoteIndexList(): getRemoteIndexList() cannot be called, zero entries in Map.");
    return directory_->getDirectoryEntries(GIDList, imageIDList);
  }

  template <class LocalOrdinal, class GlobalOrdinal>
  Teuchos::RCP<const Teuchos::Comm<int> >
  Map<LocalOrdinal,GlobalOrdinal>::getComm() const {
    return comm_;
  }

  template <class LocalOrdinal,class GlobalOrdinal>
  bool Map<LocalOrdinal,GlobalOrdinal>::checkIsDist() const {
    bool global = false;
    if(comm_->getSize() > 1) {
      char localRep = 0;
      if (numGlobalElements_ == Teuchos::as<global_size_t>(numLocalElements_)) {
        localRep = 1;
      }
      char allLocalRep;
      Teuchos::reduceAll<int>(*comm_,Teuchos::REDUCE_MIN,localRep,&allLocalRep);
      if (allLocalRep != 1) {
        global = true;
      }
    }
    return global;
  }

  //! Returns true if \c map is identical to this Map. Implemented in isSameAs().
  template <class LocalOrdinal, class GlobalOrdinal>
  bool operator== (const Map<LocalOrdinal,GlobalOrdinal> &map1, const Map<LocalOrdinal,GlobalOrdinal> &map2) {
    return map1.isSameAs(map2);
  }

  //! Returns true if \c map is not identical to this Map. Implemented in isSameAs().
  template <class LocalOrdinal, class GlobalOrdinal>
  bool operator!= (const Map<LocalOrdinal,GlobalOrdinal> &map1, const Map<LocalOrdinal,GlobalOrdinal> &map2) {
    return !map1.isSameAs(map2);
  }

} // Tpetra namespace

#endif // TPETRA_MAP_HPP

