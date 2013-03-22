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
//                    Jeremie Gaidamour (jngaida@sandia.gov)
//                    Jonathan Hu       (jhu@sandia.gov)
//                    Ray Tuminaro      (rstumin@sandia.gov)
//
// ***********************************************************************
//
// @HEADER

// WARNING: This code is experimental. Backwards compatibility should not be expected.

#ifndef XPETRA_STRIDEDTPETRAMAP_HPP
#define XPETRA_STRIDEDTPETRAMAP_HPP

#include "Xpetra_TpetraConfigDefs.hpp"

#include <Tpetra_Map_decl.hpp>

#include "Xpetra_StridedMap.hpp"
#include "Xpetra_TpetraMap.hpp"
//#include "Xpetra_Utils.hpp"


namespace Xpetra {

  template <class LocalOrdinal, class GlobalOrdinal = LocalOrdinal, class Node = Kokkos::DefaultNode::DefaultNodeType>
  class StridedTpetraMap
    : public virtual TpetraMap<LocalOrdinal,GlobalOrdinal,Node>, public virtual StridedMap<LocalOrdinal,GlobalOrdinal,Node> {

  public:

    //! @name Constructor/Destructor Methods
    //@{

     /** \brief Map constructor with Tpetra-defined contiguous uniform distribution.
     *
     *  Map constructor with Tpetra-defined contiguous uniform distribution.
     *  The elements are distributed among nodes so that the subsets of global
     *  elements are non-overlapping and contiguous and as evenly distributed
     *  across the nodes as possible.
     *
     *  If numGlobalElements ==
     *  Teuchos::OrdinalTraits<global_size_t>::invalid(), the number
     *  of global elements will be computed via a global
     *  communication.  Otherwise, it must be equal to the sum of the
     *  local elements across all nodes. This will only be verified if
     *  Trilinos' Teuchos package was built with debug support (CMake
     *  Boolean option TEUCHOS_ENABLE_DEBUG=ON).  If verification
     *  fails, a std::invalid_argument exception will be thrown.
     *
     *  \pre stridingInfo.size() > 0
     *  \pre numGlobalElements % getFixedBlockSize() == 0
     */
    StridedTpetraMap(global_size_t numGlobalElements, GlobalOrdinal indexBase, std::vector<size_t>& stridingInfo, const Teuchos::RCP< const Teuchos::Comm< int > > &comm, LocalOrdinal stridedBlockId=-1, GlobalOrdinal offset = 0,  LocalGlobal lg=GloballyDistributed, const Teuchos::RCP< Node > &node=Kokkos::DefaultNode::getDefaultNode())
    : Xpetra::TpetraMap<LocalOrdinal,GlobalOrdinal,Node>(numGlobalElements, indexBase, comm, lg, node), Xpetra::StridedMap<LocalOrdinal,GlobalOrdinal,Node>(numGlobalElements, indexBase, stridingInfo, comm, stridedBlockId, offset)
    {
      typedef Xpetra::StridedMap<LocalOrdinal,GlobalOrdinal,Node> StridedMapClass;
      // check input data and reorganize map

      global_size_t numGlobalNodes = Teuchos::OrdinalTraits<global_size_t>::invalid();
      if(numGlobalElements != Teuchos::OrdinalTraits<global_size_t>::invalid())
        numGlobalNodes = numGlobalElements / StridedMapClass::getFixedBlockSize();	// number of nodes (over all processors)

      // build an equally distributed node map
      RCP<Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> > nodeMap = Teuchos::rcp(new Tpetra::Map< LocalOrdinal, GlobalOrdinal, Node >(numGlobalNodes, indexBase, comm, toTpetra(lg), node));

      // translate local node ids to local dofs
      size_t nStridedOffset = 0;
      size_t nDofsPerNode = StridedMapClass::getFixedBlockSize(); // dofs per node for local striding block
      if(stridedBlockId > -1) {
        // determine nStridedOffset
        for(int j=0; j<stridedBlockId; j++) {
          nStridedOffset += stridingInfo[j];
        }
        nDofsPerNode = stridingInfo[stridedBlockId];

        numGlobalElements = nodeMap->getGlobalNumElements()*Teuchos::as<global_size_t>(nDofsPerNode);
      }
      std::vector<GlobalOrdinal> dofgids;
      for(LocalOrdinal i = 0; i<Teuchos::as<LocalOrdinal>(nodeMap->getNodeNumElements()); i++) {
        GlobalOrdinal gid = nodeMap->getGlobalElement(i);
        for(size_t dof = 0; dof < nDofsPerNode; ++dof) {
          // dofs are calculated by
          // global offset + node_GID * full size of strided map + striding offset of current striding block + dof id of current striding block
          dofgids.push_back(StridedMapClass::offset_ + gid*Teuchos::as<GlobalOrdinal>(StridedMapClass::getFixedBlockSize()) + Teuchos::as<GlobalOrdinal>(nStridedOffset + dof));
        }
      }

      const Teuchos::ArrayView<const GlobalOrdinal> dofgidsview = Teuchos::ArrayView<const GlobalOrdinal>(dofgids);
      TpetraMap<LocalOrdinal,GlobalOrdinal,Node>::map_ = Teuchos::rcp(new Tpetra::Map< LocalOrdinal, GlobalOrdinal, Node >(numGlobalElements, dofgidsview, indexBase, comm, node));

      //TEUCHOS_TEST_FOR_EXCEPTION(Xpetra::TpetraMap<LocalOrdinal,GlobalOrdinal,Node>::map_->getNodeNumElements() % nDofsPerNode != 0, Exceptions::RuntimeError, "StridedTpetraMap::StridedTpetraMap: wrong distribution of dofs among processors.");
      if(stridedBlockId == -1) {
        TEUCHOS_TEST_FOR_EXCEPTION(getNodeNumElements() != Teuchos::as<size_t>(nodeMap->getNodeNumElements()*nDofsPerNode), Exceptions::RuntimeError, "StridedTpetraMap::StridedTpetraMap: wrong distribution of dofs among processors.");
        TEUCHOS_TEST_FOR_EXCEPTION(getGlobalNumElements() != Teuchos::as<size_t>(nodeMap->getGlobalNumElements()*nDofsPerNode), Exceptions::RuntimeError, "StridedTpetraMap::StridedTpetraMap: wrong distribution of dofs among processors.");
      } else {
        TEUCHOS_TEST_FOR_EXCEPTION(stridingInfo.size() < Teuchos::as<size_t>(stridedBlockId), Exceptions::RuntimeError, "StridedTpetraMap::StridedTpetraMap: stridedBlockId > stridingInfo.size()");
        int nDofsInStridedBlock = stridingInfo[stridedBlockId];
        TEUCHOS_TEST_FOR_EXCEPTION(getNodeNumElements() != Teuchos::as<size_t>(nodeMap->getNodeNumElements()*nDofsInStridedBlock), Exceptions::RuntimeError, "StridedTpetraMap::StridedTpetraMap: wrong distribution of dofs among processors.");
        TEUCHOS_TEST_FOR_EXCEPTION(getGlobalNumElements() != Teuchos::as<size_t>(nodeMap->getGlobalNumElements()*nDofsInStridedBlock), Exceptions::RuntimeError, "StridedTpetraMap::StridedTpetraMap: wrong distribution of dofs among processors.");
      }

      TEUCHOS_TEST_FOR_EXCEPTION(CheckConsistency() == false, Exceptions::RuntimeError, "StridedTpetraMap::StridedTpetraMap: CheckConsistency() == false");
    }

    //! Map constructor with a user-defined contiguous distribution.
     /** \brief Map constructor with a user-defined contiguous distribution.
     *
     *  Map constructor with a user-defined contiguous distribution.
     *  The elements are distributed among nodes so that the subsets of global
     *  elements are non-overlapping and contiguous and as evenly distributed
     *  across the nodes as possible.
     *
     *  If numGlobalElements ==
     *  Teuchos::OrdinalTraits<global_size_t>::invalid(), the number
     *  of global elements will be computed via a global
     *  communication.  Otherwise, it must be equal to the sum of the
     *  local elements across all nodes. This will only be verified if
     *  Trilinos' Teuchos package was built with debug support (CMake
     *  Boolean option TEUCHOS_ENABLE_DEBUG=ON).  If verification
     *  fails, a std::invalid_argument exception will be thrown.
     *
     *  \pre stridingInfo.size() > 0
     *  \pre numGlobalElements % getFixedBlockSize() == 0
     *  \pre numLocalElements % getFixedBlockSize() == 0
     */
    StridedTpetraMap(global_size_t numGlobalElements, size_t numLocalElements, GlobalOrdinal indexBase, std::vector<size_t>& stridingInfo, const Teuchos::RCP< const Teuchos::Comm< int > > &comm, LocalOrdinal stridedBlockId=-1, GlobalOrdinal offset = 0, const Teuchos::RCP< Node > &node=Kokkos::DefaultNode::getDefaultNode())
    : Xpetra::TpetraMap<LocalOrdinal,GlobalOrdinal,Node>(numGlobalElements, numLocalElements, indexBase, comm, node), Xpetra::StridedMap<LocalOrdinal,GlobalOrdinal,Node>(numGlobalElements, numLocalElements, indexBase, stridingInfo, comm, stridedBlockId, offset)
    {
      typedef Xpetra::StridedMap<LocalOrdinal,GlobalOrdinal,Node> StridedMapClass;
      // check input data and reorganize map

      global_size_t numGlobalNodes = Teuchos::OrdinalTraits<global_size_t>::invalid();
      if(numGlobalElements != Teuchos::OrdinalTraits<global_size_t>::invalid())
        numGlobalNodes = numGlobalElements / StridedMapClass::getFixedBlockSize();	// number of nodes (over all processors)
      global_size_t numLocalNodes  = numLocalElements / StridedMapClass::getFixedBlockSize();	// number of nodes (on one processor)

      // build an equally distributed node map
      RCP<Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> > nodeMap = Teuchos::rcp(new Tpetra::Map< LocalOrdinal, GlobalOrdinal, Node >(numGlobalNodes, numLocalNodes, indexBase, comm, node));

      // translate local node ids to local dofs
      size_t nStridedOffset = 0;
      size_t nDofsPerNode = StridedMapClass::getFixedBlockSize(); // dofs per node for local striding block
      if(stridedBlockId > -1) {
        // determine nStridedOffset
        for(int j=0; j<stridedBlockId; j++) {
          nStridedOffset += stridingInfo[j];
        }
        nDofsPerNode = stridingInfo[stridedBlockId];

        numGlobalElements = nodeMap->getGlobalNumElements()*Teuchos::as<global_size_t>(nDofsPerNode);
      }
      std::vector<GlobalOrdinal> dofgids;
      for(LocalOrdinal i = 0; i<Teuchos::as<LocalOrdinal>(nodeMap->getNodeNumElements()); i++) {
        GlobalOrdinal gid = nodeMap->getGlobalElement(i);
        for(size_t dof = 0; dof < nDofsPerNode; ++dof) {
          // dofs are calculated by
          // global offset + node_GID * full size of strided map + striding offset of current striding block + dof id of current striding block
          dofgids.push_back(StridedMapClass::offset_ + gid*Teuchos::as<GlobalOrdinal>(StridedMapClass::getFixedBlockSize()) + Teuchos::as<GlobalOrdinal>(nStridedOffset + dof));
        }
      }

      const Teuchos::ArrayView<const GlobalOrdinal> dofgidsview = Teuchos::ArrayView<const GlobalOrdinal>(dofgids);
      TpetraMap<LocalOrdinal,GlobalOrdinal,Node>::map_ = Teuchos::rcp(new Tpetra::Map< LocalOrdinal, GlobalOrdinal, Node >(numGlobalElements, dofgidsview, indexBase, comm, node));

      //TEUCHOS_TEST_FOR_EXCEPTION(Xpetra::TpetraMap<LocalOrdinal,GlobalOrdinal,Node>::map_->getNodeNumElements() % nDofsPerNode != 0, Exceptions::RuntimeError, "StridedTpetraMap::StridedTpetraMap: wrong distribution of dofs among processors.");
      if(stridedBlockId == -1) {
        TEUCHOS_TEST_FOR_EXCEPTION(getNodeNumElements() != Teuchos::as<size_t>(nodeMap->getNodeNumElements()*nDofsPerNode), Exceptions::RuntimeError, "StridedTpetraMap::StridedTpetraMap: wrong distribution of dofs among processors.");
        TEUCHOS_TEST_FOR_EXCEPTION(getGlobalNumElements() != Teuchos::as<size_t>(nodeMap->getGlobalNumElements()*nDofsPerNode), Exceptions::RuntimeError, "StridedTpetraMap::StridedTpetraMap: wrong distribution of dofs among processors.");
      } else {
        TEUCHOS_TEST_FOR_EXCEPTION(stridingInfo.size() < Teuchos::as<size_t>(stridedBlockId), Exceptions::RuntimeError, "StridedTpetraMap::StridedTpetraMap: stridedBlockId > stridingInfo.size()");
        int nDofsInStridedBlock = stridingInfo[stridedBlockId];
        TEUCHOS_TEST_FOR_EXCEPTION(getNodeNumElements() != Teuchos::as<size_t>(nodeMap->getNodeNumElements()*nDofsInStridedBlock), Exceptions::RuntimeError, "StridedTpetraMap::StridedTpetraMap: wrong distribution of dofs among processors.");
        TEUCHOS_TEST_FOR_EXCEPTION(getGlobalNumElements() != Teuchos::as<size_t>(nodeMap->getGlobalNumElements()*nDofsInStridedBlock), Exceptions::RuntimeError, "StridedTpetraMap::StridedTpetraMap: wrong distribution of dofs among processors.");
      }

      TEUCHOS_TEST_FOR_EXCEPTION(CheckConsistency() == false, Exceptions::RuntimeError, "StridedTpetraMap::StridedTpetraMap: CheckConsistency() == false");

    }

    /** \brief Map constructor with user-defined non-contiguous (arbitrary) distribution.
     *
     *  createse a strided map using the GIDs in elementList and the striding information
     *  provided by user.
     *
     *  \pre stridingInfo.size() > 0
     *  \pre numGlobalElements % getFixedBlockSize() == 0
     *  \pre elementList.size() % getFixedBlockSize() == 0
     *  \post CheckConsistency() == true
     */
    StridedTpetraMap(global_size_t numGlobalElements, const Teuchos::ArrayView< const GlobalOrdinal > &elementList, GlobalOrdinal indexBase, std::vector<size_t>& stridingInfo, const Teuchos::RCP< const Teuchos::Comm< int > > &comm, LocalOrdinal stridedBlockId=-1, const Teuchos::RCP< Node > &node=Kokkos::DefaultNode::getDefaultNode())
    : Xpetra::TpetraMap<LocalOrdinal,GlobalOrdinal,Node>(numGlobalElements, elementList, indexBase, comm, node), Xpetra::StridedMap<LocalOrdinal,GlobalOrdinal,Node>(numGlobalElements, elementList, indexBase, stridingInfo, comm, stridedBlockId)
    {
      typedef Xpetra::StridedMap<LocalOrdinal,GlobalOrdinal,Node> StridedMapClass;

        if(stridedBlockId != -1) {
          TEUCHOS_TEST_FOR_EXCEPTION(stridingInfo.size() < Teuchos::as<size_t>(stridedBlockId), Exceptions::RuntimeError, "StridedTpetraMap::StridedTpetraMap: stridedBlockId > stridingInfo.size()");
        }

        // create TpetraMap using the dofs from ElementList
        TpetraMap<LocalOrdinal,GlobalOrdinal,Node>::map_ = Teuchos::rcp(new Tpetra::Map< LocalOrdinal, GlobalOrdinal, Node >(numGlobalElements, elementList, indexBase, comm, node));

        // set parameters for striding information
        TEUCHOS_TEST_FOR_EXCEPTION(CheckConsistency() == false, Exceptions::RuntimeError, "StridedTpetraMap::StridedTpetraMap: CheckConsistency() == false");
    }

    //! Map destructor.
    ~StridedTpetraMap() { }

    //@}

    //! @name Map Attribute Methods
    //@{

    //! Returns the number of elements in this Map.
    global_size_t getGlobalNumElements() const { return TpetraMap<LocalOrdinal,GlobalOrdinal,Node>::getGlobalNumElements(); }

    //! Returns the number of elements belonging to the calling node.
    size_t getNodeNumElements() const { return TpetraMap<LocalOrdinal,GlobalOrdinal,Node>::getNodeNumElements(); }

    //! Returns the index base for this Map.
    GlobalOrdinal getIndexBase() const { return TpetraMap<LocalOrdinal,GlobalOrdinal,Node>::getIndexBase(); }

    //! Returns minimum local index.
    LocalOrdinal getMinLocalIndex() const { return TpetraMap<LocalOrdinal,GlobalOrdinal,Node>::getMinLocalIndex(); }

    //! Returns maximum local index.
    LocalOrdinal getMaxLocalIndex() const { return TpetraMap<LocalOrdinal,GlobalOrdinal,Node>::getMaxLocalIndex(); }

    //! Returns minimum global index owned by this node.
    GlobalOrdinal getMinGlobalIndex() const { return TpetraMap<LocalOrdinal,GlobalOrdinal,Node>::getMinGlobalIndex(); }

    //! Returns maximum global index owned by this node.
    GlobalOrdinal getMaxGlobalIndex() const { return TpetraMap<LocalOrdinal,GlobalOrdinal,Node>::getMaxGlobalIndex(); }

    //! Return the minimum global index over all nodes.
    GlobalOrdinal getMinAllGlobalIndex() const { return TpetraMap<LocalOrdinal,GlobalOrdinal,Node>::getMinAllGlobalIndex(); }

    //! Return the maximum global index over all nodes.
    GlobalOrdinal getMaxAllGlobalIndex() const { return TpetraMap<LocalOrdinal,GlobalOrdinal,Node>::getMaxAllGlobalIndex(); }

    //! Return the local index for a given global index.
    LocalOrdinal getLocalElement(GlobalOrdinal globalIndex) const { return TpetraMap<LocalOrdinal,GlobalOrdinal,Node>::getLocalElement(globalIndex); }

    //! Return the global index for a given local index.
    GlobalOrdinal getGlobalElement(LocalOrdinal localIndex) const { return TpetraMap<LocalOrdinal,GlobalOrdinal,Node>::getGlobalElement(localIndex); }

    //! Returns the node IDs and corresponding local indices for a given list of global indices.
    LookupStatus getRemoteIndexList(const Teuchos::ArrayView< const GlobalOrdinal > &GIDList, const Teuchos::ArrayView< int > &nodeIDList, const Teuchos::ArrayView< LocalOrdinal > &LIDList) const { return TpetraMap<LocalOrdinal,GlobalOrdinal,Node>::getRemoteIndexList(GIDList, nodeIDList, LIDList); }

    //! Returns the node IDs for a given list of global indices.
    LookupStatus getRemoteIndexList(const Teuchos::ArrayView< const GlobalOrdinal > &GIDList, const Teuchos::ArrayView< int > &nodeIDList) const { return TpetraMap<LocalOrdinal,GlobalOrdinal,Node>::getRemoteIndexList(GIDList, nodeIDList); }

    //! Return a list of the global indices owned by this node.
    Teuchos::ArrayView< const GlobalOrdinal > getNodeElementList() const { return TpetraMap<LocalOrdinal,GlobalOrdinal,Node>::getNodeElementList(); }

    //! Returns true if the local index is valid for this Map on this node; returns false if it isn't.
    bool isNodeLocalElement(LocalOrdinal localIndex) const { return TpetraMap<LocalOrdinal,GlobalOrdinal,Node>::isNodeLocalElement(localIndex); }

    //! Returns true if the global index is found in this Map on this node; returns false if it isn't.
    bool isNodeGlobalElement(GlobalOrdinal globalIndex) const { return TpetraMap<LocalOrdinal,GlobalOrdinal,Node>::isNodeGlobalElement(globalIndex); }

    //! Returns true if this Map is distributed contiguously; returns false otherwise.
    bool isContiguous() const { return TpetraMap<LocalOrdinal,GlobalOrdinal,Node>::isContiguous(); }

    //! Returns true if this Map is distributed across more than one node; returns false otherwise.
    bool isDistributed() const { return TpetraMap<LocalOrdinal,GlobalOrdinal,Node>::isDistributed(); }

    //@}

    //! @name Boolean Tests
    //@{

    //! Returns true if map is compatible with this Map.
    bool isCompatible(const Map< LocalOrdinal, GlobalOrdinal, Node > &map) const { return TpetraMap<LocalOrdinal,GlobalOrdinal,Node>::isCompatible(map); }

    //! Returns true if map is identical to this Map.
    bool isSameAs(const Map< LocalOrdinal, GlobalOrdinal, Node > &map) const { return TpetraMap<LocalOrdinal,GlobalOrdinal,Node>::isSameAs(map); }

    //@}

    //! @name
    //@{

    //! Get the Comm object for this Map.
    const Teuchos::RCP< const Teuchos::Comm< int > >  getComm() const { return TpetraMap<LocalOrdinal,GlobalOrdinal,Node>::getComm(); }

    //! Get the Node object for this Map.
    const Teuchos::RCP< Node >  getNode() const { return TpetraMap<LocalOrdinal,GlobalOrdinal,Node>::getNode(); }

    //@}

    //! @name
    //@{

    //! Return a simple one-line description of this object.
    std::string description() const { return TpetraMap<LocalOrdinal,GlobalOrdinal,Node>::description(); }

    //! Print the object with some verbosity level to a FancyOStream object.
    void describe(Teuchos::FancyOStream &out, const Teuchos::EVerbosityLevel verbLevel=Teuchos::Describable::verbLevel_default) const { TpetraMap<LocalOrdinal,GlobalOrdinal,Node>::describe(out, verbLevel); }

    //@}

    //! @name Xpetra specific
    //@{

    //! TpetraMap constructor to wrap a Tpetra::Map object
    StridedTpetraMap(const Teuchos::RCP<const Tpetra::Map<LocalOrdinal, GlobalOrdinal, Node > > &map, std::vector<size_t>& stridingInfo, LocalOrdinal stridedBlockId=-1, GlobalOrdinal offset = 0)
      : Xpetra::TpetraMap<LocalOrdinal,GlobalOrdinal,Node>(map), Xpetra::StridedMap<LocalOrdinal,GlobalOrdinal,Node>(stridingInfo, stridedBlockId, offset) {
      typedef Xpetra::StridedMap<LocalOrdinal,GlobalOrdinal,Node> StridedMapClass;
      //size_t nDofsPerNode = Teuchos::as<size_t>(StridedMapClass::getFixedBlockSize());
      //TEUCHOS_TEST_FOR_EXCEPTION((Tpetra::Map<LocalOrdinal, GlobalOrdinal, Node >)::map_->getNodeNumElements() % nDofsPerNode != 0, Exceptions::RuntimeError, "StridedTpetraMap::StridedTpetraMap: wrong distribution of dofs among processors.");
      TEUCHOS_TEST_FOR_EXCEPTION(CheckConsistency() == false, Exceptions::RuntimeError, "StridedTpetraMap::StridedTpetraMap: CheckConsistency() == false");
    }

    //! Get the library used by this object (Tpetra or Epetra?)
    UnderlyingLib lib() const { return Xpetra::UseTpetra; }

    /*//! Get the underlying Tpetra map
    const RCP< const Tpetra::Map< LocalOrdinal, GlobalOrdinal, Node > > & getTpetra_Map() const { return map_; }*/

    //@}

  private:
    bool CheckConsistency() {
      typedef Xpetra::StridedMap<LocalOrdinal,GlobalOrdinal,Node> StridedMapClass;

      if(StridedMapClass::getStridedBlockId() == -1) {

        if(getNodeNumElements() % StridedMapClass::getFixedBlockSize() != 0) return false;
        if(getGlobalNumElements() % StridedMapClass::getFixedBlockSize() != 0) return false;
      }
      else {
        Teuchos::ArrayView< const GlobalOrdinal > dofGids = getNodeElementList();
        //std::sort(dofGids.begin(),dofGids.end());

        // determine nStridedOffset
        size_t nStridedOffset = 0;
        for(int j=0; j<StridedMapClass::stridedBlockId_; j++) {
          nStridedOffset += StridedMapClass::stridingInfo_[j];
        }
        //size_t nDofsPerNode = stridingInfo_[stridedBlockId_];

        const GlobalOrdinal goStridedOffset = Teuchos::as<GlobalOrdinal>(nStridedOffset);
        const GlobalOrdinal goZeroOffset = (dofGids[0] - nStridedOffset - StridedMapClass::offset_) / Teuchos::as<GlobalOrdinal>(StridedMapClass::getFixedBlockSize());

        GlobalOrdinal cnt = 0;
        for(size_t i = 0; i<Teuchos::as<size_t>(dofGids.size())/StridedMapClass::stridingInfo_[StridedMapClass::stridedBlockId_]; i+=StridedMapClass::stridingInfo_[StridedMapClass::stridedBlockId_]) {

          for(size_t j=0; j<StridedMapClass::stridingInfo_[StridedMapClass::stridedBlockId_]; j++) {
            const GlobalOrdinal gid = dofGids[i+j];
            if((gid - Teuchos::as<GlobalOrdinal>(j) - goStridedOffset - StridedMapClass::offset_) / Teuchos::as<GlobalOrdinal>(StridedMapClass::getFixedBlockSize()) - goZeroOffset - cnt != 0) {
              //std::cout << "gid: " << gid << " GID: " <<  (gid - Teuchos::as<GlobalOrdinal>(j) - goStridedOffset) / Teuchos::as<GlobalOrdinal>(getFixedBlockSize()) - goZeroOffset - cnt << std::endl;
              return false;
            }
          }
          cnt++;
        }
      }

      return true;
    }

  }; // StridedTpetraMap class


} // Xpetra namespace

#define XPETRA_STRIDEDTPETRAMAP_SHORT
#endif // XPETRA_STRIDEDTPETRAMAP_HPP
