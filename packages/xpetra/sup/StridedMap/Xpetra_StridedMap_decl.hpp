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
//
// ***********************************************************************
//
// @HEADER

// WARNING: This code is experimental. Backwards compatibility should not be expected.

#ifndef XPETRA_STRIDEDMAP_DECL_HPP
#define XPETRA_STRIDEDMAP_DECL_HPP

#include <Kokkos_DefaultNode.hpp>
#include <Teuchos_Describable.hpp>

#include "Xpetra_ConfigDefs.hpp"
#include "Xpetra_Map.hpp"

namespace Xpetra {

/*!
  @class StridedMap
  @brief Class that stores a strided map

  StridedMap extends the functionality of Xpetra::Map

  It derives from Xpetra::Map and adds a std::vector, which contains the striding information.
  E.g. for a strided map with 3dofs per node (2 velocity dofs, 1 pressure dof) the striding
  information looks like:
  std::vector<size_t> stridingInformation;
  stridingInformation.push_back(2); // 2 velocity dofs
  stridingInformation.push_back(1); // 1 pressure dof

  For this example the getFixedBlockSize() returns 3 (3 dofs per node).
  Providing a stridedBlockId parameter in the constructor the strided map only contains dofs of
  one strided block, e.g. with above stridingInformation the call

  StridingMap M(33,0,stridiningInformation,comm,0); // striding block 0 (velocity dofs)
  returns a map with the gids
  0, 1, 3, 4, 6, 7, ... (which contains only the velocity dofs)

  and
  StridingMap M(33,0,stridiningInformation,comm,1); // striding block 1 (pressure dofs)
  creates a map with only the pressure dofs
  2, 5, 8, ...

  @note: there's no support for global offset, yet.
*/
template<class LocalOrdinal, 
         class GlobalOrdinal, 
         class Node = KokkosClassic::DefaultNode::DefaultNodeType>
class StridedMap : public virtual Map<LocalOrdinal, GlobalOrdinal, Node>
{


  public:


    typedef LocalOrdinal  local_ordinal_type;
    typedef GlobalOrdinal global_ordinal_type;
    typedef Node          node_type;


  private:


#undef XPETRA_STRIDEDMAP_SHORT
#include "Xpetra_UseShortNamesOrdinal.hpp"


  public:


    //! @name Constructor/Destructor Methods
    //@{

    /** \brief Map constructor with contiguous uniform distribution.
     *
     *  Map constructor with contiguous uniform distribution.
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
#ifdef TPETRA_ENABLE_DEPRECATED_CODE
    TPETRA_DEPRECATED
    StridedMap(UnderlyingLib                                 xlib,
               global_size_t                                 numGlobalElements,
               GlobalOrdinal                                 indexBase,
               std::vector<size_t>&                          stridingInfo,
               const Teuchos::RCP<const Teuchos::Comm<int>>& comm,
               LocalOrdinal                                  stridedBlockId,      // FIXME (mfh 03 Sep 2014) This breaks for unsigned LocalOrdinal
               GlobalOrdinal                                 offset,
               LocalGlobal                                   lg,
               const Teuchos::RCP<Node>&                     /* node */);
#endif      // #ifdef TPETRA_ENABLE_DEPRECATED_CODE


    StridedMap(UnderlyingLib                                 xlib,
               global_size_t                                 numGlobalElements,
               GlobalOrdinal                                 indexBase,
               std::vector<size_t>&                          stridingInfo,
               const Teuchos::RCP<const Teuchos::Comm<int>>& comm,
               LocalOrdinal                                  stridedBlockId = -1,      // FIXME (mfh 03 Sep 2014) This breaks for unsigned LocalOrdinal
               GlobalOrdinal                                 offset         = 0,
               LocalGlobal                                   lg             = GloballyDistributed);


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
#ifdef TPETRA_ENABLE_DEPRECATED_CODE
    TPETRA_DEPRECATED
    StridedMap(UnderlyingLib                                 xlib,
               global_size_t                                 numGlobalElements,
               size_t                                        numLocalElements,
               GlobalOrdinal                                 indexBase,
               std::vector<size_t>&                          stridingInfo,
               const Teuchos::RCP<const Teuchos::Comm<int>>& comm,
               LocalOrdinal                                  stridedBlockId,
               GlobalOrdinal                                 offset,
               const Teuchos::RCP<Node>&                     /* node */);
#endif      // TPETRA_ENABLE_DEPRECATED_CODE


    StridedMap(UnderlyingLib                                 xlib,
               global_size_t                                 numGlobalElements,
               size_t                                        numLocalElements,
               GlobalOrdinal                                 indexBase,
               std::vector<size_t>&                          stridingInfo,
               const Teuchos::RCP<const Teuchos::Comm<int>>& comm,
               LocalOrdinal                                  stridedBlockId = -1,
               GlobalOrdinal                                 offset         = 0);


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
#ifdef TPETRA_ENABLE_DEPRECATED_CODE
    TPETRA_DEPRECATED
    StridedMap(UnderlyingLib                                  xlib,
               global_size_t                                  numGlobalElements,
               const Teuchos::ArrayView<const GlobalOrdinal>& elementList,
               GlobalOrdinal                                  indexBase,
               std::vector<size_t>&                           stridingInfo,
               const Teuchos::RCP<const Teuchos::Comm<int>>&  comm,
               LocalOrdinal                                   stridedBlockId,
               const Teuchos::RCP<Node>&                      /* node */);
#endif      // TPETRA_ENABLE_DEPRECATED_CODE


    StridedMap(UnderlyingLib                                  xlib,
               global_size_t                                  numGlobalElements,
               const Teuchos::ArrayView<const GlobalOrdinal>& elementList,
               GlobalOrdinal                                  indexBase,
               std::vector<size_t>&                           stridingInfo,
               const Teuchos::RCP<const Teuchos::Comm<int>>&  comm,
               LocalOrdinal                                   stridedBlockId = -1);


    StridedMap(const RCP<const Map>& map,
               std::vector<size_t>&  stridingInfo,
               GlobalOrdinal /* indexBase */,
               LocalOrdinal  stridedBlockId = -1,
               GlobalOrdinal offset         = 0);


    //! Destructor.
    virtual ~StridedMap();

    //@}

    //! @name Access functions for striding data
    //@{


    std::vector<size_t> getStridingData() const;


    void setStridingData(std::vector<size_t> stridingInfo);


    size_t getFixedBlockSize() const;


    /// returns strided block id of the dofs stored in this map
    /// or -1 if full strided map is stored in this map
    LocalOrdinal getStridedBlockId() const;


    /// returns true, if this is a strided map (i.e. more than 1 strided blocks)
    bool isStrided() const;


    /// returns true, if this is a blocked map (i.e. more than 1 dof per node)
    /// either strided or just 1 block per node
    bool isBlocked() const; 


    GlobalOrdinal getOffset() const;


    void setOffset(GlobalOrdinal offset);

    
    // returns number of strided block id which gid belongs to.
    size_t GID2StridingBlockId(GlobalOrdinal gid) const;


    //! @name Xpetra specific
    //@{


    RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node>> getMap() const;


#ifdef HAVE_XPETRA_KOKKOS_REFACTOR
    #ifdef HAVE_XPETRA_TPETRA
        using local_map_type = typename Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node>::local_map_type;

        /// \brief Get the local Map for Kokkos kernels.
        local_map_type getLocalMap() const 
        { 
            return map_->getLocalMap(); 
        }
    #else      // HAVE_XPETRA_TPETRA
        #ifdef __GNUC__
            #warning \
              "Xpetra Kokkos interface for CrsMatrix is enabled (HAVE_XPETRA_KOKKOS_REFACTOR) but Tpetra is disabled. The Kokkos interface needs Tpetra to be enabled, too."
        #endif      // __GNUC__
    #endif      // HAVE_XPETRA_TPETRA ELSE
#endif      // HAVE_XPETRA_KOKKOS_REFACTOR


    //@}

    /* // function currently not needed but maybe useful
    std::vector<GlobalOrdinal> NodeId2GlobalDofIds(GlobalOrdinal nodeId) const {
      TEUCHOS_TEST_FOR_EXCEPTION(stridingInfo_.size() == 0, Exceptions::RuntimeError, "StridedMap::NodeId2GlobalDofIds:
    stridingInfo not valid: stridingInfo.size() = 0?"); std::vector<GlobalOrdinal> dofs; if(stridedBlockId_ > -1) {
          TEUCHOS_TEST_FOR_EXCEPTION(stridingInfo_[stridedBlockId_] == 0, Exceptions::RuntimeError,
    "StridedMap::NodeId2GlobalDofIds: stridingInfo not valid: stridingInfo[stridedBlockId] = 0?");

          // determine nStridedOffset
          size_t nStridedOffset = 0;
          for(int j=0; j<stridedBlockId_; j++) {
            nStridedOffset += stridingInfo_[j];
          }

          for(size_t i = 0; i<stridingInfo_[stridedBlockId_]; i++) {
            GlobalOrdinal gid =
                nodeId * Teuchos::as<GlobalOrdinal>(getFixedBlockSize()) +
                offset_ +
                Teuchos::as<GlobalOrdinal>(nStridedOffset) +
                Teuchos::as<GlobalOrdinal>(i);
            dofs.push_back(gid);
          }
      } else {
        for(size_t i = 0; i<getFixedBlockSize(); i++) {
          GlobalOrdinal gid =
              nodeId * Teuchos::as<GlobalOrdinal>(getFixedBlockSize()) +
              offset_ +
              Teuchos::as<GlobalOrdinal>(i);
          dofs.push_back(gid);
        }
      }
      return dofs;
    }*/
    //@}


  private:


    virtual bool CheckConsistency();


  private:


    RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node>> map_;

    std::vector<size_t> stridingInfo_;        //!< vector with size of strided blocks (dofs)
    LocalOrdinal        stridedBlockId_;      //!< member variable denoting which dofs are stored in map
                                              //     stridedBlock == -1: the full map (with all strided block dofs)
                                              //     stridedBlock  > -1: only dofs of strided block with index "stridedBlockId" are
                                              //     stored in this map
    GlobalOrdinal offset_;         //!< offset for gids in map (default = 0)
    GlobalOrdinal indexBase_;      //!< index base for the strided map (default = 0)


  public:


    //! @name Map Attribute Methods
    //@{

    //! Returns the number of elements in this Map.
    global_size_t getGlobalNumElements() const;


    //! Returns the number of elements belonging to the calling node.
    size_t getNodeNumElements() const;


    //! Returns the index base for this Map.
    GlobalOrdinal getIndexBase() const;


    //! Returns minimum local index.
    LocalOrdinal getMinLocalIndex() const;


    //! Returns maximum local index.
    LocalOrdinal getMaxLocalIndex() const;


    //! Returns minimum global index owned by this node.
    GlobalOrdinal getMinGlobalIndex() const;


    //! Returns maximum global index owned by this node.
    GlobalOrdinal getMaxGlobalIndex() const;


    //! Return the minimum global index over all nodes.
    GlobalOrdinal getMinAllGlobalIndex() const;


    //! Return the maximum global index over all nodes.
    GlobalOrdinal getMaxAllGlobalIndex() const;


    //! Return the local index for a given global index.
    LocalOrdinal getLocalElement(GlobalOrdinal globalIndex) const;


    //! Return the global index for a given local index.
    GlobalOrdinal getGlobalElement(LocalOrdinal localIndex) const;


    //! Returns the node IDs and corresponding local indices for a given list of global indices.
    LookupStatus getRemoteIndexList(const Teuchos::ArrayView<const GlobalOrdinal>& GIDList,
                                    const Teuchos::ArrayView<int>&                 nodeIDList,
                                    const Teuchos::ArrayView<LocalOrdinal>&        LIDList) const;


    //! Returns the node IDs for a given list of global indices.
    LookupStatus getRemoteIndexList(const Teuchos::ArrayView<const GlobalOrdinal>& GIDList,
                                    const Teuchos::ArrayView<int>&                 nodeIDList) const;

    
    //! Return a list of the global indices owned by this node.
    Teuchos::ArrayView<const GlobalOrdinal> getNodeElementList() const;


    //! Returns true if the local index is valid for this Map on this node; returns false if it isn't.
    bool isNodeLocalElement(LocalOrdinal localIndex) const;


    //! Returns true if the global index is found in this Map on this node; returns false if it isn't.
    bool isNodeGlobalElement(GlobalOrdinal globalIndex) const;


    //! Returns true if this Map is distributed contiguously; returns false otherwise.
    bool isContiguous() const;


    //! Returns true if this Map is distributed across more than one node; returns false otherwise.
    bool isDistributed() const;

    //@}

    //! Returns true if map is compatible with this Map.
    bool isCompatible(const Map& map) const;


    //! Returns true if map is identical to this Map.
    bool isSameAs(const Map& map) const;


    //! Get the Comm object for this Map.
    Teuchos::RCP<const Teuchos::Comm<int>> getComm() const;


#ifdef TPETRA_ENABLE_DEPRECATED_CODE
    //! Get the Node object for this Map.
    Teuchos::RCP<Node> getNode() const;
#endif      // TPETRA_ENABLE_DEPRECATED_CODE


    RCP<const Map> removeEmptyProcesses() const;


    RCP<const Map> replaceCommWithSubset(const Teuchos::RCP<const Teuchos::Comm<int>>& newComm) const;


    //! Return a simple one-line description of this object.
    std::string description() const;


    //! Print the object with some verbosity level to a FancyOStream object.
    void describe(Teuchos::FancyOStream&         out,
                  const Teuchos::EVerbosityLevel verbLevel = Teuchos::Describable::verbLevel_default) const;


    //! Get the library used by this object (Tpetra or Epetra?)
    UnderlyingLib lib() const;


};      // StridedMap class



}      // namespace Xpetra

#define XPETRA_STRIDEDMAP_SHORT
#endif      // XPETRA_STRIDEDMAP_DECL_HPP


