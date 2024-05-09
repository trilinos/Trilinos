// @HEADER
// ***********************************************************************
//
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//                 Copyright (2011) Sandia Corporation
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
// Questions? Contact Roger P. Pawlowski (rppawlo@sandia.gov) and
// Eric C. Cyr (eccyr@sandia.gov)
// ***********************************************************************
// @HEADER

#ifndef __Panzer_ConnManager_hpp__
#define __Panzer_ConnManager_hpp__

#include <vector>

#include "Kokkos_StaticCrsGraph.hpp"

// Teuchos includes
#include "Teuchos_RCP.hpp"
#include "Shards_CellTopology.hpp"
#include "PanzerDofMgr_config.hpp"

namespace panzer {

  class FieldPattern; // from DOFManager

  /// Pure virtual base class for supplying mesh connectivity information to the DOF Manager.
  class ConnManager {
  public:

    using GlobalOrdinal = panzer::GlobalOrdinal;
    using LocalOrdinal = int;

    virtual ~ConnManager() = default;

    /** Tell the connection manager to build the connectivity assuming
     * a particular field pattern.
     *
     * \param[in] fp Field pattern to build connectivity for
     */
    virtual void buildConnectivity(const FieldPattern & fp) = 0;

    /** Build a clone of this connection manager, without any assumptions
     * about the required connectivity (i.e. <code>buildConnectivity</code>
     * has never been called).
     */
    virtual Teuchos::RCP<ConnManager> noConnectivityClone() const = 0;

    /** How many mesh IDs are associated with this element?
     *
     * \param[in] localElmtId Local element ID
     *
     * \returns Number of mesh IDs that are associated with this element.
     */
    virtual LocalOrdinal getConnectivitySize(LocalOrdinal localElmtId) const = 0;

    /** Get ID connectivity for a particular element
     *
     * \param[in] localElmtId Local element ID
     *
     * \returns Pointer to beginning of indices, with total size
     *          equal to <code>getConnectivitySize(localElmtId)</code>
     */
    virtual const GlobalOrdinal * getConnectivity(LocalOrdinal localElmtId) const = 0;

    /** Get the block ID for a particular element.
     *
     * \param[in] localElmtId Local element ID
     */
    virtual std::string getBlockId(LocalOrdinal localElmtId) const = 0;

    /** Returns the number of element blocks in this mesh */
    virtual size_t numElementBlocks() const = 0;

    /** What are the blockIds included in this connection manager */
    virtual void getElementBlockIds(std::vector<std::string> & elementBlockIds) const = 0;

    /** Returns the cellTopologies linked to element blocks in this connection manager */
    virtual void getElementBlockTopologies(std::vector<shards::CellTopology> & elementBlockTopologies) const = 0;

    /** Get the local element IDs for a paricular element
     * block.
     *
     * \param[in] blockID Block ID
     *
     * \returns Vector of local element IDs.
     */
    virtual const std::vector<LocalOrdinal> & getElementBlock(const std::string & blockID) const = 0;

    /** Get the local element IDs for all "neighbor" elements that reside in a paricular element
     * block (An element is a neighbor if it is in the one ring of owned elements).
     *
     * \param[in] blockID Block ID
     *
     * \returns Vector of local element IDs.
     */
    virtual const std::vector<LocalOrdinal> & getNeighborElementBlock(const std::string & blockID) const = 0;

    /** Get elements, if any, associated with <code>el</code>, excluding
     * <code>el</code> itself.
     */
    virtual const std::vector<LocalOrdinal>& getAssociatedNeighbors(const LocalOrdinal& el) const = 0;

    /** Return whether getAssociatedNeighbors will return true for at least one
     * input.
     */
    virtual bool hasAssociatedNeighbors() const = 0;
  };

} // namespace panzer

namespace panzer::Experimental {

  template <typename DataType, typename... Properties>
  using RaggedView = Kokkos::StaticCrsGraph<DataType, Properties...>;

  template <typename LocalOrdinal, typename GlobalOrdinal, typename DeviceType>
  class ConnManager : public virtual ::panzer::ConnManager {
    static_assert(std::is_same_v<LocalOrdinal,  int>);
    static_assert(std::is_same_v<GlobalOrdinal, panzer::GlobalOrdinal>);

  public:
    using local_ordinal_type = LocalOrdinal;
    using global_ordinal_type = GlobalOrdinal;
    using device_type = DeviceType;
    using execution_space = typename device_type::execution_space;
    using memory_space = typename device_type::memory_space;
    using node_type = Tpetra::KokkosCompat::KokkosDeviceWrapperNode<execution_space, memory_space>;

    using conn_manager_base_type = ::panzer::ConnManager;

    using connectivity_device_view_type = RaggedView<global_ordinal_type, memory_space>;
    using connectivity_host_view_type = typename connectivity_device_view_type::HostMirror;

    using block_elmids_device_view_type = Kokkos::View<local_ordinal_type*, memory_space>;
    using block_elmids_host_view_type = typename block_elmids_device_view_type::HostMirror;

  public:
    virtual ~ConnManager() = default;

  public:
    virtual connectivity_device_view_type getConnectivityDevice() const {
      const auto connectivity_h = this->getConnectivityHost();
      return connectivity_device_view_type(
        Kokkos::create_mirror_view_and_copy(Kokkos::view_alloc(memory_space{}), connectivity_h.entries),
        Kokkos::create_mirror_view_and_copy(Kokkos::view_alloc(memory_space{}), connectivity_h.row_map)
      );
    }

    virtual typename block_elmids_device_view_type::const_type getElementBlockDevice(const std::string& blockID) const {
      return Kokkos::create_mirror_view_and_copy(Kokkos::view_alloc(memory_space{}), this->getElementBlockHost(blockID));
    }

    virtual typename block_elmids_device_view_type::const_type getNeighborElementBlockDevice(const std::string& blockID) const {
      return Kokkos::create_mirror_view_and_copy(Kokkos::view_alloc(memory_space{}), this->getNeighborElementBlockHost(blockID));
    }

    virtual typename block_elmids_device_view_type::const_type getAssociatedNeighborsDevice(local_ordinal_type localElmtId) const {
      return Kokkos::create_mirror_view_and_copy(Kokkos::view_alloc(memory_space{}), this->getAssociatedNeighborsHost(localElmtId));
    }

  private:
    connectivity_host_view_type getConnectivityHost() const {

      local_ordinal_type maxLocalElmtId = 0;

      std::vector<std::string> elmtBlockIds;
      this->getElementBlockIds(elmtBlockIds);

      std::for_each(
        elmtBlockIds.cbegin(), elmtBlockIds.cend(),
        [&maxLocalElmtId, this] (const auto& elmtBlockId) {
          std::vector<std::vector<local_ordinal_type>> elmtBlocks{
            this->getElementBlock(elmtBlockId),
            this->getNeighborElementBlock(elmtBlockId)
          };
          std::for_each(
            elmtBlocks.cbegin(), elmtBlocks.cend(),
            [&maxLocalElmtId] (const auto& elmtBlock) {
              if ( ! elmtBlock.empty())
                maxLocalElmtId = std::max(maxLocalElmtId, *std::max_element(elmtBlock.cbegin(), elmtBlock.cend()));
            }
          );
        }
      );

      typename connectivity_host_view_type::row_map_type::non_const_type offsets_h(
        "offsets of ragged view representation of connectivity",
        maxLocalElmtId + 2
      );

      local_ordinal_type numEntries;

      Kokkos::parallel_scan(
        "compute offsets of ragged view representation of connectivity",
        Kokkos::RangePolicy<Kokkos::DefaultHostExecutionSpace>(0, maxLocalElmtId + 1),
        [=] (local_ordinal_type localElmtId, local_ordinal_type& partialSum, bool isFinal) {
          partialSum += this->getConnectivitySize(localElmtId);
          if (isFinal) offsets_h(localElmtId + 1) = partialSum;
        }, numEntries
      );

      typename connectivity_host_view_type::entries_type entries_h(
        Kokkos::view_alloc("entries of ragged view representation of connectivity", Kokkos::WithoutInitializing),
        numEntries
      );

      std::for_each(
        elmtBlockIds.cbegin(), elmtBlockIds.cend(),
        [=] (const auto& elmtBlockId) {
          std::vector<std::vector<local_ordinal_type>> elmtBlocks{
            this->getElementBlock(elmtBlockId),
            this->getNeighborElementBlock(elmtBlockId)
          };
      
          std::for_each(
            elmtBlocks.cbegin(), elmtBlocks.cend(),
            [=] (const auto& elmtBlock) {
              Kokkos::parallel_for(
                "fill entries of ragged view representation of connectivity",
                Kokkos::RangePolicy<Kokkos::DefaultHostExecutionSpace>(0, elmtBlock.size()),
                [=] (typename std::vector<local_ordinal_type>::size_type idx) {
                  const local_ordinal_type localElmtId = elmtBlock[idx];
                  const auto offset = offsets_h(localElmtId);
                  for (typename connectivity_host_view_type::row_map_type::non_const_value_type jdx = 0; jdx < offsets_h(localElmtId + 1) - offsets_h(localElmtId); ++jdx) {
                    entries_h(offset + jdx) = this->getConnectivity(localElmtId)[jdx];
                  }
                }
              );
            }
          );
        }
      );

      return connectivity_host_view_type(std::move(entries_h), std::move(offsets_h));
    }

    auto getElementBlockHost(const std::string& blockId) const {
      const auto& elmtBlock = this->getElementBlock(blockId);
      return typename block_elmids_host_view_type::const_type(elmtBlock.data(), elmtBlock.size());
    }

    auto getNeighborElementBlockHost(const std::string& blockId) const {
      const auto& elmtBlock = this->getNeighborElementBlock(blockId);
      return typename block_elmids_host_view_type::const_type(elmtBlock.data(), elmtBlock.size());
    }

    auto getAssociatedNeighborsHost(local_ordinal_type localElmtId) const {
      const auto& elmtBlock = this->getAssociatedNeighbors(localElmtId);
      return typename block_elmids_host_view_type::const_type(elmtBlock.data(), elmtBlock.size());
    }
  };

} // namespace panzer::Experimental

#endif
