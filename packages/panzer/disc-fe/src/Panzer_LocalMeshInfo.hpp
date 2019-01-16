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

#ifndef PANZER_LOCAL_MESH_INFO_HPP
#define PANZER_LOCAL_MESH_INFO_HPP

// Trilinos includes
#include "Kokkos_View.hpp"
#include "Phalanx_KokkosDeviceTypes.hpp"
#include "Shards_CellTopology.hpp"
#include "Teuchos_RCP.hpp"

// Panzer includes
#include "PanzerCore_config.hpp"

// STL includes
#include <string>

namespace panzer
{

class ConnManager;

/**
 * \class ElementSets
 *
 * \brief Class used to define block/sideset/nodeset ownership of elements
 *
 * Eventually this will be embedded in some kind of partitionable connectivity manager.
 */
class
ElementSets
{
public:

  using SetIndex = int;

  ElementSets() = default;
  ElementSets(const ElementSets & sets) = default;

  ElementSets(const std::map<std::string, SetIndex> & set_id_map,
              const PHX::View<const SetIndex*> & element_set_ids)
  {
    set_to_index_map_ = set_id_map;
    elements_ = element_set_ids;
    index_to_set_map_.clear();
    for(const auto & pr : set_id_map)
      index_to_set_map_[pr.second] = pr.first;
  }

  ~ElementSets() = default;

  int
  getNumElements() const
  {return elements_.extent(0);}

  SetIndex
  getElementSetIndex(const int & i) const
  {return elements_(i);}

  const std::string &
  getElementSet(const SetIndex & id) const
  {return getSetName(getElementSetIndex(id));}

  const std::string &
  getSetName(const SetIndex & set_index) const
  {return index_to_set_map_.at(set_index);}

  int
  getSetIndex(const std::string & set) const
  {return set_to_index_map_.at(set);}

  Teuchos::RCP<ElementSets>
  buildSubsets(const std::vector<int> & element_ids) const
  {
    PHX::View<SetIndex*> new_element_block_ids("element_block_ids",element_ids.size());
    for(unsigned int i=0; i<element_ids.size(); ++i)
      new_element_block_ids(i) = elements_(element_ids[i]);
    return Teuchos::rcp(new ElementSets(set_to_index_map_, new_element_block_ids));
  }

  std::vector<std::string>
  getSetNames() const
  {
    std::vector<std::string> sets;
    for(const auto & pr : set_to_index_map_)
      sets.push_back(pr.first);
    return sets;
  }

  bool
  hasSet(const std::string & set) const
  {return set_to_index_map_.find(set) != set_to_index_map_.end();}

protected:

  std::map<std::string, SetIndex> set_to_index_map_;
  std::map<SetIndex,std::string> index_to_set_map_;
  PHX::View<const SetIndex *> elements_;

};

/** Base class for LocalMeshInfo structures
 *
 * TODO: Replace with Connectivity manager
 */
struct LocalMeshInfoBase
{

  LocalMeshInfoBase():
    subcell_dimension(-1),
    subcell_index(-1),
    num_owned_cells(0),
    num_ghost_cells(0),
    num_virtual_cells(0),
    has_connectivity(false)
  {

  }

  int subcell_dimension;
  int subcell_index;

  panzer::LocalOrdinal num_owned_cells;
  panzer::LocalOrdinal num_ghost_cells;
  panzer::LocalOrdinal num_virtual_cells;

  // Global cell indexes -> [owned] then [ghosted] then [virtual]
  Kokkos::View<panzer::GlobalOrdinal*> global_cells;

  // These are the cell indexes in the LocalMeshInfo class
  Kokkos::View<panzer::LocalOrdinal*> local_cells;

  // Vertices
  Kokkos::View<double***,PHX::Device> cell_vertices;

  // Has the connectivity been constructed
  bool has_connectivity;

  // Connectivity information

  // Face to neighbors
  Kokkos::View<panzer::LocalOrdinal*[2]> face_to_cells;
  Kokkos::View<panzer::LocalOrdinal*[2]> face_to_lidx;
  Kokkos::View<panzer::LocalOrdinal**> cell_to_faces;

  // Used to define element block for each cell
  // TODO: Eventually we will have one for nodesets, edgesets, sidesets, and cellsets (element blocks)
  Teuchos::RCP<ElementSets> cell_sets;

};

/** Partition of LocalMeshInfo
 *
 * Used for generating worksets
 *
 * TODO: this is literally just a rename of the 'LocalMeshSidesetInfo' class
 *
 */
struct LocalMeshPartition:
    public LocalMeshInfoBase
{

  std::string element_block_name;
  Teuchos::RCP<const shards::CellTopology> cell_topology;

  // In case this is a sideset
  std::string sideset_name;

};

/** Portion of LocalMeshInfo associated with element block
 *
 * Used to represent an element block found on the local process
 *
 */
struct LocalMeshBlockInfo:
    public LocalMeshInfoBase
{
  std::string element_block_name;

  Teuchos::RCP<const shards::CellTopology> cell_topology;

};

/** Portion of LocalMeshInfo associated with sidesets
 *
 * Used to represent a sideset found on the local process
 *
 */
struct LocalMeshSidesetInfo:
    public LocalMeshBlockInfo
{
  std::string sideset_name;
};

/** Entire mesh found on a local process
 *
 */
struct LocalMeshInfo:
    public LocalMeshInfoBase
{

  /// Default constructor
  LocalMeshInfo() = default;

  /**
   * \brief Used to initialize cell/element block data
   *
   * \note THIS DOES NOT SETUP 'sidesets' IN THE LocalMeshInfo OBJECT - currently relies on STK for sideset construction
   *
   * \param[in] comm Communicator on which mesh is distributed
   * \param[in] conn Connectivity manager storing topological information
   * \param[in] owned_cell_vertices Cell vertices for OWNED cells only
   *
   */
  void
  initialize(const Teuchos::RCP<const Teuchos::Comm<int> > & comm,
             ConnManager & conn,
             const PHX::View<const double ***> & owned_cell_vertices);

  /// Default destructor
  ~LocalMeshInfo() = default;

  // Element block -> block info
  std::map<std::string, LocalMeshBlockInfo > element_blocks;

  // Element block, sideset -> sideset info
  std::map<std::string, std::map<std::string, LocalMeshSidesetInfo > > sidesets;

};

}

#endif
