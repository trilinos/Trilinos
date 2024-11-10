// Copyright(C) 1999-2020, 2022 National Technology & Engineering Solutions
// of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
// NTESS, the U.S. Government retains certain rights in this software.
//
// See packages/seacas/LICENSE for details

#pragma once

#include "iotm_export.h"

#include "Ioss_CodeTypes.h"
#include "Ioss_EntityType.h" // for EntityType

#include <cstddef> // for size_t
#include <cstdint> // for int64_t
#include <map>     // for map, etc
#include <string>  // for string
#include <unordered_map>
#include <utility> // for pair
#include <vector>  // for vector

#include <assert.h>
#include <fmt/ostream.h>

#include "Ioss_ElementPermutation.h"
#include "Ioss_ElementTopology.h"
#include "Ioss_StandardElementTypes.h"
#include "Ioss_Utils.h"

#include "Iotm_TextMeshDataTypes.h"

namespace Iotm {
  class IOTM_EXPORT TopologyMapEntry
  {
  public:
    using Ordinal     = uint16_t;
    using Permutation = uint8_t;

    static constexpr Ordinal     InvalidOrdinal     = 65535;
    static constexpr Permutation InvalidPermutation = 128;

    using DimensionArray = bool[4];

    TopologyMapEntry()
        : id(Ioss::ElementTopology::get_unique_id(std::string(Ioss::Unknown::name))),
          topology(Ioss::ElementTopology::factory(std::string(Ioss::Unknown::name))),
          initialized(false)
    {
      set_valid_spatial_dimensions({false, false, false, false});
    }

    explicit TopologyMapEntry(const std::string &name)
        : id(Ioss::ElementTopology::get_unique_id(name)),
          topology(Ioss::ElementTopology::factory(name)), initialized(false)
    {
      set_valid_spatial_dimensions({false, false, false, false});
    }

    TopologyMapEntry(const TopologyMapEntry &topo)            = default;
    TopologyMapEntry &operator=(const TopologyMapEntry &topo) = default;

    bool operator==(const Ioss::ElementTopology *topo) const { return topo == topology; }

    bool defined_on_spatial_dimension(const unsigned spatialDim) const
    {
      if (spatialDim > 3) {
        return false;
      }
      return validSpatialDimensions[spatialDim];
    }

    const std::string name() const { return topology->name(); }

    int num_nodes() const { return topology->number_nodes(); }

    bool operator==(const TopologyMapEntry &rhs) const
    {
      return id == rhs.id && topology == rhs.topology &&
             equivalent_valid_spatial_dimensions(rhs.validSpatialDimensions);
    }

    bool operator!=(const TopologyMapEntry &rhs) const { return !(*this == rhs); }

    int num_sides() const
    {
      if (topology->is_shell()) {
        // Only interested in face boundaries, not edges
        if (topology->parametric_dimension() == 2) {
          return topology->number_faces();
        }
      }

      return topology->number_boundaries();
    }

    // Side references are one-based
    bool valid_side(unsigned side) const
    {
      unsigned numSides = num_sides();
      if (side > 0 && side <= numSides)
        return true;
      return false;
    }

    std::string side_topology_name(unsigned side) const
    {
      if (!valid_side(side))
        return "";

      Ioss::ElementTopology *sideTopology = topology->boundary_type(side);
      return sideTopology->name();
    }

    const TopologyMapEntry &side_topology(unsigned side) const
    {
      if (!valid_side(side))
        return *(invalid_topology_factory());
      return *sideTopologies[side - 1];
    }

    unsigned side_topology_num_nodes(unsigned side) const
    {
      if (!valid_side(side))
        return 0;

      Ioss::ElementTopology *sideTopology = topology->boundary_type(side);
      return sideTopology->number_nodes();
    }

    std::vector<Ordinal> side_topology_node_indices(unsigned side) const
    {
      if (!valid_side(side))
        return std::vector<Ordinal>();

      Ioss::ElementTopology *sideTopology = topology->boundary_type(side);
      std::vector<Ordinal>   elementNodeOrdinalVector(sideTopology->number_nodes());

      Ioss::IntVector connectivity = topology->boundary_connectivity(side);

      for (int i = 0; i < sideTopology->number_nodes(); i++) {
        elementNodeOrdinalVector[i] = connectivity[i];
      }

      return elementNodeOrdinalVector;
    }

    bool is_shell() const { return topology->is_shell(); }

    unsigned num_permutations() const { return topology->permutation()->num_permutations(); }

    unsigned num_positive_permutations() const
    {
      return topology->permutation()->num_positive_permutations();
    }

    bool is_positive_polarity(Permutation permutation) const
    {
      return topology->permutation()->is_positive_polarity(permutation);
    }

    bool valid_permutation(Permutation permutation) const
    {
      return topology->permutation()->valid_permutation(permutation);
    }

    bool fill_permutation_indices(Permutation           permutation,
                                  std::vector<Ordinal> &nodeOrdinalVector) const
    {
      return topology->permutation()->fill_permutation_indices(permutation, nodeOrdinalVector);
    }

    std::vector<Ordinal> permutation_indices(Permutation permutation) const
    {
      return topology->permutation()->permutation_indices(permutation);
    }

    static TopologyMapEntry *invalid_topology_factory()
    {
      static TopologyMapEntry entry;

      entry.initialized = true;

      return &entry;
    }

    // Node with no permutation and no sides
    static TopologyMapEntry *node_factory()
    {
      static TopologyMapEntry entry(Ioss::Node::name);

      if (!entry.initialized) {
        entry.set_valid_spatial_dimensions({false, true, true, true});
        entry.set_side_topologies({});
        entry.initialized = true;
      }

      return &entry;
    }

    //***************************************************************************
    // topology::LINE -- topology::EDGE_RANK
    // 2 nodes with no sides defined in 2D/3D
    //***************************************************************************
    static TopologyMapEntry *line_2_factory()
    {
      static TopologyMapEntry entry(Ioss::Edge2::name);

      if (!entry.initialized) {
        entry.set_valid_spatial_dimensions({false, false, true, true});
        entry.set_side_topologies({});
        entry.initialized = true;
      }

      return &entry;
    }

    //***************************************************************************
    // topology::LINE 1D -- topology::ELEMENT_RANK
    // 2 nodes with no sides only defined on 1d problems
    //***************************************************************************
    static TopologyMapEntry *line_2_1d_factory()
    {
      static TopologyMapEntry entry(Ioss::Edge2::name);

      if (!entry.initialized) {
        entry.set_valid_spatial_dimensions({false, true, false, false});
        entry.set_side_topologies({});
        entry.initialized = true;
      }

      return &entry;
    }

    //***************************************************************************
    // topology::LINE -- topology::EDGE_RANK
    // 3 nodes with no sides defined in 2D/3D
    //***************************************************************************
    static TopologyMapEntry *line_3_factory()
    {
      static TopologyMapEntry entry(Ioss::Edge3::name);

      if (!entry.initialized) {
        entry.set_valid_spatial_dimensions({false, false, true, true});
        entry.set_side_topologies({});
        entry.initialized = true;
      }

      return &entry;
    }

    //***************************************************************************
    // topology::LINE 1D -- topology::ELEMENT_RANK
    // 3 nodes with no sides only defined on 1d problems
    //***************************************************************************
    static TopologyMapEntry *line_3_1d_factory()
    {
      static TopologyMapEntry entry(Ioss::Edge3::name);

      if (!entry.initialized) {
        entry.set_valid_spatial_dimensions({false, true, false, false});
        entry.set_side_topologies({});
        entry.initialized = true;
      }

      return &entry;
    }

    //***************************************************************************
    // topology::TRIANGLE -- topology::FACE_RANK
    // defined on spatial dimension 3d
    // 3, 4, or 6 nodes with 3 edges and no sides
    //***************************************************************************
    static TopologyMapEntry *tri_3_factory()
    {
      static TopologyMapEntry entry(Ioss::Tri3::name);

      if (!entry.initialized) {
        entry.set_valid_spatial_dimensions({false, false, false, true});
        entry.set_side_topologies({});
        entry.initialized = true;
      }

      return &entry;
    }

    static TopologyMapEntry *tri_4_factory()
    {
      static TopologyMapEntry entry(Ioss::Tri4::name);

      if (!entry.initialized) {
        entry.set_valid_spatial_dimensions({false, false, false, true});
        entry.set_side_topologies({});
        entry.initialized = true;
      }

      return &entry;
    }

    static TopologyMapEntry *tri_6_factory()
    {
      static TopologyMapEntry entry(Ioss::Tri6::name);

      if (!entry.initialized) {
        entry.set_valid_spatial_dimensions({false, false, false, true});
        entry.set_side_topologies({});
        entry.initialized = true;
      }

      return &entry;
    }

    //***************************************************************************
    // topology::QUADRILATERAL -- topology::FACE_RANK
    // defined on spatial dimension 3d
    // 4, 8, or 9 nodes with 4 edges and no sides
    //***************************************************************************
    static TopologyMapEntry *quad_4_factory()
    {
      static TopologyMapEntry entry(Ioss::Quad4::name);

      if (!entry.initialized) {
        entry.set_valid_spatial_dimensions({false, false, false, true});
        entry.set_side_topologies({});
        entry.initialized = true;
      }

      return &entry;
    }

    static TopologyMapEntry *quad_6_factory()
    {
      static TopologyMapEntry entry(Ioss::Quad6::name);

      if (!entry.initialized) {
        entry.set_valid_spatial_dimensions({false, false, false, true});
        entry.set_side_topologies({});
        entry.initialized = true;
      }

      return &entry;
    }

    static TopologyMapEntry *quad_8_factory()
    {
      static TopologyMapEntry entry(Ioss::Quad8::name);

      if (!entry.initialized) {
        entry.set_valid_spatial_dimensions({false, false, false, true});
        entry.set_side_topologies({});
        entry.initialized = true;
      }

      return &entry;
    }

    static TopologyMapEntry *quad_9_factory()
    {
      static TopologyMapEntry entry(Ioss::Quad9::name);

      if (!entry.initialized) {
        entry.set_valid_spatial_dimensions({false, false, false, true});
        entry.set_side_topologies({});
        entry.initialized = true;
      }

      return &entry;
    }

    //***************************************************************************
    // PARTICLE -- topology::ELEMENT_RANK
    // one node with no sides
    //***************************************************************************
    static TopologyMapEntry *particle_factory()
    {
      static TopologyMapEntry entry(Ioss::Sphere::name);

      if (!entry.initialized) {
        entry.set_valid_spatial_dimensions({false, true, true, true});
        entry.set_side_topologies({});
        entry.initialized = true;
      }

      return &entry;
    }

    //***************************************************************************
    // topology::BEAM_2 -- topology::ELEMENT_RANK
    // 2 nodes with 2 sides defined in 2D/3D
    //***************************************************************************
    static TopologyMapEntry *beam_2_factory()
    {
      static TopologyMapEntry entry(Ioss::Beam2::name);

      if (!entry.initialized) {
        entry.set_valid_spatial_dimensions({false, false, true, true});
        entry.set_side_topologies({line_2_factory(), line_2_factory()});
        entry.initialized = true;
      }

      return &entry;
    }

    //***************************************************************************
    // topology::BEAM_3 -- topology::ELEMENT_RANK
    // 3 nodes with 2 sides defined in 2D/3D
    //***************************************************************************
    static TopologyMapEntry *beam_3_factory()
    {
      static TopologyMapEntry entry(Ioss::Beam3::name);

      if (!entry.initialized) {
        entry.set_valid_spatial_dimensions({false, false, true, true});
        entry.set_side_topologies({line_3_factory(), line_3_factory()});
        entry.initialized = true;
      }

      return &entry;
    }

    //***************************************************************************
    // topology::SHELL_LINE -- topology::ELEMENT_RANK
    // only defined on 2d problems
    // 2 or 3 nodes with two edges and 2 sides
    //***************************************************************************

    static TopologyMapEntry *shell_line_2_factory()
    {
      static TopologyMapEntry entry(Ioss::ShellLine2D2::name);

      if (!entry.initialized) {
        entry.set_valid_spatial_dimensions({false, false, true, false});
        entry.set_side_topologies({line_2_factory(), line_2_factory()});
        entry.initialized = true;
      }

      return &entry;
    }

    static TopologyMapEntry *shell_line_3_factory()
    {
      static TopologyMapEntry entry(Ioss::ShellLine2D3::name);

      if (!entry.initialized) {
        entry.set_valid_spatial_dimensions({false, false, true, false});
        entry.set_side_topologies({line_3_factory(), line_3_factory()});
        entry.initialized = true;
      }

      return &entry;
    }

    //***************************************************************************
    // topology::SPRING -- topology::ELEM_RANK
    // 2 or 3 nodes with no sides
    //***************************************************************************

    static TopologyMapEntry *spring_2_factory()
    {
      static TopologyMapEntry entry(Ioss::Spring2::name);

      if (!entry.initialized) {
        entry.set_valid_spatial_dimensions({false, true, true, true});
        entry.set_side_topologies({});
        entry.initialized = true;
      }

      return &entry;
    }

    static TopologyMapEntry *spring_3_factory()
    {
      static TopologyMapEntry entry(Ioss::Spring3::name);

      if (!entry.initialized) {
        entry.set_valid_spatial_dimensions({false, true, true, true});
        entry.set_side_topologies({});
        entry.initialized = true;
      }

      return &entry;
    }

    //***************************************************************************
    // topology::TRIANGLE 2D -- topology::ELEMENT_RANK
    // defined on spatial dimension 2d
    // 3, 4, or 6 nodes with 3 edges and 3 sides
    //***************************************************************************
    static TopologyMapEntry *tri_3_2d_factory()
    {
      static TopologyMapEntry entry(Ioss::Tri3::name);

      if (!entry.initialized) {
        entry.set_valid_spatial_dimensions({false, false, true, false});
        entry.set_side_topologies({line_2_factory(), line_2_factory(), line_2_factory()});
        entry.initialized = true;
      }

      return &entry;
    }

    static TopologyMapEntry *tri_4_2d_factory()
    {
      static TopologyMapEntry entry(Ioss::Tri4::name);

      if (!entry.initialized) {
        entry.set_valid_spatial_dimensions({false, false, true, false});
        entry.set_side_topologies({line_2_factory(), line_2_factory(), line_2_factory()});
        entry.initialized = true;
      }

      return &entry;
    }

    static TopologyMapEntry *tri_6_2d_factory()
    {
      static TopologyMapEntry entry(Ioss::Tri6::name);

      if (!entry.initialized) {
        entry.set_valid_spatial_dimensions({false, false, true, false});
        entry.set_side_topologies({line_3_factory(), line_3_factory(), line_3_factory()});
        entry.initialized = true;
      }

      return &entry;
    }

    //***************************************************************************
    // topology::QUADRILATERAL 2D -- topology::ELEMENT_RANK
    // defined on spatial dimension 2d
    // 4, 8, or 9 nodes with 4 edges and 4 sides
    //***************************************************************************

    static TopologyMapEntry *quad_4_2d_factory()
    {
      static TopologyMapEntry entry(Ioss::Quad4::name);

      if (!entry.initialized) {
        entry.set_valid_spatial_dimensions({false, false, true, false});
        entry.set_side_topologies(
            {line_2_factory(), line_2_factory(), line_2_factory(), line_2_factory()});
        entry.initialized = true;
      }

      return &entry;
    }

    static TopologyMapEntry *quad_8_2d_factory()
    {
      static TopologyMapEntry entry(Ioss::Quad8::name);

      if (!entry.initialized) {
        entry.set_valid_spatial_dimensions({false, false, true, false});
        entry.set_side_topologies(
            {line_3_factory(), line_3_factory(), line_3_factory(), line_3_factory()});
        entry.initialized = true;
      }

      return &entry;
    }

    static TopologyMapEntry *quad_9_2d_factory()
    {
      static TopologyMapEntry entry(Ioss::Quad9::name);

      if (!entry.initialized) {
        entry.set_valid_spatial_dimensions({false, false, false, true});
        entry.set_side_topologies(
            {line_3_factory(), line_3_factory(), line_3_factory(), line_3_factory()});
        entry.initialized = true;
      }

      return &entry;
    }

    //***************************************************************************
    // topology::SHELL topology::TRIANGLE -- topology::ELEMENT_RANK
    // defined on spatial dimension 3d
    // 3, 4, or 6 nodes with 3 edges and 2 sides
    //***************************************************************************

    static TopologyMapEntry *shell_tri_3_factory()
    {
      static TopologyMapEntry entry(Ioss::TriShell3::name);

      if (!entry.initialized) {
        entry.set_valid_spatial_dimensions({false, false, false, true});
        entry.set_side_topologies({tri_3_factory(), tri_3_factory()});
        entry.initialized = true;
      }

      return &entry;
    }

    static TopologyMapEntry *shell_tri_4_factory()
    {
      static TopologyMapEntry entry(Ioss::TriShell4::name);

      if (!entry.initialized) {
        entry.set_valid_spatial_dimensions({false, false, false, true});
        entry.set_side_topologies({tri_4_factory(), tri_4_factory()});
        entry.initialized = true;
      }

      return &entry;
    }

    static TopologyMapEntry *shell_tri_6_factory()
    {
      static TopologyMapEntry entry(Ioss::TriShell6::name);

      if (!entry.initialized) {
        entry.set_valid_spatial_dimensions({false, false, false, true});
        entry.set_side_topologies({tri_6_factory(), tri_6_factory()});
        entry.initialized = true;
      }

      return &entry;
    }

    //***************************************************************************
    //  topology::SHELL topology::QUADRILATERAL -- topology::ELEMENT_RANK
    // defined on spatial dimension 3d
    // 4, 8, or 9 nodes with 4 edges and 2 sides
    //***************************************************************************
    static TopologyMapEntry *shell_quad_4_factory()
    {
      static TopologyMapEntry entry(Ioss::Shell4::name);

      if (!entry.initialized) {
        entry.set_valid_spatial_dimensions({false, false, false, true});
        entry.set_side_topologies({quad_4_factory(), quad_4_factory()});
        entry.initialized = true;
      }

      return &entry;
    }

    static TopologyMapEntry *shell_quad_8_factory()
    {
      static TopologyMapEntry entry(Ioss::Shell8::name);

      if (!entry.initialized) {
        entry.set_valid_spatial_dimensions({false, false, false, true});
        entry.set_side_topologies({quad_8_factory(), quad_8_factory()});
        entry.initialized = true;
      }

      return &entry;
    }

    static TopologyMapEntry *shell_quad_9_factory()
    {
      static TopologyMapEntry entry(Ioss::Shell9::name);

      if (!entry.initialized) {
        entry.set_valid_spatial_dimensions({false, false, false, true});
        entry.set_side_topologies({quad_9_factory(), quad_9_factory()});
        entry.initialized = true;
      }

      return &entry;
    }

    //***************************************************************************
    // topology::TETRAHEDRON -- topology::ELEMENT_RANK
    // defined on spatial dimension 3d
    // 4, 8, 10 or 11 nodes with 4 sides
    //***************************************************************************
    static TopologyMapEntry *tet_4_factory()
    {
      static TopologyMapEntry entry(Ioss::Tet4::name);

      if (!entry.initialized) {
        entry.set_valid_spatial_dimensions({false, false, false, true});
        entry.set_side_topologies(
            {tri_3_factory(), tri_3_factory(), tri_3_factory(), tri_3_factory()});
        entry.initialized = true;
      }

      return &entry;
    }

    static TopologyMapEntry *tet_8_factory()
    {
      static TopologyMapEntry entry(Ioss::Tet8::name);

      if (!entry.initialized) {
        entry.set_valid_spatial_dimensions({false, false, false, true});
        entry.set_side_topologies(
            {tri_4_factory(), tri_4_factory(), tri_4_factory(), tri_4_factory()});
        entry.initialized = true;
      }

      return &entry;
    }

    static TopologyMapEntry *tet_10_factory()
    {
      static TopologyMapEntry entry(Ioss::Tet10::name);

      if (!entry.initialized) {
        entry.set_valid_spatial_dimensions({false, false, false, true});
        entry.set_side_topologies(
            {tri_6_factory(), tri_6_factory(), tri_6_factory(), tri_6_factory()});
        entry.initialized = true;
      }

      return &entry;
    }

    static TopologyMapEntry *tet_11_factory()
    {
      static TopologyMapEntry entry(Ioss::Tet11::name);

      if (!entry.initialized) {
        entry.set_valid_spatial_dimensions({false, false, false, true});
        entry.set_side_topologies(
            {tri_6_factory(), tri_6_factory(), tri_6_factory(), tri_6_factory()});
        entry.initialized = true;
      }

      return &entry;
    }

    //***************************************************************************
    // topology::PYRAMID -- topology::ELEMENT_RANK
    // defined on spatial dimension 3d
    // 5, 13 or 14 nodes with 5 sides
    //***************************************************************************
    static TopologyMapEntry *pyramid_5_factory()
    {
      static TopologyMapEntry entry(Ioss::Pyramid5::name);

      if (!entry.initialized) {
        entry.set_valid_spatial_dimensions({false, false, false, true});
        entry.set_side_topologies(
            {tri_3_factory(), tri_3_factory(), tri_3_factory(), tri_3_factory(), quad_4_factory()});
        entry.initialized = true;
      }

      return &entry;
    }

    static TopologyMapEntry *pyramid_13_factory()
    {
      static TopologyMapEntry entry(Ioss::Pyramid13::name);

      if (!entry.initialized) {
        entry.set_valid_spatial_dimensions({false, false, false, true});
        entry.set_side_topologies(
            {tri_6_factory(), tri_6_factory(), tri_6_factory(), tri_6_factory(), quad_8_factory()});
        entry.initialized = true;
      }

      return &entry;
    }

    static TopologyMapEntry *pyramid_14_factory()
    {
      static TopologyMapEntry entry(Ioss::Pyramid14::name);

      if (!entry.initialized) {
        entry.set_valid_spatial_dimensions({false, false, false, true});
        entry.set_side_topologies(
            {tri_6_factory(), tri_6_factory(), tri_6_factory(), tri_6_factory(), quad_9_factory()});
        entry.initialized = true;
      }

      return &entry;
    }

    //***************************************************************************
    // topology::WEDGE -- topology::ELEMENT_RANK
    // defined on spatial dimension 3d
    // 6, 12, 15 or 18 nodes with 5 sides
    //***************************************************************************
    static TopologyMapEntry *wedge_6_factory()
    {
      static TopologyMapEntry entry(Ioss::Wedge6::name);

      if (!entry.initialized) {
        entry.set_valid_spatial_dimensions({false, false, false, true});
        entry.set_side_topologies({quad_4_factory(), quad_4_factory(), quad_4_factory(),
                                   tri_3_factory(), tri_3_factory()});
        entry.initialized = true;
      }

      return &entry;
    }

    static TopologyMapEntry *wedge_12_factory()
    {
      static TopologyMapEntry entry(Ioss::Wedge12::name);

      if (!entry.initialized) {
        entry.set_valid_spatial_dimensions({false, false, false, true});
        entry.set_side_topologies({quad_6_factory(), quad_6_factory(), quad_6_factory(),
                                   tri_6_factory(), tri_6_factory()});
        entry.initialized = true;
      }

      return &entry;
    }

    static TopologyMapEntry *wedge_15_factory()
    {
      static TopologyMapEntry entry(Ioss::Wedge15::name);

      if (!entry.initialized) {
        entry.set_valid_spatial_dimensions({false, false, false, true});
        entry.set_side_topologies({quad_8_factory(), quad_8_factory(), quad_8_factory(),
                                   tri_6_factory(), tri_6_factory()});
        entry.initialized = true;
      }

      return &entry;
    }

    static TopologyMapEntry *wedge_18_factory()
    {
      static TopologyMapEntry entry(Ioss::Wedge18::name);

      if (!entry.initialized) {
        entry.set_valid_spatial_dimensions({false, false, false, true});
        entry.set_side_topologies({quad_9_factory(), quad_9_factory(), quad_9_factory(),
                                   tri_6_factory(), tri_6_factory()});
        entry.initialized = true;
      }

      return &entry;
    }

    //***************************************************************************
    // topology::HEXAHEDRON -- topology::ELEMENT_RANK
    // defined on spatial dimension 3d
    // 8, 20 or 27 nodes nodes with 6 sides
    //***************************************************************************
    static TopologyMapEntry *hex_8_factory()
    {
      static TopologyMapEntry entry(Ioss::Hex8::name);

      if (!entry.initialized) {
        entry.set_valid_spatial_dimensions({false, false, false, true});
        entry.set_side_topologies({quad_4_factory(), quad_4_factory(), quad_4_factory(),
                                   quad_4_factory(), quad_4_factory(), quad_4_factory()});
        entry.initialized = true;
      }

      return &entry;
    }

    static TopologyMapEntry *hex_20_factory()
    {
      static TopologyMapEntry entry(Ioss::Hex20::name);

      if (!entry.initialized) {
        entry.set_valid_spatial_dimensions({false, false, false, true});
        entry.set_side_topologies({quad_8_factory(), quad_8_factory(), quad_8_factory(),
                                   quad_8_factory(), quad_8_factory(), quad_8_factory()});
        entry.initialized = true;
      }

      return &entry;
    }

    static TopologyMapEntry *hex_27_factory()
    {
      static TopologyMapEntry entry(Ioss::Hex27::name);

      if (!entry.initialized) {
        entry.set_valid_spatial_dimensions({false, false, false, true});
        entry.set_side_topologies({quad_9_factory(), quad_9_factory(), quad_9_factory(),
                                   quad_9_factory(), quad_9_factory(), quad_9_factory()});
        entry.initialized = true;
      }

      return &entry;
    }

  private:
    bool equivalent_valid_spatial_dimensions(const DimensionArray &validSpatialDimensions_) const
    {
      return validSpatialDimensions[0] == validSpatialDimensions_[0] &&
             validSpatialDimensions[1] == validSpatialDimensions_[1] &&
             validSpatialDimensions[2] == validSpatialDimensions_[2] &&
             validSpatialDimensions[3] == validSpatialDimensions_[3];
    }

    void set_valid_spatial_dimensions(const DimensionArray &validSpatialDimensions_)
    {
      validSpatialDimensions[0] = validSpatialDimensions_[0];
      validSpatialDimensions[1] = validSpatialDimensions_[1];
      validSpatialDimensions[2] = validSpatialDimensions_[2];
      validSpatialDimensions[3] = validSpatialDimensions_[3];
    }

    // Create and set TopologyMapEntry based side topologies since some of the implementation
    // details of the class cannot be automatically constructed from the Ioss_ElementTopology
    // object .. this will not be necessary if/when they are migrated to Ioss_ElementTopology
    void set_side_topologies(const std::vector<TopologyMapEntry *> &sideTopologies_)
    {
      int numSides = sideTopologies_.size();

      for (int side = 1; side <= numSides; side++) {
        if (topology->boundary_type(side) != sideTopologies_[side - 1]->topology) {
          std::ostringstream errmsg;
          fmt::print(errmsg,
                     "ERROR: For element topology: {} on side: {}, expected topology: {} does not "
                     "match topology: {}",
                     topology->name(), side, topology->boundary_type(side)->name(),
                     sideTopologies_[side - 1]->topology->name());
          IOSS_ERROR(errmsg);
        }
      }

      sideTopologies = sideTopologies_;
    }

    unsigned num_permutation_nodes() const
    {
      return topology->permutation()->num_permutation_nodes();
    }

    unsigned int           id{0};
    Ioss::ElementTopology *topology = nullptr;

    std::vector<TopologyMapEntry *> sideTopologies{};

    // Defines what spatial dimension the topology is valid on
    DimensionArray validSpatialDimensions;

    bool initialized{false};
  };

  inline std::ostream &operator<<(std::ostream &out, const TopologyMapEntry &t)
  {
    return out << t.name();
  }

  class IOTM_EXPORT IossTopologyMapping : public text_mesh::TopologyMapping<TopologyMapEntry>
  {
  public:
    TopologyMapEntry invalid_topology() const override { return TopologyMapEntry(); }

    // clang-format off
    void initialize_topology_map() override
    {
      m_nameToTopology = {
          {"NODE",         *TopologyMapEntry::node_factory()},
          {"LINE_2",       *TopologyMapEntry::line_2_factory()},
          {"LINE_3",       *TopologyMapEntry::line_3_factory()},
          {"TRI_3",        *TopologyMapEntry::tri_3_factory()},
          {"TRI_4",        *TopologyMapEntry::tri_4_factory()},
          {"TRI_6",        *TopologyMapEntry::tri_6_factory()},
          {"QUAD_4",       *TopologyMapEntry::quad_4_factory()},
          {"QUAD_6",       *TopologyMapEntry::quad_6_factory()},
          {"QUAD_8",       *TopologyMapEntry::quad_8_factory()},
          {"QUAD_9",       *TopologyMapEntry::quad_9_factory()},
          {"PARTICLE",     *TopologyMapEntry::particle_factory()},
          {"LINE_2_1D",    *TopologyMapEntry::line_2_1d_factory()},
          {"LINE_3_1D",    *TopologyMapEntry::line_3_1d_factory()},
          {"BEAM_2",       *TopologyMapEntry::beam_2_factory()},
          {"BEAM_3",       *TopologyMapEntry::beam_3_factory()},
          {"SHELL_LINE_2", *TopologyMapEntry::shell_line_2_factory()},
          {"SHELL_LINE_3", *TopologyMapEntry::shell_line_3_factory()},
          {"SPRING_2",     *TopologyMapEntry::spring_2_factory()},
          {"SPRING_3",     *TopologyMapEntry::spring_3_factory()},
          {"TRI_3_2D",     *TopologyMapEntry::tri_3_2d_factory()},
          {"TRI_4_2D",     *TopologyMapEntry::tri_4_2d_factory()},
          {"TRI_6_2D",     *TopologyMapEntry::tri_6_2d_factory()},
          {"QUAD_4_2D",    *TopologyMapEntry::quad_4_2d_factory()},
          {"QUAD_8_2D",    *TopologyMapEntry::quad_8_2d_factory()},
          {"QUAD_9_2D",    *TopologyMapEntry::quad_9_2d_factory()},
          {"SHELL_TRI_3",  *TopologyMapEntry::shell_tri_3_factory()},
          {"SHELL_TRI_4",  *TopologyMapEntry::shell_tri_4_factory()},
          {"SHELL_TRI_6",  *TopologyMapEntry::shell_tri_6_factory()},
          {"SHELL_QUAD_4", *TopologyMapEntry::shell_quad_4_factory()},
          {"SHELL_QUAD_8", *TopologyMapEntry::shell_quad_8_factory()},
          {"SHELL_QUAD_9", *TopologyMapEntry::shell_quad_9_factory()},
          {"TET_4",        *TopologyMapEntry::tet_4_factory()},
          {"TET_8",        *TopologyMapEntry::tet_8_factory()},
          {"TET_10",       *TopologyMapEntry::tet_10_factory()},
          {"TET_11",       *TopologyMapEntry::tet_11_factory()},
          {"PYRAMID_5",    *TopologyMapEntry::pyramid_5_factory()},
          {"PYRAMID_13",   *TopologyMapEntry::pyramid_13_factory()},
          {"PYRAMID_14",   *TopologyMapEntry::pyramid_14_factory()},
          {"WEDGE_6",      *TopologyMapEntry::wedge_6_factory()},
          {"WEDGE_12",     *TopologyMapEntry::wedge_12_factory()},
          {"WEDGE_15",     *TopologyMapEntry::wedge_15_factory()},
          {"WEDGE_18",     *TopologyMapEntry::wedge_18_factory()},
          {"HEX_8",        *TopologyMapEntry::hex_8_factory()},
          {"HEX_20",       *TopologyMapEntry::hex_20_factory()},
          {"HEX_27",       *TopologyMapEntry::hex_27_factory()}
      };
    }
    // clang-format on
  };
} // namespace Iotm
