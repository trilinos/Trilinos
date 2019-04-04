// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef percept_AdaptHelperFunctions_hpp
#define percept_AdaptHelperFunctions_hpp

#include "Teuchos_RCP.hpp"

#include <percept/PerceptMesh.hpp>
#include <adapt/UniformRefinerPattern.hpp>

namespace percept
{

inline
stk::mesh::Selector
make_active_part_selector(stk::mesh::MetaData & meta_data)
{
  stk::mesh::Selector adaptActive;
  stk::mesh::EntityRank part_ranks[] = {stk::topology::ELEMENT_RANK, meta_data.side_rank()};
  for (unsigned irank=0; irank < 2; irank++) {
    std::ostringstream active_part_name;
    active_part_name << "refine_active_elements_part_" << static_cast<unsigned int>(part_ranks[irank]);
    stk::mesh::Part* active_part = meta_data.get_part(active_part_name.str());
    if (!active_part)
      throw std::runtime_error("error - no active part can be found");
    adaptActive = adaptActive | *active_part;
  }

  return adaptActive;
}

inline
stk::mesh::Selector
make_inactive_part_selector(stk::mesh::MetaData & meta_data)
{
  stk::mesh::Selector adaptInactive;
  stk::mesh::EntityRank part_ranks[] = {stk::topology::ELEMENT_RANK, meta_data.side_rank()};
  for (unsigned irank=0; irank < 2; irank++) {
    std::ostringstream inactive_part_name;
    inactive_part_name << "refine_inactive_elements_part_" << static_cast<unsigned int>(part_ranks[irank]);
    stk::mesh::Part* inactive_part = meta_data.get_part(inactive_part_name.str());
    if (!inactive_part)
      throw std::runtime_error("error - no inactive part can be found");
    adaptInactive = adaptInactive | *inactive_part;
  }

  return adaptInactive;
}

inline
stk::mesh::Selector
make_active_part_selector(stk::mesh::MetaData & meta_data,
                          stk::mesh::Selector inputSelector)
{
  inputSelector &= make_active_part_selector(meta_data);

  return inputSelector;
}

inline
stk::mesh::Selector
make_inactive_part_selector(stk::mesh::MetaData & meta_data,
                          stk::mesh::Selector inputSelector)
{
  inputSelector &= make_inactive_part_selector(meta_data);

  return inputSelector;
}

inline
Teuchos::RCP<UniformRefinerPatternBase>
make_local_break_pattern(percept::PerceptMesh &eMesh, BlockNamesType block_names = BlockNamesType()) {
  Teuchos::RCP<UniformRefinerPatternBase> localBreakPattern = Teuchos::null;

  bool has_tri = false, has_quad = false, has_tet = false, has_wedge = false, has_hex = false, has_pyr = false;
  int ntype = 0;

  stk::mesh::PartVector pv = eMesh.get_fem_meta_data()->get_parts();
  for (unsigned ii=0; ii < pv.size(); ++ii) {
    if (stk::mesh::is_auto_declared_part(*pv[ii]))
      continue;
    stk::topology::topology_t topo = pv[ii]->topology();
    switch(topo)
      {
      case stk::topology::TRI_3_2D:
        if (!has_tri) ++ntype;
        has_tri = true;
        break;
      case stk::topology::QUAD_4_2D:
        if (!has_quad) ++ntype;
        has_quad = true;
        break;
      case stk::topology::TET_4:
        if (!has_tet) ++ntype;
        has_tet = true;
        break;
      case stk::topology::WEDGE_6:
        if (!has_wedge) ++ntype;
        has_wedge = true;
        break;
      case stk::topology::HEX_8:
        if (!has_hex) ++ntype;
        has_hex = true;
        break;
      case stk::topology::PYRAMID_5:
        if (!has_pyr) ++ntype;
        has_pyr = true;
        break;
      default:
        break;
      }
  }

  if (eMesh.get_spatial_dim() == 2)
    {
      if      (has_quad && !has_tri)
        localBreakPattern = Teuchos::rcp<UniformRefinerPatternBase>(new Local_Quad4_Quad4_N_Transition(eMesh, block_names));
      else if (!has_quad && has_tri)
        localBreakPattern = Teuchos::rcp<UniformRefinerPatternBase>(new Local_Tri3_Tri3_N_HangingNode(eMesh, block_names));
      else if (has_quad && has_tri)
        localBreakPattern = Teuchos::rcp<UniformRefinerPatternBase>(new Local_Quad4_Tri3_Hybrid_Transition(eMesh, block_names));
    }
  else
    {
      if (ntype == 1)
        {
          if (has_tet)
            localBreakPattern = Teuchos::rcp<UniformRefinerPatternBase>(new Local_Tet4_Tet4_N_HangingNode(eMesh, block_names));
          else if (has_wedge)
            localBreakPattern = Teuchos::rcp<UniformRefinerPatternBase>(new Local_Wedge6_Wedge6_N_Transition(eMesh, block_names));
          else if (has_hex)
            localBreakPattern = Teuchos::rcp<UniformRefinerPatternBase>(new Local_Hex8_Hex8_N_Transition(eMesh, block_names));
          else if (has_pyr)
            localBreakPattern = Teuchos::rcp<UniformRefinerPatternBase>(new Local_Pyr5_Pyr5_N_Transition(eMesh, block_names));
        }
      else
        {
          localBreakPattern = Teuchos::rcp<UniformRefinerPatternBase>(new Local_Hybrid_3D(eMesh, block_names));
        }
    }

  if (Teuchos::is_null(localBreakPattern))
    {
      std::string error_msg = "make_local_break_pattern: unsupported combination of element types:";
      if (has_quad) error_msg += " quad";
      if (has_tri) error_msg += " tri";
      if (has_tet) error_msg += " tet";
      if (has_wedge) error_msg += " wedge";
      if (has_hex) error_msg += " hex";
      if (has_pyr) error_msg += " pyr";
      error_msg += "\n";
      VERIFY_MSG(error_msg);
    }
  //std::cout << "localBreakPattern= " << eMesh.demangle(typeid(*localBreakPattern).name()) << "\n" << eMesh.demangled_stacktrace(20) << std::endl;
  return localBreakPattern;
}

}

#endif
