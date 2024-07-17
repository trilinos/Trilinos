// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PANZER_STK_SETUP_UTILITIES_IMPL_HPP
#define PANZER_STK_SETUP_UTILITIES_IMPL_HPP

namespace panzer_stk {
namespace workset_utils {

/////////// TO BE DEPRECATED....
template<typename ArrayT>
void getIdsAndVertices(const panzer_stk::STK_Interface& mesh,
			 std::string blockId,
			 std::vector<std::size_t>& localIds,
			 ArrayT & vertices) {
  
  std::vector<stk::mesh::Entity> elements;
  mesh.getMyElements(blockId,elements);
  
  // loop over elements of this block
  for(std::size_t elm=0;elm<elements.size();++elm) {
    stk::mesh::Entity element = elements[elm];
    
    localIds.push_back(mesh.elementLocalId(element));
  }

  // get vertices (this is slightly faster then the local id version)
  mesh.getElementVertices(elements,blockId,vertices);
}
///////// END TO BE DEPRECATED

template<typename ArrayT>
void getIdsAndNodes(const panzer_stk::STK_Interface& mesh,
			 std::string blockId,
			 std::vector<std::size_t>& localIds,
			 ArrayT & nodes) {
  
  std::vector<stk::mesh::Entity> elements;
  mesh.getMyElements(blockId,elements);
  
  // loop over elements of this block
  for(std::size_t elm=0;elm<elements.size();++elm) {
    stk::mesh::Entity element = elements[elm];
    
    localIds.push_back(mesh.elementLocalId(element));
  }

  // get nodes (this is slightly faster then the local id version)
  mesh.getElementNodes(elements,blockId,nodes);
}
}
}

#endif
