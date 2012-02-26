#ifndef PANZER_STK_SETUP_UTILITIES_IMPL_HPP
#define PANZER_STK_SETUP_UTILITIES_IMPL_HPP

namespace panzer_stk {
namespace workset_utils {

template<typename ArrayT>
void getIdsAndVertices(const panzer_stk::STK_Interface& mesh,
			 std::string blockId,
			 std::vector<std::size_t>& localIds,
			 ArrayT & vertices) {
  
  std::vector<stk::mesh::Entity*> elements;
  mesh.getMyElements(blockId,elements);
  
  // loop over elements of this block
  for(std::size_t elm=0;elm<elements.size();++elm) {
    stk::mesh::Entity * element = elements[elm];
    
    localIds.push_back(mesh.elementLocalId(element));
  }

  // get vertices
  mesh.getElementVertices(localIds,vertices);
}

}
}

#endif
