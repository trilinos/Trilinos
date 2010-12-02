
#ifndef CHARON_WORKSET_BUILDER_HPP
#define CHARON_WORKSET_BUILDER_HPP

#include <vector>
#include <map>
#include "Teuchos_RCP.hpp"

namespace shards {
  class CellTopology;
}

namespace panzer {
  
  struct Workset;
  class MeshData;
  class BoundaryCondition;
  class InputPhysicsBlock;
  class PhysicsBlock;
  class BC;

  template<typename ArrayT>
  Teuchos::RCP<std::vector<panzer::Workset> > 
  buildWorksets(const std::string& block_id,
		const std::vector<std::size_t>& local_cell_ids,
		const ArrayT& vertex_coordinates, 
		const panzer::InputPhysicsBlock& ipb,
		std::size_t workset_size,
		int base_cell_dimension);
  
  template<typename ArrayT>
  Teuchos::RCP<std::map<unsigned,panzer::Workset> >
  buildBCWorkset(const panzer::BC& bc,
		 const std::vector<std::size_t>& local_cell_ids,
		 const std::vector<std::size_t>& local_side_ids,
		 const ArrayT& vertex_coordinates, 
		 const panzer::InputPhysicsBlock& ipb,
		 unsigned base_cell_dim);

}

#include "Panzer_Workset_BuilderT.hpp"

#endif
