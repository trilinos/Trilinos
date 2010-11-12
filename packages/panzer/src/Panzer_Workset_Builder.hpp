
#ifndef CHARON_WORKSET_BUILDER_HPP
#define CHARON_WORKSET_BUILDER_HPP

#include <vector>
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

  template<typename ArrayT>
  Teuchos::RCP<std::vector<panzer::Workset> > 
  buildWorksets(const std::string& block_id,
		const std::vector<std::size_t>& local_cell_ids,
		const ArrayT& vertex_coordinates, 
		const panzer::InputPhysicsBlock& ipb,
		std::size_t workset_size,
		int base_cell_dimension);
  
//   Teuchos::RCP<std::map<unsigned,panzer::Workset> > 
//     buildBCWorkset(const SBC_Set* sideset,
// 		   const panzer::BoundaryCondition& bc,
// 		   Node_Vector_Index* node_coordinates,
// 		   const panzer::InputPhysicsBlock& ipb,
// 		   const shards::CellTopology& base_cell_topology);

}

#include "Panzer_Workset_BuilderT.hpp"

#endif
