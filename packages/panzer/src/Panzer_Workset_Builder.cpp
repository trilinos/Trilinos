#include "Panzer_config.hpp"

#ifdef HAVE_PANZER_EXPLICIT_INSTANTIATION

#include "Panzer_ExplicitTemplateInstantiation.hpp"

#include "Panzer_Workset_Builder_decl.hpp"
#include "Panzer_Workset_Builder_impl.hpp"

template
Teuchos::RCP<std::vector<panzer::Workset> >
panzer::buildWorksets<Intrepid::FieldContainer<double> >(const std::string& block_id,
							 const Teuchos::RCP<const shards::CellTopology> & blockTopo,
							 const std::vector<std::size_t>& local_cell_ids,
							 const Intrepid::FieldContainer<double>& vertex_coordinates, 
							 const panzer::InputPhysicsBlock& ipb,
							 std::size_t workset_size,
							 int base_cell_dimension);

template
Teuchos::RCP<std::vector<panzer::Workset> > 
panzer::buildWorksets(const panzer::PhysicsBlock & physBlk,
		      const std::vector<std::size_t>& local_cell_ids,
		      const Intrepid::FieldContainer<double>& vertex_coordinates, 
		      std::size_t workset_size);

template
Teuchos::RCP<std::map<unsigned,panzer::Workset> >
panzer::buildBCWorkset(const panzer::BC& bc,
		       const Teuchos::RCP<const shards::CellTopology> & blockTopo,
		       const std::vector<std::size_t>& local_cell_ids,
		       const std::vector<std::size_t>& local_side_ids,
		       const Intrepid::FieldContainer<double>& vertex_coordinates, 
		       const panzer::InputPhysicsBlock& ipb,
		       unsigned base_cell_dim);

#endif
