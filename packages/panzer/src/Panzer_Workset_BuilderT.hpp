
#include <iostream>
#include <vector>
#include <map>
#include <algorithm>
#include "Panzer_Workset.hpp"
#include "Panzer_CellData.hpp"
//#include "Charon_BoundaryCondition.h"
#include "Panzer_InputPhysicsBlock.hpp"
//#include "Charon_Shards_Utilities.h"

#include "Phalanx_DataLayout_MDALayout.hpp"

// Intrepid
#include "Shards_CellTopology.hpp"
#include "Intrepid_DefaultCubatureFactory.hpp"
#include "Intrepid_CellTools.hpp"
#include "Intrepid_FunctionSpaceTools.hpp"
#include "Intrepid_Basis.hpp"
#include "Intrepid_HGRAD_QUAD_C1_FEM.hpp"
#include "Intrepid_HGRAD_QUAD_C2_FEM.hpp"
#include "Intrepid_HGRAD_HEX_C1_FEM.hpp"
#include "Intrepid_HGRAD_HEX_C2_FEM.hpp"


using namespace std;
using namespace Teuchos;

template<typename ArrayT>
Teuchos::RCP< std::vector<panzer::Workset> > 
panzer::buildWorksets(const std::string& block_id,
		      const std::vector<std::size_t>& local_cell_ids,
		      const ArrayT& vertex_coordinates, 
		      const panzer::InputPhysicsBlock& ipb,
		      std::size_t workset_size,
		      int base_cell_dimension)
{

  std::size_t total_num_cells = vertex_coordinates.dimension(0);

  Teuchos::RCP< std::vector<panzer::Workset> > worksets_ptr = 
    Teuchos::rcp(new std::vector<panzer::Workset>);
  std::vector<panzer::Workset>& worksets = *worksets_ptr;
  {
    std::size_t num_worksets = total_num_cells / workset_size;
    bool last_set_is_full = true;
    std::size_t last_workset_size = total_num_cells % workset_size;
    if (last_workset_size != 0) {
      num_worksets += 1;
      last_set_is_full = false;
    }    

    worksets.resize(num_worksets);
    std::vector<panzer::Workset>::iterator i;
    for (i = worksets.begin(); i != worksets.end(); ++i)
      i->num_cells = workset_size;
	 
    if (!last_set_is_full)
      worksets.back().num_cells = last_workset_size;

  }

  // assign workset cell local ids
  std::vector<std::size_t>::const_iterator local_begin = local_cell_ids.begin();
  for (std::vector<panzer::Workset>::iterator wkst = worksets.begin(); wkst != worksets.end(); ++wkst) {
    std::vector<std::size_t>::const_iterator begin_iter = local_begin;
    std::vector<std::size_t>::const_iterator end_iter = begin_iter + wkst->num_cells;
    local_begin = end_iter;
    wkst->cell_local_ids.assign(begin_iter,end_iter);
  }
  
  TEUCHOS_ASSERT(local_begin == local_cell_ids.end());

  // Copy cell vertex coordinates into local workset arrays
  std::size_t offset = 0;
  for (std::vector<panzer::Workset>::iterator wkst = worksets.begin(); wkst != worksets.end(); ++wkst) {
    for (std::size_t cell = 0; cell < wkst->num_cells; ++cell)
      for (std::size_t vertex = 0; vertex < Teuchos::as<std::size_t>(vertex_coordinates.dimension(1)); ++ vertex)
	for (std::size_t dim = 0; dim < Teuchos::as<std::size_t>(vertex_coordinates.dimension(2)); ++ dim)
	  wkst->cell_vertex_coordinates(cell,vertex,dim) = vertex_coordinates(cell + offset,vertex,dim);

    offset += wkst->num_cells;
  }

  TEUCHOS_ASSERT(offset == Teuchos::as<std::size_t>(vertex_coordinates.dimension(0)));


  // Setup integration rules and basis
  RCP< vector<int> > rcp_ir_degrees = rcp(new vector<int>(0));
  RCP< vector<string> > rcp_basis_names = rcp(new vector<string>(0));

  vector<int>& ir_degrees = *rcp_ir_degrees;
  vector<string>& basis_names = *rcp_basis_names;

  std::map<std::string,int> basis_to_int_order;

  for (std::size_t eq = 0; eq < ipb.eq_sets.size(); ++eq) {
    ir_degrees.push_back(ipb.eq_sets[eq].integration_order);
    basis_names.push_back(ipb.eq_sets[eq].basis);
    basis_to_int_order[ipb.eq_sets[eq].basis] = 
      ipb.eq_sets[eq].integration_order;
  }
  std::sort(ir_degrees.begin(), ir_degrees.end());
  std::unique(ir_degrees.begin(), ir_degrees.end());
  std::sort(basis_names.begin(), basis_names.end());
  std::unique(basis_names.begin(),basis_names.end());

  for (std::size_t wkst = 0; wkst < worksets.size(); ++wkst)
    worksets[wkst].int_rules.resize(ir_degrees.size());

  std::map<int,RCP<panzer::IntegrationRule> > degree_to_int_rule;
  for (std::size_t i = 0; i < ir_degrees.size(); ++i) {
    
    const panzer::CellData volume_cell_data(workset_size, base_cell_dimension);
      
    RCP<panzer::IntegrationRule> ir = 
      rcp(new panzer::IntegrationRule(ir_degrees[i], volume_cell_data));
    
    degree_to_int_rule[ir_degrees[i]] = ir;

    for (std::size_t wkst = 0; wkst < worksets.size(); ++wkst) {
      
      worksets[wkst].ir_degrees = rcp_ir_degrees;

      worksets[wkst].int_rules[i] = 
	rcp(new panzer::IntegrationValues<double,Intrepid::FieldContainer<double> >);
      
      worksets[wkst].int_rules[i]->setupArrays(ir);

      worksets[wkst].int_rules[i]->
	evaluateValues(vertex_coordinates);
    }
  }

  for (std::size_t wkst = 0; wkst < worksets.size(); ++wkst)
    worksets[wkst].bases.resize(basis_names.size());

  for (std::size_t i = 0; i < basis_names.size(); ++i) {
    
    RCP<panzer::Basis> cb = 
      rcp(new panzer::Basis(basis_names[i], *(degree_to_int_rule[basis_to_int_order[basis_names[i]]])));
    
    for (std::size_t wkst = 0; wkst < worksets.size(); ++wkst) {

      worksets[wkst].basis_names = rcp_basis_names;

      worksets[wkst].bases[i] = 
	rcp(new panzer::BasisValues<double,Intrepid::FieldContainer<double> >);
      
      worksets[wkst].bases[i]->setupArrays(cb);

      std::size_t int_degree_index = 
	std::distance(ir_degrees.begin(), 
		      std::find(ir_degrees.begin(), 
				ir_degrees.end(), 
				basis_to_int_order[basis_names[i]]));

      worksets[wkst].bases[i]->evaluateValues(worksets[wkst].int_rules[int_degree_index]->cub_points,
					      worksets[wkst].int_rules[int_degree_index]->jac_inv,
					      worksets[wkst].int_rules[int_degree_index]->weighted_measure);
    }
  }

  return worksets_ptr;
}

// ****************************************************************
// ****************************************************************
/*
Teuchos::RCP<std::map<unsigned,charon::Workset> > 
charon::buildBCWorkset(const SBC_Set* sideset,
		       const charon::BoundaryCondition& bc,
		       Node_Vector_Index* node_coordinates,
		       const charon::InputPhysicsBlock& ipb,
		       const shards::CellTopology& base_cell_topology)
{
  // All elements of boundary condition should go into one workset.
  // However due to design of Intrepid (requires same basis for all
  // cells), we have to separate the workset based on the local side
  // index.
  
  // key is local face index, value is vector of element pointers
  std::map<unsigned,std::vector<SBC_Element*> > element_list;
  
  SBC_Set& nc_sideset = const_cast<SBC_Set&>(*sideset);
 
  for (SBC_Element_Ref side = nc_sideset.Head_Local_Element(); 
       side.notNull(); side = nc_sideset.Next_Local_Element()) {
    
    // Nevada sideset lists have all elements associated with the
    // sides in the sideset. Loop over elements in the side set and
    // add SBC_Elements to the list only if they are owned by the
    // corresponding element block (otherwise sideset contributions
    // would cancel each other out).
    
    if ( (*side)->Cell()->myBlockId() == bc.elementBlockID() &&
	 (*side)->Cell()->Processor_Ownership() == TopoEntity::OWNED) {
      
      std::vector<unsigned> element_lids;
      for (int local_node=0; local_node < (*side)->Cell()->Num_Nodes(); 
	   ++local_node) {
	element_lids.push_back((*side)->Cell()->Nodes(local_node)->globalIndex());
      }
      
      std::vector<unsigned> side_lids;
      for (int local_node=0; local_node < (*side)->Side()->Num_Nodes(); 
	   ++local_node) {
	side_lids.push_back((*side)->Side()->Nodes(local_node)->globalIndex());
      }

      unsigned face_index = 
	charon::getLocalSideIndexFromGlobalNodeList(element_lids, side_lids, 
						    base_cell_topology);
      
      element_list[face_index].push_back(*side);
      
    }	
    
  }

  std::cout << "BC:\n" << bc << std::endl;
  for (std::map<unsigned,std::vector<SBC_Element*> >::const_iterator side =  element_list.begin(); 
       side != element_list.end(); ++side) {
    TASK_0_cout << "  Workset for side " << side->first << " has " << side->second.size() << " elements." << std::endl;
  }

  // *****************
  // New BCs only work for Intrepid elements!!
  // *****************
  using std::vector;
  using std::map;
  using std::string;
  using Teuchos::RCP;
  using Teuchos::rcp;
   
  RCP< vector<int> > rcp_ir_degrees = rcp(new vector<int>(0));
  RCP< vector<string> > rcp_basis_names = rcp(new vector<string>(0));

  vector<int>& ir_degrees = *rcp_ir_degrees;
  vector<string>& basis_names = *rcp_basis_names;

  map<string,int> basis_to_int_order;

  for (std::size_t eq = 0; eq < ipb.eq_sets.size(); ++eq) {
    ir_degrees.push_back(ipb.eq_sets[eq].integration_order);
    basis_names.push_back(ipb.eq_sets[eq].basis);
    basis_to_int_order[ipb.eq_sets[eq].basis] = 
      ipb.eq_sets[eq].integration_order;
  }
  std::sort(ir_degrees.begin(), ir_degrees.end());
  std::unique(ir_degrees.begin(), ir_degrees.end());
  std::sort(basis_names.begin(), basis_names.end());
  std::unique(basis_names.begin(),basis_names.end());

  int dim = static_cast<int>(base_cell_topology.getDimension());

  // key is local face index, value is workset with all elements
  // for that local face
  Teuchos::RCP<std::map<unsigned,charon::Workset> > worksets_ptr = 
    Teuchos::rcp(new std::map<unsigned,charon::Workset>);

  std::map<unsigned,charon::Workset>& worksets = *worksets_ptr;

  // create worksets 
  for (std::map<unsigned,std::vector<SBC_Element*> >::iterator side =
	 element_list.begin(); side != element_list.end(); ++side) {
    THashList elements;
    for (std::vector<SBC_Element*>::iterator element = side->second.begin();
	 element != side->second.end(); ++element) {
      elements.insert( (*element)->Cell() );
    }
    worksets[side->first].sideset_elements = elements;
    worksets[side->first].begin = 
      worksets[side->first].sideset_elements.begin();
    worksets[side->first].end = 
      worksets[side->first].sideset_elements.end();
    worksets[side->first].num_cells = 
      worksets[side->first].sideset_elements.size();
    
    worksets[side->first].coordinate_index = node_coordinates;
  }

  // setup the integration rules
  for (std::map<unsigned,charon::Workset>::iterator wkst = worksets.begin(); 
       wkst != worksets.end(); ++wkst)
    wkst->second.int_rules.resize(ir_degrees.size());

  for (int i = 0; i < ir_degrees.size(); ++i) {
    
    for (std::map<unsigned,charon::Workset>::iterator wkst = worksets.begin();
	 wkst != worksets.end(); ++wkst) {
      
      const charon::CellData side_cell_data(wkst->second.num_cells, dim,
					    static_cast<int>(wkst->first));

      RCP<charon::IntegrationRule> ir = 
	rcp(new charon::IntegrationRule(ir_degrees[i], side_cell_data));
      
      wkst->second.ir_degrees = rcp_ir_degrees;

      wkst->second.int_rules[i] = 
	rcp(new charon::IntegrationValues<double,Intrepid::FieldContainer<double> >);
      
      wkst->second.int_rules[i]->setupArrays(ir);

      wkst->second.int_rules[i]->
	evaluateValues(wkst->second.begin, wkst->second.end, 
		       wkst->second.coordinate_index);
    }
  }



  // setup the basis functions
  for (std::map<unsigned,charon::Workset>::iterator wkst = worksets.begin(); 
       wkst != worksets.end(); ++wkst)
    wkst->second.bases.resize(basis_names.size());


  for (int i = 0; i < basis_names.size(); ++i) {
    
    for (std::map<unsigned,charon::Workset>::iterator wkst = worksets.begin(); 
	 wkst != worksets.end(); ++wkst) {
      
      std::size_t int_degree_index = 
	std::distance(ir_degrees.begin(), 
		      std::find(ir_degrees.begin(), 
				ir_degrees.end(), 
				basis_to_int_order[basis_names[i]]));


      RCP<charon::Basis> cb = 
	rcp(new charon::Basis(basis_names[i], *(wkst->second.int_rules[int_degree_index]->int_rule)));

      wkst->second.basis_names = rcp_basis_names;

      wkst->second.bases[i] = 
	rcp(new charon::BasisValues<double,Intrepid::FieldContainer<double> >);
      
      wkst->second.bases[i]->setupArrays(cb);

      wkst->second.bases[i]->evaluateValues(wkst->second.int_rules[int_degree_index]->cub_points,
					    wkst->second.int_rules[int_degree_index]->jac_inv,
					    wkst->second.int_rules[int_degree_index]->weighted_measure);

    }

  }

  return worksets_ptr;
}
*/
