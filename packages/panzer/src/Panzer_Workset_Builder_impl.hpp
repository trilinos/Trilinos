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

#ifndef PANZER_WORKSET_BUILDER_IMPL_HPP
#define PANZER_WORKSET_BUILDER_IMPL_HPP

#include <iostream>
#include <vector>
#include <map>
#include <algorithm>
#include "Panzer_Workset.hpp"
#include "Panzer_CellData.hpp"
#include "Panzer_BC.hpp"
#include "Panzer_InputPhysicsBlock.hpp"
#include "Panzer_PhysicsBlock.hpp"
#include "Panzer_Shards_Utilities.hpp"
#include "Panzer_CommonArrayFactories.hpp"

#include "Phalanx_DataLayout_MDALayout.hpp"

// Intrepid
#include "Shards_CellTopology.hpp"
#include "Intrepid_DefaultCubatureFactory.hpp"
#include "Intrepid_CellTools.hpp"
#include "Intrepid_FunctionSpaceTools.hpp"
#include "Intrepid_Basis.hpp"

template<typename ArrayT>
Teuchos::RCP< std::vector<panzer::Workset> > 
panzer::buildWorksets(const panzer::PhysicsBlock& physBlk,
		      const std::vector<std::size_t>& local_cell_ids,
		      const ArrayT& vertex_coordinates, 
		      std::size_t workset_size)
{
   int base_cell_dimension = physBlk.getBaseCellTopology().getDimension();

   // for now a simple call out to old version
   return buildWorksets(physBlk.elementBlockID(),Teuchos::rcp(new shards::CellTopology(physBlk.getBaseCellTopology())),
                        local_cell_ids,vertex_coordinates,
                        physBlk.getInputPhysicsBlock(),workset_size,base_cell_dimension);
}

template<typename ArrayT>
Teuchos::RCP< std::vector<panzer::Workset> > 
panzer::buildWorksets(const std::string& block_id,
                      const Teuchos::RCP<const shards::CellTopology> & blockTopo,
		      const std::vector<std::size_t>& local_cell_ids,
		      const ArrayT& vertex_coordinates, 
		      const panzer::InputPhysicsBlock& ipb,
		      std::size_t workset_size,
		      int base_cell_dimension)
{
  using std::vector;
  using std::string;
  using Teuchos::RCP;
  using Teuchos::rcp;

  panzer::IntrepidFieldContainerFactory arrayFactory;

  // std::size_t total_num_cells = vertex_coordinates.dimension(0);
  std::size_t total_num_cells = local_cell_ids.size();

  Teuchos::RCP< std::vector<panzer::Workset> > worksets_ptr = 
    Teuchos::rcp(new std::vector<panzer::Workset>);
  std::vector<panzer::Workset>& worksets = *worksets_ptr;
   
  // special case for 0 elements!
  if(total_num_cells==0) {
     // Setup integration rules and basis
     RCP< vector<int> > rcp_ir_degrees = rcp(new vector<int>(0));
     RCP< vector<string> > rcp_basis_names = rcp(new vector<string>(0));
      
     vector<int>& ir_degrees = *rcp_ir_degrees;
     vector<string>& basis_names = *rcp_basis_names;
     std::map<std::string,int> basis_to_int_order;
     {
        std::set<int> ir_degrees_set;
        std::set<std::string> basis_names_set;
        for (std::size_t eq = 0; eq < ipb.eq_sets.size(); ++eq) {
      
          basis_to_int_order[ipb.eq_sets[eq].basis] = 
            ipb.eq_sets[eq].integration_order;
      
          ir_degrees_set.insert(ipb.eq_sets[eq].integration_order);
          basis_names_set.insert(ipb.eq_sets[eq].basis);
        }
        ir_degrees.insert(ir_degrees.begin(),ir_degrees_set.begin(),ir_degrees_set.end());
        basis_names.insert(basis_names.begin(),basis_names_set.begin(),basis_names_set.end());
     }

     worksets.resize(1);
     std::vector<panzer::Workset>::iterator i = worksets.begin();
     i->num_cells = 0;
     i->block_id = block_id;
     i->ir_degrees = rcp_ir_degrees;
     i->basis_names = rcp_basis_names;
     i->int_rules.resize(ir_degrees.size());

     std::map<int,RCP<panzer::IntegrationRule> > degree_to_int_rule;
     for (std::size_t j = 0; j < ir_degrees.size(); ++j) {
       const panzer::CellData volume_cell_data(workset_size, base_cell_dimension,blockTopo);
         
       RCP<panzer::IntegrationRule> ir = 
         rcp(new panzer::IntegrationRule(ir_degrees[j], volume_cell_data));

       degree_to_int_rule[ir_degrees[j]] = ir;
       
       i->int_rules[j] = 
   	  rcp(new panzer::IntegrationValues<double,Intrepid::FieldContainer<double> >);
       i->int_rules[j]->setupArrays(ir);
     }

     i->bases.resize(basis_names.size());

     for (std::size_t j = 0; j < basis_names.size(); ++j) {
    
       RCP<panzer::BasisIRLayout> cb = 
          rcp(new panzer::BasisIRLayout(basis_names[j], *(degree_to_int_rule[basis_to_int_order[basis_names[j]]])));
    
       i->bases[j] = 
  	  rcp(new panzer::BasisValues<double,Intrepid::FieldContainer<double> >);
       i->bases[j]->setupArrays(cb,arrayFactory);
     }
     
     return worksets_ptr;
  } // end special case

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
	 
    if (!last_set_is_full) {
      worksets.back().num_cells = last_workset_size;
    }
  }

  // assign workset cell local ids
  std::vector<std::size_t>::const_iterator local_begin = local_cell_ids.begin();
  for (std::vector<panzer::Workset>::iterator wkst = worksets.begin(); wkst != worksets.end(); ++wkst) {
    std::vector<std::size_t>::const_iterator begin_iter = local_begin;
    std::vector<std::size_t>::const_iterator end_iter = begin_iter + wkst->num_cells;
    local_begin = end_iter;
    wkst->cell_local_ids.assign(begin_iter,end_iter);
    wkst->cell_vertex_coordinates.resize(workset_size,
					 vertex_coordinates.dimension(1),
					 vertex_coordinates.dimension(2));
    wkst->block_id = block_id;
    wkst->subcell_dim = base_cell_dimension;
    wkst->subcell_index = 0;
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

  std::set<int> ir_degrees_set;
  std::set<std::string> basis_names_set;
  for (std::size_t eq = 0; eq < ipb.eq_sets.size(); ++eq) {

    basis_to_int_order[ipb.eq_sets[eq].basis] = 
      ipb.eq_sets[eq].integration_order;

    ir_degrees_set.insert(ipb.eq_sets[eq].integration_order);
    basis_names_set.insert(ipb.eq_sets[eq].basis);
  }
  ir_degrees.insert(ir_degrees.begin(),ir_degrees_set.begin(),ir_degrees_set.end());
  basis_names.insert(basis_names.begin(),basis_names_set.begin(),basis_names_set.end());

  for (std::size_t wkst = 0; wkst < worksets.size(); ++wkst)
    worksets[wkst].int_rules.resize(ir_degrees.size());

  std::map<int,RCP<panzer::IntegrationRule> > degree_to_int_rule;
  for (std::size_t i = 0; i < ir_degrees.size(); ++i) {
    
    const panzer::CellData volume_cell_data(workset_size, base_cell_dimension,blockTopo);
      
    RCP<panzer::IntegrationRule> ir = 
      rcp(new panzer::IntegrationRule(ir_degrees[i], volume_cell_data));
    
    degree_to_int_rule[ir_degrees[i]] = ir;

    for (std::size_t wkst = 0; wkst < worksets.size(); ++wkst) {
      
      worksets[wkst].ir_degrees = rcp_ir_degrees;

      worksets[wkst].int_rules[i] = 
	rcp(new panzer::IntegrationValues<double,Intrepid::FieldContainer<double> >);
      
      worksets[wkst].int_rules[i]->setupArrays(ir);

      worksets[wkst].int_rules[i]->
	evaluateValues(worksets[wkst].cell_vertex_coordinates);
    }
  }

  for (std::size_t wkst = 0; wkst < worksets.size(); ++wkst)
    worksets[wkst].bases.resize(basis_names.size());

  for (std::size_t i = 0; i < basis_names.size(); ++i) {
    
    RCP<panzer::BasisIRLayout> cb = 
      rcp(new panzer::BasisIRLayout(basis_names[i], *(degree_to_int_rule[basis_to_int_order[basis_names[i]]])));
    
    for (std::size_t wkst = 0; wkst < worksets.size(); ++wkst) {

      worksets[wkst].basis_names = rcp_basis_names;

      worksets[wkst].bases[i] = 
	rcp(new panzer::BasisValues<double,Intrepid::FieldContainer<double> >);
      
      worksets[wkst].bases[i]->setupArrays(cb,arrayFactory);

      std::size_t int_degree_index = 
	std::distance(ir_degrees.begin(), 
		      std::find(ir_degrees.begin(), 
				ir_degrees.end(), 
				basis_to_int_order[basis_names[i]]));

      worksets[wkst].bases[i]->evaluateValues(worksets[wkst].int_rules[int_degree_index]->cub_points,
					      worksets[wkst].int_rules[int_degree_index]->jac,
					      worksets[wkst].int_rules[int_degree_index]->jac_det,
					      worksets[wkst].int_rules[int_degree_index]->jac_inv,
					      worksets[wkst].int_rules[int_degree_index]->weighted_measure,
					      worksets[wkst].cell_vertex_coordinates);
    }
  }

  return worksets_ptr;
}

// ****************************************************************
// ****************************************************************

template<typename ArrayT>
Teuchos::RCP<std::map<unsigned,panzer::Workset> >
panzer::buildBCWorkset(const panzer::BC& bc,
                       const Teuchos::RCP<const shards::CellTopology> & blockTopo,
		       const std::vector<std::size_t>& local_cell_ids,
		       const std::vector<std::size_t>& local_side_ids,
		       const ArrayT& vertex_coordinates, 
		       const panzer::InputPhysicsBlock& ipb,
		       unsigned base_cell_dim)
{
  using std::vector;
  using std::string;
  using Teuchos::RCP;
  using Teuchos::rcp;

  panzer::IntrepidFieldContainerFactory arrayFactory;

  // key is local face index, value is workset with all elements
  // for that local face
  Teuchos::RCP<std::map<unsigned,panzer::Workset> > worksets_ptr = 
    Teuchos::rcp(new std::map<unsigned,panzer::Workset>);

  // All elements of boundary condition should go into one workset.
  // However due to design of Intrepid (requires same basis for all
  // cells), we have to separate the workset based on the local side
  // index.  Each workset for a boundary condition is associated with
  // a local side for the element
  
  // std::cout << "local_side_ids.size() = " << local_side_ids.size() 
  // 	    << std::endl;

  TEUCHOS_ASSERT(local_side_ids.size() == local_cell_ids.size());
  TEUCHOS_ASSERT(local_side_ids.size() == static_cast<std::size_t>(vertex_coordinates.dimension(0)));

  // key is local face index, value is a pair of cell index and vector of element local ids
  std::map<unsigned,std::vector<std::pair<std::size_t,std::size_t> > > element_list;
  for (std::size_t cell=0; cell < local_cell_ids.size(); ++cell)
    element_list[local_side_ids[cell]].push_back(std::make_pair(cell,local_cell_ids[cell])); 

  // *****************
  // New BCs only work for Intrepid elements!!
  // *****************
  using std::map;
  using std::string;
  using Teuchos::rcp;
   
  RCP< vector<int> > rcp_ir_degrees = rcp(new vector<int>(0));
  RCP< vector<string> > rcp_basis_names = rcp(new vector<string>(0));

  vector<int>& ir_degrees = *rcp_ir_degrees;
  vector<string>& basis_names = *rcp_basis_names;

  map<string,int> basis_to_int_order;

  std::set<int> ir_degrees_set;
  std::set<std::string> basis_names_set;
  for (std::size_t eq = 0; eq < ipb.eq_sets.size(); ++eq) {

    basis_to_int_order[ipb.eq_sets[eq].basis] = 
      ipb.eq_sets[eq].integration_order;

    ir_degrees_set.insert(ipb.eq_sets[eq].integration_order);
    basis_names_set.insert(ipb.eq_sets[eq].basis);
  }
  ir_degrees.insert(ir_degrees.begin(),ir_degrees_set.begin(),ir_degrees_set.end());
  basis_names.insert(basis_names.begin(),basis_names_set.begin(),basis_names_set.end());

  std::map<unsigned,panzer::Workset>& worksets = *worksets_ptr;

  // create worksets 
  std::map<unsigned,std::vector<std::pair<std::size_t,std::size_t> > >::const_iterator side;
  for (side = element_list.begin(); side != element_list.end(); ++side) {

    std::vector<std::size_t>& cell_local_ids = worksets[side->first].cell_local_ids;
    Intrepid::FieldContainer<double> & coords = worksets[side->first].cell_vertex_coordinates;
    coords.resize(side->second.size(),vertex_coordinates.dimension(1),vertex_coordinates.dimension(2));
    for (std::size_t cell = 0; cell < side->second.size(); ++cell) {
      cell_local_ids.push_back(side->second[cell].second);

      for (std::size_t vertex = 0; vertex < Teuchos::as<std::size_t>(vertex_coordinates.dimension(1)); ++ vertex)
	for (std::size_t dim = 0; dim < Teuchos::as<std::size_t>(vertex_coordinates.dimension(2)); ++ dim)
	  coords(cell,vertex,dim) = vertex_coordinates(side->second[cell].first,vertex,dim);
    }
    worksets[side->first].num_cells = worksets[side->first].cell_local_ids.size();
    worksets[side->first].block_id = bc.elementBlockID();
    worksets[side->first].subcell_dim = base_cell_dim -1;
    worksets[side->first].subcell_index = side->first;
  }

  // setup the integration rules
  for (std::map<unsigned,panzer::Workset>::iterator wkst = worksets.begin(); 
       wkst != worksets.end(); ++wkst)
    wkst->second.int_rules.resize(ir_degrees.size());

  for (std::size_t i = 0; i < ir_degrees.size(); ++i) {
    
    for (std::map<unsigned,panzer::Workset>::iterator wkst = worksets.begin();
	 wkst != worksets.end(); ++wkst) {
      
      const panzer::CellData side_cell_data(wkst->second.num_cells, base_cell_dim,
					    static_cast<int>(wkst->first),blockTopo);

      RCP<panzer::IntegrationRule> ir = 
	rcp(new panzer::IntegrationRule(ir_degrees[i], side_cell_data));
      
      wkst->second.ir_degrees = rcp_ir_degrees;

      wkst->second.int_rules[i] = 
	rcp(new panzer::IntegrationValues<double,Intrepid::FieldContainer<double> >);
      
      wkst->second.int_rules[i]->setupArrays(ir);

      wkst->second.int_rules[i]->evaluateValues(wkst->second.cell_vertex_coordinates);
    }
  }


  // setup the basis functions
  for (std::map<unsigned,panzer::Workset>::iterator wkst = worksets.begin(); 
       wkst != worksets.end(); ++wkst)
    wkst->second.bases.resize(basis_names.size());

  for (std::size_t i = 0; i < basis_names.size(); ++i) {
    
    for (std::map<unsigned,panzer::Workset>::iterator wkst = worksets.begin(); 
	 wkst != worksets.end(); ++wkst) {
      
      std::size_t int_degree_index = 
	std::distance(ir_degrees.begin(), 
		      std::find(ir_degrees.begin(), 
				ir_degrees.end(), 
				basis_to_int_order[basis_names[i]]));

      
      RCP<panzer::BasisIRLayout> cb = 
	rcp(new panzer::BasisIRLayout(basis_names[i], *(wkst->second.int_rules[int_degree_index]->int_rule)));
      
      wkst->second.basis_names = rcp_basis_names;

      wkst->second.bases[i] = 
	rcp(new panzer::BasisValues<double,Intrepid::FieldContainer<double> >);
      
      wkst->second.bases[i]->setupArrays(cb,arrayFactory);

      wkst->second.bases[i]->evaluateValues(wkst->second.int_rules[int_degree_index]->cub_points,
					    wkst->second.int_rules[int_degree_index]->jac,
					    wkst->second.int_rules[int_degree_index]->jac_det,
					    wkst->second.int_rules[int_degree_index]->jac_inv,
					    wkst->second.int_rules[int_degree_index]->weighted_measure,
					    wkst->second.cell_vertex_coordinates);

    }

  }

  return worksets_ptr;
}

#endif
