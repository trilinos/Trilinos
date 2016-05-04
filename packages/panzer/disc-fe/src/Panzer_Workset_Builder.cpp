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

#include "PanzerDiscFE_config.hpp"

#include "Panzer_Workset_Builder_decl.hpp"

#ifdef HAVE_PANZER_EXPLICIT_INSTANTIATION

#include "Panzer_ExplicitTemplateInstantiation.hpp"
#include "Panzer_Workset_Builder_impl.hpp"

template
Teuchos::RCP<std::vector<panzer::Workset> > 
panzer::buildWorksets(const panzer::PhysicsBlock & pb,
		      const std::vector<std::size_t>& local_cell_ids,
		      const Kokkos::DynRankView<double,PHX::Device>& vertex_coordinates);

template
Teuchos::RCP<std::vector<panzer::Workset> > 
panzer::buildWorksets(const WorksetNeeds & needs,
                      const std::string & elementBlock,
		      const std::vector<std::size_t>& local_cell_ids,
		      const Kokkos::DynRankView<double,PHX::Device>& vertex_coordinates);

template
Teuchos::RCP<std::map<unsigned,panzer::Workset> >
panzer::buildBCWorkset(const panzer::PhysicsBlock & volume_pb,
		       const std::vector<std::size_t>& local_cell_ids,
		       const std::vector<std::size_t>& local_side_ids,
		       const Kokkos::DynRankView<double,PHX::Device>& vertex_coordinates);

template
Teuchos::RCP<std::map<unsigned,panzer::Workset> >
panzer::buildBCWorkset(const WorksetNeeds& needs,
                       const std::string& elementBlock,
                       const std::vector<std::size_t>& local_cell_ids,
                       const std::vector<std::size_t>& local_side_ids,
                       const Kokkos::DynRankView<double,PHX::Device>& vertex_coordinates,
                       const bool populate_value_arrays);

template
Teuchos::RCP<std::vector<panzer::Workset> > 
panzer::buildEdgeWorksets(const panzer::PhysicsBlock &,
	  	          const std::vector<std::size_t>&,
		          const std::vector<std::size_t>&,
		          const Kokkos::DynRankView<double,PHX::Device>&,
                          const panzer::PhysicsBlock &,
		          const std::vector<std::size_t>&,
		          const std::vector<std::size_t>&,
		          const Kokkos::DynRankView<double,PHX::Device>&);

template
Teuchos::RCP<std::map<unsigned,panzer::Workset> >
panzer::buildBCWorkset(const panzer::PhysicsBlock& pb_a,
                       const std::vector<std::size_t>& local_cell_ids_a,
                       const std::vector<std::size_t>& local_side_ids_a,
                       const Kokkos::DynRankView<double,PHX::Device>& vertex_coordinates_a,
                       const panzer::PhysicsBlock& pb_b,
                       const std::vector<std::size_t>& local_cell_ids_b,
                       const std::vector<std::size_t>& local_side_ids_b,
                       const Kokkos::DynRankView<double,PHX::Device>& vertex_coordinates_b);

#endif

namespace panzer {

void populateValueArrays(std::size_t num_cells,bool isSide,const WorksetNeeds & needs,
                         WorksetDetails & details,const Teuchos::RCP<WorksetDetails> other_details)
{
  using Teuchos::RCP;
  using Teuchos::rcp;

  panzer::MDFieldArrayFactory mdArrayFactory("",true);

  // setup the integration rules and bases
      
  std::vector<RCP<const panzer::PureBasis> > bases;
  std::vector<RCP<const panzer::IntegrationRule> > int_rules;
  if(isSide) {
    const panzer::CellData side_cell_data(num_cells,
                                          details.subcell_index,
                                          needs.cellData.getCellTopology());

    for(std::size_t i=0;i<needs.int_rules.size();i++) 
      int_rules.push_back(rcp(new IntegrationRule(needs.int_rules[i]->cubature_degree,side_cell_data)));

    for(std::size_t i=0;i<needs.bases.size();i++) 
      bases.push_back(rcp(new PureBasis(needs.bases[i]->type(),needs.bases[i]->order(),side_cell_data)));
  }
  else {
    int_rules = needs.int_rules;
    bases     = needs.bases;
  }

  details.ir_degrees = rcp(new std::vector<int>(0));
  details.basis_names = rcp(new std::vector<std::string>(0));

  for(std::size_t i=0;i<int_rules.size();i++) {

    details.ir_degrees->push_back(int_rules[i]->cubature_degree);
      
    RCP<panzer::IntegrationValues2<double> > iv2 = 
        rcp(new panzer::IntegrationValues2<double>("",true));
    iv2->setupArrays(int_rules[i]);
    if (Teuchos::nonnull(other_details))
      iv2->evaluateValues(details.cell_vertex_coordinates, other_details->int_rules[i]->ip_coordinates);
    else
      iv2->evaluateValues(details.cell_vertex_coordinates);
      
    details.int_rules.push_back(iv2);
      
    // Need to create all combinations of basis/ir pairings 
    for(std::size_t b=0;b<bases.size();b++) {
      RCP<panzer::BasisIRLayout> b_layout = rcp(new panzer::BasisIRLayout(bases[b],*int_rules[i]));
      details.basis_names->push_back(b_layout->name());

      std::size_t int_degree_index = std::distance(details.ir_degrees->begin(), 
                                                   std::find(details.ir_degrees->begin(), 
                                                             details.ir_degrees->end(), 
				                             int_rules[i]->order()));
      RCP<panzer::BasisValues2<double> > bv2 = 
          rcp(new panzer::BasisValues2<double>("",true,true));
      bv2->setupArrays(b_layout);
      bv2->evaluateValues(details.int_rules[int_degree_index]->cub_points,
                         details.int_rules[int_degree_index]->jac,
                         details.int_rules[int_degree_index]->jac_det,
                         details.int_rules[int_degree_index]->jac_inv,
                         details.int_rules[int_degree_index]->weighted_measure,
                         details.cell_vertex_coordinates);

      details.bases.push_back(bv2);
    }
  }

}

void populateValueArrays(std::size_t num_cells,bool isSide,const panzer::PhysicsBlock & pb,
                         WorksetDetails & details,const Teuchos::RCP<WorksetDetails> other_details)
{
  using Teuchos::RCP;

  WorksetNeeds needs;
  needs.cellData = pb.cellData();

  const std::map<int,RCP<panzer::IntegrationRule> >& int_rules = pb.getIntegrationRules();
  for(std::map<int,RCP<panzer::IntegrationRule> >::const_iterator ir_itr = int_rules.begin();
      ir_itr != int_rules.end(); ++ir_itr)
    needs.int_rules.push_back(ir_itr->second);
  
 const std::map<std::string,Teuchos::RCP<panzer::PureBasis> >& bases= pb.getBases();
 for(std::map<std::string,Teuchos::RCP<panzer::PureBasis> >::const_iterator b_itr = bases.begin();
     b_itr != bases.end(); ++b_itr)
   needs.bases.push_back(b_itr->second);
 
 return populateValueArrays(num_cells,isSide,needs,details,other_details);
}

}
