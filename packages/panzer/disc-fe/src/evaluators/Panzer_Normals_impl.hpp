// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PANZER_NORMALS_IMPL_HPP
#define PANZER_NORMALS_IMPL_HPP

#include <algorithm>
#include "Panzer_IntegrationRule.hpp"
#include "Panzer_Workset_Utilities.hpp"
#include "Intrepid2_FunctionSpaceTools.hpp"
#include "Intrepid2_CellTools.hpp"
#include "Phalanx_Scratch_Utilities.hpp"

namespace panzer {

//**********************************************************************
template<typename EvalT, typename Traits>
Normals<EvalT, Traits>::
Normals(
  const Teuchos::ParameterList& p)
   : normalize(true)
{
  // Read from parameters
  const std::string name = p.get<std::string>("Name");
  side_id = p.get<int>("Side ID");
  Teuchos::RCP<panzer::IntegrationRule> quadRule
     = p.get< Teuchos::RCP<panzer::IntegrationRule> >("IR");
  if(p.isParameter("Normalize")) // set default
     normalize = p.get<bool>("Normalize");

  // grab information from quadrature rule
  Teuchos::RCP<PHX::DataLayout> vector_dl = quadRule->dl_vector;
  quad_order = quadRule->cubature_degree;

  // build field, set as evaluated type
  normals = PHX::MDField<ScalarT,Cell,Point,Dim>(name, vector_dl);
  this->addEvaluatedField(normals);
  
  std::string n = "Normals: " + name;
  this->setName(n);
}

//**********************************************************************
template<typename EvalT, typename Traits>
void
Normals<EvalT, Traits>::
postRegistrationSetup(
  typename Traits::SetupData sd,
  PHX::FieldManager<Traits>& /* fm */)
{
  num_qp  = normals.extent(1);
  num_dim = normals.extent(2);
  
  quad_index =  panzer::getIntegrationRuleIndex(quad_order,(*sd.worksets_)[0], this->wda);
}

//**********************************************************************
template<typename EvalT, typename Traits>
void
Normals<EvalT, Traits>::
evaluateFields(
  typename Traits::EvalData workset)
{ 
  // ECC Fix: Get Physical Side Normals

  if(workset.num_cells>0) {
    Intrepid2::CellTools<PHX::exec_space>::getPhysicalSideNormals(normals.get_view(),
                                                                  this->wda(workset).int_rules[quad_index]->jac.get_view(),
                                                                  side_id, *this->wda(workset).int_rules[quad_index]->int_rule->topology);
      
    if(normalize) {
      // normalize vector: getPhysicalSideNormals does not 
      // return normalized vectors
      auto local_normals = normals;
      auto local_num_qp = num_qp;
      auto local_num_dim = num_dim;      

      if (Sacado::IsADType<ScalarT>::value) {
	scratch = PHX::View<ScalarT*>("normals:scratch",workset.num_cells,PHX::getFadSize(local_normals.get_static_view()));
      }
      else
	scratch = PHX::View<ScalarT*>("normals:scratch",workset.num_cells);
  
      auto norm = scratch;

      Kokkos::parallel_for("normalize", workset.num_cells, KOKKOS_LAMBDA (index_t c) {
        for(std::size_t q=0;q<local_num_qp;q++) {
	  norm(c) = 0.0;
   
          // compute squared norm
          for(std::size_t d=0;d<local_num_dim;d++)
            norm(c) += local_normals(c,q,d)*local_normals(c,q,d);
    
          // adjust for length of vector, now unit vectors
          norm(c) = sqrt(norm(c));
          for(std::size_t d=0;d<local_num_dim;d++)
            local_normals(c,q,d) /= norm(c);
        }
      });
    }
    // else normals correspond to differential
  }
}

//**********************************************************************

}

#endif
