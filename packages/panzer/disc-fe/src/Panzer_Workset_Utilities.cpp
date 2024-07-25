// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PANZER_WORKSET_UTILITIES_HPP
#define PANZER_WORKSET_UTILITIES_HPP

#include "Panzer_Traits.hpp"
#include "Panzer_Workset.hpp"
#include "Teuchos_Assert.hpp"
#include <vector>
#include <string>
#include <algorithm>
#include <iterator>

namespace panzer {

  std::vector<std::string>::size_type 
  getPureBasisIndex(std::string basis_name, const panzer::Workset& workset, WorksetDetailsAccessor& wda)
  {
    std::vector<std::string>::iterator basis = wda(workset).basis_names->begin();
    std::vector<std::string>::const_iterator last = wda(workset).basis_names->end();
    
    while (basis != last) {

      std::vector<std::string>::size_type index = std::distance(wda(workset).basis_names->begin(), basis);
      if (wda(workset).bases[index]->basis_layout->getBasis()->name() == basis_name)
	break;

      ++basis;
    }

    TEUCHOS_TEST_FOR_EXCEPTION(basis == wda(workset).basis_names->end(),
			       std::logic_error,
			       "Could not find the basis named \"" 
                               << basis_name << "\" in the workset!");

    return std::distance(wda(workset).basis_names->begin(), basis);
  }

  std::vector<std::string>::size_type 
  getBasisIndex(std::string basis_name, const panzer::Workset& workset, WorksetDetailsAccessor& wda)
  {
    std::vector<std::string>::iterator basis;

    basis = std::find(wda(workset).basis_names->begin(),
		      wda(workset).basis_names->end(),
		      basis_name);

    TEUCHOS_TEST_FOR_EXCEPTION(basis == wda(workset).basis_names->end(),
			       std::logic_error,
			       "Could not find the basis named \"" 
                               << basis_name << "\" in the workset!");

    return std::distance(wda(workset).basis_names->begin(), basis);
  }

  std::vector<int>::size_type
  getIntegrationRuleIndex(int ir_degree, const panzer::Workset& workset, WorksetDetailsAccessor& wda)
  {
    std::vector<int>::iterator ir;

    ir = std::find(wda(workset).ir_degrees->begin(),
		   wda(workset).ir_degrees->end(),
		   ir_degree);
    
    TEUCHOS_TEST_FOR_EXCEPTION(ir == wda(workset).ir_degrees->end(),
			       std::logic_error,
			       "Could not find the integration rule degree \"" 
                               << ir_degree << "\" in the workset!");

    return std::distance(wda(workset).ir_degrees->begin(), ir);
  }

  void printWorkset(std::ostream& os, const panzer::Workset & workset, WorksetDetailsAccessor& wda)
  {
     os << "WORKSET"
        << " block_id = \"" << wda(workset).block_id << "\""
        << " num_cells = \"" << workset.num_cells << "\"\n";
     os << "   cell_local_ids = [ "; 
     for(index_t i=0;i<workset.num_cells;i++) 
        os << wda(workset).cell_local_ids[i] << " ";
     os << "]\n";
     os << "   ir_degrees = [ "; 
     for(std::size_t i=0;i<wda(workset).ir_degrees->size();i++) 
        os << (*wda(workset).ir_degrees)[i] << " ";
     os << "]\n";
     os << "   basis_names = [ "; 
     for(std::size_t i=0;i<wda(workset).basis_names->size();i++) 
        os << (*wda(workset).basis_names)[i] << " ";
     os << "]\n";
     /*
     os << "   int rule = "; wda(workset).int_rules[0]->int_rule->print(os); os << "\n";
     os << "   basis = "; wda(workset).bases[0]->panzer_basis->print(os); os << "\n";

     for(index_t i=0;i<workset.num_cells;i++) {
        os << "   cell " << i << " vertices =\n";
        for(int j=0;j<wda(workset).cell_vertex_coordinates.extent(1);j++) {
           os << "      ";
           for(int k=0;k<wda(workset).cell_vertex_coordinates.extent(2);k++)
              os << wda(workset).cell_vertex_coordinates(i,j,k) << " ";
           os << "\n";
        }
     }

     os << "   integration rule points =\n";
     for(int j=0;j<wda(workset).int_rules[0]->cub_points.extent(0);j++) {
        os << "      ";
        for(int k=0;k<wda(workset).int_rules[0]->cub_points.extent(1);k++)
           os << wda(workset).int_rules[0]->cub_points(j,k) << " ";
        os << "\n";
     }
     os << "   integration weights = [ ";
     for(int j=0;j<wda(workset).int_rules[0]->cub_weights.extent(0);j++) {
        os << wda(workset).int_rules[0]->cub_weights(j) << " ";
     }
     os << "]\n";

     os << "   jac = [ ";
     for(int i=0;i<wda(workset).int_rules[0]->jac.size();i++) {
        os << wda(workset).int_rules[0]->jac[i] << " ";
     }
     os << "]\n";

     os << "   jac_inv = [ ";
     for(int i=0;i<wda(workset).int_rules[0]->jac_inv.size();i++) {
        os << wda(workset).int_rules[0]->jac_inv[i] << " ";
     }
     os << "]\n";

     os << "   jac_det = [ ";
     for(int i=0;i<wda(workset).int_rules[0]->jac_det.size();i++) {
        os << wda(workset).int_rules[0]->jac_det[i] << " ";
     }
     os << "]\n";

     os << "   node_coordinates = [ ";
     for(int i=0;i<wda(workset).int_rules[0]->node_coordinates.size();i++) {
        os << wda(workset).int_rules[0]->node_coordinates[i] << " ";
     }
     os << "]\n";

     os << "   weighted_basis = [ ";
     for(int i=0;i<wda(workset).bases[0]->weighted_basis.size();i++)
        os << wda(workset).bases[0]->weighted_basis[i] << " ";
     os << "]\n";
     
     os << "   weighted_grad_basis = [ ";
     for(int i=0;i<wda(workset).bases[0]->weighted_grad_basis.size();i++)
        os << wda(workset).bases[0]->weighted_grad_basis[i] << " ";
     os << "]\n";

     os << "   basis = [ ";
     for(int i=0;i<wda(workset).bases[0]->basis.size();i++)
        os << wda(workset).bases[0]->basis[i] << " ";
     os << "]\n";
     
     os << "   grad_basis = [ ";
     for(int i=0;i<wda(workset).bases[0]->grad_basis.size();i++)
        os << wda(workset).bases[0]->grad_basis[i] << " ";
     os << "]\n";
     */
  }

  std::vector<std::string>::size_type 
  getPureBasisIndex(std::string basis_name, const panzer::Workset& workset) {
    WorksetDetailsAccessor wda;
    return getPureBasisIndex(basis_name, workset, wda);
  }
  std::vector<std::string>::size_type 
  getBasisIndex(std::string basis_name, const panzer::Workset& workset) {
    WorksetDetailsAccessor wda;
    return getBasisIndex(basis_name, workset, wda);
  }
  std::vector<int>::size_type
  getIntegrationRuleIndex(int ir_degree, const panzer::Workset& workset) {
    WorksetDetailsAccessor wda;
    return getIntegrationRuleIndex(ir_degree, workset, wda);
  }
  void printWorkset(std::ostream& os, const panzer::Workset & workset) {
    WorksetDetailsAccessor wda;
    printWorkset(os, workset, wda);
  }
}

#endif
