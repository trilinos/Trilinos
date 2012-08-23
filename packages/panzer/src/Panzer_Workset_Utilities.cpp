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
  getBasisIndex(std::string basis_name, panzer::Workset& workset)
  {
    std::vector<std::string>::iterator basis;

    basis = std::find(workset.basis_names->begin(),
		      workset.basis_names->end(),
		      basis_name);

    TEUCHOS_TEST_FOR_EXCEPTION(basis == workset.basis_names->end(),
			       std::logic_error,
			       "Could not find the basis named \"" 
                               << basis_name << "\" in the workset!");

    return std::distance(workset.basis_names->begin(), basis);
  }

  std::vector<std::string>::size_type
  getIntegrationRuleIndex(int ir_degree, panzer::Workset& workset)
  {
    std::vector<int>::iterator ir;

    ir = std::find(workset.ir_degrees->begin(),
		   workset.ir_degrees->end(),
		   ir_degree);
    
    TEUCHOS_TEST_FOR_EXCEPTION(ir == workset.ir_degrees->end(),
			       std::logic_error,
			       "Could not find the integration rule degree \"" 
                               << ir_degree << "\" in the workset!");

    return std::distance(workset.ir_degrees->begin(), ir);
  }

  void printWorkset(std::ostream& os, const panzer::Workset & workset)
  {
     os << "WORKSET"
        << " block_id = \"" << workset.block_id << "\""
        << " num_cells = \"" << workset.num_cells << "\"\n";
     os << "   cell_local_ids = [ "; 
     for(std::size_t i=0;i<workset.num_cells;i++) 
        os << workset.cell_local_ids[i] << " ";
     os << "]\n";
     os << "   ir_degrees = [ "; 
     for(std::size_t i=0;i<workset.ir_degrees->size();i++) 
        os << (*workset.ir_degrees)[i] << " ";
     os << "]\n";
     os << "   basis_names = [ "; 
     for(std::size_t i=0;i<workset.basis_names->size();i++) 
        os << (*workset.basis_names)[i] << " ";
     os << "]\n";
     /*
     os << "   int rule = "; workset.int_rules[0]->int_rule->print(os); os << "\n";
     os << "   basis = "; workset.bases[0]->panzer_basis->print(os); os << "\n";

     for(std::size_t i=0;i<workset.num_cells;i++) {
        os << "   cell " << i << " vertices =\n";
        for(int j=0;j<workset.cell_vertex_coordinates.dimension(1);j++) {
           os << "      ";
           for(int k=0;k<workset.cell_vertex_coordinates.dimension(2);k++)
              os << workset.cell_vertex_coordinates(i,j,k) << " ";
           os << "\n";
        }
     }

     os << "   integration rule points =\n";
     for(int j=0;j<workset.int_rules[0]->cub_points.dimension(0);j++) {
        os << "      ";
        for(int k=0;k<workset.int_rules[0]->cub_points.dimension(1);k++)
           os << workset.int_rules[0]->cub_points(j,k) << " ";
        os << "\n";
     }
     os << "   integration weights = [ ";
     for(int j=0;j<workset.int_rules[0]->cub_weights.dimension(0);j++) {
        os << workset.int_rules[0]->cub_weights(j) << " ";
     }
     os << "]\n";

     os << "   jac = [ ";
     for(int i=0;i<workset.int_rules[0]->jac.size();i++) {
        os << workset.int_rules[0]->jac[i] << " ";
     }
     os << "]\n";

     os << "   jac_inv = [ ";
     for(int i=0;i<workset.int_rules[0]->jac_inv.size();i++) {
        os << workset.int_rules[0]->jac_inv[i] << " ";
     }
     os << "]\n";

     os << "   jac_det = [ ";
     for(int i=0;i<workset.int_rules[0]->jac_det.size();i++) {
        os << workset.int_rules[0]->jac_det[i] << " ";
     }
     os << "]\n";

     os << "   node_coordinates = [ ";
     for(int i=0;i<workset.int_rules[0]->node_coordinates.size();i++) {
        os << workset.int_rules[0]->node_coordinates[i] << " ";
     }
     os << "]\n";

     os << "   weighted_basis = [ ";
     for(int i=0;i<workset.bases[0]->weighted_basis.size();i++)
        os << workset.bases[0]->weighted_basis[i] << " ";
     os << "]\n";
     
     os << "   weighted_grad_basis = [ ";
     for(int i=0;i<workset.bases[0]->weighted_grad_basis.size();i++)
        os << workset.bases[0]->weighted_grad_basis[i] << " ";
     os << "]\n";

     os << "   basis = [ ";
     for(int i=0;i<workset.bases[0]->basis.size();i++)
        os << workset.bases[0]->basis[i] << " ";
     os << "]\n";
     
     os << "   grad_basis = [ ";
     for(int i=0;i<workset.bases[0]->grad_basis.size();i++)
        os << workset.bases[0]->grad_basis[i] << " ";
     os << "]\n";
     */
  }

}

#endif
