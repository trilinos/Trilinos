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

#ifndef PANZER_NORMALS_IMPL_HPP
#define PANZER_NORMALS_IMPL_HPP

#include <algorithm>
#include "Panzer_IntegrationRule.hpp"
#include "Panzer_Workset_Utilities.hpp"
#include "Intrepid_FunctionSpaceTools.hpp"
#include "Intrepid_CellTools.hpp"

namespace panzer {

//**********************************************************************
PHX_EVALUATOR_CTOR(Normals,p)
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
PHX_POST_REGISTRATION_SETUP(Normals,sd,fm)
{
  this->utils.setFieldData(normals,fm);

  num_qp  = normals.dimension(1);
  num_dim = normals.dimension(2);
  
  quad_index =  panzer::getIntegrationRuleIndex(quad_order,(*sd.worksets_)[0]);
}

//**********************************************************************
PHX_EVALUATE_FIELDS(Normals,workset)
{ 
   if(workset.num_cells>0) {
      Intrepid::CellTools<ScalarT>::getPhysicalSideNormals(normals,
                                        workset.int_rules[quad_index]->jac,
                                        side_id, *workset.int_rules[quad_index]->int_rule->topology);
      
      if(normalize) {
         // normalize vector: getPhysicalSideNormals does not 
         // return normalized vectors
         for(std::size_t c=0;c<workset.num_cells;c++) {
            for(std::size_t q=0;q<num_qp;q++) {
               ScalarT norm = 0.0;
   
               // compute squared norm
               for(std::size_t d=0;d<num_dim;d++)
                  norm += normals(c,q,d)*normals(c,q,d);
    
               // adjust for length of vector, now unit vectors
               norm = sqrt(norm);
               for(std::size_t d=0;d<num_dim;d++)
                  normals(c,q,d) /= norm;
            }
         }
      }
      // else normals correspond to differential
   }
 
}

//**********************************************************************

}

#endif
