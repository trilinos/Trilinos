
#include <algorithm>
#include "Panzer_IntegrationRule.hpp"
#include "Panzer_Basis.hpp"
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
  
  quad_index =  
    std::distance((*sd.worksets_)[0].ir_degrees->begin(),
		  std::find((*sd.worksets_)[0].ir_degrees->begin(),
			    (*sd.worksets_)[0].ir_degrees->end(),
			    quad_order));
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
