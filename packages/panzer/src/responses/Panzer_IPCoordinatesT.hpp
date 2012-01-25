#ifndef PANZER_EVALUATOR_IP_COORDINATES_T_HPP
#define PANZER_EVALUATOR_IP_COORDINATES_T_HPP

#include "Panzer_Workset_Utilities.hpp"
#include "Phalanx_DataLayout_MDALayout.hpp"
#include "Teuchos_Assert.hpp"
#include <cstddef>
#include <string>
#include <vector>

namespace panzer {

//**********************************************************************
template<typename EvalT, typename Traits>
IPCoordinates<EvalT, Traits>::
IPCoordinates(int in_ir_order,
	      const Teuchos::RCP<std::vector<Coordinate> >& in_coords) :
  first_evaluation(true),
  ir_order(in_ir_order),
  coords(in_coords)
{ 
  using Teuchos::RCP;
  using Teuchos::rcp;
  using PHX::MDALayout;
  using PHX::MDField;

  TEUCHOS_ASSERT(nonnull(coords));

  RCP<MDALayout<Cell,IP,Dim> > dl = rcp(new MDALayout<Cell,IP,Dim>(0,0,0));
  dummy_field = PHX::MDField<ScalarT, Cell,IP,Dim>("IP Coordinates", dl);
  
  this->addEvaluatedField(dummy_field);
 
  this->setName("IPCoordinates Evaluator");
}

//**********************************************************************
template<typename EvalT, typename Traits>
void IPCoordinates<EvalT, Traits>::
postRegistrationSetup(typename Traits::SetupData sd,
		      PHX::FieldManager<Traits>& fm)
{
  this->utils.setFieldData(dummy_field,fm);

  ir_index = panzer::getIntegrationRuleIndex(ir_order,(*sd.worksets_)[0]);
}

//**********************************************************************
template<typename EvalT, typename Traits>
void IPCoordinates<EvalT, Traits>::
evaluateFields(typename Traits::EvalData workset)
{ 
  
  if (first_evaluation) {

    Intrepid::FieldContainer<double>& workset_coords = (workset.int_rules[ir_index])->ip_coordinates;

    for (std::size_t cell = 0; cell < workset.num_cells; ++cell) {

      for (int ip = 0; ip < workset_coords.dimension(1); ++ip) {

	Coordinate pt;
	pt.lid = workset.cell_local_ids[cell];
	int num_dims = workset_coords.dimension(2);
	pt.val.resize(num_dims);

	for (int dim = 0; dim < num_dims; ++dim)
	  pt.val[dim] = workset_coords(static_cast<int>(cell),ip,dim); 

	coords->push_back(pt);
      }

    }

  }

}

//**********************************************************************
template<typename EvalT, typename Traits>
void IPCoordinates<EvalT, Traits>::preEvaluate(typename Traits::PreEvalData data)
{
  
}

//**********************************************************************
template<typename EvalT, typename Traits>
void IPCoordinates<EvalT, Traits>::postEvaluate(typename Traits::PostEvalData data)
{
  first_evaluation = false;
}

//**********************************************************************
template<typename EvalT, typename Traits>
const PHX::MDField<typename IPCoordinates<EvalT, Traits>::ScalarT,Cell,IP,Dim> 
IPCoordinates<EvalT, Traits>::getEvaluatedField() const
{
  return dummy_field;
}

//**********************************************************************

}

#endif
