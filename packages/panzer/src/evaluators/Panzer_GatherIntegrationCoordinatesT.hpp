#include "Teuchos_Assert.hpp"
#include "Phalanx_DataLayout.hpp"

#include "Panzer_IntegrationRule.hpp"
#include "Panzer_Workset_Utilities.hpp"

#include "Teuchos_FancyOStream.hpp"

template<typename EvalT,typename Traits>
std::string 
panzer::GatherIntegrationCoordinates<EvalT, Traits>::
fieldName(int degree)
{
   std::stringstream ss; 
   ss << "IR_" << degree << " IntegrationCoordinates";
   return ss.str();
}

template<typename EvalT,typename Traits>
panzer::GatherIntegrationCoordinates<EvalT, Traits>::
GatherIntegrationCoordinates(const panzer::IntegrationRule & quad)
{ 
  quadDegree_ = quad.cubature_degree;

  quadCoordinates_ = PHX::MDField<ScalarT,Cell,Point,Dim>(fieldName(quadDegree_),quad.dl_vector);

  this->addEvaluatedField(quadCoordinates_);

  this->setName("Gather "+fieldName(quadDegree_));
}

// **********************************************************************
template<typename EvalT,typename Traits>
void panzer::GatherIntegrationCoordinates<EvalT, Traits>::
postRegistrationSetup(typename Traits::SetupData sd, 
		      PHX::FieldManager<Traits>& fm)
{
  this->utils.setFieldData(quadCoordinates_,fm);

  quadIndex_ = panzer::getIntegrationRuleIndex(quadDegree_, (*sd.worksets_)[0]);
}

// **********************************************************************
template<typename EvalT,typename Traits> 
void panzer::GatherIntegrationCoordinates<EvalT, Traits>::
evaluateFields(typename Traits::EvalData workset)
{ 
  const Intrepid::FieldContainer<double> & quadCoords = workset.int_rules[quadIndex_]->ip_coordinates;  

  // just copy the array
  for(int i=0;i<quadCoords.size();i++)
     quadCoordinates_[i] = quadCoords[i];
}
