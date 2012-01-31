#ifndef PANZER_STK_SCATTER_CELL_AVG_QUANTITY_IMPL_HPP
#define PANZER_STK_SCATTER_CELL_AVG_QUANTITY_IMPL_HPP

#include "Teuchos_Assert.hpp"

#include "Phalanx_ConfigDefs.hpp"
#include "Phalanx_Evaluator_Macros.hpp"
#include "Phalanx_MDField.hpp"
#include "Phalanx_DataLayout.hpp"
#include "Phalanx_DataLayout_MDALayout.hpp"

#include "Panzer_IntegrationRule.hpp"

#include "Teuchos_FancyOStream.hpp"
#include "Teuchos_ArrayRCP.hpp"

namespace panzer_stk {

PHX_EVALUATOR_CTOR(ScatterCellAvgQuantity,p) :
   mesh_(p.get<Teuchos::RCP<STK_Interface> >("Mesh"))
{
  using panzer::Cell;
  using panzer::Point;

  std::string scatterName = p.get<std::string>("Scatter Name");
 
  const std::vector<std::string> & names = 
    *(p.get< Teuchos::RCP< std::vector<std::string> > >("Field Names"));

  Teuchos::RCP<panzer::IntegrationRule> intRule = 
    p.get< Teuchos::RCP<panzer::IntegrationRule> >("IR");

  // build dependent fields
  scatterFields_.resize(names.size());
  scatterFields_dbl_.resize(names.size());
  stkFields_.resize(names.size());
  for (std::size_t fd = 0; fd < names.size(); ++fd) {
    scatterFields_[fd] = 
      PHX::MDField<ScalarT,Cell,Point>(names[fd],intRule->dl_scalar);
    this->addDependentField(scatterFields_[fd]);

    scatterFields_dbl_[fd] =
      PHX::MDField<double,Cell,Point>(names[fd],intRule->dl_scalar);
    scatterFields_dbl_[fd].setFieldData(Teuchos::arcp<double>(intRule->dl_scalar->size()));
  }

  // setup a dummy field to evaluate
  PHX::Tag<ScalarT> scatterHolder(scatterName,Teuchos::rcp(new PHX::MDALayout<panzer::Dummy>(0)));
  this->addEvaluatedField(scatterHolder);

  this->setName(scatterName+": STK-Scatter Cell Fields");
}

PHX_POST_REGISTRATION_SETUP(ScatterCellAvgQuantity,d,fm)
{
  for (std::size_t fd = 0; fd < scatterFields_.size(); ++fd) {
    std::string fieldName = scatterFields_[fd].fieldTag().name();

    stkFields_[fd] = mesh_->getMetaData()->get_field<VariableField>(fieldName);

    // setup the field data object
    this->utils.setFieldData(scatterFields_[fd],fm);
  }
}

PHX_EVALUATE_FIELDS(ScatterCellAvgQuantity,workset)
{
   // for convenience pull out some objects from workset
   const std::vector<std::size_t> & localCellIds = workset.cell_local_ids;
   std::string blockId = workset.block_id;

   for(std::size_t fieldIndex=0; fieldIndex<scatterFields_.size();fieldIndex++) {
      PHX::MDField<ScalarT,panzer::Cell,panzer::Point> & field = scatterFields_[fieldIndex];
      std::vector<double> average(field.dimension(0),0.0);
      // write to double field
      for(int i=0; i<field.dimension(0);i++) {
         for(int j=0; j<field.dimension(1);j++) 
            average[i] += Sacado::ScalarValue<ScalarT>::eval(field(i,j));
         average[i] /= field.dimension(0);
      }

      mesh_->setCellFieldData(field.fieldTag().name(),blockId,localCellIds,average);
   }
}

} // end panzer_stk

#endif
