#include "Panzer_IntegrationRule.hpp"
#include "Panzer_Workset_Utilities.hpp"

#include "Intrepid_FunctionSpaceTools.hpp"

#include "Phalanx_DataLayout_MDALayout.hpp"

namespace panzer {

//**********************************************************************
template <typename EvalT, typename TraitsT>
void DOF_PointField<EvalT,TraitsT>::initialize(const std::string & fieldName,
                                                const PureBasis & fieldBasis,
                                                const std::string & coordinateName,
                                                const Teuchos::RCP<PHX::DataLayout> & coordLayout,
                                                const std::string & postfixFieldName)
{
  intrepidBasis = fieldBasis.getIntrepidBasis();

  int cellCount = fieldBasis.functional->dimension(0);
  int coeffCount = fieldBasis.functional->dimension(1);
  int pointCount = coordLayout->dimension(0);
  int dimCount = coordLayout->dimension(1);

  Teuchos::RCP<PHX::DataLayout> basisLayout = fieldBasis.functional;
  Teuchos::RCP<PHX::DataLayout> quadLayout = Teuchos::rcp(new PHX::MDALayout<Cell,Point>(cellCount,pointCount));

  coordinates = PHX::MDField<ScalarT,Point,Dim>(coordinateName,coordLayout);
  dof_coeff = PHX::MDField<ScalarT,Cell,Point>(fieldName,basisLayout);
  dof_field = PHX::MDField<ScalarT,Cell,Point>(fieldName+postfixFieldName,quadLayout);

  this->addDependentField(coordinates);
  this->addDependentField(dof_coeff);
  this->addEvaluatedField(dof_field);

  // build data storage for temporary conversion
  basisRef    = Intrepid::FieldContainer<double>(coeffCount,pointCount);
  basis       = Intrepid::FieldContainer<double>(cellCount,coeffCount,pointCount);
  intrpCoords = Intrepid::FieldContainer<double>(pointCount,dimCount);
  
  std::string n = "DOF_PointField: " + dof_field.fieldTag().name();
  this->setName(n);
}

//**********************************************************************
template <typename EvalT, typename TraitsT>
void DOF_PointField<EvalT,TraitsT>::postRegistrationSetup(typename TraitsT::SetupData d,
			                                  PHX::FieldManager<TraitsT>& fm)
{
  this->utils.setFieldData(coordinates,fm);
  this->utils.setFieldData(dof_coeff,fm);
  this->utils.setFieldData(dof_field,fm);
}

//**********************************************************************
template <typename EvalT, typename TraitsT>
void DOF_PointField<EvalT,TraitsT>::evaluateFields(typename TraitsT::EvalData workset)
{ 
  // Zero out arrays (intrepid does a sum! 1/17/2012)
  for (int i = 0; i < dof_field.size(); ++i)
    dof_field[i] = 0.0;

  // copy coordinates
  for (int i = 0; i < coordinates.size(); ++i)
    intrpCoords[i] = Sacado::ScalarValue<ScalarT>::eval(coordinates[i]);

  if(workset.num_cells>0) {
    // evaluate at reference points
    intrepidBasis->getValues(basisRef, intrpCoords, Intrepid::OPERATOR_VALUE);

    // transfer reference basis values to physical frame values
    Intrepid::FunctionSpaceTools::
      HGRADtransformVALUE<ScalarT>(basis,
				  basisRef);

    // evaluate function at specified points
    Intrepid::FunctionSpaceTools::
      evaluate<ScalarT>(dof_field,dof_coeff,basis);
  }
}

//**********************************************************************

}
