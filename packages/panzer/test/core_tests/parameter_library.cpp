#include <Teuchos_ConfigDefs.hpp>
#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_ScalarTraits.hpp>

#include "Panzer_ParameterLibraryEvaluationTraits.hpp"
#include "Panzer_ScalarParameterEntry.hpp"

namespace panzer {

  TEUCHOS_UNIT_TEST(parameter_library, scalar_parameter_entry)
  {
    using Teuchos::RCP;
    using Teuchos::rcp;

    // Residual
    { 
      const RCP<ScalarParameterEntry<panzer::Traits::Residual> > entry = 
	rcp(new ScalarParameterEntry<panzer::Traits::Residual>);
      
      const ScalarParameterEntry<panzer::Traits::Residual>::ScalarT value(5.0);

      entry->setValue(value);
      
      const panzer::Traits::Residual::ScalarT tol = 
	10.0*std::numeric_limits<panzer::Traits::Residual::ScalarT>::epsilon();
      const double double_tol = 10.0*std::numeric_limits<double>::epsilon();
      
      TEST_FLOATING_EQUALITY(entry->getValue(), value, tol);
      TEST_FLOATING_EQUALITY(entry->getRealValue(), 5.0, double_tol);

      const double double_value = 6.0;

      entry->setRealValue(double_value);
      
      TEST_FLOATING_EQUALITY(entry->getRealValue(), double_value, double_tol);
    }

    // Jacobian
    {
      RCP<ScalarParameterEntry<panzer::Traits::Jacobian> > entry = 
	rcp(new ScalarParameterEntry<panzer::Traits::Jacobian>);
      
      typedef Sacado::ScalarType<panzer::Traits::Jacobian::ScalarT>::type ValueType;

      ScalarParameterEntry<panzer::Traits::Jacobian>::ScalarT param(5.0);
      param.resize(1);
      ValueType d_val = ValueType(1.0);
      param.fastAccessDx(0) = d_val;

      entry->setValue(param);
      
      ValueType tol(10.0*std::numeric_limits<ValueType>::epsilon());
      const double double_tol = 10.0*std::numeric_limits<double>::epsilon();
      
      TEST_FLOATING_EQUALITY(entry->getValue().val(), param.val(), tol);
      TEST_FLOATING_EQUALITY(entry->getValue().fastAccessDx(0), d_val, tol);
      TEST_FLOATING_EQUALITY(entry->getRealValue(), 5.0, double_tol);

      const double double_value = 6.0;

      entry->setRealValue(double_value);
      
      TEST_FLOATING_EQUALITY(entry->getRealValue(), double_value, double_tol);
    }

  }

  TEUCHOS_UNIT_TEST(parameter_library, scalar_parameter_library)
  {
    using Teuchos::RCP;
    using Teuchos::rcp;
    using panzer::ParamLib;

    RCP<ParamLib> pl = rcp(new ParamLib);

    {
      RCP<ScalarParameterEntry<panzer::Traits::Residual> > entry = 
	rcp(new ScalarParameterEntry<panzer::Traits::Residual>);
      
      std::string name = "viscosity";
      
      if (!pl->isParameter(name))
	pl->addParameterFamily(name,true,false);
      
      pl->addEntry<panzer::Traits::Residual>("viscosity",entry);
    }

    {
      RCP<ScalarParameterEntry<panzer::Traits::Jacobian> > entry = 
	rcp(new ScalarParameterEntry<panzer::Traits::Jacobian>);

      std::string name = "viscosity";
      
      if (!pl->isParameter(name))
	pl->addParameterFamily(name,true,false);
      
      pl->addEntry<panzer::Traits::Jacobian>("viscosity",entry);
    }

  }

}
