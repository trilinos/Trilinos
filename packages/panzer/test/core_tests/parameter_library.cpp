#include <Teuchos_ConfigDefs.hpp>
#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_RCP.hpp>

#include "Panzer_ParameterLibraryEvaluationTraits.hpp"

namespace panzer {

  template <typename EvalType>
  class SimpleEntry : public  Sacado::ScalarParameterEntry<EvalType,panzer::EvaluationTraits> {

    typedef typename Sacado::ScalarParameterEntry<EvalType,panzer::EvaluationTraits>::ScalarT ScalarT;

    ScalarT m_parameter;

  public:

    void setRealValue(double value)
    {
      m_parameter = ScalarT(value);
    }

    void setValue(const ScalarT& value)
    {
      m_parameter = value;
    }

    const ScalarT& getValue() const
    {
      return m_parameter;
    }
    

  };

  TEUCHOS_UNIT_TEST(global_data_accessor, default_impl)
  {
    using Teuchos::RCP;
    using Teuchos::rcp;
    using panzer::ParamLib;

    RCP<ParamLib> pl = rcp(new ParamLib);

    RCP<SimpleEntry<panzer::Traits::Residual> > entry = 
      rcp(new SimpleEntry<panzer::Traits::Residual>);

    std::string name = "viscosity";

    if (!pl->isParameter(name))
      pl->addParameterFamily(name,true,false);
    
    pl->addEntry<panzer::Traits::Residual>("viscosity",entry);

  }

}
