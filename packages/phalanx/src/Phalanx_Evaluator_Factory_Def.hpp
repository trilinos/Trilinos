#ifndef PHX_FIELD_EVALUATOR_FACTORY_DEF_HPP
#define PHX_FIELD_EVALUATOR_FACTORY_DEF_HPP

#include <sstream>
#include "Sacado_mpl_size.hpp"
#include "boost/mpl/for_each.hpp"
#include "boost/mpl/range_c.hpp"
#include "Teuchos_TestForException.hpp"
#include "Phalanx_Evaluator_Factory_UFO.hpp"

//**********************************************************************
template<typename Traits, typename FactoryTraits>
Teuchos::RCP< std::vector< Teuchos::RCP<PHX::Evaluator_TemplateManager<Traits> > > > 
PHX::EvaluatorFactory<Traits, FactoryTraits>::
buildEvaluators(const std::map<std::string, Teuchos::RCP<Teuchos::ParameterList> >& data)
{
  using namespace std;
  using namespace Teuchos;
  using namespace PHX;

  RCP< std::vector< RCP<Evaluator_TemplateManager<Traits> > > > vector_tm = 
    rcp(new std::vector< RCP<Evaluator_TemplateManager<Traits> > >);

  std::map<std::string, RCP<ParameterList> >::const_iterator it = data.begin();
  
  for (; it != data.end(); ++it) {
    
    RCP< Evaluator_TemplateManager<Traits> > tm = 
      rcp(new Evaluator_TemplateManager<Traits>);
    
    RCP<ParameterList> p = it->second;

    bool found_type = false;
    int object_type = p->get<int>("Type");
    static const int size = Sacado::mpl::size<typename FactoryTraits::EvaluatorTypes>::value;
    boost::mpl::for_each< boost::mpl::range_c<int,0,size> >( UFO<Traits,FactoryTraits>(object_type, p, tm, found_type) );

    if (!found_type) {
      ostringstream msg;
      msg << "Unable to find model in EvaluatorFactory for "
	  << it->first
	  << ".  Please make sure you have set a valid integer "
	  << "for \"Type\" in the parameter list!";
      TEST_FOR_EXCEPTION(!found_type, std::logic_error, msg.str());
    }

    vector_tm->push_back(tm);

  }
  
  return vector_tm;

}
//**********************************************************************
template<typename Traits>
void PHX::registerEvaluators(const Teuchos::RCP< std::vector< Teuchos::RCP<PHX::Evaluator_TemplateManager<Traits> > > >& providers, PHX::FieldManager<Traits>& fm)
{
  using namespace std;
  using namespace Teuchos;
  // Loop over each provider template manager
  typename vector< RCP<Evaluator_TemplateManager<Traits> > >::iterator 
    tm = providers->begin();
  for (; tm != providers->end(); ++tm) {
    
    // Loop over Evaluation Types
    typename PHX::FieldManager<Traits>::iterator vmit = fm.begin();
    typename Evaluator_TemplateManager<Traits>::iterator vpit = 
      (*tm)->begin();
    for (; vpit != (*tm)->end(); ++vpit) {
      RCP<PHX::Evaluator<Traits> > vp =
	rcp_dynamic_cast<PHX::Evaluator<Traits> >(vpit.rcp());
      fm.registerEvaluator(vmit, vp);
      ++vmit;
    } 
  }

}
//**********************************************************************

#endif
