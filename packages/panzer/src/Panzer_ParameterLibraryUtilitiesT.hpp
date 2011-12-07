#ifndef PANZER_PARAMETER_LIBRARY_UTILITIES_T_HPP
#define PANZER_PARAMETER_LIBRARY_UTILITIES_T_HPP

namespace panzer {

  template<typename EvaluationType>
  Teuchos::RCP<panzer::ScalarParameterEntry<EvaluationType> >
  createAndRegisterScalarParameter(const std::string name, 
				   panzer::ParamLib& pl)
  {
    Teuchos::RCP<panzer::ScalarParameterEntry<EvaluationType> > entry = 
      Teuchos::rcp(new panzer::ScalarParameterEntry<EvaluationType>);
    
    if (!pl.isParameter(name))
      pl.addParameterFamily(name,true,false);
    
    pl.addEntry<EvaluationType>("viscosity",entry);

    return entry;
  }

}

#endif
