#ifndef PANZER_PARAMETER_LIBRARY_UTILITIES_HPP
#define PANZER_PARAMETER_LIBRARY_UTILITIES_HPP

#include "Panzer_ParameterLibrary.hpp"
#include "Panzer_ScalarParameterEntry.hpp"
#include "Teuchos_RCP.hpp"

namespace panzer {

  /** \brief Allocates a parameter entry and registers with parameter library
      
  \relates ParameterLibraryAcceptor
  */
  template<typename EvaluationType>
  Teuchos::RCP<panzer::ScalarParameterEntry<EvaluationType> >
  createAndRegisterScalarParameter(const std::string name,
				   panzer::ParamLib& pl);

}

#include "Panzer_ParameterLibraryUtilitiesT.hpp"

#endif
