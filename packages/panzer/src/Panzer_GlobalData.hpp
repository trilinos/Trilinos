#ifndef PANZER_GLOBAL_DATA_HPP
#define PANZER_GLOBAL_DATA_HPP

#include "Teuchos_RCP.hpp"
#include "Teuchos_FancyOStream.hpp"
#include "Panzer_ParameterLibraryEvaluationTraits.hpp"

namespace panzer {
  
  /** \brief Struct for global data to be stored.

     This object is unique for each instantiation of a panzer model
     evaluator.  It is intended to store data that is usually global
     to the code, but must be protected to allow for multiple
     instantiations of the application (say for multiphysics coupling
     of two panzer applications).  It stores things like the parameter
     library and a default ostream to redirect all output to.
  */
  struct GlobalData {
    
    /** \brief ostream for redirecting all panzer output for a particular instantiation.  */
    Teuchos::RCP<Teuchos::FancyOStream> os;
    
    /** \brief Sacado scalar parameter library */
    Teuchos::RCP<panzer::ParamLib> pl;
    
    /** \brief Sacado scalar parameter vector library */
    Teuchos::RCP<panzer::ParamVec> pv;
    
  };
  
  /** \brief Nonmember constructor
      
    Allocates a new global data object.  Automatically allocates the
    sacado scalar and vector parameter libraries.

    \relates GlobalData
  */
  Teuchos::RCP<panzer::GlobalData> createGlobalData();

}

#endif
