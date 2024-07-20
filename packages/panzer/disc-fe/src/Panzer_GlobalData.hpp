// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PANZER_GLOBAL_DATA_HPP
#define PANZER_GLOBAL_DATA_HPP

#include "Teuchos_RCP.hpp"
#include "Teuchos_FancyOStream.hpp"
#include "Panzer_ParameterLibrary.hpp"

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

    /** \brief Get the output stream */
    std::ostream & out() 
    { return *os; }
    
  };
  
  /** \brief Nonmember constructor
      
    Allocates a new global data object.  Automatically allocates the
    sacado parameter libraries.

    \param [in] build_default_os If set to true, the os object will be
    allocated with a pointer to cout.

    \param [in] print_process Sets the print process if the os object
    is built by this method.

    \relates GlobalData
  */
  Teuchos::RCP<panzer::GlobalData> createGlobalData(bool build_default_os = true,
						    int print_process = 0);

}

#endif
