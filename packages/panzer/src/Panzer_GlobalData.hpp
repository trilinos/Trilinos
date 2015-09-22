// @HEADER
// ***********************************************************************
//
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//                 Copyright (2011) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Roger P. Pawlowski (rppawlo@sandia.gov) and
// Eric C. Cyr (eccyr@sandia.gov)
// ***********************************************************************
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
