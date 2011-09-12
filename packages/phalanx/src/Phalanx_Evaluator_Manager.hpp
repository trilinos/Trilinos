// @HEADER
// ************************************************************************
//
//        Phalanx: A Partial Differential Equation Field Evaluation 
//       Kernel for Flexible Management of Complex Dependency Chains
//                    Copyright 2008 Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
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
// Questions? Contact Roger Pawlowski (rppawlo@sandia.gov), Sandia
// National Laboratories.
//
// ************************************************************************
// @HEADER


#ifndef PHX_FIELD_EVALUATOR_MANAGER_HPP
#define PHX_FIELD_EVALUATOR_MANAGER_HPP

#include <string>
#include <vector>
#include <iostream>
#include "Teuchos_RCP.hpp"
#include "Phalanx_ConfigDefs.hpp"
#include "Phalanx_FieldTag.hpp"
#include "Phalanx_FieldTag_STL_Functors.hpp"
#include "Phalanx_Evaluator.hpp"
#include "Phalanx_TypeStrings.hpp"
#include "Teuchos_TimeMonitor.hpp"

namespace PHX {
  
  template<typename Traits> class FieldManager;

  /*! @brief Class to sort which Evaluators should be called and the order in which to call them such that all dependencies are met.
   */
  template<typename Traits>
  class EvaluatorManager {

  public:

    EvaluatorManager(const std::string& evaluator_type_name = "???");
    
    ~EvaluatorManager();
    
    //! Require a variable to be evaluated.
    void requireField(const PHX::FieldTag& v);
    
    //! Registers a variable provider with the manager.
    void 
    registerEvaluator(const Teuchos::RCP<PHX::Evaluator<Traits> >& p);
    
    /*! Sets up all field dependencies.  This should only be called
      once all variables and DOFs have been added and all providers
      have been registered.  Sorts variable and creates dependency
      lists and evaluation order
    */
    void sortAndOrderEvaluators();
    
    /*! Calls post registration setup on all variable providers.
    */
    void postRegistrationSetup(typename Traits::SetupData d,
			       PHX::FieldManager<Traits>& vm);
    
    //! Compute the required variables for the fill on the specific element.
    void evaluateFields(typename Traits::EvalData d);
    
    /*! \brief This routine is called before each residual/Jacobian fill.
      
        This routine is called ONCE on the provider before the fill
        loop over elements is started.  This allows us to reset global
        objects between each fill.  An example is to reset a provider
        that monitors the maximum grid peclet number in a cell.  This
        call would zero out the maximum for a new fill.
    */
    void preEvaluate(typename Traits::PreEvalData d);
    
    /*! \brief This routine is called after each residual/Jacobian fill.
      
        This routine is called ONCE on the provider after the fill
        loop over elements is completed.  This allows us to evaluate
        any post fill data.  An example is to print out some
        statistics such as the maximum grid peclet number in a cell.
    */
    void postEvaluate(typename Traits::PostEvalData d);
    
    void setEvaluationTypeName(const std::string& evaluation_type_name);
    
    const std::vector< Teuchos::RCP<PHX::FieldTag> >& getFieldTags();

    bool sortingCalled() const;

    void writeGraphvizFile(const std::string filename,
			   bool writeEvaluatedFields,
			   bool writeDependentFields,
			   bool debugRegisteredEvaluators) const;

    //! Printing
    void print(std::ostream& os) const;
    
  protected:
    
    /*! @brief Create and arrange the dependency list in the correct
     *  order it should be evaluated.
     */
    void createProviderEvaluationOrder();
    
  protected:
    
    //! Fields required by the user.
    std::vector< Teuchos::RCP<PHX::FieldTag> > fields_;
    
    //@{
    /*!
      @name Evaluator Objects
      @brief Stores information about variable provider objects.
    */    
    std::vector< Teuchos::RCP<PHX::Evaluator<Traits> > > 
    varProviders;
    
    std::vector< std::vector< Teuchos::RCP<PHX::FieldTag> > > 
    providerVariables;

    std::vector< std::vector< Teuchos::RCP<PHX::FieldTag> > > 
    providerRequirements;

    std::vector<std::string> providerNames;

    // Timers used when configured with Phalanx_ENABLE_TEUCHOS_TIME_MONITOR.
    std::vector<Teuchos::RCP<Teuchos::Time> > evalTimers;

    //@}

    
    //@{
    /*! @name Evaluation Order Objects
      
        Stores information about the order that providers need to be
        called to evaluate fields correctly.
    */
    std::vector<int> providerEvalOrderIndex;
    //@}
    
    std::string evaluation_type_name_;

    //! Flag to tell the setup has been called.
    bool sorting_called_;
    
  };
  
  template<typename Traits>
  std::ostream& operator<<(std::ostream& os, 
			   const PHX::EvaluatorManager<Traits>& m);

}

#include "Phalanx_Evaluator_Manager_Def.hpp"

#endif
