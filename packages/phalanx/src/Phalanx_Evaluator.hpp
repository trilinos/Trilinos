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


#ifndef PHX_FIELDEVALUATOR_H
#define PHX_FIELDEVALUATOR_H

#include <vector>

#include "Teuchos_RCP.hpp"
#include "Phalanx_FieldTag.hpp"

namespace PHX {

  template<typename Traits> class FieldManager;

  /*! Pure virtual base class that provides field evaluation
      routines to the FieldManager.
  */
  template <typename Traits>
  class Evaluator {
    
  public:
    
    //! Ctor
    Evaluator() {};
    
    //! Dtor
    virtual ~Evaluator() {};

    /*! \brief Allows providers to grab pointers to data arrays.
      
        Called once all providers are registered with the manager.
	
	Once the field manager has allocated all data arrays, this
	method passes the field manager to the providers to allow each
	provider to grab and store pointers to the field data arrays.
	Grabbing the data arrays from the varible manager during an
	actual call to evaluateFields call is too slow due to the map
	lookup and FieldTag comparison (which uses a string compare).
	So lookups on field data are only allowed during this setup
	phase.
	
    */
    virtual void postRegistrationSetup(typename Traits::SetupData d,
				       PHX::FieldManager<Traits>& vm) = 0;

    //! Returns vector of fields that this object evaluates.
    virtual const std::vector< Teuchos::RCP<FieldTag> >& 
    evaluatedFields() const = 0;

    //! Returns vector of fields needed to compute the evaluated fields.
    virtual const std::vector< Teuchos::RCP<FieldTag> >& 
    dependentFields() const = 0;

    //! Evaluate all fields that the provider supplies.
    /*!
        Input:
	@param d - user defined data object defined by the EvalData typedef in the traits class.
    */ 
    virtual void evaluateFields(typename Traits::EvalData d) = 0;
    
    /*! \brief This routine is called before each residual/Jacobian fill.

        This routine is called ONCE on the provider before the fill
        loop over cells is started.  This allows us to reset global
        objects between each fill.  An example is to reset a provider
        that monitors the maximum grid peclet number in a cell.  This
        call would zero out the maximum for a new fill.
    */
    virtual void preEvaluate(typename Traits::PreEvalData d) = 0;

    /*! \brief This routine is called after each residual/Jacobian fill.

        This routine is called ONCE on the provider after the fill
        loop over cells is completed.  This allows us to evaluate any
        post fill data.  An example is to print out some statistics
        such as the maximum grid peclet number in a cell.
    */
    virtual void postEvaluate(typename Traits::PostEvalData d) = 0;

    //! Returns the name/identifier of this provider.
    virtual const std::string& getName() const = 0;

  };

} 

#endif 
