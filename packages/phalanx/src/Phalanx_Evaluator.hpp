// @HEADER
// ************************************************************************
// 
//        Phalanx: A Partial Differential Equation Field Evaluation 
//       Kernel for Flexible Management of Complex Dependency Chains
//                  Copyright (2008) Sandia Corporation
// 
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
// 
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
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
