// $Id$ 
// $Source$ 

//@HEADER
// ************************************************************************
// 
//            NOX: An Object-Oriented Nonlinear Solver Package
//                 Copyright (2002) Sandia Corporation
// 
//            LOCA: Library of Continuation Algorithms Package
//                 Copyright (2005) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
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
// Questions? Contact Roger Pawlowski (rppawlo@sandia.gov) or 
// Eric Phipps (etphipp@sandia.gov), Sandia National Laboratories.
// ************************************************************************
//  CVS Information
//  $Source$
//  $Author$
//  $Date$
//  $Revision$
// ************************************************************************
//@HEADER

#ifndef NOX_ROWSUMSCALING_UTILITIES_H
#define NOX_ROWSUMSCALING_UTILITIES_H

#include "NOX_Common.H"  // for NOX_Config.h
#include "Thyra_ModelEvaluator.hpp"
#include "Thyra_ScaledModelEvaluator.hpp"
#include "NOX_Abstract_PrePostOperator.H"
#include "NOX_PrePostOperator_Vector.H"
#include "NOX_PrePostOperator_RowSumScaling.H"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_TestForException.hpp"

namespace NOX {

  namespace Thyra {

  /** Nonmember function that creates objects required for row sum scaling.

    To support row sum scaling, two objects must be constructed.  The
    first object is a Thyra::ScaledModelEvaluator that decorates the
    application model evaluator and explicity applies the scaling to
    the underlying objects.  The second object is a NOX::RowSumScaling
    is a NOX::PrePostOperator that controls when the scaling values
    should be updated with row sum values.  These objects share a
    pointer to the scaling vector that must be registered with both
    objects.  This class creates all objects and takes care of the
    registrations.

    \param[in] application_model User supplied model evaluator 
    \param[out] scaled_model Newly allocated scaled model evaluator
    \param[out] scaled_ppo Newly created PrePostOperator that *must* be registered with the solver

  */
    template<typename Scalar>
    void buildRowSumScaledObjects(const Teuchos::RCP<const ::Thyra::ModelEvaluator<Scalar> >& application_model,
				  Teuchos::RCP< ::Thyra::ScaledModelEvaluator<Scalar> >& scaled_model,
				  Teuchos::RCP<NOX::Abstract::PrePostOperator>& scaled_ppo,
				  Teuchos::ParameterList& nox_parameters)
    {
      using Teuchos::RCP;
      using Teuchos::rcp;
      
      // Create f_scaling vector
      RCP< ::Thyra::VectorBase<double> > scale_vec = 
	::Thyra::createMember(application_model->get_f_space());

      // Initialize to 1
      ::Thyra::V_S(scale_vec.ptr(),1.0);

      // Create PPO
      std::string when_to_scale = nox_parameters.sublist("Solver Options").get("Update Row Sum Scaling", "Before Each Nolinear Iteration");
      if (when_to_scale == "Before Each Nolinear Iteration")
	scaled_ppo = rcp(new NOX::RowSumScaling(scale_vec, NOX::RowSumScaling::NOX_RowSumScaleAtBeginningOfIteration));
      else if (when_to_scale == "Before Each Nolinear Solve")
	scaled_ppo = rcp(new NOX::RowSumScaling(scale_vec, NOX::RowSumScaling::NOX_RowSumScaleAtBeginningOfSolve));
      else 
	TEST_FOR_EXCEPTION(true, std::logic_error,
			   "Error the option for \"Update Row Sum Scaling\" input as \"" << when_to_scale << "\" is invalid.  Please use one of the following: \"Before Each Nolinear Iteration\" or \"Before Each Nolinear Solve\".");

      // Create thyra scaled model evaluator
      scaled_model = rcp(new ::Thyra::ScaledModelEvaluator<double>);
      scaled_model->initialize(application_model);
      scaled_model->set_f_scaling(scale_vec);
    }

    void registerScaledPPOWithNOXParams(Teuchos::ParameterList& nox_parameters,
					const Teuchos::RCP<NOX::Abstract::PrePostOperator>& scaled_ppo)
    {
      if (nox_parameters.sublist("Solver Options").
	  isType<Teuchos::RCP<NOX::Abstract::PrePostOperator> >("User Defined Pre/Post Operator")) {
	    
	Teuchos::RCP<NOX::PrePostOperatorVector> ppo_vec = Teuchos::rcp(new NOX::PrePostOperatorVector);
	Teuchos::RCP<NOX::Abstract::PrePostOperator> user_ppo = 
	  nox_parameters.sublist("Solver Options").get<Teuchos::RCP<NOX::Abstract::PrePostOperator> >("User Defined Pre/Post Operator");

	ppo_vec->pushBack(scaled_ppo);
	ppo_vec->pushBack(user_ppo);

	nox_parameters.sublist("Solver Options").set("User Defined Pre/Post Operator", ppo_vec);
      }
      else
	nox_parameters.sublist("Solver Options").set("User Defined Pre/Post Operator", scaled_ppo);

    }

    template<typename Scalar>
    void setupProblemForRowSumScaling(const Teuchos::RCP<const ::Thyra::ModelEvaluator<Scalar> >& application_model,
				      Teuchos::RCP< ::Thyra::ScaledModelEvaluator<Scalar> >& scaled_model,
				      Teuchos::ParameterList& nox_parameters)
    {
      Teuchos::RCP<NOX::Abstract::PrePostOperator> scaled_ppo;
      buildRowSumScaledObjects(application_model, scaled_model, scaled_ppo, nox_parameters);
      registerScaledPPOWithNOXParams(nox_parameters, scaled_ppo);
    }

  }
}

#endif
