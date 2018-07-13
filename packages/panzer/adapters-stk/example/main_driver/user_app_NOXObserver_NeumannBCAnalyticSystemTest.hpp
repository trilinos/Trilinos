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

#ifndef USER_APP_NOX_OBSERVER_NEUMANN_BC_ANALYTIC_SYSTEM_TEST_HPP
#define USER_APP_NOX_OBSERVER_NEUMANN_BC_ANALYTIC_SYSTEM_TEST_HPP

#include "PanzerAdaptersSTK_config.hpp"

#include "NOX_Abstract_PrePostOperator.H"

#include "Teuchos_RCP.hpp"
#include "NOX_Thyra_Vector.H"
#include "Thyra_DefaultProductVector.hpp"

namespace user_app {
  
  class NOXObserver_NeumannBCAnalyticSystemTest : public NOX::Abstract::PrePostOperator {
    
  public:
    
    NOXObserver_NeumannBCAnalyticSystemTest()
    { 

    }
      
    void runPreIterate(const NOX::Solver::Generic& /* solver */)
    {

    }
    
    void runPostIterate(const NOX::Solver::Generic& /* solver */)
    {

    }
    
    void runPreSolve(const NOX::Solver::Generic& /* solver */)
    {

    }
    
    void runPostSolve(const NOX::Solver::Generic& solver)
    {
      /*
        The Neumann BC Analytic System Test is simple heat conduction
        in the x-direction.  By fixing the flux and thermal
        conductivity, there is an analytic solution with a value of
        2.0 for temperature on the left surface where the neumann
        condition is applied.  We will check exactly one node on this
        boundary and make sure it is within the expected error
        (accounting for discretization error).  

	IMPORTANT: If the global dof numbering scheme changes in the
	future, this test could fail by grabbing the wrong node.
      */

      const NOX::Abstract::Vector& x = solver.getSolutionGroup().getX();
      const NOX::Thyra::Vector* n_th_x = dynamic_cast<const NOX::Thyra::Vector*>(&x);
      TEUCHOS_TEST_FOR_EXCEPTION(n_th_x == NULL, std::runtime_error, "Failed to dynamic_cast to NOX::Thyra::Vector!")
      const ::Thyra::VectorBase<double>& th_x = n_th_x->getThyraVector(); 

      // Assume a single process job (enforced in CMakeLists.txt file)

      const ::Thyra::SpmdVectorBase<double>* th_spmd_x = dynamic_cast<const ::Thyra::SpmdVectorBase<double>* >(&th_x);

      TEUCHOS_TEST_FOR_EXCEPTION(th_spmd_x == NULL, std::runtime_error, "Failed to dynamic_cast to Thyra SPMD vector!");

      Teuchos::ArrayRCP<const double> local_values;
     
      th_spmd_x->getLocalData(outArg(local_values));

//       for (std::size_t i=0; i < local_values.size(); ++i)
// 	std::cout << "local_values[" << i << "] = " << local_values[i] << std::endl;

      // linear solve tolerance is 1e-10 and this is a linear problem with a linear basis
      double tol = 1.0e-9;

      TEUCHOS_TEST_FOR_EXCEPTION( std::fabs(local_values[0] - 2.0) > tol, std::runtime_error, "Solution value for Neumann Condition is " << local_values[0] << " and should be 2.0.  The tolerance condition, std::fabs(local_values[0] - 2.0) > tol, where tol is " << tol << " has been violated, causing the test to fail.");

      local_values = Teuchos::null;
   
    }
    
  };
}

#endif
