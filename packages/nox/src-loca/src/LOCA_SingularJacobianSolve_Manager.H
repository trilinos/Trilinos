// $Id$ 
// $Source$ 

//@HEADER
// ************************************************************************
// 
//            LOCA: Library of Continuation Algorithms Package
//                 Copyright (2005) Sandia Corporation
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

#ifndef LOCA_SINGULARJACOBIANSOLVE_MANAGER_H
#define LOCA_SINGULARJACOBIANSOLVE_MANAGER_H

#include "LOCA_SingularJacobianSolve_Generic.H"  // base class
#include "Teuchos_ParameterList.hpp"

namespace LOCA {

  namespace SingularJacobianSolve {

    //! %Manager for all singular Jacobian solve computations
    /*! 
      The parameters passed to the constructor or reset should specify
      the "Method", as described below, as well as any additional
      parameters for that particular method.
 
      <ul>
      <li> "Method" - Name of the singular jacobian solve method. 
           Valid choices are
      <ul> 
      <li> "Default" (LOCA::SingularJacobianSolve::Default) [Default]
      <li> "Nic" (LOCA::SingularJacobianSolve::Nic) 
      <li> "Nic-Day" (LOCA::SingularJacobianSolve::NicDay) 
      <li> "Iterative Refinement" (LOCA::SingularJacobianSolve::ItRef) 
      </ul>
      </ul>
    */  
    class Manager : public LOCA::SingularJacobianSolve::Generic {

    public:

      //! Constructor. 
      Manager(Teuchos::ParameterList& params);

      //! Constructor. 
      Manager(const Teuchos::ParameterList& params = Teuchos::ParameterList());

      //! Copy constructor
      Manager(const Manager& source);

      //! Destructor
      virtual ~Manager();

      //! Clone function
      virtual LOCA::SingularJacobianSolve::Generic* clone() const;

      //! Assignment operator
      virtual LOCA::SingularJacobianSolve::Generic& 
      operator = (const LOCA::SingularJacobianSolve::Generic& source);

      //! Assignment operator
      virtual Manager&
      operator = (const Manager& source);

      //! Reset parameters
      virtual NOX::Abstract::Group::ReturnType 
      reset(Teuchos::ParameterList& params);
  
      //! Computes solution based on method parameter
      virtual NOX::Abstract::Group::ReturnType 
      compute(Teuchos::ParameterList& params,
	      LOCA::Continuation::AbstractGroup& grp,
	      const NOX::Abstract::Vector& input,
	      const NOX::Abstract::Vector& approxNullVec,
	      const NOX::Abstract::Vector& jacApproxNullVec,
	      NOX::Abstract::Vector& result);

      //! Computes solution based on method parameter for multiple RHS
      virtual NOX::Abstract::Group::ReturnType 
      computeMulti(Teuchos::ParameterList& params,
		   LOCA::Continuation::AbstractGroup& grp,
		   const NOX::Abstract::Vector*const* inputs,
		   const NOX::Abstract::Vector& approxNullVec,
		   const NOX::Abstract::Vector& jacApproxNullVec,
		   NOX::Abstract::Vector** results,
		   int nVecs);

      //! Returns current method
      const std::string& getMethod() const;

    protected:
      //! Stores predictor method name
      std::string method;

      //! Stores pointer to current singular solver object
      LOCA::SingularJacobianSolve::Generic *singularSolverPtr;

    };
  } // namespace SingularJacobianSolve
} // namespace LOCA

#endif
