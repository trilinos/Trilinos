// $Id$ 
// $Source$ 

//@HEADER
// ************************************************************************
// 
//            NOX: An Object-Oriented Nonlinear Solver Package
//                 Copyright (2002) Sandia Corporation
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
// Questions? Contact Tammy Kolda (tgkolda@sandia.gov) or Roger Pawlowski
// (rppawlo@sandia.gov), Sandia National Laboratories.
// 
// ************************************************************************
//@HEADER

#ifndef LOCA_SINGULARJACOBIANSOLVE_MANAGER_H
#define LOCA_SINGULARJACOBIANSOLVE_MANAGER_H

#include "LOCA_SingularJacobianSolve_Generic.H"  // base class
#include "NOX_Parameter_List.H"

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
      <li> "%Default" (LOCA::SingularJacobianSolve::Default) [Default]
      <li> "%Nic" (LOCA::SingularJacobianSolve::Nic) 
      <li> "%Nic-Day" (LOCA::SingularJacobianSolve::NicDay) 
      <li> "%Iterative Refinement" (LOCA::SingularJacobianSolve::ItRef) 
      </ul>
      </ul>
    */  
    class Manager : public LOCA::SingularJacobianSolve::Generic {

    public:

      //! Constructor. 
      Manager(NOX::Parameter::List& params);

      //! Constructor. 
      Manager(const NOX::Parameter::List& params = NOX::Parameter::List());

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
      reset(NOX::Parameter::List& params);
  
      //! Computes solution based on method parameter
      virtual NOX::Abstract::Group::ReturnType 
      compute(NOX::Parameter::List& params,
	      LOCA::Continuation::AbstractGroup& grp,
	      const NOX::Abstract::Vector& input,
	      const NOX::Abstract::Vector& approxNullVec,
	      const NOX::Abstract::Vector& jacApproxNullVec,
	      NOX::Abstract::Vector& result);

      //! Computes solution based on method parameter for multiple RHS
      virtual NOX::Abstract::Group::ReturnType 
      computeMulti(NOX::Parameter::List& params,
		   LOCA::Continuation::AbstractGroup& grp,
		   const NOX::Abstract::Vector*const* inputs,
		   const NOX::Abstract::Vector& approxNullVec,
		   const NOX::Abstract::Vector& jacApproxNullVec,
		   NOX::Abstract::Vector** results,
		   int nVecs);

      //! Returns current method
      const string& getMethod() const;

    protected:
      //! Stores predictor method name
      string method;

      //! Stores pointer to current singular solver object
      LOCA::SingularJacobianSolve::Generic *singularSolverPtr;

    };
  } // namespace SingularJacobianSolve
} // namespace LOCA

#endif
