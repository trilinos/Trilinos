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

#ifndef LOCA_SINGULARJACOBIANSOLVE_DEFAULT_H
#define LOCA_SINGULARJACOBIANSOLVE_DEFAULT_H

#include "LOCA_SingularJacobianSolve_Generic.H"  // base class

namespace LOCA {

  namespace SingularJacobianSolve {

    //! %Default singular Jacobian solve computation class
    /*! 
     * This class computes the solution to \f$J x = b\f$ using the 
     * applyJacobianInverse method of the underlying group ignoring the
     * null vector data.
    */
    class Default : public LOCA::SingularJacobianSolve::Generic {

    public:

      //! Constructor. 
      Default(NOX::Parameter::List& params);

      //! Copy constructor
      Default(const Default& source);

      //! Destructor
      virtual ~Default();

      //! Clone function
      virtual LOCA::SingularJacobianSolve::Generic* clone() const;

      //! Assignment operator
      virtual LOCA::SingularJacobianSolve::Generic& 
      operator = (const LOCA::SingularJacobianSolve::Generic& source);

      //! Assignment operator
      virtual Default&
      operator = (const Default& source);

      //! Reset parameters
      /*!
	There are no additional parameters for the default calculation.
      */
      virtual NOX::Abstract::Group::ReturnType 
      reset(NOX::Parameter::List& params);
  
      //! Computes the solution as described above.
      virtual NOX::Abstract::Group::ReturnType 
      compute(NOX::Parameter::List& params,
	      LOCA::Continuation::AbstractGroup& grp,
	      const NOX::Abstract::Vector& input,
	      const NOX::Abstract::Vector& approxNullVec,
	      const NOX::Abstract::Vector& jacApproxNullVec,
	      NOX::Abstract::Vector& result);

      //! Computes solution for multiple RHS using applyJacobianInverseMulti
      virtual NOX::Abstract::Group::ReturnType 
      computeMulti(NOX::Parameter::List& params,
		   LOCA::Continuation::AbstractGroup& grp,
		   const NOX::Abstract::Vector*const* inputs,
		   const NOX::Abstract::Vector& approxNullVec,
		   const NOX::Abstract::Vector& jacApproxNullVec,
		   NOX::Abstract::Vector** results,
		   int nVecs);

    };
  } // namespace SingularJacobianSolve
} // namespace LOCA

#endif
