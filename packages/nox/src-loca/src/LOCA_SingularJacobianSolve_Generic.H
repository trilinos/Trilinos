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

#ifndef LOCA_SINGULARJACOBIANSOLVE_GENERIC_H
#define LOCA_SINGULARJACOBIANSOLVE_GENERIC_H

#include "NOX_Abstract_Group.H"  // for ReturnType, Parameter::List

// Forward declarations
namespace LOCA {
  namespace Continuation {
    class AbstractGroup;
  }
}

namespace LOCA {

  //! %Strategies for solving \f$Jx=b\f$ when \f$J\f$ is (nearly) singular
  namespace SingularJacobianSolve {
    
    //! %Generic singular jacobian solve interface.
    /*! 
     * %Generic interface for solving \f$Jx=b\f$ when \f$J\f$ is (nearly) 
     * singular.  All classes the implement a method for computing solutions
     * to nearly singular systems should be derived from this class.
    */
    class Generic {

    public:

      //! Constructor. 
      /*! 
	Constructors of derived objects should look like reset.
      */
      Generic() {};

      //! Copy constructor
      Generic(const Generic& source) {}

      //! Destructor
      virtual ~Generic() {};

      //! Clone function
      virtual Generic* clone() const = 0;

      //! Assignment operator
      virtual Generic& operator = (const Generic& source) = 0;

      //! Reset parameters
      virtual NOX::Abstract::Group::ReturnType 
      reset(NOX::Parameter::List& params) = 0;
  
      /*! 
	\brief Compute solution to singular system.

	The passed parameters are assumed be the (nonsingular) linear solver 
	parameters.
      */
      virtual NOX::Abstract::Group::ReturnType 
      compute(NOX::Parameter::List& params,
	      LOCA::Continuation::AbstractGroup& grp,
	      const NOX::Abstract::Vector& input,
	      const NOX::Abstract::Vector& approxNullVec,
	      const NOX::Abstract::Vector& jacApproxNullVec,
	      NOX::Abstract::Vector& result) = 0;

      /*! 
	\brief Compute solution to singular system with multiple RHS.
	
	The passed parameters are assumed be the (nonsingular) linear solver 
	parameters.
      */
      virtual NOX::Abstract::Group::ReturnType 
      computeMulti(NOX::Parameter::List& params,
		   LOCA::Continuation::AbstractGroup& grp,
		   const NOX::Abstract::Vector*const* inputs,
		   const NOX::Abstract::Vector& approxNullVec,
		   const NOX::Abstract::Vector& jacApproxNullVec,
		   NOX::Abstract::Vector** results,
		   int nVecs) = 0;

    };
  } // namespace SingularJacobianSolve
} // namespace LOCA

#endif
