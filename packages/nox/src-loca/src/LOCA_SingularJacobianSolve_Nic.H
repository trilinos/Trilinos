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
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
// any later version.
//   
// This program is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// General Public License for more details.
//   
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
// 
// Questions? Contact Tammy Kolda (tgkolda@sandia.gov) or Roger Pawlowski
// (rppawlo@sandia.gov).
// 
// ************************************************************************
//@HEADER

#ifndef LOCA_SINGULARJACOBIANSOLVE_NIC_H
#define LOCA_SINGULARJACOBIANSOLVE_NIC_H

#include "LOCA_SingularJacobianSolve_Generic.H"  // base class

namespace LOCA {

  namespace SingularJacobianSolve {

    //! This class computes the solution to \f$J x = b\f$ using the Nic method.
    class Nic : public LOCA::SingularJacobianSolve::Generic {

    public:

      //! Constructor. 
      Nic(NOX::Parameter::List& params);

      //! Copy constructor
      Nic(const Nic& source);

      //! Destructor
      virtual ~Nic();

      //! Clone function
      virtual LOCA::SingularJacobianSolve::Generic* clone() const;

      //! Assignment operator
      virtual LOCA::SingularJacobianSolve::Generic& 
      operator = (const LOCA::SingularJacobianSolve::Generic& source);

      //! Assignment operator
      virtual Nic&
      operator = (const Nic& source);

      //! Reset parameters
      /*!
	There are no additional parameters for the Nic calculation.
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

      //! Computes solution for multiple RHS 
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
