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

#ifndef LOCA_ANASAZIOPERATOR_GENERIC_H
#define LOCA_ANASAZIOPERATOR_GENERIC_H

#include "NOX_Abstract_Group.H"  // for ReturnType, Parameter::List

namespace LOCA {

  //! %Various group based operators for computing eigenvalues with Anasazi
  namespace AnasaziOperator {
    
    //! %Generic Anasazi operator interface
    /*! 
     * %Generic interface for computing eigenvalues for some operator using
     * Anasazi.  In additions, this class provides an interface for
     * transforming the eigenvalues as defined by the operator, and computing
     * the Rayleigh quotient.
     */
    class Generic {

    public:

      //! Constructor. 
      /*! 
       * Constructors of derived objects should look like reset.
       */
      Generic() {};

      //! Destructor
      virtual ~Generic() {};

      //! Reset parameters
      virtual NOX::Abstract::Group::ReturnType 
      reset(NOX::Parameter::List& eigenParams,
	    NOX::Parameter::List& solverParams,
	    NOX::Abstract::Group& grp) = 0;

      /*!
       * \brief Return name of this operator
       */
      virtual const string& label() const = 0;
  
      /*! 
       * \brief Apply the operator.
       */
      virtual NOX::Abstract::Group::ReturnType 
      apply(const NOX::Abstract::Vector& input, 
	    NOX::Abstract::Vector& output) const = 0;

      /*! 
       * \brief Transform eigenvalue.
       */
      virtual void
      transformEigenvalue(double& ev_r, double& ev_i) const = 0;

      /*!
       * \brief Compute Rayleigh quotient.
       */
      virtual NOX::Abstract::Group::ReturnType 
      rayleighQuotient(const NOX::Abstract::Vector& evec_r,
		       const NOX::Abstract::Vector& evec_i,
		       double& rq_r, double& rq_i) const = 0;

    };
  } // namespace AnasaziOperator
} // namespace LOCA

#endif
