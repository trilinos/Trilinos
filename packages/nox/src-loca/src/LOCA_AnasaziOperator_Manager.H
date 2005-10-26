// $Id$ 
// $Source$ 

//@HEADER
// ************************************************************************
// 
//                  LOCA Continuation Algorithm Package
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
// Questions? Contact Tammy Kolda (tgkolda@sandia.gov) or Roger Pawlowski
// (rppawlo@sandia.gov), Sandia National Laboratories.
// 
// ************************************************************************
//@HEADER

#ifndef LOCA_ANASAZIOPERATOR_MANAGER_H
#define LOCA_ANASAZIOPERATOR_MANAGER_H

#include "LOCA_AnasaziOperator_Generic.H"  // base class

namespace LOCA {

  namespace AnasaziOperator {
    
    //! Manager for all Anasazi operators
    /*! 
     * The parameters passed to the constructor or reset should specify the
     * "Operator", as described below, as well as any additional parameters
     * for that particular operator.
     * <ul>
     * <li> "Operator" - Name of Anasazi operator to use.  Valid choices are:
     *     <ul>
     *     <li> "%Jacobian Inverse" (LOCA::AnasaziOperator::JacobianInverse) 
     *          [Default]
     *     <li> "%Shift-Invert" (LOCA::AnasaziOperator::ShiftInvert)
     *     <li> "%Cayley" (LOCA::AnasaziOperator::Cayley)
     *     </ul>
     * </ul>
     */
    class Manager : public LOCA::AnasaziOperator::Generic {

    public:

      //! Constructor. 
      Manager(const Teuchos::RefCountPtr<NOX::Parameter::List>& eigenParams,
	      const Teuchos::RefCountPtr<NOX::Parameter::List>& solverParams,
	      const Teuchos::RefCountPtr<NOX::Abstract::Group>& grp);

      //! Destructor
      virtual ~Manager();

      //! Reset parameters
      virtual NOX::Abstract::Group::ReturnType 
      reset(const Teuchos::RefCountPtr<NOX::Parameter::List>& eigenParams,
	    const Teuchos::RefCountPtr<NOX::Parameter::List>& solverParams,
	    const Teuchos::RefCountPtr<NOX::Abstract::Group>& grp);

      /*!
       * \brief Return name of this operator
       */
      virtual const string& label() const;
  
      /*! 
       * \brief Apply the operator.
       */
      virtual NOX::Abstract::Group::ReturnType 
      apply(const NOX::Abstract::Vector& input, 
	    NOX::Abstract::Vector& output) const;

      /*! 
       * \brief Transform eigenvalue.
       */
      virtual void
      transformEigenvalue(double& ev_r, double& ev_i) const;

      /*!
       * \brief Compute Rayleigh quotient
       */
      virtual NOX::Abstract::Group::ReturnType 
      rayleighQuotient(const NOX::Abstract::Vector& evec_r,
		       const NOX::Abstract::Vector& evec_i,
		       double& rq_r, double& rq_i) const;

      //! Returns underlying operator
      LOCA::AnasaziOperator::Generic& getOperator();

      //! Returns underlying operator
      const LOCA::AnasaziOperator::Generic& getOperator() const;

    protected:

      //! Stores name of requested operator
      string opName;

      //! Stores operator
      LOCA::AnasaziOperator::Generic *op;

    };
  } // namespace AnasaziOperator
} // namespace LOCA

#endif
