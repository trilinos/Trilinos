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

#ifndef LOCA_ANASAZIOPERATOR_ABSTRACTSTRATEGY_H
#define LOCA_ANASAZIOPERATOR_ABSTRACTSTRATEGY_H

#include "NOX_Abstract_Group.H" // for ReturnType

namespace LOCA {

  //! Various group based operators for computing eigenvalues with Anasazi
  namespace AnasaziOperator {
    
    //! Abstract interface class for Anasazi operator strategies
    /*!
     * AbstractStrategy defines an abstract interface for anasazi operators. 
     * It is used by LOCA::Eigensolver::AnasaziStrategy to compute 
     * different kinds of eigenvalues of 
     * the steady-state solution after each continuation step.  
     *
     * The interface currently defines several pure virtual methods, 
     * apply(), to apply the operator, transformEigenvalues() to transform
     * the computed eigenvalues back to eigenvalues of untransformed state,
     * rayleighQuotient to compute the rayleighQuotient for the operator, and
     * label() to return the name of the operator.
     * Derived classes should implement these method for a particular operator.
     * Constructors for derived classes should be of the form:
     *
     * \code
     * class Derived : public AbstractStrategy {
     * public:
     *   Derived(
     *   const Teuchos::RCP<LOCA::GlobalData>& global_data,
     *   const Teuchos::RCP<LOCA::Parameter::SublistParser>& topParams,
     *   const Teuchos::RCP<Teuchos::ParameterList>& eigenParams,
     *   const Teuchos::RCP<Teuchos::ParameterList>& solverParams,
     *   const Teuchos::RCP<NOX::Abstract::Group>& grp);
     *   ...
     * };
     * \endcode
     *
     * where \c global_data is the LOCA global data object, \c topParams is
     * the parsed top-level parameter list, \c eigenParams 
     * is a parameter list of eigensolver parameters, \c solverParams is a 
     * parameter list of linear solver parameters, and \c grp is the 
     * group representing the Jacobian and mass matrices.
     *
     * This class and its children follow the Strategy pattern as defined
     * in Erich Gamma, et al. "Design Patterns:  Elements of Reusable 
     * Object-Oriented Software." Addison Wesley, Boston, MA, 1995.
     */
    class AbstractStrategy {

    public:

      //! Constructor
      AbstractStrategy() {}

      //! Destructor
      virtual ~AbstractStrategy() {}

      /*!
       * \brief Return name of this operator
       */
      virtual const string& label() const = 0;

      /*! 
       * \brief Apply the operator to \c input with the result in \c output
       */
      virtual void
      apply(const NOX::Abstract::MultiVector& input,
	    NOX::Abstract::MultiVector& output) const = 0;

      /*! 
       * \brief Give strategy an opportunit to massage the random seed vector
       */
      virtual void
      preProcessSeedVector(NOX::Abstract::MultiVector& ivec) {};

      /*! 
       * \brief Transform eigenvalue in place
       */
      virtual void
      transformEigenvalue(double& ev_r, double& ev_i) const = 0;

      /*!
       * \brief Compute Rayleigh quotient
       */
      virtual NOX::Abstract::Group::ReturnType 
      rayleighQuotient(NOX::Abstract::Vector& evec_r,
		       NOX::Abstract::Vector& evec_i,
		       double& rq_r, double& rq_i) const = 0;

    private:

      //! Private to prohibit copying
      AbstractStrategy(const AbstractStrategy&);

      //! Private to prohibit copying
      AbstractStrategy& operator = (const AbstractStrategy&);

    }; // Class AbstractStrategy

  } // Namespace AnasaziOperator

} // Namespace LOCA

#endif // LOCA_ANASAZIOPERATOR_ABSTRACTSTRATEGY_H
