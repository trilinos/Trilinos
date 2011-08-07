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

#ifndef LOCA_ANASAZIOPERATOR_JACOBIANINVERSE_H
#define LOCA_ANASAZIOPERATOR_JACOBIANINVERSE_H

#include "Teuchos_RCP.hpp" 

#include "LOCA_AnasaziOperator_AbstractStrategy.H"  // base class

// forward declarations
namespace LOCA {
  class GlobalData;
  namespace Parameter {
    class SublistParser;
  }
}

namespace LOCA {

  namespace AnasaziOperator {
    
    //! Anasazi operator for computing eigenvalues of the inverse-Jacobian
    /*! 
     * This class implements the LOCA::AnasaziOperator::AbstractStrategy 
     * interface for computing eigenvalues of the inverse-Jacobian.
     */
    class JacobianInverse : public LOCA::AnasaziOperator::AbstractStrategy {

    public:

      //! Constructor. 
      JacobianInverse(
	const Teuchos::RCP<LOCA::GlobalData>& global_data,
	const Teuchos::RCP<LOCA::Parameter::SublistParser>& topParams,
	const Teuchos::RCP<Teuchos::ParameterList>& eigenParams,
	const Teuchos::RCP<Teuchos::ParameterList>& solverParams,
	const Teuchos::RCP<NOX::Abstract::Group>& grp);

      //! Destructor
      virtual ~JacobianInverse();

      /*!
       * \brief Return name of this operator
       */
      virtual const string& label() const;
  
      /*! 
       * \brief Apply the operator.
       *
       * Computes \f$\mbox{output} = J^{-1}\mbox{input}\f$.
       */
      virtual void
      apply(const NOX::Abstract::MultiVector& input, 
	    NOX::Abstract::MultiVector& output) const;

      /*! 
       * \brief Transform eigenvalue.
       *
       * Transforms the given eigenvalue to the eigenvalue of the Jacobian
       * by inverting it.
       */
      virtual void
      transformEigenvalue(double& ev_r, double& ev_i) const;

      /*!
       * \brief Compute Rayleigh quotient
       *
       * Computes the Rayleigh quotient \f$z^T J z\f$ for
       * the eigenvector \f$z\f$.
       */
      virtual NOX::Abstract::Group::ReturnType 
      rayleighQuotient(NOX::Abstract::Vector& evec_r,
		       NOX::Abstract::Vector& evec_i,
		       double& rq_r, double& rq_i) const;

    protected:

      //! Global data
      Teuchos::RCP<LOCA::GlobalData> globalData;

      //! Name of this operator
      string myLabel;

      //! Stores parameters relating to the operator
      Teuchos::RCP<Teuchos::ParameterList> eigenParams;

      //! Stores linear solver parameters
      Teuchos::RCP<Teuchos::ParameterList> solverParams;

      //! Stores group representing Jacobian
      Teuchos::RCP<NOX::Abstract::Group> grp;

      //! Stores a temporary vector for computing Rayleigh quotients
      mutable Teuchos::RCP<NOX::Abstract::Vector> tmp_r;

      //! Stores a temporary vector for computing Rayleigh quotients
      mutable Teuchos::RCP<NOX::Abstract::Vector> tmp_i;

    };
  } // namespace AnasaziOperator
} // namespace LOCA

#endif
