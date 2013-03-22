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
      virtual const std::string& label() const;
  
      /*! 
       * \brief Apply the operator.
       *
       * Computes \f$\mbox{output} = J^{-1}\mbox{input}\f$.
       */
      virtual void
      apply(const NOX::Abstract::MultiVector& input, 
	    NOX::Abstract::MultiVector& output) const;

      /*! 
       * \brief  Begin PostProcessing of eigenvalues
       *
       * Compute Jacobian matrix once, for use in subsequent 
       * repeated calls to rayleighQuotient
       */
      virtual void
      beginPostProcessing();

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
      std::string myLabel;

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
