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

#ifndef LOCA_ANASAZIOPERATOR_SHIFTINVERT_H
#define LOCA_ANASAZIOPERATOR_SHIFTINVERT_H

#include "Teuchos_RCP.hpp" 

#include "LOCA_AnasaziOperator_AbstractStrategy.H"  // base class
#include "LOCA_TimeDependent_AbstractGroup.H"

// forward declarations
namespace LOCA {
  class GlobalData;
  namespace Parameter {
    class SublistParser;
  }
}

namespace LOCA {

  namespace AnasaziOperator {
    
    /*! \brief Anasazi operator for computing generalized eigenvalues 
     * using shift-invert.
     */
    /*! 
     * This class implements the LOCA::AnasaziOperator::AbstractStrategy 
     * interface for computing generalized eigenvalues \f$\lambda\f$ and 
     * eigenvectors \f$z\f$ of the system
     *\f[
     *    J z = \lambda M z
     *\f]
     * where \f$J\f$ is the Jacobian matrix and \f$M\f$ is the mass matrix.
     * The right-most eigenvalues are computed using shift-invert, i.e.
     * solving
     * \f[
     *     (J - \omega M) z = \lambda M z
     * \f]
     * where \f$\omega\f$ is a real scalar.  The resulting eigenvalue is
     * \f$\lambda + \omega\f$.
     *
     * The parameters used by this class supplied in the constructor are:
     * <ul>
     * <li> "Shift" - \f$\omega\f$ as defined above (Default 0.0)
     * </ul>
     *
     * Also the \em grp argument to the constructor must
     * be a child of LOCA::TimeDependent::AbstractGroup for the shift-invert
     * operations.
     */
    class ShiftInvert : public LOCA::AnasaziOperator::AbstractStrategy {

    public:

      //! Constructor. 
      /*!
       * Argument \em grp must be of type LOCA::TimeDependent::AbstractGroup.
       * See class description for a list of eigenParams.
       */
      ShiftInvert(
	const Teuchos::RCP<LOCA::GlobalData>& global_data,
	const Teuchos::RCP<LOCA::Parameter::SublistParser>& topParams,
	const Teuchos::RCP<Teuchos::ParameterList>& eigenParams,
	const Teuchos::RCP<Teuchos::ParameterList>& solverParams,
	const Teuchos::RCP<LOCA::TimeDependent::AbstractGroup>& grp);

      //! Destructor
      virtual ~ShiftInvert();

      /*!
       * \brief Return name of this operator
       */
      virtual const std::string& label() const;
  
      /*! 
       * \brief Apply the operator.
       *
       * Applies the inverse of the shifted operator, i.e., solves
       * \f[
       *      (J-\omega I)z = M r
       * \f]
       * for \f$z\f$, where \f$r = \mbox{input}\f$ and \f$z = \mbox{output}\f$.
       */
      virtual void
      apply(const NOX::Abstract::MultiVector& input, 
	    NOX::Abstract::MultiVector& output) const;

      /*! 
       * \brief Transform eigenvalue.
       *
       * Transforms the given eigenvalue to the eigenvalue of the 
       * Jacobian-mass matrix system by shifting and inverting it.
       */
      virtual void
      transformEigenvalue(double& ev_r, double& ev_i) const;

      /*!
       * \brief Compute Rayleigh quotient
       *
       * Computes the Rayleigh quotient \f$z^T J z / z^T M z\f$ for
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
      mutable Teuchos::RCP<Teuchos::ParameterList> eigenParams;

      //! Stores linear solver parameters
      Teuchos::RCP<Teuchos::ParameterList> solverParams;

      //! Stores group representing Jacobian and Mass matrix
      Teuchos::RCP<LOCA::TimeDependent::AbstractGroup> grp;

      //! Stores a temporary vector for computing Rayleigh quotients
      mutable Teuchos::RCP<NOX::Abstract::MultiVector> tmp_r;

      //! Stores a temporary vector for computing Rayleigh quotients
      mutable Teuchos::RCP<NOX::Abstract::MultiVector> tmp_i;

      //! Stores shift value
      double shift;

    };
  } // namespace AnasaziOperator
} // namespace LOCA

#endif
