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

#ifndef LOCA_ANASAZIOPERATOR_SHIFTINVERT_H
#define LOCA_ANASAZIOPERATOR_SHIFTINVERT_H

#include "LOCA_AnasaziOperator_Generic.H"  // base class
#include "LOCA_TimeDependent_AbstractGroup.H"

namespace LOCA {

  namespace AnasaziOperator {
    
    /*! \brief Anasazi operator for computing generalized eigenvalues 
     * using shift-invert.
     */
    /*! 
     * This class implements the LOCA::AnasaziOperator::Generic interface
     * for computing generalized eigenvalues \f$\lambda\f$ and eigenvectors
     * \f$z\f$ of the system
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
     * The parameters used by this class supplied in the constructor or reset
     * method are:
     * <ul>
     * <li> "Shift" - \f$\omega\f$ as defined above (Default 0.0)
     * </ul>
     *
     * Also the \em grp argument to the constructor and reset method must
     * be a child of LOCA::TimeDependent::AbstractGroup for the shift-invert
     * operations.
     */
    class ShiftInvert : public LOCA::AnasaziOperator::Generic {

    public:

      //! Constructor. 
      /*!
       * Argument \em grp must be of type LOCA::TimeDependent::AbstractGroup.
       * See class description for a list of eigenParams.
       */
      ShiftInvert(NOX::Parameter::List& eigenParams,
		  NOX::Parameter::List& solverParams,
		  NOX::Abstract::Group& grp);

      //! Destructor
      virtual ~ShiftInvert();

      //! Reset parameters
      /*!
       * Argument \em grp must be of type LOCA::TimeDependent::AbstractGroup.
       * See class description for a list of eigenParams.
       */
      virtual NOX::Abstract::Group::ReturnType 
      reset(NOX::Parameter::List& eigenParams,
	    NOX::Parameter::List& solverParams,
	    NOX::Abstract::Group& grp);

      /*!
       * \brief Return name of this operator
       */
      virtual const string& label() const;
  
      /*! 
       * \brief Apply the operator.
       *
       * Applies the inverse of the shifted operator, i.e., solves
       * \f[
       *      (J-\omega I)z = M r
       * \f]
       * for \f$z\f$, where \f$r = \mbox{input}\f$ and \f$z = \mbox{output}\f$.
       */
      virtual NOX::Abstract::Group::ReturnType 
      apply(const NOX::Abstract::Vector& input, 
	    NOX::Abstract::Vector& output) const;

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
      rayleighQuotient(const NOX::Abstract::Vector& evec_r,
		       const NOX::Abstract::Vector& evec_i,
		       double& rq_r, double& rq_i) const;

    protected:

      //! Name of this operator
      string myLabel;

      //! Stores parameters relating to the operator
      NOX::Parameter::List* eigenParams;

      //! Stores linear solver parameters
      NOX::Parameter::List* solverParams;

      //! Stores group representing Jacobian and Mass matrix
      LOCA::TimeDependent::AbstractGroup* grp;

      //! Stores a temporary vector for computing Rayleigh quotients
      mutable NOX::Abstract::Vector* tmp_r;

      //! Stores a temporary vector for computing Rayleigh quotients
      mutable NOX::Abstract::Vector* tmp_i;

      //! Stores shift value
      double shift;

    };
  } // namespace AnasaziOperator
} // namespace LOCA

#endif
