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
//@HEADER

#ifndef LOCA_EPETRA_ANASAZIOPERATOR_FLOQUET_H
#define LOCA_EPETRA_ANASAZIOPERATOR_FLOQUET_H

#include "Teuchos_RCP.hpp" 

#include "LOCA_AnasaziOperator_AbstractStrategy.H"  // base class
#include "LOCA_TimeDependent_AbstractGroup.H"
#include "LOCA_Epetra_Interface_xyzt.H"
#include "NOX_Epetra_Group.H"

// forward declarations
namespace LOCA {
  class GlobalData;
  namespace Parameter {
    class SublistParser;
  }
}

namespace LOCA {

  namespace Epetra {
    namespace AnasaziOperator {
    
    /*! \brief Anasazi operator for computing generalized eigenvalues 
     * using Cayley transformations.
     */
    /*! 
     * This class implements the LOCA::AnasaziOperator::AbstractStrategy 
     * interface for computing generalized eigenvalues \f$\lambda\f$ and 
     * eigenvectors \f$z\f$ of the system
     *\f[
     *    J z = \lambda M z
     *\f]
     * where \f$J\f$ is the Jacobian matrix and \f$M\f$ is the mass matrix.
     * The  eigenvalues are computed using a Cayley transformation, i.e.
     * solving
     * \f[
     *     (J - \sigma M) z = (J - \mu M) r
     * \f]
     * where \f$\sigma\f$ is the Cayley pole and \f$\mu\f$ is the Cayley
     * zero.  
     *
     * The parameters used by this class supplied in the constructor are:
     * <ul>
     * <li> "Cayley Pole" - \f$\sigma\f$ as defined above (Default 0.0)
     * <li> "Cayley Zero" - \f$\mu\f$ as defined above (Default 0.0)
     * </ul>
     *
     * Also the \em grp argument to the constructor must
     * be a child of LOCA::TimeDependent::AbstractGroup for the shift-invert
     * operations.
     */
    class Floquet : public LOCA::AnasaziOperator::AbstractStrategy {

    public:

      //! Constructor. 
      /*!
       * Argument \em grp must be of type LOCA::TimeDependent::AbstractGroup.
       * See class description for a list of eigenParams.
       */
      Floquet(
	const Teuchos::RCP<LOCA::GlobalData>& global_data,
	const Teuchos::RCP<LOCA::Parameter::SublistParser>& topParams,
	const Teuchos::RCP<Teuchos::ParameterList>& eigenParams,
	const Teuchos::RCP<Teuchos::ParameterList>& solverParams,
	const Teuchos::RCP<NOX::Abstract::Group>& grp);

      //! Destructor
      virtual ~Floquet();

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
      string myLabel;

      //! Stores parameters relating to the operator
      Teuchos::RCP<Teuchos::ParameterList> eigenParams;

      //! Stores linear solver parameters
      Teuchos::RCP<Teuchos::ParameterList> solverParams;

      //! Stores group representing Jacobian and Mass matrix
      Teuchos::RCP<NOX::Abstract::Group> grp;

      //! Stores interface for space-time periodic formulation
      Teuchos::RCP<LOCA::Epetra::Interface::xyzt> xyztInterface;

    };
  } // namespace AnasaziOperator
  } // namespace Epetra
} // namespace LOCA

#endif
