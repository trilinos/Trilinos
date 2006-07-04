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
// Questions? Contact Andy Salinger (agsalin@sandia.gov) or Eric Phipps
// (etphipp@sandia.gov), Sandia National Laboratories.
//
// ************************************************************************
//@HEADER

#ifndef LOCA_EIGENSOLVER_ABSTRACTSTRATEGY_H
#define LOCA_EIGENSOLVER_ABSTRACTSTRATEGY_H

#include "Teuchos_RefCountPtr.hpp" 

#include "NOX_Abstract_Group.H" // for ReturnType

namespace LOCA {

  //! Collection of eigensolver strategies
  namespace Eigensolver {
    
    //! Abstract interface class for Eigensolver strategies
    /*!
     * AbstractStrategy defines an abstract interface for eigensolver 
     * strategies. It is used by LOCA::Stepper to compute eigenvalues of 
     * the steady-state solution after each continuation step.  
     *
     * The interface currently defines one pure virtual method, 
     * computeEigenvalues(), to compute the eigenvalues.  Derived classes
     * should implement this method for a particular eigensolver strategy.
     * Constructors for derived classes should be of the form:
     *
     * \code
     * class Derived : public AbstractStrategy {
     * public:
     *   Derived(
     *   const Teuchos::RefCountPtr<LOCA::GlobalData>& global_data,
     *   const Teuchos::RefCountPtr<LOCA::Parameter::SublistParser>& topParams,
	 const Teuchos::RefCountPtr<Teuchos::ParameterList>& eigenParams);
     *   ...
     * };
     * \endcode
     *
     * where \c global_data is the LOCA global data object, \c topParams is
     * the parsed top-level parameter list, and \c eigenParams 
     * is a parameter list of eigensolver parameters.
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

      //! Compute eigenvalues/eigenvectors in group \c group
      /*!
       * \param group [in] NOX Group to compute eigenvalues of
       * \param evals_r [out] Real eigenvalues
       * \param evals_i [out] Imaginary eigenvalues
       * \param evecs_r [out] Real eigenvectors
       * \param evecs_i [out] Imaginary eigenvectors
       * \returns ReturnType code indicating success or failure
       */
      virtual NOX::Abstract::Group::ReturnType 
      computeEigenvalues(
	      NOX::Abstract::Group& group,
	      Teuchos::RefCountPtr< std::vector<double> >& evals_r,
	      Teuchos::RefCountPtr< std::vector<double> >& evals_i,
	      Teuchos::RefCountPtr< NOX::Abstract::MultiVector >& evecs_r,
	      Teuchos::RefCountPtr< NOX::Abstract::MultiVector >& evecs_i) = 0;

    private:

      //! Private to prohibit copying
      AbstractStrategy(const AbstractStrategy&);

      //! Private to prohibit copying
      AbstractStrategy& operator = (const AbstractStrategy&);

    }; // Class AbstractStrategy

  } // Namespace Eigensolver

} // Namespace LOCA

#endif // LOCA_EIGENSOLVER_ABSTRACTSTRATEGY_H
