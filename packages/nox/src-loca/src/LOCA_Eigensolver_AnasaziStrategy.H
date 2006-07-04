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

#ifndef LOCA_EIGENSOLVER_ANASAZISTRATEGY_H
#define LOCA_EIGENSOLVER_ANASAZISTRATEGY_H

#include "LOCA_Eigensolver_AbstractStrategy.H" // base class

#include "Teuchos_ParameterList.hpp"           // class element
#include "AnasaziOutputManager.hpp"            // class element
#include "Anasazi_LOCA_Sort.H"                 // class element

// Forward declarations
namespace Teuchos {
  class ParameterList;
}
namespace LOCA {
  class GlobalData;
  namespace Parameter {
    class SublistParser;
  }
}

namespace LOCA {

  namespace Eigensolver {

    //! Anasazi eigensolver strategy
    /*!
     * This class implements an eigensolver strategy using the generic
     * Trilinos eigensolver package Anasazi.  Since Anasazi is a generic
     * solver, this strategy will work with any group implementation.  This
     * strategy references the following parameters passed through the
     * \c eigenParams argument to the constructor:
     * <ul>
     * <li> "Operator" -- [string] (default: "Jacobian inverse") Operator
     *      to compute eigenvalues of.
     * <li> "Block Size" -- [int] (default: 1) Block size
     * <li> "Arnoldi Size" -- [int] (default: 30) Maximum length of the 
     *       Arnoldi factorization
     * <li> "NEV" -- [int] (default: 4) Number of requested eigenvalues
     * <li> "Tol" -- [double] (default: 1.0e-7) Tolerance for the converged 
     *       eigenvalues
     * <li> "Convergence Check" -- [int] (default: 1) Checks convergence 
     *       every so many steps
     * <li> "Restarts" -- [int] (default: 1) Number of restarts allowed
     * <li> "Debug Level" -- [int] (default: 1) Anasazi Debug level
     * <li> "Sorting Order" -- [string" (default: "LM") Sorting order of 
     *       printed eigenvalues
     * <li> "Save Eigenvectors" -- [int] (default: 0) How many eigenvectors
     *      to save
     * <li> "Cayley Pole" -- [double] (default: 0.0) Cayley transformation pole
     * <li> "Cayley Zero" -- [double] (default: 0.0) Cayley tranformation zero
     * </ul>
     */
    class AnasaziStrategy : public LOCA::Eigensolver::AbstractStrategy {

    public:

      //! Constructor
      /*!
       * \param global_data [in] Global data object
       * \param topParams [in] Parsed top-level parameter list.  Used to
       * obtain linear-solver parameters and to instantiate sorting strategy.
       * \param eigenParams [in] %Eigensolver parameters as described above.
       * solver.
       */
      AnasaziStrategy(
	const Teuchos::RefCountPtr<LOCA::GlobalData>& global_data,
	const Teuchos::RefCountPtr<LOCA::Parameter::SublistParser>& topParams,
	const Teuchos::RefCountPtr<Teuchos::ParameterList>& eigenParams);
      
      //! Destructor
      virtual ~AnasaziStrategy();

      //! Compute eigenvalues/eigenvectors
      /*!
       * The implementation here the sets up and calls the Anasazi 
       * BlockKrylovSchur solver for computing eigenvalues.
       */
      virtual NOX::Abstract::Group::ReturnType
      computeEigenvalues(
	      NOX::Abstract::Group& group,
	      Teuchos::RefCountPtr< std::vector<double> >& evals_r,
	      Teuchos::RefCountPtr< std::vector<double> >& evals_i,
	      Teuchos::RefCountPtr< NOX::Abstract::MultiVector >& evecs_r,
	      Teuchos::RefCountPtr< NOX::Abstract::MultiVector >& evecs_i);

    private:

      //! Private to prohibit copying
      AnasaziStrategy(const AnasaziStrategy&);

      //! Private to prohibit copying
      AnasaziStrategy& operator = (const AnasaziStrategy&);

    protected:

      //! Global data
      Teuchos::RefCountPtr<LOCA::GlobalData> globalData;

      //! Parsed parameter list
      Teuchos::RefCountPtr<LOCA::Parameter::SublistParser> topParams;

      //! Typedef to simplify templating
      typedef NOX::Abstract::MultiVector MV;

      //! Typedef to simplify templating
      typedef LOCA::AnasaziOperator::AbstractStrategy OP;

      //! Eigensolver parameter list
      Teuchos::RefCountPtr<Teuchos::ParameterList> eigenParams;

      //! Linear solver parameter list
      Teuchos::RefCountPtr<Teuchos::ParameterList> solverParams;

      //! The block size
      int blksz;

      //! The maximum length of the Arnoldi factorization
      int length;

      //! The number of requested eigenvalues
      int nev;   
      
      //! Tolerance for the converged eigenvalues
      double tol;
      
      //! This checks convergence every so many steps
      int step;    
      
      //! This is the number of restarts allowed
      int restart;
      
      //! Anasazi Debug level
      int debug;

      //! Which eigenvalues are of interest
      string which;

      //! How many eigenvectors/eigenvalues to we save
      int saveEV;

      //! Pole for Cayley transformation
      double cayleyPole;

      //! Zero for Cayley transformation
      double cayleyZero;

      //! Parameter list to pass into eigensolver
      Teuchos::ParameterList LOCA_PL;

      //! Output manager to handle the I/O from the eigensolver
      Teuchos::RefCountPtr<Anasazi::OutputManager<double> > LOCA_OM;

      //! Sorting manager to handle the sorting of eigenvalues
      Teuchos::RefCountPtr<Anasazi::LOCASort> LOCASort;

    }; // Class AnasaziStrategy

  } // Namespace Eigensolver

} // Namespace LOCA

#endif
