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

#ifdef HAVE_LOCA_ANASAZI
#include "Teuchos_ParameterList.hpp"           // class element
#include "AnasaziOutputManager.hpp"            // class element
#include "Anasazi_LOCA_Sort.H"                 // class element
#include "AnasaziMultiVec.hpp"                 // for typedef
#include "AnasaziOperator.hpp"                 // for typedef
#endif

namespace LOCA {

  namespace Eigensolver {

    //! Anasazi eigensolver strategy
    /*!
     * This class implements an eigensolver strategy using the generic
     * Trilinos eigensolver package Anasazi.  Since Anasazi is a generic
     * solver, this strategy will work with any group implementation.
     */
    class AnasaziStrategy : public LOCA::Eigensolver::AbstractStrategy {

    public:

      //! Constructor
      AnasaziStrategy(
	       const Teuchos::RefCountPtr<NOX::Parameter::List>& eigenParams,
	       const Teuchos::RefCountPtr<NOX::Parameter::List>& solverParams);
      
      //! Destructor
      virtual ~AnasaziStrategy();

      //! Compute eigenvalues/eigenvectors
      /*!
       * The implementation here the sets up and calls the Anasazi 
       * BlockKrylovSchur solver for computing eigenvalues.
       */
      virtual NOX::Abstract::Group::ReturnType
      computeEigenvalues(
		    const Teuchos::RefCountPtr<NOX::Abstract::Group>& group);

    protected:

#ifdef HAVE_LOCA_ANASAZI
      //! Typedef to simplify templating
      typedef Anasazi::MultiVec<double> MV;

      //! Typedef to simplify templating
      typedef Anasazi::Operator<double> OP;

      //! Eigensolver parameter list
      Teuchos::RefCountPtr<NOX::Parameter::List> eigenParams;

      //! Linear solver parameter list
      Teuchos::RefCountPtr<NOX::Parameter::List> solverParams;

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
      Teuchos::RefCountPtr<Anasazi::LOCASort<double, MV, OP> > LOCASort;
#endif

    }; // Class AnasaziStrategy

  } // Namespace Eigensolver

} // Namespace LOCA

#endif
