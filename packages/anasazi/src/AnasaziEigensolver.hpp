
// @HEADER
// ***********************************************************************
//
//                 Anasazi: Block Eigensolvers Package
//                 Copyright (2004) Sandia Corporation
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
// @HEADER

#ifndef ANASAZI_EIGENSOLVER_HPP
#define ANASAZI_EIGENSOLVER_HPP

/*! \file AnasaziEigensolver.hpp
    \brief Pure virtual base class which describes the basic interface to the iterative eigensolver.
*/

#include "AnasaziEigenproblem.hpp"
#include "Teuchos_RefCountPtr.hpp"

/*! \class Anasazi::Eigensolver
  \brief The Anasazi::Eigensolver is a templated virtual base class that defines the
	basic interface that any eigensolver will support.

	The Anasazi::Eigensolver class is responsible for providing the 
	current solver information to the Anasazi::StatusTest object.
*/

namespace Anasazi {

template<class ScalarType, class MV, class OP>
class Eigensolver {
    
  public:

  //@{ \name Constructor/Destructor.

  //! Default Constructor.
  Eigensolver(void) {};

  //! Destructor.
  virtual ~Eigensolver(void) {};
  //@}
  
  //@{ \name Accessor methods

  //! Get the current iteration count.
  virtual int GetNumIters() const = 0;

  //! Get the current restart count of the iteration method.
  /*! Some eigensolvers can perform restarts (i.e. Arnoldi) to reduce memory
	and orthogonalization costs.  For other eigensolvers that don't
	perform restarts (i.e. LOBPCG), this is not a valid stopping criteria.
  */
  virtual int GetNumRestarts() const = 0;

  //! Get the blocksize to be used by the iterative solver in solving this eigenproblem.
  virtual int GetBlockSize() const = 0;

  /*! \brief Get a constant reference to the current linear problem, 
    	which may include a current solution.
  */
  virtual Eigenproblem<ScalarType,MV,OP>& GetEigenproblem() const = 0;

  //@}

  //@{ \name Solver application methods.
    
  /*! \brief This method uses information given to the eigensolver
    to compute approximate solutions to the specified eigenproblem.
    
    \return Status of the solver on completion:
    <ul>
    <li> Ok - Eigensolver computed requested number of eigenvalues
    <li> Unconverged - Eigensolver reached maximum number of iterations/restarts before computing all requested eigenvalues
    <li> Failed - Numerical failure in eigensolver or bad input parameters
    </ul>    
  */
  virtual ReturnType solve() = 0;
  //@}
  
};

} // end Anasazi namespace

#endif /* ANASAZI_EIGENSOLVER_HPP */
