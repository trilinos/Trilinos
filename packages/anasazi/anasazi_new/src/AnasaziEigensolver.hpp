
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
#include "AnasaziMultiVec.hpp"

/*! \class Anasazi::Eigensolver
  \brief The Anasazi::Eigensolver is a templated virtual base class that defines the
	basic interface that any eigensolver will support.

	The Anasazi::Eigensolver class is responsible for providing the 
	current solver information to the Anasazi::StatusTest object.
*/

namespace Anasazi {

template<class TYPE>
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

  //! Get the solvers native residuals for the current eigenpairs. 
  /*! This is not be the same as the true residuals for most solvers. Sometimes the native
    residuals are not in multivector form, so the norm type is solver dependent.  

    \note
    <ol>
      <li> If the native residual is in multivector form then a non-null pointer will be
      returned, else the normvec will be populated with the current residual norms. 
      <li> If the native residual is returned in multivector form, the memory is managed
      by the calling routine.
    </ol>
  */
  virtual MultiVec<TYPE>* GetNativeResiduals( TYPE* normvec ) const = 0;

  /*! \brief Get a constant reference to the current linear problem, 
    	which may include a current solution.
  */
  virtual Eigenproblem<TYPE>& GetEigenproblem() const = 0;

  //@}

};

} // end Anasazi namespace

#endif /* ANASAZI_EIGENSOLVER_HPP */
