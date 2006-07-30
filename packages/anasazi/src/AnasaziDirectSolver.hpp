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

#ifndef ANASAZI_DIRECT_SOLVER_HPP
#define ANASAZI_DIRECT_SOLVER_HPP

/*!     \file AnasaziDirectSolver.hpp
        \brief Abstract base class providing solver capabilities for projected eigenproblems.
*/

/*!    \class Anasazi::DirectSolver
       \brief Anasazi's templated abstract base class providing 
        solver capabilities for projected eigenproblems.

       This class provides concrete, templated implementations of solvers for projected
       eigenproblems.

       \author Chris Baker, Ulrich Hetmaniuk, Rich Lehoucq, and Heidi Thornquist
*/

#include "AnasaziConfigDefs.hpp"
#include "AnasaziTypes.hpp"

#include "Teuchos_ScalarTraits.hpp"
#include "Teuchos_BLAS.hpp"
#include "Teuchos_LAPACK.hpp"

namespace Anasazi {

  template<class ScalarType>
  class DirectSolver 
  {  
  public:
    
    //! @name Constructor/Destructor
    //@{ 

    //! Basic constructor.  
    DirectSolver() {}
    
    //! Destructor.
    virtual ~DirectSolver() {};
    
    //@}
    
    //! @name Eigensolver Projection Methods
    //@{ 

    //! Routine for computing the generalized eigenpairs of the symmetric pencil <tt>(KK, MM)</tt>
    /*!
      @param size [in] Dimension of the eigenproblem (KK, MM)
      @param KK [in] Symmetric "stiffness" matrix 
      @param MM [in] Symmetric Positive "mass" matrix
      @param EV [in] Dense matrix to store the nev eigenvectors 
      @param theta [out] Array to store the eigenvalues 

      \note The code accesses only the upper triangular part of KK and MM.

      \return Integer indicating the number of computed eigenpairs.
    */
    virtual int directSolver(int size, 
                             const Teuchos::SerialDenseMatrix<int,ScalarType> &KK, 
                             const Teuchos::SerialDenseMatrix<int,ScalarType> *MM,
                             Teuchos::SerialDenseMatrix<int,ScalarType> &EV,
                             std::vector<MagnitudeType> &theta) const = 0;
    //@}
    
  };
  
} // end namespace Anasazi

#endif // ANASAZI_DIRECT_SOLVER_HPP

