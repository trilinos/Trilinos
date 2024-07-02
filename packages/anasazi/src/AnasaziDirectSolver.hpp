// @HEADER
// *****************************************************************************
//                 Anasazi: Block Eigensolvers Package
//
// Copyright 2004 NTESS and the Anasazi contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
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
                             std::vector<typename Teuchos::ScalarTraits<ScalarType>::magnitudeType> &theta) const = 0;
    //@}
    
  };
  
} // end namespace Anasazi

#endif // ANASAZI_DIRECT_SOLVER_HPP

