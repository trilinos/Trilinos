// @HEADER
// ***********************************************************************
//
//                 Anasazi: Block Eigensolvers Package
//                 Copyright 2004 Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
                             std::vector<typename Teuchos::ScalarTraits<ScalarType>::magnitudeType> &theta) const = 0;
    //@}
    
  };
  
} // end namespace Anasazi

#endif // ANASAZI_DIRECT_SOLVER_HPP

