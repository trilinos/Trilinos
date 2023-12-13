// @HEADER
// ***********************************************************************
//
//                 Belos: Block Linear Solvers Package
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
// Questions? Contact Jennifer A. Loe (jloe@sandia.gov)
//                or Heidi K. Thornquist (hkthorn@sandia.gov)
//
// ***********************************************************************
// @HEADER
//
#ifndef BELOS_DENSE_SOLVER_HPP
#define BELOS_DENSE_SOLVER_HPP

/*! \file BelosDenseSolver.hpp
  \brief Virtual base class which defines solvers for dense matrix type
*/

#include "Teuchos_RCP.hpp"
#include "Teuchos_BLAS_types.hpp"

namespace Belos {

  template<class ScalarType, class DM>
  class DenseSolver 
  {
  public:
 
    //! @name Constructor/Destructor Methods
    //@{
    //! Default constructor; matrix should be set using setMatrix(), LHS and RHS set with setVectors().
    DenseSolver()
    : spd_(false),
      equilibrate_(false),
      TRANS_(Teuchos::NO_TRANS),
      newMatrix_(false)
    {}

    //! SerialDenseSolver destructor.
    virtual ~DenseSolver() {}
    //@}

    //! @name Set Methods
    //@{

    //! Sets the pointers for coefficient matrix
    virtual int setMatrix(const Teuchos::RCP<DM>& A) { A_ = A; newMatrix_ = true; return 0; }

    //! Sets the pointers for left and right hand side vector(s).
    /*! Row dimension of X must match column dimension of matrix A, row dimension of B
      must match row dimension of A.  X and B must have the same dimensions.
    */
    virtual int setVectors(const Teuchos::RCP<DM>& X, const Teuchos::RCP<DM>& B) { X_ = X; B_ = B; return 0; }
    //@}

    //! @name Strategy Modifying Methods
    //@{

    //! Set if dense matrix is symmetric positive definite
    //! \note This method must be called before the factorization is performed, otherwise it will have no effect.
    void setSPD(bool flag) {spd_ = flag; return;}

    //! Causes equilibration to be called just before the matrix factorization as part of the call to \c factor.
    /*! \note This method must be called before the factorization is performed, otherwise it will have no effect.
    */
    void factorWithEquilibration(bool flag) {equilibrate_ = flag; return;}

    //! All subsequent function calls will work with the transpose-type set by this method (\c Teuchos::NO_TRANS, \c Teuchos::TRANS, and \c Teuchos::CONJ_TRANS).
    /*! \note This interface will set correct transpose flag for matrix, including complex-valued linear systems.
    */
    void solveWithTransposeFlag(Teuchos::ETransp trans) {TRANS_ = trans; return;}
    //@}

    //! @name Factor/Solve Methods
    //@{

    //! Computes the in-place LU factorization of the matrix.
    /*!
      \return Integer error code, set to 0 if successful.
    */
    virtual int factor() = 0;

    //! Computes the solution X to AX = B for the \e this matrix and the B provided.
    /*!
      \return Integer error code, set to 0 if successful.
    */
    virtual int solve() = 0;
    //@}

  protected:
  
    bool spd_, equilibrate_;
    Teuchos::ETransp TRANS_;

    bool newMatrix_; 
    Teuchos::RCP<DM> A_, X_, B_;
  };

} // namespace Belos

#endif // end file BELOS_DENSE_SOLVER_HPP
