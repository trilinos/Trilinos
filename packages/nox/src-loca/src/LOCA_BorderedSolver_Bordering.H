// $Id$ 
// $Source$ 

//@HEADER
// ************************************************************************
// 
//            LOCA: Library of Continuation Algorithms Package
//                 Copyright (2005) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
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
// Questions? Contact Roger Pawlowski (rppawlo@sandia.gov) or 
// Eric Phipps (etphipp@sandia.gov), Sandia National Laboratories.
// ************************************************************************
//  CVS Information
//  $Source$
//  $Author$
//  $Date$
//  $Revision$
// ************************************************************************
//@HEADER

#ifndef LOCA_BORDEREDSOLVER_BORDERING_H
#define LOCA_BORDEREDSOLVER_BORDERING_H

#include "LOCA_BorderedSolver_AbstractStrategy.H"  // base class
#include <vector>

// forward declarations
namespace LOCA {
  class GlobalData;
  namespace Parameter {
    class SublistParser;
  }
}

namespace LOCA {

  namespace BorderedSolver {
    
    //! Bordered system solver strategy based on bordering
    /*! 
     * This class solves the extended system of equations
     * \f[
     *     \begin{bmatrix}
     *          J & A    \\
     *        B^T & C
     *     \end{bmatrix}
     *     \begin{bmatrix}
     *        X \\
     *        Y 
     *     \end{bmatrix} = 
     *     \begin{bmatrix}
     *        F \\
     *        G
     *     \end{bmatrix}
     * \f]
     * via bordering (block elimination):
     * \f[
     *     \begin{aligned}
     *        X_1 &= J^{-1} F \\
     *        X_2 &= J^{-1} A \\
     *        Y   &= (C-B^T X_2)^{-1}(G-B^T X_1) \\
     *        X   &= X_1 - X_2 Y
     *     \end{aligned}
     * \f]
     * It takes advantage of any of the matrix blocks being zero and
     * concatenates \f$F\f$ and \f$A\f$ into a contiguous multivector
     * to compute \f$X_1\f$ and \f$X_2\f$ in one block solve.
     *
     * To solve the transpose of the system, a similar bordering algorithm
     * is implemented.  Note however that for the transpose, the constraint
     * object representing \f$B\f$ must implement the 
     * LOCA::MultiContinuation::ConstraintInterfaceMVDX since \f$B\f$ 
     * appears on the right-hand-side of a linear system.
    */
    class Bordering : public LOCA::BorderedSolver::AbstractStrategy {

    public:

      //! Constructor. 
      /*!
       * \param global_data [in] Global data object
       * \param topParams [in] Parsed top-level parameter list 
       * \param solverParams [in] Bordered solver parameters.  Currently none
       * are referenced.
       */
      Bordering(
	 const Teuchos::RCP<LOCA::GlobalData>& global_data,
	 const Teuchos::RCP<LOCA::Parameter::SublistParser>& topParams,
	 const Teuchos::RCP<Teuchos::ParameterList>& solverParams);

      //! Destructor
      virtual ~Bordering();

      //! Set blocks
      /*!
       * The \c blockA or \c blockC pointer may be null if either is zero.
       * Whether block B is zero will be determined by querying \c blockB
       * via ConstraintInterface::isConstraintDerivativesXZero.
       */
      virtual void setMatrixBlocks(
	  const Teuchos::RCP<const LOCA::BorderedSolver::AbstractOperator>& op,
	  const Teuchos::RCP<const NOX::Abstract::MultiVector>& blockA,
	  const Teuchos::RCP<const LOCA::MultiContinuation::ConstraintInterface>& blockB,
	  const Teuchos::RCP<const NOX::Abstract::MultiVector::DenseMatrix>& blockC);

       //! Intialize solver for a solve
      /*!
       * This should be called after setMatrixBlocks(), but before 
       * applyInverse().
       */
      virtual NOX::Abstract::Group::ReturnType 
      initForSolve();

      //! Intialize solver for a transpose solve
      /*!
       * This should be called after setMatrixBlocks(), but before 
       * applyInverseTranspose().
       */
      virtual NOX::Abstract::Group::ReturnType 
      initForTransposeSolve();
  
      /*! 
       * \brief Computed extended matrix-multivector product
       */
      /*!
       * Computes
       * \f[
       *     \begin{bmatrix}
       *        U \\
       *        V
       *     \end{bmatrix} = 
       *     \begin{bmatrix}
       *          J & A    \\
       *        B^T & C
       *     \end{bmatrix}
       *     \begin{bmatrix}
       *        X \\
       *        Y 
       *     \end{bmatrix} = 
       *     \begin{bmatrix}
       *          J*X + A*Y \\
       *        B^T*X + C*Y
       *     \end{bmatrix}.
       * \f]
       */
      virtual NOX::Abstract::Group::ReturnType 
      apply(const NOX::Abstract::MultiVector& X,
	    const NOX::Abstract::MultiVector::DenseMatrix& Y,
	    NOX::Abstract::MultiVector& U,
	    NOX::Abstract::MultiVector::DenseMatrix& V) const;

      /*! 
       * \brief Computed extended matrix transpose-multivector product
       */
      /*!
       * Computes
       * \f[
       *     \begin{bmatrix}
       *        U \\
       *        V
       *     \end{bmatrix} = 
       *     \begin{bmatrix}
       *        J^T & B    \\
       *        A^T & C
       *     \end{bmatrix}
       *     \begin{bmatrix}
       *        X \\
       *        Y 
       *     \end{bmatrix} = 
       *     \begin{bmatrix}
       *        J^T*X + B*Y \\
       *        A^T*X + C^T*Y
       *     \end{bmatrix}.
       * \f]
       */
      virtual NOX::Abstract::Group::ReturnType 
      applyTranspose(const NOX::Abstract::MultiVector& X,
		     const NOX::Abstract::MultiVector::DenseMatrix& Y,
		     NOX::Abstract::MultiVector& U,
		     NOX::Abstract::MultiVector::DenseMatrix& V) const;

      /*! 
       * \brief Solves the extended system as defined above using bordering
       */
      /*!
       * The \em params argument is the linear solver parameters. If
       * \em isZeroF or \em isZeroG is true, than the corresponding
       * \em F or \em G pointers may be NULL.
       */
      virtual NOX::Abstract::Group::ReturnType 
      applyInverse(Teuchos::ParameterList& params,
		   const NOX::Abstract::MultiVector* F,
		   const NOX::Abstract::MultiVector::DenseMatrix* G,
		   NOX::Abstract::MultiVector& X,
		   NOX::Abstract::MultiVector::DenseMatrix& Y) const;

      /*! 
       * \brief Solves the transpose of the extended system as defined above 
       * using bordering
       */
      /*!
       * The \em params argument is the linear solver parameters. If
       * \em isZeroF or \em isZeroG is true, than the corresponding
       * \em F or \em G pointers may be NULL.  Note that for the transpose
       * solve B must be of type 
       * LOCA::MultiContinuation::ConstraintInterfaceMVDX.
       */
      virtual NOX::Abstract::Group::ReturnType 
      applyInverseTranspose(Teuchos::ParameterList& params,
			    const NOX::Abstract::MultiVector* F,
			    const NOX::Abstract::MultiVector::DenseMatrix* G,
			    NOX::Abstract::MultiVector& X,
			    NOX::Abstract::MultiVector::DenseMatrix& Y) const;

    protected:

      /*! 
       * \brief Solves the extended system when F = 0
       */
      NOX::Abstract::Group::ReturnType  
      solveFZero(Teuchos::ParameterList& params,
		 const NOX::Abstract::MultiVector* AA,
		 const LOCA::MultiContinuation::ConstraintInterface* BB,
		 const NOX::Abstract::MultiVector::DenseMatrix* CC,
		 const NOX::Abstract::MultiVector::DenseMatrix* G,
		 NOX::Abstract::MultiVector& X,
		 NOX::Abstract::MultiVector::DenseMatrix& Y) const;

      /*! 
       * \brief Solves the extended system when F and A are contiguous
       */
      NOX::Abstract::Group::ReturnType  
      solveContiguous(Teuchos::ParameterList& params,
		      const NOX::Abstract::MultiVector* AA,
		      const LOCA::MultiContinuation::ConstraintInterface* BB,
		      const NOX::Abstract::MultiVector::DenseMatrix* CC,
		      std::vector<int>& indexF,
		      std::vector<int>& indexA,
		      const NOX::Abstract::MultiVector* F,
		      const NOX::Abstract::MultiVector::DenseMatrix* G,
		      NOX::Abstract::MultiVector& X,
		      NOX::Abstract::MultiVector::DenseMatrix& Y) const;

      /*! 
       * \brief Solves the transpose of the extended system when F = 0
       */
      NOX::Abstract::Group::ReturnType  
      solveFZeroTrans(Teuchos::ParameterList& params,
		      const NOX::Abstract::MultiVector* AA,
		      const NOX::Abstract::MultiVector* BB,
		      const NOX::Abstract::MultiVector::DenseMatrix* CC,
		      const NOX::Abstract::MultiVector::DenseMatrix* G,
		      NOX::Abstract::MultiVector& X,
		      NOX::Abstract::MultiVector::DenseMatrix& Y) const;

      /*! 
       * \brief Solves the transpose of the extended system when F and B are 
       * contiguous
       */
      NOX::Abstract::Group::ReturnType  
      solveContiguousTrans(
		      Teuchos::ParameterList& params,
		      const NOX::Abstract::MultiVector* AA,
		      const NOX::Abstract::MultiVector* BB,
		      const NOX::Abstract::MultiVector::DenseMatrix* CC,
		      std::vector<int>& indexF,
		      std::vector<int>& indexA,
		      const NOX::Abstract::MultiVector* F,
		      const NOX::Abstract::MultiVector::DenseMatrix* G,
		      NOX::Abstract::MultiVector& X,
		      NOX::Abstract::MultiVector::DenseMatrix& Y) const;

    private:

      //! Private to prohibit copying
      Bordering(const Bordering&);

      //! Private to prohibit copying
      Bordering& operator = (const Bordering&);

    protected:

      //! Global data object
      Teuchos::RCP<LOCA::GlobalData> globalData;

      //! Solver parameters
      Teuchos::RCP<Teuchos::ParameterList> solverParams;

      //! Pointer to operator
      Teuchos::RCP<const LOCA::BorderedSolver::AbstractOperator> op;

      //! Pointer to A block
      Teuchos::RCP<const NOX::Abstract::MultiVector> A;

      //! Pointer to B block
      Teuchos::RCP<const LOCA::MultiContinuation::ConstraintInterface> B;

      //! Pointer to C block
      Teuchos::RCP<const NOX::Abstract::MultiVector::DenseMatrix> C;

      //! flag indicating whether A block is zero
      bool isZeroA;

      //! flag indicating whether B block is zero
      bool isZeroB;

      //! flag indicating whether C block is zero
      bool isZeroC;

      //! flag indicating whether F block is zero
      mutable bool isZeroF;

      //! flag indicating whether G block is zero
      mutable bool isZeroG;

    };
  } // namespace BorderedSolver
} // namespace LOCA

#endif
