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

#ifndef LOCA_BORDEREDSOLVER_TPETRAHOUSEHOLDER_H
#define LOCA_BORDEREDSOLVER_TPETRAHOUSEHOLDER_H

#include "LOCA_BorderedSolver_AbstractStrategy.H"  // base class
#include "NOX_TpetraTypedefs.hpp"                  // class data element
#include "LOCA_BorderedSolver_HouseholderQR.H"     // class data element
#include "Teuchos_BLAS.hpp"                        // class data element

// forward declarations
namespace Thyra {
  template<typename T> class DefaultLinearOpSource;
}
namespace LOCA {
  class GlobalData;
  namespace Parameter {
    class SublistParser;
  }
  namespace MultiContinuation {
    class ConstraintInterfaceMVDX;
  }
  namespace Thyra {
    class Group;
  }
}

namespace LOCA {

  namespace BorderedSolver {

    //! Bordered system solver strategy based on Householder transformations
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
     * using Householder tranformations.  The algorithm works as follows:
     * First consider a slightly rearranged version of the extended system
     * of equations:
     * \f[
     *     \begin{bmatrix}
     *          C & B^T    \\
     *          A & J
     *     \end{bmatrix}
     *     \begin{bmatrix}
     *        Y \\
     *        X
     *     \end{bmatrix} =
     *     \begin{bmatrix}
     *        G \\
     *        F
     *     \end{bmatrix}.
     * \f]
     * Let
     * \f[
     *     Q^T
     *     \begin{bmatrix}
     *        C^T \\
     *        B
     *     \end{bmatrix} =
     *     \begin{bmatrix}
     *        R \\
     *        0
     *     \end{bmatrix}
     * \f]
     * be the QR decomposition of the constraints matrix where
     * \f$Q\in\Re^{n+m\times n+m}\f$ and \f$R\in\Re^{m\times m}\f$.
     * Define
     * \f[
     *     \begin{bmatrix}
     *        Z_Y \\
     *        Z_X
     *     \end{bmatrix} = Q^T
     *     \begin{bmatrix}
     *        Y \\
     *        X
     *     \end{bmatrix},
     * \f]
     * then the extended system of equations is equivalent to
     * \f[
     *     \begin{bmatrix}
     *          R^T & 0    \\
     *          [A & J] Q
     *     \end{bmatrix}
     *     \begin{bmatrix}
     *        Z_Y \\
     *        Z_X
     *     \end{bmatrix} =
     *     \begin{bmatrix}
     *        G \\
     *        F
     *     \end{bmatrix}
     * \f]
     * and hence
     * \f[
     *   \begin{split}
     *     Z_Y &= R^{-T} G \\
     *     [A \;\; J] Q
     *     \begin{bmatrix}
     *        0 \\
     *        Z_X
     *     \end{bmatrix} &= F - [A \;\; J] Q
     *     \begin{bmatrix}
     *        Z_Y \\
     *        0
     *     \end{bmatrix}.
     *   \end{split}
     * \f]
     * This last equation equation can be written
     * \f[
     *     P Z_X = \tilde{F}
     * \f]
     * where \f$P\in\Re^{n\times n}\f$ is given by
     * \f[
     *     P Z_X = [A \;\; J] Q
     *     \begin{bmatrix}
     *        0 \\
     *        Z_X
     *     \end{bmatrix}
     * \f]
     * and
     * \f[
     *     \tilde{F} = F - [A \;\; J] Q
     *     \begin{bmatrix}
     *        Z_Y \\
     *        0
     *     \end{bmatrix}.
     * \f]
     * We then recover \f$X\f$ and \f$Y\f$ by
     * \f[
     *     \begin{bmatrix}
     *       Y \\
     *       X
     *     \end{bmatrix} = Q
     *     \begin{bmatrix}
     *        Z_Y \\
     *        Z_X
     *     \end{bmatrix}.
     * \f]
     * It can be further shown that the \f$P\f$ operator above can be
     * written
     * \f[
     *     P = J + U V^T
     * \f]
     * where \f$U = A*Y_1 + J*Y_2\f$, \f$V = Y_2*T^T\f$ and
     * \f$Y = [Y_1 ; Y_2]\f$.
     * The equation \f$P Z_X = \tilde{F}\f$ is solved using an iterative solver
     * using the definition of \f$P Z_X\f$ above, in this case AztecOO.  The
     * system is preconditioned using the preconditioner for \f$J\f$.  The
     * operator \f$Q\f$ is generated using the standard Householder QR
     * algorithm (Algorithm 5.2.1, G. Golub and C. Van Loan, "Matrix
     * Computations," 3rd Edition, Johns Hopkins, Baltimore, 1996) and is
     * stored using the compact WY representation:  \f$Q = I + Y T Y^T\f$
     * (see R. Schreiber and C. Van Loan, "A Storage-Efficient WY Representation
     * for Products of Householder Transformations," SIAM J. Sci. Stat.
     * Comput., Vol. 10, No. 1, pp. 53-57, January 1989).
     *
     * The operator representing \f$P\f$ is encapsulated in the class
     * LOCA::Tpetra::LowRankUpdateRowMatrix if \f$J\f$ is an Tpetra::RowMatrix
     * and LOCA::Tpetra::LowRankUpdateOp otherwise.  If the row matrix
     * version is available \f$P\f$ can be scaled and also used to
     * construct a preconditioner.  If "Include UV In Preconditioner"
     * is true as discussed below, the \f$U\f$ and \f$V\f$ terms will
     * be included when computing this preconditioner, which can help
     * stability when \f$J\f$ is nearly singular.
     *
     * The class is intialized via the \c solverParams parameter list
     * argument to the constructor.  The parameters this class recognizes
     * are:
     * <ul>
     * <li> "Preconditioner Method" -- [string] (default: "Jacobian") -
     *      Method for preconditioning the \f$P\f$ operator.  Choices are:
     *      <ul>
     *      <li> "Jacobian" (default) -- Use the preconditioner for \f$J\f$
     *      <li> "SMW" -- Use the Sherman-Morrison-Woodbury formula for the
     *           inverse of \f$P\f$, replacing the inverse of \f$J\f$ with
     *           the preconditioner for \f$J\f$.
     *      </ul>
     * <li> "Scale Augmented Rows" -- [bool] (default: true) -
     *      Scale augmented rows to unit 2-norm before computing QR
     *      factorization.
     * <li> "Include UV In Preconditioner" -- [bool] (default: false) -
     *      Flag indicating whether to use the \f$U\f$ and \f$V\f$ terms
     *      in the preconditioner for \f$P\f$ when using the "Jacobian"
     *      preconditioner method.
     * <li> "Use P For Preconditioner" -- [bool] (default: false) -
     *      Flag indicating whether to use the representation of \f$P\f$ as
     *      a LOCA::Tpetra::LowRankUpdateRowMatrix for computing the
     *      preconditioner when using the "Jacobian" preconditioner method.
     *      This is valid only for preconditioners that accept an
     *      Tpetra::RowMatrix interface.
     * <li> "Transpose Solver Method" -- [string]
     *      (default: "Transpose Preconditioner") Method for preconditioning
     *      the transpose linear system.  See
     *      LOCA::Tpetra::TransposeLinearSystem::Factory for available choices.
     * </ul>
     */
    class TpetraHouseholder : public LOCA::BorderedSolver::AbstractStrategy {

    public:

      //! Constructor.
      /*!
       * \param global_data [in] Global data object
       * \param topParams [in] Parsed top-level parameter list
       * \param solverParams [in] Bordered solver parameters as described above
       */
      TpetraHouseholder(const Teuchos::RCP<LOCA::GlobalData>& global_data,
                        const Teuchos::RCP<LOCA::Parameter::SublistParser>& topParams,
                        const Teuchos::RCP<Teuchos::ParameterList>& solverParams);

      //! Destructor
      virtual ~TpetraHouseholder();

      //! Set blocks
      /*!
       * The \c blockA or \c blockC pointer may be null if either is zero.
       * Whether block B is zero will be determined by querying \c blockB
       * via ConstraintInterface::isConstraintDerivativesXZero.
       */
      virtual void
      setMatrixBlocks(const Teuchos::RCP<const LOCA::BorderedSolver::AbstractOperator>& op,
                      const Teuchos::RCP<const NOX::Abstract::MultiVector>& blockA,
                      const Teuchos::RCP<const LOCA::MultiContinuation::ConstraintInterface>& blockB,
                      const Teuchos::RCP<const NOX::Abstract::MultiVector::DenseMatrix>& blockC);

      //! Intialize solver for a solve
      /*!
       * This should be called after setMatrixBlocks(), but before
       * applyInverse().
       */
      virtual NOX::Abstract::Group::ReturnType initForSolve();

      //! Intialize solver for a transpose solve
      /*!
       * This should be called after setMatrixBlocks(), but before
       * applyInverseTranspose().
       */
      virtual NOX::Abstract::Group::ReturnType initForTransposeSolve();

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
       * \brief Solves the extended system using the technique described
       * above.
       */
      /*!
       * The \em params argument is the linear solver parameters. If
       * \em isZeroF or \em isZeroG is true, than the corresponding
       * \em F or \em G pointers may be NULL.
       *
       * Note that if either the A or B blocks are zero, the system is
       * solved using a simple block elimination scheme instead of the
       * Householder scheme.
       */
      virtual NOX::Abstract::Group::ReturnType
      applyInverse(Teuchos::ParameterList& params,
                   const NOX::Abstract::MultiVector* F,
                   const NOX::Abstract::MultiVector::DenseMatrix* G,
                   NOX::Abstract::MultiVector& X,
                   NOX::Abstract::MultiVector::DenseMatrix& Y) const;

      /*!
       * \brief Solves the transpose of the extended system as defined above
       */
      /*!
       * The \em params argument is the linear solver parameters.
       */
      virtual NOX::Abstract::Group::ReturnType
      applyInverseTranspose(Teuchos::ParameterList& params,
                            const NOX::Abstract::MultiVector* F,
                            const NOX::Abstract::MultiVector::DenseMatrix* G,
                            NOX::Abstract::MultiVector& X,
                            NOX::Abstract::MultiVector::DenseMatrix& Y) const;

    protected:

      /*!
       * \brief Solves the extended system using the technique described
       * above.
       */
      virtual NOX::Abstract::Group::ReturnType
      solve(Teuchos::ParameterList& params,
            const NOX::Abstract::MultiVector* F,
            const NOX::Abstract::MultiVector::DenseMatrix* G,
            NOX::Abstract::MultiVector& X,
            NOX::Abstract::MultiVector::DenseMatrix& Y) const;

      /*!
       * \brief Solves the transpose of the extended system as defined above
       */
      virtual NOX::Abstract::Group::ReturnType
      solveTranspose(Teuchos::ParameterList& params,
                     const NOX::Abstract::MultiVector* F,
                     const NOX::Abstract::MultiVector::DenseMatrix* G,
                     NOX::Abstract::MultiVector& X,
                     NOX::Abstract::MultiVector::DenseMatrix& Y) const;

      //! Compute \f$U\f$ and \f$V\f$ multivectors in \f$P = J+U V^T\f$.
      NOX::Abstract::Group::ReturnType
      computeUV(const NOX::Abstract::MultiVector::DenseMatrix& Y1,
                const NOX::Abstract::MultiVector& Y2,
                const NOX::Abstract::MultiVector::DenseMatrix& T,
                const NOX::Abstract::MultiVector& A,
                NOX::Abstract::MultiVector& U,
                NOX::Abstract::MultiVector& V,
                bool use_jac_transpose);

    public:
      /*!
       * \brief Overwrites the Jacobian \f$J\f$ with \f$J + U V^T\f$
       * for computing the preconditioner of \f$P\f$.
       *
       * NOTE: This should be a protected method, but cuda lambda forces this to be public!
       */
      void updateCrsMatrixForPreconditioner(const NOX::Abstract::MultiVector& U,
                                            const NOX::Abstract::MultiVector& V,
                                            NOX::TCrsMatrix& mat) const;

    protected:
      Teuchos::RCP<NOX::Abstract::MultiVector>
      createBlockMV(const NOX::Abstract::MultiVector& v) const;

      void
      setBlockMV(const NOX::Abstract::MultiVector& bv,
                 NOX::Abstract::MultiVector& v) const;

    private:

      //! Private to prohibit copying
      TpetraHouseholder(const TpetraHouseholder&);

      //! Private to prohibit copying
      TpetraHouseholder& operator = (const TpetraHouseholder&);

    protected:

      //! Enumerated type indicating preconditioner method
      enum PRECONDITIONER_METHOD {
        JACOBIAN,
        SMW
      };

      //! Global data object
      Teuchos::RCP<LOCA::GlobalData> globalData;

      //! Solver parameters
      Teuchos::RCP<Teuchos::ParameterList> solverParams;

      //! Pointer to group storing J
      Teuchos::RCP<LOCA::Thyra::Group> grp;

      // Operator
      Teuchos::RCP<const LOCA::BorderedSolver::AbstractOperator> op;

      //! Pointer to A block
      Teuchos::RCP<const NOX::Abstract::MultiVector> A;

      //! Pointer to B block
      Teuchos::RCP<const NOX::Abstract::MultiVector> B;

      //! Pointer to C block
      Teuchos::RCP<const NOX::Abstract::MultiVector::DenseMatrix> C;

      //! Pointer to constraint interface
      Teuchos::RCP<const LOCA::MultiContinuation::ConstraintInterfaceMVDX> constraints;

      //! QR Factorization object
      LOCA::BorderedSolver::HouseholderQR qrFact;

      //! Solution component of Householder multivec
      Teuchos::RCP<NOX::Abstract::MultiVector> house_x;

      //! Parameter component of Householder multivec
      NOX::Abstract::MultiVector::DenseMatrix house_p;

      //! T matrix in compact WY representation
      NOX::Abstract::MultiVector::DenseMatrix T;

      //! R matrix in QR factorization
      NOX::Abstract::MultiVector::DenseMatrix R;

      //! U matrix in low-rank update form P = J + U*V^T
      Teuchos::RCP<NOX::Abstract::MultiVector> U;

      //! V matrix in low-rank update form P = J + U*V^T
      Teuchos::RCP<NOX::Abstract::MultiVector> V;

      //! Solution component of Householder multivec for transposed system
      Teuchos::RCP<NOX::Abstract::MultiVector> house_x_trans;

      //! Parameter component of Householder multivec for transposed system
      NOX::Abstract::MultiVector::DenseMatrix house_p_trans;

      //! T matrix in compact WY representation for transposed system
      NOX::Abstract::MultiVector::DenseMatrix T_trans;

      //! R matrix in QR factorization for transposed system
      NOX::Abstract::MultiVector::DenseMatrix R_trans;

      //! U matrix in low-rank update form P = J + U*V^T for transposed system
      Teuchos::RCP<NOX::Abstract::MultiVector> U_trans;

      //! V matrix in low-rank update form P = J + U*V^T for transposed system
      Teuchos::RCP<NOX::Abstract::MultiVector> V_trans;

      //! Pointer to A block as an Tpetra multivector
      Teuchos::RCP<const NOX::Abstract::MultiVector> Ablock;

      //! Pointer to B block as an Tpetra multivector
      Teuchos::RCP<const NOX::Abstract::MultiVector> Bblock;

      //! Pointer to scaled A block
      Teuchos::RCP<NOX::Abstract::MultiVector> Ascaled;

      //! Pointer to scaled B block
      Teuchos::RCP<NOX::Abstract::MultiVector> Bscaled;

      //! Pointer to scaled C block
      Teuchos::RCP<NOX::Abstract::MultiVector::DenseMatrix> Cscaled;

      //! Pointer to Tpetra J operator
      Teuchos::RCP<NOX::TOperator> tpetraOp;

      //! Pointer to Tpetra Preconditioner operator
      mutable Teuchos::RCP<NOX::TCrsMatrix> tpetraPrecMatrix;

      //! Thyra wrapped preconditioner matrix (tpetraPrecMatrix) for when includeUV is true and use_P_for_Prec is false
      mutable Teuchos::RCP<::Thyra::DefaultLinearOpSource<double>> prec_losb;

      //! Number of constraint equations
      int numConstraints;

      //! flag indicating whether A block is zero
      bool isZeroA;

      //! flag indicating whether B block is zero
      bool isZeroB;

      //! flag indicating whether C block is zero
      bool isZeroC;

      /*!
       * \brief Flag indicating whether constraint factorization for solve
       * has been computed
       */
      bool isValidForSolve;

      /*!
       * \brief Flag indicating whether constraint factorization for transpoe
       * solve has been computed
       */
      bool isValidForTransposeSolve;

      //! BLAS Wrappers
      Teuchos::BLAS<int,double> dblas;

      //! Whether we should scale augmented rows to have unit 2-norm
      bool scale_rows;

      //! Scale values for each row
      std::vector<double> scale_vals;

      //! Preconditioner method
      PRECONDITIONER_METHOD precMethod;

      //! Flag indicating whether to include U*V^T terms in preconditioner
      bool includeUV;

      //! Flag indicating whether to use P = J + U*V^T in preconditioner
      bool use_P_For_Prec;

      //! Flag indicating whether we are doing a complex solve
      bool isComplex;

      //! Frequency for complex systems
      double omega;
    };
  } // namespace BorderedSolver
} // namespace LOCA

#endif
