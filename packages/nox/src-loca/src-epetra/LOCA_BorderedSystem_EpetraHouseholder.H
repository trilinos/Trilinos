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

#ifndef LOCA_BORDEREDSYSTEM_EPETRAHOUSEHOLDER_H
#define LOCA_BORDEREDSYSTEM_EPETRAHOUSEHOLDER_H

#include "LOCA_BorderedSystem_AbstractStrategy.H"  // base class

#include "Teuchos_BLAS.hpp" // class data element

// forward declarations
namespace LOCA {
  class GlobalData;
  namespace Parameter {
    class SublistParser;
  }
  namespace MultiContinuation {
    class ConstraintInterfaceMVDX;
  }
  namespace EpetraNew {
    class Group;
  }
}

namespace LOCA {

  namespace BorderedSystem {
    
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
     *     [A J] Q
     *     \begin{bmatrix}
     *        0 \\
     *        Z_X
     *     \end{bmatrix} &= F - [A J] Q 
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
     *     P Z_X = [A J] Q
     *     \begin{bmatrix}
     *        0 \\
     *        Z_X
     *     \end{bmatrix}
     * \f]
     * and 
     * \f[
     *     \tilde{F} = F - [A J] Q 
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
     * The equation \f$P Z_X = \tilde{F}\f$ is solved using an iterative solver
     * using the definition of \f$P Z_X\f$ above, in this case AztecOO.  The
     * system is preconditioned using the preconditioner for \f$J\f$.  The
     * operator \f$Q\f$ is generated using the standard Householder QR
     * algorithm (Algorithm 5.2.1, G. Golub and C. Van Loan, "Matrix
     * Computations," 3rd Edition, Johns Hopkins, Baltimore, 1996) and is
     * stored using the compact WY representation:  \f$Q = I + Y^T T Y\f$
     * (see R. Schreiver and C. Van Loan, "A Storage-Efficient WY Represntation
     * for Products of Householder Transformations," SIAM J. Sci. Stat. 
     * Comput., Vol. 10, No. 1, pp. 53-57, January 1989).  The operator
     * representing \f$P\f$ is encapsulated in the class 
     * LOCA::Epetra::CompactWYOp.
     */
    class EpetraHouseholder : public LOCA::BorderedSystem::AbstractStrategy {

    public:

      //! Constructor. 
      /*!
       * \param global_data [in] Global data object
       * \param topParams [in] Parsed top-level parameter list 
       * \param solverParams [in] Bordered solver parameters.  Currently none
       * are referenced.
       */
      EpetraHouseholder(
	 const Teuchos::RefCountPtr<LOCA::GlobalData>& global_data,
	 const Teuchos::RefCountPtr<LOCA::Parameter::SublistParser>& topParams,
	 const Teuchos::RefCountPtr<NOX::Parameter::List>& solverParams);

      //! Destructor
      virtual ~EpetraHouseholder();

      /*! 
       * \brief Set flag indicating whether F and A are stored in a 
       * continguous multivector.
       */
      virtual void setIsContiguous(bool flag);

      //! Set blocks
      /*!
       * The \c blockA or \c blockC pointer may be null if either is zero.
       * Whether block B is zero will be determined by querying \c blockB
       * via ConstraintInterface::isConstraintDerivativesXZero.
       */
      virtual void setMatrixBlocks(
	  const Teuchos::RefCountPtr<const NOX::Abstract::Group>& group,
	  const Teuchos::RefCountPtr<const NOX::Abstract::MultiVector>& blockA,
	  const Teuchos::RefCountPtr<const LOCA::MultiContinuation::ConstraintInterface>& blockB,
	  const Teuchos::RefCountPtr<const NOX::Abstract::MultiVector::DenseMatrix>& blockC);
  
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
      applyInverse(NOX::Parameter::List& params,
		   const NOX::Abstract::MultiVector* F,
		   const NOX::Abstract::MultiVector::DenseMatrix* G,
		   NOX::Abstract::MultiVector& X,
		   NOX::Abstract::MultiVector::DenseMatrix& Y) const;

    protected:

      /*!
       * \brief Factors the constraint derivatives using the QR-factorization
       * and generates the Y, T, and R matrices as described above.
       */
      void factorConstraints();

      /*!
       * \brief Computes the Householder vector \c V_x, \c V_p for column 
       * \c col of \c A_x and \c A_p, starting at row \c col of A_p.
       */
      /*!
       * The algorithm implemented here is essentially algorithm 5.1.1
       * of Golub and Van Loan.
       */
      void computeHouseholderVector(
			  int col,
			  const NOX::Abstract::MultiVector& A_x,
			  const NOX::Abstract::MultiVector::DenseMatrix& A_p,
			  NOX::Abstract::MultiVector& V_x,
			  NOX::Abstract::MultiVector::DenseMatrix& V_p,
			  double& beta);

      /*!
       * \brief Applies the Householder vector \c V_x, \c V_p to the
       * matrix sub-block represented by \c A_x and A_p.
       */
      void applyHouseholderVector(
			   const NOX::Abstract::MultiVector& V_x,
			   const NOX::Abstract::MultiVector::DenseMatrix& V_p,
			   double beta,
			   NOX::Abstract::MultiVector& A_x,
			   NOX::Abstract::MultiVector::DenseMatrix& A_p);

      /*!
       * \brief Applies the operator Q as described above overwriting
       * \c x and \c y.  If either of \c x or \c y are zero on input,
       * set the corresponding \c isZeroX or \c isZeroY flags.  Set\
       * \c useTranspose to \c true to instead apply the transpose of Q.
       */
      void applyCompactWY(NOX::Abstract::MultiVector& x,
			  NOX::Abstract::MultiVector::DenseMatrix& y,
			  bool isZeroX, bool isZeroY, bool useTranspose) const;

      /*!
       * \brief Another version of applyCompactWY() that does not overwrite
       * its inputs.  If either input is zero, set the corresponding pointer
       * to NULL.
       */
      void applyCompactWY(
		     const NOX::Abstract::MultiVector* input_x,
		     const NOX::Abstract::MultiVector::DenseMatrix* input_y,
		     NOX::Abstract::MultiVector& result_x,
		     NOX::Abstract::MultiVector::DenseMatrix& result_y,
		     bool useTranspose) const;

      /*!
       * \brief Solves the extended system of equations when the A block
       * is zero using a simple block elimination.
       */
      NOX::Abstract::Group::ReturnType 
      solveAZero(NOX::Parameter::List& params,
		 const LOCA::MultiContinuation::ConstraintInterface* BB,
		 const NOX::Abstract::MultiVector::DenseMatrix* CC,
		 const NOX::Abstract::MultiVector* F,
		 const NOX::Abstract::MultiVector::DenseMatrix* G,
		 NOX::Abstract::MultiVector& X,
		 NOX::Abstract::MultiVector::DenseMatrix& Y) const;

      /*!
       * \brief Solves the extended system of equations when the B block
       * is zero using a simple block elimination.
       */
      NOX::Abstract::Group::ReturnType 
      solveBZero(NOX::Parameter::List& params,
		 const NOX::Abstract::MultiVector* AA,
		 const NOX::Abstract::MultiVector::DenseMatrix* CC,
		 const NOX::Abstract::MultiVector* F,
		 const NOX::Abstract::MultiVector::DenseMatrix* G,
		 NOX::Abstract::MultiVector& X,
		 NOX::Abstract::MultiVector::DenseMatrix& Y) const;

      

    private:

      //! Private to prohibit copying
      EpetraHouseholder(const EpetraHouseholder&);

      //! Private to prohibit copying
      EpetraHouseholder& operator = (const EpetraHouseholder&);

    protected:

      //! Global data object
      Teuchos::RefCountPtr<LOCA::GlobalData> globalData;

      //! Solver parameters
      Teuchos::RefCountPtr<NOX::Parameter::List> solverParams;

      //! Pointer to group storing J
      Teuchos::RefCountPtr<const LOCA::EpetraNew::Group> grp;

      //! Pointer to A block
      Teuchos::RefCountPtr<const NOX::Abstract::MultiVector> A;

      //! Pointer to B block
      Teuchos::RefCountPtr<const NOX::Abstract::MultiVector> B;

      //! Pointer to C block
      Teuchos::RefCountPtr<const NOX::Abstract::MultiVector::DenseMatrix> C;

      //! Pointer to constraint interface
      Teuchos::RefCountPtr<const LOCA::MultiContinuation::ConstraintInterfaceMVDX> constraints;

      //! Solution component of Householder multivec
      Teuchos::RefCountPtr<NOX::Abstract::MultiVector> house_x;

      //! Parameter component of Householder multivec
      NOX::Abstract::MultiVector::DenseMatrix house_p;

      //! "T" matrix in compact WY representation
      NOX::Abstract::MultiVector::DenseMatrix T;

      //! R matrix in QR factorization
      NOX::Abstract::MultiVector::DenseMatrix R;

      //! Vector store 1 single Householder vector
      Teuchos::RefCountPtr<NOX::Abstract::MultiVector> v_x;

      //! Number of constraint equations
      int numConstraints;

      //! flag indicating whether A block is zero
      bool isZeroA;

      //! flag indicating whether B block is zero
      bool isZeroB;

      //! flag indicating whether C block is zero
      bool isZeroC;

      //! flag indicating whether F and A are contiguous
      bool isContiguous;

      //! BLAS Wrappers
      Teuchos::BLAS<int,double> dblas;

    };
  } // namespace BorderedSystem
} // namespace LOCA

#endif
