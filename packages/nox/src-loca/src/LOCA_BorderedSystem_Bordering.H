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

#ifndef LOCA_BORDEREDSYSTEM_BORDERING_H
#define LOCA_BORDEREDSYSTEM_BORDERING_H

#include "LOCA_BorderedSystem_AbstractStrategy.H"  // base class

// forward declarations
namespace LOCA {
  class GlobalData;
  namespace Parameter {
    class SublistParser;
  }
}

namespace LOCA {

  namespace BorderedSystem {
    
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
     *        Y   &= (C-B^T X_2)^{-1}(G-B^T X_1)
     *        X   &= X_1 - X_2 Y
     *     \end{aligned}
     * \f]
     * It takes advantage of any of the matrix blocks being zero and
     * concatenates \f$F\f$ and \f$A\f$ into a contiguous multivector
     * (if they are not already) to compute \f$X_1\f$ and \f$X_2\f$ in
     * one block solve.
    */
    class Bordering : public LOCA::BorderedSystem::AbstractStrategy {

    public:

      //! Constructor. 
      /*!
       * \param global_data [in] Global data object
       * \param topParams [in] Parsed top-level parameter list 
       * \param predParams [in] Bordered solver parameters.  Currently none
       * are referenced.
       */
      Bordering(
	 const Teuchos::RefCountPtr<LOCA::GlobalData>& global_data,
	 const Teuchos::RefCountPtr<LOCA::Parameter::SublistParser>& topParams,
	 const Teuchos::RefCountPtr<NOX::Parameter::List>& solverParams);

      //! Destructor
      virtual ~Bordering();

      //! Set flags indicating if any blocks are zero
      virtual void setIsZero(bool flagA, bool flagB, bool flagC, bool flagF,
			     bool flagG);

      /*! 
       * \brief Set flag indicating whether F and A are stored in a 
       * continguous multivector.
       */
      virtual void setIsContiguous(bool flag);

      //! Set blocks
      /*!
       * The \c blockA pointer may be null if this block is zero
       */
      virtual void setMatrixBlocks(
	  const Teuchos::RefCountPtr<const NOX::Abstract::Group>& group,
	  const Teuchos::RefCountPtr<const NOX::Abstract::MultiVector>& blockA,
	  const Teuchos::RefCountPtr<const LOCA::MultiContinuation::ConstraintInterface>& blockBC);
  
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
      applyInverse(NOX::Parameter::List& params,
		   const NOX::Abstract::MultiVector* F,
		   const NOX::Abstract::MultiVector::DenseMatrix* G,
		   NOX::Abstract::MultiVector& X,
		   NOX::Abstract::MultiVector::DenseMatrix& Y) const;

    protected:

      /*! 
       * \brief Solves the extended system when B = 0 and F and A are
       * not contiguous
       */
      NOX::Abstract::Group::ReturnType  
      solveBZeroNoncontiguous(
			   NOX::Parameter::List& params,
			   const NOX::Abstract::MultiVector* AA,
			   const NOX::Abstract::MultiVector::DenseMatrix* CC,
			   const NOX::Abstract::MultiVector* F,
			   const NOX::Abstract::MultiVector::DenseMatrix* G,
			   NOX::Abstract::MultiVector& X,
			   NOX::Abstract::MultiVector::DenseMatrix& Y) const;

      /*! 
       * \brief Solves the extended system when B is not 0 and F and A are
       * not contiguous
       */
      NOX::Abstract::Group::ReturnType  
      solveBNonZeroNoncontiguous(
			NOX::Parameter::List& params,
			const NOX::Abstract::MultiVector* AA,
			const LOCA::MultiContinuation::ConstraintInterface* BB,
			const NOX::Abstract::MultiVector::DenseMatrix* CC,
			const NOX::Abstract::MultiVector* F,
			const NOX::Abstract::MultiVector::DenseMatrix* G,
			NOX::Abstract::MultiVector& X,
			NOX::Abstract::MultiVector::DenseMatrix& Y) const;

      /*! 
       * \brief Solves the extended system when B is not 0 and F and A are
       * contiguous
       */
      NOX::Abstract::Group::ReturnType  
      solveBNonZeroContiguous(
			NOX::Parameter::List& params,
			const NOX::Abstract::MultiVector* AA,
			const LOCA::MultiContinuation::ConstraintInterface* BB,
			const NOX::Abstract::MultiVector::DenseMatrix* CC,
			vector<int>& indexF,
			vector<int>& indexA,
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
      Teuchos::RefCountPtr<LOCA::GlobalData> globalData;

      //! Solver parameters
      Teuchos::RefCountPtr<NOX::Parameter::List> solverParams;

      //! Pointer to group storing J
      Teuchos::RefCountPtr<const NOX::Abstract::Group>  grp;

      //! Pointer to A block
      Teuchos::RefCountPtr<const NOX::Abstract::MultiVector> A;

      //! Pointer to B block
      Teuchos::RefCountPtr<const LOCA::MultiContinuation::ConstraintInterface> B;

      //! Pointer to C block
      const NOX::Abstract::MultiVector::DenseMatrix *C;

      //! flag indicating whether A block is zero
      bool isZeroA;

      //! flag indicating whether B block is zero
      bool isZeroB;

      //! flag indicating whether C block is zero
      bool isZeroC;

      //! flag indicating whether F block is zero
      bool isZeroF;

      //! flag indicating whether G block is zero
      bool isZeroG;

      //! flag indicating whether F and A are contiguous
      bool isContiguous;

      //! flag indicating whether X block is zero
      bool isZeroX;

      //! flag indicating whether Y block is zero
      bool isZeroY;

      //! flag indicating whether B^T*X1 block is zero
      bool isZeroT1;

      //! flag indicating whether B^T*X2 block is zero
      bool isZeroT2;

    };
  } // namespace BorderedSystem
} // namespace LOCA

#endif
