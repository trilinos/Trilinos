// $Id$ 
// $Source$ 

//@HEADER
// ************************************************************************
// 
//            NOX: An Object-Oriented Nonlinear Solver Package
//                 Copyright (2002) Sandia Corporation
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
// Questions? Contact Tammy Kolda (tgkolda@sandia.gov) or Roger Pawlowski
// (rppawlo@sandia.gov), Sandia National Laboratories.
// 
// ************************************************************************
//@HEADER

#ifndef LOCA_BORDEREDSYSTEM_BORDERING_H
#define LOCA_BORDEREDSYSTEM_BORDERING_H

#include "LOCA_BorderedSystem_Generic.H"  // base class

namespace LOCA {

  namespace BorderedSystem {
    
    //! Solved bordered system of equations using bordering
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
    */
    class Bordering : public LOCA::BorderedSystem::Generic {

    public:

      //! Constructor. 
      Bordering(NOX::Parameter::List& params);

      //! Copy constructor
      Bordering(const Bordering& source);

      //! Destructor
      virtual ~Bordering();

      //! Assignment operator
      virtual Bordering& 
      operator = (const Bordering& source);

      //! Clone function
      virtual LOCA::BorderedSystem::Generic* clone() const;

      //! Assignment operator
      virtual LOCA::BorderedSystem::Generic& 
      operator = (const LOCA::BorderedSystem::Generic& source);

      //! Reset parameters
      virtual NOX::Abstract::Group::ReturnType 
      reset(NOX::Parameter::List& params);

      //! Set flags indicating if any blocks are zero
      virtual void setIsZero(bool flagA, bool flagB, bool flagC, bool flagF,
			     bool flagG);

      //! Set flag indicating whether F and A are continguous
      virtual void setIsContiguous(bool flag);

      //! Set blocks
      /*!
       * Any of these block pointers (execpt the group pointer) may be NULL 
       * if the corresponding flag is set indicating the block is zero.
       */
      virtual void setMatrixBlocks(
			const NOX::Abstract::Group* group,
			const NOX::Abstract::MultiVector* blockA,
			const NOX::Abstract::MultiVector* blockB,
			const NOX::Abstract::MultiVector::DenseMatrix* blockC);
  
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
       * \em isZeroF or \em \isZeroG is true, than the corresponding
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
			   const NOX::Abstract::MultiVector* A,
			   const NOX::Abstract::MultiVector::DenseMatrix* C,
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
			   const NOX::Abstract::MultiVector* A,
			   const NOX::Abstract::MultiVector* B,
			   const NOX::Abstract::MultiVector::DenseMatrix* C,
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
			  const NOX::Abstract::MultiVector* A,
			  const NOX::Abstract::MultiVector* B,
			  const NOX::Abstract::MultiVector::DenseMatrix* C,
			  vector<int>& indexF,
			  vector<int>& indexA,
			  const NOX::Abstract::MultiVector* F,
			  const NOX::Abstract::MultiVector::DenseMatrix* G,
			  NOX::Abstract::MultiVector& X,
			  NOX::Abstract::MultiVector::DenseMatrix& Y) const;

    protected:

      //! Pointer to group storing J
      const NOX::Abstract::Group* grp;

      //! Pointer to A block
      const NOX::Abstract::MultiVector* A;

      //! Pointer to B block
      const NOX::Abstract::MultiVector* B;

      //! Pointer to C block
      const NOX::Abstract::MultiVector::DenseMatrix* C;

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
