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

#ifndef ANASAZI_EIGENPROBLEM_H
#define ANASAZI_EIGENPROBLEM_H

#include "AnasaziMultiVec.hpp"
#include "AnasaziOperator.hpp"
#include "AnasaziReturnType.hpp"
#include "Teuchos_SerialDenseMatrix.hpp"
#include "Teuchos_SerialDenseVector.hpp"

/*! \class Anasazi::Eigenproblem
    \brief This class defines the interface required by an eigensolver and status
    test class to compute solutions to an eigenproblem.
*/


namespace Anasazi {
  
  template<class TYPE>
  class Eigenproblem {
    
  public:
    
    //@{ \name Constructors/Destructor.
    
    //! Empty constructor 
    Eigenproblem() {};
    
    //! Destructor.
    virtual ~Eigenproblem() {};
    //@}
    
    //@{ \name Set Methods.
    
    /*! \brief Set the operator for which eigenvalues will be computed.  

    NOTE:  This may be different from \c A if a spectral transformation is employed, for example.      
    */
    virtual void SetOperator( Operator<TYPE>* Op ) = 0;
    
    /*! \brief Set the operator A of the eigenvalue problem AX = BX\lambda.
    */
    virtual void SetA( Operator<TYPE>* A ) = 0;
    
    /*! \brief Set the operator B of the eigenvalue problem AX = BX\lambda.
    */
    virtual void SetB( Operator<TYPE>* B ) = 0;
    
    /*! \brief Set the preconditioner for this eigenvalue problem AX = BX\lambda.
     */
    virtual void SetPrec( Operator<TYPE>* Prec ) = 0;
    
    /*! \brief Set the initial guess.  

    NOTE:  This multivector should have the same number of columns as the blocksize.
    */
    virtual void SetInitVec( MultiVec<TYPE>* InitVec ) = 0; 
    
    /*! \brief Set auxilliary vectors.

    NOTE:  This multivector can have any number of columns, an most likely will contain vectors that
    will be used by the eigensolver to orthogonalize against.
    */
    virtual void SetAuxVec( MultiVec<TYPE>* AuxVec ) = 0;

    //! The number of eigenvalues (NEV) that are requested.
    virtual void SetNEV( const int nev ) = 0;

    //! Set the blocksize to be used by the iterative solver in solving this eigenproblem.
    virtual void SetBlockSize( const int blocksize ) = 0;

    //! Inform the eigenproblem that this problem is symmetric.
    /*! This knowledge may allow the solver to take advantage of the eigenproblems' symmetry.
      Some computational work can be avoided by setting this properly.
    */
    virtual void SetSymmetric( const bool isSym ) = 0;
    
    //! Inform the eigenproblem that is has all the information it needs to define the eigenproblem.
    /*! \note The user MUST call this routine before they send the eigenproblem to any solver!
     */
    virtual ReturnType SetProblem() = 0;

    //@}
    
    //@{ \name Accessor Methods.
    
    //! Get a pointer to the Operator.
    virtual Operator<TYPE>* GetOperator() const = 0;
    
    //! Get a pointer to the operator A of the eigenproblem AX = \lambda BX.
    virtual Operator<TYPE>* GetA() const = 0;
    
    //! Get a pointer to the operator B of the eigenproblem AX = \lambda BX.
    virtual Operator<TYPE>* GetB() const = 0;
    
    //! Get a pointer to the preconditioner.
    virtual Operator<TYPE>* GetPrec() const = 0;
    
    //! Get a pointer to the initial vector
    virtual MultiVec<TYPE>* GetInitVec() const = 0;
    
    //! Get a pointer to the auxilliary vector
    virtual MultiVec<TYPE>* GetAuxVec() const = 0;
    
    /*! \brief Get a pointer to the eigenvalues of the operator.
    
    NOTE:  If the operator is nonsymmetric, the length of this vector is 2*NEV where the 
    real part of eigenvalue \c j is entry \c j and the imaginary part is entry \c j+NEV .
    */
    virtual TYPE* GetEvals() = 0;
    
    /*! \brief Get a pointer to the eigenvectors of the operator.

    NOTE:  If the operator is nonsymmetric, this multivector has 2*NEV columns where the 
    real part of eigenvector \c j is column \c j and the imaginary part is column \c j+NEV .
    */
    virtual MultiVec<TYPE>* GetEvecs() = 0;
    
    //! Get the number of eigenvalues (NEV) that are required by this eigenproblem.
    virtual int GetNEV() const = 0;
    
    //! Get the blocksize to be used by the iterative solver in solving this eigenproblem.
    virtual int GetBlockSize() const = 0;
    
    //! Get the symmetry information for this eigenproblem.
    virtual bool IsSymmetric() const = 0;
    
    //@}	
    
    //@{ \name Inner Product Methods.
    /*! \brief Computes inner product as needed by the eigensolver, for orthogonalization purposes.
     */
    virtual ReturnType InnerProd( const MultiVec<TYPE>& X, const MultiVec<TYPE>& Y,
				  Teuchos::SerialDenseMatrix<int,TYPE>& Z ) = 0;
    //@}

    //@{ \name Norm Methods.
    /*! \brief Computes the multivector norm as needed by the eigensolver, for orthogonalization purposes.

    NOTE:  This can be different than the MvNorm method for the multivector class, which is 
    assumed to be the euclidean norm of each column.
     */
    virtual ReturnType MvNorm( MultiVec<TYPE>& X, TYPE* normvec ) = 0;
    
    //@}	
  };
   
} // end Anasazi namespace
#endif

// end AnasaziEigenproblem.hpp
