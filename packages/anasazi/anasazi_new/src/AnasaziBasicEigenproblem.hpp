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

#ifndef ANASAZI_BASIC_EIGENPROBLEM_H
#define ANASAZI_BASIC_EIGENPROBLEM_H

#include "AnasaziEigenproblem.hpp"
#include "AnasaziMultiVecTraits.hpp"
#include "AnasaziOperatorTraits.hpp"
#include "Teuchos_SerialDenseMatrix.hpp"
#include "Teuchos_SerialDenseVector.hpp"

/*! \class Anasazi::BasicEigenproblem
  \brief This provides a basic implementation for defining standard or 
  generalized eigenvalue problems.
*/

namespace Anasazi {
  
  template<class ScalarType, class MV, class OP>
  class BasicEigenproblem : public virtual Eigenproblem<ScalarType, MV, OP> {
    
  public:
    
    //@{ \name Constructors/Destructor.
    
    //! Empty constructor - allows Anasazi::BasicEigenproblem to be described at a later time through "Set Methods".
    BasicEigenproblem();
    
    //! Standard Eigenvalue Problem Constructor.
    BasicEigenproblem( const Teuchos::RefCountPtr<OP>& Op, const Teuchos::RefCountPtr<MV>& InitVec );
    
    //! Generalized Eigenvalue Problem Constructor.
    BasicEigenproblem( const Teuchos::RefCountPtr<OP>& Op, const Teuchos::RefCountPtr<OP>& B, const Teuchos::RefCountPtr<MV>& InitVec );
    
    //! Copy Constructor.
    BasicEigenproblem( const BasicEigenproblem<ScalarType, MV, OP>& Problem );	
    
    //! Destructor.
    virtual ~BasicEigenproblem() {};
    //@}
    
    //@{ \name Set Methods.
    
    /*! \brief Set the operator for which eigenvalues will be computed.  This may be different
      from the matrix \c A if a spectral transformation is employed, for example.  
    */
    void SetOperator( const Teuchos::RefCountPtr<OP>& Op ) { _Op = Op; };
    
    /*! \brief Set the operator A of the eigenvalue problem AX = BX\lambda.
    */
    void SetA( const Teuchos::RefCountPtr<OP>& A ) { _AOp = A; };
    
    /*! \brief Set the operator B of the eigenvalue problem AX = BX\lambda.
     */
    void SetB( const Teuchos::RefCountPtr<OP>& B ) { _BOp = B; }
    
    /*! \brief Set the preconditioner for this eigenvalue problem AX = BX\lambda.
     */
    void SetPrec( const Teuchos::RefCountPtr<OP>& Prec ) { _Prec = Prec; }
    
    /*! \brief Set the initial guess.  

    This vector is required to create all the space needed 
    by Anasazi to solve the eigenvalue problem.  Even if an initial guess is not known
    by the user, an initial vector must be passed in.  Sets the pointer to the input
    %Anasazi::MultiVec, so no copy is made.  
    */
    void SetInitVec( const Teuchos::RefCountPtr<MV>& InitVec ) { _InitVec = InitVec; };
    
    /*! \brief Set auxilliary vectors.

    NOTE:  This multivector can have any number of columns, an most likely will contain vectors that
    will be used by the eigensolver to orthogonalize against.
    */
    void SetAuxVec( const Teuchos::RefCountPtr<MV>& AuxVec ) { _AuxVec = AuxVec; };

    //! Inform the eigenproblem of the number of eigenvalues (NEV) that are required.
    void SetNEV( const int nev ){ _nev = nev; };

    //! Inform the eigenproblem that this problem is symmetric.
    /*! This knowledge may allow the solver to take advantage of the eigenproblems' symmetry.
      Some computational work can be avoided by setting this properly.
    */
    void SetSymmetric( const bool isSym ){ _isSym = isSym; };
    
    //! Inform the eigenproblem that is has all the information it needs to define the eigenproblem.
    /*! \note The user MUST call this routine before they send the eigenproblem to any solver!
     */
    ReturnType SetProblem();

    //@}
    
    //@{ \name Accessor Methods.
    
    //! Get a pointer to the Operator.
    Teuchos::RefCountPtr<OP> GetOperator() const { return( _Op ); };
    
    //! Get a pointer to the stiffness matrix.
    Teuchos::RefCountPtr<OP> GetA() const { return( _AOp ); };
    
    //! Get a pointer to the mass matrix.
    Teuchos::RefCountPtr<OP> GetB() const { return( _BOp ); };
    
    //! Get a pointer to the preconditioner.
    Teuchos::RefCountPtr<OP> GetPrec() const { return( _Prec ); };
    
    //! Get a pointer to the initial vector
    Teuchos::RefCountPtr<MV> GetInitVec() const { return( _InitVec ); };
    
    //! Get a pointer to the initial vector
    Teuchos::RefCountPtr<MV> GetAuxVec() const { return( _AuxVec ); };
    
    /*! \brief Get a pointer to the eigenvalues of the operator.
    
    NOTE:  If the operator is nonsymmetric, the length of this vector is 2*NEV where the 
    real part of eigenvalue \c j is entry \c j and the imaginary part is entry \c j+NEV .
    */
    Teuchos::RefCountPtr<std::vector<ScalarType> > GetEvals() { return( _Evals ); };
    
    /*! \brief Get a pointer to the eigenvectors of the operator.

    NOTE:  If the operator is nonsymmetric, this multivector has 2*NEV columns where the 
    real part of eigenvector \c j is column \c j and the imaginary part is column \c j+NEV .
    */
    Teuchos::RefCountPtr<MV> GetEvecs() { return( _Evecs ); };
    
    //! Get the number of eigenvalues (NEV) that are required by this eigenproblem.
    int GetNEV() const { return( _nev ); }

    //! Get the symmetry information for this eigenproblem.
    bool IsSymmetric() const { return( _isSym ); }
    
    //@}	
    
    //@{ \name Inner Product Methods.
    /*! \brief Computes inner product as needed by the eigensolver, for orthogonalization purposes.
     */
    ReturnType InnerProd( const MV& X, const MV& Y,
			  Teuchos::SerialDenseMatrix<int,ScalarType>& Z ) const;
    //@}

    //@{ \name Norm Methods.
    /*! \brief Computes the multivector norm as needed by the eigensolver.
      
    NOTE:  This can be different than the MvNorm method for the multivector class, which is 
    assumed to be the euclidean norm of each column.
    */
    ReturnType MvNorm( const MV& X, ScalarType* normvec ) const;
    
    //@}	
    
  protected:
    
    Teuchos::RefCountPtr<OP> _AOp, _BOp; 
    Teuchos::RefCountPtr<OP> _Op, _Prec;
    Teuchos::RefCountPtr<MV> _InitVec, _AuxVec;
    Teuchos::RefCountPtr<MV> _Evecs;
    Teuchos::RefCountPtr<std::vector<ScalarType> > _Evals;
    int _nev;
    bool _isSym;

    typedef MultiVecTraits<ScalarType,MV> MVT;
    typedef OperatorTraits<ScalarType,MV,OP> OPT;
  };		
  
  //=============================================================================
  //	Implementations (Constructors / Destructors)
  //=============================================================================
  
  template <class ScalarType, class MV, class OP>
  BasicEigenproblem<ScalarType, MV, OP>::BasicEigenproblem(void) : 
    _nev(1), 
    _isSym(false)
  {
  }
  
  //=============================================================================
    
  template <class ScalarType, class MV, class OP>
  BasicEigenproblem<ScalarType, MV, OP>::BasicEigenproblem( const Teuchos::RefCountPtr<OP>& Op, const Teuchos::RefCountPtr<MV>& InitVec ) :    
    _Op(Op), 
    _InitVec(InitVec), 
    _nev(1), 
    _isSym(false)
  {
  }
  
  //=============================================================================
  
  template <class ScalarType, class MV, class OP>
  BasicEigenproblem<ScalarType, MV, OP>::BasicEigenproblem( const Teuchos::RefCountPtr<OP>& Op, const Teuchos::RefCountPtr<OP>& B,
						      const Teuchos::RefCountPtr<MV>& InitVec ) :
    _BOp(B), 
    _Op(Op), 
    _InitVec(InitVec), 
    _nev(1), 
    _isSym(false)
  {
  }

  //=============================================================================
    
  template <class ScalarType, class MV, class OP>
  BasicEigenproblem<ScalarType, MV, OP>::BasicEigenproblem( const BasicEigenproblem<ScalarType,MV,OP>& Problem ) :
    _AOp(Problem._AOp), 
    _BOp(Problem._BOp), 
    _Op(Problem._Op), 
    _Prec(Problem._Prec), 
    _InitVec(Problem._InitVec), 
    _Evecs(Problem._Evecs),
    _Evals(Problem._Evals),
    _nev(Problem._nev), 
    _isSym(Problem._isSym)
  {
  }
  
  //=============================================================================
  //	SetProblem (sanity check method)
  //=============================================================================

  template <class ScalarType, class MV, class OP>
  ReturnType BasicEigenproblem<ScalarType, MV, OP>::SetProblem() 
  {
    //----------------------------------------------------------------
    // Sanity Checks
    //----------------------------------------------------------------
    // If there is no operator, then we can't proceed.
    if ( !_AOp.get() && !_Op.get() ) { return Failed; }

    // If there is no initial vector, then we don't have anything to clone workspace from.
    if ( !_InitVec.get() ) { return Failed; }

    // If we don't need any eigenvalues, we don't need to continue.
    if (_nev == 0) { return Failed; }

    // If there is an A, but no operator, we can set them equal.
    if (_AOp.get() && !_Op.get()) { _Op = _AOp; }

    // If this eigenproblem is being reused, then we may need to increase the space
    // for the eigenvalues / eigenvectors
    if (_Evecs.get()) {
      int old_nev = MVT::GetNumberVecs( *_Evecs );
      //
      // If the size of the old eigenproblem is larger than the new one, then
      // make sure all the storage required for the eigenproblem exists (check symmetry)
      //
      if ( _isSym ) {
	if ( _nev <= old_nev )
	  return Ok;
      } else {
	if ( 2*_nev <= old_nev )
	  return Ok;
      }
    }	
    
    //----------------------------------------------------------------
    // Allocate Memory
    // ( we need twice the storage if the problem is non-symmetric )
    //----------------------------------------------------------------
    if ( _isSym ) {      
      _Evecs = MVT::Clone( *_InitVec, _nev );
      _Evals = Teuchos::rcp( new std::vector<ScalarType>( _nev ) );
    } else {
      _Evecs = MVT::Clone( *_InitVec, 2*_nev );
      _Evals = Teuchos::rcp( new std::vector<ScalarType>( 2*_nev ) );
    }
    return Ok;
  }      
  
  
  //=============================================================================
  //	Implementations (Inner Product Methods)
  //=============================================================================
  
  template <class ScalarType, class MV, class OP>
  ReturnType BasicEigenproblem<ScalarType, MV, OP>::InnerProd( const MV& X, 
							       const MV& Y,
							       Teuchos::SerialDenseMatrix<int,ScalarType>& Z ) const
  {
    if ( _BOp.get() ) {
      Teuchos::RefCountPtr<MV> BY = MVT::CloneCopy( Y );
      
      // Apply B and check that it returned Ok.
      ReturnType ret = OPT::Apply( *_BOp, Y, *BY );
      if ( ret != Ok ) { return ret; }
      
      // Now perform inner product.  Result is stored in Z.
      MVT::MvTransMv( Teuchos::ScalarTraits<ScalarType>::one(), X, *BY, Z );
      
      return Ok;				
    } else {
      // Perform the inner product, assume B=I.
      MVT::MvTransMv( Teuchos::ScalarTraits<ScalarType>::one(), X, Y, Z );
      return Ok;				
    }
  }
  
  //=============================================================================
  //	Implementations (Norm Methods)
  //=============================================================================
  
  template <class ScalarType, class MV, class OP>
  ReturnType BasicEigenproblem<ScalarType, MV, OP>::MvNorm( const MV& X, ScalarType* normvec ) const
  {
    int IntOne = 1;
    int numvecs = MVT::GetNumberVecs( X );
    Teuchos::SerialDenseVector<int,ScalarType> DenseOne(IntOne);
    Teuchos::RefCountPtr<const MV> Xj;
    int *index = new int[IntOne];	
    ReturnType ret;
    
    for (int i=0; i<numvecs; i++) {
      index[0] = i;
      Xj = MVT::CloneView( X, index, IntOne );
      ret = InnerProd( *Xj, *Xj, DenseOne );
      if ( ret != Ok ) { return ret; }
      normvec[i] = sqrt(DenseOne(0));
    }
    
    delete [] index;
    
    return Ok;
  }

} // end Anasazi namespace
#endif

// end AnasaziBasicEigenproblem.hpp
