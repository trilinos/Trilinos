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

#include "AnasaziEpetraAdapter.hpp"

/*! \file AnasaziEpetraAdapter.cpp
 *   \brief Implementations of Anasazi multi-vector and operator classes using Epetra_MultiVector and Epetra_Operator classes
 *   */

namespace Anasazi {

  ///////////////////////////////////////////////////////////////////////////////
  //
  //--------Anasazi::EpetraMultiVec Implementation-------------------------------
  //
  ///////////////////////////////////////////////////////////////////////////////

  // Construction/Destruction
  
  EpetraMultiVec::EpetraMultiVec(const Epetra_BlockMap& Map, double * array, 
				 const int numvecs, const int stride)
    : Epetra_MultiVector(Copy, Map, array, stride, numvecs) 
  {
  }
  
  
  EpetraMultiVec::EpetraMultiVec(const Epetra_BlockMap& Map, const int numvecs)
    : Epetra_MultiVector(Map, numvecs) 
  {
  }
  
  
  EpetraMultiVec::EpetraMultiVec(Epetra_DataAccess CV, const Epetra_MultiVector& P_vec, 				
				 const std::vector<int>& index )
    : Epetra_MultiVector(CV, P_vec, &(const_cast<std::vector<int> &>(index))[0], index.size())
  {
  }
  
  
  EpetraMultiVec::EpetraMultiVec(const Epetra_MultiVector& P_vec)
    : Epetra_MultiVector(P_vec) 
  {
  }
  
  
  //
  //  member functions inherited from Anasazi::MultiVec
  //
  //
  //  Simulating a virtual copy constructor. If we could rely on the co-variance
  //  of virtual functions, we could return a pointer to EpetraMultiVec
  //  (the derived type) instead of a pointer to the pure virtual base class.
  //
  
  MultiVec<double>* EpetraMultiVec::Clone ( const int numvecs ) const
  {
    EpetraMultiVec * ptr_apv = new EpetraMultiVec(Map(), numvecs);
    return ptr_apv; // safe upcast.
  }
  //
  //  the following is a virtual copy constructor returning
  //  a pointer to the pure virtual class. vector values are
  //  copied.
  //
  
  MultiVec<double>* EpetraMultiVec::CloneCopy() const
  {
    EpetraMultiVec *ptr_apv = new EpetraMultiVec(*this);
    return ptr_apv; // safe upcast
  }
  
  
  MultiVec<double>* EpetraMultiVec::CloneCopy ( const std::vector<int>& index ) const
  {
    EpetraMultiVec * ptr_apv = new EpetraMultiVec(Copy, *this, index);
    return ptr_apv; // safe upcast.
  }
  
  
  MultiVec<double>* EpetraMultiVec::CloneView ( const std::vector<int>& index ) 
  {
    EpetraMultiVec * ptr_apv = new EpetraMultiVec(View, *this, index);
    return ptr_apv; // safe upcast.
  }
  
  
  void EpetraMultiVec::SetBlock( const MultiVec<double>& A, const std::vector<int>& index ) 
  {	
    EpetraMultiVec temp_vec(View, *this, index);

    int numvecs = index.size();
    if ( A.GetNumberVecs() != numvecs ) {
      std::vector<int> index2( numvecs );
      for(int i=0; i<numvecs; i++)
	index2[i] = i;
      EpetraMultiVec *tmp_vec = dynamic_cast<EpetraMultiVec *>(&const_cast<MultiVec<double> &>(A)); 
      assert(tmp_vec!=NULL);
      EpetraMultiVec A_vec(View, *tmp_vec, index2);
      temp_vec.MvAddMv( 1.0, A_vec, 0.0, A_vec );
    }
    else {
      temp_vec.MvAddMv( 1.0, A, 0.0, A );
    }
  }								

  //-------------------------------------------------------------
  //
  // *this <- alpha * A * B + beta * (*this)
  //
  //-------------------------------------------------------------
  
  void EpetraMultiVec::MvTimesMatAddMv ( const double alpha, const MultiVec<double>& A, 
					 const Teuchos::SerialDenseMatrix<int,double>& B, const double beta ) 
  {
    Epetra_LocalMap LocalMap(B.numRows(), 0, Map().Comm());
    Epetra_MultiVector B_Pvec(Copy, LocalMap, B.values(), B.stride(), B.numCols());
    
    EpetraMultiVec *A_vec = dynamic_cast<EpetraMultiVec *>(&const_cast<MultiVec<double> &>(A)); 
    assert(A_vec!=NULL);
    
    int ret = Multiply( 'N', 'N', alpha, *A_vec, B_Pvec, beta );
    assert( ret == 0 );
  }

  //-------------------------------------------------------------
  //
  // *this <- alpha * A + beta * B
  //
  //-------------------------------------------------------------
  
  void EpetraMultiVec::MvAddMv ( const double alpha , const MultiVec<double>& A, 
				 const double beta, const MultiVec<double>& B) 
  {
    EpetraMultiVec *A_vec = dynamic_cast<EpetraMultiVec *>(&const_cast<MultiVec<double> &>(A)); 
    assert(A_vec!=NULL);
    EpetraMultiVec *B_vec = dynamic_cast<EpetraMultiVec *>(&const_cast<MultiVec<double> &>(B)); 
    assert(B_vec!=NULL);
    
    int ret = Update( alpha, *A_vec, beta, *B_vec, 0.0 );
    assert( ret == 0 ); 
  }

  //-------------------------------------------------------------
  //
  // dense B <- alpha * A^T * (*this)
  //
  //-------------------------------------------------------------
  
  void EpetraMultiVec::MvTransMv ( const double alpha, const MultiVec<double>& A,
				   Teuchos::SerialDenseMatrix<int,double>& B
#ifdef HAVE_ANASAZI_EXPERIMENTAL
				   , ConjType conj
#endif
				   ) const
  {    
    EpetraMultiVec *A_vec = dynamic_cast<EpetraMultiVec *>(&const_cast<MultiVec<double> &>(A));
    
    if (A_vec) {
      Epetra_LocalMap LocalMap(B.numRows(), 0, Map().Comm());
      Epetra_MultiVector B_Pvec(View, LocalMap, B.values(), B.stride(), B.numCols());
      
      int ret = B_Pvec.Multiply( 'T', 'N', alpha, *A_vec, *this, 0.0 );
      assert( ret == 0 ); 
    }
  }
  
  //-------------------------------------------------------------
  //
  // b[i] = A[i]^T * this[i]
  // 
  //-------------------------------------------------------------
  
  void EpetraMultiVec::MvDot ( const MultiVec<double>& A, std::vector<double>* b
#ifdef HAVE_ANASAZI_EXPERIMENTAL
			       , ConjType conj
#endif
			       ) const
  {
    EpetraMultiVec *A_vec = dynamic_cast<EpetraMultiVec *>(&const_cast<MultiVec<double> &>(A)); 
    assert(A_vec!=NULL);
    if ((A_vec!=NULL) && (b!=NULL) && ( (int)b->size() >= A_vec->NumVectors() ) ) {
      int ret = this->Dot( *A_vec, &(*b)[0] );
      assert( ret == 0 );
    }
  }

  ///////////////////////////////////////////////////////////////////////////////
  //
  //--------Anasazi::EpetraOp Implementation-------------------------------------
  //
  ///////////////////////////////////////////////////////////////////////////////

  //
  // AnasaziOperator constructors
  //  
  EpetraOp::EpetraOp(const Teuchos::RefCountPtr<Epetra_Operator> &Op) 
    : Epetra_Op(Op)
  {
  }
  
  EpetraOp::~EpetraOp() 
  {
  }
  //
  // AnasaziOperator applications
  //
  void EpetraOp::Apply ( const MultiVec<double>& X, 
			       MultiVec<double>& Y ) const 
  {
    //
    // This standard operator computes Y = A*X
    //
    MultiVec<double> & temp_X = const_cast<MultiVec<double> &>(X);
    Epetra_MultiVector* vec_X = dynamic_cast<Epetra_MultiVector* >(&temp_X);
    Epetra_MultiVector* vec_Y = dynamic_cast<Epetra_MultiVector* >(&Y);
    
    assert( vec_X!=NULL && vec_Y!=NULL );

    int info = Epetra_Op->Apply( *vec_X, *vec_Y );
    TEST_FOR_EXCEPTION( info != 0, OperatorError, 
                        "Anasazi::EpetraOp::Apply(): Error returned from Epetra_Operator::Apply()" );
  }

  ///////////////////////////////////////////////////////////////////////////////
  //
  //--------Anasazi::EpetraGenOp Implementation----------------------------------
  //
  ///////////////////////////////////////////////////////////////////////////////

  //
  // AnasaziOperator constructors
  //
  
  EpetraGenOp::EpetraGenOp(const Teuchos::RefCountPtr<Epetra_Operator> &AOp,
			   const Teuchos::RefCountPtr<Epetra_Operator> &MOp,
			   bool isAInverse_) 
    : isAInverse( isAInverse_ ), Epetra_AOp(AOp), Epetra_MOp(MOp) 
  {
  }
    
  EpetraGenOp::~EpetraGenOp() 
  {
  }
  //
  // AnasaziOperator applications
  //
  void EpetraGenOp::Apply ( const MultiVec<double>& X, MultiVec<double>& Y ) const 
  {
    //
    // This generalized operator computes Y = A^{-1}*M*X
    //
    int info=0;
    MultiVec<double> & temp_X = const_cast<MultiVec<double> &>(X);
    Epetra_MultiVector* vec_X = dynamic_cast<Epetra_MultiVector* >(&temp_X);
    Epetra_MultiVector* vec_Y = dynamic_cast<Epetra_MultiVector* >(&Y);
    Epetra_MultiVector temp_Y(*vec_Y); 
    
    assert( vec_X!=NULL && vec_Y!=NULL );
    //
    // Need to cast away constness because the member function Apply is not declared const.  
    // Change the transpose setting for the operator if necessary and change it back when done.
    //
    // Apply M
    info = Epetra_MOp->Apply( *vec_X, temp_Y );
    TEST_FOR_EXCEPTION( info != 0, OperatorError, 
                        "Anasazi::EpetraGenOp::Apply(): Error returned from Epetra_Operator::Apply()" );
    // Apply A or A^{-1}
    if (isAInverse) {
      info = Epetra_AOp->ApplyInverse( temp_Y, *vec_Y );
    }
    else {
      info = Epetra_AOp->Apply( temp_Y, *vec_Y );
    }
    TEST_FOR_EXCEPTION( info != 0, OperatorError, 
                        "Anasazi::EpetraGenOp::Apply(): Error returned from Epetra_Operator::Apply()" );
  }
  
  int EpetraGenOp::Apply(const Epetra_MultiVector &X, Epetra_MultiVector &Y) const
  {
    //
    // This generalized operator computes Y = A^{-1}*M*X 
    //
    int info=0;
    Epetra_MultiVector temp_Y(OperatorDomainMap(), Y.NumVectors()); 
    
    // Apply M
    info = Epetra_MOp->Apply( X, temp_Y );
    if (info!=0) return info;
    
    // Apply A or A^{-1}
    if (isAInverse)
      info = Epetra_AOp->ApplyInverse( temp_Y, Y );
    else
      info = Epetra_AOp->Apply( temp_Y, Y );

    return info;
  }
  
  int EpetraGenOp::ApplyInverse(const Epetra_MultiVector &X, Epetra_MultiVector &Y) const
  {
    //
    // This generalized operator computes Y = M^{-1}*A*X 
    //
    int info=0;
    Epetra_MultiVector temp_Y(OperatorDomainMap(), Y.NumVectors()); 
    
    // Apply A or A^{-1}
    if (isAInverse)
      info = Epetra_AOp->Apply( X, temp_Y );
    else
      info = Epetra_AOp->ApplyInverse( X, temp_Y );
    if (info!=0) return info;
    
    // Apply M^{-1}
    info = Epetra_MOp->ApplyInverse( temp_Y, Y );
    
    return info;
  }

  ///////////////////////////////////////////////////////////////////////////////
  //
  //--------Anasazi::EpetraSymOp Implementation----------------------------------
  //
  ///////////////////////////////////////////////////////////////////////////////

  //
  // AnasaziOperator constructors
  //
  EpetraSymOp::EpetraSymOp(const Teuchos::RefCountPtr<Epetra_Operator> &Op, 
			   const bool isTrans) 
    : Epetra_Op(Op), isTrans_(isTrans)
  {
  }
  
  EpetraSymOp::~EpetraSymOp() 
  {
  }
  //
  // AnasaziOperator applications
  //
  void EpetraSymOp::Apply ( const MultiVec<double>& X, 
				  MultiVec<double>& Y ) const 
  {
    int info=0;
    MultiVec<double> & temp_X = const_cast<MultiVec<double> &>(X);
    Epetra_MultiVector* vec_X = dynamic_cast<Epetra_MultiVector* >(&temp_X);
    Epetra_MultiVector* vec_Y = dynamic_cast<Epetra_MultiVector* >(&Y);
    Epetra_MultiVector* temp_vec = new Epetra_MultiVector( 
							  (isTrans_) ? Epetra_Op->OperatorDomainMap() 
							  : Epetra_Op->OperatorRangeMap(), 
							  vec_X->NumVectors() );
    
    assert( vec_X!=NULL && vec_Y!=NULL && temp_vec!=NULL );
    //
    // Need to cast away constness because the member function Apply
    // is not declared const.
    //
    // Transpose the operator (if isTrans_ = true)
    if (isTrans_) {
      info = Epetra_Op->SetUseTranspose( isTrans_ );
      if (info != 0) {
        delete temp_vec;
        TEST_FOR_EXCEPTION( true, OperatorError, 
                            "Anasazi::EpetraSymOp::Apply(): Error returned from Epetra_Operator::Apply()" );
      }
    }
    //
    // Compute A*X or A'*X 
    //
    info=Epetra_Op->Apply( *vec_X, *temp_vec );
    if (info!=0) { 
      delete temp_vec; 
      TEST_FOR_EXCEPTION( true, OperatorError, 
                          "Anasazi::EpetraSymOp::Apply(): Error returned from Epetra_Operator::Apply()" );
    }
    //
    // Transpose/Un-transpose the operator based on value of isTrans_
    info=Epetra_Op->SetUseTranspose( !isTrans_ );
    if (info!=0) { 
      delete temp_vec; 
      TEST_FOR_EXCEPTION( true, OperatorError, 
                          "Anasazi::EpetraSymOp::Apply(): Error returned from Epetra_Operator::Apply()" );
    }
    
    // Compute A^T*(A*X) or A*A^T
    info=Epetra_Op->Apply( *temp_vec, *vec_Y );
    if (info!=0) { 
      delete temp_vec; 
      TEST_FOR_EXCEPTION( true, OperatorError, 
                          "Anasazi::EpetraSymOp::Apply(): Error returned from Epetra_Operator::Apply()" );
    }
    
    // Un-transpose the operator
    info=Epetra_Op->SetUseTranspose( false );
    delete temp_vec;
    TEST_FOR_EXCEPTION( info != 0, OperatorError, 
                        "Anasazi::EpetraSymOp::Apply(): Error returned from Epetra_Operator::Apply()" );
  }
  
  int EpetraSymOp::Apply(const Epetra_MultiVector &X, Epetra_MultiVector &Y) const
  {
    int info=0;
    Epetra_MultiVector temp_vec(OperatorDomainMap(), Y.NumVectors()); 
    //
    // This generalized operator computes Y = A^T*A*X or Y = A*A^T*X
    //
    // Transpose the operator (if isTrans_ = true)
    if (isTrans_) {
      info=Epetra_Op->SetUseTranspose( isTrans_ );
      if (info!=0) { return info; }
    }
    //
    // Compute A*X or A^T*X 
    //
    info=Epetra_Op->Apply( X, temp_vec );
    if (info!=0) { return info; }
    //
    // Transpose/Un-transpose the operator based on value of isTrans_
    info=Epetra_Op->SetUseTranspose( !isTrans_ );
    if (info!=0) { return info; }
    
    // Compute A^T*(A*X) or A*A^T
    info=Epetra_Op->Apply( temp_vec, Y );
    if (info!=0) { return info; }
    
    // Un-transpose the operator
    info=Epetra_Op->SetUseTranspose( false );
    return info;
  }
  
  int EpetraSymOp::ApplyInverse(const Epetra_MultiVector &X, Epetra_MultiVector &Y) const
  {
    int info=0;
    Epetra_MultiVector temp_vec(OperatorDomainMap(), Y.NumVectors()); 
    //
    // This generalized operator computes Y = (A^T*A)^{-1}*X or Y = (A*A^T)^{-1}*X
    //
    // Transpose the operator (if isTrans_ = true)
    if (!isTrans_) {
      info=Epetra_Op->SetUseTranspose( !isTrans_ );
      if (info!=0) { return info; }
    }
    //
    // Compute A^{-1}*X or A^{-T}*X 
    //
    info=Epetra_Op->ApplyInverse( X, temp_vec );
    if (info!=0) { return info; }
    //
    // Transpose/Un-transpose the operator based on value of isTrans_
    info=Epetra_Op->SetUseTranspose( isTrans_ );
    if (info!=0) { return info; }
    
    // Compute A^{-T}*(A^{-1}*X) or A^{-1}*A^{-T}
    info=Epetra_Op->Apply( temp_vec, Y );
    if (info!=0) { return info; }
    
    // Un-transpose the operator
    info=Epetra_Op->SetUseTranspose( false );
    return info;
  }

  ///////////////////////////////////////////////////////////////////////////////
  //
  //--------Anasazi::EpetraSymMVOp Implementation--------------------------------
  //
  ///////////////////////////////////////////////////////////////////////////////

  //
  // Anasazi::Operator constructors
  //
  EpetraSymMVOp::EpetraSymMVOp(const Teuchos::RefCountPtr<Epetra_MultiVector> &MV, const bool isTrans) 
    : Epetra_MV(MV), isTrans_(isTrans)
  {
    if (isTrans)
      MV_localmap = Teuchos::rcp( new Epetra_LocalMap( Epetra_MV->NumVectors(), 0, Epetra_MV->Map().Comm() ) );
    else
      MV_blockmap = Teuchos::rcp( &Epetra_MV->Map(), false );
  }
  
  //
  // AnasaziOperator applications
  //
  void EpetraSymMVOp::Apply ( const MultiVec<double>& X, MultiVec<double>& Y ) const 
  {
    int info=0;
    MultiVec<double> & temp_X = const_cast<MultiVec<double> &>(X);
    Epetra_MultiVector* vec_X = dynamic_cast<Epetra_MultiVector* >(&temp_X);
    Epetra_MultiVector* vec_Y = dynamic_cast<Epetra_MultiVector* >(&Y);
    
    if (isTrans_) {
      
      Epetra_MultiVector temp_vec( *MV_localmap, temp_X.GetNumberVecs() );
      
      /* A'*X */
      info = temp_vec.Multiply( 'T', 'N', 1.0, *Epetra_MV, *vec_X, 0.0 );
      TEST_FOR_EXCEPTION( info != 0, OperatorError, 
                          "Anasazi::EpetraSymMVOp::Apply(): Error returned from Epetra_Operator::Apply()" );
      
      /* A*(A'*X) */
      info = vec_Y->Multiply( 'N', 'N', 1.0, *Epetra_MV, temp_vec, 0.0 );      
      TEST_FOR_EXCEPTION( info != 0, OperatorError, 
                          "Anasazi::EpetraSymMVOp::Apply(): Error returned from Epetra_Operator::Apply()" );
    } 
    else {
      
      Epetra_MultiVector temp_vec( *MV_blockmap, temp_X.GetNumberVecs() );
      
      /* A*X */
      info = temp_vec.Multiply( 'N', 'N', 1.0, *Epetra_MV, *vec_X, 0.0 );
      TEST_FOR_EXCEPTION( info != 0, OperatorError, 
                          "Anasazi::EpetraSymMVOp::Apply(): Error returned from Epetra_Operator::Apply()" );
      
      /* A'*(A*X) */
      info = vec_Y->Multiply( 'T', 'N', 1.0, *Epetra_MV, temp_vec, 0.0 );
      TEST_FOR_EXCEPTION( info != 0, OperatorError, 
                          "Anasazi::EpetraSymMVOp::Apply(): Error returned from Epetra_Operator::Apply()" );
    }
  }

} // end namespace Anasazi
