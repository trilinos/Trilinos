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

#ifndef ANASAZI_EPETRA_ADAPTER_HPP
#define ANASAZI_EPETRA_ADAPTER_HPP

#include "AnasaziMultiVec.hpp"
#include "AnasaziOperator.hpp"
#include "AnasaziConfigDefs.hpp"
#include "AnasaziReturnType.hpp"

#include "Teuchos_SerialDenseMatrix.hpp"
#include "Epetra_MultiVector.h"
#include "Epetra_Operator.h"
#include "Epetra_Map.h"
#include "Epetra_LocalMap.h"

namespace Anasazi {
  
  //--------template class AnasaziEpetraMultiVec-------------------------------------
  class EpetraMultiVec : public MultiVec<double>, public Epetra_MultiVector {
  public:
    // constructors
    EpetraMultiVec(const Epetra_BlockMap& Map, double * array, const int numvecs, const int stride=0);
    EpetraMultiVec(const Epetra_BlockMap& Map, const int numvecs);
    EpetraMultiVec(Epetra_DataAccess CV, const Epetra_MultiVector& P_vec, int index[], int numvecs );
    EpetraMultiVec(const Epetra_MultiVector & P_vec);
    ~EpetraMultiVec();
    //
    //  member functions inherited from Anasazi::MultiVec
    //
    //  the following is a virtual copy constructor returning
    //  a pointer to the pure virtual class. vector values are
    //  not copied; instead a new MultiVec is created containing
    //  a non-zero amount of columns.
    //
    MultiVec<double> * Clone ( const int numvecs ) const;
    //
    //  the following is a virtual copy constructor returning
    //  a pointer to the pure virtual class. vector values are
    //  copied and a new stand-alone MultiVector is created.
    //  (deep copy).
    //
    MultiVec<double> * CloneCopy () const;
    //
    //  the following is a virtual copy constructor returning
    //  a pointer to the pure virtual class. vector values are
    //  copied and a new stand-alone MultiVector is created
    //  where only selected columns are chosen.  (deep copy).
    //
    MultiVec<double> * CloneCopy ( int index[], int numvecs ) const;
    //
    //  the following is a virtual view constructor returning
    //  a pointer to the pure virtual class. vector values are 
    //  shared and hence no memory is allocated for the columns.
    //
    MultiVec<double> * CloneView ( int index[], int numvecs);
    //
    //  this routine sets a subblock of the multivector, which
    //  need not be contiguous, and is given by the indices.
    //
    void SetBlock ( const MultiVec<double>& A, int index[], int numvecs );
    //
    int GetNumberVecs () const { return NumVectors(); }
    int GetVecLength () const { return GlobalLength(); }
    //
    // *this <- alpha * A * B + beta * (*this)
    //
    void MvTimesMatAddMv ( double alpha, const MultiVec<double>& A, 
			   const Teuchos::SerialDenseMatrix<int,double>& B, double beta );
    //
    // *this <- alpha * A + beta * B
    //
    void MvAddMv ( double alpha , const MultiVec<double>& A, double beta,
		   const MultiVec<double>& B);
    //
    // B <- alpha * A^T * (*this)
    //
    void MvTransMv ( double alpha, const MultiVec<double>& A, Teuchos::SerialDenseMatrix<int,double>& B ) const;
    //
    // b[i] = A[i]^T * this[i]
    // 
    void MvDot ( const MultiVec<double>& A, double b[] ) const;
    //
    // alpha[i] = norm of i-th column of (*this)
    //	
    void MvNorm ( double * normvec ) const {
      if (normvec)
	assert( Norm2(normvec) == 0 );
    };
    //
    // random vectors in i-th column of (*this)
    //
    void MvRandom() { assert( Random() == 0 ); };
    //
    // initializes each element of (*this) with alpha
    //
    void MvInit ( double alpha ) { assert( PutScalar( alpha ) == 0 ); };
    //
    // print (*this)
    //
    void MvPrint() const { cout<< *this << endl; };
  private:
  };
  //-------------------------------------------------------------
  
  //////////////////////////////////////////////////////////////////////
  // Construction/Destruction
  //////////////////////////////////////////////////////////////////////
  
  
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
				 int index[], int numvecs )
    : Epetra_MultiVector(CV, P_vec, index, numvecs) 
  {
  }
  
  
  EpetraMultiVec::EpetraMultiVec(const Epetra_MultiVector& P_vec)
    : Epetra_MultiVector(P_vec) 
  {
  }
  
  
  EpetraMultiVec::~EpetraMultiVec() 
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
  
  
  MultiVec<double>* EpetraMultiVec::CloneCopy ( int index[], int numvecs ) const
  {
    EpetraMultiVec * ptr_apv = new EpetraMultiVec(Copy, *this, index, numvecs );
    return ptr_apv; // safe upcast.
  }
  
  
  MultiVec<double>* EpetraMultiVec::CloneView ( int index[], int numvecs ) 
  {
    EpetraMultiVec * ptr_apv = new EpetraMultiVec(View, *this, index, numvecs );
    return ptr_apv; // safe upcast.
  }
  
  
  void EpetraMultiVec::SetBlock(const MultiVec<double>& A, int index[], int numvecs ) 
  {	
    EpetraMultiVec temp_vec(View, *this, index, numvecs);

    if ( A.GetNumberVecs() != numvecs ) {
      std::vector<int> index2( numvecs );
      for(int i=0; i<numvecs; i++)
	index2[i] = i;
      EpetraMultiVec *tmp_vec = dynamic_cast<EpetraMultiVec *>(&const_cast<MultiVec<double> &>(A)); assert(tmp_vec!=NULL);
      EpetraMultiVec A_vec(View, *tmp_vec, &index2[0], numvecs );
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
  
  void EpetraMultiVec::MvTimesMatAddMv ( double alpha, const MultiVec<double>& A, 
				    const Teuchos::SerialDenseMatrix<int,double>& B, double beta ) 
  {
    Epetra_LocalMap LocalMap(B.numRows(), 0, Map().Comm());
    Epetra_MultiVector B_Pvec(Copy, LocalMap, B.values(), B.stride(), B.numCols());
    
    EpetraMultiVec *A_vec = dynamic_cast<EpetraMultiVec *>(&const_cast<MultiVec<double> &>(A)); assert(A_vec!=NULL);
    
    assert( Multiply( 'N', 'N', alpha, *A_vec, B_Pvec, beta ) == 0 );
  }

  //-------------------------------------------------------------
  //
  // *this <- alpha * A + beta * B
  //
  //-------------------------------------------------------------
  
  void EpetraMultiVec::MvAddMv ( double alpha , const MultiVec<double>& A, 
			    double beta, const MultiVec<double>& B) 
  {
    EpetraMultiVec *A_vec = dynamic_cast<EpetraMultiVec *>(&const_cast<MultiVec<double> &>(A)); assert(A_vec!=NULL);
    EpetraMultiVec *B_vec = dynamic_cast<EpetraMultiVec *>(&const_cast<MultiVec<double> &>(B)); assert(B_vec!=NULL);
    
    assert ( Update( alpha, *A_vec, beta, *B_vec, 0.0 ) == 0 ); 
  }

  //-------------------------------------------------------------
  //
  // dense B <- alpha * A^T * (*this)
  //
  //-------------------------------------------------------------
  
  void EpetraMultiVec::MvTransMv ( double alpha, const MultiVec<double>& A,
			      Teuchos::SerialDenseMatrix<int,double>& B) const
  {    
    EpetraMultiVec *A_vec = dynamic_cast<EpetraMultiVec *>(&const_cast<MultiVec<double> &>(A));
    
    if (A_vec) {
      Epetra_LocalMap LocalMap(B.numRows(), 0, Map().Comm());
      Epetra_MultiVector B_Pvec(View, LocalMap, B.values(), B.stride(), B.numCols());
      
      assert ( B_Pvec.Multiply( 'T', 'N', alpha, *A_vec, *this, 0.0 ) == 0 ); 
    }
  }
  
  //-------------------------------------------------------------
  //
  // b[i] = A[i]^T * this[i]
  // 
  //-------------------------------------------------------------
  
  void EpetraMultiVec::MvDot ( const MultiVec<double>& A, double b[] ) const
  {
    EpetraMultiVec *A_vec = dynamic_cast<EpetraMultiVec *>(&const_cast<MultiVec<double> &>(A)); assert(A_vec!=NULL);
    if (A_vec && b) {
      assert( this->Dot( *A_vec, b ) == 0 );
    }
  }
  
  ///////////////////////////////////////////////////////////////
  //--------template class AnasaziEpetraOp---------------------
  
  class EpetraOp : public virtual Operator<double> {
  public:
    EpetraOp(const Epetra_Operator& );
    ~EpetraOp();
    ReturnType Apply ( const MultiVec<double>& x, 
		       MultiVec<double>& y ) const;
  private:
    const Epetra_Operator & Epetra_Op;
  };
  //-------------------------------------------------------------
  //
  // implementation of the AnasaziEpetraOp class.
  //
  ////////////////////////////////////////////////////////////////////
  //
  // AnasaziOperator constructors
  //
  
  EpetraOp::EpetraOp(const Epetra_Operator& Op) 
    : Epetra_Op(Op)
  {
  }
  
  EpetraOp::~EpetraOp() 
  {
  }
  //
  // AnasaziOperator applications
  //
  ReturnType EpetraOp::Apply ( const MultiVec<double>& x, 
			       MultiVec<double>& y ) const 
  {
    //
    // This standard operator computes y = A*x
    //
    MultiVec<double> & temp_x = const_cast<MultiVec<double> &>(x);
    Epetra_MultiVector* vec_x = dynamic_cast<Epetra_MultiVector* >(&temp_x);
    Epetra_MultiVector* vec_y = dynamic_cast<Epetra_MultiVector* >(&y);
    
    assert( vec_x!=NULL && vec_y!=NULL );

    int info = Epetra_Op.Apply( *vec_x, *vec_y );
    
    if (info==0) { 
      return Ok; 
    } else { 
      return Failed; 
    }	
  }
  
  ///////////////////////////////////////////////////////////////
  //--------template class AnasaziEpetraGenOp---------------------
  
  class EpetraGenOp : public virtual Operator<double> {
  public:
    EpetraGenOp(const Epetra_Operator&, const Epetra_Operator& );
    ~EpetraGenOp();
    ReturnType Apply ( const MultiVec<double>& x, 
		       MultiVec<double>& y ) const;
  private:
    const Epetra_Operator & Epetra_AOp;
    const Epetra_Operator & Epetra_BOp;
  };
  //-------------------------------------------------------------
  //
  // implementation of the AnasaziEpetraGenOp class.
  //
  ////////////////////////////////////////////////////////////////////
  //
  // AnasaziOperator constructors
  //
  
  EpetraGenOp::EpetraGenOp(const Epetra_Operator& AOp,
			   const Epetra_Operator& BOp) 
    : Epetra_AOp(AOp), Epetra_BOp(BOp) 
  {
  }
  
  
  EpetraGenOp::~EpetraGenOp() 
  {
  }
  //
  // AnasaziOperator applications
  //
  ReturnType EpetraGenOp::Apply ( const MultiVec<double>& x, 
				  MultiVec<double>& y ) const 
  {
    //
    // This generalized operator computes y = A*B*x of y = (A*B)^T*x
    //
    int info=0;
    MultiVec<double> & temp_x = const_cast<MultiVec<double> &>(x);
    Epetra_MultiVector* vec_x = dynamic_cast<Epetra_MultiVector* >(&temp_x);
    Epetra_MultiVector* vec_y = dynamic_cast<Epetra_MultiVector* >(&y);
    Epetra_MultiVector temp_y(*vec_y); 
    
    assert( vec_x!=NULL && vec_y!=NULL );
    //
    // Need to cast away constness because the member function Apply is not declared const.  
    // Change the transpose setting for the operator if necessary and change it back when done.
    //
    // Apply B
    info = Epetra_BOp.Apply( *vec_x, temp_y );
    assert(info==0);
    // Apply A
    info = Epetra_AOp.Apply( temp_y, *vec_y );
    if (info==0) { 
      return Ok; 
    } else { 
      return Failed; 
    }	
  }
  
  ///////////////////////////////////////////////////////////////
  //--------template class AnasaziEpetraSymOp---------------------
  class EpetraSymOp : public virtual Operator<double> {
  public:
    EpetraSymOp(const Epetra_Operator& Op );
    ~EpetraSymOp();
    ReturnType Apply ( const MultiVec<double>& x, 
		       MultiVec<double>& y ) const;
  private:
    const Epetra_Operator& Epetra_Op;
  };
  //-------------------------------------------------------------
  //
  // implementation of the AnasaziEpetraSymOp class.
  //
  ////////////////////////////////////////////////////////////////////
  //
  // AnasaziOperator constructors
  //
  EpetraSymOp::EpetraSymOp(const Epetra_Operator& Op) 
    : Epetra_Op(Op)
  {
  }
  
  EpetraSymOp::~EpetraSymOp() 
  {
  }
  //
  // AnasaziOperator applications
  //
  ReturnType EpetraSymOp::Apply ( const MultiVec<double>& x, 
				  MultiVec<double>& y ) const 
  {
    int info=0;
    MultiVec<double> & temp_x = const_cast<MultiVec<double> &>(x);
    Epetra_MultiVector* vec_x = dynamic_cast<Epetra_MultiVector* >(&temp_x);
    Epetra_MultiVector* vec_y = dynamic_cast<Epetra_MultiVector* >(&y);
    Epetra_MultiVector temp_vec( Epetra_Op.OperatorRangeMap(), vec_x->NumVectors() );
    
    assert( vec_x!=NULL && vec_y!=NULL );
    //
    // Need to cast away constness because the member function Apply
    // is not declared const.
    //
    // Compute A*x
    info = Epetra_Op.Apply( *vec_x, temp_vec );
    if (info!=0) { return Failed; }
    
    // Transpose the operator
    info = const_cast<Epetra_Operator&>(Epetra_Op).SetUseTranspose( true );
    if (info!=0) { return Failed; }
    
    // Compute A^T*(A*x)
    info = Epetra_Op.Apply( temp_vec, *vec_y );
    if (info!=0) { return Failed; }
    
    // Un-transpose the operator
    info = const_cast<Epetra_Operator&>(Epetra_Op).SetUseTranspose( false );
    
    if (info==0)
      return Ok; 
    else
      return Failed; 
  }
  
  ////////////////////////////////////////////////////////////////////
  //
  // Implementation of the Anasazi::MultiVecTraits for Epetra::MultiVector.
  //
  ////////////////////////////////////////////////////////////////////

  template<>
  class MultiVecTraits<double, Epetra_MultiVector>
  {
  public:
    ///
    static Teuchos::RefCountPtr<Epetra_MultiVector> Clone( const Epetra_MultiVector& mv, const int numvecs )
    { return Teuchos::rcp( new Epetra_MultiVector(mv.Map(), numvecs) ); }
    ///
    static Teuchos::RefCountPtr<Epetra_MultiVector> CloneCopy( const Epetra_MultiVector& mv )
    { return Teuchos::rcp( new Epetra_MultiVector( mv ) ); }
    ///
    static Teuchos::RefCountPtr<Epetra_MultiVector> CloneCopy( const Epetra_MultiVector& mv, int index[], int numvecs )
    { return Teuchos::rcp( new Epetra_MultiVector(Copy, mv, index, numvecs) ); }
    ///
    static Teuchos::RefCountPtr<Epetra_MultiVector> CloneView( Epetra_MultiVector& mv, int index[], int numvecs )
    { return Teuchos::rcp( new Epetra_MultiVector(View, mv, index, numvecs) ); }
    ///
    static Teuchos::RefCountPtr<const Epetra_MultiVector> CloneView( const Epetra_MultiVector& mv, int index[], int numvecs )
    { return Teuchos::rcp( new Epetra_MultiVector(View, mv, index, numvecs) ); }
    ///
    static int GetVecLength( const Epetra_MultiVector& mv )
    { return mv.GlobalLength(); }
    ///
    static int GetNumberVecs( const Epetra_MultiVector& mv )
    { return mv.NumVectors(); }
    ///
    static void MvTimesMatAddMv( double alpha, const Epetra_MultiVector& A, 
				 const Teuchos::SerialDenseMatrix<int,double>& B, 
				 double beta, Epetra_MultiVector& mv )
    { 
      Epetra_LocalMap LocalMap(B.numRows(), 0, mv.Map().Comm());
      Epetra_MultiVector B_Pvec(Copy, LocalMap, B.values(), B.stride(), B.numCols());

      assert( mv.Multiply( 'N', 'N', alpha, A, B_Pvec, beta ) == 0 );   
    }
    ///
    static void MvAddMv( double alpha, const Epetra_MultiVector& A, double beta, const Epetra_MultiVector& B, Epetra_MultiVector& mv )
    { 
      assert( mv.Update( alpha, A, beta, B, 0.0 ) == 0 );
    }
    ///
    static void MvTransMv( double alpha, const Epetra_MultiVector& A, const Epetra_MultiVector& mv, Teuchos::SerialDenseMatrix<int,double>& B )
    { 
      Epetra_LocalMap LocalMap(B.numRows(), 0, mv.Map().Comm());
      Epetra_MultiVector B_Pvec(View, LocalMap, B.values(), B.stride(), B.numCols());
      
      assert( B_Pvec.Multiply( 'T', 'N', alpha, A, mv, 0.0 ) == 0 );
    }
    ///
    static void MvDot( const Epetra_MultiVector& mv, const Epetra_MultiVector& A, double b[] )
    {
      if (b) 
	assert( mv.Dot( A, b ) == 0 );
    }
    ///
    static void MvNorm( const Epetra_MultiVector& mv, double *normvec )
    { 
      if (normvec)
	assert( mv.Norm2(normvec) == 0 );
    }
    ///
    static void SetBlock( const Epetra_MultiVector& A, int index[], int numvecs, Epetra_MultiVector& mv )
    { 
      // Extract the "numvecs" columns of mv indicated by the index vector.
      Epetra_MultiVector temp_vec(View, mv, index, numvecs);

      if ( A.NumVectors() != numvecs ) {
        std::vector<int> index2( numvecs );
        for(int i=0; i<numvecs; i++)
	  index2[i] = i;
        Epetra_MultiVector A_vec(View, A, &index2[0], numvecs);      
        assert( temp_vec.Update( 1.0, A_vec, 0.0, A_vec, 0.0 ) == 0 );
      }
      else
        assert( temp_vec.Update( 1.0, A, 0.0, A, 0.0 ) == 0 );
    }
    ///
    static void MvRandom( Epetra_MultiVector& mv )
    { assert( mv.Random() == 0 ); }
    ///
    static void MvInit( Epetra_MultiVector& mv, double alpha = Teuchos::ScalarTraits<double>::zero() )
    { assert( mv.PutScalar(alpha) == 0); }
    ///
    static void MvPrint( const Epetra_MultiVector& mv, ostream& os )
    { os << mv << endl; }
    
  };        

  ////////////////////////////////////////////////////////////////////
  //
  // Implementation of the Anasazi::OperatorTraits for Epetra::Operator.
  //
  ////////////////////////////////////////////////////////////////////

  template <> 
  class OperatorTraits < double, Epetra_MultiVector, Epetra_Operator >
  {
  public:
    
    ///
    static ReturnType Apply ( const Epetra_Operator& Op, 
			      const Epetra_MultiVector& x, 
			      Epetra_MultiVector& y )
    { return ( Op.Apply( x, y ) == 0 ? Ok : Failed ); }
    
  };
  
} // end of Anasazi namespace 

#endif 
// end of file ANASAZI_EPETRA_ADAPTER_HPP
