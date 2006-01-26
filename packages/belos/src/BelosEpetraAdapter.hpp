// @HEADER
// ***********************************************************************
//
//                 Belos: Block Linear Solvers Package
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

#ifndef BELOS_EPETRA_ADAPTER_HPP
#define BELOS_EPETRA_ADAPTER_HPP

/*! \file BelosEpetraAdapter.hpp
    \brief Provides several interfaces between Belos virtual classes and Epetra concrete classes.
*/

#include "Epetra_MultiVector.h"
#include "Epetra_Operator.h"
#include "Epetra_Map.h"
#include "Epetra_LocalMap.h"
#include "Teuchos_SerialDenseMatrix.hpp"

#include "BelosConfigDefs.hpp"
#include "BelosMultiVec.hpp"
#include "BelosOperator.hpp"

namespace Belos {
  
  //--------template class BelosEpetraMultiVec-------------------------------------
  class EpetraMultiVec : public MultiVec<double>, public Epetra_MultiVector {
  public:
    // constructors
    EpetraMultiVec(const Epetra_BlockMap& Map, double * array, const int numvecs, const int stride=0);
    EpetraMultiVec(const Epetra_BlockMap& Map, const int numvecs);
    EpetraMultiVec(Epetra_DataAccess CV, const Epetra_MultiVector& P_vec, const std::vector<int>& index);
    EpetraMultiVec& operator=(const EpetraMultiVec& pv) { Epetra_MultiVector::operator=(pv); return *this; }
    EpetraMultiVec(const Epetra_MultiVector & P_vec);
    ~EpetraMultiVec();
    //
    //  member functions inherited from Belos::MultiVec
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
    MultiVec<double> * CloneCopy ( const std::vector<int>& index ) const;
    //
    //  the following is a virtual view constructor returning
    //  a pointer to the pure virtual class. vector values are 
    //  shared and hence no memory is allocated for the columns.
    //
    MultiVec<double> * CloneView ( const std::vector<int>& index );
    //
    //  this routine sets a subblock of the multivector, which
    //  need not be contiguous, and is given by the indices.
    //
    void SetBlock ( const MultiVec<double>& A, const std::vector<int>& index );
    //
    int GetNumberVecs () const { return NumVectors(); }
    int GetVecLength () const { return GlobalLength(); }
    //
    // *this <- alpha * A * B + beta * (*this)
    //
    void MvTimesMatAddMv ( const double alpha, const MultiVec<double>& A, 
			   const Teuchos::SerialDenseMatrix<int,double>& B, const double beta );
    //
    // *this <- alpha * A + beta * B
    //
    void MvAddMv ( const double alpha, const MultiVec<double>& A, const double beta,
		   const MultiVec<double>& B);
    //
    // B <- alpha * A^T * (*this)
    //
    void MvTransMv ( const double alpha, const MultiVec<double>& A, Teuchos::SerialDenseMatrix<int,double>& B ) const;
    //
    // b[i] = A[i]^T * this[i]
    // 
    void MvDot ( const MultiVec<double>& A, std::vector<double>* b ) const;
    //
    // alpha[i] = norm of i-th column of (*this)
    //	
    void MvNorm ( std::vector<double>* normvec, NormType norm_type = TwoNorm ) const;
    //
    // random vectors in i-th column of (*this)
    //
    void MvRandom() { assert( Random() == 0 ); };
    //
    // initializes each element of (*this) with alpha
    //
    void MvInit ( const double alpha ) { assert( PutScalar( alpha ) == 0 ); };
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
				 const std::vector<int>& index )
    : Epetra_MultiVector(CV, P_vec, &(const_cast<std::vector<int> &>(index))[0], index.size())
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
  //  member functions inherited from Belos::MultiVec
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
      EpetraMultiVec *tmp_vec = dynamic_cast<EpetraMultiVec *>(&const_cast<MultiVec<double> &>(A)); assert(tmp_vec!=NULL);
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
    
    EpetraMultiVec *A_vec = dynamic_cast<EpetraMultiVec *>(&const_cast<MultiVec<double> &>(A)); assert(A_vec!=NULL);
    
    assert( Multiply( 'N', 'N', alpha, *A_vec, B_Pvec, beta ) == 0 );
  }

  //-------------------------------------------------------------
  //
  // *this <- alpha * A + beta * B
  //
  //-------------------------------------------------------------
  
  void EpetraMultiVec::MvAddMv ( const double alpha , const MultiVec<double>& A, 
				 const double beta, const MultiVec<double>& B) 
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
  
  void EpetraMultiVec::MvTransMv ( const double alpha, const MultiVec<double>& A,
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
  
  void EpetraMultiVec::MvDot ( const MultiVec<double>& A, std::vector<double>* b ) const
  {
    EpetraMultiVec *A_vec = dynamic_cast<EpetraMultiVec *>(&const_cast<MultiVec<double> &>(A)); assert(A_vec!=NULL);
    if (A_vec && b && ( (int)b->size() >= A_vec->NumVectors() ) ) {
      assert( this->Dot( *A_vec, &(*b)[0] ) == 0 );
    }
  }

  //-------------------------------------------------------------
  //
  // alpha[i] = norm of i-th column of (*this)
  //
  //-------------------------------------------------------------
  
  void EpetraMultiVec::MvNorm ( std::vector<double>* normvec, NormType norm_type ) const {
    if (normvec && ((int)normvec->size() >= GetNumberVecs())) {
      switch( norm_type ) {
      case ( OneNorm ) :
	assert( Norm1(&(*normvec)[0]) == 0 );
	break;
      case ( TwoNorm ) :
	assert( Norm2(&(*normvec)[0]) == 0 );
	break;
      case ( InfNorm ) :	
	assert( NormInf(&(*normvec)[0]) == 0 );
	break;
      default:
	break;
      }
    }
  }
  
  ///////////////////////////////////////////////////////////////
  //--------template class BelosEpetraOp---------------------
  
  class EpetraOp : public virtual Operator<double> {
  public:
    EpetraOp( const Teuchos::RefCountPtr<Epetra_Operator> &Op );
    ~EpetraOp() {};
    ReturnType Apply ( const MultiVec<double>& x, MultiVec<double>& y, ETrans trans=NOTRANS ) const;
  private:
    Teuchos::RefCountPtr<Epetra_Operator> Epetra_Op;
  };
  //-------------------------------------------------------------
  //
  // implementation of the BelosEpetraOp class.
  //
  ////////////////////////////////////////////////////////////////////
  //
  // BelosOperator constructors
  //
  
  EpetraOp::EpetraOp( const Teuchos::RefCountPtr<Epetra_Operator> &Op ) 
    : Epetra_Op(Op)
  {
  }
  
  //
  // BelosOperator applications
  //
  ReturnType EpetraOp::Apply ( const MultiVec<double>& x, 
			       MultiVec<double>& y, ETrans trans ) const {
    int info=0;
    MultiVec<double> & temp_x = const_cast<MultiVec<double> &>(x);
    Epetra_MultiVector* vec_x = dynamic_cast<Epetra_MultiVector* >(&temp_x);
    Epetra_MultiVector* vec_y = dynamic_cast<Epetra_MultiVector* >(&y);
    //
    assert( vec_x!=NULL && vec_y!=NULL );
    //
    // Set the operator to apply the transpose if that is requested.
    // (TO DO:  insert a throw if the application returns a nonzero )
    //
    if ( trans ) { 
      info = Epetra_Op->SetUseTranspose( true );
      if (info != 0) { return Undefined; }
    }
    info = Epetra_Op->Apply( *vec_x, *vec_y );
    
    if ( trans ) { 
      info = Epetra_Op->SetUseTranspose( false );
      if (info != 0) { return Undefined; }
    }
    
    if (info != 0) { return Error; }
    return Ok;	
  }
  
  ///////////////////////////////////////////////////////////////
  //--------template class BelosEpetraPrecOp---------------------
  
  class EpetraPrecOp : public virtual Operator<double> {
  public:
    EpetraPrecOp( const Teuchos::RefCountPtr<Epetra_Operator> &Op );
    ~EpetraPrecOp() {};
    ReturnType Apply ( const MultiVec<double>& x, MultiVec<double>& y, ETrans trans=NOTRANS ) const;
  private:
    Teuchos::RefCountPtr<Epetra_Operator> Epetra_Op;
  };
  //-------------------------------------------------------------
  //
  // implementation of the BelosEpetraPrecOp class.
  //
  ////////////////////////////////////////////////////////////////////
  //
  // BelosOperator constructors
  //
  
  EpetraPrecOp::EpetraPrecOp( const Teuchos::RefCountPtr<Epetra_Operator> &Op ) 
    : Epetra_Op(Op)
  {
  }
  
  //
  // BelosOperator applications
  //
  ReturnType EpetraPrecOp::Apply ( const MultiVec<double>& x, 
			           MultiVec<double>& y, ETrans trans ) const {
    int info=0;
    MultiVec<double> & temp_x = const_cast<MultiVec<double> &>(x);
    Epetra_MultiVector* vec_x = dynamic_cast<Epetra_MultiVector* >(&temp_x);
    Epetra_MultiVector* vec_y = dynamic_cast<Epetra_MultiVector* >(&y);
    //
    assert( vec_x!=NULL && vec_y!=NULL );
    //
    // Set the operator to apply the transpose if that is requested.
    // (TO DO:  insert a throw if the application returns a nonzero )
    //
    if ( trans ) { 
      info = Epetra_Op->SetUseTranspose( true );
      if (info != 0) { return Undefined; }
    }
    info = Epetra_Op->ApplyInverse( *vec_x, *vec_y );
    
    if ( trans ) { 
      info = Epetra_Op->SetUseTranspose( false );
      if (info != 0) { return Undefined; }
    }
    
    if (info != 0) { return Error; }
    return Ok;	
  }
  
  ////////////////////////////////////////////////////////////////////
  //
  // Implementation of the Belos::MultiVecTraits for Epetra::MultiVector.
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
    static Teuchos::RefCountPtr<Epetra_MultiVector> CloneCopy( const Epetra_MultiVector& mv, const std::vector<int>& index )
    { 
      std::vector<int>& tmp_index = const_cast<std::vector<int> &>( index );
      return Teuchos::rcp( new Epetra_MultiVector(Copy, mv, &tmp_index[0], index.size()) ); 
    }
    ///
    static Teuchos::RefCountPtr<Epetra_MultiVector> CloneView( Epetra_MultiVector& mv, const std::vector<int>& index )
    { 
      std::vector<int>& tmp_index = const_cast<std::vector<int> &>( index );
      return Teuchos::rcp( new Epetra_MultiVector(View, mv, &tmp_index[0], index.size()) ); 
    }
    ///
    static Teuchos::RefCountPtr<const Epetra_MultiVector> CloneView( const Epetra_MultiVector& mv, const std::vector<int>& index )
    { 
      std::vector<int>& tmp_index = const_cast<std::vector<int> &>( index );
      return Teuchos::rcp( new Epetra_MultiVector(View, mv, &tmp_index[0], index.size()) ); 
    }
    ///
    static int GetVecLength( const Epetra_MultiVector& mv )
    { return mv.GlobalLength(); }
    ///
    static int GetNumberVecs( const Epetra_MultiVector& mv )
    { return mv.NumVectors(); }
    ///
    static void MvTimesMatAddMv( const double alpha, const Epetra_MultiVector& A, 
				 const Teuchos::SerialDenseMatrix<int,double>& B, 
				 const double beta, Epetra_MultiVector& mv )
    { 
      Epetra_LocalMap LocalMap(B.numRows(), 0, mv.Map().Comm());
      Epetra_MultiVector B_Pvec(Copy, LocalMap, B.values(), B.stride(), B.numCols());
      
      assert( mv.Multiply( 'N', 'N', alpha, A, B_Pvec, beta ) == 0 );   
    }
    ///
    static void MvAddMv( const double alpha, const Epetra_MultiVector& A, const double beta, const Epetra_MultiVector& B, Epetra_MultiVector& mv )
    { 
      assert( mv.Update( alpha, A, beta, B, 0.0 ) == 0 );
    }
    ///
    static void MvTransMv( const double alpha, const Epetra_MultiVector& A, const Epetra_MultiVector& mv, Teuchos::SerialDenseMatrix<int,double>& B )
    { 
      Epetra_LocalMap LocalMap(B.numRows(), 0, mv.Map().Comm());
      Epetra_MultiVector B_Pvec(View, LocalMap, B.values(), B.stride(), B.numCols());
      
      assert( B_Pvec.Multiply( 'T', 'N', alpha, A, mv, 0.0 ) == 0 );
    }
    ///
    static void MvDot( const Epetra_MultiVector& mv, const Epetra_MultiVector& A, std::vector<double>* b )
    {
      assert( mv.Dot( A, &(*b)[0] ) == 0 );
    }
    ///
    static void MvNorm( const Epetra_MultiVector& mv, std::vector<double>* normvec, NormType type = TwoNorm )
    { 
      if (normvec && ((int)normvec->size() >= mv.NumVectors())) {
	switch( type ) {
	case ( OneNorm ) :
	  assert( mv.Norm1(&(*normvec)[0]) == 0 );
	  break;
	case ( TwoNorm ) :
	  assert( mv.Norm2(&(*normvec)[0]) == 0 );
	  break;
	case ( InfNorm ) :	
	  assert( mv.NormInf(&(*normvec)[0]) == 0 );
	  break;
	default:
	  break;
	}
      }
    }
    ///
    static void SetBlock( const Epetra_MultiVector& A, const std::vector<int>& index, Epetra_MultiVector& mv )
    { 
      // Extract the "numvecs" columns of mv indicated by the index vector.
      int numvecs = index.size();
      std::vector<int>& tmp_index = const_cast<std::vector<int> &>( index );
      Epetra_MultiVector temp_vec(View, mv, &tmp_index[0], numvecs);
      
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
  // Implementation of the Belos::OperatorTraits for Epetra::Operator.
  //
  ////////////////////////////////////////////////////////////////////
  
  template <> 
  class OperatorTraits < double, Epetra_MultiVector, Epetra_Operator >
  {
  public:
    
    ///
    static ReturnType Apply ( const Epetra_Operator& Op, 
			      const Epetra_MultiVector& x, 
			      Epetra_MultiVector& y,
			      ETrans trans=NOTRANS )
    { 
      int info = 0;
      if ( trans )
	const_cast<Epetra_Operator &>(Op).SetUseTranspose( true );
      info = Op.Apply( x, y );
      if ( trans )
	const_cast<Epetra_Operator &>(Op).SetUseTranspose( false );      
      
      return ( info == 0 ? Ok : Error ); 
    }
    
  };
  
} // end of Belos namespace 

#endif 
// end of file BELOS_EPETRA_ADAPTER_HPP
