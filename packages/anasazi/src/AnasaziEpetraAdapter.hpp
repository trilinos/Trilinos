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

/*! \file AnasaziEpetraAdapter.hpp
  \brief Implementations of Anasazi multi-vector and operator classes using Epetra_MultiVector and Epetra_Operator classes
*/

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
  
  /*! \class EpetraMultiVec
    \brief Basic adapter class for Anasazi::MultiVec that uses Epetra_MultiVector.

    \note The Epetra package performs double-precision arithmetic, so the use of Epetra with Anasazi will
    only provide a double-precision eigensolver.
  */
  class EpetraMultiVec : public MultiVec<double>, public Epetra_MultiVector {
  public:
    //@{ \name Constructors/Destructors

    //! Basic EpetraMultiVec constructor.
    /*! @param Map [in] An Epetra_LocalMap, Epetra_Map or Epetra_BlockMap.
      @param numvecs [in] Number of vectors in multi-vector.

      \returns Pointer to an EpetraMultiVec
    */
    EpetraMultiVec(const Epetra_BlockMap& Map, const int numvecs);

    //! Copy constructor.
    EpetraMultiVec(const Epetra_MultiVector & P_vec);
    
    //! Create multi-vector with values from two dimensional array.
    /*! @param Map [in] An Epetra_LocalMap, Epetra_Map or Epetra_BlockMap
      @param array [in] Pointer to an array of double precision numbers.  The first vector starts at \c array, the
      second at \c array+stride, and so on.  This array is copied.
      @param numvecs [in] Number of vectors in the multi-vector.
      @param stride [in] The stride between vectors in memory of \c array.

      \returns Pointer to an EpetraMultiVec
    */
    EpetraMultiVec(const Epetra_BlockMap& Map, double * array, const int numvecs, const int stride=0);

    //! Create multi-vector from list of vectors in an existing EpetraMultiVec.
    /*! @param CV [in] Enumerated type set to Copy or View.
      @param P_vec [in] An existing fully constructed Epetra_MultiVector.
      @param index [in] A integer vector containing the indices of the vectors to copy out of \c P_vec.

      \returns Pointer to an EpetraMultiVec
    */
    EpetraMultiVec(Epetra_DataAccess CV, const Epetra_MultiVector& P_vec, const std::vector<int>& index);

    //! Destructor
    virtual ~EpetraMultiVec() {};

    //@}

    //@{ \name Creation methods

    /*! \brief Creates a new empty EpetraMultiVec containing \c numvecs columns.
      
    \returns Pointer to an EpetraMultiVec
    */
    MultiVec<double> * Clone ( const int numvecs ) const;

    /*! \brief Creates a new EpetraMultiVec and copies contents of \c *this into
      the new vector (deep copy).
      
      \returns Pointer to an EpetraMultiVec
    */	
    MultiVec<double> * CloneCopy () const;

    /*! \brief Creates a new EpetraMultiVec and copies the selected contents of \c *this 
      into the new vector (deep copy).  
      
      The copied vectors from \c *this are indicated by the \c index.size() indices in \c index.
      
      \returns Pointer to an EpetraMultiVec
    */
    MultiVec<double> * CloneCopy ( const std::vector<int>& index ) const;
    
    /*! \brief Creates a new EpetraMultiVec that shares the selected contents of \c *this.
      
    The index of the \c numvecs vectors shallow copied from \c *this are indicated by the
    indices given in \c index.
    
    \returns Pointer to an EpetraMultiVec
    */
    MultiVec<double> * CloneView ( const std::vector<int>& index );

    //@}

    //@{ \name Attribute methods	

    //! Obtain the vector length of *this.
    int GetNumberVecs () const { return NumVectors(); }

    //! Obtain the number of vectors in *this.
    int GetVecLength () const { return GlobalLength(); }

    //@}

    //@{ \name Update methods
    /*! \brief Update \c *this with \f$\alpha AB + \beta (*this)\f$.
     */
    void MvTimesMatAddMv ( const double alpha, const MultiVec<double>& A, 
			   const Teuchos::SerialDenseMatrix<int,double>& B, const double beta );

    /*! \brief Replace \c *this with \f$\alpha A + \beta B\f$.
     */
    void MvAddMv ( const double alpha, const MultiVec<double>& A, const double beta,
		   const MultiVec<double>& B);

    /*! \brief Compute a dense matrix \c B through the matrix-matrix multiply \f$\alpha A^T(*this)\f$.
    */
    void MvTransMv ( const double alpha, const MultiVec<double>& A, Teuchos::SerialDenseMatrix<int,double>& B ) const;

    /*! \brief Compute a vector \c b where the components are the individual dot-products, i.e. \f$ b[i] = A[i]^T(this[i])\f$ where \c A[i] is the i-th column of \c A.
	*/
    void MvDot ( const MultiVec<double>& A, std::vector<double>* b ) const;

    //@}
    //@{ \name Norm method
    
    /*! \brief Compute the 2-norm of each individual vector of \c *this.  
      Upon return, \c normvec[i] holds the 2-norm of the \c i-th vector of \c *this
    */
    void MvNorm ( std::vector<double>* normvec ) const {
      if (normvec && (normvec->size() >= GetNumberVecs()) ) {
	int ret = Norm2(&(*normvec)[0]);
	assert( ret == 0 );
      }
    };
    //@}

    //@{ \name Initialization methods
    /*! \brief Copy the vectors in \c A to a set of vectors in \c *this.  

    The \c numvecs vectors in \c A are copied to a subset of vectors in \c *this
    indicated by the indices given in \c index.
    */
    void SetBlock ( const MultiVec<double>& A, const std::vector<int>& index );

    /*! \brief Fill the vectors in \c *this with random numbers.
     */
    void MvRandom() { int ret = Random(); assert( ret == 0 ); };

    /*! \brief Replace each element of the vectors in \c *this with \c alpha.
     */
    void MvInit ( const double alpha ) { int ret = PutScalar( alpha ); assert( ret == 0 ); };

    //@}
    //@{ \name Print method.
    /*! \brief Print \c *this EpetraMultiVec.
     */
    void MvPrint() const { std::cout<< *this << endl; };
    //@}

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
				   Teuchos::SerialDenseMatrix<int,double>& B) const
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
  
  void EpetraMultiVec::MvDot ( const MultiVec<double>& A, std::vector<double>* b ) const
  {
    EpetraMultiVec *A_vec = dynamic_cast<EpetraMultiVec *>(&const_cast<MultiVec<double> &>(A)); 
    assert(A_vec!=NULL);
    if (A_vec && b && ( b->size() >= A_vec->NumVectors() ) ) {
      int ret = this->Dot( *A_vec, &(*b)[0] );
      assert( ret == 0 );
    }
  }
  
  ///////////////////////////////////////////////////////////////
  //--------template class AnasaziEpetraOp---------------------
  
  /*! \class EpetraOp
    \brief Basic adapter class for Anasazi::Operator that uses Epetra_Operator.

    \note The Epetra package performs double-precision arithmetic, so the use of Epetra with Anasazi will
    only provide a double-precision eigensolver.
  */
  class EpetraOp : public virtual Operator<double> {
  public:
    //@{ \name Constructor/Destructor
    
    //! Basic constructor.  Accepts reference-counted pointer to an Epetra_Operator.
    EpetraOp(const Teuchos::RefCountPtr<Epetra_Operator> &Op );
    
    //! Destructor
    ~EpetraOp();
    //@}
    
    //@{ \name Operator application method.
    
    /*! \brief This method takes the Anasazi::MultiVec \c x and
      applies the operator to it resulting in the Anasazi::MultiVec \c y.
    */
    ReturnType Apply ( const MultiVec<double>& x, 
		       MultiVec<double>& y ) const;
    //@} 
    
  private:
    Teuchos::RefCountPtr<Epetra_Operator> Epetra_Op;
  };
  //-------------------------------------------------------------
  //
  // implementation of the AnasaziEpetraOp class.
  //
  ////////////////////////////////////////////////////////////////////
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

    int info = Epetra_Op->Apply( *vec_x, *vec_y );
    
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
    EpetraGenOp(const Teuchos::RefCountPtr<Epetra_Operator> &AOp, 
                const Teuchos::RefCountPtr<Epetra_Operator> &MOp );
    ~EpetraGenOp();
    ReturnType Apply ( const MultiVec<double>& x, 
		       MultiVec<double>& y ) const;
  private:
    Teuchos::RefCountPtr<Epetra_Operator> Epetra_AOp;
    Teuchos::RefCountPtr<Epetra_Operator> Epetra_MOp;
  };
  //-------------------------------------------------------------
  //
  // implementation of the AnasaziEpetraGenOp class.
  //
  ////////////////////////////////////////////////////////////////////
  //
  // AnasaziOperator constructors
  //
  
  EpetraGenOp::EpetraGenOp(const Teuchos::RefCountPtr<Epetra_Operator> &AOp,
			   const Teuchos::RefCountPtr<Epetra_Operator> &MOp) 
    : Epetra_AOp(AOp), Epetra_MOp(MOp) 
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
    // This generalized operator computes y = A*M*x of y = (A*M)^T*x
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
    // Apply M
    info = Epetra_MOp->Apply( *vec_x, temp_y );
    assert(info==0);
    // Apply A
    info = Epetra_AOp->Apply( temp_y, *vec_y );
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
    EpetraSymOp(const Teuchos::RefCountPtr<Epetra_Operator> &Op, const bool isTrans = false );
    ~EpetraSymOp();
    ReturnType Apply ( const MultiVec<double>& x, 
		       MultiVec<double>& y ) const;
  private:
    Teuchos::RefCountPtr<Epetra_Operator> Epetra_Op;
    bool isTrans_;
  };
  //-------------------------------------------------------------
  //
  // implementation of the AnasaziEpetraSymOp class.
  //
  ////////////////////////////////////////////////////////////////////
  //
  // AnasaziOperator constructors
  //
  EpetraSymOp::EpetraSymOp(const Teuchos::RefCountPtr<Epetra_Operator> &Op, const bool isTrans) 
    : Epetra_Op(Op), isTrans_(isTrans)
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
    Epetra_MultiVector* temp_vec = new Epetra_MultiVector( 
							  (isTrans_) ? Epetra_Op->OperatorDomainMap() 
							  : Epetra_Op->OperatorRangeMap(), 
							  vec_x->NumVectors() );
    
    assert( vec_x!=NULL && vec_y!=NULL && temp_vec!=NULL );
    //
    // Need to cast away constness because the member function Apply
    // is not declared const.
    //
    // Transpose the operator (if isTrans_ = true)
    if (isTrans_) {
      info=Epetra_Op->SetUseTranspose( isTrans_ );
      if (info!=0) { delete temp_vec; return Failed; }
    }
    //
    // Compute A*x or A'*x 
    //
    info=Epetra_Op->Apply( *vec_x, *temp_vec );
    if (info!=0) { delete temp_vec; return Failed; }
    //
    // Transpose/Un-transpose the operator based on value of isTrans_
    info=Epetra_Op->SetUseTranspose( !isTrans_ );
    if (info!=0) { delete temp_vec; return Failed; }
    
    // Compute A^T*(A*x) or A*A^T
    info=Epetra_Op->Apply( *temp_vec, *vec_y );
    if (info!=0) { delete temp_vec; return Failed; }
    
    // Un-transpose the operator
    info=Epetra_Op->SetUseTranspose( false );
    delete temp_vec;
    
    if (info==0)
      return Ok; 
    else
      return Failed; 
  }
  

  ///////////////////////////////////////////////////////////////
  //--------template class AnasaziEpetraSymMVOp---------------------
  class EpetraSymMVOp : public virtual Operator<double> {
  public:
    EpetraSymMVOp(const Teuchos::RefCountPtr<Epetra_MultiVector> &MV, const bool isTrans = false );
    ~EpetraSymMVOp();
    ReturnType Apply ( const MultiVec<double>& x, 
		       MultiVec<double>& y ) const;
  private:
    Teuchos::RefCountPtr<Epetra_MultiVector> Epetra_MV;
    bool isTrans_;
  };
  //-------------------------------------------------------------
  //
  // implementation of the Anasazi::EpetraSymMVOp class.
  //
  ////////////////////////////////////////////////////////////////////
  //
  // Anasazi::Operator constructors
  //
  EpetraSymMVOp::EpetraSymMVOp(const Teuchos::RefCountPtr<Epetra_MultiVector> &MV, const bool isTrans) 
    : Epetra_MV(MV), isTrans_(isTrans)
  {
  }
  
  EpetraSymMVOp::~EpetraSymMVOp() 
  {
  }
  //
  // AnasaziOperator applications
  //
  ReturnType EpetraSymMVOp::Apply ( const MultiVec<double>& x, 
				    MultiVec<double>& y ) const 
  {
    int info=0;
    MultiVec<double> & temp_x = const_cast<MultiVec<double> &>(x);
    Epetra_MultiVector* vec_x = dynamic_cast<Epetra_MultiVector* >(&temp_x);
    Epetra_MultiVector* vec_y = dynamic_cast<Epetra_MultiVector* >(&y);
    
    if (isTrans_) {
      Epetra_LocalMap localMap( Epetra_MV->NumVectors(), 0, Epetra_MV->Map().Comm() );
      Epetra_MultiVector Pvec( localMap, temp_x.GetNumberVecs() );
      
      /* A'*x */
      info = Pvec.Multiply( 'T', 'N', 1.0, *Epetra_MV, *vec_x, 0.0 );
      
      /* A*(A'*x) */
      info = vec_y->Multiply( 'N', 'N', 1.0, *Epetra_MV, Pvec, 0.0 );
      
    } 
    else {
      Epetra_MultiVector temp( Epetra_MV->Map(), temp_x.GetNumberVecs() );
      
      /* A*x */
      info = temp.Multiply( 'N', 'N', 1.0, *Epetra_MV, *vec_x, 0.0 );
      
      /* A'*(A*x) */
      info = vec_y->Multiply( 'T', 'N', 1.0, *Epetra_MV, temp, 0.0 );
    }
    
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

  /*! \class MultiVecTraits< double, Epetra_MultiVector >
    \brief Template specialization of Anasazi::MultiVecTraits class using the Epetra_MultiVector class.

    This interface will ensure that any Epetra_MultiVector will be accepted by the Anasazi
    templated solvers.  

    \note The Epetra package performs double-precision arithmetic, so the use of Epetra with Anasazi will
    only provide a double-precision eigensolver.
  */

  template<>
  class MultiVecTraits<double, Epetra_MultiVector>
  {
  public:

    //@{ \name Creation methods

    /*! \brief Creates a new empty Epetra_MultiVector containing \c numvecs columns.
      
    \return Reference-counted pointer to the new Epetra_MultiVector.
    */
    static Teuchos::RefCountPtr<Epetra_MultiVector> Clone( const Epetra_MultiVector& mv, const int numvecs )
    { return Teuchos::rcp( new Epetra_MultiVector(mv.Map(), numvecs) ); }

    /*! \brief Creates a new Epetra_MultiVector and copies contents of \c mv into the new vector (deep copy).
      
      \return Reference-counted pointer to the new Epetra_MultiVector.
    */
    static Teuchos::RefCountPtr<Epetra_MultiVector> CloneCopy( const Epetra_MultiVector& mv )
    { return Teuchos::rcp( new Epetra_MultiVector( mv ) ); }

    /*! \brief Creates a new Epetra_MultiVector and copies the selected contents of \c mv into the new vector (deep copy).  

      The copied vectors from \c mv are indicated by the \c index.size() indices in \c index.      
      \return Reference-counted pointer to the new Epetra_MultiVector.
    */
    static Teuchos::RefCountPtr<Epetra_MultiVector> CloneCopy( const Epetra_MultiVector& mv, const std::vector<int>& index )
    { 
      std::vector<int>& tmp_index = const_cast<std::vector<int> &>( index );
      return Teuchos::rcp( new Epetra_MultiVector(Copy, mv, &tmp_index[0], index.size()) ); 
    }

    /*! \brief Creates a new Epetra_MultiVector that shares the selected contents of \c mv (shallow copy).

    The index of the \c numvecs vectors shallow copied from \c mv are indicated by the indices given in \c index.
    \return Reference-counted pointer to the new Epetra_MultiVector.
    */      
    static Teuchos::RefCountPtr<Epetra_MultiVector> CloneView( Epetra_MultiVector& mv, const std::vector<int>& index )
    { 
      std::vector<int>& tmp_index = const_cast<std::vector<int> &>( index );
      return Teuchos::rcp( new Epetra_MultiVector(View, mv, &tmp_index[0], index.size()) ); 
    }

    /*! \brief Creates a new const Epetra_MultiVector that shares the selected contents of \c mv (shallow copy).

    The index of the \c numvecs vectors shallow copied from \c mv are indicated by the indices given in \c index.
    \return Reference-counted pointer to the new const Epetra_MultiVector.
    */      
    static Teuchos::RefCountPtr<const Epetra_MultiVector> CloneView( const Epetra_MultiVector& mv, const std::vector<int>& index )
    { 
      std::vector<int>& tmp_index = const_cast<std::vector<int> &>( index );
      return Teuchos::rcp( new Epetra_MultiVector(View, mv, &tmp_index[0], index.size()) ); 
    }

    //@}

    //@{ \name Attribute methods

    //! Obtain the vector length of \c mv.
    static int GetVecLength( const Epetra_MultiVector& mv )
    { return mv.GlobalLength(); }

    //! Obtain the number of vectors in \c mv
    static int GetNumberVecs( const Epetra_MultiVector& mv )
    { return mv.NumVectors(); }
    //@}

    //@{ \name Update methods

    /*! \brief Update \c mv with \f$ \alpha AB + \beta mv \f$.
     */
    static void MvTimesMatAddMv( const double alpha, const Epetra_MultiVector& A, 
				 const Teuchos::SerialDenseMatrix<int,double>& B, 
				 const double beta, Epetra_MultiVector& mv )
    { 
      Epetra_LocalMap LocalMap(B.numRows(), 0, mv.Map().Comm());
      Epetra_MultiVector B_Pvec(Copy, LocalMap, B.values(), B.stride(), B.numCols());

      int ret = mv.Multiply( 'N', 'N', alpha, A, B_Pvec, beta );
      assert( ret == 0 );   
    }

    /*! \brief Replace \c mv with \f$\alpha A + \beta B\f$.
     */
    static void MvAddMv( const double alpha, const Epetra_MultiVector& A, const double beta, const Epetra_MultiVector& B, Epetra_MultiVector& mv )
    { 
      int ret = mv.Update( alpha, A, beta, B, 0.0 );
      assert( ret == 0 );
    }

    /*! \brief Compute a dense matrix \c B through the matrix-matrix multiply \f$ \alpha A^Tmv \f$.
    */
    static void MvTransMv( const double alpha, const Epetra_MultiVector& A, const Epetra_MultiVector& mv, Teuchos::SerialDenseMatrix<int,double>& B )
    { 
      Epetra_LocalMap LocalMap(B.numRows(), 0, mv.Map().Comm());
      Epetra_MultiVector B_Pvec(View, LocalMap, B.values(), B.stride(), B.numCols());
      
      int ret = B_Pvec.Multiply( 'T', 'N', alpha, A, mv, 0.0 );
      assert( ret == 0 );
    }
    
    /*! \brief Compute a vector \c b where the components are the individual dot-products of the \c i-th columns of \c A and \c mv, i.e.\f$b[i] = A[i]^Tmv[i]\f$.
     */
    static void MvDot( const Epetra_MultiVector& mv, const Epetra_MultiVector& A, std::vector<double>* b )
    {
      int ret = mv.Dot( A, &(*b)[0] );
      assert( ret == 0 );
    }

    //@}
    //@{ \name Norm method

    /*! \brief Compute the 2-norm of each individual vector of \c mv.  
      Upon return, \c normvec[i] holds the value of \f$||mv_i||_2\f$, the \c i-th column of \c mv.
    */
    static void MvNorm( const Epetra_MultiVector& mv, std::vector<double>* normvec )
    { 
      int ret = mv.Norm2(&(*normvec)[0]);
      assert( ret == 0 );
    }

    //@}

    //@{ \name Initialization methods
    /*! \brief Copy the vectors in \c A to a set of vectors in \c mv indicated by the indices given in \c index.
     */
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
        int ret = temp_vec.Update( 1.0, A_vec, 0.0, A_vec, 0.0 );
	assert( ret == 0 );
      }
      else {
        int ret = temp_vec.Update( 1.0, A, 0.0, A, 0.0 );
	assert( ret == 0 );
      }
    }

    /*! \brief Replace the vectors in \c mv with random vectors.
     */
    static void MvRandom( Epetra_MultiVector& mv )
    { int ret = mv.Random(); assert( ret == 0 ); }

    /*! \brief Replace each element of the vectors in \c mv with \c alpha.
     */
    static void MvInit( Epetra_MultiVector& mv, double alpha = Teuchos::ScalarTraits<double>::zero() )
    { int ret = mv.PutScalar(alpha); assert( ret == 0 ); }

    //@}

    //@{ \name Print method

    /*! \brief Print the \c mv multi-vector to the \c os output stream.
     */
    static void MvPrint( const Epetra_MultiVector& mv, ostream& os )
    { os << mv << endl; }

    //@}
  };        

  ////////////////////////////////////////////////////////////////////
  //
  // Implementation of the Anasazi::OperatorTraits for Epetra::Operator.
  //
  ////////////////////////////////////////////////////////////////////

  /*! \class OperatorTraits< double, Epetra_MultiVector, Epetra_Operator >
    \brief Template specialization of Anasazi::OperatorTraits class using the Epetra_Operator virtual base class and 
    Epetra_MultiVector class.

    This interface will ensure that any Epetra_Operator and Epetra_MultiVector will be accepted by the Anasazi
    templated solvers.

    \note The Epetra package performs double-precision arithmetic, so the use of Epetra with Anasazi will
    only provide a double-precision eigensolver.
  */

  template <> 
  class OperatorTraits < double, Epetra_MultiVector, Epetra_Operator >
  {
  public:
    
    /*! \brief This method takes the Epetra_MultiVector \c x and
      applies the Epetra_Operator \c Op to it resulting in the Epetra_MultiVector \c y.
    */    
    static ReturnType Apply ( const Epetra_Operator& Op, 
			      const Epetra_MultiVector& x, 
			      Epetra_MultiVector& y )
    { return ( Op.Apply( x, y ) == 0 ? Ok : Failed ); }
    
  };
  
} // end of Anasazi namespace 

#endif 
// end of file ANASAZI_EPETRA_ADAPTER_HPP
