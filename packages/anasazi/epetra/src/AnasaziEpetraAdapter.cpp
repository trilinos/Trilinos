// @HEADER
// *****************************************************************************
//                 Anasazi: Block Eigensolvers Package
//
// Copyright 2004 NTESS and the Anasazi contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
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
  
  EpetraMultiVec::EpetraMultiVec(const Epetra_BlockMap& Map_in, double * array, 
                                 const int numvecs, const int stride)
    : Epetra_MultiVector(Epetra_DataAccess::Copy, Map_in, array, stride, numvecs) 
  {
  }
  
  
  EpetraMultiVec::EpetraMultiVec(const Epetra_BlockMap& Map_in, const int numvecs)
    : Epetra_MultiVector(Map_in, numvecs) 
  {
  }
  
  
  EpetraMultiVec::EpetraMultiVec(Epetra_DataAccess CV, 
                                 const Epetra_MultiVector& P_vec, 
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
    EpetraMultiVec * ptr_apv = new EpetraMultiVec(Epetra_DataAccess::Copy, *this, index);
    return ptr_apv; // safe upcast.
  }
  
  
  MultiVec<double>* EpetraMultiVec::CloneViewNonConst ( const std::vector<int>& index ) 
  {
    EpetraMultiVec * ptr_apv = new EpetraMultiVec(Epetra_DataAccess::View, *this, index);
    return ptr_apv; // safe upcast.
  }
  
  const MultiVec<double>* EpetraMultiVec::CloneView ( const std::vector<int>& index ) const
  {
    EpetraMultiVec * ptr_apv = new EpetraMultiVec(Epetra_DataAccess::View, *this, index);
    return ptr_apv; // safe upcast.
  }
  
  
  void EpetraMultiVec::SetBlock( const MultiVec<double>& A, const std::vector<int>& index ) 
  {
    // this should be revisited to e
    EpetraMultiVec temp_vec(Epetra_DataAccess::View, *this, index);

    int numvecs = index.size();
    if ( A.GetNumberVecs() != numvecs ) {
      std::vector<int> index2( numvecs );
      for(int i=0; i<numvecs; i++)
        index2[i] = i;
      EpetraMultiVec *tmp_vec = dynamic_cast<EpetraMultiVec *>(&const_cast<MultiVec<double> &>(A)); 
      TEUCHOS_TEST_FOR_EXCEPTION( tmp_vec==NULL, std::invalid_argument, "Anasazi::EpetraMultiVec::SetBlocks() cast of MultiVec<double> to EpetraMultiVec failed.");
      EpetraMultiVec A_vec(Epetra_DataAccess::View, *tmp_vec, index2);
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
    Epetra_MultiVector B_Pvec(Epetra_DataAccess::View, LocalMap, B.values(), B.stride(), B.numCols());
    
    EpetraMultiVec *A_vec = dynamic_cast<EpetraMultiVec *>(&const_cast<MultiVec<double> &>(A)); 
    TEUCHOS_TEST_FOR_EXCEPTION( A_vec==NULL,  std::invalid_argument, "Anasazi::EpetraMultiVec::SetBlocks() cast of MultiVec<double> to EpetraMultiVec failed.");
    
    TEUCHOS_TEST_FOR_EXCEPTION( 
        Multiply( 'N', 'N', alpha, *A_vec, B_Pvec, beta ) != 0,
        EpetraMultiVecFailure, "Anasazi::EpetraMultiVec::MvTimesMatAddMv() call to Epetra_MultiVec::Multiply() returned a nonzero value.");
  }

  //-------------------------------------------------------------
  //
  // *this <- alpha * A + beta * B
  //
  //-------------------------------------------------------------
  
  void EpetraMultiVec::MvAddMv ( double alpha , const MultiVec<double>& A, 
                                 double beta, const MultiVec<double>& B) 
  {
    EpetraMultiVec *A_vec = dynamic_cast<EpetraMultiVec *>(&const_cast<MultiVec<double> &>(A)); 
    TEUCHOS_TEST_FOR_EXCEPTION( A_vec==NULL,  std::invalid_argument, "Anasazi::EpetraMultiVec::MvAddMv() cast of MultiVec<double> to EpetraMultiVec failed.");
    EpetraMultiVec *B_vec = dynamic_cast<EpetraMultiVec *>(&const_cast<MultiVec<double> &>(B)); 
    TEUCHOS_TEST_FOR_EXCEPTION( B_vec==NULL,  std::invalid_argument, "Anasazi::EpetraMultiVec::MvAddMv() cast of MultiVec<double> to EpetraMultiVec failed.");
    
    TEUCHOS_TEST_FOR_EXCEPTION( 
        Update( alpha, *A_vec, beta, *B_vec, 0.0 ) != 0,
        EpetraMultiVecFailure, "Anasazi::EpetraMultiVec::MvAddMv() call to Epetra_MultiVec::Update() returned a nonzero value.");
  }

  //-------------------------------------------------------------
  //
  // dense B <- alpha * A^T * (*this)
  //
  //-------------------------------------------------------------
  
  void EpetraMultiVec::MvTransMv ( double alpha, const MultiVec<double>& A,
                                   Teuchos::SerialDenseMatrix<int,double>& B
#ifdef HAVE_ANASAZI_EXPERIMENTAL
                                   , ConjType conj
#endif
                                  ) const
  {    
    EpetraMultiVec *A_vec = dynamic_cast<EpetraMultiVec *>(&const_cast<MultiVec<double> &>(A));
    
    if (A_vec) {
      Epetra_LocalMap LocalMap(B.numRows(), 0, Map().Comm());
      Epetra_MultiVector B_Pvec(Epetra_DataAccess::View, LocalMap, B.values(), B.stride(), B.numCols());
      
    TEUCHOS_TEST_FOR_EXCEPTION( 
        B_Pvec.Multiply( 'T', 'N', alpha, *A_vec, *this, 0.0 ) != 0,
        EpetraMultiVecFailure, "Anasazi::EpetraMultiVec::MvTransMv() call to Epetra_MultiVec::Multiply() returned a nonzero value.");
    }
  }
  
  //-------------------------------------------------------------
  //
  // b[i] = A[i]^T * this[i]
  // 
  //-------------------------------------------------------------
  
  void EpetraMultiVec::MvDot ( const MultiVec<double>& A, std::vector<double> & b
#ifdef HAVE_ANASAZI_EXPERIMENTAL
                               , ConjType conj
#endif
                             ) const
  {
    EpetraMultiVec *A_vec = dynamic_cast<EpetraMultiVec *>(&const_cast<MultiVec<double> &>(A)); 
    TEUCHOS_TEST_FOR_EXCEPTION( A_vec==NULL,  std::invalid_argument, "Anasazi::EpetraMultiVec::MvDot() cast of MultiVec<double> to EpetraMultiVec failed.");

    if (( (int)b.size() >= A_vec->NumVectors() ) ) {
      TEUCHOS_TEST_FOR_EXCEPTION( 
          this->Dot( *A_vec, &b[0] ) != 0,
          EpetraMultiVecFailure, "Anasazi::EpetraMultiVec::MvDot() call to Epetra_MultiVec::Dot() returned a nonzero value.");
    }
  }

  //-------------------------------------------------------------
  //
  // this[i] = alpha[i] * this[i]
  // 
  //-------------------------------------------------------------
  void EpetraMultiVec::MvScale ( const std::vector<double>& alpha )
  {
    // Check to make sure the vector is as long as the multivector has columns.
    int numvecs = this->NumVectors();
    TEUCHOS_TEST_FOR_EXCEPTION( (int)alpha.size() != numvecs, std::invalid_argument, 
        "Anasazi::EpetraMultiVec::MvScale() alpha argument size was inconsistent with number of vectors in mv.");
    
    std::vector<int> tmp_index( 1, 0 );
    for (int i=0; i<numvecs; i++) {
      Epetra_MultiVector temp_vec(Epetra_DataAccess::View, *this, &tmp_index[0], 1);
      TEUCHOS_TEST_FOR_EXCEPTION( 
          temp_vec.Scale( alpha[i] ) != 0,
          EpetraMultiVecFailure, "Anasazi::EpetraMultiVec::MvScale() call to Epetra_MultiVec::Scale() returned a nonzero value.");
      tmp_index[0]++;
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
  EpetraOp::EpetraOp(const Teuchos::RCP<Epetra_Operator> &Op) 
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
    Epetra_MultiVector* vec_X = dynamic_cast<EpetraMultiVecAccessor*>(&temp_X)->GetEpetraMultiVec();
    Epetra_MultiVector* vec_Y = dynamic_cast<EpetraMultiVecAccessor*>(&Y)->GetEpetraMultiVec();
    
    TEUCHOS_TEST_FOR_EXCEPTION( vec_X==NULL, std::invalid_argument, "Anasazi::EpetraOp::Apply() cast of MultiVec<double> to Epetra_MultiVector failed.");
    TEUCHOS_TEST_FOR_EXCEPTION( vec_Y==NULL, std::invalid_argument, "Anasazi::EpetraOp::Apply() cast of MultiVec<double> to Epetra_MultiVector failed.");

    int info = Epetra_Op->Apply( *vec_X, *vec_Y );
    TEUCHOS_TEST_FOR_EXCEPTION( info != 0, OperatorError, 
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
  
  EpetraGenOp::EpetraGenOp(const Teuchos::RCP<Epetra_Operator> &AOp,
                           const Teuchos::RCP<Epetra_Operator> &MOp,
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
    Epetra_MultiVector* vec_X = dynamic_cast<EpetraMultiVecAccessor*>(&temp_X)->GetEpetraMultiVec();
    Epetra_MultiVector* vec_Y = dynamic_cast<EpetraMultiVecAccessor*>(&Y)->GetEpetraMultiVec();
    Epetra_MultiVector temp_Y(*vec_Y); 
    
    TEUCHOS_TEST_FOR_EXCEPTION( vec_X==NULL, std::invalid_argument, "Anasazi::EpetraGenOp::Apply() cast of MultiVec<double> to Epetra_MultiVector failed.");
    TEUCHOS_TEST_FOR_EXCEPTION( vec_Y==NULL, std::invalid_argument, "Anasazi::EpetraGenOp::Apply() cast of MultiVec<double> to Epetra_MultiVector failed.");
    //
    // Need to cast away constness because the member function Apply is not declared const.  
    // Change the transpose setting for the operator if necessary and change it back when done.
    //
    // Apply M
    info = Epetra_MOp->Apply( *vec_X, temp_Y );
    TEUCHOS_TEST_FOR_EXCEPTION( info != 0, OperatorError, 
                        "Anasazi::EpetraGenOp::Apply(): Error returned from Epetra_Operator::Apply()" );
    // Apply A or A^{-1}
    if (isAInverse) {
      info = Epetra_AOp->ApplyInverse( temp_Y, *vec_Y );
    }
    else {
      info = Epetra_AOp->Apply( temp_Y, *vec_Y );
    }
    TEUCHOS_TEST_FOR_EXCEPTION( info != 0, OperatorError, 
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
  EpetraSymOp::EpetraSymOp(const Teuchos::RCP<Epetra_Operator> &Op, 
                           bool isTrans) 
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
    Epetra_MultiVector* vec_X = dynamic_cast<EpetraMultiVecAccessor*>(&temp_X)->GetEpetraMultiVec();
    Epetra_MultiVector* vec_Y = dynamic_cast<EpetraMultiVecAccessor*>(&Y)->GetEpetraMultiVec();
    Epetra_MultiVector* temp_vec = new Epetra_MultiVector( 
        (isTrans_) ? Epetra_Op->OperatorDomainMap() 
        : Epetra_Op->OperatorRangeMap(), 
        vec_X->NumVectors() );
    
    TEUCHOS_TEST_FOR_EXCEPTION( vec_X==NULL   , std::invalid_argument, "Anasazi::EpetraSymOp::Apply() cast of MultiVec<double> to Epetra_MultiVector failed.");
    TEUCHOS_TEST_FOR_EXCEPTION( vec_Y==NULL   , std::invalid_argument, "Anasazi::EpetraSymOp::Apply() cast of MultiVec<double> to Epetra_MultiVector failed.");
    TEUCHOS_TEST_FOR_EXCEPTION( temp_vec==NULL, std::invalid_argument, "Anasazi::EpetraSymOp::Apply() allocation Epetra_MultiVector failed.");
    //
    // Need to cast away constness because the member function Apply
    // is not declared const.
    //
    // Transpose the operator (if isTrans_ = true)
    if (isTrans_) {
      info = Epetra_Op->SetUseTranspose( isTrans_ );
      if (info != 0) {
        delete temp_vec;
        TEUCHOS_TEST_FOR_EXCEPTION( true, OperatorError, 
                            "Anasazi::EpetraSymOp::Apply(): Error returned from Epetra_Operator::Apply()" );
      }
    }
    //
    // Compute A*X or A'*X 
    //
    info=Epetra_Op->Apply( *vec_X, *temp_vec );
    if (info!=0) { 
      delete temp_vec; 
      TEUCHOS_TEST_FOR_EXCEPTION( true, OperatorError, 
                          "Anasazi::EpetraSymOp::Apply(): Error returned from Epetra_Operator::Apply()" );
    }
    //
    // Transpose/Un-transpose the operator based on value of isTrans_
    info=Epetra_Op->SetUseTranspose( !isTrans_ );
    if (info!=0) { 
      delete temp_vec; 
      TEUCHOS_TEST_FOR_EXCEPTION( true, OperatorError, 
                          "Anasazi::EpetraSymOp::Apply(): Error returned from Epetra_Operator::Apply()" );
    }
    
    // Compute A^T*(A*X) or A*A^T
    info=Epetra_Op->Apply( *temp_vec, *vec_Y );
    if (info!=0) { 
      delete temp_vec; 
      TEUCHOS_TEST_FOR_EXCEPTION( true, OperatorError, 
                          "Anasazi::EpetraSymOp::Apply(): Error returned from Epetra_Operator::Apply()" );
    }
    
    // Un-transpose the operator
    info=Epetra_Op->SetUseTranspose( false );
    delete temp_vec;
    TEUCHOS_TEST_FOR_EXCEPTION( info != 0, OperatorError, 
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
  EpetraSymMVOp::EpetraSymMVOp(const Teuchos::RCP<const Epetra_MultiVector> &MV, 
                               bool isTrans) 
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
    Epetra_MultiVector* vec_X = dynamic_cast<EpetraMultiVecAccessor*>(&temp_X)->GetEpetraMultiVec();
    Epetra_MultiVector* vec_Y = dynamic_cast<EpetraMultiVecAccessor*>(&Y)->GetEpetraMultiVec();
    
    if (isTrans_) {
      
      Epetra_MultiVector temp_vec( *MV_localmap, temp_X.GetNumberVecs() );
      
      /* A'*X */
      info = temp_vec.Multiply( 'T', 'N', 1.0, *Epetra_MV, *vec_X, 0.0 );
      TEUCHOS_TEST_FOR_EXCEPTION( info != 0, OperatorError, 
                          "Anasazi::EpetraSymMVOp::Apply(): Error returned from Epetra_MultiVector::Multiply()" );
      
      /* A*(A'*X) */
      info = vec_Y->Multiply( 'N', 'N', 1.0, *Epetra_MV, temp_vec, 0.0 );      
      TEUCHOS_TEST_FOR_EXCEPTION( info != 0, OperatorError, 
                          "Anasazi::EpetraSymMVOp::Apply(): Error returned from Epetra_MultiVector::Multiply()" );
    } 
    else {
      
      Epetra_MultiVector temp_vec( *MV_blockmap, temp_X.GetNumberVecs() );
      
      /* A*X */
      info = temp_vec.Multiply( 'N', 'N', 1.0, *Epetra_MV, *vec_X, 0.0 );
      TEUCHOS_TEST_FOR_EXCEPTION( info != 0, OperatorError, 
                          "Anasazi::EpetraSymMVOp::Apply(): Error returned from Epetra_MultiVector::Multiply()" );
      
      /* A'*(A*X) */
      info = vec_Y->Multiply( 'T', 'N', 1.0, *Epetra_MV, temp_vec, 0.0 );
      TEUCHOS_TEST_FOR_EXCEPTION( info != 0, OperatorError, 
                          "Anasazi::EpetraSymMVOp::Apply(): Error returned from Epetra_MultiVector::Multiply()" );
    }
  }

  ///////////////////////////////////////////////////////////////////////////////
  //
  //--------Anasazi::EpetraWSymMVOp Implementation--------------------------------
  //
  ///////////////////////////////////////////////////////////////////////////////

  //
  // Anasazi::Operator constructors
  //
  EpetraWSymMVOp::EpetraWSymMVOp(const Teuchos::RCP<const Epetra_MultiVector> &MV, 
                                 const Teuchos::RCP<Epetra_Operator> &OP ) 
    : Epetra_MV(MV), Epetra_OP(OP)
  {
      MV_blockmap = Teuchos::rcp( &Epetra_MV->Map(), false );
      Epetra_WMV = Teuchos::rcp( new Epetra_MultiVector( *MV_blockmap, Epetra_MV->NumVectors() ) ); 
      Epetra_OP->Apply( *Epetra_MV, *Epetra_WMV );
  }
  
  //
  // AnasaziOperator applications
  //
  void EpetraWSymMVOp::Apply ( const MultiVec<double>& X, MultiVec<double>& Y ) const 
  {
    int info=0;
    MultiVec<double> & temp_X = const_cast<MultiVec<double> &>(X);
    Epetra_MultiVector* vec_X = dynamic_cast<EpetraMultiVecAccessor*>(&temp_X)->GetEpetraMultiVec();
    Epetra_MultiVector* vec_Y = dynamic_cast<EpetraMultiVecAccessor*>(&Y)->GetEpetraMultiVec();
    
    Epetra_MultiVector temp_vec( *MV_blockmap, temp_X.GetNumberVecs() );
      
    /* WA*X */
    info = temp_vec.Multiply( 'N', 'N', 1.0, *Epetra_WMV, *vec_X, 0.0 );
    TEUCHOS_TEST_FOR_EXCEPTION( info != 0, OperatorError, 
                        "Anasazi::EpetraWSymMVOp::Apply(): Error returned from Epetra_MultiVector::Multiply()" );
      
    /* A'*(WA*X) */
    info = vec_Y->Multiply( 'T', 'N', 1.0, *Epetra_MV, temp_vec, 0.0 );
    TEUCHOS_TEST_FOR_EXCEPTION( info != 0, OperatorError, 
                        "Anasazi::EpetraWSymMVOp::Apply(): Error returned from Epetra_MultiVector::Multiply()" );
  }
  
  ///////////////////////////////////////////////////////////////////////////////
  //
  //--------Anasazi::EpetraW2SymMVOp Implementation--------------------------------
  //
  ///////////////////////////////////////////////////////////////////////////////

  //
  // Anasazi::Operator constructors
  //
  EpetraW2SymMVOp::EpetraW2SymMVOp(const Teuchos::RCP<const Epetra_MultiVector> &MV, 
                                 const Teuchos::RCP<Epetra_Operator> &OP ) 
    : Epetra_MV(MV), Epetra_OP(OP)
  {
      MV_blockmap = Teuchos::rcp( &Epetra_MV->Map(), false );
      Epetra_WMV = Teuchos::rcp( new Epetra_MultiVector( *MV_blockmap, Epetra_MV->NumVectors() ) ); 
      Epetra_OP->Apply( *Epetra_MV, *Epetra_WMV );
  }
  
  //
  // AnasaziOperator applications
  //
  void EpetraW2SymMVOp::Apply ( const MultiVec<double>& X, MultiVec<double>& Y ) const 
  {
    int info=0;
    MultiVec<double> & temp_X = const_cast<MultiVec<double> &>(X);
    Epetra_MultiVector* vec_X = dynamic_cast<EpetraMultiVecAccessor*>(&temp_X)->GetEpetraMultiVec();
    Epetra_MultiVector* vec_Y = dynamic_cast<EpetraMultiVecAccessor*>(&Y)->GetEpetraMultiVec();
    
    Epetra_MultiVector temp_vec( *MV_blockmap, temp_X.GetNumberVecs() );
      
    /* WA*X */
    info = temp_vec.Multiply( 'N', 'N', 1.0, *Epetra_WMV, *vec_X, 0.0 );
    TEUCHOS_TEST_FOR_EXCEPTION( info != 0, OperatorError, 
                        "Anasazi::EpetraW2SymMVOp::Apply(): Error returned from Epetra_MultiVector::Multiply()" );
      
    /* (WA)'*(WA*X) */
    info = vec_Y->Multiply( 'T', 'N', 1.0, *Epetra_WMV, temp_vec, 0.0 );
    TEUCHOS_TEST_FOR_EXCEPTION( info != 0, OperatorError, 
                        "Anasazi::EpetraW2SymMVOp::Apply(): Error returned from Epetra_MultiVector::Multiply()" );
    
  }
} // end namespace Anasazi
