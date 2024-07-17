// @HEADER
// *****************************************************************************
//                 Anasazi: Block Eigensolvers Package
//
// Copyright 2004 NTESS and the Anasazi contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "AnasaziSpecializedEpetraAdapter.hpp"
#include "Teuchos_ScalarTraits.hpp"

/*! \file AnasaziSpecializedEpetraAdapter.cpp
  \brief Implementations of specialized Anasazi multi-vector and operator classes using Epetra_MultiVector and Epetra_Operator
*/

namespace Anasazi {

  ///////////////////////////////////////////////////////////////////////////////
  //
  //--------Anasazi::EpetraOpMultiVec Implementation-------------------------------
  //
  ///////////////////////////////////////////////////////////////////////////////

  // Construction/Destruction
 
   EpetraOpMultiVec::EpetraOpMultiVec(const Teuchos::RCP<Epetra_Operator> &Op, const Epetra_BlockMap& Map_in, const int numvecs)
     : Epetra_OP( Op )
   {
     Epetra_MV = Teuchos::rcp( new Epetra_MultiVector(Map_in, numvecs) );
     Epetra_MV_Temp = Teuchos::rcp( new Epetra_MultiVector(Map_in, numvecs) );
   }    
   
   EpetraOpMultiVec::EpetraOpMultiVec(const Teuchos::RCP<Epetra_Operator> &Op, 
                                      const Epetra_BlockMap& Map_in, double * array, const int numvecs, const int stride)
     : Epetra_OP( Op )
   {
     Epetra_MV = Teuchos::rcp( new Epetra_MultiVector(Epetra_DataAccess::Copy, Map_in, array, stride, numvecs) ); 
     Epetra_MV_Temp = Teuchos::rcp( new Epetra_MultiVector(Map_in, numvecs) ); 
   }
    
   EpetraOpMultiVec::EpetraOpMultiVec(const Teuchos::RCP<Epetra_Operator> &Op, 
                                      Epetra_DataAccess CV, const Epetra_MultiVector& P_vec, const std::vector<int>& index)
     : Epetra_OP( Op )
   {
     Epetra_MV = Teuchos::rcp( new Epetra_MultiVector(CV, P_vec, &(const_cast<std::vector<int> &>(index))[0], index.size()) );
     Epetra_MV_Temp = Teuchos::rcp( new Epetra_MultiVector( P_vec.Map(), index.size()) );
   }

   EpetraOpMultiVec::EpetraOpMultiVec(const EpetraOpMultiVec& P_vec)
     : Epetra_OP( P_vec.Epetra_OP )
   {
     Epetra_MV = Teuchos::rcp( new Epetra_MultiVector( *(P_vec.Epetra_MV) ) );
     Epetra_MV_Temp = Teuchos::rcp( new Epetra_MultiVector( *(P_vec.Epetra_MV_Temp) ) );
   }

  //
  //  member functions inherited from Anasazi::MultiVec
  //
  //
  //  Simulating a virtual copy constructor. If we could rely on the co-variance
  //  of virtual functions, we could return a pointer to EpetraOpMultiVec
  //  (the derived type) instead of a pointer to the pure virtual base class.
  //
  
  MultiVec<double>* EpetraOpMultiVec::Clone ( const int numvecs ) const
  {
    EpetraOpMultiVec *ptr_apv = new EpetraOpMultiVec( Epetra_OP, Epetra_MV->Map(), numvecs );
    return ptr_apv; // safe upcast.
  }
  //
  //  the following is a virtual copy constructor returning
  //  a pointer to the pure virtual class. vector values are
  //  copied.
  //
  
  MultiVec<double>* EpetraOpMultiVec::CloneCopy() const
  {
    EpetraOpMultiVec *ptr_apv = new EpetraOpMultiVec(*this);
    return ptr_apv; // safe upcast
  }
  
  
  MultiVec<double>* EpetraOpMultiVec::CloneCopy ( const std::vector<int>& index ) const
  {
    EpetraOpMultiVec * ptr_apv = new EpetraOpMultiVec( Epetra_OP, Copy, *Epetra_MV, index);
    return ptr_apv; // safe upcast.
  }
  
  
  MultiVec<double>* EpetraOpMultiVec::CloneViewNonConst ( const std::vector<int>& index ) 
  {
    EpetraOpMultiVec * ptr_apv = new EpetraOpMultiVec( Epetra_OP, Epetra_DataAccess::View, *Epetra_MV, index);
    return ptr_apv; // safe upcast.
  }
  
  const MultiVec<double>* EpetraOpMultiVec::CloneView ( const std::vector<int>& index ) const
  {
    EpetraOpMultiVec * ptr_apv = new EpetraOpMultiVec( Epetra_OP, Epetra_DataAccess::View, *Epetra_MV, index);
    return ptr_apv; // safe upcast.
  }
  
  
  void EpetraOpMultiVec::SetBlock( const MultiVec<double>& A, const std::vector<int>& index ) 
  {
    // this should be revisited to e
    EpetraOpMultiVec temp_vec(Epetra_OP, Epetra_DataAccess::View, *Epetra_MV, index);

    int numvecs = index.size();
    if ( A.GetNumberVecs() != numvecs ) {
      std::vector<int> index2( numvecs );
      for(int i=0; i<numvecs; i++)
        index2[i] = i;
      EpetraOpMultiVec *tmp_vec = dynamic_cast<EpetraOpMultiVec *>(&const_cast<MultiVec<double> &>(A)); 
      TEUCHOS_TEST_FOR_EXCEPTION( tmp_vec==NULL, std::invalid_argument, "Anasazi::EpetraOpMultiVec::SetBlocks() cast of MultiVec<double> to EpetraOpMultiVec failed.");
      EpetraOpMultiVec A_vec(Epetra_OP, Epetra_DataAccess::View, *(tmp_vec->GetEpetraMultiVector()), index2);
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
  
  void EpetraOpMultiVec::MvTimesMatAddMv ( double alpha, const MultiVec<double>& A, 
      const Teuchos::SerialDenseMatrix<int,double>& B, double beta ) 
  {
    Epetra_LocalMap LocalMap(B.numRows(), 0, Epetra_MV->Map().Comm());
    Epetra_MultiVector B_Pvec(Epetra_DataAccess::View, LocalMap, B.values(), B.stride(), B.numCols());
    
    EpetraOpMultiVec *A_vec = dynamic_cast<EpetraOpMultiVec *>(&const_cast<MultiVec<double> &>(A)); 
    TEUCHOS_TEST_FOR_EXCEPTION( A_vec==NULL,  std::invalid_argument, "Anasazi::EpetraOpMultiVec::SetBlocks() cast of MultiVec<double> to EpetraOpMultiVec failed.");
    
    TEUCHOS_TEST_FOR_EXCEPTION( 
        Epetra_MV->Multiply( 'N', 'N', alpha, *(A_vec->GetEpetraMultiVector()), B_Pvec, beta ) != 0,
        EpetraSpecializedMultiVecFailure, "Anasazi::EpetraOpMultiVec::MvTimesMatAddMv() call to Epetra_MultiVec::Multiply() returned a nonzero value.");
  }

  //-------------------------------------------------------------
  //
  // *this <- alpha * A + beta * B
  //
  //-------------------------------------------------------------
  
  void EpetraOpMultiVec::MvAddMv ( double alpha , const MultiVec<double>& A, 
                                 double beta, const MultiVec<double>& B) 
  {
    EpetraOpMultiVec *A_vec = dynamic_cast<EpetraOpMultiVec *>(&const_cast<MultiVec<double> &>(A)); 
    TEUCHOS_TEST_FOR_EXCEPTION( A_vec==NULL,  std::invalid_argument, "Anasazi::EpetraOpMultiVec::MvAddMv() cast of MultiVec<double> to EpetraOpMultiVec failed.");
    EpetraOpMultiVec *B_vec = dynamic_cast<EpetraOpMultiVec *>(&const_cast<MultiVec<double> &>(B)); 
    TEUCHOS_TEST_FOR_EXCEPTION( B_vec==NULL,  std::invalid_argument, "Anasazi::EpetraOpMultiVec::MvAddMv() cast of MultiVec<double> to EpetraOpMultiVec failed.");
    
    TEUCHOS_TEST_FOR_EXCEPTION( 
        Epetra_MV->Update( alpha, *(A_vec->GetEpetraMultiVector()), beta, *(B_vec->GetEpetraMultiVector()), 0.0 ) != 0,
        EpetraSpecializedMultiVecFailure, "Anasazi::EpetraOpMultiVec::MvAddMv() call to Epetra_MultiVec::Update() returned a nonzero value.");
  }

  //-------------------------------------------------------------
  //
  // dense B <- alpha * A^T * OP * (*this)
  //
  //-------------------------------------------------------------
  
  void EpetraOpMultiVec::MvTransMv ( double alpha, const MultiVec<double>& A,
                                   Teuchos::SerialDenseMatrix<int,double>& B
#ifdef HAVE_ANASAZI_EXPERIMENTAL
                                   , ConjType conj
#endif
                                  ) const
  {    
    EpetraOpMultiVec *A_vec = dynamic_cast<EpetraOpMultiVec *>(&const_cast<MultiVec<double> &>(A));
    
    if (A_vec) {
      Epetra_LocalMap LocalMap(B.numRows(), 0, Epetra_MV->Map().Comm());
      Epetra_MultiVector B_Pvec(Epetra_DataAccess::View, LocalMap, B.values(), B.stride(), B.numCols());
     
      int info = Epetra_OP->Apply( *Epetra_MV, *Epetra_MV_Temp );
      TEUCHOS_TEST_FOR_EXCEPTION( info != 0, EpetraSpecializedMultiVecFailure, 
        "Anasazi::EpetraOpMultiVec::MvTransMv(): Error returned from Epetra_Operator::Apply()" );

      TEUCHOS_TEST_FOR_EXCEPTION( 
        B_Pvec.Multiply( 'T', 'N', alpha, *(A_vec->GetEpetraMultiVector()), *Epetra_MV_Temp, 0.0 ) != 0,
        EpetraSpecializedMultiVecFailure, "Anasazi::EpetraOpMultiVec::MvTransMv() call to Epetra_MultiVector::Multiply() returned a nonzero value.");
    }
  }
  
  //-------------------------------------------------------------
  //
  // b[i] = A[i]^T * OP * this[i]
  // 
  //-------------------------------------------------------------
  
  void EpetraOpMultiVec::MvDot ( const MultiVec<double>& A, std::vector<double> & b
#ifdef HAVE_ANASAZI_EXPERIMENTAL
                               , ConjType conj
#endif
                             ) const
  {
    EpetraOpMultiVec *A_vec = dynamic_cast<EpetraOpMultiVec *>(&const_cast<MultiVec<double> &>(A)); 
    TEUCHOS_TEST_FOR_EXCEPTION( A_vec==NULL,  std::invalid_argument, "Anasazi::EpetraOpMultiVec::MvDot() cast of MultiVec<double> to EpetraOpMultiVec failed.");

    int info = Epetra_OP->Apply( *Epetra_MV, *Epetra_MV_Temp );
    TEUCHOS_TEST_FOR_EXCEPTION( info != 0, EpetraSpecializedMultiVecFailure, 
      "Anasazi::EpetraOpMultiVec::MvDot(): Error returned from Epetra_Operator::Apply()" );

    if (( (int)b.size() >= A_vec->GetNumberVecs() ) ) {
      TEUCHOS_TEST_FOR_EXCEPTION( 
          Epetra_MV_Temp->Dot( *(A_vec->GetEpetraMultiVector()), &b[0] ) != 0,
          EpetraSpecializedMultiVecFailure, "Anasazi::EpetraOpMultiVec::MvDot() call to Epetra_MultiVector::Dot() returned an error.");
    }
  }

  //-------------------------------------------------------------
  //
  // normvec[i] = || this[i] ||_OP
  // 
  //-------------------------------------------------------------

  void EpetraOpMultiVec::MvNorm ( std::vector<double> & normvec ) const
  {
    int info = Epetra_OP->Apply( *Epetra_MV, *Epetra_MV_Temp );
    TEUCHOS_TEST_FOR_EXCEPTION( info != 0, EpetraSpecializedMultiVecFailure,
      "Anasazi::EpetraOpMultiVec::MvNorm(): Error returned from Epetra_Operator::Apply()" );

    if (( (int)normvec.size() >= Epetra_MV->NumVectors() ) ) {
      TEUCHOS_TEST_FOR_EXCEPTION(
          Epetra_MV_Temp->Dot( *Epetra_MV, &normvec[0] ) != 0,
          EpetraSpecializedMultiVecFailure, "Anasazi::EpetraOpMultiVec::MvNorm() call to Epetra_MultiVector::Dot() returned an error.");
    }

    for (int i=0; i<Epetra_MV->NumVectors(); ++i)
      normvec[i] = Teuchos::ScalarTraits<double>::squareroot( normvec[i] );
  }

  //-------------------------------------------------------------
  //
  // this[i] = alpha[i] * this[i]
  // 
  //-------------------------------------------------------------
  void EpetraOpMultiVec::MvScale ( const std::vector<double>& alpha )
  {
    // Check to make sure the vector is as long as the multivector has columns.
    int numvecs = this->GetNumberVecs();
    TEUCHOS_TEST_FOR_EXCEPTION( (int)alpha.size() != numvecs, std::invalid_argument, 
        "Anasazi::EpetraOpMultiVec::MvScale() alpha argument size was inconsistent with number of vectors in mv.");
    
    std::vector<int> tmp_index( 1, 0 );
    for (int i=0; i<numvecs; i++) {
      Epetra_MultiVector temp_vec(Epetra_DataAccess::View, *Epetra_MV, &tmp_index[0], 1);
      TEUCHOS_TEST_FOR_EXCEPTION( 
          temp_vec.Scale( alpha[i] ) != 0,
          EpetraSpecializedMultiVecFailure, "Anasazi::EpetraOpMultiVec::MvScale() call to Epetra_MultiVector::Scale() returned a nonzero value.");
      tmp_index[0]++;
    }
  }

} // namespace Anasazi
