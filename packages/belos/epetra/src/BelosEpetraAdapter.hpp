//@HEADER
// ************************************************************************
//
//                 Belos: Block Linear Solvers Package
//                  Copyright 2004 Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ************************************************************************
//@HEADER

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
#include "BelosTypes.hpp"

#ifdef HAVE_BELOS_TSQR
#  include <Epetra_TsqrAdaptor.hpp>
#endif // HAVE_BELOS_TSQR

namespace Belos {
 
  //! @name Epetra Adapter Exceptions
  //@{

  /** \brief EpetraMultiVecFailure is thrown when a return value from an Epetra
   * call on an Epetra_MultiVector is non-zero.
   */
  class EpetraMultiVecFailure : public BelosError {public:
    EpetraMultiVecFailure(const std::string& what_arg) : BelosError(what_arg)
    {}};

  /** \brief EpetraOpFailure is thrown when a return value from an Epetra
   * call on an Epetra_Operator is non-zero.
   */
  class EpetraOpFailure : public BelosError {public:
    EpetraOpFailure(const std::string& what_arg) : BelosError(what_arg)
    {}};

  //@}
 
  //--------template class BelosEpetraMultiVec-------------------------------------
  class EpetraMultiVec : public MultiVec<double>, public Epetra_MultiVector {
  public:
    // constructors
    EpetraMultiVec(const Epetra_BlockMap& Map_in, double * array, const int numvecs, const int stride=0);
    EpetraMultiVec(const Epetra_BlockMap& Map_in, const int numvecs, bool zeroOut=true);
    EpetraMultiVec(Epetra_DataAccess CV_in, const Epetra_MultiVector& P_vec, const std::vector<int>& index);
    EpetraMultiVec& operator=(const EpetraMultiVec& pv) { Epetra_MultiVector::operator=(pv); return *this; }
    EpetraMultiVec(const Epetra_MultiVector & P_vec);
    ~EpetraMultiVec();
    //
    //  member functions inherited from Belos::MultiVec
    //
    //  the following is a virtual copy constructor returning
    //  a pointer to the pure virtual class. std::vector values are
    //  not copied; instead a new MultiVec is created containing
    //  a non-zero amount of columns.  
    //
    MultiVec<double> * Clone ( const int numvecs ) const;
    //
    //  the following is a virtual copy constructor returning
    //  a pointer to the pure virtual class. std::vector values are
    //  copied and a new stand-alone MultiVector is created.
    //  (deep copy).
    //
    MultiVec<double> * CloneCopy () const;
    //
    //  the following is a virtual copy constructor returning
    //  a pointer to the pure virtual class. std::vector values are
    //  copied and a new stand-alone MultiVector is created
    //  where only selected columns are chosen.  (deep copy).
    //
    MultiVec<double> * CloneCopy ( const std::vector<int>& index ) const;
    //
    //  the following is a virtual view constructor returning
    //  a pointer to the pure virtual class. std::vector values are 
    //  shared and hence no memory is allocated for the columns.
    //
    MultiVec<double> * CloneViewNonConst ( const std::vector<int>& index );
    //
    //  the following is a virtual view constructor returning
    //  a pointer to the pure virtual class. std::vector values are 
    //  shared and hence no memory is allocated for the columns.
    //
    const MultiVec<double> * CloneView ( const std::vector<int>& index ) const;
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

    /*! \brief Scale each element of the vectors in \c *this with \c alpha.
     */
    void MvScale ( const double alpha ) { 
      TEST_FOR_EXCEPTION( this->Scale( alpha )!=0, EpetraMultiVecFailure, 
			  "Belos::EpetraMultiVec::MvScale() call to Scale() returned a nonzero value."); }

    /*! \brief Scale each element of the \c i-th vector in \c *this with \c alpha[i].
     */
    void MvScale ( const std::vector<double>& alpha );
    //
    // B <- alpha * A^T * (*this)
    //
    void MvTransMv ( const double alpha, const MultiVec<double>& A, Teuchos::SerialDenseMatrix<int,double>& B ) const;
    //
    // b[i] = A[i]^T * this[i]
    // 
    void MvDot ( const MultiVec<double>& A, std::vector<double>& b ) const;
    //
    // alpha[i] = norm of i-th column of (*this)
    //	
    void MvNorm ( std::vector<double>& normvec, NormType norm_type = TwoNorm ) const;
    //
    // random vectors in i-th column of (*this)
    //
    void MvRandom() { 
      TEST_FOR_EXCEPTION( Random()!=0, EpetraMultiVecFailure, 
			  "Belos::EpetraMultiVec::MvRandom() call to Random() returned a nonzero value."); }
    //
    // initializes each element of (*this) with alpha
    //
    void MvInit ( const double alpha ) { 
      TEST_FOR_EXCEPTION( PutScalar(alpha)!=0, EpetraMultiVecFailure, 
			  "Belos::EpetraMultiVec::MvInit() call to PutScalar() returned a nonzero value."); }
    //
    // print (*this)
    //
    void MvPrint( std::ostream& os ) const { os << *this << std::endl; };
  private:
  };
  
  
  ///////////////////////////////////////////////////////////////
  //--------template class BelosEpetraOp---------------------
  
  class EpetraOp : public virtual Operator<double> {
  public:
    EpetraOp( const Teuchos::RCP<Epetra_Operator> &Op );
    ~EpetraOp() {};
    void Apply ( const MultiVec<double>& x, MultiVec<double>& y, ETrans trans=NOTRANS ) const;
  private:
    Teuchos::RCP<Epetra_Operator> Epetra_Op;
  };
  
  ///////////////////////////////////////////////////////////////
  //--------template class BelosEpetraPrecOp---------------------
  
  class EpetraPrecOp : public virtual Operator<double>, public virtual Epetra_Operator {
  public:
    //! Basic constructor for applying the operator as its inverse.
    EpetraPrecOp( const Teuchos::RCP<Epetra_Operator> &Op );

    //! Destructor
    virtual ~EpetraPrecOp() {};

    //! Apply method for a Belos::MultiVec [inherited from Belos::Operator class]
    void Apply ( const MultiVec<double>& x, MultiVec<double>& y, ETrans trans=NOTRANS ) const;

    //! Apply method for an Epetra_MultiVector [inherited from Epetra_Operator class]
    int Apply( const Epetra_MultiVector &X, Epetra_MultiVector &Y ) const;

    //! Apply inverse method for an Epetra_MultiVector [inherited from Epetra_Operator class]
    int ApplyInverse( const Epetra_MultiVector &X, Epetra_MultiVector &Y ) const;

    //! Returns a character std::string describing the operator.
    const char* Label() const { return "Epetra_Operator applying A^{-1} as A"; };

    //! Returns the current UseTranspose setting [always false for this operator].
    bool UseTranspose() const { return (false); };

    //! If set true, the transpose of this operator will be applied [not functional for this operator].
    int SetUseTranspose(bool UseTranspose_in) { return 0; };

    //! Returns true if this object can provide an approximate inf-norm [always false for this operator].
    bool HasNormInf() const { return (false); };

    //! Returns the infinity norm of the global matrix [not functional for this operator].
    double NormInf() const  { return (-1.0); };

    //! Returns the Epetra_Comm communicator associated with this operator.
    const Epetra_Comm& Comm() const { return Epetra_Op->Comm(); };

    //! Returns the Epetra_Map object associated with the domain of this operator.
    const Epetra_Map& OperatorDomainMap() const { return Epetra_Op->OperatorDomainMap(); };

    //! Returns the Epetra_Map object associated with the range of this operator.
    const Epetra_Map& OperatorRangeMap() const { return Epetra_Op->OperatorRangeMap(); };
    
  private:

    Teuchos::RCP<Epetra_Operator> Epetra_Op;

  };

  
  ////////////////////////////////////////////////////////////////////
  //
  // Implementation of the Belos::MultiVecTraits for Epetra_MultiVector.
  //
  ////////////////////////////////////////////////////////////////////

  template<>
  class MultiVecTraits<double, Epetra_MultiVector>
  {
  public:
    ///
    static Teuchos::RCP<Epetra_MultiVector> Clone( const Epetra_MultiVector& mv, const int numvecs )
    { return Teuchos::rcp( new Epetra_MultiVector(mv.Map(), numvecs, false) ); }
    ///
    static Teuchos::RCP<Epetra_MultiVector> CloneCopy( const Epetra_MultiVector& mv )
    { return Teuchos::rcp( new Epetra_MultiVector( mv ) ); }
    ///
    static Teuchos::RCP<Epetra_MultiVector> CloneCopy( const Epetra_MultiVector& mv, const std::vector<int>& index )
    { 
      std::vector<int>& tmp_index = const_cast<std::vector<int> &>( index );
      return Teuchos::rcp( new Epetra_MultiVector(Copy, mv, &tmp_index[0], index.size()) ); 
    }
    ///
    static Teuchos::RCP<Epetra_MultiVector> CloneViewNonConst( Epetra_MultiVector& mv, const std::vector<int>& index )
    { 
      std::vector<int>& tmp_index = const_cast<std::vector<int> &>( index );
      return Teuchos::rcp( new Epetra_MultiVector(View, mv, &tmp_index[0], index.size()) ); 
    }
    ///
    static Teuchos::RCP<const Epetra_MultiVector> CloneView( const Epetra_MultiVector& mv, const std::vector<int>& index )
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
    static bool HasConstantStride( const Epetra_MultiVector& mv )
    { return mv.ConstantStride(); }
    ///
    static void MvTimesMatAddMv( const double alpha, const Epetra_MultiVector& A, 
				 const Teuchos::SerialDenseMatrix<int,double>& B, 
				 const double beta, Epetra_MultiVector& mv )
    { 
      Epetra_LocalMap LocalMap(B.numRows(), 0, mv.Map().Comm());
      Epetra_MultiVector B_Pvec(View, LocalMap, B.values(), B.stride(), B.numCols());
      
      int info = mv.Multiply( 'N', 'N', alpha, A, B_Pvec, beta );
      TEST_FOR_EXCEPTION(info!=0, EpetraMultiVecFailure, 
			 "Belos::MultiVecTraits<double,Epetra_MultiVector>::MvTimesMatAddMv call to Multiply() returned a nonzero value.");
    }
    ///
    static void MvAddMv( const double alpha, const Epetra_MultiVector& A, const double beta, const Epetra_MultiVector& B, Epetra_MultiVector& mv )
    { 
      int info = mv.Update( alpha, A, beta, B, 0.0 );
      TEST_FOR_EXCEPTION(info!=0, EpetraMultiVecFailure, 
			 "Belos::MultiVecTraits<double,Epetra_MultiVec>::MvAddMv call to Update() returned a nonzero value.");
    }

    /*! \brief Scale each element of the vectors in \c mv with \c alpha.
     */
    static void MvScale ( Epetra_MultiVector& mv, const double alpha )
    { int ret = mv.Scale( alpha );
      TEST_FOR_EXCEPTION(ret!=0, EpetraMultiVecFailure, 
			 "Belos::MultiVecTraits<double,Epetra_MultiVec>::MvScale call to Scale() returned a nonzero value.");
    }

    /*! \brief Scale each element of the \c i-th vector in \c mv with \c alpha[i].
     */
    static void MvScale ( Epetra_MultiVector& mv, const std::vector<double>& alpha )
    {
      // Check to make sure the vector is as long as the multivector has columns.
      int numvecs = mv.NumVectors();
      TEST_FOR_EXCEPTION((int)alpha.size() != numvecs, EpetraMultiVecFailure, 
			 "Belos::MultiVecTraits<double,Epetra_MultiVec>::MvScale scaling vector (alpha) not same size as number of input vectors (mv).");

      int ret = 0;
      std::vector<int> tmp_index( 1, 0 );
      for (int i=0; i<numvecs; i++) {
        Epetra_MultiVector temp_vec(::View, mv, &tmp_index[0], 1);
        ret = temp_vec.Scale( alpha[i] );
        TEST_FOR_EXCEPTION(ret!=0, EpetraMultiVecFailure, 
	                   "Belos::MultiVecTraits<double,Epetra_MultiVec>::MvScale call to Scale() returned a nonzero value.");
        tmp_index[0]++;
      }
    }

    ///
    static void MvTransMv( const double alpha, const Epetra_MultiVector& A, const Epetra_MultiVector& mv, Teuchos::SerialDenseMatrix<int,double>& B )
    { 
      Epetra_LocalMap LocalMap(B.numRows(), 0, mv.Map().Comm());
      Epetra_MultiVector B_Pvec(View, LocalMap, B.values(), B.stride(), B.numCols());
      
      int info = B_Pvec.Multiply( 'T', 'N', alpha, A, mv, 0.0 );
      TEST_FOR_EXCEPTION(info!=0, EpetraMultiVecFailure, 
			 "Belos::MultiVecTraits<double,Epetra_MultiVector>::MvTransMv call to Multiply() returned a nonzero value.");
    }
    ///
    static void MvDot( const Epetra_MultiVector& mv, const Epetra_MultiVector& A, std::vector<double>& b )
    {
      int info = mv.Dot( A, &b[0] );
      TEST_FOR_EXCEPTION(info!=0, EpetraMultiVecFailure, 
			 "Belos::MultiVecTraits<double,Epetra_MultiVector>::MvDot call to Dot() returned a nonzero value.");   
    }
    ///
    static void MvNorm( const Epetra_MultiVector& mv, std::vector<double>& normvec, NormType type = TwoNorm )
    { 
      if ((int)normvec.size() >= mv.NumVectors()) {
        int info = 0;
	switch( type ) {
	case ( OneNorm ) :
	  info = mv.Norm1(&normvec[0]);
	  break;
	case ( TwoNorm ) :
	  info = mv.Norm2(&normvec[0]);
	  break;
	case ( InfNorm ) :	
	  info = mv.NormInf(&normvec[0]);
	  break;
	default:
	  break;
	}
        TEST_FOR_EXCEPTION(info!=0, EpetraMultiVecFailure, 
			   "Belos::MultiVecTraits<double,Epetra_MultiVector>::MvNorm call to Norm() returned a nonzero value.");
      }
    }
    ///
    static void SetBlock( const Epetra_MultiVector& A, const std::vector<int>& index, Epetra_MultiVector& mv )
    { 
      // Extract the "numvecs" columns of mv indicated by the index std::vector.
      int numvecs = index.size(), info = 0;
      std::vector<int>& tmp_index = const_cast<std::vector<int> &>( index );
      Epetra_MultiVector temp_vec(View, mv, &tmp_index[0], numvecs);
      
      if ( A.NumVectors() != numvecs ) {
        std::vector<int> index2( numvecs );
        for(int i=0; i<numvecs; i++)
	  index2[i] = i;
        Epetra_MultiVector A_vec(View, A, &index2[0], numvecs);      
        info = temp_vec.Update( 1.0, A_vec, 0.0, A_vec, 0.0 );
      }
      else {
        info = temp_vec.Update( 1.0, A, 0.0, A, 0.0 );
      }
      TEST_FOR_EXCEPTION(info!=0, EpetraMultiVecFailure, 
			 "Belos::MultiVecTraits<double,Epetra_MultiVector>::SetBlock call to Update() returned a nonzero value.");
    }
    ///
    static void MvRandom( Epetra_MultiVector& mv )
    { TEST_FOR_EXCEPTION( mv.Random()!=0, EpetraMultiVecFailure, 
			  "Belos::MultiVecTraits<double,Epetra_MultiVector>::MvRandom() call to Random() returned a nonzero value.");
    }
    ///
    static void MvInit( Epetra_MultiVector& mv, double alpha = Teuchos::ScalarTraits<double>::zero() )
    { TEST_FOR_EXCEPTION( mv.PutScalar(alpha)!=0, EpetraMultiVecFailure, 
			  "Belos::MultiVecTraits<double,Epetra_MultiVector>::MvInit() call to PutScalar() returned a nonzero value.");
    }
    ///
    static void MvPrint( const Epetra_MultiVector& mv, std::ostream& os )
    { os << mv << std::endl; }

#ifdef HAVE_BELOS_TSQR
    /// \typedef tsqr_adaptor_type
    /// \brief TsqrAdaptor specialization for Epetra_MultiVector
    ///
    typedef Epetra::TsqrAdaptor tsqr_adaptor_type;
#endif // HAVE_BELOS_TSQR
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
    static void Apply ( const Epetra_Operator& Op, 
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
      TEST_FOR_EXCEPTION(info!=0, EpetraOpFailure, 
			 "Belos::OperatorTraits<double,Epetra_MultiVector,Epetra_Operator>::Apply call to Apply() returned a nonzero value.");
    }
    
  };
  
} // end of Belos namespace 

#endif 
// end of file BELOS_EPETRA_ADAPTER_HPP
