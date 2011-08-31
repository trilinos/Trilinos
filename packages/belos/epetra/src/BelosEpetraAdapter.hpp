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

  /// \class EpetraMultiVec
  /// \brief Implementation of Belos::MultiVec using Epetra_MultiVector.
  ///
  /// Belos::MultiVec offers a simple abstract interface for
  /// multivector operations in Belos solver algorithms.  This class
  /// implements Belos::MultiVec by extending Epetra_MultiVector.
  class EpetraMultiVec : public MultiVec<double>, public Epetra_MultiVector {
  public:
    // constructors
    EpetraMultiVec(const Epetra_BlockMap& Map_in, double * array, const int numvecs, const int stride=0);
    EpetraMultiVec(const Epetra_BlockMap& Map_in, const int numvecs, bool zeroOut=true);
    EpetraMultiVec(Epetra_DataAccess CV_in, const Epetra_MultiVector& P_vec, const std::vector<int>& index);
    EpetraMultiVec& operator=(const EpetraMultiVec& pv) { Epetra_MultiVector::operator=(pv); return *this; }
    EpetraMultiVec(const Epetra_MultiVector & P_vec);
    ~EpetraMultiVec();

    //! @name Member functions inherited from Belos::MultiVec
    //@{

    /// A virtual "copy constructor" that returns a pointer to a new
    /// object of the pure virtual class.  This vector's entries are
    /// not copied; instead, a new MultiVec is created with the same
    /// data distribution, but with numvecs columns (numvecs > 0).
    /// 
    /// \param numvecs [in] The number of columns in the output
    ///   multivector.  Must be positive.
    MultiVec<double> * Clone ( const int numvecs ) const;

    /// A virtual "copy constructor" returning a pointer to a new
    /// object of the pure virtual class.  This vector's entries are
    /// copied and a new stand-alone multivector is created.  (deep
    /// copy).
    MultiVec<double> * CloneCopy () const;

    /// A virtual "copy constructor" returning a pointer to the pure
    /// virtual class.  This vector's entries are copied and a new
    /// stand-alone MultiVector is created where only selected columns
    /// are chosen.  (deep copy).
    MultiVec<double> * CloneCopy ( const std::vector<int>& index ) const;

    /// A virtual view "constructor" returning a pointer to the pure
    /// virtual class.  This vector's entries are shared and hence no
    /// memory is allocated for the columns.
    MultiVec<double> * CloneViewNonConst ( const std::vector<int>& index );

    /// A virtual view constructor returning a pointer to the pure
    /// virtual class.  This vector's entries are shared and hence no
    /// memory is allocated for the columns.
    const MultiVec<double> * CloneView ( const std::vector<int>& index ) const;

    /// Set a subblock of the multivector, which need not be
    /// contiguous, and is given by the indices.
    void SetBlock ( const MultiVec<double>& A, const std::vector<int>& index );
    
    //! The number of columns in the multivector.
    int GetNumberVecs () const { return NumVectors(); }

    //! The (global) number of rows in the multivector.
    int GetVecLength () const { return GlobalLength(); }

    //! *this <- alpha * A * B + beta * (*this)
    void MvTimesMatAddMv ( const double alpha, const MultiVec<double>& A, 
			   const Teuchos::SerialDenseMatrix<int,double>& B, const double beta );
    //! *this <- alpha * A + beta * B
    void MvAddMv ( const double alpha, const MultiVec<double>& A, const double beta,
		   const MultiVec<double>& B);

    //! Scale each element of the vectors in \c *this with \c alpha.
    void MvScale ( const double alpha ) { 
      TEST_FOR_EXCEPTION( this->Scale( alpha )!=0, EpetraMultiVecFailure, 
			  "Belos::EpetraMultiVec::MvScale() call to Scale() returned a nonzero value."); }

    //! Scale each element of the \c i-th vector in \c *this with \c alpha[i].
    void MvScale ( const std::vector<double>& alpha );

    //! B <- alpha * A^T * (*this)
    void MvTransMv ( const double alpha, const MultiVec<double>& A, Teuchos::SerialDenseMatrix<int,double>& B ) const;

    //! b[i] = A[i]^T * this[i]
    void MvDot ( const MultiVec<double>& A, std::vector<double>& b ) const;

    //! alpha[i] = norm of i-th column of (*this)
    void MvNorm ( std::vector<double>& normvec, NormType norm_type = TwoNorm ) const;

    //! Fill all columns of *this with random values.
    void MvRandom() { 
      TEST_FOR_EXCEPTION( Random()!=0, EpetraMultiVecFailure, 
			  "Belos::EpetraMultiVec::MvRandom() call to Random() returned a nonzero value."); }

    //! Initialize each element of (*this) to the scalar value alpha. 
    void MvInit ( const double alpha ) { 
      TEST_FOR_EXCEPTION( PutScalar(alpha)!=0, EpetraMultiVecFailure, 
			  "Belos::EpetraMultiVec::MvInit() call to PutScalar() returned a nonzero value."); }

    //! Print (*this) to the given output stream.
    void MvPrint( std::ostream& os ) const { os << *this << std::endl; };
  private:
  };
  
  /// \class EpetraOp
  /// \brief Implementation of Belos::Operator using Epetra_Operator.
  ///
  class EpetraOp : public virtual Operator<double> {
  public:
    EpetraOp( const Teuchos::RCP<Epetra_Operator> &Op );
    ~EpetraOp() {};
    void Apply ( const MultiVec<double>& x, MultiVec<double>& y, ETrans trans=NOTRANS ) const;
  private:
    Teuchos::RCP<Epetra_Operator> Epetra_Op;
  };
  

  /// \class EpetraPrecOp
  /// \brief Implementation of Belos::Operator using Epetra_Operator as a preconditioner.
  ///
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

  //! Full specialization of Belos::MultiVecTraits for Epetra_MultiVector.
  template<>
  class MultiVecTraits<double, Epetra_MultiVector>
  {
  public:

    static Teuchos::RCP<Epetra_MultiVector> 
    Clone (const Epetra_MultiVector& mv, const int outNumVecs)
    { 
      TEST_FOR_EXCEPTION(outNumVecs <= 0, std::invalid_argument,
			 "Belos::MultiVecTraits<double, Epetra_MultiVector>::"
			 "Clone(mv, outNumVecs = " << outNumVecs << "): "
			 "outNumVecs must be positive.");
      // NOTE (mfh 13 Jan 2011) Anasazi currently lets Epetra fill in
      // the entries of the returned multivector with zeros, but Belos
      // does not.  We retain this different behavior for now, but the
      // two versions should be reconciled.
      return Teuchos::rcp (new Epetra_MultiVector (mv.Map(), outNumVecs, false)); 
    }

    static Teuchos::RCP<Epetra_MultiVector> 
    CloneCopy (const Epetra_MultiVector& mv)
    { 
      return Teuchos::rcp (new Epetra_MultiVector (mv)); 
    }

    static Teuchos::RCP<Epetra_MultiVector> 
    CloneCopy (const Epetra_MultiVector& mv, const std::vector<int>& index)
    { 
      const int inNumVecs = GetNumberVecs (mv);
      const int outNumVecs = index.size();

      // Simple, inexpensive tests of the index vector.
      TEST_FOR_EXCEPTION(outNumVecs == 0, std::invalid_argument,
			 "Belos::MultiVecTraits<double,Epetra_MultiVector>::"
			 "CloneCopy(mv, index = {}): At least one vector must be"
			 " cloned from mv.");
      if (outNumVecs > inNumVecs)
	{
	  std::ostringstream os;
	  os << "Belos::MultiVecTraits<double,Epetra_MultiVector>::"
	    "CloneCopy(mv, index = {";
	  for (int k = 0; k < outNumVecs - 1; ++k)
	    os << index[k] << ", ";
	  os << index[outNumVecs-1] << "}): There are " << outNumVecs 
	     << " indices to copy, but only " << inNumVecs << " columns of mv.";
	  TEST_FOR_EXCEPTION(true, std::invalid_argument, os.str());
	}
#ifdef TEUCHOS_DEBUG
      // In debug mode, we perform more expensive tests of the index
      // vector, to ensure all the elements are in range.
      // Dereferencing the iterator is valid because index has length
      // > 0.
      const int minIndex = *std::min_element (index.begin(), index.end());
      const int maxIndex = *std::max_element (index.begin(), index.end());

      if (minIndex < 0)
	{
	  std::ostringstream os;
	  os << "Belos::MultiVecTraits<double,Epetra_MultiVector>::"
	    "CloneCopy(mv, index = {";
	  for (int k = 0; k < outNumVecs - 1; ++k)
	    os << index[k] << ", ";
	  os << index[outNumVecs-1] << "}): Indices must be nonnegative, but "
	    "the smallest index " << minIndex << " is negative.";
	  TEST_FOR_EXCEPTION(true, std::invalid_argument, os.str());
	}
      if (maxIndex >= inNumVecs)
	{
	  std::ostringstream os;
	  os << "Belos::MultiVecTraits<double,Epetra_MultiVector>::"
	    "CloneCopy(mv, index = {";
	  for (int k = 0; k < outNumVecs - 1; ++k)
	    os << index[k] << ", ";
	  os << index[outNumVecs-1] << "}): Indices must be strictly less than "
	    "the number of vectors " << inNumVecs << " in mv; the largest index " 
	     << maxIndex << " is out of bounds.";
	  TEST_FOR_EXCEPTION(true, std::invalid_argument, os.str());
	}
#endif // TEUCHOS_DEBUG
      // Cast to nonconst, because Epetra_MultiVector's constructor
      // wants a nonconst int array argument.  It doesn't actually
      // change the entries of the array.
      std::vector<int>& tmpind = const_cast< std::vector<int>& > (index);
      return Teuchos::rcp (new Epetra_MultiVector (Copy, mv, &tmpind[0], index.size())); 
      // return Teuchos::rcp (new Epetra_MultiVector (::Copy, mv, &tmpind[0], index.size())); 
    }

    static Teuchos::RCP<Epetra_MultiVector> 
    CloneCopy (const Epetra_MultiVector& mv, const Teuchos::Range1D& index)
    { 
      const int inNumVecs = GetNumberVecs (mv);
      const int outNumVecs = index.size();
      const bool validRange = outNumVecs > 0 && index.lbound() >= 0 && 
	index.ubound() < inNumVecs;
      if (! validRange)
	{
	  std::ostringstream os;
	  os <<	"Belos::MultiVecTraits<double,Epetra_MultiVector>::Clone(mv,"
	    "index=[" << index.lbound() << ", " << index.ubound() << "]): ";
	  TEST_FOR_EXCEPTION(outNumVecs == 0, std::invalid_argument,
			     os.str() << "Column index range must be nonempty.");
	  TEST_FOR_EXCEPTION(index.lbound() < 0, std::invalid_argument,
			     os.str() << "Column index range must be nonnegative.");
	  TEST_FOR_EXCEPTION(index.ubound() >= inNumVecs, std::invalid_argument,
			     os.str() << "Column index range must not exceed "
			     "number of vectors " << inNumVecs << " in the "
			     "input multivector.");
	}
      return Teuchos::rcp (new Epetra_MultiVector (Copy, mv, index.lbound(), index.size()));
    }

    static Teuchos::RCP<Epetra_MultiVector> 
    CloneViewNonConst (Epetra_MultiVector& mv, const std::vector<int>& index)
    { 
      const int inNumVecs = GetNumberVecs (mv);
      const int outNumVecs = index.size();

      // Simple, inexpensive tests of the index vector.
      TEST_FOR_EXCEPTION(outNumVecs == 0, std::invalid_argument,
			 "Belos::MultiVecTraits<double,Epetra_MultiVector>::"
			 "CloneViewNonConst(mv, index = {}): The output view "
			 "must have at least one column.");
      if (outNumVecs > inNumVecs)
	{
	  std::ostringstream os;
	  os << "Belos::MultiVecTraits<double,Epetra_MultiVector>::"
	    "CloneViewNonConst(mv, index = {";
	  for (int k = 0; k < outNumVecs - 1; ++k)
	    os << index[k] << ", ";
	  os << index[outNumVecs-1] << "}): There are " << outNumVecs 
	     << " indices to view, but only " << inNumVecs << " columns of mv.";
	  TEST_FOR_EXCEPTION(true, std::invalid_argument, os.str());
	}
#ifdef TEUCHOS_DEBUG
      // In debug mode, we perform more expensive tests of the index
      // vector, to ensure all the elements are in range.
      // Dereferencing the iterator is valid because index has length
      // > 0.
      const int minIndex = *std::min_element (index.begin(), index.end());
      const int maxIndex = *std::max_element (index.begin(), index.end());

      if (minIndex < 0)
	{
	  std::ostringstream os;
	  os << "Belos::MultiVecTraits<double,Epetra_MultiVector>::"
	    "CloneViewNonConst(mv, index = {";
	  for (int k = 0; k < outNumVecs - 1; ++k)
	    os << index[k] << ", ";
	  os << index[outNumVecs-1] << "}): Indices must be nonnegative, but "
	    "the smallest index " << minIndex << " is negative.";
	  TEST_FOR_EXCEPTION(true, std::invalid_argument, os.str());
	}
      if (maxIndex >= inNumVecs)
	{
	  std::ostringstream os;
	  os << "Belos::MultiVecTraits<double,Epetra_MultiVector>::"
	    "CloneViewNonConst(mv, index = {";
	  for (int k = 0; k < outNumVecs - 1; ++k)
	    os << index[k] << ", ";
	  os << index[outNumVecs-1] << "}): Indices must be strictly less than "
	    "the number of vectors " << inNumVecs << " in mv; the largest index " 
	     << maxIndex << " is out of bounds.";
	  TEST_FOR_EXCEPTION(true, std::invalid_argument, os.str());
	}
#endif // TEUCHOS_DEBUG
      // Cast to nonconst, because Epetra_MultiVector's constructor
      // wants a nonconst int array argument.  It doesn't actually
      // change the entries of the array.
      std::vector<int>& tmpind = const_cast< std::vector<int>& > (index);
      return Teuchos::rcp (new Epetra_MultiVector (View, mv, &tmpind[0], index.size()));
    }

    static Teuchos::RCP<Epetra_MultiVector> 
    CloneViewNonConst (Epetra_MultiVector& mv, const Teuchos::Range1D& index)
    { 
      const bool validRange = index.size() > 0 && 
	index.lbound() >= 0 && 
	index.ubound() < mv.NumVectors();
      if (! validRange)
	{
	  std::ostringstream os;
	  os <<	"Belos::MultiVecTraits<double,Epetra_MultiVector>::CloneView"
	    "NonConst(mv,index=[" << index.lbound() << ", " << index.ubound() 
	     << "]): ";
	  TEST_FOR_EXCEPTION(index.size() == 0, std::invalid_argument,
			     os.str() << "Column index range must be nonempty.");
	  TEST_FOR_EXCEPTION(index.lbound() < 0, std::invalid_argument,
			     os.str() << "Column index range must be nonnegative.");
	  TEST_FOR_EXCEPTION(index.ubound() >= mv.NumVectors(), 
			     std::invalid_argument,
			     os.str() << "Column index range must not exceed "
			     "number of vectors " << mv.NumVectors() << " in "
			     "the input multivector.");
	}
      return Teuchos::rcp (new Epetra_MultiVector (View, mv, index.lbound(), index.size()));
    }

    static Teuchos::RCP<const Epetra_MultiVector> 
    CloneView (const Epetra_MultiVector& mv, const std::vector<int>& index)
    { 
      const int inNumVecs = GetNumberVecs (mv);
      const int outNumVecs = index.size();

      // Simple, inexpensive tests of the index vector.
      TEST_FOR_EXCEPTION(outNumVecs == 0, std::invalid_argument,
			 "Belos::MultiVecTraits<double,Epetra_MultiVector>::"
			 "CloneView(mv, index = {}): The output view "
			 "must have at least one column.");
      if (outNumVecs > inNumVecs)
	{
	  std::ostringstream os;
	  os << "Belos::MultiVecTraits<double,Epetra_MultiVector>::"
	    "CloneView(mv, index = {";
	  for (int k = 0; k < outNumVecs - 1; ++k)
	    os << index[k] << ", ";
	  os << index[outNumVecs-1] << "}): There are " << outNumVecs 
	     << " indices to view, but only " << inNumVecs << " columns of mv.";
	  TEST_FOR_EXCEPTION(true, std::invalid_argument, os.str());
	}
#ifdef TEUCHOS_DEBUG
      // In debug mode, we perform more expensive tests of the index
      // vector, to ensure all the elements are in range.
      // Dereferencing the iterator is valid because index has length
      // > 0.
      const int minIndex = *std::min_element (index.begin(), index.end());
      const int maxIndex = *std::max_element (index.begin(), index.end());

      if (minIndex < 0)
	{
	  std::ostringstream os;
	  os << "Belos::MultiVecTraits<double,Epetra_MultiVector>::"
	    "CloneView(mv, index = {";
	  for (int k = 0; k < outNumVecs - 1; ++k)
	    os << index[k] << ", ";
	  os << index[outNumVecs-1] << "}): Indices must be nonnegative, but "
	    "the smallest index " << minIndex << " is negative.";
	  TEST_FOR_EXCEPTION(true, std::invalid_argument, os.str());
	}
      if (maxIndex >= inNumVecs)
	{
	  std::ostringstream os;
	  os << "Belos::MultiVecTraits<double,Epetra_MultiVector>::"
	    "CloneView(mv, index = {";
	  for (int k = 0; k < outNumVecs - 1; ++k)
	    os << index[k] << ", ";
	  os << index[outNumVecs-1] << "}): Indices must be strictly less than "
	    "the number of vectors " << inNumVecs << " in mv; the largest index " 
	     << maxIndex << " is out of bounds.";
	  TEST_FOR_EXCEPTION(true, std::invalid_argument, os.str());
	}
#endif // TEUCHOS_DEBUG
      // Cast to nonconst, because Epetra_MultiVector's constructor
      // wants a nonconst int array argument.  It doesn't actually
      // change the entries of the array.
      std::vector<int>& tmpind = const_cast< std::vector<int>& > (index);
      return Teuchos::rcp (new Epetra_MultiVector (View, mv, &tmpind[0], index.size()));
    }

    static Teuchos::RCP<Epetra_MultiVector> 
    CloneView (const Epetra_MultiVector& mv, const Teuchos::Range1D& index)
    { 
      const bool validRange = index.size() > 0 && 
	index.lbound() >= 0 && 
	index.ubound() < mv.NumVectors();
      if (! validRange)
	{
	  std::ostringstream os;
	  os <<	"Belos::MultiVecTraits<double,Epetra_MultiVector>::CloneView"
	    "(mv,index=[" << index.lbound() << ", " << index.ubound() 
	     << "]): ";
	  TEST_FOR_EXCEPTION(index.size() == 0, std::invalid_argument,
			     os.str() << "Column index range must be nonempty.");
	  TEST_FOR_EXCEPTION(index.lbound() < 0, std::invalid_argument,
			     os.str() << "Column index range must be nonnegative.");
	  TEST_FOR_EXCEPTION(index.ubound() >= mv.NumVectors(), 
			     std::invalid_argument,
			     os.str() << "Column index range must not exceed "
			     "number of vectors " << mv.NumVectors() << " in "
			     "the input multivector.");
	}
      return Teuchos::rcp (new Epetra_MultiVector(View, mv, index.lbound(), index.size()));
    }

    static int GetVecLength( const Epetra_MultiVector& mv )
    { return mv.GlobalLength(); }

    static int GetNumberVecs( const Epetra_MultiVector& mv )
    { return mv.NumVectors(); }

    static bool HasConstantStride( const Epetra_MultiVector& mv )
    { return mv.ConstantStride(); }

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

    static void MvAddMv( const double alpha, const Epetra_MultiVector& A, const double beta, const Epetra_MultiVector& B, Epetra_MultiVector& mv )
    { 
      int info = mv.Update( alpha, A, beta, B, 0.0 );
      TEST_FOR_EXCEPTION(info!=0, EpetraMultiVecFailure, 
			 "Belos::MultiVecTraits<double,Epetra_MultiVec>::MvAddMv call to Update() returned a nonzero value.");
    }

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

    static void 
    SetBlock (const Epetra_MultiVector& A, 
	      const std::vector<int>& index, 
	      Epetra_MultiVector& mv)
    { 
      const int inNumVecs = GetNumberVecs (A);
      const int outNumVecs = index.size();

      // NOTE (mfh 13 Jan 2011) Belos allows A to have more columns
      // than index.size(), in which case we just take the first
      // index.size() columns of A.  Anasazi requires that A have the
      // same number of columns as index.size().  Changing Anasazi's
      // behavior should not break existing Anasazi solvers, but the
      // tests need to be done.
      if (inNumVecs < outNumVecs)
	{
	  std::ostringstream os;
	  os << "Belos::MultiVecTraits<double,Epetra_MultiVector>::"
	    "SetBlock(A, mv, index = {";
	  if (outNumVecs > 0)
	    {
	      for (int k = 0; k < outNumVecs - 1; ++k)
		os << index[k] << ", ";
	      os << index[outNumVecs-1];
	    }
	  os << "}): A has only " << inNumVecs << " columns, but there are "
	     << outNumVecs << " indices in the index vector.";
	  TEST_FOR_EXCEPTION(true, std::invalid_argument, os.str());
	}
      // Make a view of the columns of mv indicated by the index std::vector.
      Teuchos::RCP<Epetra_MultiVector> mv_view = CloneViewNonConst (mv, index);

      // View of columns [0, outNumVecs-1] of the source multivector A.
      // If A has fewer columns than mv_view, then create a view of
      // the first outNumVecs columns of A.
      Teuchos::RCP<const Epetra_MultiVector> A_view;
      if (outNumVecs == inNumVecs)
	A_view = Teuchos::rcpFromRef (A); // Const, non-owning RCP
      else
	A_view = CloneView (A, Teuchos::Range1D(0, outNumVecs - 1));

      // Assignment calls Epetra_MultiVector::Assign(), which deeply
      // copies the data directly, ignoring the underlying
      // Epetra_Map(s).  If A and mv don't have the same data
      // distribution (Epetra_Map), this may result in incorrect or
      // undefined behavior.  Epetra_MultiVector::Update() also
      // ignores the Epetra_Maps, so we might as well just use the
      // (perhaps slightly cheaper) Assign() method via operator=().
      *mv_view = *A_view;
    }

    static void 
    SetBlock (const Epetra_MultiVector& A, 
	      const Teuchos::Range1D& index, 
	      Epetra_MultiVector& mv)
    { 
      const int numColsA = A.NumVectors();
      const int numColsMv = mv.NumVectors();
      // 'index' indexes into mv; it's the index set of the target.
      const bool validIndex = index.lbound() >= 0 && index.ubound() < numColsMv;
      // We can't take more columns out of A than A has.
      const bool validSource = index.size() <= numColsA;

      if (! validIndex || ! validSource)
	{
	  std::ostringstream os;
	  os <<	"Belos::MultiVecTraits<double, Epetra_MultiVector>::SetBlock"
	    "(A, index=[" << index.lbound() << ", " << index.ubound() << "], "
	    "mv): ";
	  TEST_FOR_EXCEPTION(index.lbound() < 0, std::invalid_argument,
			     os.str() << "Range lower bound must be nonnegative.");
	  TEST_FOR_EXCEPTION(index.ubound() >= numColsMv, std::invalid_argument,
			     os.str() << "Range upper bound must be less than "
			     "the number of columns " << numColsA << " in the "
			     "'mv' output argument.");
	  TEST_FOR_EXCEPTION(index.size() > numColsA, std::invalid_argument,
			     os.str() << "Range must have no more elements than"
			     " the number of columns " << numColsA << " in the "
			     "'A' input argument.");
	  TEST_FOR_EXCEPTION(true, std::logic_error, "Should never get here!");
	}

      // View of columns [index.lbound(), index.ubound()] of the
      // target multivector mv.  We avoid view creation overhead by
      // only creating a view if the index range is different than [0,
      // (# columns in mv) - 1].
      Teuchos::RCP<Epetra_MultiVector> mv_view;
      if (index.lbound() == 0 && index.ubound()+1 == numColsMv)
	mv_view = Teuchos::rcpFromRef (mv); // Non-const, non-owning RCP
      else
	mv_view = CloneViewNonConst (mv, index);

      // View of columns [0, index.size()-1] of the source multivector
      // A.  If A has fewer columns than mv_view, then create a view
      // of the first index.size() columns of A.
      Teuchos::RCP<const Epetra_MultiVector> A_view;
      if (index.size() == numColsA)
	A_view = Teuchos::rcpFromRef (A); // Const, non-owning RCP
      else
	A_view = CloneView (A, Teuchos::Range1D(0, index.size()-1));

      // Assignment calls Epetra_MultiVector::Assign(), which deeply
      // copies the data directly, ignoring the underlying
      // Epetra_Map(s).  If A and mv don't have the same data
      // distribution (Epetra_Map), this may result in incorrect or
      // undefined behavior.  Epetra_MultiVector::Update() also
      // ignores the Epetra_Maps, so we might as well just use the
      // (perhaps slightly cheaper) Assign() method via operator=().
      *mv_view = *A_view; 
    }

    static void 
    Assign (const Epetra_MultiVector& A, 
	    Epetra_MultiVector& mv)
    { 
      const int numColsA = GetNumberVecs (A);
      const int numColsMv = GetNumberVecs (mv);
      if (numColsA > numColsMv)
	{
	  std::ostringstream os;
	  os <<	"Belos::MultiVecTraits<double, Epetra_MultiVector>::Assign"
	    "(A, mv): ";
	  TEST_FOR_EXCEPTION(numColsA > numColsMv, std::invalid_argument,
			     os.str() << "Input multivector 'A' has " 
			     << numColsA << " columns, but output multivector "
			     "'mv' has only " << numColsMv << " columns.");
	  TEST_FOR_EXCEPTION(true, std::logic_error, "Should never get here!");
	}
      // View of the first [0, numColsA-1] columns of mv.
      Teuchos::RCP<Epetra_MultiVector> mv_view;
      if (numColsMv == numColsA)
	mv_view = Teuchos::rcpFromRef (mv); // Non-const, non-owning RCP
      else // numColsMv > numColsA
	mv_view = CloneView (mv, Teuchos::Range1D(0, numColsA - 1));
      
      // Assignment calls Epetra_MultiVector::Assign(), which deeply
      // copies the data directly, ignoring the underlying
      // Epetra_Map(s).  If A and mv don't have the same data
      // distribution (Epetra_Map), this may result in incorrect or
      // undefined behavior.  Epetra_MultiVector::Update() also
      // ignores the Epetra_Maps, so we might as well just use the
      // (perhaps slightly cheaper) Assign() method via operator=().
      *mv_view = A;
    }

    static void MvRandom( Epetra_MultiVector& mv )
    { TEST_FOR_EXCEPTION( mv.Random()!=0, EpetraMultiVecFailure, 
			  "Belos::MultiVecTraits<double,Epetra_MultiVector>::MvRandom() call to Random() returned a nonzero value.");
    }

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
  // Implementation of the Belos::OperatorTraits for Epetra_Operator.
  //
  ////////////////////////////////////////////////////////////////////
  
  //! Specialization of OperatorTraits for Epetra_Operator.
  template <> 
  class OperatorTraits <double, Epetra_MultiVector, Epetra_Operator>
  {
  public:
    
    //! Specialization of Apply() for Epetra_Operator.
    static void 
    Apply (const Epetra_Operator& Op, 
	   const Epetra_MultiVector& x, 
	   Epetra_MultiVector& y,
	   ETrans trans=NOTRANS)
    { 
      int info = 0;
      // Applying the transpose of an Epetra_Operator requires
      // temporarily setting its "use transpose" flag.  We unset the
      // flag after we are done.  Hopefully you haven't been depending
      // on that flag staying set...
      if (trans)
	const_cast<Epetra_Operator &>(Op).SetUseTranspose (true);
      info = Op.Apply (x, y);
      if (trans)
	const_cast<Epetra_Operator &>(Op).SetUseTranspose (false);      
      TEST_FOR_EXCEPTION(info!=0, EpetraOpFailure, 
			 "Belos::OperatorTraits<double, Epetra_MultiVector, "
			 "Epetra_Operator>::Apply() " << 
			 (trans!=NOTRANS ? "(applying the transpose)" : "")
			 << " returned a nonzero value info=" << info << ", "
			 "indicating an error in applying the operator.");
    }
  };

} // end of Belos namespace 

#endif 
// end of file BELOS_EPETRA_ADAPTER_HPP
