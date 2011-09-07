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
  /// This class wraps an \c Epetra_Operator, when the wrapped
  /// operator is a preconditioner or other object that is normally
  /// applied using ApplyInverse().  EpetraPrecOp's \c Apply() method
  /// thus invokes the underlying object's
  /// Epetra_Operator::ApplyInverse() method, and its \c
  /// ApplyInverse() method invokes the underlying object's
  /// Epetra_Operator::Apply() method.
  ///
  /// EpetraPrecOp implements both \c Belos::Operator and
  /// Epetra_Operator.  Thus, you can use it in Belos' templated
  /// solvers with either the OP = Belos::Operator or OP =
  /// Epetra_Operator specializations.
  ///
  class EpetraPrecOp : 
    public virtual Operator<double>, 
    public virtual Epetra_Operator 
  {
  public:
    /// Basic constructor.
    ///
    /// \param Op [in/out] The operator to wrap.  EpetraPrecOp's
    ///   Apply() method will invoke this operator's ApplyInverse()
    ///   method, and vice versa.
    EpetraPrecOp (const Teuchos::RCP<Epetra_Operator>& Op);

    //! Virtual destructor, for memory safety of derived classes.
    virtual ~EpetraPrecOp() {};

    /// \brief Apply the operator (or its transpose) to x, putting the result in y.
    ///
    /// This method is part of the Belos::MultiVec implementation.
    /// EpetraPrecOp's Apply() methods invoke the underlying
    /// operator's ApplyInverse() method.  
    ///
    /// All Epetra operators are real-valued, never complex-valued, so
    /// setting either trans=CONJTRANS or trans=TRANS means that this
    /// method will attempt to apply the transpose.  If you attempt to
    /// apply the transpose, EpetraPrecOp will invoke the underlying
    /// operator's SetUseTranspose() method.  This may or may not
    /// succeed.  If it does <i>not</i> succeed, this method will
    /// throw an \c EpetraOpFailure exception.  Furthermore, if the
    /// underlying operator's ApplyInverse() method does not succeed
    /// (i.e., returns a nonzero error code), this method will throw
    /// an \c EpetraOpFailure exception.
    ///
    /// \note The trans argument will always have its literal meaning,
    ///   even if the underlying Epetra_Operator's transpose flag is
    ///   set (i.e., if UseTranspose()==true).  Implementing this
    ///   requires temporarily changing the transpose flag of the
    ///   underlying operator.  However, the flag will be changed back
    ///   to its original value before this method returns, whatever
    ///   that original value was.  This behavior is <i>different</i>
    ///   than that of the two-argument version of Apply() below,
    ///   because the three-argument version of Apply() here
    ///   implements the Belos::Operator and Belos::OperatorTraits
    ///   interfaces.  Those interfaces expect the transpose-ness of
    ///   the operator to be stateless.
    void 
    Apply (const MultiVec<double>& x, 
	   MultiVec<double>& y, 
	   ETrans trans=NOTRANS) const;

    /// \brief Apply the operator to x, putting the result in y.
    ///
    /// This method is part of the Epetra_MultiVector implementation.
    /// EpetraPrecOp's Apply() methods invoke the underlying
    /// operator's ApplyInverse() method.  This version of Apply()
    /// does <i>not</i> attempt to check for errors, but it returns
    /// the error code that the underlying operator's ApplyInverse()
    /// method returns.
    /// 
    /// \note If the underlying operator's transpose flag is set
    ///   (i.e., if UseTranspose() returns true), then EpetraPrecOp
    ///   will apply the transpose of the inverse.  This behavior is
    ///   <i>different</i> than that of the three-argument version of
    ///   Apply() above.  This is so because the two-argument version
    ///   of Apply() implements Epetra_Operator, which expects the
    ///   effect of SetUseTranspose() to persist for all subsequent
    ///   calls to Apply() and ApplyInverse().
    ///
    /// \return Zero if successful, else nonzero.  The value of the
    ///   error code is the same as would be returned by the
    ///   underlying operator's ApplyInverse() method.
    int 
    Apply (const Epetra_MultiVector &X, 
	   Epetra_MultiVector &Y) const;

    /// \brief Apply inverse method for an Epetra_MultiVector.
    ///
    /// This method is part of the Epetra_MultiVector implementation.
    /// EpetraPrecOp's ApplyInverse() method invokes the underlying
    /// operator's Apply() method.  This method does <i>not</i>
    /// attempt to check for errors, but it returns the error code
    /// that the underlying operator's Apply() method returns.
    ///
    /// \note If the underlying operator's transpose flag is set
    ///   (i.e., if UseTranspose() returns true), then EpetraPrecOp
    ///   will apply the transpose of the underlying operator.  This
    ///   behavior is the same as that of the two-argument version of
    ///   Apply(), for the same reason as discussed in that method's
    ///   documentation.
    ///
    /// \return Zero if successful, else nonzero.  The value of the
    ///   error code is the same as would be returned by the
    ///   underlying operator's Apply() method.
    int 
    ApplyInverse (const Epetra_MultiVector &X, 
		  Epetra_MultiVector &Y) const;

    //! Return a human-readable string describing the operator.
    const char* Label() const { 
      return "Epetra_Operator applying A^{-1} as A"; 
    }

    //! Return the current UseTranspose setting.
    bool UseTranspose() const { 
      return Epetra_Op->UseTranspose (); 
    }

    /// \brief If set true, the transpose of this operator will be applied.
    ///
    /// \param UseTranspose_in [in] True if you want to apply the
    ///   transpose, else false.  Epetra operators are real-valued,
    ///   never complex-valued, so applying the conjugate transpose is
    ///   the same as applying the transpose.
    ///
    /// \return Zero if setting was successful, else nonzero.
    ///
    /// This method invokes the underlying operator's
    /// SetUseTranspose().  That operator returns zero if successful,
    /// else nonzero.  Note that SetUseTranspose() affects all
    /// subsequent applications of the operator, until the next
    /// SetUseTranspose() call.
    int SetUseTranspose (bool UseTranspose_in) { 
      return Epetra_Op->SetUseTranspose (UseTranspose_in);
    }

    /// \brief Return true if this object can provide an approximate inf-norm.
    ///
    /// If this method returns false, then the \c NormInf() method
    /// should not be used.
    bool HasNormInf () const { 
      return Epetra_Op->HasNormInf ();
    }

    /// \brief Return the infinity norm of the global matrix.
    ///
    /// The returned value only makes sense if HasNormInf() == true.
    double NormInf() const  { 
      return Epetra_Op->NormInf ();
    }

    //! Return the Epetra_Comm communicator associated with this operator.
    const Epetra_Comm& Comm() const { return Epetra_Op->Comm(); };

    //! Return the Epetra_Map object representing the domain of this operator.
    const Epetra_Map& OperatorDomainMap() const { 
      return Epetra_Op->OperatorDomainMap(); 
    }

    //! Return the Epetra_Map object representing the range of this operator.
    const Epetra_Map& OperatorRangeMap() const { 
      return Epetra_Op->OperatorRangeMap();
    }
    
  private:
    //! The underlying operator that this EpetraPrecOp instance wraps.
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

    static void MvRandom (Epetra_MultiVector& mv)
    { 
      TEST_FOR_EXCEPTION( mv.Random()!=0, EpetraMultiVecFailure, 
			  "Belos::MultiVecTraits<double,Epetra_MultiVector>::"
			  "MvRandom() call to Random() returned a nonzero value.");
    }

    static void MvInit (Epetra_MultiVector& mv, 
			double alpha = Teuchos::ScalarTraits<double>::zero())
    { 
      TEST_FOR_EXCEPTION( mv.PutScalar(alpha)!=0, EpetraMultiVecFailure, 
			  "Belos::MultiVecTraits<double,Epetra_MultiVector>::"
			  "MvInit() call to PutScalar() returned a nonzero value.");
    }

    static void MvPrint (const Epetra_MultiVector& mv, std::ostream& os)
    { 
      os << mv << std::endl; 
    }

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
    
    /// \brief Specialization of Apply() for Epetra_Operator.
    ///
    /// This method throws an EpetraOpFailure if you attempt to apply
    /// the transpose and the underlying operator does not implement
    /// the transpose.
    static void 
    Apply (const Epetra_Operator& Op, 
	   const Epetra_MultiVector& x, 
	   Epetra_MultiVector& y,
	   ETrans trans=NOTRANS)
    { 
      int info = 0;

      // Remember the original value of the transpose flag, so that we
      // can restore it on exit.
      //
      // Do NOT change this behavior just because you don't want an
      // exception to be thrown.  The semantics of OperatorTraits
      // demand that the Apply() method be STATELESS, and that the
      // 'trans' argument have its literal meaning.  If the underlying
      // operator doesn't implement applying the transpose, then it is
      // RIGHT for this method to throw an exception immediately,
      // rather than to silently compute the wrong thing.  However, if
      // you don't try to apply the transpose, then no exception will
      // be thrown.
      const bool originalTransposeFlag = Epetra_Op->UseTranspose ();

      // Whether we want to apply the transpose of the operator.
      // Recall that Epetra operators are always real-valued, never
      // complex-valued, so the conjugate transpose means the same
      // thing as the transpose.
      const bool newTransposeFlag = (trans != NOTRANS);

      // If necessary, set (or unset) the transpose flag to the value
      // corresponding to 'trans'.
      if (newTransposeFlag != originalTransposeFlag) {
	// Toggle the transpose flag temporarily.  We will restore the
	// original flag's value before this method exits.
	info = const_cast<Epetra_Operator &>(Op).SetUseTranspose (newTransposeFlag);
	TEST_FOR_EXCEPTION(info != 0, EpetraOpFailure,
			   "Belos::OperatorTraits::Apply: Toggling the "
			   "transpose flag of the operator failed, returning a "
			   "nonzero error code of " << info << ".  That probably"
			   " means that the Epetra_Operator object doesn't know "
			   "how to apply its transpose.  Are you perhaps using a"
			   " subclass of Epetra_Operator for which applying the "
			   "transpose is not implemented?");
      }

      info = Op.Apply (x, y);
      TEST_FOR_EXCEPTION(info != 0, EpetraOpFailure, 
			 "Belos::OperatorTraits::Apply: Calling Apply() on the "
			 "Epetra_Operator object failed, returning a nonzero "
			 "error code of " << info << ".");

      // SetUseTranspose() changes the state of the operator, so if
      // applicable, we have to change the state back.
      if (newTransposeFlag != originalTransposeFlag) {
	info = const_cast<Epetra_Operator &>(Op).SetUseTranspose (originalTransposeFlag);
	TEST_FOR_EXCEPTION(info != 0, EpetraOpFailure,
			   "Belos::OperatorTraits::Apply: Resetting the original "
			   "transpose flag value of the Epetra_Operator failed, "
			   "returning a nonzero error code of " << info << ".  "
			   "This is an unusual error, since we were able to call "
			   "its SetUseTranspose method successfully before.  "
			   "This suggests a bug in the subclass of Epetra_Operator "
			   "which you are currently using.  This is probably not a "
			   "Belos bug.");
      }

      // Make sure that the transpose flag has its original value.  If
      // not, we throw std::logic_error instead of EpetraOpFailure,
      // since that's definitely a bug.  It's safe to do this check
      // whether or not we actually had to toggle the transpose flag.
      // Any reasonable implementation of Epetra_Operator should make
      // calling UseTranspose() cheap.
      //
      // Note to code maintainers: The reason we capture the value of
      // UseTranspose() instead of just calling it twice, is that if
      // the exception test really does trigger, then there is
      // something fundamentally wrong.  If something is that wrong,
      // then we want to avoid further calls to the operator's
      // methods.  For example, it could be that the UseTranspose()
      // method is erroneously changing the value of the flag, so if
      // we call that method twice, it might have the right value on
      // the second call.  This would make the resulting exception
      // message confusing.
      const bool finalTransposeFlag = Epetra_Op->UseTranspose ();
      TEST_FOR_EXCEPTION(originalTransposeFlag != finalTransposeFlag,
			 std::logic_error,
			 "Belos::OperatorTraits::Apply: The value of the "
			 "Epetra_Operator's transpose flag changed unexpectedly!"
			 "  The original value at the top of this method was "
			 << originalTransposeFlag << ", and its new value is "
			 << finalTransposeFlag << ".  This suggests either a "
			 "bug in Belos (in the implementation of this routine), "
			 "or a bug in the operator.");
    }
  };

} // end of Belos namespace 

#endif 
// end of file BELOS_EPETRA_ADAPTER_HPP
