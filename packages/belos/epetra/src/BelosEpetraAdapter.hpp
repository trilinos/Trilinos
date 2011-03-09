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
#include "BelosInnerSolver.hpp"
#include "BelosVectorSpaceTraits.hpp"

#ifdef HAVE_BELOS_TSQR
#  include <Epetra_TsqrAdaptor.hpp>
#endif // HAVE_BELOS_TSQR

namespace Belos {

  /// \brief Specialization of VectorSpaceTraits for Epetra objects.
  ///
  /// This specialization lets Belos reason about vector spaces
  /// (called "maps" in Epetra) for Epetra objects: Epetra_MultiVector
  /// for multivectors, and Epetra_Operator for operators.  We use
  /// Epetra_BlockMap as the vector space type, with dynamic casting
  /// to Epetra_Map as necessary if we need to make a deep copy of the
  /// map.
  /// 
  /// \note Since this is a full specialization of the template class,
  ///   the implementation can live in the .cpp file.
  template<>
  class VectorSpaceTraits<Epetra_BlockMap> {
  public:
    typedef Epetra_BlockMap vector_space_type;

    /// \brief Are vectors from the two vector spaces compatible?
    ///
    /// \note This method may require communication, so it may be wise
    ///   to cache the result.
    static bool 
    compatible (const vector_space_type& first, 
		const vector_space_type& second);

    /// \brief Return a persistent view of the given vector space.
    ///
    /// "Persistent" means that the view will outlast the Epetra
    /// distributed object from which the vector space object comes.
    /// This is particularly important with Epetra, since its objects'
    /// methods return vector spaces as const references (which are
    /// not guaranteed to survive destruction of their host
    /// distributed object). 
    ///
    /// The cost of persistence for Epetra objects is that if space is
    /// a weak RCP (likely coming from Teuchos::rcpFromRef() called on
    /// a const reference coming from the Epetra object), a deep copy
    /// will be made.  A deep copy is _not_ made if space is a strong
    /// RCP, so you could avoid unnecessary deep copies by making your
    /// own deep copy once and reusing it as necessary.  Tpetra does
    /// not have this problem.
    static Teuchos::RCP<const vector_space_type>
    persistentView (const Teuchos::RCP<const vector_space_type>& space);

    /// \brief Return a persistent view of the given vector space.
    ///
    /// \note Epetra's VectorSpaceTraits gets an overload of \c
    ///   persistentView() for "const vector_space_type&", since
    ///   passing vector spaces by const reference (rather than by
    ///   RCP) is a common Epetra idiom.
    static Teuchos::RCP<const vector_space_type>
    persistentView (const vector_space_type& space);
  };
 
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

    //! @name Vector space typedefs and methods
    //@{
    
    /// \typedef vector_space_type
    ///
    /// \note If you ask for the domain or range map of an
    ///   Epetra_Operator, you'll get an Epetra_Map, but if you ask
    ///   for the map of an Epetra distributed object (such as an
    ///   Epetra_MultiVector), you'll get an Epetra_BlockMap.  An
    ///   Epetra_Map is-an Epetra_BlockMap, so we use the more generic
    ///   interface, except for when we need to make a deep copy of
    ///   the map (in which case we attempt a dynamic cast to
    ///   Epetra_Map and do the right thing, depending on what we
    ///   get).
    typedef Epetra_BlockMap vector_space_type;

    /// Return a persistent view to the vector space in which x lives.
    ///
    /// "Persistent" means that the vector space object will persist
    /// beyond the scope of x.  For the Epetra specialization, this
    /// generally means that a deep copy of the vector space is made.
    /// The Tpetra specialization does not need to make a deep copy.
    ///
    /// \note The term "range" comes from Thyra; an
    ///   Epetra_MultiVector's Epetra_Map and a Tpetra::MultiVector's
    ///   Tpetra::Map both correspond to the "range" of the
    ///   multivector, i.e., the distribution of its rows.
    static Teuchos::RCP<const vector_space_type> 
    getRange (const Epetra_MultiVector& x);
    //@}

    static Teuchos::RCP<Epetra_MultiVector> 
    Clone (const Epetra_MultiVector& mv, const int outNumVecs)
    { 
      TEST_FOR_EXCEPTION(outNumVecs <= 0, std::invalid_argument,
			 "Belos::MultiVecTraits<double, Epetra_MultiVector>::"
			 "Clone(mv, outNumVecs = " << outNumVecs << "): "
			 "outNumVecs must be positive.");
      // FIXME (mfh 13 Jan 2011) Anasazi currently lets Epetra fill in
      // the entries of the returned multivector with zeros, but Belos
      // does not.  We retain this different behavior for now, but the
      // two versions will need to be reconciled.
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

      // FIXME (mfh 13 Jan 2011) Belos allows A to have more columns
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

    //! @name Vector space typedefs and methods
    //@{
    
    /// \typedef vector_space_type
    ///
    /// OP objects have a domain and range "vector space," which may
    /// or may not be different.  OP objects take MV objects from the
    /// domain space as input, and produce OP objects from the range
    /// space as input.  "Vector space" includes the idea of
    /// distributed-memory data distribution, among other things.
    ///  
    /// \note If you ask for the domain or range map of an
    ///   Epetra_Operator, you'll get an Epetra_Map, but if you ask
    ///   for the map of an Epetra distributed object (such as an
    ///   Epetra_MultiVector), you'll get an Epetra_BlockMap.  An
    ///   Epetra_Map is-an Epetra_BlockMap, so we use the more generic
    ///   interface, except for when we need to make a deep copy of
    ///   the map (in which case we attempt a dynamic cast to
    ///   Epetra_Map and do the right thing, depending on what we
    ///   get).
    typedef Epetra_BlockMap vector_space_type;

    /// Return a persistent view to the domain vector space of A.
    ///
    /// \note The Epetra specialization requires copying the vector
    ///   space (Map, in Epetra terms), since Epetra distributed
    ///   objects only return a const vector space reference that is
    ///   not guaranteed to persist beyond the scope of the
    ///   distributed object.
    static Teuchos::RCP<const vector_space_type> 
    getDomain (const Epetra_Operator& A);

    /// Return a persistent view to the range vector space of A.
    ///
    /// \note The Epetra specialization requires copying the vector
    ///   space (Map, in Epetra terms), since Epetra distributed
    ///   objects only return a const vector space reference that is
    ///   not guaranteed to persist beyond the scope of the
    ///   distributed object.
    static Teuchos::RCP<const vector_space_type> 
    getRange (const Epetra_Operator& A);
    //@}
    
  };

  // Forward declaration of InnerSolver class for EpetraInnerSolver.
  template<class Scalar, class MV, class OP>
  class InnerSolver;
  
  /// \class EpetraInnerSolver
  /// \brief Adaptor between InnerSolver and Epetra_Operator.
  /// 
  /// This wrapper lets you use as an Epetra_Operator any
  /// implementation of Belos::InnerSolver<Scalar, MV, OP> with MV =
  /// Epetra_MultiVector and OP = Epetra_Operator.
  class EpetraInnerSolver : public Epetra_Operator {
  public:
    typedef double scalar_type;
    typedef Epetra_MultiVector multivector_type;
    typedef Epetra_Operator operator_type;
    typedef OperatorTraits<scalar_type, multivector_type, operator_type>::vector_space_type vector_space_type;
    typedef InnerSolver<scalar_type, multivector_type, operator_type> inner_solver_type;

    /// \brief Constructor.
    ///
    /// \param solver [in/out] The actual inner solver implementation.
    EpetraInnerSolver (const Teuchos::RCP<inner_solver_type>& solver);

    /// \brief Return the underlying inner solver object.
    ///
    /// This breach of encapsulation makes EpetraInnerSolver into an
    /// "envelope."  First, the inner solver hides inside an
    /// Epetra_Operator until it gets inside a Belos solver that
    /// recognizes the Epetra_Operator as an EpetraInnerSolver.  Then,
    /// the Belos solver can take the InnerSolver out of the
    /// "envelope," destroy the envelope (by setting its RCP to null)
    /// if it wants, and work directly with the (more feature-rich)
    /// InnerSolver.
    ///
    /// \note This method is declared const in order to cheat
    ///   Belos::LinearProblem into letting the operator act like an
    ///   envelope.  It's technically correct to call this method
    ///   const, since it doesn't let the caller assign to the pointer
    ///   (even though it lets the caller call nonconst methods on the
    ///   InnerSolver).
    Teuchos::RCP<inner_solver_type> getInnerSolver() const {
      return solver_;
    }

    /// \brief Compute Y := solver(Y,X).
    ///
    /// \warning X and Y may not alias one another.
    ///
    /// \note The contents of Y on input may be relevant, depending on
    ///   the inner solver implementation.  For example, Y on input
    ///   may be treated as the initial guess of an iterative solver.
    ///
    /// \return Zero on success, else nonzero.
    int Apply (const Epetra_MultiVector& X, Epetra_MultiVector& Y) const;

    /// \brief Apply the inverse of an InnerSolver (not supported).
    /// 
    /// InnerSolver does not support applying the inverse of the
    /// solver.  It seems like the inverse of solving AX=B should be
    /// applying A.  However, this is not necessarily true, since the
    /// solver may only be approximate.
    ///
    /// \return Zero if the inverse was successfully applied, else nonzero.
    int ApplyInverse (const Epetra_MultiVector& X, Epetra_MultiVector& Y) const;

    /// \brief Set Apply() to apply the transpose if possible (it's not).
    ///
    /// If UseTranspose is true, all subsequent calls to Apply() will
    /// attempt to apply the transpose operator, until \c
    /// SetUseTranspose(false) is called.  If applying the transpose
    /// operator is not supported, this method returns nonzero,
    /// otherwise it returns zero.
    int SetUseTranspose (bool UseTranspose);

    /// Computation of the infinity norm (of the whole inner solver)
    /// is not supported.  This method always throws an exception.
    double NormInf() const;

    //! Label for the Epetra_Operator implementation.
    const char* Label() const;

    //! Can you apply the transpose of the operator? (No, never.)
    bool UseTranspose() const { return false; }

    //! Can you compute the infinity norm of the operator? (No, never.)
    bool HasNormInf() const { return false; }

    /// \brief Temporary reference to the domain vector space of the operator.
    ///
    /// The temporary reference should be valid until destruction of
    /// *this.  For a persistent view, call persistentDomainMap() instead.
    const Epetra_Map& OperatorDomainMap() const;

    /// \brief Temporary reference to the range vector space of the operator.
    ///
    /// The temporary reference should be valid until destruction of
    /// *this.  For a persistent view, call persistentRangeMap() instead.
    const Epetra_Map& OperatorRangeMap() const;

    /// \brief Persistent view of the domain vector space of the operator.
    const Teuchos::RCP<const vector_space_type>& persistentDomainMap() const {
      return domain_;
    }

    /// \brief Persistent view of the range vector space of the operator.
    const Teuchos::RCP<const vector_space_type>& persistentRangeMap() const {
      return range_;
    }

    /// \brief The communicator associated with this Epetra object.
    ///
    /// \return A temporary reference to the communicator, which is
    ///   valid within the scope of *this.
    ///
    /// \note I have interpreted this to mean the communicator
    ///   associated with the range of the operator.  Epetra_CrsMatrix
    ///   returns the communicator associated with the row Map (since
    ///   its Epetra_DistObject parent constructor is called with
    ///   rowMap as an argument), which may be different than the
    ///   range Map.  I would imagine that all of an
    ///   Epetra_CrsMatrix's maps should use the same communicator,
    ///   otherwise sparse matrix-vector multiply wouldn't make sense.
    ///   So it shouldn't really matter which map's communicator to
    ///   return.
    const Epetra_Comm& Comm() const {
      return range_->Comm();
    }

  private:
    /// \brief Default construction is not allowed.
    ///
    /// It's not sensible to construct this object without an
    /// instantiated InnerSolver implementation, so we forbid it
    /// syntactically.
    EpetraInnerSolver ();

    //! The inner solver implementation.
    Teuchos::RCP<inner_solver_type> solver_;

    /// \brief The domain vector space.
    ///
    /// \note We have to keep RCPs of the domain and range vector
    ///   spaces around, because the Tpetra::Operator interface
    ///   requires that we return "const RCP&" instead of "RCP" for
    ///   the domain and range.  Instantiating a new RCP in the method
    ///   and returning a const reference to it would result in a
    ///   dangling reference.
    Teuchos::RCP<const vector_space_type> domain_;

    /// \brief The range vector space.
    /// 
    /// See note on \c domain_.
    Teuchos::RCP<const vector_space_type> range_;    
  };

  /// \brief Specialization of makeInnerSolverOperator() for Epetra objects.
  ///
  /// Take an InnerSolver instance, and wrap it in an implementation
  /// of the Epetra_Operator interface.  That way you can use it
  /// alongside any other implementation of the Epetra_Operator
  /// interface.
  ///
  /// \note This is necessary because Belos' solvers require that the
  ///   preconditioner(s) and the matrix all have the same type (OP).
  template<>
  class InnerSolverTraits<double, Epetra_MultiVector, Epetra_Operator> {
  public:
    typedef double scalar_type;
    typedef Epetra_MultiVector multivector_type;
    typedef Epetra_Operator operator_type;
    typedef InnerSolver<scalar_type, multivector_type, operator_type> inner_solver_type;
    typedef EpetraInnerSolver wrapper_type;

    /// \brief Wrap the given inner solver in a wrapper_type.
    ///
    /// The wrapper_type class implements the operator_type interface,
    /// which can be used directly in Belos.
    static Teuchos::RCP<operator_type>
    makeInnerSolverOperator (const Teuchos::RCP<inner_solver_type>& solver)
    {
      using Teuchos::rcp;
      using Teuchos::rcp_implicit_cast;
      return rcp_implicit_cast<operator_type> (rcp (new wrapper_type (solver)));
    }

    /// \brief Return the given wrapper's inner solver object.
    ///
    /// If op is an inner solver wrapper instance, return the inner
    /// solver object.  Otherwise, throw an std::bad_cast exception.
    ///
    /// \note After calling this method, the inner solver object will
    ///   persist beyond the scope of op.  Thus, if you don't care
    ///   about the wrapper that implements the operator_type
    ///   interface, you can get rid of the wrapper (by setting the
    ///   RCP to null) and keep the inner solver.
    static Teuchos::RCP<inner_solver_type>
    getInnerSolver (const Teuchos::RCP<operator_type>& op)
    {
      using Teuchos::RCP;
      using Teuchos::rcp_dynamic_cast;
      RCP<wrapper_type> wrapper = rcp_dynamic_cast<wrapper_type> (op, true);
      return wrapper->getInnerSolver();
    }
  };
  
} // end of Belos namespace 

#endif 
// end of file BELOS_EPETRA_ADAPTER_HPP
