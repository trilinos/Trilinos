// @HEADER
// ***********************************************************************
//
//                 Belos: Block Eigensolvers Package
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

/*! \file BelosThyraAdapter.hpp
  \brief Specializations of the Belos multi-std::vector and operator traits 
  classes using Thyra base classes LinearOpBase and MultiVectorBase.
*/

#ifndef BELOS_THYRA_ADAPTER_HPP
#define BELOS_THYRA_ADAPTER_HPP

#include "BelosConfigDefs.hpp"
#include "BelosMultiVecTraits.hpp"
#include "BelosOperatorTraits.hpp"
#include "BelosVectorSpaceTraits.hpp"
#include "BelosInnerSolver.hpp"

#include <Thyra_DetachedMultiVectorView.hpp>
#include <Thyra_MultiVectorBase.hpp>
#include <Thyra_MultiVectorStdOps.hpp>
#ifdef HAVE_BELOS_TSQR
#  include <Thyra_TsqrAdaptor.hpp>
#endif // HAVE_BELOS_TSQR

namespace Belos {

  /// \brief Specialization of VectorSpaceTraits for Thyra objects.
  /// 
  /// \note Belos' Thyra adaptor uses LinearOpBase for the OP
  ///   (operator) type.  MultiVectorBase "is-a" LinearOpBase, but
  ///   that doesn't matter here.
  template<class Scalar>
  class VectorSpaceTraits<Thyra::VectorSpaceBase<Scalar> > {
  public:
    typedef Thyra::VectorSpaceBase<Scalar> vector_space_type;

    static bool 
    compatible (const vector_space_type& first, 
		const vector_space_type& second) {
      return first.isCompatible (second);
    }

    static Teuchos::RCP<const vector_space_type>
    persistentView (const Teuchos::RCP<const vector_space_type>& space)
    {
      using Teuchos::RCP;

      // Thyra distributed objects return RCP<const
      // VectorSpaceBase<Scalar> > when asked for their data
      // distributions, as do Thyra::LinearOpBase objects when asked
      // for their domain and range.  Those RCPs are persistent
      // (strong references) by construction, so we can just let them
      // pass through.  For safety, we still check whether it's a weak
      // reference, and return a deep copy if it is.
      if (space.strength() == Teuchos::RCP_WEAK)
	{
	  // VectorSpaceBase::clone() is allowed to return null,
	  // indicating that cloning is not supported.
	  RCP<const vector_space_type> copy = space->clone();
	  TEST_FOR_EXCEPTION(copy.is_null(), std::logic_error,
			     "Making a deep copy of the given Thyra::"
			     "VectorSpaceBase object is not supported.");
	  return copy;
	}
      else
	return space;
    }
  };
  
  ////////////////////////////////////////////////////////////////////////////
  //
  // Implementation of the Belos::MultiVecTraits for Thyra::MultiVectorBase
  //
  ////////////////////////////////////////////////////////////////////////////

  /*!  \brief Template specialization of Belos::MultiVecTraits class using the
    Thyra::MultiVectorBase class.

    This interface will ensure that any implementation of MultiVectorBaseClass 
    will be accepted by the Belos templated solvers.  

  */

  template<class ScalarType>
  class MultiVecTraits< ScalarType, Thyra::MultiVectorBase<ScalarType> >
  {
  private:
    typedef Thyra::MultiVectorBase<ScalarType> TMVB;
    typedef Teuchos::ScalarTraits<ScalarType> ST;
    typedef typename ST::magnitudeType magType;
  public:

    //! @name Vector space typedefs and methods
    //@{
    
    /// \typedef vector_space_type
    ///
    /// MV objects live in a "vector space."  Two objects of the same
    /// MV type might live in different vector spaces.  "Vector space"
    /// includes the idea of distributed-memory data distribution,
    /// among other things.
    typedef Thyra::VectorSpaceBase<ScalarType> vector_space_type;

    /// Return a persistent view to the vector space in which x lives.
    ///
    /// "Persistent" means that the vector space object will persist
    /// beyond the scope of x.  "Range" refers to the data
    /// distribution of the rows of the multivector, with the
    /// assumption that it is distributed in block row fashion among
    /// processor(s).
    static Teuchos::RCP<const vector_space_type> 
    getRange (const Thyra::MultiVectorBase<ScalarType>& x) {
      return x.range ();
    }
    //@}

    /** \name Creation methods */
    //@{

    /*! \brief Creates a new empty MultiVectorBase containing \c numvecs columns.
      
    \return Reference-counted pointer to the new MultiVectorBase.
    */
    static Teuchos::RCP<TMVB> Clone( const TMVB& mv, const int numvecs )
    { 
      Teuchos::RCP<TMVB> c = Thyra::createMembers( mv.range(), numvecs ); 
      return c;
    }

    /*! \brief Creates a new MultiVectorBase and copies contents of \c mv into the new std::vector (deep copy).
      
      \return Reference-counted pointer to the new MultiVectorBase.
    */
    static Teuchos::RCP<TMVB> CloneCopy( const TMVB& mv )
    { 
      int numvecs = mv.domain()->dim();
      // create the new multivector
      Teuchos::RCP< TMVB > cc = Thyra::createMembers( mv.range(), numvecs );
      // copy the data from the source multivector to the new multivector
      Thyra::assign(&*cc, mv);
      return cc;
    }

    /*! \brief Creates a new MultiVectorBase and copies the selected contents of \c mv into the new std::vector (deep copy).  

      The copied vectors from \c mv are indicated by the \c indeX.size() indices in \c index.      
      \return Reference-counted pointer to the new MultiVectorBase.
    */
    static Teuchos::RCP<TMVB> CloneCopy( const TMVB& mv, const std::vector<int>& index )
    { 
      int numvecs = index.size();
      // create the new multivector
      Teuchos::RCP<TMVB> cc = Thyra::createMembers( mv.range(), numvecs );
      // create a view to the relevant part of the source multivector
      Teuchos::RCP<const TMVB> view = mv.subView(index);
      // copy the data from the relevant view to the new multivector
      Thyra::assign(&*cc, *view);
      return cc;
    }

    static Teuchos::RCP<TMVB> 
    CloneCopy (const TMVB& mv, const Teuchos::Range1D& index)
    { 
      const int numVecs = index.size();      
      // Create the new multivector
      Teuchos::RCP<TMVB> cc = Thyra::createMembers (mv.range(), numVecs);
      // Create a view to the relevant part of the source multivector
      Teuchos::RCP<const TMVB> view = mv.subView (index);
      // Copy the data from the view to the new multivector.
      // &*cc converts cc from an RCP to a Ptr.
      Thyra::assign (&*cc, *view);
      return cc;
    }

    /*! \brief Creates a new MultiVectorBase that shares the selected contents of \c mv (shallow copy).

    The index of the \c numvecs vectors shallow copied from \c mv are indicated by the indices given in \c index.
    \return Reference-counted pointer to the new MultiVectorBase.
    */      
    static Teuchos::RCP<TMVB> CloneViewNonConst( TMVB& mv, const std::vector<int>& index )
    {
      int numvecs = index.size();

      // We do not assume that the indices are sorted, nor do we check that
      // index.size() > 0. This code is fail-safe, in the sense that a zero
      // length index std::vector will pass the error on the Thyra.

      // Thyra has two ways to create an indexed View:
      // * contiguous (via a range of columns)
      // * indexed (via a std::vector of column indices)
      // The former is significantly more efficient than the latter, in terms of
      // computations performed with/against the created view.
      // We will therefore check to see if the given indices are contiguous, and
      // if so, we will use the contiguous view creation method.

      int lb = index[0];
      bool contig = true;
      for (int i=0; i<numvecs; i++) {
        if (lb+i != index[i]) contig = false;
      }

      Teuchos::RCP< TMVB > cc;
      if (contig) {
        const Thyra::Range1D rng(lb,lb+numvecs-1);
        // create a contiguous view to the relevant part of the source multivector
        cc = mv.subView(rng);
      }
      else {
        // create an indexed view to the relevant part of the source multivector
        cc = mv.subView(index);
      }
      return cc;
    }

    static Teuchos::RCP<TMVB> 
    CloneViewNonConst (TMVB& mv, const Teuchos::Range1D& index)
    {
      // We let Thyra be responsible for checking that the index range
      // is nonempty.
      //
      // Create and return a contiguous view to the relevant part of
      // the source multivector.
      return mv.subView (index);
    }


    /*! \brief Creates a new const MultiVectorBase that shares the selected contents of \c mv (shallow copy).

    The index of the \c numvecs vectors shallow copied from \c mv are indicated by the indices given in \c index.
    \return Reference-counted pointer to the new const MultiVectorBase.
    */      
    static Teuchos::RCP<const TMVB> CloneView( const TMVB& mv, const std::vector<int>& index )
    {
      int numvecs = index.size();

      // We do not assume that the indices are sorted, nor do we check that
      // index.size() > 0. This code is fail-safe, in the sense that a zero
      // length index std::vector will pass the error on the Thyra.

      // Thyra has two ways to create an indexed View:
      // * contiguous (via a range of columns)
      // * indexed (via a std::vector of column indices)
      // The former is significantly more efficient than the latter, in terms of
      // computations performed with/against the created view.
      // We will therefore check to see if the given indices are contiguous, and
      // if so, we will use the contiguous view creation method.

      int lb = index[0];
      bool contig = true;
      for (int i=0; i<numvecs; i++) {
        if (lb+i != index[i]) contig = false;
      }

      Teuchos::RCP< const TMVB > cc;
      if (contig) {
        const Thyra::Range1D rng(lb,lb+numvecs-1);
        // create a contiguous view to the relevant part of the source multivector
        cc = mv.subView(rng);
      }
      else {
        // create an indexed view to the relevant part of the source multivector
        cc = mv.subView(index);
      }
      return cc;
    }

    static Teuchos::RCP<const TMVB> 
    CloneView (const TMVB& mv, const Teuchos::Range1D& index)
    {
      // We let Thyra be responsible for checking that the index range
      // is nonempty.
      //
      // Create and return a contiguous view to the relevant part of
      // the source multivector.
      return mv.subView (index);
    }

    //@}

    /** \name Attribute methods */
    //@{

    //! Obtain the std::vector length of \c mv.
    static int GetVecLength( const TMVB& mv )
    { return mv.range()->dim(); }

    //! Obtain the number of vectors in \c mv
    static int GetNumberVecs( const TMVB& mv )
    { return mv.domain()->dim(); }

    //@}

    /** \name Update methods */
    //@{

    /*! \brief Update \c mv with \f$ \alpha AB + \beta mv \f$.
     */
    static void MvTimesMatAddMv( const ScalarType alpha, const TMVB& A, 
         const Teuchos::SerialDenseMatrix<int,ScalarType>& B, 
         const ScalarType beta, TMVB& mv )
    {
      using Teuchos::arrayView; using Teuchos::arcpFromArrayView;
      const int m = B.numRows();
      const int n = B.numCols();
      // Create a view of the B object!
      Teuchos::RCP< const TMVB >
        B_thyra = Thyra::createMembersView(
          A.domain(),
          RTOpPack::ConstSubMultiVectorView<ScalarType>(
            0, m, 0, n,
            arcpFromArrayView(arrayView(&B(0,0), B.stride()*B.numCols())), B.stride()
            )
          );
      // perform the operation via A: mv <- alpha*A*B_thyra + beta*mv
      Thyra::apply<ScalarType>(A, Thyra::NOTRANS, *B_thyra, Teuchos::outArg(mv), alpha, beta);
    }

    /*! \brief Replace \c mv with \f$\alpha A + \beta B\f$.
     */
    static void MvAddMv( const ScalarType alpha, const TMVB& A, 
                         const ScalarType beta,  const TMVB& B, TMVB& mv )
    { 
      ScalarType coef[2], zero = Teuchos::ScalarTraits<ScalarType>::zero();
      const TMVB * in[2];

      in[0] = &A;
      in[1] = &B;
      coef[0] = alpha;
      coef[1] = beta;

      Thyra::linear_combination(2,coef,in,zero,&mv);
    }

    /*! \brief Scale each element of the vectors in \c *this with \c alpha.
     */
    static void MvScale ( TMVB& mv, const ScalarType alpha )
    { Thyra::scale(alpha,&mv); }

    /*! \brief Scale each element of the \c i-th vector in \c *this with \c alpha[i].
     */
    static void MvScale ( TMVB& mv, const std::vector<ScalarType>& alpha )
    {
      for (unsigned int i=0; i<alpha.size(); i++) {
        Thyra::scale<ScalarType>(alpha[i], mv.col(i).ptr());
      }
    }

    /*! \brief Compute a dense matrix \c B through the matrix-matrix multiply \f$ \alpha A^Tmv \f$.
    */
    static void MvTransMv( const ScalarType alpha, const TMVB& A, const TMVB& mv, 
                           Teuchos::SerialDenseMatrix<int,ScalarType>& B )
    { 
      using Teuchos::arrayView; using Teuchos::arcpFromArrayView;
      // Create a multivector to hold the result (m by n)
      int m = A.domain()->dim();
      int n = mv.domain()->dim();
      // Create a view of the B object!
      Teuchos::RCP< TMVB >
        B_thyra = Thyra::createMembersView(
          A.domain(),
          RTOpPack::SubMultiVectorView<ScalarType>(
            0, m, 0, n,
            arcpFromArrayView(arrayView(&B(0,0), B.stride()*B.numCols())), B.stride()
            )
          );
      Thyra::apply<ScalarType>(A, Thyra::CONJTRANS, mv, B_thyra.ptr(), alpha);
    }

    /*! \brief Compute a std::vector \c b where the components are the individual dot-products of the 
        \c i-th columns of \c A and \c mv, i.e.\f$b[i] = A[i]^Tmv[i]\f$.
     */
    static void MvDot( const TMVB& mv, const TMVB& A, std::vector<ScalarType>& b )
    { Thyra::dots(mv,A,&(b[0])); }

    //@}

    /** \name Norm method */
    //@{

    /*! \brief Compute the 2-norm of each individual std::vector of \c mv.  
      Upon return, \c normvec[i] holds the value of \f$||mv_i||_2\f$, the \c i-th column of \c mv.
    */
    static void MvNorm( const TMVB& mv, std::vector<magType>& normvec, NormType type = TwoNorm )
    { Thyra::norms_2(mv,&(normvec[0])); }

    //@}

    /** \name Initialization methods */
    //@{

    /*! \brief Copy the vectors in \c A to a set of vectors in \c mv indicated by the indices given in \c index.
     */
    static void SetBlock( const TMVB& A, const std::vector<int>& index, TMVB& mv )
    { 
      // Extract the "numvecs" columns of mv indicated by the index std::vector.
      int numvecs = index.size();
      std::vector<int> indexA(numvecs);
      int numAcols = A.domain()->dim();
      for (int i=0; i<numvecs; i++) {
        indexA[i] = i;
      }
      // Thyra::assign requires that both arguments have the same number of 
      // vectors. Enforce this, by shrinking one to match the other.
      if ( numAcols < numvecs ) {
        // A does not have enough columns to satisfy index_plus. Shrink
        // index_plus.
        numvecs = numAcols;
      }
      else if ( numAcols > numvecs ) {
        numAcols = numvecs;
        indexA.resize( numAcols );
      }
      // create a view to the relevant part of the source multivector
      Teuchos::RCP< const TMVB > relsource = A.subView(indexA);
      // create a view to the relevant part of the destination multivector
      Teuchos::RCP< TMVB > reldest = mv.subView(index);
      // copy the data to the destination multivector subview
      Thyra::assign(&*reldest, *relsource);
    }

    static void 
    SetBlock (const TMVB& A, const Teuchos::Range1D& index, TMVB& mv)
    { 
      const int numColsA = A.domain()->dim();
      const int numColsMv = mv.domain()->dim();
      // 'index' indexes into mv; it's the index set of the target.
      const bool validIndex = index.lbound() >= 0 && index.ubound() < numColsMv;
      // We can't take more columns out of A than A has.
      const bool validSource = index.size() <= numColsA;

      if (! validIndex || ! validSource)
	{
	  std::ostringstream os;
	  os << "Belos::MultiVecTraits<Scalar, Thyra::MultiVectorBase<Scalar> "
	    ">::SetBlock(A, [" << index.lbound() << ", " << index.ubound() 
	     << "], mv): ";
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

      // View of the relevant column(s) of the target multivector mv.
      // We avoid view creation overhead by only creating a view if
      // the index range is different than [0, (# columns in mv) - 1].
      Teuchos::RCP<TMVB> mv_view;
      if (index.lbound() == 0 && index.ubound()+1 == numColsMv)
	mv_view = Teuchos::rcpFromRef (mv); // Non-const, non-owning RCP
      else
	mv_view = mv.subView (index);

      // View of the relevant column(s) of the source multivector A.
      // If A has fewer columns than mv_view, then create a view of
      // the first index.size() columns of A.
      Teuchos::RCP<const TMVB> A_view;
      if (index.size() == numColsA)
	A_view = Teuchos::rcpFromRef (A); // Const, non-owning RCP
      else
	A_view = A.subView (Teuchos::Range1D(0, index.size()-1));

      // Copy the data to the destination multivector.  "&*" is a
      // typical Thyra idiom for turning an RCP into a raw pointer.
      Thyra::assign (&*mv_view, *A_view);
    }

    static void 
    Assign (const TMVB& A, TMVB& mv)
    { 
      const int numColsA = A.domain()->dim();
      const int numColsMv = mv.domain()->dim();
      if (numColsA > numColsMv)
	{
	  std::ostringstream os;
	  os << "Belos::MultiVecTraits<Scalar, Thyra::MultiVectorBase<Scalar>"
	    " >::Assign(A, mv): ";
	  TEST_FOR_EXCEPTION(numColsA > numColsMv, std::invalid_argument,
			     os.str() << "Input multivector 'A' has " 
			     << numColsA << " columns, but output multivector "
			     "'mv' has only " << numColsMv << " columns.");
	  TEST_FOR_EXCEPTION(true, std::logic_error, "Should never get here!");
	}
      // Copy the data to the destination multivector.
      if (numColsA == numColsMv)
	Thyra::assign (&mv, A);
      else
	{
	  Teuchos::RCP<TMVB> mv_view = 
	    CloneViewNonConst (mv, Teuchos::Range1D(0, numColsA-1));
	  // "&*" is a typical Thyra idiom for turning an RCP into a
	  // raw pointer.
	  Thyra::assign (&*mv_view, A);
	}
    }

    /*! \brief Replace the vectors in \c mv with random vectors.
     */
    static void MvRandom( TMVB& mv )
    { 
      // Thyra::randomize generates via a uniform distribution on [l,u]
      // We will use this to generate on [-1,1]
      Thyra::randomize(-Teuchos::ScalarTraits<ScalarType>::one(),
                        Teuchos::ScalarTraits<ScalarType>::one(),
                       &mv);
    }

    /*! \brief Replace each element of the vectors in \c mv with \c alpha.
     */
    static void MvInit( TMVB& mv, ScalarType alpha = Teuchos::ScalarTraits<ScalarType>::zero() )
    { Thyra::assign(&mv,alpha); }

    //@}

    /** \name Print method */
    //@{

    /*! \brief Print the \c mv multi-std::vector to the \c os output stream.
     */
    static void MvPrint( const TMVB& mv, std::ostream& os )
      { os << describe(mv,Teuchos::VERB_EXTREME); }

    //@}

#ifdef HAVE_BELOS_TSQR
    /// \typedef tsqr_adaptor_type
    /// \brief TsqrAdaptor specialization for Epetra_MultiVector
    ///
    /// \warning This is a STUB for now; it doesn't work.
    ///
    typedef Thyra::TsqrAdaptor< ScalarType > tsqr_adaptor_type;
#endif // HAVE_BELOS_TSQR
  };        

  ///////////////////////////////////////////////////////////////////////// 
  //
  // Implementation of the Belos::OperatorTraits for Thyra::LinearOpBase
  //
  ///////////////////////////////////////////////////////////////////////// 

  /*!  \brief Template specialization of Belos::OperatorTraits class using the
    Thyra::LinearOpBase virtual base class and Thyra::MultiVectorBase class.

    This interface will ensure that any LinearOpBase and MultiVectorBase
    implementations will be accepted by the Belos templated solvers.

  */

  template <class ScalarType> 
  class OperatorTraits < ScalarType, Thyra::MultiVectorBase<ScalarType>, Thyra::LinearOpBase<ScalarType> >
  {
  private:
    typedef Thyra::MultiVectorBase<ScalarType> TMVB;
    typedef Thyra::LinearOpBase<ScalarType>    TLOB;
  public:
    
    /*! \brief This method takes the MultiVectorBase \c x and
      applies the LinearOpBase \c Op to it resulting in the MultiVectorBase \c y.
    */    
    static void Apply ( const TLOB& Op, const TMVB& x, TMVB& y )
    { 
      Thyra::apply<ScalarType>(Op, Thyra::NOTRANS, x, Teuchos::outArg(y));
    }

    //! @name Vector space typedefs and methods
    //@{
    
    /// \typedef vector_space_type
    ///
    /// Operator objects have a domain and range "vector space," which
    /// may or may not be different.  OP objects take MV objects from
    /// the domain space as input, and produce OP objects from the
    /// range space as input.  "Vector space" includes the idea of
    /// distributed-memory data distribution, among other things.
    typedef Thyra::VectorSpaceBase<ScalarType> vector_space_type;

    /// Return a persistent view to the domain vector space of A.
    ///
    /// "Persistent" means that the vector space object will persist
    /// beyond the scope of A.
    static Teuchos::RCP<const vector_space_type> 
    getDomain (const Thyra::LinearOpBase<ScalarType>& A) {
      return A.domain ();
    }

    /// Return a persistent view to the range vector space of A.
    ///
    /// "Persistent" means that the vector space object will persist
    /// beyond the scope of A.
    static Teuchos::RCP<const vector_space_type> 
    getRange (const Thyra::LinearOpBase<ScalarType>& A) {
      return A.range ();
    }
    //@}
  };


  /// \class ThyraInnerSolver
  /// \brief Adaptor between InnerSolver and Thyra::LinearOpBase.
  /// 
  /// This wrapper lets you use as a Thyra::LinearOpBase any
  /// implementation of Belos::InnerSolver<Scalar, MV, OP> with MV =
  /// Thyra::MultiVectorBase<Scalar> and OP =
  /// Thyra::LinearOpBase<Scalar>.
  /// 
  template <class Scalar>
  class ThyraInnerSolver : public Thyra::LinearOpBase<Scalar> {
  public:
    typedef Scalar scalar_type;
    typedef Thyra::MultiVectorBase<Scalar> multivector_type;
    typedef Thyra::LinearOpBase<Scalar> operator_type;
    typedef typename OperatorTraits<scalar_type, multivector_type, operator_type>::vector_space_type vector_space_type;
    typedef InnerSolver<Scalar, multivector_type, operator_type> inner_solver_type;

    /// \brief Constructor.
    ///
    /// \param solver [in/out] The actual inner solver implementation.
    ThyraInnerSolver (const Teuchos::RCP<inner_solver_type>& solver) :
      solver_ (solver)
    {}

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

    //! A persistent view of the domain vector space of the operator
    Teuchos::RCP<const vector_space_type> domain() const {
      return solver_->getDomain();
    }
    //! A persistent view of the range vector space of the operator
    Teuchos::RCP<const vector_space_type> range() const {
      return solver_->getRange();
    }

  protected:
    //! Whether apply() with the given transpose option is supported.
    bool opSupportedImpl (Thyra::EOpTransp M_trans) const {
      // Applying the (conjugate) transpose is not supported.
      return (M_trans == Thyra::NOTRANS);
    }

    /// \brief Compute Y := alpha*solver(Y,X) + beta*Y.
    ///
    /// \note The contents of Y on input may be relevant, depending on
    ///   the inner solver implementation.  For example, Y on input
    ///   may be treated as the initial guess of an iterative solver.
    void 
    applyImpl(Thyra::EOpTransp M_trans,
	      const multivector_type& X,
	      const Teuchos::Ptr<multivector_type>& Y,
	      const Scalar alpha,
	      const Scalar beta) const
    {
      using Teuchos::RCP;
      using Teuchos::rcpFromRef;
      typedef multivector_type MV;
      typedef Teuchos::ScalarTraits<Scalar> STS;

      TEST_FOR_EXCEPTION(! opSupportedImpl(M_trans), std::invalid_argument,
			 "ThyraInnerSolver only supports applying the operator"
			 " itself, not its transpose or conjugate transpose.");
      TEST_FOR_EXCEPTION(Y.is_null(), std::invalid_argument,
			 "ThyraInnerSolver does not support callign apply() "
			 "with Y null.");
      const Scalar zero = STS::zero();
      const Scalar one = STS::one();

      // "X" is the right-hand side in this case, and "Y" is the
      // "left-hand side."
      RCP<const MV> X_ptr = rcpFromRef (X);
      RCP<MV> Y_ptr = rcpFromRef (*Y);

      // Compute Y := alpha*solver(Y,X) + beta*Y.
      //
      // There are special cases for beta==0 or alpha==0; these
      // special cases ignore NaNs in Y or in the result of
      // solver(Y,X).  Note also that the inner solver is not invoked
      // at all if alpha==0.
      if (beta == 0)
	{
	  if (alpha != zero)
	    solver_->solve (Y_ptr, X_ptr);
	  if (alpha != one)
	    Thyra::scale (alpha, Y_ptr.ptr());
	}
      else if (alpha == 0)
	Thyra::scale (beta, Y_ptr.ptr());
      else 
	{ // The solver overwrites Y with the result of the solve, so
	  // we have to make a copy of Y first before invoking the
	  // solver.
	  RCP<MV> Y_copy = Y_ptr->clone_mv ();
	  solver_->solve (Y_copy, X_ptr);
	  // Y := alpha*Y_copy + beta*Y.
	  if (beta != one)
	    Thyra::scale (beta, Y_ptr.ptr());
	  Thyra::update (alpha, *Y_copy, Y_ptr.ptr());
	}
    }

  private:
    //! Default construction is not allowed.
    ThyraInnerSolver ();

    //! The inner solver implementation.
    Teuchos::RCP<inner_solver_type> solver_;
  };

  
  /// \brief Specialization of makeInnerSolverOperator() for Thyra objects.
  ///
  /// Take an InnerSolver instance, and wrap it in an implementation
  /// of the Thyra::LinearOpBase interface.  That way you can use it
  /// alongside any other implementation of the Thyra::LinearOpBase
  /// interface.
  template<class Scalar>
  class InnerSolverTraits<Scalar, Thyra::MultiVectorBase<Scalar>, Thyra::LinearOpBase<Scalar> > {
  public:
    typedef Scalar scalar_type;
    typedef Thyra::MultiVectorBase<Scalar> multivector_type;
    typedef Thyra::LinearOpBase<Scalar> operator_type;
    typedef InnerSolver<scalar_type, multivector_type, operator_type> inner_solver_type;
    typedef ThyraInnerSolver<Scalar> wrapper_type;

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
// end of file BELOS_THYRA_ADAPTER_HPP
