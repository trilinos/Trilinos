// @HEADER
// *****************************************************************************
//                 Anasazi: Block Eigensolvers Package
//
// Copyright 2004 NTESS and the Anasazi contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ANASAZI_MULTI_VEC_HPP
#define ANASAZI_MULTI_VEC_HPP

/// \file AnasaziMultiVec.hpp
/// \brief Interface for multivectors used by Anasazi' linear solvers.
///
/// We provide two options for letting Anasazi' linear solvers use
/// arbitrary multivector types.  One is via compile-time
/// polymorphism, by specializing Anasazi::MultiVecTraits.  The other
/// is via run-time polymorphism, by implementing Anasazi::MultiVec
/// (the interface defined in this header file).  Anasazi ultimately
/// only uses Anasazi::MultiVecTraits (it uses Anasazi::MultiVec via a
/// specialization of Anasazi::MultiVecTraits for Anasazi::MultiVec),
/// so the preferred way to tell Anasazi how to use your multivector
/// class is via an Anasazi:: MultiVecTraits specialization.  However,
/// some users find a run-time polymorphic interface useful, so we
/// provide it as a service to them.

#include "AnasaziConfigDefs.hpp"
#include "AnasaziMultiVecTraits.hpp"

namespace Anasazi {

/// \class MultiVec
/// \brief Interface for multivectors used by Anasazi's linear solvers.
/// \author Ulrich Hetmaniuk, Rich Lehoucq, and Heidi Thornquist
///
/// \tparam ScalarType The type of entries of the multivector.
///
/// Anasazi accesses multivectors through a traits interface called
/// MultiVecTraits.  If you want to use Anasazi with your own
/// multivector class MV, you may either specialize MultiVecTraits for
/// MV, or you may wrap MV in your own class that implements MultiVec.
/// Specializing MultiVecTraits works via compile-time polymorphism,
/// whereas implementing the MultiVec interface works via run-time
/// polymorphism.  You may pick whichever option you like.  However,
/// specializing MultiVecTraits is the preferred method.  This is
/// because Anasazi's linear solvers always use a specialization of
/// MultiVecTraits to access multivector operations.  They only use
/// MultiVec through a specialization of the MultiVecTraits traits
/// class, which is implemented below in this header file.
///
/// If you want your multivector class (or a wrapper thereof) to
/// implement the MultiVec interface, you should inherit from
/// MultiVec<ScalarType>, where ScalarType is the type of entries in
/// the multivector.  For example, a multivector with entries of type
/// double would inherit from MultiVec<double>.
template <class ScalarType>
class MultiVec {
public:
  //! @name Constructor/Destructor
  //@{ 
  //! Default constructor.
  MultiVec() {}

  //! Destructor (virtual for memory safety of derived classes).
  virtual ~MultiVec () {}

  //@}
  //! @name Creation methods
  //@{ 

  /// \brief Create a new MultiVec with \c numvecs columns.
  /// \return Pointer to the new multivector with uninitialized values.
  virtual MultiVec<ScalarType> * Clone ( const int numvecs ) const = 0;

  /// \brief Create a new MultiVec and copy contents of \c *this into it (deep copy).
  /// \return Pointer to the new multivector	
  virtual MultiVec<ScalarType> * CloneCopy () const = 0;
	
  /*! \brief Creates a new Anasazi::MultiVec and copies the selected contents of \c *this 
    into the new vector (deep copy).  The copied 
    vectors from \c *this are indicated by the \c index.size() indices in \c index.

    \return Pointer to the new multivector	
  */
  virtual MultiVec<ScalarType> * CloneCopy ( const std::vector<int>& index ) const = 0;

  /*! \brief Creates a new Anasazi::MultiVec that shares the selected contents of \c *this.
    The index of the \c numvecs vectors shallow copied from \c *this are indicated by the
    indices given in \c index.

    \return Pointer to the new multivector	
  */
  virtual MultiVec<ScalarType> * CloneViewNonConst ( const std::vector<int>& index ) = 0;
	
  /*! \brief Creates a new Anasazi::MultiVec that shares the selected contents of \c *this.
    The index of the \c numvecs vectors shallow copied from \c *this are indicated by the
    indices given in \c index.

    \return Pointer to the new multivector	
  */
  virtual const MultiVec<ScalarType> * CloneView ( const std::vector<int>& index ) const = 0;

  //@}
  //! @name Dimension information methods	
  //@{ 

  //! The number of rows in the multivector.
  virtual ptrdiff_t GetGlobalLength () const = 0;

  //! The number of vectors (i.e., columns) in the multivector.
  virtual int GetNumberVecs () const = 0;

  //@}
  //! @name Update methods
  //@{ 

  //! Update \c *this with \c alpha * \c A * \c B + \c beta * (\c *this).
  virtual void 
  MvTimesMatAddMv (ScalarType alpha, 
		   const MultiVec<ScalarType>& A, 
		   const Teuchos::SerialDenseMatrix<int,ScalarType>& B, ScalarType beta) = 0;

  //! Replace \c *this with \c alpha * \c A + \c beta * \c B.
  virtual void MvAddMv ( ScalarType alpha, const MultiVec<ScalarType>& A, ScalarType beta, const MultiVec<ScalarType>& B ) = 0;

  //! Scale each element of the vectors in \c *this with \c alpha.
  virtual void MvScale ( ScalarType alpha ) = 0;

  //! Scale each element of the <tt>i</tt>-th vector in \c *this with <tt>alpha[i]</tt>.
  virtual void MvScale ( const std::vector<ScalarType>& alpha ) = 0;

  /*! \brief Compute a dense matrix \c B through the matrix-matrix multiply 
    \c alpha * \c A^T * (\c *this).
  */
  virtual void MvTransMv ( ScalarType alpha, const MultiVec<ScalarType>& A, Teuchos::SerialDenseMatrix<int,ScalarType>& B
#ifdef HAVE_ANASAZI_EXPERIMENTAL
			   , ConjType conj = Anasazi::CONJ
#endif				 
			   ) const = 0;

  /// \brief Compute the dot product of each column of *this with the corresponding column of A.
  ///
  /// Compute a vector \c b whose entries are the individual
  /// dot-products.  That is, <tt>b[i] = A[i]^H * (*this)[i]</tt>
  /// where <tt>A[i]</tt> is the i-th column of A.
  virtual void MvDot ( const MultiVec<ScalarType>& A, std::vector<ScalarType> & b 
#ifdef HAVE_ANASAZI_EXPERIMENTAL
		       , ConjType conj = Anasazi::CONJ
#endif
		       ) const = 0;
  
  //@}
  //! @name Norm method
  //@{ 

  /// \brief Compute the 2-norm of each vector in \c *this.  
  ///
  /// \param normvec [out] On output, normvec[i] holds the 2-norm of the
  ///   \c i-th vector of \c *this.
  virtual void MvNorm ( std::vector<typename Teuchos::ScalarTraits<ScalarType>::magnitudeType> & normvec ) const = 0;

  //@}
  //! @name Initialization methods
  //@{ 

  /// \brief Copy the vectors in \c A to a set of vectors in \c *this.  
  ///
  /// The \c numvecs vectors in \c A are copied to a subset of vectors
  /// in \c *this indicated by the indices given in \c index.
  virtual void SetBlock ( const MultiVec<ScalarType>& A, const std::vector<int>& index ) = 0;
	
  //! Fill all the vectors in \c *this with random numbers.
  virtual void MvRandom () = 0;

  //! Replace each element of the vectors in \c *this with \c alpha.
  virtual void MvInit ( ScalarType alpha ) = 0;

  //@}
  //! @name Print method
  //@{ 

  //! Print \c *this multivector to the \c os output stream.
  virtual void MvPrint ( std::ostream& os ) const = 0;
  //@}

#ifdef HAVE_ANASAZI_TSQR
  //! @name TSQR-related methods
  //@{ 

  /// \brief Compute the QR factorization *this = QR, using TSQR.
  ///
  /// The *this multivector on input is the multivector A to factor.
  /// It is overwritten with garbage on output.
  ///
  /// \param Q [out] On input: a multivector with the same number of
  ///   rows and columns as A (the *this multivector).  Its contents
  ///   are overwritten on output with the (explicitly stored) Q
  ///   factor in the QR factorization of A.
  ///
  /// \param R [out] On output: the R factor in the QR factorization
  ///   of the (input) multivector A.
  ///
  /// \param forceNonnegativeDiagonal [in] If true, then (if
  ///   necessary) do extra work (modifying both the Q and R
  ///   factors) in order to force the R factor to have a
  ///   nonnegative diagonal.
  ///
  /// For syntax's sake, we provide a default implementation of this
  /// method that throws std::logic_error.  You should implement this
  /// method if you intend to use TsqrOrthoManager or
  /// TsqrMatOrthoManager with your subclass of MultiVec.
  virtual void 
  factorExplicit (MultiVec<ScalarType>& Q, 
		  Teuchos::SerialDenseMatrix<int, ScalarType>& R,
		  const bool forceNonnegativeDiagonal=false)
  {
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "The Anasazi::MultiVec<" 
      << Teuchos::TypeNameTraits<ScalarType>::name() << "> subclass which you "
      "are using does not implement the TSQR-related method factorExplicit().");
  }

  /// \brief Use result of factorExplicit() to compute rank-revealing decomposition.
  ///
  /// When calling this method, the *this multivector should be the Q
  /// factor output of factorExplicit().  Using that Q factor and the
  /// R factor from factorExplicit(), compute the singular value
  /// decomposition (SVD) of R (\f$R = U \Sigma V^*\f$).  If R is full
  /// rank (with respect to the given relative tolerance tol), don't
  /// change Q (= *this) or R.  Otherwise, compute \f$Q := Q \cdot
  /// U\f$ and \f$R := \Sigma V^*\f$ in place (the latter may be no
  /// longer upper triangular).
  ///
  /// The *this multivector on input must be the explicit Q factor
  /// output of a previous call to factorExplicit().  On output: On
  /// output: If R is of full numerical rank with respect to the
  /// tolerance tol, Q is unmodified.  Otherwise, Q is updated so that
  /// the first rank columns of Q are a basis for the column space of
  /// A (the original matrix whose QR factorization was computed by
  /// factorExplicit()).  The remaining columns of Q are a basis for
  /// the null space of A.
  ///
  /// \param R [in/out] On input: N by N upper triangular matrix with
  ///   leading dimension LDR >= N.  On output: if input is full rank,
  ///   R is unchanged on output.  Otherwise, if \f$R = U \Sigma
  ///   V^*\f$ is the SVD of R, on output R is overwritten with
  ///   \f$\Sigma \cdot V^*\f$.  This is also an N by N matrix, but
  ///   may not necessarily be upper triangular.
  ///
  /// \param tol [in] Relative tolerance for computing the numerical
  ///   rank of the matrix R.
  ///
  /// For syntax's sake, we provide a default implementation of this
  /// method that throws std::logic_error.  You should implement this
  /// method if you intend to use TsqrOrthoManager or
  /// TsqrMatOrthoManager with your subclass of MultiVec.
  virtual int
  revealRank (Teuchos::SerialDenseMatrix<int, ScalarType>& R,
	      const typename Teuchos::ScalarTraits<ScalarType>::magnitudeType& tol)
  {
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "The Anasazi::MultiVec<" 
      << Teuchos::TypeNameTraits<ScalarType>::name() << "> subclass which you "
      "are using does not implement the TSQR-related method revealRank().");
  }

  //@}
#endif // HAVE_ANASAZI_TSQR
};

namespace details {
/// \class MultiVecTsqrAdapter
/// \brief TSQR adapter for MultiVec.
///
/// TSQR (Tall Skinny QR factorization) is an orthogonalization
/// kernel that is as accurate as Householder QR, yet requires only
/// \f$2 \log P\f$ messages between $P$ MPI processes, independently
/// of the number of columns in the multivector.  
///
/// TSQR works independently of the particular multivector
/// implementation, and interfaces to the latter via an adapter
/// class.  Each multivector type MV needs its own adapter class.
/// The specialization of MultiVecTraits for MV refers to its
/// corresponding adapter class as its \c tsqr_adaptor_type [sic;
/// sorry about the lack of standard spelling of "adapter"] typedef.
///
/// This class is the TSQR adapter for MultiVec.  It merely calls 
/// MultiVec's corresponding methods for TSQR functionality.
template<class ScalarType>
class MultiVecTsqrAdapter {
public:
  typedef MultiVec<ScalarType> MV;
  typedef ScalarType scalar_type; 
  typedef int ordinal_type; // This doesn't matter either
  typedef int node_type; // Nor does this
  typedef Teuchos::SerialDenseMatrix<ordinal_type, scalar_type> dense_matrix_type;
  typedef typename Teuchos::ScalarTraits<scalar_type>::magnitudeType magnitude_type;
  
  //! Compute QR factorization A = QR, using TSQR.
  void
  factorExplicit (MV& A,
		  MV& Q,
		  dense_matrix_type& R,
		  const bool forceNonnegativeDiagonal=false)
  {
    A.factorExplicit (Q, R, forceNonnegativeDiagonal);
  }

  //! Compute rank-revealing decomposition using results of factorExplicit().
  int
  revealRank (MV& Q,
	      dense_matrix_type& R,
	      const magnitude_type& tol)
  {
    return Q.revealRank (R, tol);
  }
};
} // namespace details

  /// \brief Specialization of MultiVecTraits for Belos::MultiVec.
  ///
  /// Anasazi interfaces to every multivector implementation through a
  /// specialization of MultiVecTraits.  Thus, we provide a
  /// specialization of MultiVecTraits for the MultiVec run-time
  /// polymorphic interface above.
  ///
  /// \tparam ScalarType The type of entries in the multivector; the
  ///   template parameter of MultiVec.
  template<class ScalarType>
  class MultiVecTraits<ScalarType,MultiVec<ScalarType> > {
  public:
    //! @name Creation methods
    //@{ 

    /// \brief Create a new empty \c MultiVec containing \c numvecs columns.
    /// \return Reference-counted pointer to the new \c MultiVec.
    static Teuchos::RCP<MultiVec<ScalarType> > 
    Clone (const MultiVec<ScalarType>& mv, const int numvecs) {
      return Teuchos::rcp (const_cast<MultiVec<ScalarType>&> (mv).Clone (numvecs)); 
    }

    /*! \brief Creates a new \c Anasazi::MultiVec and copies contents of \c mv into the new vector (deep copy).
      
      \return Reference-counted pointer to the new \c Anasazi::MultiVec.
    */
    static Teuchos::RCP<MultiVec<ScalarType> > CloneCopy( const MultiVec<ScalarType>& mv )
    { return Teuchos::rcp( const_cast<MultiVec<ScalarType>&>(mv).CloneCopy() ); }

    /*! \brief Creates a new \c Anasazi::MultiVec and copies the selected contents of \c mv into the new vector (deep copy).  

      The copied vectors from \c mv are indicated by the \c index.size() indices in \c index.      
      \return Reference-counted pointer to the new \c Anasazi::MultiVec.
    */
    static Teuchos::RCP<MultiVec<ScalarType> > CloneCopy( const MultiVec<ScalarType>& mv, const std::vector<int>& index )
    { return Teuchos::rcp( const_cast<MultiVec<ScalarType>&>(mv).CloneCopy(index) ); }

    /*! \brief Creates a new \c Anasazi::MultiVec that shares the selected contents of \c mv (shallow copy).

    The index of the \c numvecs vectors shallow copied from \c mv are indicated by the indices given in \c index.
    \return Reference-counted pointer to the new \c Anasazi::MultiVec.
    */    
    static Teuchos::RCP<MultiVec<ScalarType> > CloneViewNonConst( MultiVec<ScalarType>& mv, const std::vector<int>& index )
    { return Teuchos::rcp( mv.CloneViewNonConst(index) ); }

    /*! \brief Creates a new const \c Anasazi::MultiVec that shares the selected contents of \c mv (shallow copy).

    The index of the \c numvecs vectors shallow copied from \c mv are indicated by the indices given in \c index.
    \return Reference-counted pointer to the new const \c Anasazi::MultiVec.
    */      
    static Teuchos::RCP<const MultiVec<ScalarType> > CloneView( const MultiVec<ScalarType>& mv, const std::vector<int>& index )
    { return Teuchos::rcp( const_cast<MultiVec<ScalarType>&>(mv).CloneView(index) ); }

    //@}

    //! @name Attribute methods
    //@{ 

    //! Obtain the vector length of \c mv.
    static ptrdiff_t GetGlobalLength( const MultiVec<ScalarType>& mv )
    { return mv.GetGlobalLength(); }

    //! Obtain the number of vectors in \c mv
    static int GetNumberVecs( const MultiVec<ScalarType>& mv )
    { return mv.GetNumberVecs(); }

    //@}

    //! @name Update methods
    //@{ 

    /*! \brief Update \c mv with \f$ \alpha AB + \beta mv \f$.
     */
    static void MvTimesMatAddMv( ScalarType alpha, const MultiVec<ScalarType>& A, 
				 const Teuchos::SerialDenseMatrix<int,ScalarType>& B, 
				 ScalarType beta, MultiVec<ScalarType>& mv )
    { mv.MvTimesMatAddMv(alpha, A, B, beta); }

    /*! \brief Replace \c mv with \f$\alpha A + \beta B\f$.
     */
    static void MvAddMv( ScalarType alpha, const MultiVec<ScalarType>& A, ScalarType beta, const MultiVec<ScalarType>& B, MultiVec<ScalarType>& mv )
    { mv.MvAddMv(alpha, A, beta, B); }

    /*! \brief Compute a dense matrix \c B through the matrix-matrix multiply \f$ \alpha A^Tmv \f$.
     */
    static void MvTransMv( ScalarType alpha, const MultiVec<ScalarType>& A, const MultiVec<ScalarType>& mv, Teuchos::SerialDenseMatrix<int,ScalarType>& B
#ifdef HAVE_ANASAZI_EXPERIMENTAL
			   , ConjType conj = Anasazi::CONJ
#endif
			   )
    { mv.MvTransMv(alpha, A, B
#ifdef HAVE_ANASAZI_EXPERIMENTAL
		   , conj
#endif
		   ); }
    
    /*! \brief Compute a vector \c b where the components are the individual dot-products of the \c i-th columns of \c A and \c mv, i.e.\f$b[i] = A[i]^H mv[i]\f$.
     */
    static void MvDot( const MultiVec<ScalarType>& mv, const MultiVec<ScalarType>& A, std::vector<ScalarType> & b
#ifdef HAVE_ANASAZI_EXPERIMENTAL
		       , ConjType conj = Anasazi::CONJ
#endif
		       )
    { mv.MvDot( A, b
#ifdef HAVE_ANASAZI_EXPERIMENTAL
		, conj
#endif
		); }

    //! Scale each element of the vectors in \c *this with \c alpha.
    static void MvScale ( MultiVec<ScalarType>& mv, ScalarType alpha )
    { mv.MvScale( alpha ); }
    
    //! Scale each element of the \c i-th vector in \c *this with \c alpha[i].
    static void MvScale ( MultiVec<ScalarType>& mv, const std::vector<ScalarType>& alpha )
    { mv.MvScale( alpha ); }
    
    //@}
    //! @name Norm method
    //@{ 

    /*! \brief Compute the 2-norm of each individual vector of \c mv.  
      Upon return, \c normvec[i] holds the value of \f$||mv_i||_2\f$, the \c i-th column of \c mv.
    */
    static void MvNorm( const MultiVec<ScalarType>& mv, std::vector<typename Teuchos::ScalarTraits<ScalarType>::magnitudeType> & normvec )
    { mv.MvNorm(normvec); }

    //@}
    //! @name Initialization methods
    //@{ 
    /*! \brief Copy the vectors in \c A to a set of vectors in \c mv indicated by the indices given in \c index.

    The \c numvecs vectors in \c A are copied to a subset of vectors in \c mv indicated by the indices given in \c index,
    i.e.<tt> mv[index[i]] = A[i]</tt>.
    */
    static void SetBlock( const MultiVec<ScalarType>& A, const std::vector<int>& index, MultiVec<ScalarType>& mv )
    { mv.SetBlock(A, index); }

    /*! \brief Replace the vectors in \c mv with random vectors.
     */
    static void MvRandom( MultiVec<ScalarType>& mv )
    { mv.MvRandom(); }

    /*! \brief Replace each element of the vectors in \c mv with \c alpha.
     */
    static void MvInit( MultiVec<ScalarType>& mv, ScalarType alpha = Teuchos::ScalarTraits<ScalarType>::zero() )
    { mv.MvInit(alpha); }

    //@}
    //! @name Print method
    //@{ 

    //! Print the \c mv multi-vector to the \c os output stream.
    static void MvPrint( const MultiVec<ScalarType>& mv, std::ostream& os )
    { mv.MvPrint(os); }

    //@}

#ifdef HAVE_ANASAZI_TSQR
    /// \typedef tsqr_adaptor_type
    /// \brief TSQR adapter for MultiVec.
    ///
    /// Our TSQR adapter for MultiVec calls MultiVec's virtual
    /// methods.  If you want to use TSQR with your MultiVec subclass,
    /// you must implement these methods yourself, as the default
    /// implementations throw std::logic_error.
    typedef details::MultiVecTsqrAdapter<ScalarType> tsqr_adaptor_type;
#endif // HAVE_ANASAZI_TSQR
  };

} // namespace Anasazi

#endif

// end of file AnasaziMultiVec.hpp
