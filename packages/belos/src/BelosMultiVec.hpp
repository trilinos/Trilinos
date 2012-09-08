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

#ifndef BELOS_MULTI_VEC_HPP
#define BELOS_MULTI_VEC_HPP

/// \file BelosMultiVec.hpp
/// \brief Interface for multivectors used by Belos' linear solvers.
///
/// We provide two options for letting Belos' linear solvers use
/// arbitrary multivector types.  One is via compile-time
/// polymorphism, by specializing MultiVecTraits.  The other is via
/// run-time polymorphism, by implementing MultiVec (the interface
/// defined in this header file).  Belos ultimately only uses
/// MultiVecTraits (it uses MultiVec via a specialization of
/// MultiVecTraits for MultiVec), so the preferred way to tell Belos
/// how to use your multivector class is via a MultiVecTraits
/// specialization.  However, some users find a run-time polymorphic
/// interface useful, so we provide it as a service to them.

#include "BelosMultiVecTraits.hpp"
#include "BelosTypes.hpp"
#include "BelosConfigDefs.hpp"

namespace Belos {

/// \class MultiVec
/// \brief Interface for multivectors used by Belos' linear solvers.
/// \author Michael Heroux, Rich Lehoucq, and Heidi Thornquist
///
/// \tparam ScalarType The type of entries of the multivector.
///
/// Belos accesses multivectors through a traits interface called
/// MultiVecTraits.  If you want to use Belos with your own
/// multivector class MV, you may either specialize MultiVecTraits for
/// MV, or you may wrap MV in your own class that implements MultiVec.
/// Specializing MultiVecTraits works via compile-time polymorphism,
/// whereas implementing the MultiVec interface works via run-time
/// polymorphism.  You may pick whichever option you like.  However,
/// specializing MultiVecTraits is the preferred method.  This is
/// because Belos' linear solvers always use a specialization of
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
  MultiVec() {};
  
  //! Destructor (virtual for memory safety of derived classes).
  virtual ~MultiVec () {};
  
  //@}
  //! @name Creation methods for new multivectors
  //@{ 
  
  /// \brief Create a new MultiVec with \c numvecs columns.
  /// \return Pointer to the new multivector with uninitialized values.
  virtual MultiVec<ScalarType> * Clone ( const int numvecs ) const = 0;
  
  /// \brief Create a new MultiVec and copy contents of \c *this into it (deep copy).
  /// \return Pointer to the new multivector	
  virtual MultiVec<ScalarType> * CloneCopy () const = 0;
  
  /*! \brief Creates a new %Belos::MultiVec and copies the selected contents of \c *this 
    into the new multivector (deep copy).  The copied 
    vectors from \c *this are indicated by the \c index.size() indices in \c index.
    
    \return Pointer to the new multivector	
  */
  
  virtual MultiVec<ScalarType> * CloneCopy ( const std::vector<int>& index ) const = 0;
  
  /*! \brief Creates a new %Belos::MultiVec that shares the selected contents of \c *this.
    The index of the \c numvecs vectors copied from \c *this are indicated by the
    indices given in \c index.
    
    \return Pointer to the new multivector	
  */
  
  virtual MultiVec<ScalarType> * CloneViewNonConst ( const std::vector<int>& index ) = 0;
  
  /*! \brief Creates a new %Belos::MultiVec that shares the selected contents of \c *this.
    The index of the \c numvecs vectors copied from \c *this are indicated by the
    indices given in \c index.
    
    \return Pointer to the new multivector	
  */
  
  virtual const MultiVec<ScalarType> * CloneView ( const std::vector<int>& index ) const = 0;
  //@}
  
  //! @name Dimension information methods	
  //@{ 
  //! Obtain the multivector length of *this multivector block.
  
  virtual int GetVecLength () const = 0;
  
  //! Obtain the number of vectors in *this multivector block.
  
  virtual int GetNumberVecs () const = 0;
  
  //@}
  //! @name Update methods
  //@{ 
  /*! \brief Update \c *this with \c alpha * \c A * \c B + \c beta * (\c *this).
   */
  
  virtual void MvTimesMatAddMv ( const ScalarType alpha, const MultiVec<ScalarType>& A, 
				 const Teuchos::SerialDenseMatrix<int,ScalarType>& B, const ScalarType beta ) = 0;
  
  /*! \brief Replace \c *this with \c alpha * \c A + \c beta * \c B.
   */
  
  virtual void MvAddMv ( const ScalarType alpha, const MultiVec<ScalarType>& A, const ScalarType beta, const MultiVec<ScalarType>& B ) = 0;
  
  /*! \brief Scale each element of the vectors in \c *this with \c alpha.
   */
  
  virtual void MvScale ( const ScalarType alpha ) = 0;
  
  /*! \brief Scale each element of the \c i-th vector in \c *this with \c alpha[i].
   */
  
  virtual void MvScale ( const std::vector<ScalarType>& alpha ) = 0;
  
  /*! \brief Compute a dense matrix \c B through the matrix-matrix multiply 
    \c alpha * \c A^T * (\c *this).
  */
  
  virtual void MvTransMv ( const ScalarType alpha, const MultiVec<ScalarType>& A, Teuchos::SerialDenseMatrix<int,ScalarType>& B) const = 0;
  
  /*! \brief Compute a multivector \c b where the components are the individual dot-products, i.e.\c b[i] = \c A[i]^T*\c this[i] where \c A[i] is the i-th column of A.
   */
  
  virtual void MvDot ( const MultiVec<ScalarType>& A, std::vector<ScalarType>& b ) const = 0;
  
  //@}
  //! @name Norm method
  //@{ 
  
  /*! \brief Compute the 2-norm of each vector of \c *this.  
    Upon return, \c normvec[i] holds the 2-norm of the \c i-th vector of \c *this
  */
  
  virtual void MvNorm ( std::vector<typename Teuchos::ScalarTraits<ScalarType>::magnitudeType>& normvec, NormType type = TwoNorm ) const = 0;
  
  //@}
  //! @name Initialization methods
  //@{ 
  /*! \brief Copy the vectors in \c A to a set of vectors in \c *this.  The \c 
    numvecs vectors in \c A are copied to a subset of vectors in \c *this
    indicated by the indices given in \c index.
  */
  
  virtual void SetBlock ( const MultiVec<ScalarType>& A, const std::vector<int>& index ) = 0;
  
  /*! \brief Replace the vectors in \c *this with random vectors.
   */
  
  virtual void MvRandom () = 0;
  
  /*! \brief Replace each element of the vectors in \c *this with \c alpha.
   */
  
  virtual void MvInit ( const ScalarType alpha ) = 0;
  
  //@}
  //! @name Print method
  //@{ 
  /*! \brief Print the \c *this multivector.
   */
  virtual void MvPrint ( std::ostream& os ) const = 0;
  //@}
};


  ////////////////////////////////////////////////////////////////////
  //
  // Implementation of the Belos::MultiVecTraits for Belos::MultiVec.
  //
  ////////////////////////////////////////////////////////////////////


  /// \brief Specialization of MultiVecTraits for Belos::MultiVec.
  ///
  /// This is a partial specialization of MultiVecTraits for
  /// Belos::MultiVec, which is a generic interface that users may
  /// implement in order to wrap their own multivector
  /// implementations.
  template<class ScalarType>
  class MultiVecTraits<ScalarType,MultiVec<ScalarType> >
  {
  public:

    ///
    static Teuchos::RCP<MultiVec<ScalarType> > Clone( const MultiVec<ScalarType>& mv, const int numvecs )
    { return Teuchos::rcp( const_cast<MultiVec<ScalarType>&>(mv).Clone(numvecs) ); }
    ///
    static Teuchos::RCP<MultiVec<ScalarType> > CloneCopy( const MultiVec<ScalarType>& mv )
    { return Teuchos::rcp( const_cast<MultiVec<ScalarType>&>(mv).CloneCopy() ); }
    ///
    static Teuchos::RCP<MultiVec<ScalarType> > CloneCopy( const MultiVec<ScalarType>& mv, const std::vector<int>& index )
    { return Teuchos::rcp( const_cast<MultiVec<ScalarType>&>(mv).CloneCopy(index) ); }
    ///
    static Teuchos::RCP<MultiVec<ScalarType> > CloneViewNonConst( MultiVec<ScalarType>& mv, const std::vector<int>& index )
    { return Teuchos::rcp( mv.CloneViewNonConst(index) ); }
    ///
    static Teuchos::RCP<const MultiVec<ScalarType> > CloneView( const MultiVec<ScalarType>& mv, const std::vector<int>& index )
    { return Teuchos::rcp( const_cast<MultiVec<ScalarType>&>(mv).CloneView(index) ); }
    ///
    static int GetVecLength( const MultiVec<ScalarType>& mv )
    { return mv.GetVecLength(); }
    ///
    static int GetNumberVecs( const MultiVec<ScalarType>& mv )
    { return mv.GetNumberVecs(); }
    ///
    static void MvTimesMatAddMv( ScalarType alpha, const MultiVec<ScalarType>& A, 
				 const Teuchos::SerialDenseMatrix<int,ScalarType>& B, 
				 ScalarType beta, MultiVec<ScalarType>& mv )
    { mv.MvTimesMatAddMv(alpha, A, B, beta); }
    ///
    static void MvAddMv( ScalarType alpha, const MultiVec<ScalarType>& A, ScalarType beta, const MultiVec<ScalarType>& B, MultiVec<ScalarType>& mv )
    { mv.MvAddMv(alpha, A, beta, B); }
    ///
    static void MvScale ( MultiVec<ScalarType>& mv, const ScalarType alpha )
    { mv.MvScale( alpha ); } 

    static void MvScale ( MultiVec<ScalarType>& mv, const std::vector<ScalarType>& alpha )
    { mv.MvScale(alpha); }
    ///
    static void MvTransMv( const ScalarType alpha, const MultiVec<ScalarType>& A, const MultiVec<ScalarType>& mv, Teuchos::SerialDenseMatrix<int,ScalarType>& B )
    { mv.MvTransMv(alpha, A, B); }
    ///
    static void MvDot( const MultiVec<ScalarType>& mv, const MultiVec<ScalarType>& A, std::vector<ScalarType>& b )
    { mv.MvDot( A, b ); }
    ///
    static void MvNorm( const MultiVec<ScalarType>& mv, std::vector<typename Teuchos::ScalarTraits<ScalarType>::magnitudeType>& normvec, NormType type = TwoNorm )
    { mv.MvNorm(normvec,type); }
    ///
    static void SetBlock( const MultiVec<ScalarType>& A, const std::vector<int>& index, MultiVec<ScalarType>& mv )
    { mv.SetBlock(A, index); }
    ///
    static void MvRandom( MultiVec<ScalarType>& mv )
    { mv.MvRandom(); }
    ///
    static void MvInit( MultiVec<ScalarType>& mv, ScalarType alpha = Teuchos::ScalarTraits<ScalarType>::zero() )
    { mv.MvInit(alpha); }
    ///
    static void MvPrint( const MultiVec<ScalarType>& mv, std::ostream& os )
    { mv.MvPrint(os); }

#ifdef HAVE_BELOS_TSQR
    /// \typedef tsqr_adaptor_type
    /// \brief TsqrAdaptor specialization for the multivector type MV.
    ///
    /// By default, we provide a "stub" implementation.  It has the
    /// right methods and typedefs, but its constructors and methods
    /// all throw std::logic_error.  Later, we will extend MultiVec
    /// with a TSQR interface itself, and write an adapter for
    /// MultiVec that calls MultiVec's methods.  This will allow users
    /// of MultiVec to extend it to implement a TSQR adapter if they
    /// wish.
    typedef Belos::details::StubTsqrAdapter<MultiVec<ScalarType> > tsqr_adaptor_type;
#endif // HAVE_BELOS_TSQR
  };


} // namespace Belos

#endif

// end of file BelosMultiVec.hpp
