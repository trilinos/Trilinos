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
//
#ifndef ANASAZI_MULTI_VEC_TRAITS_HPP
#define ANASAZI_MULTI_VEC_TRAITS_HPP

/// \file AnasaziMultiVecTraits.hpp
/// \brief Declaration of basic traits for the multivector type
///
/// Anasazi::MultiVecTraits declares basic traits for the multivector
/// type MV used in Anasazi's orthogonalizations and solvers.  A
/// specialization of MultiVecTraits that defines all the traits must
/// be made for each specific multivector type.  Here, we only provide
/// default definitions that fail at compile time if no specialization
/// of MultiVecTraits exists for the given combination of scalar type
/// (ScalarType) and multivector type (MV).

#include "AnasaziTypes.hpp"
#include "Teuchos_Range1D.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_SerialDenseMatrix.hpp"

namespace Anasazi {

  /*! \brief This is the default struct used by MultiVecTraits<ScalarType, MV> class to produce a
      compile time error when the specialization does not exist for multivector type <tt>MV</tt>.
  */
  template< class ScalarType, class MV >
  struct UndefinedMultiVecTraits
  {
    //! This function should not compile if there is an attempt to instantiate!
    /*! \note Any attempt to compile this function results in a compile time error.  This means
      that the template specialization of Anasazi::MultiVecTraits class for type <tt>MV</tt> does 
      not exist, or is not complete.
    */
    static inline ScalarType notDefined() { return MV::this_type_is_missing_a_specialization(); };
  };


  /// \brief Traits class which defines basic operations on multivectors.
  /// \ingroup anasazi_opvec_interfaces
  ///
  /// \tparam ScalarType The type of the entries in the multivectors.
  /// \tparam MV The type of the multivectors themselves.
  ///
  /// This traits class tells Anasazi's solvers how to perform
  /// multivector operations for the multivector type MV.  These
  /// operations include creating copies or views, finding the number
  /// of rows or columns (i.e., vectors) in a given multivector, and
  /// computing inner products, norms, and vector sums.  (Anasazi's
  /// solvers use the OperatorTraits traits class to apply operators
  /// to multivectors.)
  ///
  /// Anasazi gives users two different ways to tell its solvers how
  /// to compute with multivectors of a given type MV.  The first and
  /// preferred way is for users to specialize MultiVecTraits, this
  /// traits class, for their given MV type.  Anasazi provides
  /// specializations for MV = Epetra_MultiVector,
  /// Tpetra::MultiVector, and Thyra::MultiVectorBase.  The second way
  /// is for users to make their multivector type (or a wrapper
  /// thereof) inherit from MultiVec.  This works because Anasazi
  /// provides a specialization of MultiVecTraits for MultiVec.
  /// Specializing MultiVecTraits is more flexible because it does not
  /// require a multivector type to inherit from MultiVec; this is
  /// possible even if you do not have control over the interface of a
  /// class.
  ///
  /// If you have a different multivector type MV that you would like
  /// to use with Anasazi, and if that type does not inherit from
  /// MultiVec, then you must implement a specialization of
  /// MultiVecTraits for MV.  Otherwise, this traits class will report
  /// a compile-time error (relating to UndefinedMultiVecTraits).
  template<class ScalarType, class MV>
  class MultiVecTraits 
  {
  public:

    //! @name Creation methods
    //@{

    /*! \brief Creates a new empty \c MV containing \c numvecs columns.

    \return Reference-counted pointer to the new multivector of type \c MV.
    */
    static Teuchos::RCP<MV> Clone( const MV& mv, const int numvecs )
    { UndefinedMultiVecTraits<ScalarType, MV>::notDefined(); return Teuchos::null; }     

    /*! \brief Creates a new \c MV and copies contents of \c mv into the new vector (deep copy).

      \return Reference-counted pointer to the new multivector of type \c MV.
    */
    static Teuchos::RCP<MV> CloneCopy( const MV& mv )
    { UndefinedMultiVecTraits<ScalarType, MV>::notDefined(); return Teuchos::null; }     

    /*! \brief Creates a new \c MV and copies the selected contents of \c mv into the new vector (deep copy).

      The copied vectors from \c mv are indicated by the \c index.size() indices in \c index.
      \return Reference-counted pointer to the new multivector of type \c MV.
    */
    static Teuchos::RCP<MV> CloneCopy( const MV& mv, const std::vector<int>& index )
    { UndefinedMultiVecTraits<ScalarType, MV>::notDefined(); return Teuchos::null; }     

    /// \brief Deep copy of specified columns of mv
    ///
    /// Create a new MV, and copy (deep copy) the columns of mv
    /// specified by the given inclusive index range into the new
    /// multivector.
    ///
    /// \param mv [in] Multivector to copy
    /// \param index [in] Inclusive index range of columns of mv
    /// \return Reference-counted pointer to the new multivector of type \c MV.
    static Teuchos::RCP<MV> CloneCopy( const MV& mv, const Teuchos::Range1D& index )
    { UndefinedMultiVecTraits<ScalarType, MV>::notDefined(); return Teuchos::null; }     

    /*! \brief Creates a new \c MV that shares the selected contents of \c mv (shallow copy).

    The index of the \c numvecs vectors shallow copied from \c mv are indicated by the indices given in \c index.
    \return Reference-counted pointer to the new multivector of type \c MV.
    */
    static Teuchos::RCP<MV> CloneViewNonConst( MV& mv, const std::vector<int>& index )
    { UndefinedMultiVecTraits<ScalarType, MV>::notDefined(); return Teuchos::null; }     

    /// \brief Non-const view of specified columns of mv
    ///
    /// Return a non-const view of the columns of mv specified by the
    /// given inclusive index range.
    ///
    /// \param mv [in] Multivector to view (shallow non-const copy)
    /// \param index [in] Inclusive index range of columns of mv
    /// \return Reference-counted pointer to the non-const view of specified columns of mv
    static Teuchos::RCP<MV> CloneViewNonConst( MV& mv, const Teuchos::Range1D& index )
    { UndefinedMultiVecTraits<ScalarType, MV>::notDefined(); return Teuchos::null; }     

    /*! \brief Creates a new const \c MV that shares the selected contents of \c mv (shallow copy).

    The index of the \c numvecs vectors shallow copied from \c mv are indicated by the indices given in \c index.
    \return Reference-counted pointer to the new const multivector of type \c MV.
    */
    static Teuchos::RCP<const MV> CloneView( const MV& mv, const std::vector<int>& index )
    { UndefinedMultiVecTraits<ScalarType, MV>::notDefined(); return Teuchos::null; }     

    /// \brief Const view of specified columns of mv
    ///
    /// Return a const view of the columns of mv specified by the
    /// given inclusive index range.
    ///
    /// \param mv [in] Multivector to view (shallow const copy)
    /// \param index [in] Inclusive index range of columns of mv
    /// \return Reference-counted pointer to the const view of specified columns of mv
    static Teuchos::RCP<MV> CloneView( MV& mv, const Teuchos::Range1D& index )
    { UndefinedMultiVecTraits<ScalarType, MV>::notDefined(); return Teuchos::null; }     

    //@}

    //! @name Attribute methods
    //@{

    //! Obtain the vector length of \c mv.
    static int GetVecLength( const MV& mv )
    { UndefinedMultiVecTraits<ScalarType, MV>::notDefined(); return 0; }     

    //! Obtain the number of vectors in \c mv
    static int GetNumberVecs( const MV& mv )
    { UndefinedMultiVecTraits<ScalarType, MV>::notDefined(); return 0; }     

    //@}

    //! @name Update methods
    //@{

    /*! \brief Update \c mv with \f$ \alpha AB + \beta mv \f$.
     */
    static void MvTimesMatAddMv( const ScalarType alpha, const MV& A, 
                                 const Teuchos::SerialDenseMatrix<int,ScalarType>& B, 
                                 const ScalarType beta, MV& mv )
    { UndefinedMultiVecTraits<ScalarType, MV>::notDefined(); }     

    /*! \brief Replace \c mv with \f$\alpha A + \beta B\f$.
     */
    static void MvAddMv( const ScalarType alpha, const MV& A, const ScalarType beta, const MV& B, MV& mv )
    { UndefinedMultiVecTraits<ScalarType, MV>::notDefined(); }     

    /*! \brief Scale each element of the vectors in \c mv with \c alpha.
     */
    static void MvScale ( MV& mv, const ScalarType alpha )
    { UndefinedMultiVecTraits<ScalarType, MV>::notDefined(); }
    
    /*! \brief Scale each element of the \c i-th vector in \c mv with \c alpha[i].
     */
    static void MvScale ( MV& mv, const std::vector<ScalarType>& alpha )
    { UndefinedMultiVecTraits<ScalarType, MV>::notDefined(); }

    /*! \brief Compute a dense matrix \c B through the matrix-matrix multiply \f$ \alpha A^Hmv \f$.
    */
    static void MvTransMv( const ScalarType alpha, const MV& A, const MV& mv, Teuchos::SerialDenseMatrix<int,ScalarType>& B)
    { UndefinedMultiVecTraits<ScalarType, MV>::notDefined(); }     

    /*! \brief Compute a vector \c b where the components are the individual dot-products of the \c i-th columns of \c A and \c mv, i.e.\f$b[i] = A[i]^Hmv[i]\f$.
     */
    static void MvDot ( const MV& mv, const MV& A, std::vector<ScalarType> &b) 
    { UndefinedMultiVecTraits<ScalarType, MV>::notDefined(); }     

    //@}
    //! @name Norm method
    //@{

    /*! \brief Compute the 2-norm of each individual vector of \c mv.
      Upon return, \c normvec[i] holds the value of \f$||mv_i||_2\f$, the \c i-th column of \c mv.
    */
    static void MvNorm( const MV& mv, std::vector<typename Teuchos::ScalarTraits<ScalarType>::magnitudeType> &normvec )
    { UndefinedMultiVecTraits<ScalarType, MV>::notDefined(); }     

    //@}

    //! @name Initialization methods
    //@{
    /*! \brief Copy the vectors in \c A to a set of vectors in \c mv indicated by the indices given in \c index.

    The \c numvecs vectors in \c A are copied to a subset of vectors in \c mv indicated by the indices given in \c index,
    i.e.<tt> mv[index[i]] = A[i]</tt>.
    */
    static void SetBlock( const MV& A, const std::vector<int>& index, MV& mv )
    { UndefinedMultiVecTraits<ScalarType, MV>::notDefined(); }     

    /// \brief Deep copy of A into specified columns of mv
    ///
    /// (Deeply) copy the first <tt>index.size()</tt> columns of \c A
    /// into the columns of \c mv specified by the given index range.
    ///
    /// Postcondition: <tt>mv[i] = A[i - index.lbound()]</tt>
    /// for all <tt>i</tt> in <tt>[index.lbound(), index.ubound()]</tt>
    ///
    /// \param A [in] Source multivector
    /// \param index [in] Inclusive index range of columns of mv;
    ///   index set of the target
    /// \param mv [out] Target multivector
    static void SetBlock( const MV& A, const Teuchos::Range1D& index, MV& mv )
    { UndefinedMultiVecTraits<ScalarType, MV>::notDefined(); }     

    /// \brief mv := A
    /// 
    /// Assign (deep copy) A into mv.
    static void Assign( const MV& A, MV& mv )
    { UndefinedMultiVecTraits<ScalarType, MV>::notDefined(); }     

    /*! \brief Replace the vectors in \c mv with random vectors.
     */
    static void MvRandom( MV& mv )
    { UndefinedMultiVecTraits<ScalarType, MV>::notDefined(); }     

    /*! \brief Replace each element of the vectors in \c mv with \c alpha.
     */
    static void MvInit( MV& mv, const ScalarType alpha = Teuchos::ScalarTraits<ScalarType>::zero() )
    { UndefinedMultiVecTraits<ScalarType, MV>::notDefined(); }     

    //@}

    //! @name Print method
    //@{

    /*! \brief Print the \c mv multi-vector to the \c os output stream.
     */
    static void MvPrint( const MV& mv, std::ostream& os )
    { UndefinedMultiVecTraits<ScalarType, MV>::notDefined(); }     

    //@}
  };
  
} // namespace Anasazi

#endif // ANASAZI_MULTI_VEC_TRAITS_HPP
