// @HEADER
// *****************************************************************************
//                 Belos: Block Linear Solvers Package
//
// Copyright 2004-2016 NTESS and the Belos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/*! \file BelosOrthoManager.hpp
  \brief  Templated virtual class for providing orthogonalization/orthonormalization methods.
*/

#ifndef BELOS_ORTHOMANAGER_HPP
#define BELOS_ORTHOMANAGER_HPP

/*! \class Belos::OrthoManager

  \brief Belos's templated virtual class for providing routines for orthogonalization and
  orthonormzalition of multivectors.

  This class defines concepts of orthogonality through the definition of an
  inner product. It also provides computational routines for orthogonalization.

  A concrete implementation of this class is necessary. The user can create
  their own implementation if those supplied are not suitable for their needs.

  Note: The OrthoManager class now inherits from 
  Teuchos::ParameterListAcceptorDefaultBase. New derived classes must implement
  <tt>setParameterList()</tt> and <tt>getValidParameters()</tt>.

  \author Chris Baker, Teri Barth, and Heidi Thornquist
*/

#include "BelosConfigDefs.hpp"
#include "BelosTypes.hpp"
#include "Teuchos_ScalarTraits.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_SerialDenseMatrix.hpp"
#include "Teuchos_Array.hpp"


namespace Belos {


  //! @name OrthoManager Exceptions
  //@{

  /** \brief Exception thrown to signal error in an orthogonalization manager method.
   */
  class OrthoError : public BelosError
  {public: OrthoError(const std::string& what_arg) : BelosError(what_arg) {}};

  //@}

  template <class ScalarType, class MV>
  class OrthoManager : 
    public Teuchos::ParameterListAcceptorDefaultBase 
  {
  public:
    //! @name Constructor/Destructor
    //@{
    //! Default constructor.
    OrthoManager() {};

    //! Destructor.
    virtual ~OrthoManager() {};
    //@}

    //! @name Orthogonalization methods
    //@{

    /*! \brief Provides the inner product defining the orthogonality concepts.

    All concepts of orthogonality discussed in this class are with respect to this inner product.

    \note This can be different than the MvTransMv method from the multivector class. For example,
    if there is a mass matrix \c M, then this might be the \c M inner product (\f$x^HMx\f$).

     */
    virtual void innerProd( const MV &X, const MV &Y, Teuchos::SerialDenseMatrix<int,ScalarType>& Z ) const = 0;


    /// \brief Compute the norm(s) of the column(s) of X.
    ///
    /// The norm computed is the norm induced by the inner product
    /// defined by \c innerProd().
    ///
    /// \param X [in] The multivector whose columns this method will
    ///   compute norms.
    ///
    /// \param normvec [out] On output, normvec[j] is the norm of
    ///   column j of X.  This method reserves the right to resize
    ///   normvec if it does not have enough entries, but it may not
    ///   necessarily resize normvec if it has too many entries.
    virtual void norm (const MV& X, std::vector<typename Teuchos::ScalarTraits<ScalarType>::magnitudeType>& normvec) const = 0;

    /// \brief Project X against the (orthogonal) entries of Q
    ///
    /// Given a list of (mutually and internally) orthonormal bases \c
    /// Q, this method takes a multivector \c X and projects it onto
    /// the space orthogonal to the individual <tt>Q[i]</tt>,
    /// optionally returning the coefficients of \c X for the
    /// individual <tt>Q[i]</tt>.  All of this is done with respect to
    /// the inner product innerProd().
    ///
    /// After calling this routine, \c X will be orthogonal to each of
    /// the <tt>Q[i]</tt>.
    ///
    /// \param X [in/out] The multivector to be modified.  On output,
    ///   \c X will be orthogonal to <tt>Q[i]</tt> with respect to
    ///   innerProd().
    ///
    /// \param C [out] The coefficients of \c X in the \c *Q[i], with
    ///   respect to innerProd(). If <tt>C[i]</tt> is a non-null
    ///   pointer and \c *C[i] matches the dimensions of \c X and \c
    ///   *Q[i], then the coefficients computed during the
    ///   orthogonalization routine will be stored in the matrix \c
    ///   *C[i]. If <tt>C[i]</tt> is a nnon-null pointer whose size
    ///   does not match the dimensions of \c X and \c *Q[i], then a
    ///   std::invalid_argument std::exception will be
    ///   thrown. Otherwise, if <tt>C.size() < i</tt> or <tt>C[i]</tt>
    ///   is a null pointer, then the orthogonalization manager will
    ///   declare storage for the coefficients and the user will not
    ///   have access to them.
    ///
    /// \param Q [in] A list of multivector bases specifying the
    ///   subspaces to be orthogonalized against. Each <tt>Q[i]</tt>
    ///   is assumed to have orthonormal columns, and the
    ///   <tt>Q[i]</tt> are assumed to be mutually orthogonal.
    ///
    virtual void
    project (MV &X,
             Teuchos::Array<Teuchos::RCP<Teuchos::SerialDenseMatrix<int,ScalarType> > > C,
             Teuchos::ArrayView<Teuchos::RCP<const MV> > Q) const = 0;

    /*! \brief This method takes a multivector \c X and attempts to compute an orthonormal basis for \f$colspan(X)\f$, with respect to innerProd().
     *
     * This routine returns an integer \c rank stating the rank of the computed basis. If \c X does not have full rank and the normalize() routine does
     * not attempt to augment the subspace, then \c rank may be smaller than the number of columns in \c X. In this case, only the first \c rank columns of
     * output \c X and first \c rank rows of \c B will be valid.
     *
     @param X [in/out] The multivector to the modified.
       On output, \c X will have some number of orthonormal columns (with respect to innerProd()).

     @param B [out] The coefficients of the original \c X with respect to the computed basis. This matrix is not necessarily triangular; see the documentation
       for specific orthogonalization managers.

     @return Rank of the basis computed by this method.
    */
    virtual int normalize ( MV &X, Teuchos::RCP<Teuchos::SerialDenseMatrix<int,ScalarType> > B ) const = 0;

  protected:
    virtual int
    projectAndNormalizeImpl (MV &X,
                             Teuchos::Array<Teuchos::RCP<Teuchos::SerialDenseMatrix<int,ScalarType> > > C,
                             Teuchos::RCP<Teuchos::SerialDenseMatrix<int,ScalarType> > B,
                             Teuchos::ArrayView<Teuchos::RCP<const MV> > Q) const = 0;

  public:

    /// \brief Project X against the Q[i] and normalize X
    ///
    /// Given a set of bases <tt>Q[i]</tt> and a multivector \c X,
    /// this method computes an orthonormal basis for \f$colspan(X) -
    /// \sum_i colspan(Q[i])\f$.
    ///
    /// This routine returns an integer \c rank stating the rank of
    /// the computed basis. If the subspace \f$colspan(X) - \sum_i
    /// colspan(Q[i])\f$ does not have dimension as large as the
    /// number of columns of \c X and the orthogonalization manager
    /// doe not attempt to augment the subspace, then \c rank may be
    /// smaller than the number of columns of \c X. In this case, only
    /// the first \c rank columns of output \c X and first \c rank
    /// rows of \c B will be valid.
    ///
    /// \note This routine guarantees both the orthgonality
    ///   constraints against the <tt>Q[i]</tt> as well as the
    ///   orthonormality constraints. Therefore, this method is not
    ///   necessarily equivalent to calling project() followed by a
    ///   call to normalize().  See the documentation for specific
    ///   orthogonalization managers.
    ///
    /// \param X [in/out] The multivector to the modified. On output,
    ///   the relevant rows of \c X will be orthogonal to the
    ///   <tt>Q[i]</tt> and will have orthonormal columns (with
    ///   respect to innerProd()).
    ///
    /// \param C [out] The coefficients of the original \c X in the \c
    ///   *Q[i], with respect to innerProd(). If <tt>C[i]</tt> is a
    ///   non-null pointer and \c *C[i] matches the dimensions of \c X
    ///   and \c *Q[i], then the coefficients computed during the
    ///   orthogonalization routine will be stored in the matrix \c
    ///   *C[i]. If <tt>C[i]</tt> is a non-null pointer whose size does
    ///   not match the dimensions of \c X and \c *Q[i], then a
    ///   std::invalid_argument std::exception will be
    ///   thrown. Otherwise, if <tt>C.size() < i</tt> or <tt>C[i]</tt>
    ///   is a null pointer, then the orthogonalization manager will
    ///   declare storage for the coefficients and the user will not
    ///   have access to them.
    ///
    /// \param B [out] The coefficients of the original \c X with
    ///   respect to the computed basis. This matrix is not
    ///   necessarily upper triangular (as it would be for textbook
    ///   Gram-Schmidt orthogonalization of a full-rank matrix, for
    ///   example).  See the documentation for specific
    ///   orthogonalization managers.
    ///
    /// \param Q [in] A list of multivector bases specifying the
    ///   subspaces to be orthogonalized against. Each <tt>Q[i]</tt>
    ///   is assumed to have orthonormal columns, and the
    ///   <tt>Q[i]</tt> are assumed to be mutually orthogonal.
    ///
    /// \return Rank of the basis computed by this method.
    ///
    int
    projectAndNormalize (MV &X,
                         Teuchos::Array<Teuchos::RCP<Teuchos::SerialDenseMatrix<int,ScalarType> > > C,
                         Teuchos::RCP<Teuchos::SerialDenseMatrix<int,ScalarType> > B,
                         Teuchos::ArrayView<Teuchos::RCP<const MV> > Q) const
    {
      return this->projectAndNormalizeImpl (X, C, B, Q);
    }

    //@}

    //! @name Error methods
    //@{

    /*! \brief This method computes the error in orthonormality of a multivector.
     */
    virtual typename Teuchos::ScalarTraits< ScalarType >::magnitudeType
    orthonormError(const MV &X) const = 0;

    /*! \brief This method computes the error in orthogonality of two multivectors.
     */
    virtual typename Teuchos::ScalarTraits<ScalarType>::magnitudeType
    orthogError(const MV &X1, const MV &X2) const = 0;

    //@}


    //! @name Label methods
    //@{

    /*! \brief This method sets the label used by the timers in the orthogonalization manager.
     */
    virtual void setLabel(const std::string& label) = 0;

    /*! \brief This method returns the label being used by the timers in the orthogonalization manager.
     */
    virtual const std::string& getLabel() const = 0;

    //@}

  };

} // end of Belos namespace


#endif

// end of file BelosOrthoManager.hpp
