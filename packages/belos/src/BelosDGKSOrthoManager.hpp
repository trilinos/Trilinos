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


/*! \file BelosDGKSOrthoManager.hpp
  \brief Classical Gram-Schmidt (with DGKS correction) implementation of the Belos::OrthoManager class
*/

#ifndef BELOS_DGKS_ORTHOMANAGER_HPP
#define BELOS_DGKS_ORTHOMANAGER_HPP

/*!   \class Belos::DGKSOrthoManager
      \brief An implementation of the Belos::MatOrthoManager that performs orthogonalization
      using (potentially) multiple steps of classical Gram-Schmidt.

      \author Chris Baker, Ulrich Hetmaniuk, Rich Lehoucq, and Heidi Thornquist
*/

// #define ORTHO_DEBUG

#include "BelosConfigDefs.hpp"
#include "BelosMultiVecTraits.hpp"
#include "BelosOperatorTraits.hpp"
#include "BelosMatOrthoManager.hpp"

#include "Teuchos_as.hpp"
#include "Teuchos_ParameterListAcceptorDefaultBase.hpp"
#ifdef BELOS_TEUCHOS_TIME_MONITOR
#include "Teuchos_TimeMonitor.hpp"
#endif // BELOS_TEUCHOS_TIME_MONITOR

namespace Belos {

  /// \brief "Default" parameters for robustness and accuracy.
  template<class ScalarType, class MV, class OP>
  Teuchos::RCP<Teuchos::ParameterList> getDGKSDefaultParameters ();

  /// \brief "Fast" but possibly unsafe or less accurate parameters.
  template<class ScalarType, class MV, class OP>
  Teuchos::RCP<Teuchos::ParameterList> getDGKSFastParameters();

  template<class ScalarType, class MV, class OP>
  class DGKSOrthoManager :
    public MatOrthoManager<ScalarType,MV,OP>,
    public Teuchos::ParameterListAcceptorDefaultBase
  {
  private:
    typedef typename Teuchos::ScalarTraits<ScalarType>::magnitudeType MagnitudeType;
    typedef typename Teuchos::ScalarTraits<MagnitudeType> MGT;
    typedef Teuchos::ScalarTraits<ScalarType>  SCT;
    typedef MultiVecTraits<ScalarType,MV>      MVT;
    typedef OperatorTraits<ScalarType,MV,OP>   OPT;

  public:
    //! @name Constructor/Destructor
    //@{

    //! Constructor specifying re-orthogonalization tolerance.
    DGKSOrthoManager( const std::string& label = "Belos",
                      Teuchos::RCP<const OP> Op = Teuchos::null,
                      const int max_blk_ortho = max_blk_ortho_default_,
                      const MagnitudeType blk_tol = blk_tol_default_,
                      const MagnitudeType dep_tol = dep_tol_default_,
                      const MagnitudeType sing_tol = sing_tol_default_ )
      : MatOrthoManager<ScalarType,MV,OP>(Op),
        max_blk_ortho_( max_blk_ortho ),
        blk_tol_( blk_tol ),
        dep_tol_( dep_tol ),
        sing_tol_( sing_tol ),
        label_( label )
    {
#ifdef BELOS_TEUCHOS_TIME_MONITOR
      std::stringstream ss;
      ss << label_ + ": DGKS[" << max_blk_ortho_ << "]";

      std::string orthoLabel = ss.str() + ": Orthogonalization";
      timerOrtho_ = Teuchos::TimeMonitor::getNewCounter(orthoLabel);

      std::string updateLabel = ss.str() + ": Ortho (Update)";
      timerUpdate_ = Teuchos::TimeMonitor::getNewCounter(updateLabel);

      std::string normLabel = ss.str() + ": Ortho (Norm)";
      timerNorm_ = Teuchos::TimeMonitor::getNewCounter(normLabel);

      std::string ipLabel = ss.str() + ": Ortho (Inner Product)";
      timerInnerProd_ = Teuchos::TimeMonitor::getNewCounter(ipLabel);
#endif
    }

    //! Constructor that takes a list of parameters.
    DGKSOrthoManager (const Teuchos::RCP<Teuchos::ParameterList>& plist,
                      const std::string& label = "Belos",
                      Teuchos::RCP<const OP> Op = Teuchos::null)
      : MatOrthoManager<ScalarType,MV,OP>(Op),
        max_blk_ortho_ ( max_blk_ortho_default_ ),
        blk_tol_ ( blk_tol_default_ ),
        dep_tol_ ( dep_tol_default_ ),
        sing_tol_ ( sing_tol_default_ ),
        label_( label )
    {
      setParameterList (plist);

#ifdef BELOS_TEUCHOS_TIME_MONITOR
      std::stringstream ss;
      ss << label_ + ": DGKS[" << max_blk_ortho_ << "]";

      std::string orthoLabel = ss.str() + ": Orthogonalization";
      timerOrtho_ = Teuchos::TimeMonitor::getNewCounter(orthoLabel);

      std::string updateLabel = ss.str() + ": Ortho (Update)";
      timerUpdate_ = Teuchos::TimeMonitor::getNewCounter(updateLabel);

      std::string normLabel = ss.str() + ": Ortho (Norm)";
      timerNorm_ = Teuchos::TimeMonitor::getNewCounter(normLabel);

      std::string ipLabel = ss.str() + ": Ortho (Inner Product)";
      timerInnerProd_ = Teuchos::TimeMonitor::getNewCounter(ipLabel);
#endif
    }

    //! Destructor
    ~DGKSOrthoManager() {}
    //@}

    //! @name Implementation of Teuchos::ParameterListAcceptorDefaultBase interface
    //@{

    void
    setParameterList (const Teuchos::RCP<Teuchos::ParameterList>& plist)
    {
      using Teuchos::ParameterList;
      using Teuchos::parameterList;
      using Teuchos::RCP;

      RCP<const ParameterList> defaultParams = getValidParameters();
      RCP<ParameterList> params;
      if (plist.is_null()) {
        // No need to validate default parameters.
        params = parameterList (*defaultParams);
      } else {
        params = plist;
        params->validateParametersAndSetDefaults (*defaultParams);
      }

      // Using temporary variables and fetching all values before
      // setting the output arguments ensures the strong exception
      // guarantee for this function: if an exception is thrown, no
      // externally visible side effects (in this case, setting the
      // output arguments) have taken place.
      const int maxNumOrthogPasses = params->get<int> ("maxNumOrthogPasses");
      const MagnitudeType blkTol = params->get<MagnitudeType> ("blkTol");
      const MagnitudeType depTol = params->get<MagnitudeType> ("depTol");
      const MagnitudeType singTol = params->get<MagnitudeType> ("singTol");

      max_blk_ortho_ = maxNumOrthogPasses;
      blk_tol_ = blkTol;
      dep_tol_ = depTol;
      sing_tol_ = singTol;

      setMyParamList (params);
    }

    Teuchos::RCP<const Teuchos::ParameterList>
    getValidParameters () const
    {
      if (defaultParams_.is_null()) {
        defaultParams_ = Belos::getDGKSDefaultParameters<ScalarType, MV, OP>();
      }

      return defaultParams_;
    }

    //@}

    //! @name Accessor routines
    //@{

    //! Set parameter for block re-orthogonalization threshhold.
    void setBlkTol( const MagnitudeType blk_tol ) {
      // Update the parameter list as well.
      Teuchos::RCP<Teuchos::ParameterList> params = getNonconstParameterList();
      if (! params.is_null()) {
        // If it's null, then we haven't called setParameterList()
        // yet.  It's entirely possible to construct the parameter
        // list on demand, so we don't try to create the parameter
        // list here.
        params->set ("blkTol", blk_tol);
      }
      blk_tol_ = blk_tol;
    }

    //! Set parameter for re-orthogonalization threshhold.
    void setDepTol( const MagnitudeType dep_tol ) {
      // Update the parameter list as well.
      Teuchos::RCP<Teuchos::ParameterList> params = getNonconstParameterList();
      if (! params.is_null()) {
        params->set ("depTol", dep_tol);
      }
      dep_tol_ = dep_tol;
    }

    //! Set parameter for singular block detection.
    void setSingTol( const MagnitudeType sing_tol ) {
      // Update the parameter list as well.
      Teuchos::RCP<Teuchos::ParameterList> params = getNonconstParameterList();
      if (! params.is_null()) {
        params->set ("singTol", sing_tol);
      }
      sing_tol_ = sing_tol;
    }

    //! Return parameter for block re-orthogonalization threshhold.
    MagnitudeType getBlkTol() const { return blk_tol_; }

    //! Return parameter for re-orthogonalization threshhold.
    MagnitudeType getDepTol() const { return dep_tol_; }

    //! Return parameter for singular block detection.
    MagnitudeType getSingTol() const { return sing_tol_; }

    //@}


    //! @name Orthogonalization methods
    //@{

    /*! \brief Given a list of (mutually and internally) orthonormal bases \c Q, this method
     * takes a multivector \c X and projects it onto the space orthogonal to the individual <tt>Q[i]</tt>,
     * optionally returning the coefficients of \c X for the individual <tt>Q[i]</tt>. All of this is done with respect
     * to the inner product innerProd().
     *
     * After calling this routine, \c X will be orthogonal to each of the \c <tt>Q[i]</tt>.
     *
     * The method uses either one or two steps of classical Gram-Schmidt. The algebraically
     * equivalent projection matrix is \f$P_Q = I - Q Q^H Op\f$, if \c Op is the matrix specified for
     * use in the inner product. Note, this is not an orthogonal projector.
     *
     @param X [in/out] The multivector to be modified.
       On output, \c X will be orthogonal to <tt>Q[i]</tt> with respect to innerProd().

     @param MX [in/out] The image of \c X under the operator \c Op.
       If \f$ MX != 0\f$: On input, this is expected to be consistent with \c X. On output, this is updated consistent with updates to \c X.
       If \f$ MX == 0\f$ or \f$ Op == 0\f$: \c MX is not referenced.

     @param C [out] The coefficients of \c X in the \c *Q[i], with respect to innerProd(). If <tt>C[i]</tt> is a non-null pointer
       and \c *C[i] matches the dimensions of \c X and \c *Q[i], then the coefficients computed during the orthogonalization
       routine will be stored in the matrix \c *C[i]. If <tt>C[i]</tt> is a non-null pointer whose size does not match the dimensions of
       \c X and \c *Q[i], then a std::invalid_argument std::exception will be thrown. Otherwise, if <tt>C.size() < i</tt> or <tt>C[i]</tt> is a null
       pointer, then the orthogonalization manager will declare storage for the coefficients and the user will not have access to them.

     @param Q [in] A list of multivector bases specifying the subspaces to be orthogonalized against. Each <tt>Q[i]</tt> is assumed to have
     orthonormal columns, and the <tt>Q[i]</tt> are assumed to be mutually orthogonal.
    */
    void project ( MV &X, Teuchos::RCP<MV> MX,
                   Teuchos::Array<Teuchos::RCP<Teuchos::SerialDenseMatrix<int,ScalarType> > > C,
                   Teuchos::ArrayView<Teuchos::RCP<const MV> > Q) const;


    /*! \brief This method calls project(X,Teuchos::null,C,Q); see documentation for that function.
    */
    void project ( MV &X,
                   Teuchos::Array<Teuchos::RCP<Teuchos::SerialDenseMatrix<int,ScalarType> > > C,
                   Teuchos::ArrayView<Teuchos::RCP<const MV> > Q) const {
      project(X,Teuchos::null,C,Q);
    }



    /*! \brief This method takes a multivector \c X and attempts to compute an orthonormal basis for \f$colspan(X)\f$, with respect to innerProd().
     *
     * The method uses classical Gram-Schmidt, so that the coefficient matrix \c B is upper triangular.
     *
     * This routine returns an integer \c rank stating the rank of the computed basis. If \c X does not have full rank and the normalize() routine does
     * not attempt to augment the subspace, then \c rank may be smaller than the number of columns in \c X. In this case, only the first \c rank columns of
     * output \c X and first \c rank rows of \c B will be valid.
     *
     * The method attempts to find a basis with dimension the same as the number of columns in \c X. It does this by augmenting linearly dependant
     * vectors in \c X with random directions. A finite number of these attempts will be made; therefore, it is possible that the dimension of the
     * computed basis is less than the number of vectors in \c X.
     *
     @param X [in/out] The multivector to the modified.
       On output, \c X will have some number of orthonormal columns (with respect to innerProd()).

     @param MX [in/out] The image of \c X under the operator \c Op.
       If \f$ MX != 0\f$: On input, this is expected to be consistent with \c X. On output, this is updated consistent with updates to \c X.
       If \f$ MX == 0\f$ or \f$ Op == 0\f$: \c MX is not referenced.

     @param B [out] The coefficients of the original \c X with respect to the computed basis. The first rows in \c B
            corresponding to the valid columns in \c X will be upper triangular.

     @return Rank of the basis computed by this method.
    */
    int normalize ( MV &X, Teuchos::RCP<MV> MX,
                    Teuchos::RCP<Teuchos::SerialDenseMatrix<int,ScalarType> > B) const;


    /*! \brief This method calls normalize(X,Teuchos::null,B); see documentation for that function.
    */
    int normalize ( MV &X, Teuchos::RCP<Teuchos::SerialDenseMatrix<int,ScalarType> > B ) const {
      return normalize(X,Teuchos::null,B);
    }

  protected:

    /*! \brief Given a set of bases <tt>Q[i]</tt> and a multivector \c X, this method computes an orthonormal basis for \f$colspan(X) - \sum_i colspan(Q[i])\f$.
     *
     *  This routine returns an integer \c rank stating the rank of
     *  the computed basis. If the subspace \f$colspan(X) - \sum_i
     *  colspan(Q[i])\f$ does not have dimension as large as the
     *  number of columns of \c X and the orthogonalization manager
     *  doe not attempt to augment the subspace, then \c rank may be
     *  smaller than the number of columns of \c X. In this case, only
     *  the first \c rank columns of output \c X and first \c rank
     *  rows of \c B will be valid.
     *
     * The method attempts to find a basis with dimension the same as
     * the number of columns in \c X. It does this by augmenting
     * linearly dependant vectors with random directions. A finite
     * number of these attempts will be made; therefore, it is
     * possible that the dimension of the computed basis is less than
     * the number of vectors in \c X.
     *
     @param X [in/out] The multivector to the modified.  On output,
       the relevant rows of \c X will be orthogonal to the
       <tt>Q[i]</tt> and will have orthonormal columns (with respect
       to innerProd()).

     @param MX [in/out] The image of \c X under the operator \c Op.
       If \f$ MX != 0\f$: On input, this is expected to be consistent
       with \c X. On output, this is updated consistent with updates
       to \c X.  If \f$ MX == 0\f$ or \f$ Op == 0\f$: \c MX is not
       referenced.

     @param C [out] The coefficients of the original \c X in the \c
     *Q[i], with respect to innerProd(). If <tt>C[i]</tt> is a
     non-null pointer and \c *C[i] matches the dimensions of \c X and
     \c *Q[i], then the coefficients computed during the
     orthogonalization routine will be stored in the matrix \c
     *C[i]. If <tt>C[i]</tt> is a non-null pointer whose size does not
     match the dimensions of \c X and \c *Q[i], then *C[i] will first
     be resized to the correct size.  This will destroy the original
     contents of the matrix.  (This is a change from previous
     behavior, in which a std::invalid_argument exception was thrown
     if *C[i] was of the wrong size.)  Otherwise, if <tt>C.size() <
     i<\tt> or <tt>C[i]</tt> is a null pointer, then the
     orthogonalization manager will declare storage for the
     coefficients and the user will not have access to them.

     @param B [out] The coefficients of the original \c X with respect
       to the computed basis. The first rows in \c B corresponding to
       the valid columns in \c X will be upper triangular.

     @param Q [in] A list of multivector bases specifying the
       subspaces to be orthogonalized against. Each <tt>Q[i]</tt> is
       assumed to have orthonormal columns, and the <tt>Q[i]</tt> are
       assumed to be mutually orthogonal.

     @return Rank of the basis computed by this method.
    */
    virtual int
    projectAndNormalizeWithMxImpl (MV &X,
                                   Teuchos::RCP<MV> MX,
                                   Teuchos::Array<Teuchos::RCP<Teuchos::SerialDenseMatrix<int,ScalarType> > > C,
                                   Teuchos::RCP<Teuchos::SerialDenseMatrix<int,ScalarType> > B,
                                   Teuchos::ArrayView<Teuchos::RCP<const MV> > Q) const;

  public:
    //@}
    //! @name Error methods
    //@{

    /// \brief Compute \fn$\| X^* M X - I \|_F\fn$
    ///
    /// This method computes the error in orthonormality of a
    /// multivector, measured as the Frobenius norm of the difference
    /// <tt>innerProd(X,X) - I</tt>.
    typename Teuchos::ScalarTraits<ScalarType>::magnitudeType
    orthonormError(const MV &X) const {
      return orthonormError(X,Teuchos::null);
    }

    /// \brief Compute \fn$\| X^* M X - I \|_F\fn$
    ///
    /// This method computes the error in orthonormality of a
    /// multivector, measured as the Frobenius norm of the difference
    /// <tt>innerProd(X,X) - I</tt>.  The method has the option of
    /// exploiting a caller-provided \c MX, which is used if not null.
    typename Teuchos::ScalarTraits<ScalarType>::magnitudeType
    orthonormError(const MV &X, Teuchos::RCP<const MV> MX) const;

    /*! \brief This method computes the error in orthogonality of two multivectors, measured
     * as the Frobenius norm of <tt>innerProd(X,Y)</tt>.
     */
    typename Teuchos::ScalarTraits<ScalarType>::magnitudeType
    orthogError(const MV &X1, const MV &X2) const {
      return orthogError(X1,Teuchos::null,X2);
    }

    /*! \brief This method computes the error in orthogonality of two multivectors, measured
     * as the Frobenius norm of <tt>innerProd(X,Y)</tt>.
     *  The method has the option of exploiting a caller-provided \c MX.
     */
    typename Teuchos::ScalarTraits<ScalarType>::magnitudeType
    orthogError(const MV &X1, Teuchos::RCP<const MV> MX1, const MV &X2) const;

    //@}

    //! @name Label methods
    //@{

    /*! \brief This method sets the label used by the timers in the orthogonalization manager.
     */
    void setLabel(const std::string& label);

    /*! \brief This method returns the label being used by the timers in the orthogonalization manager.
     */
    const std::string& getLabel() const { return label_; }

    //@}

    //! @name Default orthogonalization constants
    //@{

    //! Max number of (re)orthogonalization steps, including the first (default).
    static const int max_blk_ortho_default_;
    //! Block reorthogonalization threshold (default).
    static const MagnitudeType blk_tol_default_;
    //! (Non-block) reorthogonalization threshold (default).
    static const MagnitudeType dep_tol_default_;
    //! Singular block detection threshold (default).
    static const MagnitudeType sing_tol_default_;

    //! Max number of (re)orthogonalization steps, including the first (fast).
    static const int max_blk_ortho_fast_;
    //! Block reorthogonalization threshold (fast).
    static const MagnitudeType blk_tol_fast_;
    //! (Non-block) reorthogonalization threshold (fast).
    static const MagnitudeType dep_tol_fast_;
    //! Singular block detection threshold (fast).
    static const MagnitudeType sing_tol_fast_;

    //@}

  private:

    //! Max number of (re)orthogonalization steps, including the first.
    int max_blk_ortho_;
    //! Block reorthogonalization threshold.
    MagnitudeType blk_tol_;
    //! (Non-block) reorthogonalization threshold.
    MagnitudeType dep_tol_;
    //! Singular block detection threshold.
    MagnitudeType sing_tol_;

    //! Label for timer(s).
    std::string label_;
#ifdef BELOS_TEUCHOS_TIME_MONITOR
    Teuchos::RCP<Teuchos::Time> timerOrtho_, timerUpdate_, timerNorm_, timerInnerProd_;
#endif // BELOS_TEUCHOS_TIME_MONITOR

    //! Default parameter list.
    mutable Teuchos::RCP<Teuchos::ParameterList> defaultParams_;

    //! Routine to find an orthonormal basis for X
    int findBasis(MV &X, Teuchos::RCP<MV> MX,
                  Teuchos::RCP<Teuchos::SerialDenseMatrix<int,ScalarType> > C,
                  bool completeBasis, int howMany = -1 ) const;

    //! Routine to compute the block orthogonalization
    bool blkOrtho1 ( MV &X, Teuchos::RCP<MV> MX,
                    Teuchos::Array<Teuchos::RCP<Teuchos::SerialDenseMatrix<int,ScalarType> > > C,
                    Teuchos::ArrayView<Teuchos::RCP<const MV> > Q) const;

    //! Routine to compute the block orthogonalization
    bool blkOrtho ( MV &X, Teuchos::RCP<MV> MX,
                    Teuchos::Array<Teuchos::RCP<Teuchos::SerialDenseMatrix<int,ScalarType> > > C,
                    Teuchos::ArrayView<Teuchos::RCP<const MV> > Q) const;

    /// Project X against QQ and normalize X, one vector at a time
    ///
    /// \note QQ is called QQ, rather than Q, because we convert it
    ///   internally from an ArrayView to an Array (named Q inside).
    ///   This is because the C++ compiler doesn't know how to do type
    ///   inference (Array has a constructor that takes an ArrayView
    ///   input).  This routine wants an Array rather than an
    ///   ArrayView internally, because it likes to add (via
    ///   push_back()) and remove (via resize()) elements to the Q
    ///   array.  Remember that Arrays can be passed by value, just
    ///   like std::vector objects, so this routine can add whatever
    ///   it likes to the Q array without changing it from the
    ///   caller's perspective.
    int blkOrthoSing ( MV &X, Teuchos::RCP<MV> MX,
                       Teuchos::Array<Teuchos::RCP<Teuchos::SerialDenseMatrix<int,ScalarType> > > C,
                       Teuchos::RCP<Teuchos::SerialDenseMatrix<int,ScalarType> > B,
                       Teuchos::ArrayView<Teuchos::RCP<const MV> > QQ) const;
  };

  // Set static variables.
  template<class ScalarType, class MV, class OP>
  const int DGKSOrthoManager<ScalarType,MV,OP>::max_blk_ortho_default_ = 2;

  template<class ScalarType, class MV, class OP>
  const typename DGKSOrthoManager<ScalarType,MV,OP>::MagnitudeType
  DGKSOrthoManager<ScalarType,MV,OP>::blk_tol_default_
    = 10*Teuchos::ScalarTraits<typename DGKSOrthoManager<ScalarType,MV,OP>::MagnitudeType>::squareroot(
      Teuchos::ScalarTraits<typename DGKSOrthoManager<ScalarType,MV,OP>::MagnitudeType>::eps() );

  template<class ScalarType, class MV, class OP>
  const typename DGKSOrthoManager<ScalarType,MV,OP>::MagnitudeType
  DGKSOrthoManager<ScalarType,MV,OP>::dep_tol_default_ 
    = Teuchos::ScalarTraits<typename DGKSOrthoManager<ScalarType,MV,OP>::MagnitudeType>::one()
      / Teuchos::ScalarTraits<typename DGKSOrthoManager<ScalarType,MV,OP>::MagnitudeType>::squareroot( 
      2*Teuchos::ScalarTraits<typename DGKSOrthoManager<ScalarType,MV,OP>::MagnitudeType>::one() );

  template<class ScalarType, class MV, class OP>
  const typename DGKSOrthoManager<ScalarType,MV,OP>::MagnitudeType
  DGKSOrthoManager<ScalarType,MV,OP>::sing_tol_default_ 
    = 10*Teuchos::ScalarTraits<typename DGKSOrthoManager<ScalarType,MV,OP>::MagnitudeType>::eps();

  template<class ScalarType, class MV, class OP>
  const int DGKSOrthoManager<ScalarType,MV,OP>::max_blk_ortho_fast_ = 1;

  template<class ScalarType, class MV, class OP>
  const typename DGKSOrthoManager<ScalarType,MV,OP>::MagnitudeType
  DGKSOrthoManager<ScalarType,MV,OP>::blk_tol_fast_
    = Teuchos::ScalarTraits<typename DGKSOrthoManager<ScalarType,MV,OP>::MagnitudeType>::zero();

  template<class ScalarType, class MV, class OP>
  const typename DGKSOrthoManager<ScalarType,MV,OP>::MagnitudeType
  DGKSOrthoManager<ScalarType,MV,OP>::dep_tol_fast_
    = Teuchos::ScalarTraits<typename DGKSOrthoManager<ScalarType,MV,OP>::MagnitudeType>::zero();

  template<class ScalarType, class MV, class OP>
  const typename DGKSOrthoManager<ScalarType,MV,OP>::MagnitudeType
  DGKSOrthoManager<ScalarType,MV,OP>::sing_tol_fast_
    = Teuchos::ScalarTraits<typename DGKSOrthoManager<ScalarType,MV,OP>::MagnitudeType>::zero();

  //////////////////////////////////////////////////////////////////////////////////////////////////
  // Set the label for this orthogonalization manager and create new timers if it's changed
  template<class ScalarType, class MV, class OP>
  void DGKSOrthoManager<ScalarType,MV,OP>::setLabel(const std::string& label)
  {
    if (label != label_) {
      label_ = label;
#ifdef BELOS_TEUCHOS_TIME_MONITOR
      std::stringstream ss;
      ss << label_ + ": DGKS[" << max_blk_ortho_ << "]";

      std::string orthoLabel = ss.str() + ": Orthogonalization";
      timerOrtho_ = Teuchos::TimeMonitor::getNewCounter(orthoLabel);

      std::string updateLabel = ss.str() + ": Ortho (Update)";
      timerUpdate_ = Teuchos::TimeMonitor::getNewCounter(updateLabel);

      std::string normLabel = ss.str() + ": Ortho (Norm)";
      timerNorm_ = Teuchos::TimeMonitor::getNewCounter(normLabel);

      std::string ipLabel = ss.str() + ": Ortho (Inner Product)";
      timerInnerProd_ = Teuchos::TimeMonitor::getNewCounter(ipLabel);
#endif
    }
  }

  //////////////////////////////////////////////////////////////////////////////////////////////////
  // Compute the distance from orthonormality
  template<class ScalarType, class MV, class OP>
  typename Teuchos::ScalarTraits<ScalarType>::magnitudeType
  DGKSOrthoManager<ScalarType,MV,OP>::orthonormError(const MV &X, Teuchos::RCP<const MV> MX) const {
    const ScalarType ONE = SCT::one();
    int rank = MVT::GetNumberVecs(X);
    Teuchos::SerialDenseMatrix<int,ScalarType> xTx(rank,rank);
    {
#ifdef BELOS_TEUCHOS_TIME_MONITOR
    Teuchos::TimeMonitor innerProdTimer( *timerInnerProd_ );
#endif
    MatOrthoManager<ScalarType,MV,OP>::innerProd(X,X,MX,xTx);
    }
    for (int i=0; i<rank; i++) {
      xTx(i,i) -= ONE;
    }
    return xTx.normFrobenius();
  }

  //////////////////////////////////////////////////////////////////////////////////////////////////
  // Compute the distance from orthogonality
  template<class ScalarType, class MV, class OP>
  typename Teuchos::ScalarTraits<ScalarType>::magnitudeType
  DGKSOrthoManager<ScalarType,MV,OP>::orthogError(const MV &X1, Teuchos::RCP<const MV> MX1, const MV &X2) const {
    int r1 = MVT::GetNumberVecs(X1);
    int r2  = MVT::GetNumberVecs(X2);
    Teuchos::SerialDenseMatrix<int,ScalarType> xTx(r2,r1);
    {
#ifdef BELOS_TEUCHOS_TIME_MONITOR
    Teuchos::TimeMonitor innerProdTimer( *timerInnerProd_ );
#endif
    MatOrthoManager<ScalarType,MV,OP>::innerProd(X2,X1,MX1,xTx);
    }
    return xTx.normFrobenius();
  }

  //////////////////////////////////////////////////////////////////////////////////////////////////
  // Find an Op-orthonormal basis for span(X) - span(W)
  template<class ScalarType, class MV, class OP>
  int
  DGKSOrthoManager<ScalarType, MV, OP>::
  projectAndNormalizeWithMxImpl (MV &X,
                                 Teuchos::RCP<MV> MX,
                                 Teuchos::Array<Teuchos::RCP<Teuchos::SerialDenseMatrix<int,ScalarType> > > C,
                                 Teuchos::RCP<Teuchos::SerialDenseMatrix<int,ScalarType> > B,
                                 Teuchos::ArrayView<Teuchos::RCP<const MV> > Q) const
  {
    using Teuchos::Array;
    using Teuchos::null;
    using Teuchos::is_null;
    using Teuchos::RCP;
    using Teuchos::rcp;
    using Teuchos::SerialDenseMatrix;
    typedef SerialDenseMatrix< int, ScalarType > serial_dense_matrix_type;
    typedef typename Array< RCP< const MV > >::size_type size_type;

#ifdef BELOS_TEUCHOS_TIME_MONITOR
    Teuchos::TimeMonitor orthotimer(*timerOrtho_);
#endif

    ScalarType    ONE  = SCT::one();
    const MagnitudeType ZERO = SCT::magnitude(SCT::zero());

    int nq = Q.size();
    int xc = MVT::GetNumberVecs( X );
    ptrdiff_t xr = MVT::GetGlobalLength( X );
    int rank = xc;

    // If the user doesn't want to store the normalization
    // coefficients, allocate some local memory for them.  This will
    // go away at the end of this method.
    if (is_null (B)) {
      B = rcp (new serial_dense_matrix_type (xc, xc));
    }
    // Likewise, if the user doesn't want to store the projection
    // coefficients, allocate some local memory for them.  Also make
    // sure that all the entries of C are the right size.  We're going
    // to overwrite them anyway, so we don't have to worry about the
    // contents (other than to resize them if they are the wrong
    // size).
    if (C.size() < nq)
      C.resize (nq);
    for (size_type k = 0; k < nq; ++k)
      {
        const int numRows = MVT::GetNumberVecs (*Q[k]);
        const int numCols = xc; // Number of vectors in X

        if (is_null (C[k]))
          C[k] = rcp (new serial_dense_matrix_type (numRows, numCols));
        else if (C[k]->numRows() != numRows || C[k]->numCols() != numCols)
          {
            int err = C[k]->reshape (numRows, numCols);
            TEUCHOS_TEST_FOR_EXCEPTION(err != 0, std::runtime_error,
                               "DGKS orthogonalization: failed to reshape "
                               "C[" << k << "] (the array of block "
                               "coefficients resulting from projecting X "
                               "against Q[1:" << nq << "]).");
          }
      }

    /******   DO NO MODIFY *MX IF _hasOp == false   ******/
    if (this->_hasOp) {
      if (MX == Teuchos::null) {
        // we need to allocate space for MX
        MX = MVT::Clone(X,MVT::GetNumberVecs(X));
        OPT::Apply(*(this->_Op),X,*MX);
      }
    }
    else {
      // Op == I  -->  MX = X (ignore it if the user passed it in)
      MX = Teuchos::rcp( &X, false );
    }

    int mxc = MVT::GetNumberVecs( *MX );
    ptrdiff_t mxr = MVT::GetGlobalLength( *MX );

    // short-circuit
    TEUCHOS_TEST_FOR_EXCEPTION( xc == 0 || xr == 0, std::invalid_argument, "Belos::DGKSOrthoManager::projectAndNormalize(): X must be non-empty" );

    int numbas = 0;
    for (int i=0; i<nq; i++) {
      numbas += MVT::GetNumberVecs( *Q[i] );
    }

    // check size of B
    TEUCHOS_TEST_FOR_EXCEPTION( B->numRows() != xc || B->numCols() != xc, std::invalid_argument,
                        "Belos::DGKSOrthoManager::projectAndNormalize(): Size of X must be consistant with size of B" );
    // check size of X and MX
    TEUCHOS_TEST_FOR_EXCEPTION( xc<0 || xr<0 || mxc<0 || mxr<0, std::invalid_argument,
                        "Belos::DGKSOrthoManager::projectAndNormalize(): MVT returned negative dimensions for X,MX" );
    // check size of X w.r.t. MX
    TEUCHOS_TEST_FOR_EXCEPTION( xc!=mxc || xr!=mxr, std::invalid_argument,
                        "Belos::DGKSOrthoManager::projectAndNormalize(): Size of X must be consistant with size of MX" );
    // check feasibility
    //TEUCHOS_TEST_FOR_EXCEPTION( numbas+xc > xr, std::invalid_argument,
    //                    "Belos::DGKSOrthoManager::projectAndNormalize(): Orthogonality constraints not feasible" );

    // Some flags for checking dependency returns from the internal orthogonalization methods
    bool dep_flg = false;

    if (xc == 1) {

      // Use the cheaper block orthogonalization.
      // NOTE: Don't check for dependencies because the update has one vector.
      dep_flg = blkOrtho1( X, MX, C, Q );

      // Normalize the new block X
      if ( B == Teuchos::null ) {
        B = Teuchos::rcp( new Teuchos::SerialDenseMatrix<int,ScalarType>(xc,xc) );
      }
      std::vector<ScalarType> diag(1);
      {
#ifdef BELOS_TEUCHOS_TIME_MONITOR
        Teuchos::TimeMonitor normTimer( *timerNorm_ );
#endif
        MVT::MvDot( X, *MX, diag );
      }
      (*B)(0,0) = SCT::squareroot(SCT::magnitude(diag[0]));

      if (SCT::magnitude((*B)(0,0)) > ZERO) {
        rank = 1;
        MVT::MvScale( X, ONE/(*B)(0,0) );
        if (this->_hasOp) {
          // Update MXj.
          MVT::MvScale( *MX, ONE/(*B)(0,0) );
        }
      }
    }
    else {

      // Make a temporary copy of X and MX, just in case a block dependency is detected.
      Teuchos::RCP<MV> tmpX, tmpMX;
      tmpX = MVT::CloneCopy(X);
      if (this->_hasOp) {
        tmpMX = MVT::CloneCopy(*MX);
      }

      // Use the cheaper block orthogonalization.
      dep_flg = blkOrtho( X, MX, C, Q );

      // If a dependency has been detected in this block, then perform
      // the more expensive single-vector orthogonalization.
      if (dep_flg) {
        rank = blkOrthoSing( *tmpX, tmpMX, C, B, Q );

        // Copy tmpX back into X.
        MVT::Assign( *tmpX, X );
        if (this->_hasOp) {
          MVT::Assign( *tmpMX, *MX );
        }
      }
      else {
        // There is no dependency, so orthonormalize new block X
        rank = findBasis( X, MX, B, false );
        if (rank < xc) {
          // A dependency was found during orthonormalization of X,
          // rerun orthogonalization using more expensive single-vector orthogonalization.
          rank = blkOrthoSing( *tmpX, tmpMX, C, B, Q );

          // Copy tmpX back into X.
          MVT::Assign( *tmpX, X );
          if (this->_hasOp) {
            MVT::Assign( *tmpMX, *MX );
          }
        }
      }
    } // if (xc == 1) 

    // this should not raise an std::exception; but our post-conditions oblige us to check
    TEUCHOS_TEST_FOR_EXCEPTION( rank > xc || rank < 0, std::logic_error,
                        "Belos::DGKSOrthoManager::projectAndNormalize(): Debug error in rank variable." );

    // Return the rank of X.
    return rank;
  }



  //////////////////////////////////////////////////////////////////////////////////////////////////
  // Find an Op-orthonormal basis for span(X), with rank numvectors(X)
  template<class ScalarType, class MV, class OP>
  int DGKSOrthoManager<ScalarType, MV, OP>::normalize(
                                MV &X, Teuchos::RCP<MV> MX,
                                Teuchos::RCP<Teuchos::SerialDenseMatrix<int,ScalarType> > B ) const {

#ifdef BELOS_TEUCHOS_TIME_MONITOR
    Teuchos::TimeMonitor orthotimer(*timerOrtho_);
#endif

    // call findBasis, with the instruction to try to generate a basis of rank numvecs(X)
    return findBasis(X, MX, B, true);
    
  }



  //////////////////////////////////////////////////////////////////////////////////////////////////
  template<class ScalarType, class MV, class OP>
  void DGKSOrthoManager<ScalarType, MV, OP>::project(
                          MV &X, Teuchos::RCP<MV> MX,
                          Teuchos::Array<Teuchos::RCP<Teuchos::SerialDenseMatrix<int,ScalarType> > > C,
                          Teuchos::ArrayView<Teuchos::RCP<const MV> > Q) const {
    // For the inner product defined by the operator Op or the identity (Op == 0)
    //   -> Orthogonalize X against each Q[i]
    // Modify MX accordingly
    //
    // Note that when Op is 0, MX is not referenced
    //
    // Parameter variables
    //
    // X  : Vectors to be transformed
    //
    // MX : Image of the block vector X by the mass matrix
    //
    // Q  : Bases to orthogonalize against. These are assumed orthonormal, mutually and independently.
    //

#ifdef BELOS_TEUCHOS_TIME_MONITOR
    Teuchos::TimeMonitor orthotimer(*timerOrtho_);
#endif

    int xc = MVT::GetNumberVecs( X );
    ptrdiff_t xr = MVT::GetGlobalLength( X );
    int nq = Q.size();
    std::vector<int> qcs(nq);
    // short-circuit
    if (nq == 0 || xc == 0 || xr == 0) {
      return;
    }
    ptrdiff_t qr = MVT::GetGlobalLength ( *Q[0] );
    // if we don't have enough C, expand it with null references
    // if we have too many, resize to throw away the latter ones
    // if we have exactly as many as we have Q, this call has no effect
    C.resize(nq);


    /******   DO NO MODIFY *MX IF _hasOp == false   ******/
    if (this->_hasOp) {
      if (MX == Teuchos::null) {
        // we need to allocate space for MX
        MX = MVT::Clone(X,MVT::GetNumberVecs(X));
        OPT::Apply(*(this->_Op),X,*MX);
      }
    }
    else {
      // Op == I  -->  MX = X (ignore it if the user passed it in)
      MX = Teuchos::rcp( &X, false );
    }
    int mxc = MVT::GetNumberVecs( *MX );
    ptrdiff_t mxr = MVT::GetGlobalLength( *MX );

    // check size of X and Q w.r.t. common sense
    TEUCHOS_TEST_FOR_EXCEPTION( xc<0 || xr<0 || mxc<0 || mxr<0, std::invalid_argument,
                        "Belos::DGKSOrthoManager::project(): MVT returned negative dimensions for X,MX" );
    // check size of X w.r.t. MX and Q
    TEUCHOS_TEST_FOR_EXCEPTION( xc!=mxc || xr!=mxr || xr!=qr, std::invalid_argument,
                        "Belos::DGKSOrthoManager::project(): Size of X not consistant with MX,Q" );

    // tally up size of all Q and check/allocate C
    int baslen = 0;
    for (int i=0; i<nq; i++) {
      TEUCHOS_TEST_FOR_EXCEPTION( MVT::GetGlobalLength( *Q[i] ) != qr, std::invalid_argument,
                          "Belos::DGKSOrthoManager::project(): Q lengths not mutually consistant" );
      qcs[i] = MVT::GetNumberVecs( *Q[i] );
      TEUCHOS_TEST_FOR_EXCEPTION( qr < qcs[i], std::invalid_argument,
                          "Belos::DGKSOrthoManager::project(): Q has less rows than columns" );
      baslen += qcs[i];

      // check size of C[i]
      if ( C[i] == Teuchos::null ) {
        C[i] = Teuchos::rcp( new Teuchos::SerialDenseMatrix<int,ScalarType>(qcs[i],xc) );
      }
      else {
        TEUCHOS_TEST_FOR_EXCEPTION( C[i]->numRows() != qcs[i] || C[i]->numCols() != xc , std::invalid_argument,
                           "Belos::DGKSOrthoManager::project(): Size of Q not consistant with size of C" );
      }
    }

    // Use the cheaper block orthogonalization, don't check for rank deficiency.
    blkOrtho( X, MX, C, Q );

  }

  //////////////////////////////////////////////////////////////////////////////////////////////////
  // Find an Op-orthonormal basis for span(X), with the option of extending the subspace so that
  // the rank is numvectors(X)
  template<class ScalarType, class MV, class OP>
  int DGKSOrthoManager<ScalarType, MV, OP>::findBasis(
                                                      MV &X, Teuchos::RCP<MV> MX,
                                                      Teuchos::RCP<Teuchos::SerialDenseMatrix<int,ScalarType> > B,
                                                      bool completeBasis, int howMany ) const {
    // For the inner product defined by the operator Op or the identity (Op == 0)
    //   -> Orthonormalize X
    // Modify MX accordingly
    //
    // Note that when Op is 0, MX is not referenced
    //
    // Parameter variables
    //
    // X  : Vectors to be orthonormalized
    //
    // MX : Image of the multivector X under the operator Op
    //
    // Op  : Pointer to the operator for the inner product
    //
    //

    const ScalarType ONE  = SCT::one();
    const MagnitudeType ZERO = SCT::magnitude(SCT::zero());

    int xc = MVT::GetNumberVecs( X );
    ptrdiff_t xr = MVT::GetGlobalLength( X );

    if (howMany == -1) {
      howMany = xc;
    }

    /*******************************************************
     *  If _hasOp == false, we will not reference MX below *
     *******************************************************/

    // if Op==null, MX == X (via pointer)
    // Otherwise, either the user passed in MX or we will allocated and compute it
    if (this->_hasOp) {
      if (MX == Teuchos::null) {
        // we need to allocate space for MX
        MX = MVT::Clone(X,xc);
        OPT::Apply(*(this->_Op),X,*MX);
      }
    }

    /* if the user doesn't want to store the coefficienets,
     * allocate some local memory for them
     */
    if ( B == Teuchos::null ) {
      B = Teuchos::rcp( new Teuchos::SerialDenseMatrix<int,ScalarType>(xc,xc) );
    }

    int mxc = (this->_hasOp) ? MVT::GetNumberVecs( *MX ) : xc;
    ptrdiff_t mxr = (this->_hasOp) ? MVT::GetGlobalLength( *MX ) : xr;

    // check size of C, B
    TEUCHOS_TEST_FOR_EXCEPTION( xc == 0 || xr == 0, std::invalid_argument,
                        "Belos::DGKSOrthoManager::findBasis(): X must be non-empty" );
    TEUCHOS_TEST_FOR_EXCEPTION( B->numRows() != xc || B->numCols() != xc, std::invalid_argument,
                        "Belos::DGKSOrthoManager::findBasis(): Size of X not consistant with size of B" );
    TEUCHOS_TEST_FOR_EXCEPTION( xc != mxc || xr != mxr, std::invalid_argument,
                        "Belos::DGKSOrthoManager::findBasis(): Size of X not consistant with size of MX" );
    TEUCHOS_TEST_FOR_EXCEPTION( static_cast<ptrdiff_t>(xc) > xr, std::invalid_argument,
                        "Belos::DGKSOrthoManager::findBasis(): Size of X not feasible for normalization" );
    TEUCHOS_TEST_FOR_EXCEPTION( howMany < 0 || howMany > xc, std::invalid_argument,
                        "Belos::DGKSOrthoManager::findBasis(): Invalid howMany parameter" );

    /* xstart is which column we are starting the process with, based on howMany
     * columns before xstart are assumed to be Op-orthonormal already
     */
    int xstart = xc - howMany;

    for (int j = xstart; j < xc; j++) {

      // numX is
      // * number of currently orthonormal columns of X
      // * the index of the current column of X
      int numX = j;
      bool addVec = false;

      // Get a view of the vector currently being worked on.
      std::vector<int> index(1);
      index[0] = numX;
      Teuchos::RCP<MV> Xj = MVT::CloneViewNonConst( X, index );
      Teuchos::RCP<MV> MXj;
      if ((this->_hasOp)) {
        // MXj is a view of the current vector in MX
        MXj = MVT::CloneViewNonConst( *MX, index );
      }
      else {
        // MXj is a pointer to Xj, and MUST NOT be modified
        MXj = Xj;
      }

      // Get a view of the previous vectors.
      std::vector<int> prev_idx( numX );
      Teuchos::RCP<const MV> prevX, prevMX;
      Teuchos::RCP<MV> oldMXj;

      if (numX > 0) {
        for (int i=0; i<numX; i++) {
          prev_idx[i] = i;
        }
        prevX = MVT::CloneView( X, prev_idx );
        if (this->_hasOp) {
          prevMX = MVT::CloneView( *MX, prev_idx );
        }

        oldMXj = MVT::CloneCopy( *MXj );
      }

      // Make storage for these Gram-Schmidt iterations.
      Teuchos::SerialDenseMatrix<int,ScalarType> product(numX, 1);
      std::vector<ScalarType> oldDot( 1 ), newDot( 1 );
      //
      // Save old MXj vector and compute Op-norm
      //
      {
#ifdef BELOS_TEUCHOS_TIME_MONITOR
      Teuchos::TimeMonitor normTimer( *timerNorm_ );
#endif
      MVT::MvDot( *Xj, *MXj, oldDot );
      }
      // Xj^H Op Xj should be real and positive, by the hermitian positive definiteness of Op
      TEUCHOS_TEST_FOR_EXCEPTION( SCT::real(oldDot[0]) < ZERO, OrthoError,
                          "Belos::DGKSOrthoManager::findBasis(): Negative definiteness discovered in inner product" );

      if (numX > 0) {
        // Apply the first step of Gram-Schmidt

        // product <- prevX^T MXj
        {
#ifdef BELOS_TEUCHOS_TIME_MONITOR
        Teuchos::TimeMonitor innerProdTimer( *timerInnerProd_ );
#endif
        MatOrthoManager<ScalarType,MV,OP>::innerProd(*prevX,*Xj,MXj,product);
        }
        // Xj <- Xj - prevX prevX^T MXj
        //     = Xj - prevX product
        {
#ifdef BELOS_TEUCHOS_TIME_MONITOR
        Teuchos::TimeMonitor updateTimer( *timerUpdate_ );
#endif
        MVT::MvTimesMatAddMv( -ONE, *prevX, product, ONE, *Xj );
        }

        // Update MXj
        if (this->_hasOp) {
          // MXj <- Op*Xj_new
          //      = Op*(Xj_old - prevX prevX^T MXj)
          //      = MXj - prevMX product
#ifdef BELOS_TEUCHOS_TIME_MONITOR
          Teuchos::TimeMonitor updateTimer( *timerUpdate_ );
#endif
          MVT::MvTimesMatAddMv( -ONE, *prevMX, product, ONE, *MXj );
        }

        // Compute new Op-norm
        {
#ifdef BELOS_TEUCHOS_TIME_MONITOR
        Teuchos::TimeMonitor normTimer( *timerNorm_ );
#endif
        MVT::MvDot( *Xj, *MXj, newDot );
        }

        // Check if a correction is needed.
        if ( MGT::squareroot(SCT::magnitude(newDot[0])) < dep_tol_*MGT::squareroot(SCT::magnitude(oldDot[0])) ) {
          // Apply the second step of Gram-Schmidt
          // This is the same as above
          Teuchos::SerialDenseMatrix<int,ScalarType> P2(numX,1);
          {
#ifdef BELOS_TEUCHOS_TIME_MONITOR
          Teuchos::TimeMonitor innerProdTimer( *timerInnerProd_ );
#endif
          MatOrthoManager<ScalarType,MV,OP>::innerProd(*prevX,*Xj,MXj,P2);
          }
          product += P2;
 
          {
#ifdef BELOS_TEUCHOS_TIME_MONITOR
          Teuchos::TimeMonitor updateTimer( *timerUpdate_ );
#endif
          MVT::MvTimesMatAddMv( -ONE, *prevX, P2, ONE, *Xj );
          }
          if ((this->_hasOp)) {
#ifdef BELOS_TEUCHOS_TIME_MONITOR
            Teuchos::TimeMonitor updateTimer( *timerUpdate_ );
#endif
            MVT::MvTimesMatAddMv( -ONE, *prevMX, P2, ONE, *MXj );
          }
        } // if (newDot[0] < dep_tol_*oldDot[0])

      } // if (numX > 0)

      // Compute Op-norm with old MXj
      if (numX > 0) {
#ifdef BELOS_TEUCHOS_TIME_MONITOR
        Teuchos::TimeMonitor normTimer( *timerNorm_ );
#endif
        MVT::MvDot( *Xj, *oldMXj, newDot );
      }
      else {
        newDot[0] = oldDot[0];
      }

      // Check to see if the new vector is dependent.
      if (completeBasis) {
        //
        // We need a complete basis, so add random vectors if necessary
        //
        if ( SCT::magnitude(newDot[0]) < SCT::magnitude(sing_tol_*oldDot[0]) ) {

          // Add a random vector and orthogonalize it against previous vectors in block.
          addVec = true;
#ifdef ORTHO_DEBUG
          std::cout << "Belos::DGKSOrthoManager::findBasis() --> Random for column " << numX << std::endl;
#endif
          //
          Teuchos::RCP<MV> tempXj = MVT::Clone( X, 1 );
          Teuchos::RCP<MV> tempMXj;
          MVT::MvRandom( *tempXj );
          if (this->_hasOp) {
            tempMXj = MVT::Clone( X, 1 );
            OPT::Apply( *(this->_Op), *tempXj, *tempMXj );
          }
          else {
            tempMXj = tempXj;
          }
          {
#ifdef BELOS_TEUCHOS_TIME_MONITOR
          Teuchos::TimeMonitor normTimer( *timerNorm_ );
#endif
          MVT::MvDot( *tempXj, *tempMXj, oldDot );
          }
          //
          for (int num_orth=0; num_orth<max_blk_ortho_; num_orth++){
            {
#ifdef BELOS_TEUCHOS_TIME_MONITOR
            Teuchos::TimeMonitor innerProdTimer( *timerInnerProd_ );
#endif
            MatOrthoManager<ScalarType,MV,OP>::innerProd(*prevX,*tempXj,tempMXj,product);
            }
            {
#ifdef BELOS_TEUCHOS_TIME_MONITOR
            Teuchos::TimeMonitor updateTimer( *timerUpdate_ );
#endif
            MVT::MvTimesMatAddMv( -ONE, *prevX, product, ONE, *tempXj );
            }
            if (this->_hasOp) {
#ifdef BELOS_TEUCHOS_TIME_MONITOR
              Teuchos::TimeMonitor updateTimer( *timerUpdate_ );
#endif
              MVT::MvTimesMatAddMv( -ONE, *prevMX, product, ONE, *tempMXj );
            }
          }
          // Compute new Op-norm
          {
#ifdef BELOS_TEUCHOS_TIME_MONITOR
          Teuchos::TimeMonitor normTimer( *timerNorm_ );
#endif
          MVT::MvDot( *tempXj, *tempMXj, newDot );
          }
          //
          if ( SCT::magnitude(newDot[0]) >= SCT::magnitude(oldDot[0]*sing_tol_) ) {
            // Copy vector into current column of _basisvecs
            MVT::Assign( *tempXj, *Xj );
            if (this->_hasOp) {
              MVT::Assign( *tempMXj, *MXj );
            }
          }
          else {
            return numX;
          }
        }
      }
      else {
        //
        // We only need to detect dependencies.
        //
        if ( SCT::magnitude(newDot[0]) < SCT::magnitude(oldDot[0]*blk_tol_) ) {
          return numX;
        }
      }

      // If we haven't left this method yet, then we can normalize the new vector Xj.
      // Normalize Xj.
      // Xj <- Xj / std::sqrt(newDot)
      ScalarType diag = SCT::squareroot(SCT::magnitude(newDot[0]));

      if (SCT::magnitude(diag) > ZERO) {
        MVT::MvScale( *Xj, ONE/diag );
        if (this->_hasOp) {
          // Update MXj.
          MVT::MvScale( *MXj, ONE/diag );
        }
      }

      // If we've added a random vector, enter a zero in the j'th diagonal element.
      if (addVec) {
        (*B)(j,j) = ZERO;
      }
      else {
        (*B)(j,j) = diag;
      }

      // Save the coefficients, if we are working on the original vector and not a randomly generated one
      if (!addVec) {
        for (int i=0; i<numX; i++) {
          (*B)(i,j) = product(i,0);
        }
      }

    } // for (j = 0; j < xc; ++j)

    return xc;
  }

  //////////////////////////////////////////////////////////////////////////////////////////////////
  // Routine to compute the block orthogonalization
  template<class ScalarType, class MV, class OP>
  bool
  DGKSOrthoManager<ScalarType, MV, OP>::blkOrtho1 ( MV &X, Teuchos::RCP<MV> MX,
                                                   Teuchos::Array<Teuchos::RCP<Teuchos::SerialDenseMatrix<int,ScalarType> > > C,
                                                   Teuchos::ArrayView<Teuchos::RCP<const MV> > Q) const
  {
    int nq = Q.size();
    int xc = MVT::GetNumberVecs( X );
    const ScalarType ONE  = SCT::one();

    std::vector<int> qcs( nq );
    for (int i=0; i<nq; i++) {
      qcs[i] = MVT::GetNumberVecs( *Q[i] );
    }

    // Compute the initial Op-norms
    std::vector<ScalarType> oldDot( 1 ), newDot( 1 );
    {
#ifdef BELOS_TEUCHOS_TIME_MONITOR
    Teuchos::TimeMonitor normTimer( *timerNorm_ );
#endif
    MVT::MvDot( X, *MX, oldDot );
    }

    Teuchos::Array<Teuchos::RCP<MV> > MQ(nq);
    // Define the product Q^T * (Op*X)
    for (int i=0; i<nq; i++) {
      // Multiply Q' with MX
      {
#ifdef BELOS_TEUCHOS_TIME_MONITOR
      Teuchos::TimeMonitor innerProdTimer( *timerInnerProd_ );
#endif
      MatOrthoManager<ScalarType,MV,OP>::innerProd(*Q[i],X,MX,*C[i]);
      }
      // Multiply by Q and subtract the result in X
      {
#ifdef BELOS_TEUCHOS_TIME_MONITOR
      Teuchos::TimeMonitor updateTimer( *timerUpdate_ );
#endif
      MVT::MvTimesMatAddMv( -ONE, *Q[i], *C[i], ONE, X );
      }

      // Update MX, with the least number of applications of Op as possible
      if (this->_hasOp) {
        if (xc <= qcs[i]) {
          OPT::Apply( *(this->_Op), X, *MX);
        }
        else {
          // this will possibly be used again below; don't delete it
          MQ[i] = MVT::Clone( *Q[i], qcs[i] );
          OPT::Apply( *(this->_Op), *Q[i], *MQ[i] );
          {
#ifdef BELOS_TEUCHOS_TIME_MONITOR
          Teuchos::TimeMonitor updateTimer( *timerUpdate_ );
#endif
          MVT::MvTimesMatAddMv( -ONE, *MQ[i], *C[i], ONE, *MX );
          }
        }
      }
    }

    {
#ifdef BELOS_TEUCHOS_TIME_MONITOR
    Teuchos::TimeMonitor normTimer( *timerNorm_ );
#endif
    MVT::MvDot( X, *MX, newDot );
    }

/*  // Compute correction bound, compare with PETSc bound.
    MagnitudeType hnrm = C[0]->normFrobenius();
    for (int i=1; i<nq; i++)
    {
      hnrm += C[i]->normFrobenius();
    }

    std::cout << "newDot < 1/sqrt(2)*oldDot < hnrm = " << MGT::squareroot(SCT::magnitude(newDot[0])) << " < " << dep_tol_*MGT::squareroot(SCT::magnitude(oldDot[0])) << " < " << hnrm << std::endl;
*/

    // Check if a correction is needed.
    if ( MGT::squareroot(SCT::magnitude(newDot[0])) < dep_tol_*MGT::squareroot(SCT::magnitude(oldDot[0])) ) {
    // Apply the second step of Gram-Schmidt

      for (int i=0; i<nq; i++) {
        Teuchos::SerialDenseMatrix<int,ScalarType> C2(C[i]->numRows(), C[i]->numCols());

        // Apply another step of classical Gram-Schmidt
        {
#ifdef BELOS_TEUCHOS_TIME_MONITOR
        Teuchos::TimeMonitor innerProdTimer( *timerInnerProd_ );
#endif
        MatOrthoManager<ScalarType,MV,OP>::innerProd(*Q[i],X,MX,C2);
        }
        *C[i] += C2;

        {
#ifdef BELOS_TEUCHOS_TIME_MONITOR
        Teuchos::TimeMonitor updateTimer( *timerUpdate_ );
#endif
        MVT::MvTimesMatAddMv( -ONE, *Q[i], C2, ONE, X );
        }

        // Update MX, with the least number of applications of Op as possible
        if (this->_hasOp) {
          if (MQ[i].get()) {
#ifdef BELOS_TEUCHOS_TIME_MONITOR
            Teuchos::TimeMonitor updateTimer( *timerUpdate_ );
#endif
            // MQ was allocated and computed above; use it
            MVT::MvTimesMatAddMv( -ONE, *MQ[i], C2, ONE, *MX );
          }
          else if (xc <= qcs[i]) {
            // MQ was not allocated and computed above; it was cheaper to use X before and it still is
            OPT::Apply( *(this->_Op), X, *MX);
          }
        }
      } // for (int i=0; i<nq; i++)
    }

    return false;
  }

  //////////////////////////////////////////////////////////////////////////////////////////////////
  // Routine to compute the block orthogonalization
  template<class ScalarType, class MV, class OP>
  bool
  DGKSOrthoManager<ScalarType, MV, OP>::blkOrtho ( MV &X, Teuchos::RCP<MV> MX,
                                                   Teuchos::Array<Teuchos::RCP<Teuchos::SerialDenseMatrix<int,ScalarType> > > C,
                                                   Teuchos::ArrayView<Teuchos::RCP<const MV> > Q) const
  {
    int nq = Q.size();
    int xc = MVT::GetNumberVecs( X );
    bool dep_flg = false;
    const ScalarType ONE  = SCT::one();

    std::vector<int> qcs( nq );
    for (int i=0; i<nq; i++) {
      qcs[i] = MVT::GetNumberVecs( *Q[i] );
    }

    // Perform the Gram-Schmidt transformation for a block of vectors

    // Compute the initial Op-norms
    std::vector<ScalarType> oldDot( xc );
    {
#ifdef BELOS_TEUCHOS_TIME_MONITOR
    Teuchos::TimeMonitor normTimer( *timerNorm_ );
#endif
    MVT::MvDot( X, *MX, oldDot );
    }

    Teuchos::Array<Teuchos::RCP<MV> > MQ(nq);
    // Define the product Q^T * (Op*X)
    for (int i=0; i<nq; i++) {
      // Multiply Q' with MX
      {
#ifdef BELOS_TEUCHOS_TIME_MONITOR
      Teuchos::TimeMonitor innerProdTimer( *timerInnerProd_ );
#endif
      MatOrthoManager<ScalarType,MV,OP>::innerProd(*Q[i],X,MX,*C[i]);
      }
      // Multiply by Q and subtract the result in X
      {
#ifdef BELOS_TEUCHOS_TIME_MONITOR
      Teuchos::TimeMonitor updateTimer( *timerUpdate_ );
#endif
      MVT::MvTimesMatAddMv( -ONE, *Q[i], *C[i], ONE, X );
      }

      // Update MX, with the least number of applications of Op as possible
      if (this->_hasOp) {
        if (xc <= qcs[i]) {
          OPT::Apply( *(this->_Op), X, *MX);
        }
        else {
          // this will possibly be used again below; don't delete it
          MQ[i] = MVT::Clone( *Q[i], qcs[i] );
          OPT::Apply( *(this->_Op), *Q[i], *MQ[i] );
          {
#ifdef BELOS_TEUCHOS_TIME_MONITOR
          Teuchos::TimeMonitor updateTimer( *timerUpdate_ );
#endif
          MVT::MvTimesMatAddMv( -ONE, *MQ[i], *C[i], ONE, *MX );
          }
        }
      }
    }

    // Do as many steps of classical Gram-Schmidt as required by max_blk_ortho_
    for (int j = 1; j < max_blk_ortho_; ++j) {

      for (int i=0; i<nq; i++) {
        Teuchos::SerialDenseMatrix<int,ScalarType> C2(C[i]->numRows(), C[i]->numCols());

        // Apply another step of classical Gram-Schmidt
        {
#ifdef BELOS_TEUCHOS_TIME_MONITOR
        Teuchos::TimeMonitor innerProdTimer( *timerInnerProd_ );
#endif
        MatOrthoManager<ScalarType,MV,OP>::innerProd(*Q[i],X,MX,C2);
        }
        *C[i] += C2;

        {
#ifdef BELOS_TEUCHOS_TIME_MONITOR
        Teuchos::TimeMonitor updateTimer( *timerUpdate_ );
#endif
        MVT::MvTimesMatAddMv( -ONE, *Q[i], C2, ONE, X );
        }

        // Update MX, with the least number of applications of Op as possible
        if (this->_hasOp) {
          if (MQ[i].get()) {
#ifdef BELOS_TEUCHOS_TIME_MONITOR
            Teuchos::TimeMonitor updateTimer( *timerUpdate_ );
#endif
            // MQ was allocated and computed above; use it
            MVT::MvTimesMatAddMv( -ONE, *MQ[i], C2, ONE, *MX );
          }
          else if (xc <= qcs[i]) {
            // MQ was not allocated and computed above; it was cheaper to use X before and it still is
            OPT::Apply( *(this->_Op), X, *MX);
          }
        }
      } // for (int i=0; i<nq; i++)
    } // for (int j = 0; j < max_blk_ortho; ++j)

    // Compute new Op-norms
    std::vector<ScalarType> newDot(xc);
    {
#ifdef BELOS_TEUCHOS_TIME_MONITOR
    Teuchos::TimeMonitor normTimer( *timerNorm_ );
#endif
    MVT::MvDot( X, *MX, newDot );
    }

    // Check to make sure the new block of vectors are not dependent on previous vectors
    for (int i=0; i<xc; i++){
      if (SCT::magnitude(newDot[i]) < SCT::magnitude(oldDot[i] * blk_tol_)) {
        dep_flg = true;
        break;
      }
    } // end for (i=0;...)

    return dep_flg;
  }


  template<class ScalarType, class MV, class OP>
  int
  DGKSOrthoManager<ScalarType, MV, OP>::blkOrthoSing ( MV &X, Teuchos::RCP<MV> MX,
                                                       Teuchos::Array<Teuchos::RCP<Teuchos::SerialDenseMatrix<int,ScalarType> > > C,
                                                       Teuchos::RCP<Teuchos::SerialDenseMatrix<int,ScalarType> > B,
                                                       Teuchos::ArrayView<Teuchos::RCP<const MV> > QQ) const
  {
    Teuchos::Array<Teuchos::RCP<const MV> > Q (QQ);

    const ScalarType ONE  = SCT::one();
    const ScalarType ZERO  = SCT::zero();

    int nq = Q.size();
    int xc = MVT::GetNumberVecs( X );
    std::vector<int> indX( 1 );
    std::vector<ScalarType> oldDot( 1 ), newDot( 1 );

    std::vector<int> qcs( nq );
    for (int i=0; i<nq; i++) {
      qcs[i] = MVT::GetNumberVecs( *Q[i] );
    }

    // Create pointers for the previous vectors of X that have already been orthonormalized.
    Teuchos::RCP<const MV> lastQ;
    Teuchos::RCP<MV> Xj, MXj;
    Teuchos::RCP<Teuchos::SerialDenseMatrix<int,ScalarType> > lastC;

    // Perform the Gram-Schmidt transformation for each vector in the block of vectors.
    for (int j=0; j<xc; j++) {

      bool dep_flg = false;

      // Get a view of the previously orthogonalized vectors and B, add it to the arrays.
      if (j > 0) {
        std::vector<int> index( j );
        for (int ind=0; ind<j; ind++) {
          index[ind] = ind;
        }
        lastQ = MVT::CloneView( X, index );

        // Add these views to the Q and C arrays.
        Q.push_back( lastQ );
        C.push_back( B );
        qcs.push_back( MVT::GetNumberVecs( *lastQ ) );
      }

      // Get a view of the current vector in X to orthogonalize.
      indX[0] = j;
      Xj = MVT::CloneViewNonConst( X, indX );
      if (this->_hasOp) {
        MXj = MVT::CloneViewNonConst( *MX, indX );
      }
      else {
        MXj = Xj;
      }

      // Compute the initial Op-norms
      {
#ifdef BELOS_TEUCHOS_TIME_MONITOR
      Teuchos::TimeMonitor normTimer( *timerNorm_ );
#endif
      MVT::MvDot( *Xj, *MXj, oldDot );
      }

      Teuchos::Array<Teuchos::RCP<MV> > MQ(Q.size());
      // Define the product Q^T * (Op*X)
      for (int i=0; i<Q.size(); i++) {

        // Get a view of the current serial dense matrix
        Teuchos::SerialDenseMatrix<int,ScalarType> tempC( Teuchos::View, *C[i], qcs[i], 1, 0, j );

        // Multiply Q' with MX
        {
#ifdef BELOS_TEUCHOS_TIME_MONITOR
        Teuchos::TimeMonitor innerProdTimer( *timerInnerProd_ );
#endif
        MatOrthoManager<ScalarType,MV,OP>::innerProd(*Q[i],*Xj,MXj,tempC);
        }
        // Multiply by Q and subtract the result in Xj
        {
#ifdef BELOS_TEUCHOS_TIME_MONITOR
        Teuchos::TimeMonitor updateTimer( *timerUpdate_ );
#endif
        MVT::MvTimesMatAddMv( -ONE, *Q[i], tempC, ONE, *Xj );
        }

        // Update MXj, with the least number of applications of Op as possible
        if (this->_hasOp) {
          if (xc <= qcs[i]) {
            OPT::Apply( *(this->_Op), *Xj, *MXj);
          }
          else {
            // this will possibly be used again below; don't delete it
            MQ[i] = MVT::Clone( *Q[i], qcs[i] );
            OPT::Apply( *(this->_Op), *Q[i], *MQ[i] );
            {
#ifdef BELOS_TEUCHOS_TIME_MONITOR
            Teuchos::TimeMonitor updateTimer( *timerUpdate_ );
#endif
            MVT::MvTimesMatAddMv( -ONE, *MQ[i], tempC, ONE, *MXj );
            }
          }
        }
      }

      // Compute the Op-norms
      {
#ifdef BELOS_TEUCHOS_TIME_MONITOR
      Teuchos::TimeMonitor normTimer( *timerNorm_ );
#endif
      MVT::MvDot( *Xj, *MXj, newDot );
      }

      // Do one step of classical Gram-Schmidt orthogonalization
      // with a second correction step if needed.

      if ( SCT::magnitude(newDot[0]) < SCT::magnitude(oldDot[0]*dep_tol_) ) {

        for (int i=0; i<Q.size(); i++) {
          Teuchos::SerialDenseMatrix<int,ScalarType> tempC( Teuchos::View, *C[i], qcs[i], 1, 0, j );
          Teuchos::SerialDenseMatrix<int,ScalarType> C2( qcs[i], 1 );

          // Apply another step of classical Gram-Schmidt
          {
#ifdef BELOS_TEUCHOS_TIME_MONITOR
          Teuchos::TimeMonitor innerProdTimer( *timerInnerProd_ );
#endif
          MatOrthoManager<ScalarType,MV,OP>::innerProd(*Q[i],*Xj,MXj,C2);
          }
          tempC += C2;
          {
#ifdef BELOS_TEUCHOS_TIME_MONITOR
          Teuchos::TimeMonitor updateTimer( *timerUpdate_ );
#endif
          MVT::MvTimesMatAddMv( -ONE, *Q[i], C2, ONE, *Xj );
          }

          // Update MXj, with the least number of applications of Op as possible
          if (this->_hasOp) {
            if (MQ[i].get()) {
#ifdef BELOS_TEUCHOS_TIME_MONITOR
              Teuchos::TimeMonitor updateTimer( *timerUpdate_ );
#endif
              // MQ was allocated and computed above; use it
              MVT::MvTimesMatAddMv( -ONE, *MQ[i], C2, ONE, *MXj );
            }
            else if (xc <= qcs[i]) {
              // MQ was not allocated and computed above; it was cheaper to use X before and it still is
              OPT::Apply( *(this->_Op), *Xj, *MXj);
            }
          }
        } // for (int i=0; i<Q.size(); i++)

        // Compute the Op-norms after the correction step.
        {
#ifdef BELOS_TEUCHOS_TIME_MONITOR
        Teuchos::TimeMonitor normTimer( *timerNorm_ );
#endif
        MVT::MvDot( *Xj, *MXj, newDot );
        }
      } // if ()

      // Check for linear dependence.
      if (SCT::magnitude(newDot[0]) < SCT::magnitude(oldDot[0]*sing_tol_)) {
        dep_flg = true;
      }

      // Normalize the new vector if it's not dependent
      if (!dep_flg) {
        ScalarType diag = SCT::squareroot(SCT::magnitude(newDot[0]));

        MVT::MvScale( *Xj, ONE/diag );
        if (this->_hasOp) {
          // Update MXj.
          MVT::MvScale( *MXj, ONE/diag );
        }

        // Enter value on diagonal of B.
        (*B)(j,j) = diag;
      }
      else {
        // Create a random vector and orthogonalize it against all previous columns of Q.
        Teuchos::RCP<MV> tempXj = MVT::Clone( X, 1 );
        Teuchos::RCP<MV> tempMXj;
        MVT::MvRandom( *tempXj );
        if (this->_hasOp) {
          tempMXj = MVT::Clone( X, 1 );
          OPT::Apply( *(this->_Op), *tempXj, *tempMXj );
        }
        else {
          tempMXj = tempXj;
        }
        {
#ifdef BELOS_TEUCHOS_TIME_MONITOR
        Teuchos::TimeMonitor normTimer( *timerNorm_ );
#endif
        MVT::MvDot( *tempXj, *tempMXj, oldDot );
        }
        //
        for (int num_orth=0; num_orth<max_blk_ortho_; num_orth++) {

          for (int i=0; i<Q.size(); i++) {
            Teuchos::SerialDenseMatrix<int,ScalarType> product( qcs[i], 1 );

            // Apply another step of classical Gram-Schmidt
            {
#ifdef BELOS_TEUCHOS_TIME_MONITOR
            Teuchos::TimeMonitor innerProdTimer( *timerInnerProd_ );
#endif
            MatOrthoManager<ScalarType,MV,OP>::innerProd(*Q[i],*tempXj,tempMXj,product);
            }
            {
#ifdef BELOS_TEUCHOS_TIME_MONITOR
            Teuchos::TimeMonitor updateTimer( *timerUpdate_ );
#endif
            MVT::MvTimesMatAddMv( -ONE, *Q[i], product, ONE, *tempXj );
            }

            // Update MXj, with the least number of applications of Op as possible
            if (this->_hasOp) {
              if (MQ[i].get()) {
#ifdef BELOS_TEUCHOS_TIME_MONITOR
                Teuchos::TimeMonitor updateTimer( *timerUpdate_ );
#endif
                // MQ was allocated and computed above; use it
                MVT::MvTimesMatAddMv( -ONE, *MQ[i], product, ONE, *tempMXj );
              }
              else if (xc <= qcs[i]) {
                // MQ was not allocated and computed above; it was cheaper to use X before and it still is
                OPT::Apply( *(this->_Op), *tempXj, *tempMXj);
              }
            }
          } // for (int i=0; i<nq; i++)

        }

        // Compute the Op-norms after the correction step.
        {
#ifdef BELOS_TEUCHOS_TIME_MONITOR
        Teuchos::TimeMonitor normTimer( *timerNorm_ );
#endif
        MVT::MvDot( *tempXj, *tempMXj, newDot );
        }

        // Copy vector into current column of Xj
        if ( SCT::magnitude(newDot[0]) >= SCT::magnitude(oldDot[0]*sing_tol_) ) {
          ScalarType diag = SCT::squareroot(SCT::magnitude(newDot[0]));

          // Enter value on diagonal of B.
          (*B)(j,j) = ZERO;

          // Copy vector into current column of _basisvecs
          MVT::MvAddMv( ONE/diag, *tempXj, ZERO, *tempXj, *Xj );
          if (this->_hasOp) {
            MVT::MvAddMv( ONE/diag, *tempMXj, ZERO, *tempMXj, *MXj );
          }
        }
        else {
          return j;
        }
      } // if (!dep_flg)

      // Remove the vectors from array
      if (j > 0) {
        Q.resize( nq );
        C.resize( nq );
        qcs.resize( nq );
      }

    } // for (int j=0; j<xc; j++)

    return xc;
  }

  template<class ScalarType, class MV, class OP>
  Teuchos::RCP<Teuchos::ParameterList> getDGKSDefaultParameters ()
  {
    using Teuchos::ParameterList;
    using Teuchos::parameterList;
    using Teuchos::RCP;

    RCP<ParameterList> params = parameterList ("DGKS");

    // Default parameter values for DGKS orthogonalization.
    // Documentation will be embedded in the parameter list.
    params->set ("maxNumOrthogPasses", DGKSOrthoManager<ScalarType, MV, OP>::max_blk_ortho_default_,
                 "Maximum number of orthogonalization passes (includes the "
                 "first).  Default is 2, since \"twice is enough\" for Krylov "
                 "methods.");
    params->set ("blkTol", DGKSOrthoManager<ScalarType, MV, OP>::blk_tol_default_, 
                 "Block reorthogonalization threshold.");
    params->set ("depTol", DGKSOrthoManager<ScalarType, MV, OP>::dep_tol_default_,
                 "(Non-block) reorthogonalization threshold.");
    params->set ("singTol", DGKSOrthoManager<ScalarType, MV, OP>::sing_tol_default_, 
                 "Singular block detection threshold.");

    return params;
  }

  template<class ScalarType, class MV, class OP>
  Teuchos::RCP<Teuchos::ParameterList> getDGKSFastParameters ()
  {
    using Teuchos::ParameterList;
    using Teuchos::RCP;

    RCP<ParameterList> params = getDGKSDefaultParameters<ScalarType, MV, OP>();

    params->set ("maxNumOrthogPasses", 
                 DGKSOrthoManager<ScalarType, MV, OP>::max_blk_ortho_fast_);
    params->set ("blkTol", 
                 DGKSOrthoManager<ScalarType, MV, OP>::blk_tol_fast_);
    params->set ("depTol", 
                 DGKSOrthoManager<ScalarType, MV, OP>::dep_tol_fast_);
    params->set ("singTol", 
                 DGKSOrthoManager<ScalarType, MV, OP>::sing_tol_fast_);

    return params;
  }

} // namespace Belos

#endif // BELOS_DGKS_ORTHOMANAGER_HPP

