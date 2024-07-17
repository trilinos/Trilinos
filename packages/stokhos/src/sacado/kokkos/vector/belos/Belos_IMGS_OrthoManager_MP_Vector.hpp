// @HEADER
// *****************************************************************************
//                 Belos: Block Linear Solvers Package
//
// Copyright 2004 NTESS and the Belos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/*! 
  \file Belos_IMGS_OrthoManager_MP_Vector.hpp
  \brief Iterated Modified Gram-Schmidt (IMGS) implementation of the Belos::OrthoManager class

  Specialized for Sacado::MP::Vector scalar type to deal with lucky breakdown that can occur for one sample
  but not for the others.    
*/

#ifndef BELOS_IMGS_ORTHOMANAGER_MP_VECTOR_HPP
#define BELOS_IMGS_ORTHOMANAGER_MP_VECTOR_HPP

#include "Stokhos_Sacado_Kokkos_MP_Vector.hpp"
#include "BelosIMGSOrthoManager.hpp"

namespace Belos {

  template<class Storage, class MV, class OP>
  class IMGSOrthoManager<Sacado::MP::Vector<Storage>,MV,OP> :
    public MatOrthoManager<Sacado::MP::Vector<Storage>,MV,OP>
  {
  private:
    typedef Sacado::MP::Vector<Storage> ScalarType;
    typedef typename Teuchos::ScalarTraits<ScalarType>::magnitudeType MagnitudeType;
    typedef typename Teuchos::ScalarTraits<MagnitudeType> MGT;
    typedef Teuchos::ScalarTraits<ScalarType>  SCT;
    typedef MultiVecTraits<ScalarType,MV>      MVT;
    typedef OperatorTraits<ScalarType,MV,OP>   OPT;

  public:
    //! @name Constructor/Destructor
    //@{

    //! Constructor specifying re-orthogonalization tolerance.
    IMGSOrthoManager( const std::string& label = "Belos",
                      Teuchos::RCP<const OP> Op = Teuchos::null,
                      const int max_ortho_steps = max_ortho_steps_default_,
                      const MagnitudeType blk_tol = blk_tol_default_,
                      const MagnitudeType sing_tol = sing_tol_default_ )
      : MatOrthoManager<ScalarType,MV,OP>(Op),
        max_ortho_steps_( max_ortho_steps ),
        blk_tol_( blk_tol ),
        sing_tol_( sing_tol ),
        label_( label )
    {
#ifdef BELOS_TEUCHOS_TIME_MONITOR
        std::stringstream ss;
        ss << label_ + ": IMGS[" << max_ortho_steps_ << "]";

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
    IMGSOrthoManager (const Teuchos::RCP<Teuchos::ParameterList>& plist,
                      const std::string& label = "Belos",
                      Teuchos::RCP<const OP> Op = Teuchos::null) :
      MatOrthoManager<ScalarType,MV,OP>(Op),
      max_ortho_steps_ (max_ortho_steps_default_),
      blk_tol_ (blk_tol_default_),
      sing_tol_ (sing_tol_default_),
      label_ (label)
    {
      setParameterList (plist);

#ifdef BELOS_TEUCHOS_TIME_MONITOR
      std::stringstream ss;
      ss << label_ + ": IMGS[" << max_ortho_steps_ << "]";

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
    ~IMGSOrthoManager() {}
    //@}

    //! @name Implementation of Teuchos::ParameterListAcceptorDefaultBase interface
    //@{
    void
    setParameterList (const Teuchos::RCP<Teuchos::ParameterList>& plist)
    {
      using Teuchos::Exceptions::InvalidParameterName;
      using Teuchos::ParameterList;
      using Teuchos::parameterList;
      using Teuchos::RCP;

      RCP<const ParameterList> defaultParams = getValidParameters();
      RCP<ParameterList> params;
      if (plist.is_null()) {
        params = parameterList (*defaultParams);
      } else {
        params = plist;
        // Some users might want to specify "blkTol" as "depTol".  Due
        // to this case, we don't invoke
        // validateParametersAndSetDefaults on params.  Instead, we go
        // through the parameter list one parameter at a time and look
        // for alternatives.
      }

      // Using temporary variables and fetching all values before
      // setting the output arguments ensures the strong exception
      // guarantee for this function: if an exception is thrown, no
      // externally visible side effects (in this case, setting the
      // output arguments) have taken place.
      int maxNumOrthogPasses;
      MagnitudeType blkTol;
      MagnitudeType singTol;

      try {
        maxNumOrthogPasses = params->get<int> ("maxNumOrthogPasses");
      } catch (InvalidParameterName&) {
        maxNumOrthogPasses = defaultParams->get<int> ("maxNumOrthogPasses");
        params->set ("maxNumOrthogPasses", maxNumOrthogPasses);
      }

      // Handling of the "blkTol" parameter is a special case.  This
      // is because some users may prefer to call this parameter
      // "depTol" for consistency with DGKS.  However, our default
      // parameter list calls this "blkTol", and we don't want the
      // default list's value to override the user's value.  Thus, we
      // first check the user's parameter list for both names, and
      // only then access the default parameter list.
      try {
        blkTol = params->get<MagnitudeType> ("blkTol");
      } catch (InvalidParameterName&) {
        try {
          blkTol = params->get<MagnitudeType> ("depTol");
          // "depTol" is the wrong name, so remove it and replace with
          // "blkTol".  We'll set "blkTol" below.
          params->remove ("depTol");
        } catch (InvalidParameterName&) {
          blkTol = defaultParams->get<MagnitudeType> ("blkTol");
        }
        params->set ("blkTol", blkTol);
      }

      try {
        singTol = params->get<MagnitudeType> ("singTol");
      } catch (InvalidParameterName&) {
        singTol = defaultParams->get<MagnitudeType> ("singTol");
        params->set ("singTol", singTol);
      }

      max_ortho_steps_ = maxNumOrthogPasses;
      blk_tol_ = blkTol;
      sing_tol_ = singTol;

      this->setMyParamList (params);
    }

    Teuchos::RCP<const Teuchos::ParameterList>
    getValidParameters () const
    {
      if (defaultParams_.is_null()) {
        defaultParams_ = Belos::getIMGSDefaultParameters<ScalarType, MV, OP>();
      }

      return defaultParams_;
    }

    //@}

    //! @name Accessor routines
    //@{

    //! Set parameter for block re-orthogonalization threshhold.
    void setBlkTol( const MagnitudeType blk_tol ) { blk_tol_ = blk_tol; }

    //! Set parameter for singular block detection.
    void setSingTol( const MagnitudeType sing_tol ) { sing_tol_ = sing_tol; }

    //! Return parameter for block re-orthogonalization threshhold.
    MagnitudeType getBlkTol() const { return blk_tol_; }

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
     * The method uses either one or two steps of modified Gram-Schmidt. The algebraically
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
     * The method uses modified Gram-Schmidt, so that the coefficient matrix \c B is upper triangular.
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
     *  This routine returns an integer \c rank stating the rank of the computed basis. If the subspace \f$colspan(X) - \sum_i colspan(Q[i])\f$ does not
     *  have dimension as large as the number of columns of \c X and the orthogonalization manager doe not attempt to augment the subspace, then \c rank
     *  may be smaller than the number of columns of \c X. In this case, only the first \c rank columns of output \c X and first \c rank rows of \c B will
     *  be valid.
     *
     * The method attempts to find a basis with dimension the same as the number of columns in \c X. It does this by augmenting linearly dependant
     * vectors with random directions. A finite number of these attempts will be made; therefore, it is possible that the dimension of the
     * computed basis is less than the number of vectors in \c X.
     *
     @param X [in/out] The multivector to the modified.
       On output, the relevant rows of \c X will be orthogonal to the <tt>Q[i]</tt> and will have orthonormal columns (with respect to innerProd()).

     @param MX [in/out] The image of \c X under the operator \c Op.
       If \f$ MX != 0\f$: On input, this is expected to be consistent with \c X. On output, this is updated consistent with updates to \c X.
       If \f$ MX == 0\f$ or \f$ Op == 0\f$: \c MX is not referenced.

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

     @param B [out] The coefficients of the original \c X with respect to the computed basis. The first rows in \c B
            corresponding to the valid columns in \c X will be upper triangular.

     @param Q [in] A list of multivector bases specifying the subspaces to be orthogonalized against. Each <tt>Q[i]</tt> is assumed to have
     orthonormal columns, and the <tt>Q[i]</tt> are assumed to be mutually orthogonal.

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

    /*! \brief This method computes the error in orthonormality of a multivector, measured
     * as the Frobenius norm of the difference <tt>innerProd(X,Y) - I</tt>.
     */
    typename Teuchos::ScalarTraits<ScalarType>::magnitudeType
    orthonormError(const MV &X) const {
      return orthonormError(X,Teuchos::null);
    }

    /*! \brief This method computes the error in orthonormality of a multivector, measured
     * as the Frobenius norm of the difference <tt>innerProd(X,Y) - I</tt>.
     *  The method has the option of exploiting a caller-provided \c MX.
     */
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
    static const int max_ortho_steps_default_;
    //! Block reorthogonalization threshold (default).
    static const MagnitudeType blk_tol_default_;
    //! Singular block detection threshold (default).
    static const MagnitudeType sing_tol_default_;

    //! Max number of (re)orthogonalization steps, including the first (fast).
    static const int max_ortho_steps_fast_;
    //! Block reorthogonalization threshold (fast).
    static const MagnitudeType blk_tol_fast_;
    //! Singular block detection threshold (fast).
    static const MagnitudeType sing_tol_fast_;

    //@}

  private:

    //! Max number of (re)orthogonalization steps, including the first.
    int max_ortho_steps_;
    //! Block reorthogonalization tolerance.
    MagnitudeType blk_tol_;
    //! Singular block detection threshold.
    MagnitudeType sing_tol_;

    //! Label for timers.
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
  template<class StorageType, class MV, class OP>
  const int IMGSOrthoManager<Sacado::MP::Vector<StorageType>,MV,OP>::max_ortho_steps_default_ = 1;

  template<class StorageType, class MV, class OP>
  const typename IMGSOrthoManager<Sacado::MP::Vector<StorageType>,MV,OP>::MagnitudeType
  IMGSOrthoManager<Sacado::MP::Vector<StorageType>,MV,OP>::blk_tol_default_
    = 10*Teuchos::ScalarTraits<typename IMGSOrthoManager<Sacado::MP::Vector<StorageType>,MV,OP>::MagnitudeType>::squareroot(
      Teuchos::ScalarTraits<typename IMGSOrthoManager<Sacado::MP::Vector<StorageType>,MV,OP>::MagnitudeType>::eps() );

  template<class StorageType, class MV, class OP>
  const typename IMGSOrthoManager<Sacado::MP::Vector<StorageType>,MV,OP>::MagnitudeType
  IMGSOrthoManager<Sacado::MP::Vector<StorageType>,MV,OP>::sing_tol_default_
    = 10*Teuchos::ScalarTraits<typename IMGSOrthoManager<Sacado::MP::Vector<StorageType>,MV,OP>::MagnitudeType>::eps();

  template<class StorageType, class MV, class OP>
  const int IMGSOrthoManager<Sacado::MP::Vector<StorageType>,MV,OP>::max_ortho_steps_fast_ = 1;

  template<class StorageType, class MV, class OP>
  const typename IMGSOrthoManager<Sacado::MP::Vector<StorageType>,MV,OP>::MagnitudeType
  IMGSOrthoManager<Sacado::MP::Vector<StorageType>,MV,OP>::blk_tol_fast_
    = Teuchos::ScalarTraits<typename IMGSOrthoManager<Sacado::MP::Vector<StorageType>,MV,OP>::MagnitudeType>::zero();

  template<class StorageType, class MV, class OP>
  const typename IMGSOrthoManager<Sacado::MP::Vector<StorageType>,MV,OP>::MagnitudeType
  IMGSOrthoManager<Sacado::MP::Vector<StorageType>,MV,OP>::sing_tol_fast_
    = Teuchos::ScalarTraits<typename IMGSOrthoManager<Sacado::MP::Vector<StorageType>,MV,OP>::MagnitudeType>::zero();

  //////////////////////////////////////////////////////////////////////////////////////////////////
  // Set the label for this orthogonalization manager and create new timers if it's changed
  template<class StorageType, class MV, class OP>
  void IMGSOrthoManager<Sacado::MP::Vector<StorageType>,MV,OP>::setLabel(const std::string& label)
  {
    if (label != label_) {
      label_ = label;
#ifdef BELOS_TEUCHOS_TIME_MONITOR
      std::stringstream ss;
      ss << label_ + ": IMGS[" << max_ortho_steps_ << "]";

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
  template<class StorageType, class MV, class OP>
  typename Teuchos::ScalarTraits<Sacado::MP::Vector<StorageType> >::magnitudeType
  IMGSOrthoManager<Sacado::MP::Vector<StorageType>,MV,OP>::orthonormError(const MV &X, Teuchos::RCP<const MV> MX) const {
    const ScalarType ONE = SCT::one();
    int rank = MVT::GetNumberVecs(X);
    Teuchos::SerialDenseMatrix<int,ScalarType> xTx(rank,rank);
    MatOrthoManager<ScalarType,MV,OP>::innerProd(X,X,MX,xTx);
    for (int i=0; i<rank; i++) {
      xTx(i,i) -= ONE;
    }
    return xTx.normFrobenius();
  }

  //////////////////////////////////////////////////////////////////////////////////////////////////
  // Compute the distance from orthogonality
  template<class StorageType, class MV, class OP>
  typename Teuchos::ScalarTraits<Sacado::MP::Vector<StorageType> >::magnitudeType
  IMGSOrthoManager<Sacado::MP::Vector<StorageType>,MV,OP>::orthogError(const MV &X1, Teuchos::RCP<const MV> MX1, const MV &X2) const {
    int r1 = MVT::GetNumberVecs(X1);
    int r2  = MVT::GetNumberVecs(X2);
    Teuchos::SerialDenseMatrix<int,ScalarType> xTx(r2,r1);
    MatOrthoManager<ScalarType,MV,OP>::innerProd(X2,X1,MX1,xTx);
    return xTx.normFrobenius();
  }

  //////////////////////////////////////////////////////////////////////////////////////////////////
  // Find an Op-orthonormal basis for span(X) - span(W)
  template<class StorageType, class MV, class OP>
  int
  IMGSOrthoManager<Sacado::MP::Vector<StorageType>, MV, OP>::
  projectAndNormalizeWithMxImpl(MV &X,
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
    const MagnitudeType ZERO = MGT::zero();

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
                               "IMGS orthogonalization: failed to reshape "
                               "C[" << k << "] (the array of block "
                               "coefficients resulting from projecting X "
                               "against Q[1:" << nq << "]).");
          }
      }

    /******   DO NOT MODIFY *MX IF _hasOp == false   ******/
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
    TEUCHOS_TEST_FOR_EXCEPTION( xc == 0 || xr == 0, std::invalid_argument, "Belos::IMGSOrthoManager::projectAndNormalize(): X must be non-empty" );

    int numbas = 0;
    for (int i=0; i<nq; i++) {
      numbas += MVT::GetNumberVecs( *Q[i] );
    }

    // check size of B
    TEUCHOS_TEST_FOR_EXCEPTION( B->numRows() != xc || B->numCols() != xc, std::invalid_argument,
                        "Belos::IMGSOrthoManager::projectAndNormalize(): Size of X must be consistant with size of B" );
    // check size of X and MX
    TEUCHOS_TEST_FOR_EXCEPTION( xc<0 || xr<0 || mxc<0 || mxr<0, std::invalid_argument,
                        "Belos::IMGSOrthoManager::projectAndNormalize(): MVT returned negative dimensions for X,MX" );
    // check size of X w.r.t. MX
    TEUCHOS_TEST_FOR_EXCEPTION( xc!=mxc || xr!=mxr, std::invalid_argument,
                        "Belos::IMGSOrthoManager::projectAndNormalize(): Size of X must be consistant with size of MX" );
    // check feasibility
    //TEUCHOS_TEST_FOR_EXCEPTION( numbas+xc > xr, std::invalid_argument,
    //                    "Belos::IMGSOrthoManager::projectAndNormalize(): Orthogonality constraints not feasible" );

    // Some flags for checking dependency returns from the internal orthogonalization methods
    bool dep_flg = false;

    // Make a temporary copy of X and MX, just in case a block dependency is detected.
    Teuchos::RCP<MV> tmpX, tmpMX;
    tmpX = MVT::CloneCopy(X);
    if (this->_hasOp) {
      tmpMX = MVT::CloneCopy(*MX);
    }

    if (xc == 1) {

      // Use the cheaper block orthogonalization.
      // NOTE: Don't check for dependencies because the update has one vector.
      dep_flg = blkOrtho1( X, MX, C, Q );

      // Normalize the new block X
      if ( B == Teuchos::null ) {
        B = Teuchos::rcp( new Teuchos::SerialDenseMatrix<int,ScalarType>(xc,xc) );
      }
      std::vector<ScalarType> diag(xc);
      {
#ifdef BELOS_TEUCHOS_TIME_MONITOR
        Teuchos::TimeMonitor normTimer( *timerNorm_ );
#endif
        MVT::MvDot( X, *MX, diag );
      }
      (*B)(0,0) = SCT::squareroot(SCT::magnitude(diag[0]));

      ScalarType scale = ONE;
      mask_assign((*B)(0,0)!= ZERO, scale) /= (*B)(0,0);

      if(MaskLogic::OR((*B)(0,0)!= ZERO) )
        rank = 1;
      MVT::MvScale( X, scale );
      if (this->_hasOp) {
        // Update MXj.
        MVT::MvScale( *MX, scale );
      }
    }
    else {

      // Use the cheaper block orthogonalization.
      dep_flg = blkOrtho( X, MX, C, Q );

      // If a dependency has been detected in this block, then perform
      // the more expensive nonblock (single vector at a time)
      // orthogonalization.
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
          // rerun orthogonalization using more expensive single-
          // vector orthogonalization.
          rank = blkOrthoSing( *tmpX, tmpMX, C, B, Q );

          // Copy tmpX back into X.
          MVT::Assign( *tmpX, X );
          if (this->_hasOp) {
            MVT::Assign( *tmpMX, *MX );
          }
        }
      }
    } // if (xc == 1) {

    // this should not raise an std::exception; but our post-conditions oblige us to check
    TEUCHOS_TEST_FOR_EXCEPTION( rank > xc || rank < 0, std::logic_error,
                        "Belos::IMGSOrthoManager::projectAndNormalize(): Debug error in rank variable." );

    // Return the rank of X.
    return rank;
  }



  //////////////////////////////////////////////////////////////////////////////////////////////////
  // Find an Op-orthonormal basis for span(X), with rank numvectors(X)
  template<class StorageType, class MV, class OP>
  int IMGSOrthoManager<Sacado::MP::Vector<StorageType>, MV, OP>::normalize(
                                MV &X, Teuchos::RCP<MV> MX,
                                Teuchos::RCP<Teuchos::SerialDenseMatrix<int,ScalarType> > B ) const {

#ifdef BELOS_TEUCHOS_TIME_MONITOR
    Teuchos::TimeMonitor orthotimer(*timerOrtho_);
#endif

    // call findBasis, with the instruction to try to generate a basis of rank numvecs(X)
    return findBasis(X, MX, B, true);
  }



  //////////////////////////////////////////////////////////////////////////////////////////////////
  template<class StorageType, class MV, class OP>
  void IMGSOrthoManager<Sacado::MP::Vector<StorageType>, MV, OP>::project(
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
    // MX : Image of the block of vectors X by the mass matrix
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


    /******   DO NOT MODIFY *MX IF _hasOp == false   ******/
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
                        "Belos::IMGSOrthoManager::project(): MVT returned negative dimensions for X,MX" );
    // check size of X w.r.t. MX and Q
    TEUCHOS_TEST_FOR_EXCEPTION( xc!=mxc || xr!=mxr || xr!=qr, std::invalid_argument,
                        "Belos::IMGSOrthoManager::project(): Size of X not consistant with MX,Q" );

    // tally up size of all Q and check/allocate C
    int baslen = 0;
    for (int i=0; i<nq; i++) {
      TEUCHOS_TEST_FOR_EXCEPTION( MVT::GetGlobalLength( *Q[i] ) != qr, std::invalid_argument,
                          "Belos::IMGSOrthoManager::project(): Q lengths not mutually consistant" );
      qcs[i] = MVT::GetNumberVecs( *Q[i] );
      TEUCHOS_TEST_FOR_EXCEPTION( qr < qcs[i], std::invalid_argument,
                          "Belos::IMGSOrthoManager::project(): Q has less rows than columns" );
      baslen += qcs[i];

      // check size of C[i]
      if ( C[i] == Teuchos::null ) {
        C[i] = Teuchos::rcp( new Teuchos::SerialDenseMatrix<int,ScalarType>(qcs[i],xc) );
      }
      else {
        TEUCHOS_TEST_FOR_EXCEPTION( C[i]->numRows() != qcs[i] || C[i]->numCols() != xc , std::invalid_argument,
                           "Belos::IMGSOrthoManager::project(): Size of Q not consistant with size of C" );
      }
    }

    // Use the cheaper block orthogonalization, don't check for rank deficiency.
    blkOrtho( X, MX, C, Q );

  }

  //////////////////////////////////////////////////////////////////////////////////////////////////
  // Find an Op-orthonormal basis for span(X), with the option of extending the subspace so that
  // the rank is numvectors(X)
  template<class StorageType, class MV, class OP>
  int IMGSOrthoManager<Sacado::MP::Vector<StorageType>, MV, OP>::findBasis(
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
                        "Belos::IMGSOrthoManager::findBasis(): X must be non-empty" );
    TEUCHOS_TEST_FOR_EXCEPTION( B->numRows() != xc || B->numCols() != xc, std::invalid_argument,
                        "Belos::IMGSOrthoManager::findBasis(): Size of X not consistant with size of B" );
    TEUCHOS_TEST_FOR_EXCEPTION( xc != mxc || xr != mxr, std::invalid_argument,
                        "Belos::IMGSOrthoManager::findBasis(): Size of X not consistant with size of MX" );
    TEUCHOS_TEST_FOR_EXCEPTION( xc > xr, std::invalid_argument,
                        "Belos::IMGSOrthoManager::findBasis(): Size of X not feasible for normalization" );
    TEUCHOS_TEST_FOR_EXCEPTION( howMany < 0 || howMany > xc, std::invalid_argument,
                        "Belos::IMGSOrthoManager::findBasis(): Invalid howMany parameter" );

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

      Teuchos::RCP<MV> oldMXj; 
      if (numX > 0) {
        oldMXj = MVT::CloneCopy( *MXj );
      }
 
      // Make storage for these Gram-Schmidt iterations.
      Teuchos::SerialDenseVector<int,ScalarType> product(numX);
      Teuchos::SerialDenseVector<int,ScalarType> P2(1);
      Teuchos::RCP<const MV> prevX, prevMX;

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
                          "Belos::IMGSOrthoManager::findBasis(): Negative definiteness discovered in inner product" );

      // Perform MGS one vector at a time
      for (int ii=0; ii<numX; ii++) {

        index[0] = ii;
        prevX = MVT::CloneView( X, index );
        if (this->_hasOp) {
          prevMX = MVT::CloneView( *MX, index );
        }

        for (int i=0; i<max_ortho_steps_; ++i) {

          // product <- prevX^T MXj
          {
#ifdef BELOS_TEUCHOS_TIME_MONITOR
            Teuchos::TimeMonitor innerProdTimer( *timerInnerProd_ );
#endif
            MatOrthoManager<ScalarType,MV,OP>::innerProd(*prevX,*Xj,MXj,P2);
          }

          // Xj <- Xj - prevX prevX^T MXj
          //     = Xj - prevX product
          {
#ifdef BELOS_TEUCHOS_TIME_MONITOR
            Teuchos::TimeMonitor updateTimer( *timerUpdate_ );
#endif
            MVT::MvTimesMatAddMv( -ONE, *prevX, P2, ONE, *Xj );
          }

          // Update MXj
          if (this->_hasOp) {
            // MXj <- Op*Xj_new
            //      = Op*(Xj_old - prevX prevX^T MXj)
            //      = MXj - prevMX product
#ifdef BELOS_TEUCHOS_TIME_MONITOR
            Teuchos::TimeMonitor updateTimer( *timerUpdate_ );
#endif
            MVT::MvTimesMatAddMv( -ONE, *prevMX, P2, ONE, *MXj );
          }

          // Set coefficients
          if ( i==0 )
            product[ii] = P2[0];
          else
            product[ii] += P2[0];

        } // for (int i=0; i<max_ortho_steps_; ++i)

      } // for (int ii=0; ii<numX; ++ii)

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
          std::cout << "Belos::IMGSOrthoManager::findBasis() --> Random for column " << numX << std::endl;
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
          // Perform MGS one vector at a time
          for (int ii=0; ii<numX; ii++) {

            index[0] = ii;
            prevX = MVT::CloneView( X, index );
            if (this->_hasOp) {
              prevMX = MVT::CloneView( *MX, index );
            }

            for (int num_orth=0; num_orth<max_ortho_steps_; num_orth++){
              {
#ifdef BELOS_TEUCHOS_TIME_MONITOR
                Teuchos::TimeMonitor innerProdTimer( *timerInnerProd_ );
#endif
                MatOrthoManager<ScalarType,MV,OP>::innerProd(*prevX,*tempXj,tempMXj,P2);
              }
              {
#ifdef BELOS_TEUCHOS_TIME_MONITOR
                Teuchos::TimeMonitor updateTimer( *timerUpdate_ );
#endif
                MVT::MvTimesMatAddMv( -ONE, *prevX, P2, ONE, *tempXj );
              }
              if (this->_hasOp) {
#ifdef BELOS_TEUCHOS_TIME_MONITOR
                Teuchos::TimeMonitor updateTimer( *timerUpdate_ );
#endif
                MVT::MvTimesMatAddMv( -ONE, *prevMX, P2, ONE, *tempMXj );
              }

              // Set coefficients
              if ( num_orth==0 )
                product[ii] = P2[0];
              else
                product[ii] += P2[0];
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
      ScalarType scale = ONE;
      mask_assign(SCT::magnitude(diag) > ZERO, scale) /= diag;
      MVT::MvScale( *Xj, scale );
      if (this->_hasOp) {
        // Update MXj.
        MVT::MvScale( *MXj, scale );
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
          (*B)(i,j) = product(i);
        }
      }

    } // for (j = 0; j < xc; ++j)

    return xc;
  }

  //////////////////////////////////////////////////////////////////////////////////////////////////
  // Routine to compute the block orthogonalization
  template<class StorageType, class MV, class OP>
  bool
  IMGSOrthoManager<Sacado::MP::Vector<StorageType>, MV, OP>::blkOrtho1 ( MV &X, Teuchos::RCP<MV> MX,
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

    // Perform the Gram-Schmidt transformation for a block of vectors
    std::vector<int> index(1);
    Teuchos::RCP<const MV> tempQ;

    Teuchos::Array<Teuchos::RCP<MV> > MQ(nq);
    // Define the product Q^T * (Op*X)
    for (int i=0; i<nq; i++) {

      // Perform MGS one vector at a time
      for (int ii=0; ii<qcs[i]; ii++) {

        index[0] = ii;
        tempQ = MVT::CloneView( *Q[i], index );
        Teuchos::SerialDenseMatrix<int,ScalarType> tempC( Teuchos::View, *C[i], 1, 1, ii, 0 );

        // Multiply Q' with MX
        {
#ifdef BELOS_TEUCHOS_TIME_MONITOR
          Teuchos::TimeMonitor innerProdTimer( *timerInnerProd_ );
#endif
          MatOrthoManager<ScalarType,MV,OP>::innerProd(*tempQ,X,MX,tempC);
        }
        // Multiply by Q and subtract the result in X
        {
#ifdef BELOS_TEUCHOS_TIME_MONITOR
          Teuchos::TimeMonitor updateTimer( *timerUpdate_ );
#endif
          MVT::MvTimesMatAddMv( -ONE, *tempQ, tempC, ONE, X );
        }
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
          MVT::MvTimesMatAddMv( -ONE, *MQ[i], *C[i], ONE, *MX );
        }
      }
    }

    // Do as many steps of modified Gram-Schmidt as required by max_ortho_steps_
    for (int j = 1; j < max_ortho_steps_; ++j) {

      for (int i=0; i<nq; i++) {

        Teuchos::SerialDenseMatrix<int,ScalarType> C2(qcs[i],1);

        // Perform MGS one vector at a time
        for (int ii=0; ii<qcs[i]; ii++) {

          index[0] = ii;
          tempQ = MVT::CloneView( *Q[i], index );
          Teuchos::SerialDenseMatrix<int,ScalarType> tempC( Teuchos::View, *C[i], 1, 1, ii, 0 );
          Teuchos::SerialDenseMatrix<int,ScalarType> tempC2( Teuchos::View, C2, 1, 1, ii, 0 );

          // Apply another step of modified Gram-Schmidt
          {
#ifdef BELOS_TEUCHOS_TIME_MONITOR
            Teuchos::TimeMonitor innerProdTimer( *timerInnerProd_ );
#endif
            MatOrthoManager<ScalarType,MV,OP>::innerProd( *tempQ, X, MX, tempC2 );
          }
          tempC += tempC2;
          {
#ifdef BELOS_TEUCHOS_TIME_MONITOR
            Teuchos::TimeMonitor updateTimer( *timerUpdate_ );
#endif
            MVT::MvTimesMatAddMv( -ONE, *tempQ, tempC2, ONE, X );
          }

        }

        // Update MX, with the least number of applications of Op as possible
        if (this->_hasOp) {
          if (MQ[i].get()) {
            // MQ was allocated and computed above; use it
            MVT::MvTimesMatAddMv( -ONE, *MQ[i], C2, ONE, *MX );
          }
          else if (xc <= qcs[i]) {
            // MQ was not allocated and computed above; it was cheaper to use X before and it still is
            OPT::Apply( *(this->_Op), X, *MX);
          }
        }
      } // for (int i=0; i<nq; i++)
    } // for (int j = 0; j < max_ortho_steps; ++j)

    return false;
  }

  //////////////////////////////////////////////////////////////////////////////////////////////////
  // Routine to compute the block orthogonalization
  template<class StorageType, class MV, class OP>
  bool
  IMGSOrthoManager<Sacado::MP::Vector<StorageType>, MV, OP>::blkOrtho ( MV &X, Teuchos::RCP<MV> MX,
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

    std::vector<int> index(1);
    Teuchos::Array<Teuchos::RCP<MV> > MQ(nq);
    Teuchos::RCP<const MV> tempQ;

    // Define the product Q^T * (Op*X)
    for (int i=0; i<nq; i++) {

      // Perform MGS one vector at a time
      for (int ii=0; ii<qcs[i]; ii++) {

        index[0] = ii;
        tempQ = MVT::CloneView( *Q[i], index );
        Teuchos::SerialDenseMatrix<int,ScalarType> tempC( Teuchos::View, *C[i], 1, xc, ii, 0 );

        // Multiply Q' with MX
        {
#ifdef BELOS_TEUCHOS_TIME_MONITOR
          Teuchos::TimeMonitor innerProdTimer( *timerInnerProd_ );
#endif
          MatOrthoManager<ScalarType,MV,OP>::innerProd( *tempQ, X, MX, tempC);
        }
        // Multiply by Q and subtract the result in X
        {
#ifdef BELOS_TEUCHOS_TIME_MONITOR
          Teuchos::TimeMonitor updateTimer( *timerUpdate_ );
#endif
          MVT::MvTimesMatAddMv( -ONE, *tempQ, tempC, ONE, X );
        }
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
          MVT::MvTimesMatAddMv( -ONE, *MQ[i], *C[i], ONE, *MX );
        }
      }
    }

    // Do as many steps of modified Gram-Schmidt as required by max_ortho_steps_
    for (int j = 1; j < max_ortho_steps_; ++j) {

      for (int i=0; i<nq; i++) {
        Teuchos::SerialDenseMatrix<int,ScalarType> C2(qcs[i],xc);

        // Perform MGS one vector at a time
        for (int ii=0; ii<qcs[i]; ii++) {

          index[0] = ii;
          tempQ = MVT::CloneView( *Q[i], index );
          Teuchos::SerialDenseMatrix<int,ScalarType> tempC( Teuchos::View, *C[i], 1, xc, ii, 0 );
          Teuchos::SerialDenseMatrix<int,ScalarType> tempC2( Teuchos::View, C2, 1, xc, ii, 0 );

          // Apply another step of modified Gram-Schmidt
          {
#ifdef BELOS_TEUCHOS_TIME_MONITOR
            Teuchos::TimeMonitor innerProdTimer( *timerInnerProd_ );
#endif
            MatOrthoManager<ScalarType,MV,OP>::innerProd( *tempQ, X, MX, tempC2 );
          }
          tempC += tempC2;
          {
#ifdef BELOS_TEUCHOS_TIME_MONITOR
            Teuchos::TimeMonitor updateTimer( *timerUpdate_ );
#endif
            MVT::MvTimesMatAddMv( -ONE, *tempQ, tempC2, ONE, X );
          }
        }

        // Update MX, with the least number of applications of Op as possible
        if (this->_hasOp) {
          if (MQ[i].get()) {
            // MQ was allocated and computed above; use it
            MVT::MvTimesMatAddMv( -ONE, *MQ[i], C2, ONE, *MX );
          }
          else if (xc <= qcs[i]) {
            // MQ was not allocated and computed above; it was cheaper to use X before and it still is
            OPT::Apply( *(this->_Op), X, *MX);
          }
        }
      } // for (int i=0; i<nq; i++)
    } // for (int j = 0; j < max_ortho_steps; ++j)

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

  template<class StorageType, class MV, class OP>
  int
  IMGSOrthoManager<Sacado::MP::Vector<StorageType>, MV, OP>::blkOrthoSing ( MV &X, Teuchos::RCP<MV> MX,
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
      Teuchos::RCP<const MV> tempQ;

      // Define the product Q^T * (Op*X)
      for (int i=0; i<Q.size(); i++) {

        // Perform MGS one vector at a time
        for (int ii=0; ii<qcs[i]; ii++) {

          indX[0] = ii;
          tempQ = MVT::CloneView( *Q[i], indX );
          // Get a view of the current serial dense matrix
          Teuchos::SerialDenseMatrix<int,ScalarType> tempC( Teuchos::View, *C[i], 1, 1, ii, j );

          // Multiply Q' with MX
          MatOrthoManager<ScalarType,MV,OP>::innerProd(*tempQ,*Xj,MXj,tempC);

          // Multiply by Q and subtract the result in Xj
          MVT::MvTimesMatAddMv( -ONE, *tempQ, tempC, ONE, *Xj );
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
            Teuchos::SerialDenseMatrix<int,ScalarType> tempC( Teuchos::View, *C[i], qcs[i], 1, 0, j );
            MVT::MvTimesMatAddMv( -ONE, *MQ[i], tempC, ONE, *MXj );
          }
        }
      }

      // Do any additional steps of modified Gram-Schmidt orthogonalization
      for (int num_ortho_steps=1; num_ortho_steps < max_ortho_steps_; ++num_ortho_steps) {

        for (int i=0; i<Q.size(); i++) {
          Teuchos::SerialDenseMatrix<int,ScalarType> C2( qcs[i], 1 );

          // Perform MGS one vector at a time
          for (int ii=0; ii<qcs[i]; ii++) {

            indX[0] = ii;
            tempQ = MVT::CloneView( *Q[i], indX );
            // Get a view of the current serial dense matrix
            Teuchos::SerialDenseMatrix<int,ScalarType> tempC2( Teuchos::View, C2, 1, 1, ii );

            // Apply another step of modified Gram-Schmidt
            MatOrthoManager<ScalarType,MV,OP>::innerProd( *tempQ, *Xj, MXj, tempC2);
            MVT::MvTimesMatAddMv( -ONE, *tempQ, tempC2, ONE, *Xj );
          }

          // Add the coefficients into C[i]
          Teuchos::SerialDenseMatrix<int,ScalarType> tempC( Teuchos::View, *C[i], qcs[i], 1, 0, j );
          tempC += C2;

          // Update MXj, with the least number of applications of Op as possible
          if (this->_hasOp) {
            if (MQ[i].get()) {
              // MQ was allocated and computed above; use it
              MVT::MvTimesMatAddMv( -ONE, *MQ[i], C2, ONE, *MXj );
            }
            else if (xc <= qcs[i]) {
              // MQ was not allocated and computed above; it was cheaper to use X before and it still is
              OPT::Apply( *(this->_Op), *Xj, *MXj);
            }
          }
        } // for (int i=0; i<Q.size(); i++)

      } // for (int num_ortho_steps=1; num_ortho_steps < max_ortho_steps_; ++num_ortho_steps)

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
        for (int num_orth=0; num_orth<max_ortho_steps_; num_orth++) {

          for (int i=0; i<Q.size(); i++) {
            Teuchos::SerialDenseMatrix<int,ScalarType> product( qcs[i], 1 );

            // Perform MGS one vector at a time
            for (int ii=0; ii<qcs[i]; ii++) {

              indX[0] = ii;
              tempQ = MVT::CloneView( *Q[i], indX );
              Teuchos::SerialDenseMatrix<int,ScalarType> tempC( Teuchos::View, product, 1, 1, ii );

              // Apply another step of modified Gram-Schmidt
              MatOrthoManager<ScalarType,MV,OP>::innerProd( *tempQ, *tempXj, tempMXj, tempC );
              MVT::MvTimesMatAddMv( -ONE, *tempQ, tempC, ONE, *tempXj );

            }

            // Update MXj, with the least number of applications of Op as possible
            if (this->_hasOp) {
              if (MQ[i].get()) {
                // MQ was allocated and computed above; use it
                MVT::MvTimesMatAddMv( -ONE, *MQ[i], product, ONE, *tempMXj );
              }
              else if (xc <= qcs[i]) {
                // MQ was not allocated and computed above; it was cheaper to use X before and it still is
                OPT::Apply( *(this->_Op), *tempXj, *tempMXj);
              }
            }
          } // for (int i=0; i<nq; i++)
        } // for (int num_orth=0; num_orth<max_orth_steps_; num_orth++)

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
  
} // namespace Belos

#endif // BELOS_IMGS_ORTHOMANAGER_HPP

