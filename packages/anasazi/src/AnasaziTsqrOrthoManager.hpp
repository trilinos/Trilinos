// @HEADER
// *****************************************************************************
//                 Anasazi: Block Eigensolvers Package
//
// Copyright 2004 NTESS and the Anasazi contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/// \file AnasaziTsqrOrthoManager.hpp
/// \brief Orthogonalization manager based on Tall Skinny QR (TSQR)

#ifndef __AnasaziTsqrOrthoManager_hpp
#define __AnasaziTsqrOrthoManager_hpp

#include "AnasaziTsqrOrthoManagerImpl.hpp"
#include "AnasaziSVQBOrthoManager.hpp"


namespace Anasazi {

  /// \class OutOfPlaceNormalizerMixin
  /// \brief Mixin for out-of-place orthogonalization 
  /// \author Mark Hoemmen
  ///
  /// This class presents an abstract interface for multiple
  /// inheritance ("mixin") for special orthogonalization methods that
  /// normalize "out-of-place."  OrthoManager and MatOrthoManager both
  /// normalize (and projectAndNormalize) multivectors "in place,"
  /// meaning that the input and output multivectors are the same (X,
  /// in both cases).  Gram-Schmidt (modified or classical) is an
  /// example of an orthogonalization method that can normalize (and
  /// projectAndNormalize) in place.  TSQR (the Tall Skinny QR
  /// factorization, see \c TsqrOrthoManager.hpp for references) is an
  /// orthogonalization method which cannot normalize (or
  /// projectAndNormalize) in place.
  ///
  /// Tsqr(Mat)OrthoManager implements (Mat)OrthoManager's normalize()
  /// and projectAndNormalize() methods with scratch space and
  /// copying.  However, if you handle Tsqr(Mat)OrthoManager through
  /// this mixin, you can exploit TSQR's unique interface to avoid
  /// copying back and forth between scratch space.
  template<class Scalar, class MV>
  class OutOfPlaceNormalizerMixin {
  public:
    typedef Scalar scalar_type;
    typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType magnitude_type;
    /// \typedef multivector_type
    /// \brief Multivector type with which this class was specialized.
    typedef MV multivector_type;

    typedef Teuchos::SerialDenseMatrix<int, Scalar> mat_type;
    typedef Teuchos::RCP<mat_type>                  mat_ptr;

    /// \brief Normalize X into Q*B.
    ///
    /// \param X [in/out] On input: Multivector to normalize.  On
    ///   output: Possibly overwritten with invalid values.
    /// \param Q [out] On output: Normalized multivector.
    /// \param B [out] On output: Normalization coefficients.
    ///
    /// \return Rank of the input multivector X.
    virtual int 
    normalizeOutOfPlace (MV& X, MV& Q, mat_ptr B) const = 0;

    /// \brief Project and normalize X_in into X_out.
    ///
    /// Project X_in against Q, storing projection coefficients in C,
    /// and normalize X_in into X_out, storing normalization
    /// coefficients in B.  On output, X_out has the resulting
    /// orthogonal vectors.  X_in may be overwritten with invalid
    /// values.
    ///
    /// \param X_in [in/out] On input: The vectors to project against
    ///   Q and normalize.  On output: possibly overwritten with 
    ///   invalid values.
    /// \param X_out [out] The normalized input vectors after
    ///   projection against Q.
    /// \param C [out] Projection coefficients 
    /// \param B [out] Normalization coefficients
    /// \param Q [in] The orthogonal basis against which to project
    ///
    /// \return Rank of X_in after projection
    virtual int 
    projectAndNormalizeOutOfPlace (MV& X_in, 
				   MV& X_out, 
				   Teuchos::Array<mat_ptr> C,
				   mat_ptr B,
				   Teuchos::ArrayView<Teuchos::RCP<const MV> > Q) const = 0;

    //! Trivial virtual destructor, to silence compiler warnings.
    virtual ~OutOfPlaceNormalizerMixin () {}
  };

  /// \class TsqrOrthoManager
  /// \brief TSQR-based OrthoManager subclass
  /// \author Mark Hoemmen
  ///
  /// Subclass of OrthoManager, implemented using TsqrOrthoManagerImpl
  /// (TSQR + Block Gram-Schmidt).
  template<class Scalar, class MV>
  class TsqrOrthoManager : 
    public OrthoManager<Scalar, MV>, 
    public OutOfPlaceNormalizerMixin<Scalar, MV>,
    public Teuchos::ParameterListAcceptor
  {
  public:
    typedef Scalar scalar_type;
    typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType magnitude_type;
    //! \typedef Multivector type with which this class was specialized
    typedef MV multivector_type;

    typedef Teuchos::SerialDenseMatrix<int, Scalar> mat_type;
    typedef Teuchos::RCP<mat_type>                  mat_ptr;

    void setParameterList (const Teuchos::RCP<Teuchos::ParameterList>& params) {
      impl_.setParameterList (params);
    }

    Teuchos::RCP<Teuchos::ParameterList> getNonconstParameterList () {
      return impl_.getNonconstParameterList ();
    }

    Teuchos::RCP<Teuchos::ParameterList> unsetParameterList () {
      return impl_.unsetParameterList ();
    }

    /// \brief Default valid parameter list.
    ///
    /// Get a (pointer to a) default list of parameters for
    /// configuring a TsqrOrthoManager instance.
    ///
    /// \note TSQR implementation configuration options are stored
    ///   under "TsqrImpl" as a sublist.
    Teuchos::RCP<const Teuchos::ParameterList> getValidParameters() const {
      return impl_.getValidParameters();
    }

    /// \brief Get "fast" parameters for TsqrOrthoManager.
    ///
    /// Get a (pointer to a) list of parameters for configuring a
    /// TsqrOrthoManager instance for maximum speed, at the cost of
    /// accuracy (no block reorthogonalization) and robustness to rank
    /// deficiency (no randomization of the null space basis).
    ///
    /// \note TSQR implementation configuration options are stored
    ///   under "TsqrImpl" as a sublist.
    Teuchos::RCP<const Teuchos::ParameterList> getFastParameters() const {
      return impl_.getFastParameters();
    }

    /// \brief Constructor (that sets user-specified parameters).
    ///
    /// \param params [in/out] Configuration parameters, both for this
    ///   orthogonalization manager, and for TSQR itself (as the
    ///   "TsqrImpl" sublist).  This can be null, in which case
    ///   default parameters will be set for now; you can always call
    ///   \c setParameterList() later to change these.
    /// 
    /// \param label [in] Label for timers.  This only matters if the
    ///   compile-time option for enabling timers is set.
    ///
    /// Call \c getValidParameters() for default parameters and their
    /// documentation, including TSQR implementation parameters.  Call
    /// \c getFastParameters() to get documented parameters for faster
    /// computation, possibly at the expense of accuracy and
    /// robustness.
    TsqrOrthoManager (const Teuchos::RCP<Teuchos::ParameterList>& params, 
		      const std::string& label = "Anasazi") :
      impl_ (params, label)
    {}

    /// \brief Constructor (that sets default parameters).
    ///
    /// \param label [in] Label for timers.  This only matters if the
    ///   compile-time option for enabling timers is set.
    TsqrOrthoManager (const std::string& label) :
      impl_ (label)
    {}

    //! Destructor, declared virtual for safe inheritance.
    virtual ~TsqrOrthoManager() {}

    void innerProd (const MV &X, const MV& Y, mat_type& Z) const {
      return impl_.innerProd (X, Y, Z);
    }

    void norm (const MV& X, std::vector<magnitude_type>& normVec) const {
      return impl_.norm (X, normVec);
    }

    void 
    project (MV &X, 
	     Teuchos::Array<Teuchos::RCP<const MV> > Q,
	     Teuchos::Array<Teuchos::RCP<Teuchos::SerialDenseMatrix<int,Scalar> > > C 
	     = Teuchos::tuple(Teuchos::RCP< Teuchos::SerialDenseMatrix<int,Scalar> >(Teuchos::null))) const
    {
      return impl_.project (X, C, Q);
    }

    int 
    normalize (MV &X, mat_ptr B = Teuchos::null) const
    {
      return impl_.normalize (X, B);
    }

    int 
    projectAndNormalize (MV &X, 
			 Teuchos::Array<Teuchos::RCP<const MV> > Q,
			 Teuchos::Array<Teuchos::RCP<Teuchos::SerialDenseMatrix<int,Scalar> > > C 
			 = Teuchos::tuple(Teuchos::RCP< Teuchos::SerialDenseMatrix<int,Scalar> >(Teuchos::null)),
			 Teuchos::RCP<Teuchos::SerialDenseMatrix<int,Scalar> > B = Teuchos::null) const
    {
      return impl_.projectAndNormalize (X, C, B, Q);
    }

    /// \brief Normalize X into Q*B, overwriting X with invalid values.
    ///
    /// We expose this interface to applications because TSQR is not
    /// able to compute an orthogonal basis in place; it needs scratch
    /// space.  Applications can exploit this interface to avoid
    /// excessive copying of vectors when using TSQR for
    /// orthogonalization.
    ///
    /// \param X [in/out] Input vector(s) to normalize
    /// \param Q [out] Normalized output vector(s)
    /// \param B [out] Normalization coefficients
    ///
    /// \return Rank of X.
    ///
    /// \note Q must have at least as many columns as X.  It may have
    ///   more columns than X; those columns are ignored.
    int 
    normalizeOutOfPlace (MV& X, MV& Q, mat_ptr B) const
    {
      return impl_.normalizeOutOfPlace (X, Q, B);
    }

    /// \brief Project and normalize X_in into X_out; overwrite X_in.
    ///
    /// Project X_in against Q, storing projection coefficients in C,
    /// and normalize X_in into X_out, storing normalization
    /// coefficients in B.  On output, X_out has the resulting
    /// orthogonal vectors and X_in is overwritten with invalid
    /// values.
    ///
    /// \param X_in [in/out] On input: The vectors to project against
    ///   Q and normalize.  Overwritten with invalid values on output.
    /// \param X_out [out] The normalized input vectors after
    ///   projection against Q.
    /// \param C [out] Projection coefficients 
    /// \param B [out] Normalization coefficients
    /// \param Q [in] The orthogonal basis against which to project
    ///
    /// \return Rank of X_in after projection
    ///
    /// \note We expose this interface to applications for the same
    ///   reason that we expose \c normalizeOutOfPlace().
    int 
    projectAndNormalizeOutOfPlace (MV& X_in, 
				   MV& X_out,
				   Teuchos::Array<mat_ptr> C,
				   mat_ptr B,
				   Teuchos::ArrayView<Teuchos::RCP<const MV> > Q) const
    {
      return impl_.projectAndNormalizeOutOfPlace (X_in, X_out, C, B, Q);
    }

    magnitude_type orthonormError (const MV &X) const {
      return impl_.orthonormError (X);
    }

    magnitude_type orthogError (const MV &X1, const MV &X2) const {
      return impl_.orthogError (X1, X2);
    }

  private:
    /// \brief The implementation of TSQR.
    ///
    /// The object is delcared "mutable" because it has internal
    /// scratch space state, and because the OrthoManager interface
    /// requires most of its methods to be declared const.
    mutable TsqrOrthoManagerImpl<Scalar, MV> impl_;
  };


  /// \class TsqrMatOrthoManager
  /// \brief MatOrthoManager subclass using TSQR or SVQB
  ///
  /// When the inner product matrix has not been set, this class uses
  /// TSQR + Block Gram-Schmidt (via \c TsqrOrthoManagerImpl).  If the
  /// inner product matrix <i>has</i> been set, then this class uses
  /// the SVQB algorithm (Stathopoulos and Wu 2002: CholeskyQR + SVD)
  /// for orthogonalization.
  /// 
  /// TSQR uses multivector scratch space.  However, scratch space
  /// initialization is "lazy," so scratch space will not be allocated
  /// if TSQR is not used.
  template<class Scalar, class MV, class OP>
  class TsqrMatOrthoManager : 
    public MatOrthoManager<Scalar, MV, OP>,
    public OutOfPlaceNormalizerMixin<Scalar, MV>,
    public Teuchos::ParameterListAcceptorDefaultBase
  {
  public:
    typedef Scalar scalar_type;
    typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType magnitude_type;
    //! Multivector type with which this class was specialized
    typedef MV multivector_type;
    //! Operator type with which this class was specialized
    typedef OP operator_type;

    typedef Teuchos::SerialDenseMatrix<int, Scalar> mat_type;
    typedef Teuchos::RCP<mat_type>                  mat_ptr;

  private:
    /// \typedef base_type
    ///
    /// This will be used to help C++ resolve getOp().  We can't call
    /// getOp() directly, because C++ can't figure out that it belongs
    /// to the base class MatOrthoManager.  (Remember that at this
    /// point, we might not have specialized the specific base class
    /// yet; it's just a template at the moment and not a "real
    /// class.")
    typedef MatOrthoManager<Scalar, MV, OP> base_type;

    /// \typedef tsqr_type
    /// \brief Implementation of TSQR-based orthogonalization
    typedef TsqrOrthoManagerImpl<Scalar, MV> tsqr_type;

    /// \typedef svqb_type
    /// \brief Implementation of SVQB-based orthogonalization
    typedef SVQBOrthoManager<Scalar, MV, OP> svqb_type;

    /// \typedef MVT
    /// \brief Traits class for the multivector type
    typedef MultiVecTraits<Scalar, MV> MVT;

  public:
    /// \brief Constructor (that sets user-specified parameters).
    ///
    /// \param params [in/out] Configuration parameters, both for this
    ///   orthogonalization manager, and for TSQR itself (as the "TSQR
    ///   implementation" sublist).  This can be null, in which case
    ///   default parameters will be set for now; you can always call
    ///   \c setParameterList() later to change these.
    /// 
    /// \param label [in] Label for timers.  This only matters if the
    ///   compile-time option for enabling timers is set.
    ///
    /// \param Op [in] Inner product with respect to which to
    ///   orthogonalize vectors.  If Teuchos::null, use the Euclidean
    ///   inner product.
    ///
    /// Call \c getValidParameters() for default parameters and their
    /// documentation, including TSQR implementation parameters.  Call
    /// \c getFastParameters() to get documented parameters for faster
    /// computation, possibly at the expense of accuracy and
    /// robustness.
    TsqrMatOrthoManager (const Teuchos::RCP<Teuchos::ParameterList>& params,
			 const std::string& label = "Belos",
			 Teuchos::RCP<const OP> Op = Teuchos::null) :
      MatOrthoManager<Scalar, MV, OP> (Op),
      tsqr_ (params, label),
      pSvqb_ (Teuchos::null) // Initialized lazily
    {}

    /// \brief Constructor (that sets default parameters).
    ///
    /// \param Op [in] Inner product with respect to which to
    ///   orthogonalize vectors.  If Teuchos::null, use the Euclidean
    ///   inner product.
    ///
    /// \param label [in] Label for timers.  This only matters if the
    ///   compile-time option for enabling timers is set.
    TsqrMatOrthoManager (const std::string& label = "Belos",
			 Teuchos::RCP<const OP> Op = Teuchos::null) :
      MatOrthoManager<Scalar, MV, OP> (Op),
      tsqr_ (label),
      pSvqb_ (Teuchos::null) // Initialized lazily
    {}

    //! Destructor (declared virtual for memory safety of derived classes).
    virtual ~TsqrMatOrthoManager() {}

    /// \brief Get default parameters for TsqrMatOrthoManager.
    ///
    /// Get a (pointer to a) default list of parameters for
    /// configuring a TsqrMatOrthoManager instance.
    ///
    /// \note TSQR implementation configuration options are stored
    ///   under "TSQR implementation" as a sublist.
    Teuchos::RCP<const Teuchos::ParameterList> getValidParameters() const {
      return tsqr_.getValidParameters ();
    }

    /// \brief Get "fast" parameters for TsqrMatOrthoManager.
    ///
    /// Get a (pointer to a) list of parameters for configuring a
    /// TsqrMatOrthoManager instance for maximum speed, at the cost of
    /// accuracy (no block reorthogonalization) and robustness to rank
    /// deficiency (no randomization of the null space basis).
    ///
    /// \note TSQR implementation configuration options are stored
    ///   under "TSQR implementation" as a sublist.
    Teuchos::RCP<const Teuchos::ParameterList> getFastParameters() {
      return tsqr_.getFastParameters ();
    }

    void setParameterList (const Teuchos::RCP<Teuchos::ParameterList>& params) {
      tsqr_.setParameterList (params);
    }

    virtual void 
    setOp (Teuchos::RCP< const OP > Op) 
    {
      // We override the base class' setOp() so that the
      // SVQBOrthoManager instance gets the new op.
      //
      // The base_type notation helps C++ resolve the name for a
      // member function of a templated base class.
      base_type::setOp (Op); // base class gets a copy of the Op too

      if (! Op.is_null()) {
	ensureSvqbInit (); // Make sure the SVQB object has been initialized
	pSvqb_->setOp (Op);
      }
    }

    Teuchos::RCP<const OP> getOp () const { 
      // The base_type notation helps C++ resolve the name for a
      // member function of a templated base class.
      return base_type::getOp(); 
    }

    void 
    projectMat (MV &X, 
		Teuchos::Array<Teuchos::RCP<const MV> > Q,
		Teuchos::Array<Teuchos::RCP<mat_type> > C = 
		  Teuchos::tuple (Teuchos::RCP<mat_type> (Teuchos::null)),
		Teuchos::RCP<MV> MX = Teuchos::null,
		Teuchos::Array<Teuchos::RCP<const MV> > MQ = 
		  Teuchos::tuple (Teuchos::null)) const
    {
      if (getOp().is_null()) {
	// FIXME (mfh 12 Jan 2011): Do we need to check if C is null
	// and allocate space if so?
	tsqr_.project (X, C, Q);
	// FIXME (mfh 20 Jul 2010) What about MX and MQ?
	// 
	// FIXME (mfh 12 Jan 2011) Need to port MVT::Assign from Belos to Anasazi
#if 0
	if (! MX.is_null()) {
	  // MX gets a copy of X; M is the identity operator.
	  MVT::Assign (X, *MX);
	}
#endif // 0
      } else {
	ensureSvqbInit ();
	pSvqb_->projectMat (X, Q, C, MX, MQ);
      }
    }

    int 
    normalizeMat (MV &X, 
		  mat_ptr B = Teuchos::null,
		  Teuchos::RCP<MV> MX = Teuchos::null) const
    {
      if (getOp().is_null()) {
	// FIXME (mfh 12 Jan 2011): Do we need to check if B is null
	// and allocate space if so?
	const int rank = tsqr_.normalize (X, B);
	// FIXME (mfh 20 Jul 2010) What about MX?
	// 
	// FIXME (mfh 12 Jan 2011) Need to port MVT::Assign from Belos to Anasazi
#if 0
	if (! MX.is_null()) {
	  // MX gets a copy of X; M is the identity operator.
	  MVT::Assign (X, *MX);
	}
#endif // 0
	return rank;
      } else {
	ensureSvqbInit ();
	return pSvqb_->normalizeMat (X, B, MX);
      }
    }

    int 
    projectAndNormalizeMat (MV &X, 
			    Teuchos::Array<Teuchos::RCP<const MV> > Q,
			    Teuchos::Array<Teuchos::RCP<mat_type> > C = 
			      Teuchos::tuple (Teuchos::RCP<mat_type> (Teuchos::null)),
			    Teuchos::RCP<mat_type> B = Teuchos::null,
			    Teuchos::RCP<MV> MX = Teuchos::null,
			    Teuchos::Array<Teuchos::RCP<const MV> > MQ = 
			      Teuchos::tuple (Teuchos::RCP<const MV> (Teuchos::null))) const 
    {
      if (getOp().is_null()) {
	// FIXME (mfh 12 Jan 2011): Do we need to check if C and B
	// are null and allocate space if so?
	const int rank = tsqr_.projectAndNormalize (X, C, B, Q); 
	// FIXME (mfh 20 Jul 2010) What about MX and MQ?
	// mfh 12 Jan 2011: Ignore MQ (assume Q == MQ).
	// 
	// FIXME (mfh 12 Jan 2011) Need to port MVT::Assign from Belos to Anasazi
#if 0
	if (! MX.is_null()) {
	  // MX gets a copy of X; M is the identity operator.
	  MVT::Assign (X, *MX);
	}
#endif // 0
	return rank;
      } else {
	ensureSvqbInit ();
	return pSvqb_->projectAndNormalizeMat (X, Q, C, B, MX, MQ);
      }
    }

    int 
    normalizeOutOfPlace (MV& X, MV& Q, mat_ptr B) const
    {
      if (getOp().is_null()) {
	return tsqr_.normalizeOutOfPlace (X, Q, B);
      } else {
	// SVQB normalizes in place, so we have to copy.
	ensureSvqbInit ();
	const int rank = pSvqb_->normalize (X, B);
	MVT::Assign (X, Q);
	return rank;
      }
    }

    int 
    projectAndNormalizeOutOfPlace (MV& X_in, 
				   MV& X_out,
				   Teuchos::Array<mat_ptr> C,
				   mat_ptr B,
				   Teuchos::ArrayView<Teuchos::RCP<const MV> > Q) const
    {
      using Teuchos::null;

      if (getOp().is_null()) {
	return tsqr_.projectAndNormalizeOutOfPlace (X_in, X_out, C, B, Q);
      } else {
	ensureSvqbInit ();
	// SVQB's projectAndNormalize wants an Array, not an ArrayView.
	// Copy into a temporary Array and copy back afterwards.
	Teuchos::Array<Teuchos::RCP<const MV> > Q_array (Q);
	const int rank = pSvqb_->projectAndNormalize (X_in, Q_array, C, B);
	Q.assign (Q_array);
	// SVQB normalizes in place, so we have to copy X_in to X_out.
	MVT::Assign (X_in, X_out);
	return rank;
      }
    }

    magnitude_type
    orthonormErrorMat (const MV &X, Teuchos::RCP<const MV> MX = Teuchos::null) const
    {
      if (getOp().is_null()) {
	return tsqr_.orthonormError (X);
	// FIXME (mfh 20 Jul 2010) What about MX?
      } else {
	ensureSvqbInit ();
	return pSvqb_->orthonormErrorMat (X, MX);
      }
    }

    magnitude_type
    orthogErrorMat (const MV &X, 
		    const MV &Y,
		    Teuchos::RCP<const MV> MX = Teuchos::null, 
		    Teuchos::RCP<const MV> MY = Teuchos::null) const
    {
      if (getOp().is_null()) {
	return tsqr_.orthogError (X, Y);
	// FIXME (mfh 20 Jul 2010) What about MX and MY?
      } else {
	ensureSvqbInit ();
	return pSvqb_->orthogErrorMat (X, Y, MX, MY);
      }
    }

  private:
    //! Ensure that the SVQBOrthoManager instance is initialized.
    void 
    ensureSvqbInit () const
    {
      // NOTE (mfh 12 Oct 2011) For now, we rely on the default values
      // of SVQB parameters.
      if (pSvqb_.is_null()) {
	pSvqb_ = Teuchos::rcp (new svqb_type (getOp()));
      }
    }

    /// \brief TSQR + BGS orthogonalization manager implementation.
    ///
    /// We use this when getOp() == null (Euclidean inner product).
    /// It gets initialized in the constructor.  This must be declared
    /// mutable because of the requirements of the MatOrthoManager
    /// interface.  tsqr_type's interface is honest about what methods
    /// are really const.
    mutable tsqr_type tsqr_;

    /// \brief SVQB orthogonalization manager.
    ///
    /// We initialize and use this when getOp() != null (could be a
    /// non-Euclidean inner product, but not necessarily).
    mutable Teuchos::RCP<svqb_type> pSvqb_;
  };

} // namespace Anasazi

#endif // __AnasaziTsqrOrthoManager_hpp

