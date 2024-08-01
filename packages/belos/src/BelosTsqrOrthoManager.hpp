// @HEADER
// *****************************************************************************
//                 Belos: Block Linear Solvers Package
//
// Copyright 2004-2016 NTESS and the Belos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/// \file BelosTsqrOrthoManager.hpp
/// \brief Orthogonalization manager based on Tall Skinny QR (TSQR)

#ifndef __BelosTsqrOrthoManager_hpp
#define __BelosTsqrOrthoManager_hpp

#include "BelosTsqrOrthoManagerImpl.hpp"
// Belos doesn't have an SVQB orthomanager implemented yet, so we just
// use a default orthomanager for the case where the inner product
// matrix is nontrivial.
#include "BelosDGKSOrthoManager.hpp"


namespace Belos {

/// \class OutOfPlaceNormalizerMixin
/// \brief Mixin for out-of-place orthogonalization
/// \author Mark Hoemmen
///
/// This class presents an abstract interface for multiple inheritance
/// ("mixin") for special orthogonalization methods that normalize
/// "out-of-place."  OrthoManager and MatOrthoManager both normalize
/// (and projectAndNormalize) multivectors "in place," meaning that
/// the input and output multivectors are the same (X, in both cases).
/// Gram-Schmidt (modified or classical) is an example of an
/// orthogonalization method that can normalize (and
/// projectAndNormalize) in place.  TSQR (the Tall Skinny QR
/// factorization, see TsqrOrthoManager.hpp for references) is an
/// orthogonalization method which cannot normalize (or
/// projectAndNormalize) in place.
///
/// Tsqr(Mat)OrthoManager implements (Mat)OrthoManager's normalize()
/// and projectAndNormalize() methods with scratch space and copying.
/// However, if you handle Tsqr(Mat)OrthoManager through this mixin,
/// you can exploit TSQR's unique interface to avoid copying back and
/// forth between scratch space.
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
    public OutOfPlaceNormalizerMixin<Scalar, MV>
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
  ///   setParameterList() later to change these.
  ///
  /// \param label [in] Label for timers.  This only matters if the
  ///   compile-time option for enabling timers is set.
  ///
  /// Call getValidParameters() for default parameters and their
  /// documentation, including TSQR implementation parameters.  Call
  /// getFastParameters() to get documented parameters for faster
  /// computation, possibly at the expense of accuracy and robustness.
  TsqrOrthoManager (const Teuchos::RCP<Teuchos::ParameterList>& params,
                    const std::string& label = "Belos") :
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

  /// \brief Set callback to be invoked on reorthogonalization.
  ///
  /// This callback is invoked right after the first projection
  /// step, and only if reorthogonalization will be necessary.  It
  /// is called before actually reorthogonalizing.  The first
  /// argument gives the norms of the columns of the input
  /// multivector before the first projection pass, and the second
  /// argument gives their norms after the first projection pass.
  ///
  /// The callback is null by default.  If the callback is null, no
  /// callback will be invoked.
  ///
  /// For details and suggested uses, please refer to the
  /// documentation of ReorthogonalizationCallback.
  ///
  /// \warning Please do not rely on the interface to this method.
  ///   This method may change or go away at any time.
  ///
  /// \warning We assume that the input arguments of the callback's
  ///   operator() are only valid views within the scope of the
  ///   function.  Your callback should not keep the views.
  void
  setReorthogonalizationCallback (const Teuchos::RCP<ReorthogonalizationCallback<Scalar> >& callback)
  {
    impl_.setReorthogonalizationCallback (callback);
  }

  void innerProd (const MV &X, const MV &Y, mat_type& Z) const {
    return impl_.innerProd (X, Y, Z);
  }

  void norm (const MV& X, std::vector<magnitude_type>& normVec) const {
    return impl_.norm (X, normVec);
  }

  void
  project (MV &X,
           Teuchos::Array<mat_ptr> C,
           Teuchos::ArrayView<Teuchos::RCP<const MV> > Q) const
  {
    return impl_.project (X, C, Q);
  }

  int
  normalize (MV &X, mat_ptr B) const
  {
    return impl_.normalize (X, B);
  }

protected:
  virtual int
  projectAndNormalizeImpl (MV &X,
                           Teuchos::Array<mat_ptr> C,
                           mat_ptr B,
                           Teuchos::ArrayView<Teuchos::RCP<const MV> > Q) const
  {
    return impl_.projectAndNormalize (X, C, B, Q);
  }

public:
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
  ///   reason that we expose normalizeOutOfPlace().
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

  /// Set the label for (the timers for) this orthogonalization
  /// manager, and create new timers if the label has changed.
  ///
  /// \param label [in] New label for timers
  ///
  /// \note Belos::OrthoManager wants this virtual function to be
  ///   implemented; Anasazi::OrthoManager does not.
  void setLabel (const std::string& label) {
    impl_.setLabel (label);
  }

  const std::string& getLabel() const { return impl_.getLabel(); }

private:
  /// \brief The implementation of TSQR.
  ///
  /// The object is delcared "mutable" because it has internal
  /// scratch space state, and because the OrthoManager interface
  /// requires most of its methods to be declared const.
  mutable TsqrOrthoManagerImpl<Scalar, MV> impl_;

  //! Label for timers (if timers are enabled at compile time).
  std::string label_;
};


/// \class TsqrMatOrthoManager
/// \brief MatOrthoManager subclass using TSQR or DGKS
/// \author Mark Hoemmen
///
/// When the inner product matrix has not been set, this class uses
/// TSQR + Block Gram-Schmidt (via TsqrOrthoManagerImpl).  If the
/// inner product matrix <i>has</i> been set, then this class uses
/// classical Gram-Schmidt with reorthogonalization (via
/// DGKSOrthoManager).
///
/// TSQR uses multivector scratch space.  However, scratch space
/// initialization is "lazy," so scratch space will not be allocated
/// if TSQR is not used.
template<class Scalar, class MV, class OP>
class TsqrMatOrthoManager :
    public MatOrthoManager<Scalar, MV, OP>,
    public OutOfPlaceNormalizerMixin<Scalar, MV>
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

  /// \typedef dgks_type
  /// \brief Implementation of DGKS-based orthogonalization
  typedef DGKSOrthoManager<Scalar, MV, OP> dgks_type;

  /// \typedef MVT
  /// \brief Traits class for the multivector type
  typedef MultiVecTraits<Scalar, MV> MVT;

public:
  /// \brief Constructor (that sets user-specified parameters).
  ///
  /// \param label [in] Label for timers.  This only matters if the
  ///   compile-time option for enabling timers is set.
  ///
  /// \param params [in/out] Configuration parameters, both for this
  ///   orthogonalization manager, and for TSQR itself (as the "TSQR
  ///   implementation" sublist).  This can be null, in which case
  ///   default parameters will be set for now; you can always call
  ///   setParameterList() later to change these.
  ///
  /// \param Op [in] Inner product with respect to which to
  ///   orthogonalize vectors.  If Teuchos::null, use the Euclidean
  ///   inner product.
  ///
  /// Call getValidParameters() for default parameters and their
  /// documentation, including TSQR implementation parameters.  Call
  /// getFastParameters() to get documented parameters for faster
  /// computation, possibly at the expense of accuracy and
  /// robustness.
  TsqrMatOrthoManager (const Teuchos::RCP<Teuchos::ParameterList>& params,
                       const std::string& label = "Belos",
                       Teuchos::RCP<const OP> Op = Teuchos::null) :
    MatOrthoManager<Scalar, MV, OP> (Op),
    tsqr_ (params, label),
    pDgks_ (Teuchos::null) // Initialized lazily
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
    pDgks_ (Teuchos::null) // Initialized lazily
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

  const std::string& getLabel() const { return tsqr_.getLabel (); }

  void
  setOp (Teuchos::RCP<const OP> Op)
  {
    // We override the base class' setOp() so that the
    // DGKSOrthoManager instance gets the new op.
    //
    // The base_type notation helps C++ resolve the name for a
    // member function of a templated base class.
    base_type::setOp (Op); // base class gets a copy of the Op too

    if (! Op.is_null()) {
      ensureDgksInit (); // Make sure the DGKS object has been initialized
      pDgks_->setOp (Op);
    }
  }

  Teuchos::RCP<const OP> getOp () const {
    // The base_type notation helps C++ resolve the name for a
    // member function of a templated base class.
    return base_type::getOp();
  }

  void
  project (MV &X,
           Teuchos::RCP<MV> MX,
           Teuchos::Array<mat_ptr> C,
           Teuchos::ArrayView<Teuchos::RCP<const MV> > Q) const
  {
    if (getOp().is_null()) {
      tsqr_.project (X, C, Q);
      if (! MX.is_null()) {
        // MX gets a copy of X; M is the identity operator.
        MVT::Assign (X, *MX);
      }
    } else {
      ensureDgksInit ();
      pDgks_->project (X, MX, C, Q);
    }
  }

  void
  project (MV &X,
           Teuchos::Array<mat_ptr> C,
           Teuchos::ArrayView<Teuchos::RCP<const MV> > Q) const
  {
    project (X, Teuchos::null, C, Q);
  }

  int
  normalize (MV& X, Teuchos::RCP<MV> MX, mat_ptr B) const
  {
    if (getOp().is_null()) {
      const int rank = tsqr_.normalize (X, B);
      if (! MX.is_null()) {
        // MX gets a copy of X; M is the identity operator.
        MVT::Assign (X, *MX);
      }
      return rank;
    } else {
      ensureDgksInit ();
      return pDgks_->normalize (X, MX, B);
    }
  }

  int normalize (MV& X, mat_ptr B) const {
    return normalize (X, Teuchos::null, B);
  }

  // Attempted fix for a warning about hiding
  // OrthoManager::projectAndNormalize(). Unfortunately, this fix turns out
  // to produce a compilation error with cray++, see bug #6129,
  // <https://software.sandia.gov/bugzilla/show_bug.cgi?id=6129>.
  //using Belos::OrthoManager<Scalar, MV>::projectAndNormalize;

protected:
  virtual int
  projectAndNormalizeWithMxImpl (MV &X,
                                 Teuchos::RCP<MV> MX,
                                 Teuchos::Array<mat_ptr> C,
                                 mat_ptr B,
                                 Teuchos::ArrayView<Teuchos::RCP<const MV> > Q) const
  {
    if (getOp().is_null()) {
      const int rank = tsqr_.projectAndNormalize (X, C, B, Q);
      if (! MX.is_null()) {
        // MX gets a copy of X; M is the identity operator.
        MVT::Assign (X, *MX);
      }
      return rank;
    } else {
      ensureDgksInit ();
      return pDgks_->projectAndNormalize (X, MX, C, B, Q);
    }
  }

public:
  int
  normalizeOutOfPlace (MV& X, MV& Q, mat_ptr B) const
  {
    if (getOp().is_null()) {
      return tsqr_.normalizeOutOfPlace (X, Q, B);
    } else {
      // DGKS normalizes in place, so we have to copy.
      ensureDgksInit ();
      const int rank = pDgks_->normalize (X, B);
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
      // DGKS normalizes in place, so we have to copy.
      ensureDgksInit ();
      const int rank = pDgks_->projectAndNormalize (X_in, null, C, B, Q);
      MVT::Assign (X_in, X_out);
      return rank;
    }
  }

  magnitude_type
  orthonormError (const MV &X, Teuchos::RCP<const MV> MX) const
  {
    if (getOp().is_null()) {
      return tsqr_.orthonormError (X); // Ignore MX
    } else {
      ensureDgksInit ();
      return pDgks_->orthonormError (X, MX);
    }
  }

  magnitude_type orthonormError (const MV &X) const {
    return orthonormError (X, Teuchos::null);
  }

  magnitude_type orthogError (const MV &X1, const MV &X2) const {
    return orthogError (X1, Teuchos::null, X2);
  }

  magnitude_type
  orthogError (const MV &X1,
               Teuchos::RCP<const MV> MX1,
               const MV &X2) const
  {
    if (getOp ().is_null ()) {
      // Ignore MX1, since we don't need to write to it.
      return tsqr_.orthogError (X1, X2);
    } else {
      ensureDgksInit ();
      return pDgks_->orthogError (X1, MX1, X2);
    }
  }

  void
  setLabel (const std::string& label)
  {
    tsqr_.setLabel (label);

    // Make sure DGKS gets the new label, if it's initialized.
    // Otherwise, it will get the new label on initialization.
    if (! pDgks_.is_null ()) {
      pDgks_->setLabel (label);
    }
  }

private:
  //! Ensure that the DGKSOrthoManager instance is initialized.
  void
  ensureDgksInit () const
  {
    // NOTE (mfh 11 Jan 2011) DGKS has a parameter that needs to be
    // set.  For now, we just use the default value of the parameter.
    if (pDgks_.is_null ()) {
      pDgks_ = Teuchos::rcp (new dgks_type (getLabel (), getOp ()));
    }
  }

  /// \brief TSQR + BGS orthogonalization manager implementation.
  ///
  /// We use this when getOp() returns null (so that the inner product
  /// is the default Euclidean inner product).  It gets initialized in
  /// the constructor.  This must be declared mutable because of the
  /// requirements of the MatOrthoManager interface.  tsqr_type's
  /// interface is honest about what methods are really const.
  mutable tsqr_type tsqr_;

  /// \brief DGKS orthogonalization manager.
  ///
  /// We initialize and use this when getOp() does <i>not</i> return
  /// null (could be a non-Euclidean inner product, but not
  /// necessarily).
  mutable Teuchos::RCP<dgks_type> pDgks_;
};

} // namespace Belos

#endif // __BelosTsqrOrthoManager_hpp
