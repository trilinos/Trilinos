//@HEADER
// ************************************************************************
//
//                 Belos: Block Linear Solvers Package
//                  Copyright 2010 Sandia Corporation
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

/// \file BelosTsqrOrthoManager.hpp
/// \brief Orthogonalization manager based on Tall Skinny QR (TSQR)

#ifndef __BelosTsqrOrthoManager_hpp
#define __BelosTsqrOrthoManager_hpp

#include "BelosTsqrOrthoManagerImpl.hpp"
// Belos doesn't have an SVQB orthomanager implemented yet, so we just
// use a default orthomanager for the case where the inner product
// matrix is nontrivial.
#include "BelosDGKSOrthoManager.hpp" 

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

namespace Belos {

  /// \class TsqrOrthoManager
  /// \brief TSQR-based OrthoManager subclass
  ///
  /// This is the actual subclass of OrthoManager, implemented using
  /// TsqrOrthoManagerImpl (TSQR + Block Gram-Schmidt).
  template<class Scalar, class MV>
  class TsqrOrthoManager : public OrthoManager<Scalar, MV> {
  public:
    typedef Scalar scalar_type;
    typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType magnitude_type;
    //! \typedef Multivector type with which this class was specialized
    typedef MV multivector_type;

    typedef Teuchos::SerialDenseMatrix<int, Scalar>           serial_matrix_type;
    typedef Teuchos::RCP<serial_matrix_type>                  serial_matrix_ptr;
    typedef Teuchos::Array<Teuchos::RCP<serial_matrix_type> > prev_coeffs_type;

    /// \brief Get default parameters for TsqrOrthoManager
    ///
    /// Get a (pointer to a) default list of parameters for
    /// configuring a TsqrOrthoManager instance.
    ///
    /// \note To get nondefault behavior, a good thing to do is to
    ///   make a deep copy of the returned parameter list, and then
    ///   modify individual entries as desired.
    ///
    /// \warning This method is not reentrant.  It should only be
    ///   called by one thread at a time.
    ///
    static Teuchos::RCP<const Teuchos::ParameterList> getDefaultParameters() {
      return TsqrOrthoManagerImpl<Scalar, MV>::getDefaultParameters();
    }

    //! Constructor
    TsqrOrthoManager (const Teuchos::RCP<const Teuchos::ParameterList>& params, 
		      const std::string& label = "Belos") :
      impl_ (params, label)
    {}

    virtual ~TsqrOrthoManager() {}

    virtual void 
    innerProd (const MV &X, const MV &Y, serial_matrix_type& Z) const
    {
      return impl_.innerProd (X, Y, Z);
    }

    virtual void 
    norm (const MV& X, std::vector<magnitude_type> &normvec) const
    {
      return impl_.norm (X, normvec);
    }

    virtual void 
    project (MV &X, 
	     prev_coeffs_type C, 
	     Teuchos::Array<Teuchos::RCP<const MV> > Q) const
    {
      return impl_.project (X, C, Q);
    }

    virtual int 
    normalize (MV &X, serial_matrix_ptr B) const
    {
      return impl_.normalize (X, B);
    }

    virtual int 
    projectAndNormalize (MV &X, 
			 prev_coeffs_type C,
			 serial_matrix_ptr B,
			 Teuchos::Array<Teuchos::RCP<const MV> > Q) const
    {
      return impl_.projectAndNormalize (X, C, B, Q);
    }

    /// \brief Normalize X into Q*B, overwriting X with invalid values
    ///
    /// Normalize X into Q*B, overwriting X with invalid values.  
    ///
    /// \note We expose this interface to applications because TSQR is
    ///   not able to compute an orthogonal basis in place; it needs
    ///   scratch space.  Applications can exploit this interface to
    ///   avoid excessive copying of vectors when using TSQR for
    ///   orthogonalization.
    ///
    /// \param X [in/out] Vector(s) to orthogonalize
    /// \param B [out] Orthogonalization coefficients
    ///
    /// \return Rank of X
    ///
    /// \note Q must have at least as many columns as X.  It may have
    /// more columns than X; those columns are ignored.
    int 
    normalizeNoCopy (MV& X, MV& Q, serial_matrix_ptr B) const
    {
      return impl_.normalizeNoCopy (X, Q, B);
    }

    /// \brief Project and normalize X_in into X_out; overwrite X_in
    ///
    /// Project X_in against Q, storing projection coefficients in C,
    /// and normalize X_in into X_out, storing normalization
    /// coefficients in B.  On output, X_out has the resulting
    /// orthogonal vectors and X_in is overwritten with invalid values.
    ///
    /// \param X_in [in/out] On input: The vectors to project against
    ///   Q and normalize.  Overwritten with invalid values on output.
    /// \param X_out [out] On output: the normalized input vectors
    ///   after projection against Q.
    /// \param C [out] The projection coefficients 
    /// \param B [out] The normalization coefficients
    /// \param Q [in] The orthogonal basis against which to project
    ///
    /// \return Rank of X_in after projection
    ///
    /// \note We expose this interface to applications for the same
    ///   reason that we expose normalizeNoCopy().
    int 
    projectAndNormalizeNoCopy (MV& X_in, 
			       MV& X_out,
			       prev_coeffs_type C,
			       serial_matrix_ptr B,
			       Teuchos::Array<Teuchos::RCP<const MV> > Q) const
    {
      return impl_.projectAndNormalizeNoCopy (X_in, X_out, C, B, Q);
    }

    virtual typename Teuchos::ScalarTraits<Scalar>::magnitudeType 
    orthonormError (const MV &X) const
    {
      return impl_.orthonormError (X);
    }

    virtual magnitude_type orthogError (const MV &X1, const MV &X2) const {
      return impl_.orthogError (X1, X2);
    }

    /// Set the label for (the timers for) this orthogonalization
    /// manager, and create new timers if the label has changed.
    ///
    /// \param label [in] New label for timers
    ///
    /// \note Belos::OrthoManager wants this virtual function to be
    ///   implemented; Anasazi::OrthoManager does not.
    void setLabel (const std::string& label) { impl_.setLabel (label); }

    //! Return timers label
    const std::string& getLabel() const { return impl_.getLabel(); }

  private:
    /// "Mutable" because it has internal scratch space state.  I know
    /// it's bad, but it's the only way this class can be part of the
    /// OrthoManager hierarchy.
    mutable TsqrOrthoManagerImpl<Scalar, MV> impl_;

    //! Label for timers (if timers are enabled at compile time)
    std::string label_;
  };


  /// \class TsqrMatOrthoManager
  /// \brief MatOrthoManager subclass using TSQR or DGKS
  ///
  /// Subclass of MatOrthoManager.  When getOp() == null (Euclidean
  /// inner product), uses TSQR + Block Gram-Schmidt for
  /// orthogonalization.  When getOp() != null, uses DGKSOrthoManager
  /// (Classical Gram-Schmidt with reorthogonalization) for
  /// orthogonalization.  Avoids communication only in the TSQR case.
  /// Initialization of either orthogonalization manager is "lazy," so
  /// you don't have to pay for scratch space if you don't use it.
  ///
  template<class Scalar, class MV, class OP>
  class TsqrMatOrthoManager : public MatOrthoManager<Scalar, MV, OP> {
  public:
    typedef Scalar scalar_type;
    typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType magnitude_type;
    //! Multivector type with which this class was specialized
    typedef MV multivector_type;
    //! Operator type with which this class was specialized
    typedef OP operator_type;

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
    typedef Teuchos::RCP<MV>                        mv_ptr;
    typedef Teuchos::RCP<const MV>                  const_mv_ptr;
    typedef Teuchos::Array<const_mv_ptr>            const_prev_mvs_type;
    typedef Teuchos::SerialDenseMatrix<int, Scalar> serial_matrix_type;
    typedef Teuchos::RCP<serial_matrix_type>                  serial_matrix_ptr;
    typedef Teuchos::Array<Teuchos::RCP<serial_matrix_type> > prev_coeffs_type;

    /// \brief Get default parameters for TsqrMatOrthoManager
    ///
    /// Get a (pointer to a) default list of parameters for
    /// configuring a TsqrMatOrthoManager instance.
    ///
    /// \note To get nondefault behavior, a good thing to do is to
    ///   make a deep copy of the returned parameter list, and then
    ///   modify individual entries as desired.
    ///
    /// \warning This method is not reentrant.  It should only be
    ///   called by one thread at a time.
    ///
    static Teuchos::RCP<const Teuchos::ParameterList> getDefaultParameters() {
      // FIXME (mfh 11 Jan 2011) What about DGKS parameters?
      return TsqrOrthoManagerImpl<Scalar, MV>::getDefaultParameters();
    }

    //! Default constructor (sets Op to Teuchos::null)
    TsqrMatOrthoManager () :
      MatOrthoManager<Scalar, MV, OP>(Teuchos::null),
      pTsqr_ (Teuchos::null), // Lazy initialization
      pDgks_ (Teuchos::null)  // Lazy initialization
    {}

    /// \brief Set label for timers (if timers enabled)
    ///
    /// \note Belos::OrthoManager wants this virtual function to be
    ///   implemented; Anasazi::OrthoManager does not.
    void setLabel (const std::string& label) {
      label_ = label; 
    }
    //! Return label for timers (if timers enabled)
    const std::string& getLabel() const { return label_; }

    /// \brief Constructor
    ///
    /// \param params [in] Parameters used to set up the
    ///   orthogonalization.  Call the getDefaultParameters() class
    ///   method for default parameters and their documentation.
    ///
    /// \param label [in] Label for timers (if timers are used) 
    ///
    /// \param Op [in] Inner product with respect to which to
    ///   orthogonalize vectors.  If Teuchos::null, use the Euclidean
    ///   inner product.
    TsqrMatOrthoManager (const Teuchos::RCP<const Teuchos::ParameterList>& params,
			 const std::string& label = "Belos",
			 Teuchos::RCP< const OP > Op = Teuchos::null) :
      MatOrthoManager<Scalar, MV, OP>(Op),
      params_ (params),
      label_ (label),
      pTsqr_ (Teuchos::null), // Lazy initialization
      pDgks_ (Teuchos::null)  // Lazy initialization
    {}

    //! \brief Destructor
    virtual ~TsqrMatOrthoManager() {}

    /// \brief Return the inner product operator
    ///
    /// Return the inner product operator used for orthogonalization.
    /// If it is Teuchos::null, then the vectors are orthogonalized
    /// with respect to the Euclidean inner product.
    ///
    /// \note We override the base class' setOp() so that the
    ///   DGKSOrthoManager gets the new op.
    virtual void 
    setOp (Teuchos::RCP<const OP> Op) 
    {
      // We use this notation to help C++ resolve the name.
      // Otherwise, it won't know where to find setOp(), since this is
      // a member function of the base class which does not depend on
      // the template parameters.
      base_type::setOp (Op); // base class gets a copy of the Op too
      ensureDgksInit (); // Make sure the DGKS object has been initialized
      pDgks_->setOp (Op);
    }

    /// Return the inner product operator, if any
    ///
    /// \note We override only to help C++ do name lookup in the other
    ///   member functions.
    virtual Teuchos::RCP<const OP> getOp () const { 
      return base_type::getOp(); 
    }

    /// Project X against Q with respect to the inner product computed
    /// by \c innerProd().  Store the resulting coefficients in C.  If
    /// MX is not null, assume that MX is the result of applying the
    /// operator to X, and exploit this when computing the inner
    /// product.
    virtual void 
    project (MV &X, 
	     Teuchos::RCP<MV> MX,
	     prev_coeffs_type C,
	     const_prev_mvs_type Q) const
    {
      if (getOp().is_null())
	{
	  ensureTsqrInit ();
	  pTsqr_->project (X, C, Q);
	  if (! MX.is_null())
	    // MX gets a copy of X; M is the identity operator.
	    MVT::Assign (X, *MX);
	}
      else
	{
	  ensureDgksInit ();
	  pDgks_->project (X, MX, C, Q);
	}
    }

    /// Project X against Q with respect to the inner product computed
    /// by \c innerProd().  Store the resulting coefficients in C.
    virtual void 
    project (MV &X, 
	     prev_coeffs_type C,
	     const_prev_mvs_type Q) const
    {
      project (X, Teuchos::null, C, Q);
    }

    virtual int 
    normalize (MV& X, 
	       Teuchos::RCP<MV> MX,
	       Teuchos::RCP<Teuchos::SerialDenseMatrix<int,Scalar> > B) const
    {
      if (getOp().is_null())
	{
	  ensureTsqrInit ();
	  const int rank = pTsqr_->normalize (X, B);
	  if (! MX.is_null())
	    // MX gets a copy of X; M is the identity operator.
	    MVT::Assign (X, *MX);
	  return rank;
	}
      else
	{
	  ensureDgksInit ();
	  return pDgks_->normalize (X, MX, B);
	}
    }

    virtual int
    normalize (MV& X, 
	       Teuchos::RCP<Teuchos::SerialDenseMatrix<int,Scalar> > B) const 
    {
      return normalize (X, Teuchos::null, B);
    }


    virtual int 
    projectAndNormalize (MV &X, 
			 Teuchos::RCP<MV> MX,
			 Teuchos::Array<serial_matrix_ptr> C,
			 serial_matrix_ptr B, 
			 Teuchos::Array<Teuchos::RCP<const MV> > Q ) const
    {
      if (getOp().is_null())
	{
	  ensureTsqrInit ();
	  const int rank = pTsqr_->projectAndNormalize (X, C, B, Q); 
	  if (! MX.is_null())
	    // MX gets a copy of X; M is the identity operator.
	    MVT::Assign (X, *MX);
	  return rank;
	}
      else
	{
	  ensureDgksInit ();
	  return pDgks_->projectAndNormalize (X, MX, C, B, Q);
	}
    }

    /// \brief Normalize X into Q*B, overwriting X with invalid values
    ///
    /// Normalize X into Q*B, overwriting X with invalid values.  
    ///
    /// \note We expose this interface to applications because TSQR is
    ///   not able to compute an orthogonal basis in place; it needs
    ///   scratch space.  Applications can exploit this interface to
    ///   avoid excessive copying of vectors when using TSQR for
    ///   orthogonalization.
    ///
    /// \param X [in/out] Vector(s) to orthogonalize
    /// \param B [out] Orthogonalization coefficients
    ///
    /// \return Rank of X
    ///
    /// \note Q must have at least as many columns as X.  It may have
    /// more columns than X; those columns are ignored.
    int 
    normalizeNoCopy (MV& X, MV& Q, serial_matrix_ptr B) const
    {
      if (getOp().is_null())
	{
	  ensureTsqrInit ();
	  return pTsqr_->normalizeNoCopy (X, Q, B);
	}
      else
	{
	  // DGKS normalizes in place, so we have to copy.
	  ensureDgksInit ();
	  const int rank = pDgks_->normalize (X, B);
	  MVT::Assign (X, Q);
	  return rank;
	}
    }

    /// \brief Project and normalize X_in into X_out; overwrite X_in
    ///
    /// Project X_in against Q, storing projection coefficients in C,
    /// and normalize X_in into X_out, storing normalization
    /// coefficients in B.  On output, X_out has the resulting
    /// orthogonal vectors and X_in is overwritten with invalid values.
    ///
    /// \param X_in [in/out] On input: The vectors to project against
    ///   Q and normalize.  Overwritten with invalid values on output.
    /// \param X_out [out] On output: the normalized input vectors
    ///   after projection against Q.
    /// \param C [out] The projection coefficients 
    /// \param B [out] The normalization coefficients
    /// \param Q [in] The orthogonal basis against which to project
    ///
    /// \return Rank of X_in after projection
    ///
    /// \note We expose this interface to applications for the same
    ///   reason that we expose normalizeNoCopy().
    int 
    projectAndNormalizeNoCopy (MV& X_in, 
			       MV& X_out,
			       prev_coeffs_type C,
			       serial_matrix_ptr B,
			       const_prev_mvs_type Q) const
    {
      if (getOp().is_null())
	{
	  ensureTsqrInit ();
	  return pTsqr_->projectAndNormalizeNoCopy (X_in, X_out, C, B, Q);
	}
      else
	{
	  // DGKS normalizes in place, so we have to copy.
	  ensureDgksInit ();
	  const int rank = pDgks_->projectAndNormalize (X_in, Teuchos::null, C, B, Q);
	  MVT::Assign (X_in, X_out);
	  return rank;
	}
    }

    virtual magnitude_type
    orthonormError (const MV &X,
		    Teuchos::RCP<const MV> MX) const
    {
      if (getOp().is_null())
	{
	  ensureTsqrInit ();
	  return pTsqr_->orthonormError (X); // Ignore MX
	}
      else
	{
	  ensureDgksInit ();
	  return pDgks_->orthonormError (X, MX);
	}
    }

    virtual magnitude_type
    orthonormError (const MV &X) const
    {
      return orthonormError (X, Teuchos::null);
    }

    virtual magnitude_type
    orthogError (const MV &X1, 
		 const MV &X2) const
    {
      return orthogError (X1, Teuchos::null, X2);
    }

    virtual magnitude_type
    orthogError (const MV &X1, 
		 Teuchos::RCP<const MV> MX1,
		 const MV &X2) const
    {
      if (getOp().is_null())
	{
	  ensureTsqrInit ();
	  // Ignore MX1, since we don't need to write to it.
	  return pTsqr_->orthogError (X1, X2);
	}
      else
	{
	  ensureDgksInit ();
	  return pDgks_->orthogError (X1, MX1, X2);
	}
    }

  private:
    void
    ensureTsqrInit () const
    {
      if (pTsqr_.is_null())
	pTsqr_ = Teuchos::rcp (new tsqr_type (params_, getLabel()));
    }
    void 
    ensureDgksInit () const
    {
      // FIXME (mfh 11 Jan 2011) 
      //
      // DGKS has a parameter that needs to be set.
      if (pDgks_.is_null())
	pDgks_ = Teuchos::rcp (new dgks_type (getLabel(), getOp()));
    }

    //! Parameter list for initializing the orthogonalization
    Teuchos::RCP<const Teuchos::ParameterList> params_;
    //! Label for timers (if timers are used)
    std::string label_;
    /// TSQR + BGS orthogonalization manager implementation, used when
    /// getOp() == null (Euclidean inner product).
    mutable Teuchos::RCP<tsqr_type> pTsqr_;
    /// DGKS orthogonalization manager, used when getOp() != null
    /// (could be a non-Euclidean inner product, but not necessarily).
    mutable Teuchos::RCP<dgks_type> pDgks_;
  };

} // namespace Belos

#endif // __BelosTsqrOrthoManager_hpp

