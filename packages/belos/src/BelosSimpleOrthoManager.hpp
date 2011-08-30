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

/// \file BelosSimpleOrthoManager.hpp
/// \brief Simple OrthoManager implementation for benchmarks.
///
#ifndef __Belos_SimpleOrthoManager_hpp
#define __Belos_SimpleOrthoManager_hpp

#include <BelosConfigDefs.hpp>
#include <BelosMultiVecTraits.hpp>
#include <BelosOrthoManager.hpp>
#include <BelosOutputManager.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_StandardCatchMacros.hpp>
#include <Teuchos_TimeMonitor.hpp>

namespace Belos {

  /// \class SimpleOrthoManager
  /// \brief Simple OrthoManager implementation for benchmarks.
  /// \author Mark Hoemmen
  ///
  /// This is a very simple orthogonalization method and should only
  /// be used for benchmarks.  It performs optional unconditional
  /// reorthogonalization (no norm tests), but has no rank-revealing
  /// features.
  template<class Scalar, class MV>
  class SimpleOrthoManager : public OrthoManager<Scalar, MV> {
  public:
    typedef Scalar scalar_type;
    typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType magnitude_type;
    typedef Teuchos::SerialDenseMatrix<int, Scalar> mat_type;
    typedef Teuchos::RCP<Teuchos::SerialDenseMatrix<int, Scalar> > mat_ptr;
    
  private:
    typedef MultiVecTraits<Scalar, MV> MVT;
    typedef Teuchos::ScalarTraits<Scalar> STS;
    typedef Teuchos::ScalarTraits<magnitude_type> STM;

    //! Label for Belos timer display
    std::string label_;
    //! Output manager (used mainly for debugging)
    Teuchos::RCP<OutputManager<Scalar> > outMan_;
    //! Whether or not to do (unconditional) reorthogonalization
    bool reorthogonalize_;
    //! Whether to use MGS or CGS in the normalize() step
    bool useMgs_;

#ifdef BELOS_TEUCHOS_TIME_MONITOR
    //! Timer for all orthogonalization operations
    Teuchos::RCP<Teuchos::Time> timerOrtho_;
    //! Timer for projection operations
    Teuchos::RCP<Teuchos::Time> timerProject_;
    //! Timer for normalization operations
    Teuchos::RCP<Teuchos::Time> timerNormalize_;

    /// \brief Instantiate and return a timer with an appropriate label.
    ///
    /// \param prefix [in] Prefix for the timer label, e.g., "Belos"
    /// \param timerName [in] Name of the timer, or what the timer
    ///   is timing, e.g., "Projection" or "Normalization"
    ///
    /// \return Smart pointer to a new Teuchos::Time timer object,
    ///   to be used via Teuchos::TimeMonitor
    static Teuchos::RCP<Teuchos::Time>
    makeTimer (const std::string& prefix, 
	       const std::string& timerName)
    {
      const std::string timerLabel = 
	prefix.empty() ? timerName : (prefix + ": " + timerName);
      return Teuchos::TimeMonitor::getNewTimer (timerLabel);
    }
#endif // BELOS_TEUCHOS_TIME_MONITOR
    
  public:

    /// \brief Get a default list of parameters.
    ///
    /// The "default" parameter list sets reasonably safe options in
    /// terms of accuracy of the computed orthogonalization.  Call \c
    /// getFastParameters() if you prefer to sacrifice some accuracy
    /// for speed.
    ///
    /// \warning This class method is nonreentrant.
    ///
    static Teuchos::RCP<const Teuchos::ParameterList> 
    getDefaultParameters ()
    {
      using Teuchos::ParameterList;
      using Teuchos::RCP;

      // This part makes this class method non-reentrant.
      static RCP<ParameterList> params;
      if (! params.is_null())
	return params;

      params = Teuchos::parameterList();
      const std::string defaultNormalizationMethod ("MGS");
      params->set ("Normalization", defaultNormalizationMethod,
		   "Which normalization method to use. Valid values are \"MGS\""
		   " (for Modified Gram-Schmidt) and \"CGS\" (for Classical "
		   "Gram-Schmidt).");
      const bool defaultReorthogonalization = false;
      params->set ("Reorthogonalization", defaultReorthogonalization,
		   "Whether to perform one (unconditional) reorthogonalization "
		   "pass.");
      return params;
    }

    /// \brief Get a "fast" list of parameters.
    ///
    /// The "fast" parameter list favors speed of orthogonalization,
    /// but sacrifices some accuracy.  Call \c getDefaultParameters()
    /// for safer options in terms of accuracy.
    ///
    /// \warning This class method is nonreentrant.
    ///
    static Teuchos::RCP<const Teuchos::ParameterList> 
    getFastParameters ()
    {
      using Teuchos::ParameterList;
      using Teuchos::RCP;
      using Teuchos::rcp;

      // This part makes this class method non-reentrant.
      static RCP<ParameterList> fastParams;
      if (! fastParams.is_null())
	return fastParams;
      
      fastParams = rcp (new ParameterList (*getDefaultParameters()));
      const std::string fastNormalizationMethod ("CGS");
      fastParams->set ("Normalization", fastNormalizationMethod);
      const bool fastReorthogonalization = false;
      fastParams->set ("Reorthogonalization", fastReorthogonalization);
      return fastParams;
    }

    /// \brief Constructor
    ///
    /// \param outMan [in/out] Output manager.  If not null, use for
    ///   various kinds of status output (in particular, for debugging).
    ///
    /// \param label [in] Label for Belos timers.
    ///
    /// \param params [in] List of configuration parameters.  Call
    ///   getDefaultParameters() or getFastParameters() for sample
    ///   valid parameter lists.
    ///
    SimpleOrthoManager (const Teuchos::RCP<OutputManager<Scalar> >& outMan,
			const std::string& label,
			const Teuchos::RCP<const Teuchos::ParameterList>& params) :
      label_ (label),
      outMan_ (outMan)
    {
      using Teuchos::ParameterList;
      using Teuchos::RCP;
      using Teuchos::Exceptions::InvalidParameter;

#ifdef BELOS_TEUCHOS_TIME_MONITOR
      timerOrtho_ = makeTimer (label, "All orthogonalization");
      timerProject_ = makeTimer (label, "Projection");
      timerNormalize_ = makeTimer (label, "Normalization");
#endif // BELOS_TEUCHOS_TIME_MONITOR
      
      std::string normalizeImpl;
      if (params.is_null())
	normalizeImpl = "MGS";
      else
	{
	  try {
	    normalizeImpl = params->get<std::string>("Normalization");
	  } catch (InvalidParameter&) {
	    normalizeImpl = "MGS";
	  }
	}
      bool reorthogonalize;
      if (params.is_null())
	reorthogonalize = false;
      else
	{
	  try {
	    reorthogonalize = params->get<bool>("Reorthogonalization");
	  } catch (InvalidParameter&) {
	    reorthogonalize = false;
	  }
	}
      if (normalizeImpl == "MGS" || normalizeImpl == "Mgs" || normalizeImpl == "mgs")
	useMgs_ = true;
      else 
	useMgs_ = false;
      reorthogonalize_ = reorthogonalize;

      if (! outMan_.is_null())
	{
	  using std::endl;
	  std::ostream& dbg = outMan_->stream(Debug);
	  dbg << "Belos::SimpleOrthoManager constructor:" << endl
	      << "-- Normalization method: " 
	      << (useMgs_ ? "MGS" : "CGS") << endl
	      << "-- Reorthogonalize (unconditionally)? " 
	      << (reorthogonalize_ ? "Yes" : "No") << endl;
	}
    }
    
    //! Virtual destructor for memory safety of derived classes.
    virtual ~SimpleOrthoManager() {}
      
    void innerProd (const MV &X, const MV &Y, mat_type& Z ) const {
      MVT::MvTransMv (STS::one(), X, Y, Z);
    }

    void norm (const MV& X, std::vector<magnitude_type>& normVec) const {
      const int numCols = MVT::GetNumberVecs (X);
      // std::vector<T>::size_type is unsigned; int is signed.  Mixed
      // unsigned/signed comparisons trigger compiler warnings.
      if (normVec.size() < static_cast<size_t>(numCols))
	normVec.resize (numCols); // Resize normvec if necessary.
      MVT::MvNorm (X, normVec);
    }

    void 
    project (MV &X, 
	     Teuchos::Array<mat_ptr> C,
	     Teuchos::ArrayView<Teuchos::RCP<const MV> > Q) const
    {
#ifdef BELOS_TEUCHOS_TIME_MONITOR
      Teuchos::TimeMonitor timerMonitorOrtho(*timerOrtho_);
      Teuchos::TimeMonitor timerMonitorProject(*timerProject_);
#endif // BELOS_TEUCHOS_TIME_MONITOR

      allocateProjectionCoefficients (C, Q, X, true);
      rawProject (X, Q, C);
      if (reorthogonalize_) // Unconditional reorthogonalization
	{
	  Teuchos::Array<Teuchos::RCP<Teuchos::SerialDenseMatrix<int,Scalar> > > C2;
	  allocateProjectionCoefficients (C2, Q, X, false);
	  for (int k = 0; k < Q.size(); ++k)
	    *C[k] += *C2[k];
	}
    }

    int 
    normalize (MV &X, mat_ptr B) const
    {
#ifdef BELOS_TEUCHOS_TIME_MONITOR
      Teuchos::TimeMonitor timerMonitorOrtho(*timerOrtho_);
      Teuchos::TimeMonitor timerMonitorProject(*timerNormalize_);
#endif // BELOS_TEUCHOS_TIME_MONITOR

      if (useMgs_)
	return normalizeMgs (X, B);
      else
	return normalizeCgs (X, B);
    }

    int 
    projectAndNormalize (MV &X, 
			 Teuchos::Array<mat_ptr> C,
			 mat_ptr B,
			 Teuchos::ArrayView<Teuchos::RCP<const MV> > Q) const
    {
      // Don't need time monitors here: project() and normalize() have
      // their own.
      project (X, C, Q);
      return normalize (X, B);
    }

    magnitude_type 
    orthonormError(const MV &X) const 
    {
      const Scalar ONE = STS::one();
      const int ncols = MVT::GetNumberVecs(X);
      mat_type XTX (ncols, ncols);
      innerProd (X, X, XTX);
      for (int k = 0; k < ncols; ++k)
	XTX(k,k) -= ONE;
      return XTX.normFrobenius();
    }

    magnitude_type
    orthogError(const MV &X1, const MV &X2) const
    {
      const int ncols_X1 = MVT::GetNumberVecs (X1);
      const int ncols_X2 = MVT::GetNumberVecs (X2);
      mat_type X1_T_X2 (ncols_X1, ncols_X2);
      innerProd (X1, X2, X1_T_X2);
      return X1_T_X2.normFrobenius();
    }

    void setLabel (const std::string& label) { label_ = label; }
    const std::string& getLabel() const { return label_; }

  private:

    int
    normalizeMgs (MV &X, 
		  Teuchos::RCP<Teuchos::SerialDenseMatrix<int,Scalar> > B) const
    {
      using Teuchos::Range1D;
      using Teuchos::RCP;
      using Teuchos::rcp;
      using Teuchos::View;

      const int numCols = MVT::GetNumberVecs (X);
      if (numCols == 0)
	return 0;

      if (B.is_null())
	B = rcp (new mat_type (numCols, numCols));
      else if (B->numRows() != numCols || B->numCols() != numCols)
	B->shape (numCols, numCols);

      // Modified Gram-Schmidt orthogonalization
      std::vector<magnitude_type> normVec (1);
      for (int j = 0; j < numCols; ++j)
	{
	  RCP<MV> X_cur = MVT::CloneViewNonConst (X, Range1D(j, j));
	  MV& X_j = *X_cur;
	  for (int i = 0; i < j; ++i)
	    {
	      RCP<const MV> X_prv = MVT::CloneView (X, Range1D(i, i));
	      const MV& X_i = *X_prv;
	      mat_type B_ij (View, *B, 1, 1, i, j);
	      innerProd (X_i, X_j, B_ij);
	      MVT::MvTimesMatAddMv (-STS::one(), X_i, B_ij, STS::one(), X_j);
	      if (reorthogonalize_) // Unconditional reorthogonalization
		{		    
		  // innerProd() overwrites B(i,j), so save the
		  // first-pass projection coefficient and update
		  // B(i,j) afterwards.
		  const Scalar B_ij_first = (*B)(i, j);
		  innerProd (X_i, X_j, B_ij);
		  MVT::MvTimesMatAddMv (-STS::one(), X_i, B_ij, STS::one(), X_j);
		  (*B)(i, j) += B_ij_first;
		}
	    }
	  // Normalize column j of X
	  norm (X_j, normVec);
	  const magnitude_type theNorm = normVec[0];
	  (*B)(j, j) = theNorm;
	  if (normVec[0] != STM::zero())
	    MVT::MvScale (X_j, STS::one() / theNorm);
	  else
	    return j; // break out early
	}
      return numCols; // full rank, as far as we know
    }


    int 
    normalizeCgs (MV &X, mat_ptr B) const
    {
      using Teuchos::Range1D;
      using Teuchos::RCP;
      using Teuchos::rcp;
      using Teuchos::View;

      const int numCols = MVT::GetNumberVecs (X);
      if (numCols == 0)
	return 0;

      if (B.is_null())
	B = rcp (new mat_type (numCols, numCols));
      else if (B->numRows() != numCols || B->numCols() != numCols)
	B->shape (numCols, numCols);
      mat_type& B_ref = *B;

      // Classical Gram-Schmidt orthogonalization 
      std::vector<magnitude_type> normVec (1);

      // Space for reorthogonalization
      mat_type B2 (numCols, numCols);

      // Do the first column first.
      {
	RCP<MV> X_cur = MVT::CloneViewNonConst (X, Range1D(0, 0));
	// Normalize column 0 of X
	norm (*X_cur, normVec);
	const magnitude_type theNorm = normVec[0];
	B_ref(0,0) = theNorm;
	if (theNorm != STM::zero())
	  {
	    const Scalar invNorm = STS::one() / theNorm;
	    MVT::MvScale (*X_cur, invNorm);
	  }
	else
	  return 0; // break out early
      }

      // Orthogonalize the remaining columns of X
      for (int j = 1; j < numCols; ++j)
	{
	  RCP<MV> X_cur = MVT::CloneViewNonConst (X, Range1D(j, j));
	  RCP<const MV> X_prv = MVT::CloneView (X, Range1D(0, j-1));	      
	  mat_type B_prvcur (View, B_ref, j, 1, 0, j);

	  // Project X_cur against X_prv (first pass)
	  innerProd (*X_prv, *X_cur, B_prvcur);
	  MVT::MvTimesMatAddMv (-STS::one(), *X_prv, B_prvcur, STS::one(), *X_cur);
	  // Unconditional reorthogonalization: 
	  // project X_cur against X_prv (second pass)
	  if (reorthogonalize_) 
	    {		    
	      mat_type B2_prvcur (View, B2, j, 1, 0, j);
	      innerProd (*X_prv, *X_cur, B2_prvcur);
	      MVT::MvTimesMatAddMv (-STS::one(), *X_prv, B2_prvcur, STS::one(), *X_cur);
	      B_prvcur += B2_prvcur;
	    }
	  // Normalize column j of X
	  norm (*X_cur, normVec);
	  const magnitude_type theNorm = normVec[0];
	  B_ref(j,j) = theNorm;
	  if (theNorm != STM::zero())
	    {
	      const Scalar invNorm = STS::one() / theNorm;
	      MVT::MvScale (*X_cur, invNorm);
	    }
	  else
	    return j; // break out early
	}
      return numCols; // full rank, as far as we know
    }


    void
    allocateProjectionCoefficients (Teuchos::Array<mat_ptr>& C, 
				    Teuchos::ArrayView<Teuchos::RCP<const MV> > Q,
				    const MV& X,
				    const bool attemptToRecycle = true) const
    {
      using Teuchos::rcp;

      const int num_Q_blocks = Q.size();
      const int ncols_X = MVT::GetNumberVecs (X);
      C.resize (num_Q_blocks);
      // # of block(s) that had to be (re)allocated (either allocated
      // freshly, or resized).
      int numAllocated = 0;
      if (attemptToRecycle)
	{
	  for (int i = 0; i < num_Q_blocks; ++i) 
	    {
	      const int ncols_Qi = MVT::GetNumberVecs (*Q[i]);
	      // Create a new C[i] if necessary, otherwise resize if
	      // necessary, otherwise fill with zeros.
	      if (C[i].is_null())
		{
		  C[i] = rcp (new mat_type (ncols_Qi, ncols_X));
		  numAllocated++;
		}
	      else 
		{
		  mat_type& Ci = *C[i];
		  if (Ci.numRows() != ncols_Qi || Ci.numCols() != ncols_X)
		    {
		      Ci.shape (ncols_Qi, ncols_X);
		      numAllocated++;
		    }
		  else
		    Ci.putScalar (STS::zero());
		}
	    }
	}
      else // Just allocate; don't try to check if we can recycle
	{
	  for (int i = 0; i < num_Q_blocks; ++i) 
	    {
	      const int ncols_Qi = MVT::GetNumberVecs (*Q[i]);
	      C[i] = rcp (new mat_type (ncols_Qi, ncols_X));
	      numAllocated++;
	    }
	}
      if (! outMan_.is_null())
	{
	  using std::endl;
	  std::ostream& dbg = outMan_->stream(Debug);
	  dbg << "SimpleOrthoManager::allocateProjectionCoefficients: " 
	      << "Allocated " << numAllocated << " blocks out of " 
	      << num_Q_blocks << endl;
	}
    }

    void
    rawProject (MV& X, 
		Teuchos::ArrayView<Teuchos::RCP<const MV> > Q,
		Teuchos::ArrayView<mat_ptr> C) const
    {	
      // "Modified Gram-Schmidt" version of Block Gram-Schmidt.
      const int num_Q_blocks = Q.size();
      for (int i = 0; i < num_Q_blocks; ++i)
	{
	  mat_type& Ci = *C[i];
	  const MV& Qi = *Q[i];
	  innerProd (Qi, X, Ci);
	  MVT::MvTimesMatAddMv (-STS::one(), Qi, Ci, STS::one(), X);
	}
    }

  };
} // namespace Belos
