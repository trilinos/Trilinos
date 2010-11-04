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

#include "BelosConfigDefs.hpp" // HAVE_BELOS_TSQR
#include "BelosMultiVecTraits.hpp"
// Belos doesn't have an SVQB orthomanager implemented yet, so we just
// use a default orthomanager for the case where the inner product
// matrix is nontrivial.
#include "BelosDGKSOrthoManager.hpp" 
#include "Teuchos_LAPACK.hpp"
#include "Teuchos_ParameterList.hpp"

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

namespace Belos {

  /// \brief Get default parameters for Tsqr(Mat)OrthoManager
  ///
  /// Get a (pointer to a) default list of parameters for setting up
  /// TsqrOrthoManager or TsqrMatOrthoManager (the same parameters
  /// work for both).  
  template< class MagnitudeType >
  Teuchos::RCP< const Teuchos::ParameterList >
  getDefaultTsqrParameters ()
  {
    using Teuchos::RCP;
    typedef Teuchos::ScalarTraits< MagnitudeType > SCTM;

    // FIXME (mfh 28 Oct 2010) The use of persistent data means that
    // this routine is NOT reentrant.  That means that it is not
    // thread-safe, for example.  One way to make this routine
    // reentrant would be to use a construct like pthread_once().
    static RCP< const Teuchos::ParameterList > defaultParams_;
    if (Teuchos::is_null (defaultParams_))
      {
	RCP< Teuchos::ParameterList > params = Teuchos::parameterList();

	const bool defaultRandomizeNullSpace = false;
	params->set ("randomizeNullSpace", defaultRandomizeNullSpace, 
		     "Whether to fill in null space vectors with random data.");
	const bool defaultReorthogonalizeBlocks = false;
	params->set ("reorthogonalizeBlocks", defaultReorthogonalizeBlocks,
		     "Whether to do block reorthogonalization at all.");
	// This parameter corresponds to the "blk_tol_" parameter in
	// Belos' DGKSOrthoManager.  We choose the same default value.
	const MagnitudeType defaultBlockReorthogThreshold = 
	  MagnitudeType(10) * SCTM::squareroot (SCTM::eps());
	params->set ("blockReorthogThreshold", defaultBlockReorthogThreshold, 
		     "If reorthogonalizeBlocks==true, and if the norm of "
		     "any column within a block decreases by this much or "
		     "more after orthogonalization, we reorthogonalize the "
		     "block.");
	// This parameter corresponds to the "sing_tol_" parameter in
	// Belos' DGKSOrthoManager.  We choose the same default value.
	const MagnitudeType defaultRelativeRankTolerance = 
	  MagnitudeType(10) * SCTM::eps();
	// If zero, then we will always declare blocks to be
	// numerically full rank, as long as no singular values are
	// zero.
	params->set ("relativeRankTolerance", defaultRelativeRankTolerance,
		     "Relative tolerance to determine the numerical rank "
		     "of a block, in the normalize() method.");
	const bool defaultThrowOnReorthogFault = true;
	// See Stewart's 2008 paper on block Gram-Schmidt for an explanation.
	params->set ("throwOnReorthogFault", defaultThrowOnReorthogFault,
		     "Whether to throw an exception if a "
		     "reorthogonalization fault occurs.  Handling the "
		     "fault correctly can be expensive and the benefits "
		     "of doing so for a Krylov method aren't clear.");
	defaultParams_ = params;
      }
    return defaultParams_;
  }

  /// \class TsqrOrthoError
  /// \brief TsqrOrthoManager(Impl) error
  class TsqrOrthoError : public OrthoError
  {
  public: 
    TsqrOrthoError (const std::string& what_arg) : 
      OrthoError (what_arg) {}
  };

  /// \class TsqrOrthoFault
  /// \brief Orthogonalization fault
  ///
  /// Stewart (SISC 2008) presents a Block Gram-Schmidt (BGS)
  /// algorithm with careful reorthogonalization.  He defines an
  /// "orthogonalization fault" as happening when the second BGS pass
  /// does not succeed.  This is possible in BGS, but not possible in
  /// (non-block) Gram-Schmidt, if you use Stewart's randomization
  /// procedure for the latter.  Stewart gives an algorithm for
  /// recovering from an orthogonalization fault, but the algorithm is
  /// expensive: it involves careful reorthogonalization with
  /// non-block Gram-Schmidt.  If the "throwOnReorthogFault" option is
  /// set, we choose instead to report the orthogonalization fault as
  /// an exception.
  ///
  /// \note This is not a (subclass of) TsqrOrthoError, because the
  ///   latter is a logic or runtime bug, whereas a TsqrOrthoFault is
  ///   a property of the input and admits recovery.
  class TsqrOrthoFault : public OrthoError
  {
  public: 
    TsqrOrthoFault (const std::string& what_arg) : 
      OrthoError (what_arg) {}
  };


  class NonNullOperatorError : public OrthoError
  {
  public:
    NonNullOperatorError () : 
      OrthoError ("Sorry, TsqrOrthoManager doesn\'t work with a non-null Op "
		  "argument.  I know this is bad class design, but it will "
		  "have to do for now, since TsqrOrthoManager has to inherit "
		  "from MatOrthoManager in order to work in Belos\'s "
		  "eigensolvers.  If you want to solve this problem yourself, "
		  "the thing to do is to have TsqrOrthoManager degrade to "
		  "DGKSOrthoManager when this->getOp() != Teuchos::null.")
    {}
  };

  /// \class TsqrOrthoManagerImpl
  /// \brief TSQR-based OrthoManager subclass implementation
  /// 
  /// TsqrOrthoManagerImpl implements the interface defined by
  /// OrthoManager.  It doesn't actually inherit from OrthoManager,
  /// which gives us a bit more freedom when defining the actual
  /// subclass of OrthoManager (TsqrOrthoManager).  
  ///
  /// This class uses a combination of Tall Skinny QR (TSQR) and Block
  /// Gram-Schmidt (BGS) to orthogonalize multivectors.
  ///
  /// The Block Gram-Schmidt procedure used here is inspired by that
  /// of G. W. Stewart ("Block Gram-Schmidt Orthogonalization", SISC
  /// vol 31 #1 pp. 761--775, 2008), except that we use TSQR+SVD
  /// instead of standard Gram-Schmidt with orthogonalization to
  /// handle the current block.  "Orthogonalization faults" may still
  /// happen, but we do not handle them by default.  Rather, we make
  /// one BGS pass, do TSQR+SVD, check the resulting column norms, and
  /// make a second BGS pass (+ TSQR+SVD) if necessary.  If we then
  /// detect an orthogonalization fault, we throw TsqrOrthoFault.
  ///
  template< class ScalarType, class MV >
  class TsqrOrthoManagerImpl
  {
  private:
    typedef typename Teuchos::ScalarTraits< ScalarType >::magnitudeType MagnitudeType;
    typedef Teuchos::ScalarTraits< ScalarType >    SCT;
    typedef Teuchos::ScalarTraits< MagnitudeType > SCTM;
    typedef MultiVecTraits< ScalarType, MV >       MVT;

    typedef typename MVT::tsqr_adaptor_type tsqr_adaptor_type;
    typedef Teuchos::RCP< tsqr_adaptor_type > tsqr_adaptor_ptr;

  public:

    /// Constructor
    ///
    /// \param tsqrParams [in] Configuration parameters, both for this
    ///   orthogonalization manager, and for TSQR itself.  See the TSQR
    ///   documentation for how to set the latter.  They depend on
    ///   which multivector (MV) class you are using (since each MV
    ///   class maps to a specific TSQR implementation).
    ///
    /// \param label [in] Label for Belos timers.  Only has an effect 
    ///   if the compile-time option for enabling Belos timers is set.
    ///
    TsqrOrthoManagerImpl (const Teuchos::ParameterList& tsqrParams,
			  const std::string& label) :
      tsqrParams_ (tsqrParams),
      label_ (label),
      tsqrAdaptor_ (Teuchos::null),
      Q_ (Teuchos::null),
      eps_ (SCT::eps()),
      randomizeNullSpace_ (false),    // Set later by readParams()
      reorthogonalizeBlocks_ (false), // Set later by readParams()
      throwOnReorthogFault_ (true),   // Set later by readParams()
      blockReorthogThreshold_ (0),    // Set later by readParams()
      relativeRankTolerance_ (0)      // Set later by readParams()
    {
#ifdef BELOS_TEUCHOS_TIME_MONITOR
      std::string sublabel = label_ = ": All orthogonalization";
      timerOrtho_ = Teuchos::TimeMonitor::getNewTimer (sublabel);
      sublabel = label_ = ": Projection";
      timerProject_ = Teuchos::TimeMonitor::getNewTimer (sublabel);
      sublabel = label_ = ": Normalization";
      timerNormalize_ = Teuchos::TimeMonitor::getNewTimer (sublabel);
#endif // BELOS_TEUCHOS_TIME_MONITOR

      // Extract values for the four parameters
      // reorthogonalizeBlocks_, throwOnReorthogFault_,
      // blockReorthogThreshold_, and relativeRankTolerance_ from the
      // given parameter list.  Use default values if none are
      // provided.  Other parameters in tsqrParams (e.g.,
      // "cacheBlockSize") get passed along to the underlying TSQR
      // implementation.
      readParams (tsqrParams);
    }

    typedef Teuchos::RCP< MV >                                   mv_ptr;
    typedef Teuchos::RCP< const MV >                             const_mv_ptr;
    typedef Teuchos::Array< const_mv_ptr >                       const_prev_mvs_type;
    typedef Teuchos::SerialDenseMatrix< int, ScalarType >        serial_matrix_type;
    typedef Teuchos::RCP< serial_matrix_type >                   serial_matrix_ptr;
    typedef Teuchos::Array< Teuchos::RCP< serial_matrix_type > > prev_coeffs_type;

    void 
    innerProd (const MV& X, 
	       const MV& Y, 
	       serial_matrix_type& Z) const
    {
      MVT::MvTransMv (SCT::one(), X, Y, Z);
    }

    /// Compute the norm of each column j of X, and return it in
    /// normvec[j].
    void
    norm (const MV& X, 
	  std::vector< MagnitudeType >& normvec) const
    {
      const int ncols = MVT::GetNumberVecs (X);

      // FIXME (mfh 10 Jul 2010, mfh 28 Oct 2010)
      // The current MultiVecTraits, for Tpetra at least, computes
      // column norms one column at a time.  This results in too much
      // communication.  Instead, we compute the whole inner product.
      // This takes more computation, but is the only option
      // available, given the current MVT interface.
      serial_matrix_type XTX (ncols, ncols); 
      innerProd (X, X, XTX);

      // Record the results.  Make sure normvec is long enough.
      if (normvec.size() < ncols)
	normvec.resize (ncols);
      for (int k = 0; k < ncols; ++k)
	normvec[k] = SCTM::squareroot (XTX(k,k));
    }

    /// \brief Compute \f$C := Q^* X\f$ and \f$X := X - Q C\f$.
    ///
    /// \warning This method does not do reorthogonalization.  This is
    /// because the reorthogonalization test requires normalization as
    /// well.
    void 
    project (MV& X, 
	     prev_coeffs_type C,
	     const_prev_mvs_type Q) 
    {
#ifdef BELOS_TEUCHOS_TIME_MONITOR
      // project() is part of orthogonalization, so we time it with
      // total orthogonalization as well as with projection.  If
      // project() is called from projectAndNormalize(), the
      // TimeMonitor won't start timerOrtho_, because it is already
      // running in projectAndNormalize().
      Teuchos::TimeMonitor timerMonitorOrtho(*timerOrtho_);
      Teuchos::TimeMonitor timerMonitorProject(*timerProject_);
#endif
      int nrows_X, ncols_X, num_Q_blocks, ncols_Q_total;
      checkProjectionDims (nrows_X, ncols_X, num_Q_blocks, ncols_Q_total, X, Q);
      // Test for quick exit: any dimension of X is zero, or there are
      // zero Q blocks, or the total number of columns of the Q blocks
      // is zero.
      if (nrows_X == 0 || ncols_X == 0 || num_Q_blocks == 0 || ncols_Q_total == 0)
	return;

      // If we don't have enough C, expanding it creates null references.
      // If we have too many, resizing just throws away the later ones.
      // If we have exactly as many as we have Q, this call has no effect.
      C.resize (num_Q_blocks);
      for (int i = 0; i < num_Q_blocks; ++i) 
	{
	  const int ncols_Q_i = MVT::GetNumberVecs (*Q[i]);
	  // Create a new C[i] if necessary.
	  if (C[i] == Teuchos::null)
	    C[i] = Teuchos::rcp (new serial_matrix_type (ncols_Q_i, ncols_X));
	}
      rawProject (X, Q, C);
    }


    int 
    normalize (MV& X, serial_matrix_ptr B)
    {
#ifdef BELOS_TEUCHOS_TIME_MONITOR
      Teuchos::TimeMonitor timerMonitorNormalize(*timerNormalize_);
      // If normalize() is called internally -- i.e., called from
      // projectAndNormalize() -- the timer will not be started or 
      // stopped, because it is already running.  TimeMonitor handles
      // recursive invocation by doing nothing.
      Teuchos::TimeMonitor timerMonitorOrtho(*timerOrtho_);
#endif

      // Internal data used by this method require a specific MV
      // object for initialization (e.g., to get a Map / communicator,
      // and to initialize scratch space).  Thus, we delay (hence
      // "lazy") initialization until we get an X.
      lazyInit (X);
      
      // MVT returns int for this, even though the local_ordinal_type
      // of the MV may be some other type.
      const int ncols = MVT::GetNumberVecs (X);

      // TSQR's rank-revealing part doesn't work unless B is provided.
      // If B is not provided, allocate a temporary B for use in TSQR.
      // If it is provided, adjust dimensions as necessary.
      if (B == Teuchos::null)
	B = Teuchos::rcp (new serial_matrix_type (ncols, ncols));
      else
	B->shape (ncols, ncols);
      //
      // Compute rank-revealing decomposition (in this case, TSQR of X
      // followed by SVD of the R factor and appropriate updating of
      // the resulting Q factor) of X.  X is modified in place, and Q_
      // contains the resulting explicit Q factor.  Later, we will
      // copy this back into X.  
      //
      // The matrix *B will only be upper triangular if X is of full
      // numerical rank.
      //
      int rank;
      try {
	tsqrAdaptor_->factorExplicit (X, *Q_, *B);
	// This call will only modify *B if *B on input is not of full
	// numerical rank.
	rank = tsqrAdaptor_->revealRank (*Q_, *B, relativeRankTolerance());
      } catch (std::exception& e) {
	throw OrthoError (e.what());
      }

      // Now we should copy Q_ back into X, but don't do it yet: if we
      // want to fill the last ncols-rank columns with random data, we
      // should do so in Q_, because it's fresh in the cache.
      if (randomizeNullSpace_)
	{
	  // If X did not have full (numerical rank), augment the last
	  // ncols-rank columns of X with random data.
	  const int ncolsToFill = ncols - rank;
	  if (ncolsToFill > 0)
	    {
	      // ind: Indices of columns of X to fill with random data.
	      std::vector< int > fillIndices (ncolsToFill);
	      for (int j = 0; j < ncolsToFill; ++j)
		fillIndices[j] = j + rank;

	      mv_ptr Q_null = MVT::CloneViewNonConst (*Q_, fillIndices);
	      MVT::MvRandom (*Q_null);
	      // Done with Q_null; tell Teuchos to deallocate it.
	      Q_null = Teuchos::null;
	    }
	}

      // MultiVecTraits (MVT) doesn't have a "deep copy from one MV
      // into an existing MV in place" method, but it does have two
      // methods which may be used to this effect:
      // 
      // 1. MvAddMv() (compute X := 1*Q_ + 0*X)
      //
      // 2. SetBlock() (Copy from A to mv by setting the index vector
      //    to [0, 1, ..., GetNumberVecs(mv)-1])
      //
      // MVT doesn't state the aliasing rules for MvAddMv(), but I'm
      // guessing it's OK for B and mv to alias one another (that is
      // the way AXPY works in the BLAS).
      MVT::MvAddMv (ScalarType(1), *Q_, ScalarType(0), X, X);

      // Don't deallocate B, even if the user provided Teuchos::null
      // as the B input (or the default B input took over).  The RCP's
      // destructor should work just fine.
      return rank;
    }


    int 
    projectAndNormalize (MV &X,
			 prev_coeffs_type C,
			 serial_matrix_ptr B,
			 const_prev_mvs_type Q)
    {
#ifdef BELOS_TEUCHOS_TIME_MONITOR
      Teuchos::TimeMonitor timerMonitorOrtho(*timerOrtho_);
#endif

      // Fetch dimensions of X and Q.
      int nrows_X, ncols_X, num_Q_blocks, ncols_Q_total;
      checkProjectionDims (nrows_X, ncols_X, num_Q_blocks, ncols_Q_total, X, Q);

      // Test for quick exit: any dimension of X is zero.
      if (nrows_X == 0 || ncols_X == 0)
	return 0;

      // If there are zero Q blocks or zero Q columns, just normalize!
      if (num_Q_blocks == 0 || ncols_Q_total == 0)
	return normalize (X, B);
      
      // If we don't have enough C, expanding it creates null references.
      // If we have too many, resizing just throws away the later ones.
      // If we have exactly as many as we have Q, this call has no effect.
      C.resize (num_Q_blocks);
      prev_coeffs_type newC (num_Q_blocks);
      for (int i = 0; i < num_Q_blocks; ++i) 
	{
	  const int ncols_Q_i = MVT::GetNumberVecs (*Q[i]);
	  // Create a new C[i] if necessary.
	  if (C[i] == Teuchos::null)
	    C[i] = Teuchos::rcp (new serial_matrix_type (ncols_Q_i, ncols_X));
	  // Fill C[i] with zeros.
	  C[i]->putScalar (ScalarType(0));
	  // Create newC[i] as a clone of C[i].  (All that really
	  // matters is that is has the same dimensions and is filled
	  // with zeros; cloning C[i] accomplishes both in this case.)
	  newC[i] = Teuchos::rcp (new serial_matrix_type (*C[i]));
	}

      // Keep track of the column norms of X, both before and after
      // each orthogonalization pass.
      std::vector< MagnitudeType > normsBeforeFirstPass (ncols_X, MagnitudeType(0));
      if (reorthogonalizeBlocks_)
	MVT::MvNorm (X, normsBeforeFirstPass);

      // First BGS pass.  "Modified Gram-Schmidt" version of Block
      // Gram-Schmidt:
      //
      // \li \f$C^{\text{new}}_i := Q_i^* \cdot X\f$
      // \li \f$X := X - Q_i \cdot C^{\text{new}}_i\f$
      rawProject (X, Q, newC);

      // Update the C matrices:
      //
      // \li \f$C_i := C_i + C^{\text{new}}_i\f$
      for (int i = 0; i < num_Q_blocks; ++i)
	*C[i] += *newC[i];

      // Normalize the matrix X.
      if (B == Teuchos::null)
	B = Teuchos::rcp (new serial_matrix_type (ncols_X, ncols_X));
      int rank = normalize (X, B);

      if (reorthogonalizeBlocks_)
	{
	  std::vector< MagnitudeType > normsAfterFirstPass (ncols_X, MagnitudeType(0));
	  std::vector< MagnitudeType > normsAfterSecondPass (ncols_X, MagnitudeType(0));

	  // Compute post-first-pass (pre-normalization) norms, using B.
	  // normalize() doesn't guarantee in general that B is upper
	  // triangular, so we compute norms using the entire column of B.
	  Teuchos::BLAS< int, ScalarType > blas;
	  for (int j = 0; j < ncols_X; ++j)
	    {
	      const ScalarType* const B_j = &(*B)(0,j);
	      // Teuchos::BLAS returns a MagnitudeType result on
	      // ScalarType inputs.
	      normsAfterFirstPass[j] = blas.NRM2 (ncols_X, B_j, 1);
	    }
	  // Test whether any of the norms dropped below the
	  // reorthogonalization threshold.
	  bool reorthog = false;
	  for (int j = 0; j < ncols_X; ++j)
	    // If any column's norm decreased too much, mark this
	    // block for reorthogonalization.  Note that this test
	    // will _not_ activate reorthogonalization if a column's
	    // norm before the first project-and-normalize step was
	    // zero.  It _will_ activate reorthogonalization if the
	    // column's norm before was not zero, but is zero now.
	    if (normsAfterFirstPass[j] < blockReorthogThreshold() * normsBeforeFirstPass[j])
	      {
		reorthog = true; 
		break;
	      }

	  // Perform another Block Gram-Schmidt pass if necessary.  "Twice
	  // is enough" (Kahan's theorem) for a Krylov method, unless
	  // (using Stewart's term) there is an "orthogonalization fault"
	  // (indicated by reorthogFault).
	  bool reorthogFault = false;
	  // Indices of X at which there was an orthogonalization fault.
	  std::vector< int > faultIndices (0);
	  if (reorthog)
	    {
	      // Block Gram-Schmidt (again):
	      //
	      // \li \f$C^{\text{new}} = Q^* X\f$
	      // \li \f$X := X - Q C^{\text{new}}\f$
	      // \li \f$C := C + C^{\text{new}}\f$
	      rawProject (X, Q, newC);
	      for (int i = 0; i < num_Q_blocks; ++i)
		*C[i] += *newC[i];
	      
	      // Normalize the matrix X.
	      serial_matrix_ptr B_new (new serial_matrix_type (ncols_X, ncols_X));
	      rank = normalize (X, B_new);
	      *B += *B_new;
	      
	      // Compute post-second-pass (pre-normalization) norms, using
	      // B.  normalize() doesn't guarantee in general that B is
	      // upper triangular, so we compute norms using the entire
	      // column of B.
	      for (int j = 0; j < ncols_X; ++j)
		{
		  // FIXME Should we use B_new or B here?
		  const ScalarType* const B_j = &(*B)(0,j);
		  // Teuchos::BLAS returns a MagnitudeType result on
		  // ScalarType inputs.
		  normsAfterSecondPass[j] = blas.NRM2 (ncols_X, B_j, 1);
		}
	      // Test whether any of the norms dropped below the
	      // reorthogonalization threshold.  If so, it's an
	      // orthogonalization fault, which requires expensive
	      // recovery.
	      reorthogFault = false;
	      for (int j = 0; j < ncols_X; ++j)
		if (normsAfterSecondPass[j] / normsAfterFirstPass[j] <= blockReorthogThreshold())
		  {
		    reorthogFault = true; 
		    faultIndices.push_back (j);
		  }
	    } // if (reorthog) // reorthogonalization pass

	  if (reorthogFault)
	    {
	      if (throwOnReorthogFault_)
		{
		  using std::endl;
		  typedef std::vector<int>::size_type size_type;
		  std::ostringstream os;
		  
		  os << "Orthogonalization fault at the following column(s) of X:" << endl;
		  os << "Column\tNorm decrease factor" << endl;
		  for (size_type k = 0; k < faultIndices.size(); ++k)
		    {
		      const int index = faultIndices[k];
		      const MagnitudeType decreaseFactor = 
			normsAfterSecondPass[index] / normsAfterFirstPass[index];
		      os << index << "\t" << decreaseFactor << endl;
		    }
		  throw TsqrOrthoFault (os.str());
		}
	      else // Slowly reorthogonalize, one vector at a time, the offending vectors of X.
		{
		  throw std::logic_error ("Slow-and-careful reorthogonalization "
					  "is not yet implemented");

		  // for (int k = 0; k < faultIndices.size(); ++k)
		  // 	{
		  // 	  // Get a view of the current column of X.
		  // 	  std::vector<int> currentIndex (1, faultIndices[k]);
		  // 	  Teuchos::RCP< MV > X_j = MVT::CloneViewNonConst (X, currentIndex);
		  
		  // 	  // Reorthogonalize X_j against all columns of each Q[i].
		  // 	}
		}
	    }
	} // if (reorthogonalizeBlocks_)
      return rank;
    }

    /// \brief Return \f$ \| I - X^* \cdot X \|_F \f$
    ///
    /// Return the Frobenius norm of I - X^* X, which is an absolute
    /// measure of the orthogonality of the columns of X.
    MagnitudeType 
    orthonormError (const MV &X) const
    {
      const ScalarType ONE = SCT::one();
      const int ncols = MVT::GetNumberVecs(X);
      serial_matrix_type XTX (ncols, ncols);
      innerProd (X, X, XTX);
      for (int k = 0; k < ncols; ++k)
	XTX(k,k) -= ONE;
      return XTX.normFrobenius();
    }

    MagnitudeType 
    orthogError (const MV &X1, 
		 const MV &X2) const
    {
      const int ncols_X1 = MVT::GetNumberVecs (X1);
      const int ncols_X2 = MVT::GetNumberVecs (X2);
      serial_matrix_type X1_T_X2 (ncols_X1, ncols_X2);
      innerProd (X1, X2, X1_T_X2);
      return X1_T_X2.normFrobenius();
    }

    /// Relative tolerance for triggering a block reorthogonalization.
    /// If any column norm in a block decreases by this amount, then
    /// we reorthogonalize.
    MagnitudeType blockReorthogThreshold() const { return blockReorthogThreshold_; }
    /// Relative tolerance for determining (via the SVD) whether a
    /// block is of full numerical rank.
    MagnitudeType relativeRankTolerance() const { return relativeRankTolerance_; }

  private:
#ifdef BELOS_TEUCHOS_TIME_MONITOR
    Teuchos::RCP< Teuchos::Time > timerOrtho_, timerProject_, timerNormalize_;
#endif // BELOS_TEUCHOS_TIME_MONITOR

    /// 
    /// Parameters for initializing TSQR
    Teuchos::ParameterList tsqrParams_;
    ///
    /// Label for timers
    std::string label_;
    ///
    /// Interface between Belos and TSQR (the Tall Skinny QR
    /// factorization).
    tsqr_adaptor_ptr tsqrAdaptor_;
    ///
    /// Scratch space for TSQR
    mv_ptr Q_;
    /// 
    // Machine precision for ScalarType
    MagnitudeType eps_;

    /// Whether to fill in null space vectors (after normalization)
    /// with random data.
    bool randomizeNullSpace_;
    ///
    /// Whether to reorthogonalize blocks at all.
    bool reorthogonalizeBlocks_;
    /// Whether to throw an exception when an orthogonalization fault
    /// occurs.  Recovery is possible, but expensive.
    bool throwOnReorthogFault_;
    ///
    /// Relative reorthogonalization threshold in Block Gram-Schmidt.
    MagnitudeType blockReorthogThreshold_;
    ///
    /// Relative tolerance for measuring the numerical rank of a matrix.
    MagnitudeType relativeRankTolerance_;

    static bool 
    getBoolParamWithDefault (const Teuchos::ParameterList& params, 
			     const std::string& key)
    {
      using Teuchos::Exceptions::InvalidParameterName;
      using Teuchos::Exceptions::InvalidParameterType;

      Teuchos::RCP< const Teuchos::ParameterList > defaultParams = 
	getDefaultTsqrParameters< MagnitudeType >();
      try {
	// No validation is necessary; bool is either true or false.
	return params.get< bool >(key);
      } catch (InvalidParameterName&) {
	// Don't try to catch an exception when reading from the
	// default parameters, since they key had better be there.
	return defaultParams->get< bool >(key);
      } catch (InvalidParameterType&) {
	std::ostringstream os;
	os << "The value of parameter \"" << key << "\" to Tsqr(Mat)"
	  "OrthoManager(Impl) must be a boolean truth value.  Here "
	  "is the documentation for that parameter:" << std::endl
	   << defaultParams->getEntry(key).docString();
	throw std::invalid_argument (os.str());
      }
    }

    static MagnitudeType
    getMagParamWithDefault (const Teuchos::ParameterList& params, 
			    const std::string& key)
    {
      using Teuchos::Exceptions::InvalidParameterName;
      using Teuchos::Exceptions::InvalidParameterType;

      Teuchos::RCP< const Teuchos::ParameterList > defaultParams = 
	getDefaultTsqrParameters< MagnitudeType >();
      try {
	const MagnitudeType value = params.get< MagnitudeType >(key);
	if (value >= MagnitudeType(0)) // Validate
	  return value;
	else
	  {
	    std::ostringstream os;
	    os << "You specified " << key << " = " << value
	       << ", but that parameter must be a nonnegative real "
	       << "floating-point value.  Here is the documentation " 
	       << "for that parameter:" << std::endl
	       << defaultParams->getEntry(key).docString();
	    throw std::invalid_argument (os.str());
	  }
      } catch (InvalidParameterName&) { 
	// Don't try to catch an exception when reading from the
	// default parameters, since they key had better be there.
	return defaultParams->get< MagnitudeType >(key);
      } catch (InvalidParameterType&) {
	std::ostringstream os;
	os << "The value of parameter \"" << key << "\" to Tsqr(Mat)OrthoMa"
	  "nager(Impl) must be a nonnegative real floating-point value.  "
	  "Here is the documentation for that parameter:" << std::endl
	   << defaultParams->getEntry(key).docString();
	throw std::invalid_argument (os.str());
      }
    }

    void
    readParams (const Teuchos::ParameterList& params)
    {
      randomizeNullSpace_ = 
	getBoolParamWithDefault (params, "randomizeNullSpace");
      reorthogonalizeBlocks_ = 
	getBoolParamWithDefault (params, "reorthogonalizeBlocks");
      throwOnReorthogFault_ = 
	getBoolParamWithDefault (params, "throwOnReorthogFault");
      blockReorthogThreshold_ = 
	getMagParamWithDefault (params, "blockReorthogThreshold");
      relativeRankTolerance_ = 
	getMagParamWithDefault (params, "relativeRankTolerance");
    }

    /// \brief Initialize the TSQR adaptor and scratch space
    ///
    /// Initialize the TSQR adaptor and scratch space for TSQR.  Both
    /// require a specific MV object, so we have to delay their
    /// initialization until we get an X input (for normalize(), since
    /// only that method uses the TSQR adaptor and the scratch space).
    /// (Hence, "lazy," for delayed initialization.)
    void
    lazyInit (const MV& X)
    {
      // The TSQR adaptor object requires a specific MV object for
      // initialization.  As long as subsequent MV objects use the
      // same communicator (e.g., the same Teuchos::Comm<int>), we
      // don't need to reinitialize the adaptor.
      //
      // FIXME (mfh 15 Jul 2010) If tsqrAdaptor_ has already been
      // initialized, check to make sure that X has the same
      // communicator as the multivector previously used to initialize
      // tsqrAdaptor_.
      if (tsqrAdaptor_ == Teuchos::null)
	tsqrAdaptor_ = Teuchos::rcp (new tsqr_adaptor_type (X, tsqrParams_));

      const int nrows = MVT::GetVecLength (X);
      const int ncols = MVT::GetNumberVecs (X);

      // Q_ is temporary workspace.  It must have the same dimensions
      // as X.  If not, we have to reallocate.  We also have to
      // allocate (not "re-") if we haven't allocated Q_ before.  (We
      // can't allocate Q_ until we have some X, so we need a
      // multivector as the "prototype.")
      if (Q_ == Teuchos::null || 
	  nrows != MVT::GetVecLength(X) || 
	  ncols != MVT::GetVecLength(X))
	// Allocate Q_ to have the same dimensions as X.  Contents of
	// X are not copied.
	Q_ = MVT::Clone (X, ncols);
    }

    void
    checkProjectionDims (int& nrows_X, 
			 int& ncols_X, 
			 int& num_Q_blocks,
			 int& ncols_Q_total,
			 const MV& X, 
			 const_prev_mvs_type Q) const
    {
      // Test for quick exit
      num_Q_blocks = Q.length();
      nrows_X = MVT::GetVecLength (X);
      ncols_X = MVT::GetNumberVecs (X);

      //
      // Make sure that for each i, the dimensions of X and Q[i] are
      // compatible.
      //
      ncols_Q_total = 0; // total over all Q blocks
      for (int i = 0; i < num_Q_blocks; ++i)
	{
	  const int nrows_Q = MVT::GetVecLength (*Q[i]);
	  TEST_FOR_EXCEPTION( (nrows_Q != nrows_X), 
			      std::invalid_argument,
			      "Belos::TsqrOrthoManager::project(): "
			      "Size of X not consistant with size of Q" );
	  ncols_Q_total += MVT::GetNumberVecs (*Q[i]);
	}
    }


    /// Like project(), but does no allocation of blocks of C, and
    /// does no updating of newC (see project() for details).
    void
    rawProject (MV& X, 
		const_prev_mvs_type Q,
		prev_coeffs_type C = Teuchos::tuple (serial_matrix_ptr (Teuchos::null))) const
    {
#ifdef BELOS_TEUCHOS_TIME_MONITOR
      // rawProject() is part of orthogonalization, so we time it with
      // total orthogonalization as well as with projection.  If
      // rawProject() is called from project(), the TimeMonitor won't
      // start timerProject_, because it is already running in
      // project().  Similarly, if rawProject() is called from
      // projectAndNormalize(), the TimeMonitor won't start
      // timerProject_, because it is already running in
      // projectAndNormalize().
      Teuchos::TimeMonitor timerMonitorOrtho(*timerOrtho_);
      Teuchos::TimeMonitor timerMonitorProject(*timerProject_);
#endif

      int nrows_X, ncols_X, num_Q_blocks, ncols_Q_total;
      checkProjectionDims (nrows_X, ncols_X, num_Q_blocks, ncols_Q_total, X, Q);
      // Test for quick exit: any dimension of X is zero, or there are
      // zero Q blocks, or the total number of columns of the Q blocks
      // is zero.
      if (nrows_X == 0 || ncols_X == 0 || num_Q_blocks == 0 || ncols_Q_total == 0)
	return;

      // "Modified Gram-Schmidt" version of Block Gram-Schmidt.
      for (int i = 0; i < num_Q_blocks; ++i)
	{
	  if (C[i] == Teuchos::null)
	    {
	      std::ostringstream os;
	      os << "C[" << i << "] is null" << std::endl;
	      throw std::logic_error (os.str());
	    }
	  innerProd (*Q[i], X, *C[i]);
	  MVT::MvTimesMatAddMv (ScalarType(-1), *Q[i], *C[i], ScalarType(1), X);
	}
    }
  };


  /// \class TsqrOrthoManager
  /// \brief TSQR-based OrthoManager subclass
  ///
  /// This is the actual subclass of OrthoManager, implemented using
  /// TsqrOrthoManagerImpl (TSQR + Block Gram-Schmidt).
  template< class ScalarType, class MV >
  class TsqrOrthoManager : public OrthoManager< ScalarType, MV > {
  public:
    typedef typename Teuchos::ScalarTraits<ScalarType>::magnitudeType magnitude_type;

    TsqrOrthoManager (const Teuchos::ParameterList& tsqrParams, 
		      const std::string& label = "Belos") :
      impl_ (tsqrParams, label)
    {}

    virtual ~TsqrOrthoManager() {}

    virtual void 
    innerProd (const MV &X, 
	       const MV &Y, 
	       Teuchos::SerialDenseMatrix<int,ScalarType>& Z) const
    {
      return impl_.innerProd (X, Y, Z);
    }

    virtual void 
    norm (const MV& X, 
	  std::vector< magnitude_type > &normvec) const
    {
      return impl_.norm (X, normvec);
    }

    virtual void 
    project (MV &X, 
	     Teuchos::Array< Teuchos::RCP< Teuchos::SerialDenseMatrix< int, ScalarType > > > C,
	     Teuchos::Array<Teuchos::RCP<const MV> > Q) const
    {
      return impl_.project (X, C, Q);
    }

    virtual int 
    normalize (MV &X, 
	       Teuchos::RCP< Teuchos::SerialDenseMatrix< int, ScalarType > > B) const
    {
      return impl_.normalize (X, B);
    }

    virtual int 
    projectAndNormalize (MV &X, 
			 Teuchos::Array< Teuchos::RCP< Teuchos::SerialDenseMatrix< int, ScalarType > > > C,
			 Teuchos::RCP< Teuchos::SerialDenseMatrix< int, ScalarType > > B,
			 Teuchos::Array<Teuchos::RCP<const MV> > Q) const
    {
      return impl_.projectAndNormalize (X, C, B, Q);
    }

    virtual typename Teuchos::ScalarTraits< ScalarType >::magnitudeType 
    orthonormError (const MV &X) const
    {
      return impl_.orthonormError (X);
    }

    virtual typename Teuchos::ScalarTraits<ScalarType>::magnitudeType 
    orthogError (const MV &X1, 
		 const MV &X2) const 
    {
      return impl_.orthogError (X1, X2);
    }

    /// Set the label for (the timers for) this orthogonalization
    /// manager, and create new timers if the label has changed.
    ///
    /// \param label [in] New label for timers
    ///
    /// \note Belos::OrthoManager wants this virtual function to be
    ///   implemented; Anasazi::OrthoManager does not.
    void 
    setLabel (const std::string& label) 
    { 
      impl_.setLabel (label);
    }

    const std::string& getLabel() const { return impl_.getLabel(); }

  private:
    /// "Mutable" because it has internal scratch space state.  I know
    /// it's bad, but it's the only way this class can be part of the
    /// OrthoManager hierarchy.
    mutable TsqrOrthoManagerImpl< ScalarType, MV > impl_;

    std::string label_;
  };


  /// \class TsqrMatOrthoManager
  /// \brief MatOrthoManager subclass using TSQR or DGKS
  ///
  /// Subclass of MatOrthoManager.  When getOp() == null (Euclidean
  /// inner product), uses TSQR + Block Gram-Schmidt for
  /// orthogonalization.  When getOp() != null, uses DGKSOrthoManager
  /// for orthogonalization.  Avoids communication only in the TSQR
  /// case.  Initialization of either orthogonalization manager is
  /// "lazy," so you don't have to pay for scratch space if you don't
  /// use it.
  ///
  template< class ScalarType, class MV, class OP >
  class TsqrMatOrthoManager : public MatOrthoManager< ScalarType, MV, OP > {
  private:
    // This will be used to help C++ resolve getOp().  We can't call
    // getOp() directly, because C++ can't figure out that it belongs
    // to the base class MatOrthoManager.  (Remember that at this
    // point, we might not have specialized the specific base class
    // yet; it's just a template at the moment and not a "real
    // class.")
    typedef MatOrthoManager< ScalarType, MV, OP > base_type;

    typedef TsqrOrthoManagerImpl< ScalarType, MV > tsqr_type;
    typedef DGKSOrthoManager< ScalarType, MV, OP > dgks_type;

  public:
    typedef Teuchos::RCP< MV >       mv_ptr;
    typedef Teuchos::RCP< const MV > const_mv_ptr;
    typedef Teuchos::Array< const_mv_ptr >                       const_prev_mvs_type;
    typedef Teuchos::SerialDenseMatrix< int, ScalarType >        serial_matrix_type;
    typedef Teuchos::RCP< serial_matrix_type >                   serial_matrix_ptr;
    typedef Teuchos::Array< Teuchos::RCP< serial_matrix_type > > prev_coeffs_type;
    typedef typename Teuchos::ScalarTraits< ScalarType >::magnitudeType magnitude_type; 

    /// \brief Default constructor (sets Op to Teuchos::null)
    ///
    TsqrMatOrthoManager () :
      MatOrthoManager< ScalarType, MV, OP >(Teuchos::null),
      pTsqr_ (Teuchos::null), // Lazy initialization
      pDgks_ (Teuchos::null)  // Lazy initialization
    {}

    /// Belos::OrthoManager wants this virtual function to be
    /// implemented; Anasazi::OrthoManager does not.
    void setLabel (const std::string& label) {
      label_ = label; 
    }
    const std::string& getLabel() const { return label_; }

    /// \brief Constructor
    ///
    /// \param tsqrParams [in] Parameters used to initialize TSQR 
    ///
    /// \param Op [in] Inner product with respect to which to
    ///   orthogonalize vectors.  If Teuchos::null, use the Euclidean
    ///   inner product.
    TsqrMatOrthoManager (const Teuchos::ParameterList& tsqrParams, 
			 const std::string& label = "Belos",
			 Teuchos::RCP< const OP > Op = Teuchos::null) :
      MatOrthoManager< ScalarType, MV, OP >(Op),
      tsqrParams_ (tsqrParams),
      label_ (label),
      pTsqr_ (Teuchos::null), // Lazy initialization
      pDgks_ (Teuchos::null)  // Lazy initialization
    {}

    /// \brief Destructor
    ///
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
    setOp (Teuchos::RCP< const OP > Op) 
    {
      // We use this notation to help C++ resolve the name.
      // Otherwise, it won't know where to find setOp(), since this is
      // a member function of the base class which does not depend on
      // the template parameters.
      base_type::setOp (Op); // base class gets a copy of the Op too
      ensureDgksInit (); // Make sure the DGKS object has been initialized
      pDgks_->setOp (Op);
    }

    /// \brief Return the inner product operator, if any
    ///
    /// \note We override only to help C++ do name lookup in the other member functions.
    virtual Teuchos::RCP< const OP > getOp () const { return base_type::getOp(); }

    virtual void 
    project (MV &X, 
	     Teuchos::RCP< MV > MX,
	     prev_coeffs_type C,
	     const_prev_mvs_type Q) const
    {
      if (getOp() == Teuchos::null)
	{
	  ensureTsqrInit ();
	  // MX is not modified by MatOrthoManager::project(), so we
	  // don't need to use it here.
	  pTsqr_->project (X, C, Q);
	}
      else
	{
	  ensureDgksInit ();
	  pDgks_->project (X, MX, C, Q);
	}
    }

    virtual void 
    project (MV &X, 
	     prev_coeffs_type C,
	     const_prev_mvs_type Q) const
    {
      project (X, Teuchos::null, C, Q);
    }

    virtual int 
    normalize (MV& X, 
	       Teuchos::RCP< MV > MX,
	       Teuchos::RCP< Teuchos::SerialDenseMatrix< int, ScalarType > > B) const
    {
      if (getOp() == Teuchos::null)
	{
	  ensureTsqrInit ();
	  const int rank = pTsqr_->normalize (X, B);
	  if (MX != Teuchos::null)
	    ; // FIXME: MX should get a copy of X
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
	       Teuchos::RCP< Teuchos::SerialDenseMatrix< int, ScalarType > > B) const 
    {
      return normalize (X, Teuchos::null, B);
    }


    virtual int 
    projectAndNormalize (MV &X, 
			 Teuchos::RCP< MV > MX,
			 Teuchos::Array<Teuchos::RCP<Teuchos::SerialDenseMatrix<int,ScalarType> > > C, 
			 Teuchos::RCP<Teuchos::SerialDenseMatrix<int,ScalarType> > B, 
			 Teuchos::Array<Teuchos::RCP<const MV> > Q ) const
    {
      if (getOp() == Teuchos::null)
	{
	  ensureTsqrInit ();
	  const int rank = pTsqr_->projectAndNormalize (X, C, B, Q); 
	  if (MX != Teuchos::null)
	    ; // FIXME: MX should get a copy of X
	  return rank;
	}
      else
	{
	  ensureDgksInit ();
	  return pDgks_->projectAndNormalize (X, MX, C, B, Q);
	}
    }

    virtual magnitude_type
    orthonormError (const MV &X,
		    Teuchos::RCP< const MV > MX) const
    {
      if (getOp() == Teuchos::null)
	{
	  ensureTsqrInit ();
	  // Ignore MX.
	  return pTsqr_->orthonormError (X);
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
		 Teuchos::RCP< const MV > MX1,
		 const MV &X2) const
    {
      if (getOp() == Teuchos::null)
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
      if (pTsqr_ == Teuchos::null)
	pTsqr_ = Teuchos::rcp (new tsqr_type (tsqrParams_, getLabel()));
    }
    void 
    ensureDgksInit () const
    {
      if (pDgks_ == Teuchos::null)
	pDgks_ = Teuchos::rcp (new dgks_type (getLabel(), getOp()));
    }

    ///
    /// Parameter list for initializing TSQR
    Teuchos::ParameterList tsqrParams_;
    ///
    /// Label for Belos timers
    std::string label_;
    ///
    /// TSQR + BGS orthogonalization manager implementation, used when
    /// getOp() == null (Euclidean inner product).
    mutable Teuchos::RCP< tsqr_type > pTsqr_;
    ///
    /// DGKS orthogonalization manager, used when getOp() != null
    /// (could be a non-Euclidean inner product, but not necessarily).
    mutable Teuchos::RCP< dgks_type > pDgks_;
  };

} // namespace Belos

#endif // __BelosTsqrOrthoManager_hpp

