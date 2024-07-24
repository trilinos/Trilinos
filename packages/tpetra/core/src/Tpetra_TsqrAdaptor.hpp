// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TPETRA_TSQRADAPTOR_HPP
#define TPETRA_TSQRADAPTOR_HPP

/// \file Tpetra_TsqrAdaptor.hpp
/// \brief Adaptor from Tpetra::MultiVector to TSQR
/// \author Mark Hoemmen

#include "Tpetra_ConfigDefs.hpp"

#ifdef HAVE_TPETRA_TSQR
#  include "Tsqr_NodeTsqrFactory.hpp" // create intranode TSQR object
#  include "Tsqr.hpp" // full (internode + intranode) TSQR
#  include "Tsqr_DistTsqr.hpp" // internode TSQR
// Subclass of TSQR::MessengerBase, implemented using Teuchos
// communicator template helper functions
#  include "Tsqr_TeuchosMessenger.hpp"
#  include "Tpetra_MultiVector.hpp"
#  include "Teuchos_ParameterListAcceptorDefaultBase.hpp"
#  include <stdexcept>

namespace Tpetra {

  /// \class TsqrAdaptor
  /// \brief Adaptor from Tpetra::MultiVector to TSQR
  /// \author Mark Hoemmen
  ///
  /// \tparam MV A specialization of MultiVector.
  ///
  /// TSQR (Tall Skinny QR factorization) is an orthogonalization
  /// kernel that is as accurate as Householder QR, yet requires only
  /// \f$2 \log P\f$ messages between $P$ MPI processes, independently
  /// of the number of columns in the multivector.
  ///
  /// TSQR works independently of the particular multivector
  /// implementation, and interfaces to the latter via an adaptor
  /// class.  This class is the adaptor class for \c MultiVector.  It
  /// templates on the particular specialization of MultiVector, so
  /// that it can pick up the specialization's typedefs.  In
  /// particular, TSQR chooses its intranode implementation based on
  /// the Kokkos Node type of the multivector.
  ///
  /// \warning The current implementation of this adaptor requires
  ///   that all Tpetra::MultiVector inputs use the same Map.
  template<class MV>
  class TsqrAdaptor : public Teuchos::ParameterListAcceptorDefaultBase {
  public:
    using scalar_type = typename MV::scalar_type;
    using ordinal_type = typename MV::local_ordinal_type;
    using dense_matrix_type =
      Teuchos::SerialDenseMatrix<ordinal_type, scalar_type>;
    using magnitude_type =
      typename Teuchos::ScalarTraits<scalar_type>::magnitudeType;

  private:
    using node_tsqr_factory_type =
      TSQR::NodeTsqrFactory<scalar_type, ordinal_type,
                            typename MV::device_type>;
    using node_tsqr_type = TSQR::NodeTsqr<ordinal_type, scalar_type>;
    using dist_tsqr_type = TSQR::DistTsqr<ordinal_type, scalar_type>;
    using tsqr_type = TSQR::Tsqr<ordinal_type, scalar_type>;

    TSQR::MatView<ordinal_type, scalar_type>
    get_mat_view(MV& X)
    {
      TEUCHOS_ASSERT( ! tsqr_.is_null() );
      // FIXME (mfh 18 Oct 2010, 22 Dec 2019) Check Teuchos::Comm<int>
      // object in Q to make sure it is the same communicator as the
      // one we are using in our dist_tsqr_type implementation.

      const ordinal_type lclNumRows(X.getLocalLength());
      const ordinal_type numCols(X.getNumVectors());
      scalar_type* X_ptr = nullptr;
      // LAPACK and BLAS functions require "LDA" >= 1, even if the
      // corresponding matrix dimension is zero.
      ordinal_type X_stride = 1;
      if(tsqr_->wants_device_memory()) {
        auto X_view = X.getLocalViewDevice(Access::ReadWrite);
        X_ptr = reinterpret_cast<scalar_type*>(X_view.data());
        X_stride = static_cast<ordinal_type>(X_view.stride(1));
        if(X_stride == 0) {
          X_stride = ordinal_type(1); // see note above
        }
      }
      else {
        auto X_view = X.getLocalViewHost(Access::ReadWrite);
        X_ptr = reinterpret_cast<scalar_type*>(X_view.data());
        X_stride = static_cast<ordinal_type>(X_view.stride(1));
        if(X_stride == 0) {
          X_stride = ordinal_type(1); // see note above
        }
      }
      using mat_view_type = TSQR::MatView<ordinal_type, scalar_type>;
      return mat_view_type(lclNumRows, numCols, X_ptr, X_stride);
    }

  public:
    /// \brief Constructor that accepts a Teuchos::ParameterList.
    ///
    /// \param plist [in/out] List of parameters for configuring TSQR.
    ///   The specific parameter keys that are read depend on the TSQR
    ///   implementation.  For details, call getValidParameters() and
    ///   examine the documentation embedded therein.
    TsqrAdaptor(const Teuchos::RCP<Teuchos::ParameterList>& plist) :
      nodeTsqr_(node_tsqr_factory_type::getNodeTsqr()),
      distTsqr_(new dist_tsqr_type),
      tsqr_(new tsqr_type(nodeTsqr_, distTsqr_))
    {
      setParameterList(plist);
    }

    //! Constructor(that uses default parameters).
    TsqrAdaptor() :
      nodeTsqr_(node_tsqr_factory_type::getNodeTsqr()),
      distTsqr_(new dist_tsqr_type),
      tsqr_(new tsqr_type(nodeTsqr_, distTsqr_))
    {
      setParameterList(Teuchos::null);
    }

    //! Get all valid parameters (with default values) that TSQR understands.
    Teuchos::RCP<const Teuchos::ParameterList>
    getValidParameters() const
    {
      if(defaultParams_.is_null()) {
        auto params = Teuchos::parameterList("TSQR implementation");
        params->set("NodeTsqr", *(nodeTsqr_->getValidParameters()));
        params->set("DistTsqr", *(distTsqr_->getValidParameters()));
        defaultParams_ = params;
      }
      return defaultParams_;
    }

    /// \brief Set TSQR's parameters.
    ///
    /// \param plist [in/out] List of parameters.
    ///
    /// This method accepts the following sublists:
    ///   - "NodeTsqr": Parameters for the intra-process part of TSQR.
    ///   - "DistTsqr": Parameters for the inter-process part of TSQR.
    ///
    /// Only experts should attempt to set these parameters.  The
    /// default parameters generally perform well.
    ///
    /// The exact set of parameters valid for the "NodeTsqr" sublist
    /// depends on the intra-process TSQR implementation, which in
    /// turn is a function of the Kokkos Node type.  All
    /// implementations accept the "Cache Size Hint" parameter, which
    /// is the cache size in bytes (as a size_t) to use for the
    /// intra-process part of TSQR.  If zero, TSQR will pick a
    /// reasonable default.  The size should correspond to that of the
    /// largest cache that is private to each CPU core, if such a
    /// private cache exists; otherwise, it should correspond to the
    /// amount of shared cache, divided by the number of cores sharing
    /// that cache.  I found through experimentation that TSQR's
    /// performance is not sensitive to this parameter's value, as
    /// long as it is not too large or too small.  The default value
    /// should be fine.
    void
    setParameterList(const Teuchos::RCP<Teuchos::ParameterList>& plist)
    {
      auto params = plist.is_null() ?
        Teuchos::parameterList(*getValidParameters()) : plist;
      using Teuchos::sublist;
      nodeTsqr_->setParameterList(sublist(params, "NodeTsqr"));
      distTsqr_->setParameterList(sublist(params, "DistTsqr"));

      this->setMyParamList(params);
    }

    /// \brief Compute QR factorization [Q,R] = qr(A,0).
    ///
    /// \param A [in/out] On input: the multivector to factor.
    ///   Overwritten with garbage on output.
    ///
    /// \param Q [out] On output: the (explicitly stored) Q factor in
    ///   the QR factorization of the (input) multivector A.
    ///
    /// \param R [out] On output: the R factor in the QR factorization
    ///   of the (input) multivector A.
    ///
    /// \param forceNonnegativeDiagonal [in] If true, then (if
    ///   necessary) do extra work (modifying both the Q and R
    ///   factors) in order to force the R factor to have a
    ///   nonnegative diagonal.
    ///
    /// \warning Currently, this method only works if A and Q have the
    ///   same communicator and row distribution ("Map," in Petra
    ///   terms) as those of the multivector given to this adapter
    ///   instance's constructor.  Otherwise, the result of this
    ///   method is undefined.
    void
    factorExplicit(MV& A,
                   MV& Q,
                   dense_matrix_type& R,
                   const bool forceNonnegativeDiagonal=false)
    {
      TEUCHOS_TEST_FOR_EXCEPTION
        (! A.isConstantStride(), std::invalid_argument, "TsqrAdaptor::"
         "factorExplicit: Input MultiVector A must have constant stride.");
      TEUCHOS_TEST_FOR_EXCEPTION
        (! Q.isConstantStride(), std::invalid_argument, "TsqrAdaptor::"
         "factorExplicit: Input MultiVector Q must have constant stride.");
      prepareTsqr(Q); // Finish initializing TSQR.
      TEUCHOS_ASSERT( ! tsqr_.is_null() );

      auto A_view = get_mat_view(A);
      auto Q_view = get_mat_view(Q);
      constexpr bool contiguousCacheBlocks = false;
      tsqr_->factorExplicitRaw(A_view.extent(0),
                               A_view.extent(1),
                               A_view.data(), A_view.stride(1),
                               Q_view.data(), Q_view.stride(1),
                               R.values(), R.stride(),
                               contiguousCacheBlocks,
                               forceNonnegativeDiagonal);
    }

    /// \brief Rank-revealing decomposition
    ///
    /// Using the R factor and explicit Q factor from
    /// factorExplicit(), compute the singular value decomposition
    /// (SVD) of R: \f$R = U \Sigma V^*\f$.  If R is full rank (with
    /// respect to the given relative tolerance \c tol), do not modify
    /// Q or R.  Otherwise, compute \f$Q := Q \cdot U\f$ and \f$R :=
    /// \Sigma V^*\f$ in place.  If R was modified, then it may not
    /// necessarily be upper triangular on output.
    ///
    /// \param Q [in/out] On input: explicit Q factor computed by
    ///   factorExplicit().  (Must be an orthogonal resp. unitary
    ///   matrix.)  On output: If R is of full numerical rank with
    ///   respect to the tolerance tol, Q is unmodified.  Otherwise, Q
    ///   is updated so that the first \c rank columns of Q are a
    ///   basis for the column space of A (the original matrix whose
    ///   QR factorization was computed by factorExplicit()).  The
    ///   remaining columns of Q are a basis for the null space of A.
    ///
    /// \param R [in/out] On input: N by N upper triangular matrix
    ///   with leading dimension LDR >= N.  On output: if input is
    ///   full rank, R is unchanged on output.  Otherwise, if \f$R = U
    ///   \Sigma V^*\f$ is the SVD of R, on output R is overwritten
    ///   with \f$\Sigma \cdot V^*\f$.  This is also an N by N matrix,
    ///   but it may not necessarily be upper triangular.
    ///
    /// \param tol [in] Relative tolerance for computing the numerical
    ///   rank of the matrix R.
    ///
    /// \return Rank \f$r\f$ of R: \f$ 0 \leq r \leq N\f$.
    int
    revealRank(MV& Q,
               dense_matrix_type& R,
               const magnitude_type& tol)
    {
      TEUCHOS_TEST_FOR_EXCEPTION
        (! Q.isConstantStride(), std::invalid_argument, "TsqrAdaptor::"
         "revealRank: Input MultiVector Q must have constant stride.");
      prepareTsqr(Q); // Finish initializing TSQR.

      auto Q_view = get_mat_view(Q);
      constexpr bool contiguousCacheBlocks = false;
      return tsqr_->revealRankRaw(Q_view.extent(0),
                                  Q_view.extent(1),
                                  Q_view.data(), Q_view.stride(1),
                                  R.values(), R.stride(),
                                  tol, contiguousCacheBlocks);
    }

  private:
    //! The intranode TSQR implementation instance.
    Teuchos::RCP<node_tsqr_type> nodeTsqr_;

    //! The internode TSQR implementation instance.
    Teuchos::RCP<dist_tsqr_type> distTsqr_;

    //! The (full) TSQR implementation instance.
    Teuchos::RCP<tsqr_type> tsqr_;

    //! Default parameter list.  Initialized by getValidParameters().
    mutable Teuchos::RCP<const Teuchos::ParameterList> defaultParams_;

    //! Whether TSQR has been fully initialized.
    bool ready_ = false;

    /// \brief Finish TSQR initialization.
    ///
    /// The intranode and internode TSQR implementations both have a
    /// two-stage initialization procedure: first, setting parameters
    /// (which may happen at construction), and second, getting
    /// information they need from the multivector input in order to
    /// finish initialization.  For intranode TSQR, this includes the
    /// Kokkos Node instance; for internode TSQR, this includes the
    /// communicator.  The second stage of initialization happens in
    /// this class' computational routines; all of those routines
    /// accept one or more multivector inputs, which this method can
    /// use for finishing initialization.  Thus, users of this class
    /// never need to see the two-stage initialization.
    ///
    /// \param mv [in] Multivector object, used only to access the
    ///   underlying communicator object (in this case,
    ///   Teuchos::Comm<int>, accessed via the Tpetra::Map belonging
    ///   to the multivector) and Kokkos Node instance.  All
    ///   multivector objects used with this Adaptor instance must
    ///   have the same map, communicator, and Kokkos Node instance.
    void
    prepareTsqr(const MV& mv)
    {
      if(! ready_) {
        prepareDistTsqr(mv);
        ready_ = true;
      }
    }

    /// \brief Finish interprocess TSQR initialization.
    ///
    /// \param mv [in] A valid Tpetra::MultiVector instance whose
    ///   communicator wrapper we will use to prepare TSQR.
    ///
    /// \note It's OK to call this method more than once; it is idempotent.
    void
    prepareDistTsqr(const MV& mv)
    {
      using Teuchos::RCP;
      using Teuchos::rcp_implicit_cast;
      using mess_type = TSQR::TeuchosMessenger<scalar_type>;
      using base_mess_type = TSQR::MessengerBase<scalar_type>;

      auto comm = mv.getMap()->getComm();
      RCP<mess_type> mess(new mess_type(comm));
      auto messBase = rcp_implicit_cast<base_mess_type>(mess);
      distTsqr_->init(messBase);
    }
  };

} // namespace Tpetra

#endif // HAVE_TPETRA_TSQR

#endif // TPETRA_TSQRADAPTOR_HPP
