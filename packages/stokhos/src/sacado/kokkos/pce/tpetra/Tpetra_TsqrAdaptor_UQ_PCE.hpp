// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TPETRA_TSQR_ADAPTOR_UQ_PCE_HPP
#define TPETRA_TSQR_ADAPTOR_UQ_PCE_HPP

#include <Tpetra_ConfigDefs.hpp> // HAVE_TPETRA_TSQR, etc.

#ifdef HAVE_TPETRA_TSQR

#include "Stokhos_Sacado_Kokkos_UQ_PCE.hpp"

#  include "Tsqr_NodeTsqrFactory.hpp" // create intranode TSQR object
#  include "Tsqr.hpp" // full (internode + intranode) TSQR
#  include "Tsqr_DistTsqr.hpp" // internode TSQR
// Subclass of TSQR::MessengerBase, implemented using Teuchos
// communicator template helper functions
#  include "Tsqr_TeuchosMessenger.hpp"
#  include "Tpetra_MultiVector.hpp"
#  include "Teuchos_ParameterListAcceptorDefaultBase.hpp"
#  include <stdexcept>

// Base TsqrAdator template we will specialize
#  include "Tpetra_TsqrAdaptor.hpp"

namespace Tpetra {

  /// \class TsqrAdaptor
  /// \brief Adaptor from Tpetra::MultiVector to TSQR for UQ::PCE scalar type
  /// \author Eric Phipps
  ///
  /// This specialization works be extracting the underlying array within the
  /// multivector and converting to a standard scalar type.
  template <class Storage, class LO, class GO, class Node>
  class TsqrAdaptor< Tpetra::MultiVector< Sacado::UQ::PCE<Storage>,
                                          LO, GO, Node > > :
    public Teuchos::ParameterListAcceptorDefaultBase {
  public:
    typedef Tpetra::MultiVector< Sacado::UQ::PCE<Storage>, LO, GO, Node > MV;
    typedef typename MV::scalar_type mp_scalar_type;

    // For Sacado::UQ::PCE< Storage<Ordinal,Scalar,Device> > this is Scalar
    typedef typename mp_scalar_type::scalar_type scalar_type;
    typedef typename mp_scalar_type::ordinal_type mp_ordinal_type;
    typedef typename MV::local_ordinal_type ordinal_type;
    typedef Teuchos::SerialDenseMatrix<ordinal_type, scalar_type> dense_matrix_type;
    typedef typename Teuchos::ScalarTraits<scalar_type>::magnitudeType magnitude_type;

  private:
    using node_tsqr_factory_type =
      TSQR::NodeTsqrFactory<scalar_type, ordinal_type,
                            typename MV::device_type>;
    using node_tsqr_type = TSQR::NodeTsqr<ordinal_type, scalar_type>;
    using dist_tsqr_type = TSQR::DistTsqr<ordinal_type, scalar_type>;
    using tsqr_type = TSQR::Tsqr<ordinal_type, scalar_type>;

  public:
    /// \brief Constructor (that accepts a parameter list).
    ///
    /// \param plist [in/out] List of parameters for configuring TSQR.
    ///   The specific parameter keys that are read depend on the TSQR
    ///   implementation.  For details, call \c getValidParameters()
    ///   and examine the documentation embedded therein.
    TsqrAdaptor (const Teuchos::RCP<Teuchos::ParameterList>& plist) :
      nodeTsqr_ (node_tsqr_factory_type::getNodeTsqr ()),
      distTsqr_ (new dist_tsqr_type),
      tsqr_ (new tsqr_type (nodeTsqr_, distTsqr_)),
      ready_ (false)
    {
      setParameterList (plist);
    }

    //! Constructor (that uses default parameters).
    TsqrAdaptor () :
      nodeTsqr_ (new node_tsqr_factory_type::getNodeTsqr ()),
      distTsqr_ (new dist_tsqr_type),
      tsqr_ (new tsqr_type (nodeTsqr_, distTsqr_)),
      ready_ (false)
    {
      setParameterList (Teuchos::null);
    }

    Teuchos::RCP<const Teuchos::ParameterList>
    getValidParameters () const
    {
      using Teuchos::RCP;
      using Teuchos::rcp;
      using Teuchos::ParameterList;
      using Teuchos::parameterList;

      if (defaultParams_.is_null()) {
        RCP<ParameterList> params = parameterList ("TSQR implementation");
        params->set ("NodeTsqr", *(nodeTsqr_->getValidParameters ()));
        params->set ("DistTsqr", *(distTsqr_->getValidParameters ()));
        defaultParams_ = params;
      }
      return defaultParams_;
    }

    void
    setParameterList (const Teuchos::RCP<Teuchos::ParameterList>& plist)
    {
      using Teuchos::ParameterList;
      using Teuchos::parameterList;
      using Teuchos::RCP;
      using Teuchos::sublist;

      RCP<ParameterList> params = plist.is_null() ?
        parameterList (*getValidParameters ()) : plist;
      nodeTsqr_->setParameterList (sublist (params, "NodeTsqr"));
      distTsqr_->setParameterList (sublist (params, "DistTsqr"));

      this->setMyParamList (params);
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
    factorExplicit (MV& A,
                    MV& Q,
                    dense_matrix_type& R,
                    const bool forceNonnegativeDiagonal=false)
    {
      prepareTsqr (Q); // Finish initializing TSQR.

      ordinal_type numRows;
      ordinal_type numCols;
      ordinal_type LDA;
      ordinal_type LDQ;
      scalar_type* A_ptr;
      scalar_type* Q_ptr;

      getNonConstView (numRows, numCols, A_ptr, LDA, A);
      getNonConstView (numRows, numCols, Q_ptr, LDQ, Q);
      const bool contiguousCacheBlocks = false;
      tsqr_->factorExplicitRaw (numRows, numCols, A_ptr, LDA,
                                Q_ptr, LDQ, R.values (), R.stride (),
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
    revealRank (MV& Q,
                dense_matrix_type& R,
                const magnitude_type& tol)
    {
      prepareTsqr (Q); // Finish initializing TSQR.

      // FIXME (mfh 18 Oct 2010) Check Teuchos::Comm<int> object in Q
      // to make sure it is the same communicator as the one we are
      // using in our dist_tsqr_type implementation.

      ordinal_type numRows;
      ordinal_type numCols;
      scalar_type* Q_ptr;
      ordinal_type LDQ;
      getNonConstView (numRows, numCols, Q_ptr, LDQ, Q);
      const bool contiguousCacheBlocks = false;
      return tsqr_->revealRankRaw (numRows, numCols, Q_ptr, LDQ,
                                   R.values (), R.stride (), tol,
                                   contiguousCacheBlocks);
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
    bool ready_;

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
    prepareTsqr (const MV& mv)
    {
      if (! ready_) {
        prepareDistTsqr (mv);
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
    prepareDistTsqr (const MV& mv)
    {
      using Teuchos::RCP;
      using Teuchos::rcp_implicit_cast;
      typedef TSQR::TeuchosMessenger<scalar_type> mess_type;
      typedef TSQR::MessengerBase<scalar_type> base_mess_type;

      RCP<const Teuchos::Comm<int> > comm = mv.getMap()->getComm();
      RCP<mess_type> mess (new mess_type (comm));
      RCP<base_mess_type> messBase = rcp_implicit_cast<base_mess_type> (mess);
      distTsqr_->init (messBase);
    }

    /// \brief Extract a nonpersistent view of A's data as a
    ///   scalar_type matrix, stored as a flat column-major array.
    ///
    /// \warning TSQR does not currently support multivectors with
    ///   nonconstant stride.  If A has nonconstant stride, this
    ///   method will throw an exception.
    static void
    getNonConstView (ordinal_type& numRows,
                     ordinal_type& numCols,
                     scalar_type*& A_ptr,
                     ordinal_type& LDA,
                     const MV& A)
    {
      // FIXME (mfh 25 Oct 2010) We should be able to run TSQR even if
      // storage of A uses nonconstant stride internally.  We would
      // have to copy and pack into a matrix with constant stride, and
      // then unpack on exit.  For now we choose just to raise an
      // exception.
      TEUCHOS_TEST_FOR_EXCEPTION
        (! A.isConstantStride(), std::invalid_argument,
         "TSQR does not currently support Tpetra::MultiVector "
         "inputs that do not have constant stride.");

      // FIXME (mfh 16 Jan 2016) When I got here, I found strides[0]
      // instead of strides[1] for the stride.  I don't think this is
      // right.  However, I don't know about these Stokhos scalar
      // types so I'll just do what was here.
      //
      // STOKHOS' TYPES ARE NOT TESTED WITH TSQR REGULARLY SO IT IS
      // POSSIBLE THAT THE ORIGINAL CODE WAS WRONG.

      typedef typename MV::dual_view_type view_type;
      typedef typename view_type::t_dev::array_type flat_array_type;

      // Reinterpret the data as a longer array of the base scalar
      // type.  TSQR currently forbids MultiVector input with
      // nonconstant stride, so we need not worry about that here.

      flat_array_type flat_mv = A.getLocalViewDevice();

      numRows = static_cast<ordinal_type> (flat_mv.extent(0));
      numCols = static_cast<ordinal_type> (flat_mv.extent(1));
      A_ptr = flat_mv.data ();

      ordinal_type strides[2];
      flat_mv.stride (strides);
      LDA = strides[0];
    }
  };
} // namespace Tpetra

#endif // HAVE_TPETRA_TSQR

#endif // TPETRA_TSQR_ADAPTOR_UQ_PCE_HPP
