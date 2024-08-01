// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef EPETRA_TSQRADAPTOR_HPP
#define EPETRA_TSQRADAPTOR_HPP

/// \file Epetra_TsqrAdaptor.hpp
/// \brief Epetra_MultiVector to TSQR adaptor
///
/// \note (mfh 27 Oct 2010) This file is in Tpetra (rather than
/// Epetra, where it would seem to belong) as a temporary fix.
/// Otherwise, Epetra would need an optional package dependency on
/// Teuchos and Kokkos, which would break third-party code linking to
/// the Epetra library.  Third-party code should use FIND_PACKAGE on
/// Trilinos to get the correct list of libraries against which to
/// link, but we make this easy temporary fix now so they have time to
/// fix their build systems later.

#include "Tpetra_ConfigDefs.hpp"

#if defined(TPETRA_ENABLE_DEPRECATED_CODE)
#if defined(TPETRA_DEPRECATED_DECLARATIONS)
#warning This file is deprecated due to Epetra removal and will be removed
#endif
#else
#error This file is deprecated due to Epetra removal and will be removed
#endif

#if defined(TPETRA_ENABLE_DEPRECATED_CODE) && defined(HAVE_TPETRA_EPETRA) && defined(HAVE_TPETRA_TSQR)

#include "Tsqr_NodeTsqrFactory.hpp" // create intranode TSQR object
#include "Tsqr.hpp" // full (internode + intranode) TSQR
#include "Tsqr_DistTsqr.hpp" // internode TSQR
#include "Epetra_Comm.h"
// Subclass of TSQR::MessengerBase, implemented using Teuchos
// communicator template helper functions
#include "Epetra_TsqrMessenger.hpp"
#include "Epetra_MultiVector.h"
#include "Teuchos_ParameterListAcceptorDefaultBase.hpp"
#include <stdexcept>

namespace Epetra {

  /// \class TsqrAdaptor
  /// \brief Adaptor from Epetra_MultiVector to TSQR.
  /// \author Mark Hoemmen
  ///
  /// TSQR (Tall Skinny QR factorization) is an orthogonalization
  /// kernel that is as accurate as Householder QR, yet requires only
  /// \f$2 \log P\f$ messages between $P$ MPI processes, independently
  /// of the number of columns in the multivector.
  ///
  /// TSQR works independently of the particular multivector
  /// implementation, and interfaces to the latter via an adaptor
  /// class.  This class is the adaptor class for \c
  /// Epetra_MultiVector.  It templates on the MultiVector (MV) type
  /// so that it can pick up that class' typedefs.  In particular,
  /// TSQR chooses its intranode implementation based on the Kokkos
  /// Node type of the multivector.
  ///
  /// \note Epetra objects live in the global namespace.  TSQR
  ///   requires support for namespaces, so it's acceptable for us to
  ///   create an "Epetra" namespace to contain this adaptor.
  ///
  /// \warning The current implementation of this adaptor requires
  ///   that all Epetra_MultiVector inputs use the same communicator
  ///   object (that is, the same Epetra_Comm) and map.
  TPETRA_DEPRECATED_MSG("epetra removal")
  class TsqrAdaptor : public Teuchos::ParameterListAcceptorDefaultBase {
  public:
    typedef Epetra_MultiVector MV;

    /// \typedef magnitude_type
    ///
    /// Epetra_MultiVector's "Scalar" type is double; it is not a
    /// templated object.  TSQR supports Tpetra as well, in which the
    /// "Scalar" type is a template parameter.  In fact, TSQR supports
    /// complex arithmetic (see the magnitude_type typedef).
    typedef double scalar_type;

    /// \typedef ordinal_type
    ///
    /// In Tpetra terms, this would be the "LocalOrdinal" type.  TSQR
    /// does not depend on the "GlobalOrdinal" type.  Epetra does not
    /// distinguish between the LocalOrdinal and GlobalOrdinal types:
    /// both are int.
    typedef int ordinal_type;

    /// \typedef device_type
    ///
    /// TSQR depends on a Kokkos::Device type.  For Epetra, use a
    /// host-only type.  Typical types are Kokkos::Serial or
    /// Kokkos::OpenMP, depending on build settings.
    using device_type =
      Kokkos::Device<Kokkos::DefaultHostExecutionSpace,
                     Kokkos::HostSpace>;

    /// \typedef dense_matrix_type
    ///
    /// How we pass around small dense matrices that are either local
    /// to each MPI process, or globally replicated.
    ///
    /// \note TSQR lives in the Kokkos package, which requires the
    ///   Teuchos package, so it's acceptable for us to require
    ///   Teuchos components.
    using dense_matrix_type =
      Teuchos::SerialDenseMatrix<ordinal_type, scalar_type>;

    /// \typedef magnitude_type
    ///
    /// Epetra_MultiVector's "Scalar" type is real.  TSQR supports
    /// complex arithmetic as well, in which magnitude_type would
    /// differ from scalar_type.
    using magnitude_type = double;

  private:
    using matview_type = TSQR::MatView<ordinal_type, scalar_type>;
    using node_tsqr_factory_type =
      TSQR::NodeTsqrFactory<scalar_type, ordinal_type, device_type>;
    // Don't need a "typename" here, because there are no template
    // parameters involved in the type definition.
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
      nodeTsqr_ (node_tsqr_factory_type::getNodeTsqr ()),
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
    ///   same communicator and row distribution ("map," in Petra
    ///   terms) as those of the multivector given to this TsqrAdaptor
    ///   instance's constructor.  Otherwise, the result of this
    ///   method is undefined.
    void
    factorExplicit (MV& A,
                    MV& Q,
                    dense_matrix_type& R,
                    const bool forceNonnegativeDiagonal=false)
    {
      prepareTsqr (Q); // Finish initializing TSQR.

      scalar_type* const A_ptr = A.Values ();
      scalar_type* const Q_ptr = Q.Values ();
      scalar_type* const R_ptr = R.values ();
      const ordinal_type numRows = A.MyLength ();
      const ordinal_type numCols = A.NumVectors ();
      const ordinal_type lda = A.Stride ();
      const ordinal_type ldq = Q.Stride ();
      const ordinal_type ldr = R.stride ();

      const bool contiguousCacheBlocks = false;
      tsqr_->factorExplicitRaw (numRows, numCols, A_ptr, lda,
                                Q_ptr, ldq, R_ptr, ldr,
                                contiguousCacheBlocks,
                                forceNonnegativeDiagonal);
    }

    /// \brief Rank-revealing decomposition
    ///
    /// Using the R factor and explicit Q factor from
    /// factorExplicit(), compute the singular value decomposition
    /// (SVD) of R (\f$R = U \Sigma V^*\f$).  If R is full rank (with
    /// respect to the given relative tolerance tol), don't change Q
    /// or R.  Otherwise, compute \f$Q := Q \cdot U\f$ and \f$R :=
    /// \Sigma V^*\f$ in place (the latter may be no longer upper
    /// triangular).
    ///
    /// \param Q [in/out] On input: explicit Q factor computed by
    ///   factorExplicit().  (Must be an orthogonal resp. unitary
    ///   matrix.)  On output: If R is of full numerical rank with
    ///   respect to the tolerance tol, Q is unmodified.  Otherwise, Q
    ///   is updated so that the first rank columns of Q are a basis
    ///   for the column space of A (the original matrix whose QR
    ///   factorization was computed by factorExplicit()).  The
    ///   remaining columns of Q are a basis for the null space of A.
    ///
    /// \param R [in/out] On input: ncols by ncols upper triangular
    ///   matrix with leading dimension ldr >= ncols.  On output: if
    ///   input is full rank, R is unchanged on output.  Otherwise, if
    ///   \f$R = U \Sigma V^*\f$ is the SVD of R, on output R is
    ///   overwritten with $\Sigma \cdot V^*$.  This is also an ncols by
    ///   ncols matrix, but may not necessarily be upper triangular.
    ///
    /// \param tol [in] Relative tolerance for computing the numerical
    ///   rank of the matrix R.
    ///
    /// \return Rank \f$r\f$ of R: \f$ 0 \leq r \leq ncols\f$.
    int
    revealRank (MV& Q,
                dense_matrix_type& R,
                const magnitude_type& tol)
    {
      TEUCHOS_TEST_FOR_EXCEPTION
        (! Q.ConstantStride (), std::invalid_argument, "TsqrAdaptor::"
         "revealRank: Input MultiVector Q must have constant stride.");
      prepareTsqr (Q); // Finish initializing TSQR.
      // FIXME (mfh 25 Oct 2010) Check Epetra_Comm object in Q to make
      // sure it is the same communicator as the one we are using in
      // our dist_tsqr_type implementation.
      return tsqr_->revealRankRaw (Q.MyLength (), Q.NumVectors (),
                                   Q.Values (), Q.Stride (),
                                   R.values (), R.stride (), tol, false);
    }

  private:
    //! The intranode TSQR implementation instance.
    Teuchos::RCP<node_tsqr_type> nodeTsqr_;

    //! The internode TSQR implementation instance.
    Teuchos::RCP<dist_tsqr_type> distTsqr_;

    //! The (full) TSQR implementation instance.
    Teuchos::RCP<tsqr_type> tsqr_;

    //! Default parameter list.  Initialized by \c getValidParameters().
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
    ///   underlying communicator object (in this case, Epetra_Comm).
    ///   All multivector objects used with this Adaptor instance must
    ///   have the same map and communicator.
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
    /// \param mv [in] A multivector, from which to extract the
    ///   Epetra_Comm communicator wrapper to use to initialize TSQR.
    ///
    /// \note It's OK to call this method more than once; it is idempotent.
    void
    prepareDistTsqr (const MV& mv)
    {
      using Teuchos::RCP;
      using Teuchos::rcp;
      using TSQR::Epetra::makeTsqrMessenger;
      typedef TSQR::MessengerBase<scalar_type> base_mess_type;

      // If mv falls out of scope, its Epetra_Comm may become invalid.
      // Thus, we clone the input Epetra_Comm, so that the messenger
      // owns the object.
      RCP<const Epetra_Comm> comm = rcp (mv.Comm().Clone());
      RCP<base_mess_type> messBase = makeTsqrMessenger<scalar_type> (comm);
      distTsqr_->init (messBase);
    }
  };

} // namespace Epetra

#endif // defined(TPETRA_ENABLE_DEPRECATED_CODE) && defined(HAVE_TPETRA_EPETRA) && defined(HAVE_TPETRA_TSQR)

#endif // EPETRA_TSQRADAPTOR_HPP

