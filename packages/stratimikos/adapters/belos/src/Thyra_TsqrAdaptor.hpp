// @HEADER
// *****************************************************************************
//         Stratimikos: Thyra-based strategies for linear solvers
//
// Copyright 2006 NTESS and the Stratimikos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef __Thyra_TsqrAdaptor_hpp
#define __Thyra_TsqrAdaptor_hpp

#include "BelosConfigDefs.hpp"

// BelosThyraAdapter.hpp only includes this file if HAVE_BELOS_TSQR is
// defined.  Thus, it's OK to include TSQR header files here.

#include "Thyra_MultiVectorBase.hpp"
#include "Thyra_SpmdVectorSpaceBase.hpp"

#ifdef HAVE_MPI
#  include "Teuchos_DefaultMpiComm.hpp"
#endif // HAVE_MPI
#include "Teuchos_DefaultSerialComm.hpp"
#include "Teuchos_ParameterListAcceptorDefaultBase.hpp"

#include <stdexcept>


namespace Thyra {

  /// \class TsqrAdaptor
  /// \brief Stub adaptor from Thyra::MultiVectorBase to TSQR
  ///
  /// TSQR (Tall Skinny QR factorization) is an orthogonalization
  /// kernel that is as accurate as Householder QR, yet requires only
  /// \f$2 \log P\f$ messages between $P$ MPI processes, independently
  /// of the number of columns in the multivector.
  ///
  /// TSQR works independently of the particular multivector
  /// implementation, and interfaces to the latter via an adaptor
  /// class.  This class is the adaptor class for \c MultiVectorBase.
  /// It templates on the MultiVector (MV) type so that it can pick up
  /// that class' typedefs.  In particular, TSQR chooses its intranode
  /// implementation based on the Kokkos Node type of the multivector.
  ///
  /// \warning This is a stub adaptor that just placates the compiler
  ///   and does nothing.  It's not hard to implement a Thyra adaptor,
  ///   but in order for the adaptor to be efficient, it requires
  ///   special cases for extracting the actual multivector
  ///   implementation (e.g., Epetra_MultiVector or
  ///   Tpetra::MultiVector) out of the Thyra wrapper.
  template<class Scalar>
  class TsqrAdaptor : public Teuchos::ParameterListAcceptorDefaultBase {
  public:
    typedef Thyra::MultiVectorBase<Scalar> MV;
    typedef Scalar scalar_type;
    typedef int ordinal_type; // MultiVectorBase really does use int for this
    typedef Teuchos::SerialDenseMatrix<ordinal_type, scalar_type> dense_matrix_type;
    typedef typename Teuchos::ScalarTraits<scalar_type>::magnitudeType magnitude_type;

    /// \brief Constructor (that accepts a parameter list).
    ///
    /// \param plist [in] List of parameters for configuring TSQR.
    ///   The specific parameter keys that are read depend on the TSQR
    ///   implementation.  For details, call \c getValidParameters()
    ///   and examine the documentation embedded therein.
    TsqrAdaptor (const Teuchos::RCP<Teuchos::ParameterList>& /* plist */)
    {
      throw std::logic_error ("Thyra adaptor for TSQR not implemented");
    }

    //! Constructor (that uses default parameters).
    TsqrAdaptor ()
    {
      throw std::logic_error ("Thyra adaptor for TSQR not implemented");
    }

    Teuchos::RCP<const Teuchos::ParameterList>
    getValidParameters () const
    {
      throw std::logic_error ("Thyra adaptor for TSQR not implemented");
    }

    void
    setParameterList (const Teuchos::RCP<Teuchos::ParameterList>& /* plist */)
    {
      throw std::logic_error ("Thyra adaptor for TSQR not implemented");
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
    factorExplicit (MV& /* A */,
                    MV& /* Q */,
                    dense_matrix_type& /* R */,
                    const bool /* forceNonnegativeDiagonal */ = false)
    {
      throw std::logic_error ("Thyra adaptor for TSQR not implemented");
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
    revealRank (MV& /* Q */,
                dense_matrix_type& /* R */,
                const magnitude_type& /* tol */)
    {
      throw std::logic_error ("Thyra adaptor for TSQR not implemented");
    }

  private:
    /// \brief Attempt to get a communicator out of the given multivector.
    ///
    /// This only works if the multivector's range (VectorSpaceBase)
    /// is actually an SpmdVectorSpaceBase object, and if that
    /// object's Comm is either an MpiComm (in an MPI build) or a
    /// SerialComm (in either an MPI build or a no-MPI build).
    ///
    /// If the attempt does <i>not</i> succeed, this method throws
    /// std::runtime_error.  If it <i>does</i> succeed, it returns the
    /// (suitably wrapped) communicator.
    static Teuchos::RCP<const Teuchos::Comm<int> >
    getComm (const MV& X)
    {
      using Teuchos::RCP;
      using Teuchos::rcp;
      using Teuchos::rcp_dynamic_cast;
      using Teuchos::rcp_implicit_cast;
      typedef Thyra::VectorSpaceBase<Scalar> space_base_type;
      typedef Thyra::SpmdVectorSpaceBase<Scalar> space_type;

      // Thyra stores the communicator in the "vector space," but only
      // if that vector space is an SpmdVectorSpaceBase.
      RCP<const space_base_type> rangeBase = X.range ();
      TEUCHOS_TEST_FOR_EXCEPTION(rangeBase.is_null (), std::runtime_error, "X.range() is null.");
      RCP<const space_type> range = rcp_dynamic_cast<const space_type> (rangeBase);
      TEUCHOS_TEST_FOR_EXCEPTION(range.is_null (), std::runtime_error, "X.range() is not an SpmdVectorSpaceBase.");

      // Thyra annoyingly uses a (possibly) different template
      // parameter for its Teuchos::Comm than everybody else.  The
      // least hackish way to work around this is to convert the Comm
      // to one of two subclasses (MpiComm or SerialComm).  If it's an
      // MpiComm, we can extract the RCP<const OpaqueWrapper<MPI_Comm>
      // > and make a new MpiComm<int> from it.  If it's a SerialComm,
      // just create a new SerialComm<int>.  If it's neither of those,
      // then I have no idea what to do.  Note that MpiComm is only
      // defined if HAVE_MPI is defined.
      RCP<const Teuchos::Comm<Thyra::Ordinal> > thyraComm = range->getComm ();
#ifdef HAVE_MPI
      RCP<const Teuchos::MpiComm<Thyra::Ordinal> > thyraMpiComm =
        rcp_dynamic_cast<const Teuchos::MpiComm<Thyra::Ordinal> > (thyraComm);
      if (thyraMpiComm.is_null ()) {
        RCP<const Teuchos::SerialComm<Thyra::Ordinal> > thyraSerialComm =
          rcp_dynamic_cast<const Teuchos::SerialComm<Thyra::Ordinal> > (thyraComm);
        TEUCHOS_TEST_FOR_EXCEPTION(
          thyraSerialComm.is_null (), std::runtime_error,
          "Thyra's communicator is neither an MpiComm nor a SerialComm.  "
          "Sorry, I have no idea what to do with it in that case.");
        // It's a SerialComm.  Make a SerialComm of our own.
        // SerialComm instances are all the same, so there's no need
        // to keep the original one.
        return rcp_implicit_cast<const Teuchos::Comm<int> > (rcp (new Teuchos::SerialComm<int>));
      }
      else { // Yippie, we have an MpiComm.
        RCP<const Teuchos::OpaqueWrapper<MPI_Comm> > rawMpiComm = thyraMpiComm->getRawMpiComm ();
        // NOTE (mfh 18 Jun 2013) Since the error handler is attached
        // to the MPI_Comm, not to the Teuchos widget, we don't have
        // to set the error handler again on the new MpiComm object.
        return rcp_implicit_cast<const Teuchos::Comm<int> > (rcp (new Teuchos::MpiComm<int> (rawMpiComm)));
      }
#else // NOT HAVE_MPI
      // Either it's a SerialComm or I don't know what to do with it.
      RCP<const Teuchos::SerialComm<Thyra::Ordinal> > thyraSerialComm =
        rcp_dynamic_cast<const Teuchos::SerialComm<Thyra::Ordinal> > (thyraComm);
      TEUCHOS_TEST_FOR_EXCEPTION(
        thyraSerialComm.is_null (), std::runtime_error,
        "Thyra's communicator is not a SerialComm, and MPI is not enabled, so "
        "it can't be an MpiComm either.  That means it must be some other "
        "subclass of Comm, about which I don't know.  "
        "Sorry, I have no idea what to do with it in that case.");
      // It's a SerialComm.  Make a SerialComm of our own.
      // SerialComm instances are all the same, so there's no need
      // to keep the original one.
      return rcp_implicit_cast<const Teuchos::Comm<int> > (rcp (new Teuchos::SerialComm<int>));
#endif // HAVE_MPI
    }

    /// \brief Finish interprocess TSQR initialization.
    ///
    /// Input X is a valid Thyra::MultiVectorBase instance whose
    /// communicator wrapper we will use to prepare TSQR.  It is not
    /// modified.
    ///
    /// \note It's OK to call this method more than once; it is idempotent.
    ///
    /// This method may fail if MV is not the right kind of
    /// multivector, that is, if it does not have a communicator or if
    /// we don't know how to extract a communicator from it.  If it
    /// fails in this way, it will throw std::runtime_error.
    void
    prepareDistTsqr (const MV& /* X */) {}

    /// \brief Finish TSQR initialization.
    ///
    /// The intranode and internode TSQR implementations both have a
    /// two-stage initialization procedure: first, setting parameters
    /// (which may happen at construction), and second, getting
    /// information they need from the multivector input in order to
    /// finish initialization.  For intranode TSQR, this may include
    /// the Kokkos Node instance; for internode TSQR, this includes
    /// the communicator.  The second stage of initialization happens
    /// in this class' computational routines; all of those routines
    /// accept one or more multivector inputs, which this method can
    /// use for finishing initialization.  Thus, users of this class
    /// never need to see the two-stage initialization.
    ///
    /// \param X [in] Multivector object, used only to access the
    ///   underlying communicator object (in this case, the
    ///   Teuchos::Comm<int>) and (possibly) the Kokkos Node instance.
    ///   All multivector objects used with this adapter must have the
    ///   same communicator and Kokkos Node instance (if applicable).
    void
    prepareTsqr (const MV& /* X */) {}
  };

} // namespace Tpetra

#endif // __Thyra_TsqrAdaptor_hpp

