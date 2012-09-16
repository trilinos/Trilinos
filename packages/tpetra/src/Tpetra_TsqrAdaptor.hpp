// @HEADER
// ***********************************************************************
// 
//          Tpetra: Templated Linear Algebra Services Package
//                 Copyright (2008) Sandia Corporation
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
// @HEADER

#ifndef __Tpetra_TsqrAdaptor_hpp
#define __Tpetra_TsqrAdaptor_hpp

#include <Tpetra_ConfigDefs.hpp> // HAVE_TPETRA_TSQR, etc.

#ifdef HAVE_TPETRA_TSQR
#  include <Tsqr_NodeTsqrFactory.hpp> // create intranode TSQR object
#  include <Tsqr.hpp> // full (internode + intranode) TSQR
#  include <Tsqr_DistTsqr.hpp> // internode TSQR
// Subclass of TSQR::MessengerBase, implemented using Teuchos
// communicator template helper functions
#  include <Tsqr_TeuchosMessenger.hpp> 
#  include <Tpetra_MultiVector.hpp>
#  include <Teuchos_ParameterListAcceptorDefaultBase.hpp>
#  include <stdexcept>


namespace Tpetra {

  /// \class TsqrAdaptor
  /// \brief Adaptor from Tpetra::MultiVector to TSQR
  /// \author Mark Hoemmen
  ///
  /// \tparam MV A specialization of \c Tpetra::MultiVector.
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
  ///   that all Tpetra::MultiVector inputs use the same communicator
  ///   object (that is, the same Epetra_Comm) and map.
  template<class MV>
  class TsqrAdaptor : public Teuchos::ParameterListAcceptorDefaultBase {
  public:
    typedef typename MV::scalar_type scalar_type;
    typedef typename MV::local_ordinal_type ordinal_type;
    typedef typename MV::node_type node_type;
    typedef Teuchos::SerialDenseMatrix<ordinal_type, scalar_type> dense_matrix_type;
    typedef typename Teuchos::ScalarTraits<scalar_type>::magnitudeType magnitude_type;

  private:
    typedef TSQR::MatView<ordinal_type, scalar_type> matview_type;
    typedef TSQR::NodeTsqrFactory<node_type, scalar_type, ordinal_type> node_tsqr_factory_type;
    typedef typename node_tsqr_factory_type::node_tsqr_type node_tsqr_type;
    typedef TSQR::DistTsqr<ordinal_type, scalar_type> dist_tsqr_type;
    typedef TSQR::Tsqr<ordinal_type, scalar_type, node_tsqr_type> tsqr_type;

  public:
    /// \brief Constructor (that accepts a parameter list).
    ///
    /// \param plist [in/out] List of parameters for configuring TSQR.
    ///   The specific parameter keys that are read depend on the TSQR
    ///   implementation.  For details, call \c getValidParameters()
    ///   and examine the documentation embedded therein.
    TsqrAdaptor (const Teuchos::RCP<Teuchos::ParameterList>& plist) :
      nodeTsqr_ (new node_tsqr_type),
      distTsqr_ (new dist_tsqr_type),
      tsqr_ (new tsqr_type (nodeTsqr_, distTsqr_)),
      ready_ (false)
    {
      setParameterList (plist);
    }

    //! Constructor (that uses default parameters).
    TsqrAdaptor () : 
      nodeTsqr_ (new node_tsqr_type),
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
      typedef Kokkos::MultiVector<scalar_type, node_type> KMV;

      prepareTsqr (Q); // Finish initializing TSQR.
      KMV A_view = getNonConstView (A);
      KMV Q_view = getNonConstView (Q);
      tsqr_->factorExplicit (A_view, Q_view, R, false, 
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
    ///   overwritten with \f$\Sigma \cdot V^*\f$.  This is also an
    ///   ncols by ncols matrix, but may not necessarily be upper
    ///   triangular.
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
      typedef Kokkos::MultiVector<scalar_type, node_type> KMV;

      prepareTsqr (Q); // Finish initializing TSQR.      

      // FIXME (mfh 18 Oct 2010) Check Teuchos::Comm<int> object in Q
      // to make sure it is the same communicator as the one we are
      // using in our dist_tsqr_type implementation.
      KMV Q_view = getNonConstView (Q);
      return tsqr_->revealRank (Q_view, R, tol, false);
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
        prepareNodeTsqr (mv);
        ready_ = true;
      }
    }

    /// \brief Finish intranode TSQR initialization.
    ///
    /// \note It's OK to call this method more than once; it is idempotent.
    void
    prepareNodeTsqr (const MV& mv)
    {
      node_tsqr_factory_type::prepareNodeTsqr (nodeTsqr_, mv.getMap()->getNode());
    }

    /// \brief Finish internode TSQR initialization.
    ///
    /// \param mv [in] A valid Tpetra_MultiVector instance whose
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

    /// \brief Extract A's underlying Kokkos::MultiVector instance.
    ///
    /// TSQR represents the local (to each MPI process) part of a
    /// multivector as a Kokkos::MultiVector (KMV), which gives a
    /// nonconstant view of the original multivector's data.  This
    /// class method tells TSQR how to get the KMV from the input
    /// multivector.  The KMV is not a persistent view of the data;
    /// its scope is contained within the scope of the multivector.
    ///
    /// \warning TSQR does not currently support multivectors with
    ///   nonconstant stride.  If A has nonconstant stride, this
    ///   method will throw an exception.
    static Kokkos::MultiVector<scalar_type, node_type>
    getNonConstView (MV& A)
    {
      // FIXME (mfh 25 Oct 2010) We should be able to run TSQR even if
      // storage of A uses nonconstant stride internally.  We would
      // have to copy and pack into a matrix with constant stride, and
      // then unpack on exit.  For now we choose just to raise an
      // exception.
      TEUCHOS_TEST_FOR_EXCEPTION(! A.isConstantStride(), std::invalid_argument,
                                 "TSQR does not currently support Tpetra::MultiVector "
                                 "inputs that do not have constant stride.");
      return A.getLocalMVNonConst();
    }
  };

} // namespace Tpetra

#endif // HAVE_TPETRA_TSQR

#endif // __Tpetra_TsqrAdaptor_hpp

