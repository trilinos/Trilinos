//@HEADER
// ************************************************************************
// 
//               Epetra: Linear Algebra Services Package 
//                 Copyright 2010 Sandia Corporation
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

#ifndef __Epetra_TsqrAdaptor_hpp
#define __Epetra_TsqrAdaptor_hpp

///
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
///

#include <Tpetra_ConfigDefs.hpp>

#if defined(HAVE_TPETRA_EPETRA) && defined(HAVE_TPETRA_TSQR)

#include <Kokkos_DefaultNode.hpp> // Include minimal Kokkos Node types
#include <Tsqr_NodeTsqrFactory.hpp> // create intranode TSQR object
#include <Tsqr.hpp> // full (internode + intranode) TSQR
#include <Tsqr_DistTsqr.hpp> // internode TSQR

#include <Epetra_Comm.h>
// Subclass of TSQR::MessengerBase, implemented using Teuchos
// communicator template helper functions
#include <Epetra_TsqrMessenger.hpp>
#include <Epetra_MultiVector.h>

#include <stdexcept>

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

namespace Epetra {

  /// \class TsqrAdaptor
  /// \brief Adaptor from Epetra_MultiVector to TSQR
  ///
  /// TSQR (Tall Skinny QR factorization) is an orthogonalization
  /// kernel that is as accurate as Householder QR, yet requires only
  /// \f$2 \log P\f$ messages between $P$ MPI processes, independently
  /// of the number of columns in the multivector.  
  ///
  /// TSQR works independently of the particular multivector
  /// implementation, and interfaces to the latter via an adaptor
  /// class.  Epetra::TsqrAdaptor is the adaptor class for
  /// Epetra_MultiVector.  It templates on the MultiVector (MV) type
  /// so that it can pick up that class' typedefs.  In particular,
  /// TSQR chooses its intranode implementation based on the Kokkos
  /// Node type of the multivector.
  ///
  /// \note Epetra objects live in the global namespace.  TSQR
  ///   requires support for namespaces, so it's acceptable for us to
  ///   create an "Epetra" namespace to contain this adaptor.
  ///
  class TsqrAdaptor {
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

    /// \typedef node_type
    ///
    /// TSQR depends on a Kokkos Node type.  We could use the
    /// Kokkos::DefaultNode::DefaultNodeType typedef, but (a) we want
    /// to ensure the expected "sequential within one MPI process"
    /// semantics of Epetra, and (b) we don't have a good
    /// platform-independent automatic mechanism for determining how
    /// many threads each MPI process should use, when running
    /// multiple MPI processes on a node.  Thus, we use
    /// Kokkos::SerialNode.
    typedef Kokkos::SerialNode node_type;

    /// \typedef dense_matrix_type 
    ///
    /// How we pass around small dense matrices that are either local
    /// to each MPI process, or globally replicated.  
    ///
    /// \note TSQR lives in the Kokkos package, which requires the
    ///   Teuchos package, so it's acceptable for us to require
    ///   Teuchos components.
    typedef Teuchos::SerialDenseMatrix< ordinal_type, scalar_type > dense_matrix_type;

    /// \typedef magnitude_type
    ///
    /// Epetra_MultiVector's "Scalar" type is real.  TSQR supports
    /// complex arithmetic as well, in which magnitude_type would
    /// differ from scalar_type.
    typedef double magnitude_type;

  private:
    typedef TSQR::MatView< ordinal_type, scalar_type > matview_type;
    typedef TSQR::NodeTsqrFactory< node_type, scalar_type, ordinal_type > node_tsqr_factory_type;
    // Don't need a "typename" here, because there are no template
    // parameters involved in the type definition.
    typedef node_tsqr_factory_type::node_tsqr_type node_tsqr_type;
    typedef TSQR::DistTsqr< ordinal_type, scalar_type > dist_tsqr_type;
    typedef TSQR::Tsqr< ordinal_type, scalar_type, node_tsqr_type > tsqr_type;

  public:
    /// \brief Constructor
    ///
    /// \param mv [in] Multivector object, used only to access the
    ///   underlying communicator object (in this case, Epetra_Comm,
    ///   accessed directly from the Epetra_MultiVector input).  All
    ///   multivector objects with which this Adaptor works must use
    ///   the same map and communicator.
    ///
    /// \param plist [in] List of parameters for configuring TSQR.
    ///   The specific parameter keys that are read depend on the
    ///   TSQR implementation.  "cacheBlockSize" (cache block size
    ///   per core, in bytes) tends to be defined for all of the
    ///   non-GPU implementations.  For details, check the specific
    ///   NodeTsqrFactory implementation.
    ///
    /// \warning The current implementation of this adaptor requires
    ///   that all Epetra_MultiVector inputs use the same communicator
    ///   object (that is, the same Epetra_Comm).  This will be fixed
    ///   in the future -- it's not hard to fix.
    ///
    TsqrAdaptor (const MV& mv,
		 const Teuchos::ParameterList& plist) :
      pTsqr_ (new tsqr_type (makeNodeTsqr (plist), makeDistTsqr (mv)))
    {
      Teuchos::ParameterList emptyParams;
      pNode_ = Teuchos::rcp (new Kokkos::SerialNode (emptyParams));
    }

    /// \brief Compute QR factorization [Q,R] = qr(A,0)
    ///
    void
    factorExplicit (MV& A,
		    MV& Q,
		    dense_matrix_type& R)
    {
      typedef Kokkos::MultiVector< scalar_type, node_type > KMV;

      // FIXME (mfh 25 Oct 2010) Check Epetra_Comm objects in A and Q
      // to make sure they are the same communicator as the one we are
      // using in our dist_tsqr_type implementation.
      KMV A_view = getNonConstView (A);
      KMV Q_view = getNonConstView (Q);
      pTsqr_->factorExplicit (A_view, Q_view, R, false);
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
    ///
    int
    revealRank (MV& Q,
		dense_matrix_type& R,
		const magnitude_type& tol)
    {
      typedef Kokkos::MultiVector< scalar_type, node_type > KMV;

      // FIXME (mfh 25 Oct 2010) Check Epetra_Comm object in Q to make
      // sure it is the same communicator as the one we are using in
      // our dist_tsqr_type implementation.

      KMV Q_view = getNonConstView (Q);
      return pTsqr_->revealRank (Q_view, R, tol, false);
    }

  private:
    /// Smart pointer to the TSQR implementation object
    ///
    Teuchos::RCP< tsqr_type > pTsqr_;

    /// Smart pointer to a Kokkos::SerialNode.  It really doesn't hold
    /// or do anything, but we keep it around so that we don't have to
    /// call the SerialNode constructor more than once.
    Teuchos::RCP< Kokkos::SerialNode > pNode_;

    /// \brief Return a nonconstant view of the input MultiVector
    ///
    /// TSQR represents the local (to each MPI process) part of a
    /// MultiVector as a Kokkos::MultiVector (KMV), which gives a
    /// nonconstant view of the original MultiVector's data.  This
    /// class method tells TSQR how to get the KMV from the input.
    ///
    /// \warning TSQR does not currently support multivectors with
    ///   nonconstant stride.  This method will raise an exception
    ///   if A has nonconstant stride.
    ///
    /// \param MV [in] MultiVector object, for which to return a
    ///   Kokkos::MultiVector nonconstant view.
    ///
    /// \return Nonconstant view of the input MultiVector object. 
    ///
    /// \note The view never escapes the scope of the creating 
    Kokkos::MultiVector< scalar_type, node_type >
    getNonConstView (MV& A)
    {
      if (! A.ConstantStride())
	{
	  // FIXME (mfh 25 Oct 2010) Storage of A uses nonconstant
	  // stride internally, but that doesn't necessarily mean we
	  // can't run TSQR.  We would have to copy and pack into a
	  // matrix with constant stride, and then unpack on exit.
	  std::ostringstream os;
	  os << "TSQR does not currently support Epetra_MultiVector "
	    "inputs that do not have constant stride.";
	  throw std::runtime_error (os.str());
	}

      const int numRows = A.MyLength();
      const int numCols = A.NumVectors();
      const int stride  = A.Stride();
      // A_ptr does _not_ own the data.  TSQR only operates within the
      // scope of the multivector objects on which it operates, so it
      // doesn't need ownership of the data.
      Teuchos::ArrayRCP< double > A_ptr (A.Values(), 0, numRows*stride, false);

      typedef Kokkos::MultiVector< scalar_type, node_type > KMV;
      KMV A_kmv (pNode_);
      A_kmv.initializeValues (numRows, numCols, A_ptr, stride);
      return A_kmv;
    }

    /// Initialize and return internode TSQR implementation
    ///
    static Teuchos::RCP< dist_tsqr_type > 
    makeDistTsqr (const MV& mv)
    {
      using Teuchos::RCP;
      typedef TSQR::MessengerBase< scalar_type > base_mess_type;
      using TSQR::Epetra::makeTsqrMessenger;

      // pComm is a nonowning pointer to the reference.
      RCP< const Epetra_Comm > pComm = Teuchos::rcpFromRef (mv.Comm());
      RCP< base_mess_type > pMessBase = makeTsqrMessenger< scalar_type >(pComm);
      RCP< dist_tsqr_type > pDistTsqr (new dist_tsqr_type (pMessBase));
      return pDistTsqr;
    }

    /// Initialize and return intranode TSQR implementation
    ///
    static Teuchos::RCP< node_tsqr_type >
    makeNodeTsqr (const Teuchos::ParameterList& plist)
    {
      return node_tsqr_factory_type::makeNodeTsqr (plist);
    }
  };

} // namespace Epetra

#endif // defined(HAVE_TPETRA_EPETRA) && defined(HAVE_TPETRA_TSQR)

#endif // __Epetra_TsqrAdaptor_hpp

