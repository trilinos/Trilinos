// @HEADER
// ***********************************************************************
// 
//          Tpetra: Templated Linear Algebra Services Package
//                 Copyright (2010) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
// @HEADER

#ifndef __Tpetra_TsqrAdaptor_hpp
#define __Tpetra_TsqrAdaptor_hpp

#include <Kokkos_ConfigDefs.hpp> // HAVE_KOKKOS_TSQR, etc.

#include <Tsqr_NodeTsqrFactory.hpp> // create intranode TSQR object
#include <Tsqr.hpp> // full (internode + intranode) TSQR
#include <Tsqr_DistTsqrRB.hpp> // internode TSQR
// Subclass of TSQR::MessengerBase, implemented using Teuchos
// communicator template helper functions
#include <Tsqr_TeuchosMessenger.hpp> 

#include <Tpetra_MultiVector.hpp>
#include <Teuchos_SerialDenseMatrix.hpp>

#include <stdexcept>

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

namespace Tpetra {

  /// \class TsqrAdaptor
  /// \brief Adaptor from Tpetra::MultiVector to TSQR
  ///
  /// TSQR (Tall Skinny QR factorization) is an orthogonalization
  /// kernel that is as accurate as Householder QR, yet requires only
  /// \f$2 \log P\f$ messages between $P$ MPI processes, independently
  /// of the number of columns in the multivector.  
  ///
  /// TSQR works independently of the particular multivector
  /// implementation, and interfaces to the latter via an adaptor
  /// class.  Tpetra::TsqrAdaptor is the adaptor class for
  /// Tpetra::MultiVector.  It templates on the MultiVector (MV) type
  /// so that it can pick up that class' typedefs.  In particular,
  /// TSQR chooses its intranode implementation based on the Kokkos
  /// Node type of the multivector.
  ///
  template< class MV >
  class TsqrAdaptor {
  public:
    typedef typename MV::scalar_type scalar_type;
    typedef typename MV::local_ordinal_type ordinal_type;
    typedef typename MV::node_type node_type;
    typedef Teuchos::SerialDenseMatrix< ordinal_type, scalar_type > dense_matrix_type;
    typedef typename Teuchos::ScalarTraits< scalar_type >::magnitudeType magnitude_type;

  private:
    typedef TSQR::MatView< ordinal_type, scalar_type > matview_type;
    typedef TSQR::NodeTsqrFactory< node_type, scalar_type, ordinal_type > node_tsqr_factory_type;
    typedef typename node_tsqr_factory_type::node_tsqr_type node_tsqr_type;
    typedef TSQR::DistTsqrRB< ordinal_type, scalar_type > dist_tsqr_type;
    typedef TSQR::Tsqr< ordinal_type, scalar_type, node_tsqr_type, dist_tsqr_type > tsqr_type;

  public:
    /// \brief Constructor
    ///
    /// \param mv [in] Multivector object, used only to access the
    ///   underlying communicator object (in this case,
    ///   Teuchos::Comm<int>, accessed via the Tpetra::Map belonging
    ///   to the multivector).  All multivector objects with which
    ///   this Adaptor works must use the same map and communicator.
    ///
    /// \param plist [in] List of parameters for configuring TSQR.
    ///   The specific parameter keys that are read depend on the
    ///   TSQR implementation.  "cacheBlockSize" (cache block size
    ///   per core, in bytes) tends to be defined for all of the
    ///   non-GPU implementations.  For details, check the specific
    ///   NodeTsqrFactory implementation.
    TsqrAdaptor (const MV& mv,
		 const Teuchos::ParameterList& plist) :
      pTsqr_ (new tsqr_type (makeNodeTsqr (plist), makeDistTsqr (mv)))
    {}

    /// \brief Compute QR factorization [Q,R] = qr(A,0)
    ///
    void
    factorExplicit (MV& A,
		    MV& Q,
		    dense_matrix_type& R)
    {
      // FIXME (mfh 18 Oct 2010) Check Teuchos::Comm<int> objects in A
      // and Q to make sure they are the same communicator as the one
      // we are using in our dist_tsqr_type implementation.

      matview_type R_view (R.numRows(), R.numCols(), R.values(), R.stride());
      pTsqr_->factorExplicit (getNonConstView (A), getNonConstView (Q), R_view);
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
    void
    revealRank (MV& Q,
		dense_matrix_type& R,
		const magnitude_type& tol)
    {
      // FIXME (mfh 18 Oct 2010) Check Teuchos::Comm<int> object in Q
      // to make sure it is the same communicator as the one we are
      // using in our dist_tsqr_type implementation.

      matview_type Q_view = getNonConstView (Q);
      matview_type R_view (R.numRows(), R.numCols(), R.values(), R.stride());
      pTsqr_->revealRank (Q_view.ncols(), Q_view.ncols(), Q.get(), Q.lda(),
			  R.get(), R.lda(), tol);
    }

  private:
    /// Smart pointer to the TSQR implementation object
    ///
    Teuchos::RCP< tsqr_type > pTsqr_;

    /// Return a TSQR::MatView (with raw pointer) from the given
    /// multivector object.  TSQR does not currently support
    /// multivectors with nonconstant stride.
    static matview_type 
    getNonConstView (const MV& A)
    {
      if (! A.isConstantStride())
	{
	  // FIXME (mfh 14 June 2010) Storage of A uses nonconstant
	  // stride internally, but that doesn't necessarily mean we
	  // can't run TSQR.  It depends on what get1dViewNonConst()
	  // returns.  If it's copied and packed into a matrix with
	  // constant stride, then we are free to run TSQR.
	  std::ostringstream os;
	  os << "TSQR does not currently support Tpetra::MultiVector "
	    "inputs that do not have constant stride.";
	  throw std::runtime_error (os.str());
	}
      return matview_type (A.getLocalLength(), 
			   A.getNumVectors, 
			   A.get1dViewNonConst().getRawPtr(), 
			   A.getStride());
    }

    /// Initialize and return internode TSQR implementation
    ///
    static RCP< dist_tsqr_type > 
    makeDistTsqr (const MV& mv)
    {
      using Teuchos::Comm;
      using Teuchos::RCP;
      typedef TSQR::TeuchosMessenger< scalar_type > mess_type;
      typedef TSQR::MessengerBase< scalar_type > base_mess_type;

      RCP< Comm<int> > pComm = mv.getMap()->getComm();
      RCP< mess_type > pMess (new mess_type (pComm));
      RCP< base_mess_type > pMessBase = Teuchos::rcp_implicit_cast (pMess);
      RCP< dist_tsqr_type > pDistTsqr (new dist_tsqr_type (pMessBase));
      return pDistTsqr;
    }

    /// Initialize and return intranode TSQR implementation
    ///
    static RCP< node_tsqr_type >
    makeNodeTsqr (const Teuchos::ParameterList& plist)
    {
      return node_tsqr_factory_type::makeNodeTsqr (plist);
    }
  };

} // namespace Tpetra

#endif // __Tpetra_TsqrAdaptor_hpp

