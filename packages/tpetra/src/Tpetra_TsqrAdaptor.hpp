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

#include <Tsqr_NodeTsqrFactory.hpp>
#include <Tsqr.hpp>
#include <Tsqr_DistTsqrRB.hpp>

#include <Tpetra_MultiVector.hpp>
#include <Teuchos_SerialDenseMatrix.hpp>

#include <stdexcept>

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

namespace Tpetra {

  /// \class TsqrAdaptor
  /// \brief Adaptor from Tpetra::MultiVector to TSQR
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
    TsqrAdaptor (const multivector_type& mv,
		 const Teuchos::ParameterList& plist) :
      pTsqr_ (new tsqr_type (makeNodeTsqr (plist), makeDistTsqr (mv)))
    {}
    
    void
    factorExplicit (multivector_type& A,
		    multivector_type& Q,
		    dense_matrix_type& R)
    {
      // FIXME (mfh 18 Oct 2010) Check Teuchos::Comm<int> objects in A
      // and Q to make sure they are the same communicator as the one
      // we are using in our dist_tsqr_type implementation.

      matview_type R_view (R.numRows(), R.numCols(), R.values(), R.stride());
      pTsqr_->factorExplicit (getNonConstView (A), getNonConstView (Q), R_view);
    }

    void
    revealRank (multivector_type& Q,
		dense_matrix_type& R,
		const magnitude_type& tol)
    {
      // FIXME (mfh 18 Oct 2010) Check Teuchos::Comm<int> object in Q
      // to make sure it is the same communicator as the one we are
      // using in our dist_tsqr_type implementation.

      matview_type Q_view = getNonConstView (Q);
      matview_type R_view (R.numRows(), R.numCols(), R.values(), R.stride());
      pTsqr_->reveal_rank (Q_view.ncols(), Q_view.ncols(), Q.get(), Q.lda(),
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
    getNonConstView (const multivector_type& A)
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
    makeDistTsqr (const multivector_type& mv)
    {
      using Teuchos::Comm;
      using Teuchos::RCP;
      using TSQR::MessengerBase;
      using TSQR::Teuchos::TeuchosMessenger;

      RCP< Comm<int> > pComm = mv.getMap()->getComm();
      RCP< MessengerBase< scalar_type > > pMess (new TeuchosMessenger (pComm));
      RCP< dist_tsqr_type > pDistTsqr (new dist_tsqr_type (pMess));
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

