#ifndef __TSQR_Trilinos_TsqrAdaptor_Tpetra_MultiVector_SerialNode_hpp
#define __TSQR_Trilinos_TsqrAdaptor_Tpetra_MultiVector_SerialNode_hpp

/// \file TsqrAdaptor_Tpetra_MultiVector.hpp
///
/// \warning Users should _not_ include this file directly.  Include
///   "TsqrAdaptor.hpp" instead.  If HAVE_ANASAZI_TPETRA is defined,
///   then this file will be included automatically.  If for some
///   reason you need to include this file directly, be sure to
///   include "TsqrAdaptor.hpp" first.

#include "TsqrTypeAdaptor_Tpetra_MultiVector_SerialNode.hpp"
#ifdef HAVE_KOKKOS_TBB
#  include "TsqrTypeAdaptor_Tpetra_MultiVector_TBBNode.hpp"
#endif // HAVE_KOKKOS_TBB

#include <stdexcept>
#include <sstream>

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

namespace TSQR {
  namespace Trilinos {

    template< class S, class LO, class GO, class NodeType >
    class TsqrAdaptor< S, LO, GO, 
		       Tpetra::MultiVector< S, LO, GO, NodeType >,
		       TsqrTypeAdaptor< S, LO, GO, Tpetra::MultiVector< S, LO, GO, NodeType > > >
    {
    public:
      typedef Tpetra::MultiVector< S, LO, GO, NodeType > MV;
      typedef TsqrTypeAdaptor< S, LO, GO, MV > TsqrTypeAdaptorType;

      typedef S   scalar_type;
      typedef LO  local_ordinal_type;
      typedef GO  global_ordinal_type;
      typedef MV  multivector_type;

      typedef typename TsqrTypeAdaptorType::node_tsqr_type node_tsqr_type;
      typedef typename TsqrTypeAdaptorType::tsqr_type      tsqr_type;

      typedef Teuchos::RCP< node_tsqr_type >           node_tsqr_ptr;
      typedef Teuchos::RCP< tsqr_type >                tsqr_ptr;
      typedef Teuchos::RCP< const Teuchos::Comm<int> > comm_ptr;
      typedef Teuchos::RCP< MessengerBase< S > >       messenger_ptr;
      typedef typename tsqr_type::FactorOutput         factor_output_type;
      typedef Teuchos::SerialDenseMatrix< LO, S >      dense_matrix_type;

      TsqrAdaptor (const comm_ptr& comm,
		   const Teuchos::ParameterList& plist)
      {
	typedef TsqrFactory< LO, S, node_tsqr_type, tsqr_type > factory_type;
	factory_type::makeTsqr (plist, comm, pMessenger_, pNodeTsqr_, pTsqr_);
      }

      factor_output_type
      factor (multivector_type& A,
	      dense_matrix_type& R,
	      const bool contiguousCacheBlocks = false)
      {
	local_ordinal_type nrowsLocal, ncols, LDA;
	fetch_MV_dims (A, nrowsLocal, ncols, LDA);

	// This is guaranteed to be _correct_ for any Node type, but
	// won't necessary be efficient.  The desired model is that
	// A_local requires no copying.
	Teuchos::ArrayRCP< scalar_type > A_local = A.get1dViewNonConst();

	// Reshape R if necessary.  This operation zeros out all the
	// entries of R, which is what we want anyway.
	if (R.numRows() != ncols || R.numCols() != ncols)
	  {
	    if (0 != R.shape (ncols, ncols))
	      throw std::runtime_error ("Failed to reshape matrix R");
	  }
	return pTsqr_->factor (nrowsLocal, ncols, A_local.get(), LDA, 
			       R.values(), R.stride(), contiguousCacheBlocks);
      }

      void 
      explicitQ (const multivector_type& Q_in, 
		 const factor_output_type& factorOutput,
		 multivector_type& Q_out, 
		 const bool contiguousCacheBlocks = false)
      {
	using Teuchos::ArrayRCP;

	local_ordinal_type nrowsLocal, ncols_in, LDQ_in;
	fetch_MV_dims (Q_in, nrowsLocal, ncols_in, LDQ_in);
	local_ordinal_type nrowsLocal_out, ncols_out, LDQ_out;
	fetch_MV_dims (Q_out, nrowsLocal_out, ncols_out, LDQ_out);

	if (nrowsLocal_out != nrowsLocal)
	  {
	    std::ostringstream os;
	    os << "TSQR explicit Q: input Q factor\'s node-local part has a di"
	      "fferent number of rows (" << nrowsLocal << ") than output Q fac"
	      "tor\'s node-local part (" << nrowsLocal_out << ").";
	    throw std::runtime_error (os.str());
	  }
	ArrayRCP< const scalar_type > pQin = Q_in.get1dView();
	ArrayRCP< scalar_type > pQout = Q_out.get1dViewNonConst();
	pTsqr_->explicit_Q (nrowsLocal, 
			    ncols_in, pQin.get(), LDQ_in, 
			    factorOutput,
			    ncols_out, pQout.get(), LDQ_out,
			    contiguousCacheBlocks);
      }

      void 
      cacheBlock (const multivector_type& A_in,
		  multivector_type& A_out)
      {
	using Teuchos::ArrayRCP;

	local_ordinal_type nrowsLocal, ncols, LDA_in;
	fetch_MV_dims (A_in, nrowsLocal, ncols, LDA_in);
	local_ordinal_type nrowsLocal_out, ncols_out, LDA_out;
	fetch_MV_dims (A_out, nrowsLocal_out, ncols_out, LDA_out);

	if (nrowsLocal_out != nrowsLocal)
	  {
	    std::ostringstream os;
	    os << "TSQR cache block: the input matrix\'s node-local part has a"
	      " different number of rows (" << nrowsLocal << ") than the outpu"
	      "t matrix\'s node-local part (" << nrowsLocal_out << ").";
	    throw std::runtime_error (os.str());
	  }
	else if (ncols_out != ncols)
	  {
	    std::ostringstream os;
	    os << "TSQR cache block: the input matrix\'s node-local part has a"
	      " different number of columns (" << ncols << ") than the output "
	      "matrix\'s node-local part (" << ncols_out << ").";
	    throw std::runtime_error (os.str());
	  }
	ArrayRCP< const scalar_type > pA_in = A_in.get1dView();
	ArrayRCP< scalar_type > pA_out = A_out.get1dViewNonConst();
	pTsqr_->cache_block (nrowsLocal, ncols, pA_out.get(), 
			     pA_in.get(), LDA_in);
      }


      void 
      unCacheBlock (const multivector_type& A_in,
		    multivector_type& A_out)
      {
	using Teuchos::ArrayRCP;

	local_ordinal_type nrowsLocal, ncols, LDA_in;
	fetch_MV_dims (A_in, nrowsLocal, ncols, LDA_in);
	local_ordinal_type nrowsLocal_out, ncols_out, LDA_out;
	fetch_MV_dims (A_out, nrowsLocal_out, ncols_out, LDA_out);

	if (nrowsLocal_out != nrowsLocal)
	  {
	    std::ostringstream os;
	    os << "TSQR un-cache-block: the input matrix\'s node-local part ha"
	      "s a different number of rows (" << nrowsLocal << ") than the ou"
	      "tput matrix\'s node-local part (" << nrowsLocal_out << ").";
	    throw std::runtime_error (os.str());
	  }
	else if (ncols_out != ncols)
	  {
	    std::ostringstream os;
	    os << "TSQR cache block: the input matrix\'s node-local part has a"
	      " different number of columns (" << ncols << ") than the output "
	      "matrix\'s node-local part (" << ncols_out << ").";
	    throw std::runtime_error (os.str());
	  }
	ArrayRCP< const scalar_type > pA_in = A_in.get1dView();
	ArrayRCP< scalar_type > pA_out = A_out.get1dViewNonConst();
	pTsqr_->un_cache_block (nrowsLocal, ncols, pA_out.get(), 
				LDA_out, pA_in.get());
      }

    private:
      void
      fetch_MV_dims (const multivector_type& A,
		     local_ordinal_type& nrowsLocal, 
		     local_ordinal_type& ncols, 
		     local_ordinal_type& LDA)
      {
	nrowsLocal = A.getLocalLength();
	ncols = A.getNumVectors();
	if (nrowsLocal < ncols)
	  {
	    std::ostringstream os;
	    os << "TSQR: The local component of the input matrix has fewer row"
	      "s (" << nrowsLocal << ") than columns (" << ncols << ").  TSQR "
	      "does not support this case.";
	    throw std::runtime_error (os.str());
	  }
	if (! A.isConstantStride())
	  {
	    // FIXME (mfh 14 June 2010) Storage of A uses nonconstant
	    // stride internally, but that doesn't necessarily mean we
	    // can't run TSQR.  It depends on what get1dViewNonConst()
	    // returns.  If it's copied and packed into a matrix with
	    // constant stride, then we are free to run TSQR.
	    std::ostringstream os;
	    os << "TSQR does not support Tpetra::MultiVector inputs that do no"
	      "t have constant stride.";
	    throw std::runtime_error (os.str());
	  }
	LDA = A.getStride();
      }

      messenger_ptr pMessenger_;
      node_tsqr_ptr pNodeTsqr_;
      tsqr_ptr pTsqr_;
    };
  } // namespace Trilinos
} // namespace TSQR

#endif // __TSQR_Trilinos_TsqrAdaptor_Tpetra_MultiVector_SerialNode_hpp
