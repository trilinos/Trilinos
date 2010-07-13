#ifndef __TSQR_Trilinos_TsqrAdaptor_Epetra_MultiVector_hpp
#define __TSQR_Trilinos_TsqrAdaptor_Epetra_MultiVector_hpp

/// \file TsqrAdaptor_Epetra_MultiVector.hpp
///
/// \warning Users should _not_ include this file directly.  Include
///   "TsqrAdaptor.hpp" instead.  If HAVE_ANASAZI_EPETRA is defined,
///   then this file will be included automatically.  If for some
///   reason you need to include this file directly, be sure to
///   include "TsqrAdaptor.hpp" first.

#include "Epetra_MultiVector.h" // sic (not .hpp)
#include <stdexcept>
#include <sstream>

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

namespace TSQR {
  namespace Trilinos {

    template<>
    class TsqrAdaptor< double, int, int, Epetra_MultiVector >
    {
    public:
      typedef double             scalar_type;
      typedef int                local_ordinal_type;
      typedef int                global_ordinal_type;
      typedef Epetra_MultiVector multivector_type;

      typedef TsqrTypeAdaptor< double, int, int, Epetra_MultiVector >::
      node_tsqr_type node_tsqr_type;
      typedef TsqrTypeAdaptor< double, int, int, Epetra_MultiVector >::
      tsqr_type      tsqr_type;

      typedef Teuchos::RCP< node_tsqr_type >            node_tsqr_ptr;
      typedef Teuchos::RCP< tsqr_type >                 tsqr_ptr;
      typedef Teuchos::RCP< const Teuchos::Comm<int> >  comm_ptr;
      typedef Teuchos::RCP< MessengerBase<double> >     messenger_ptr;
      typedef tsqr_type::FactorOutput                   factor_output_type;
      typedef Teuchos::SerialDenseMatrix<int, double>   dense_matrix_type;

      TsqrAdaptor (const comm_ptr& comm,
		   const Teuchos::ParameterList& plist)
      {
	// This typedef is an implementation detail.
	typedef TsqrFactory< int, double, node_tsqr_type, tsqr_type > factory_type;
	// plist and comm are inputs.
	// Construct *pMessenger_, *pNodeTsqr_, and *pTsqr_.
	factory_type::makeTsqr (plist, comm, pMessenger_, pNodeTsqr_, pTsqr_);
      }

      factor_output_type
      factor (multivector_type& A,
	      dense_matrix_type& R, 
	      const bool contiguousCacheBlocks = false)
      {
	local_ordinal_type nrowsLocal, ncols, LDA;
	fetch_MV_dims (A, nrowsLocal, ncols, LDA);
	scalar_type* const A_local = A.Values();

	// Reshape R if necessary.  This operation zeros out all the
	// entries of R, which is what we want anyway.
	if (R.numRows() != ncols || R.numCols() != ncols)
	  {
	    if (0 != R.shape (ncols, ncols))
	      throw std::runtime_error ("Failed to reshape matrix R");
	  }
	return pTsqr_->factor (nrowsLocal, ncols, A_local, LDA, R.values(),
			       R.stride(), contiguousCacheBlocks);
      }

      void 
      explicitQ (const multivector_type& Q_in, 
		 const factor_output_type& factorOutput,
		 multivector_type& Q_out, 
		 const bool contiguousCacheBlocks = false)
      {
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
	pTsqr_->explicit_Q (nrowsLocal, 
			    ncols_in, Q_in.Values(), LDQ_in, 
			    factorOutput,
			    ncols_out, Q_out.Values(), LDQ_out, 
			    contiguousCacheBlocks);
      }

      void 
      cacheBlock (const multivector_type& A_in,
		  multivector_type& A_out)
      {
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
	pTsqr_->cache_block (nrowsLocal, ncols, A_out.Values(),
			     A_in.Values(), LDA_in);
      }


      void 
      unCacheBlock (const multivector_type& A_in,
		    multivector_type& A_out)
      {
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
	pTsqr_->un_cache_block (nrowsLocal, ncols, A_out.Values(), LDA_out,
				A_in.Values());
      }

    private:
      void
      fetch_MV_dims (const multivector_type& A,
		     local_ordinal_type& nrowsLocal,
		     local_ordinal_type& ncols,
		     local_ordinal_type& LDA)
      {
	nrowsLocal = A.MyLength();
	ncols = A.NumVectors();
	if (nrowsLocal < ncols)
	  {
	    std::ostringstream os;
	    os << "TSQR: The local component of the input matrix has fewer row"
	      "s (" << nrowsLocal << ") than columns (" << ncols << ").  TSQR "
	      "does not support this case.";
	    throw std::runtime_error (os.str());
	  }
	if (! A.ConstantStride())
	  {
	    std::ostringstream os;
	    os << "TSQR does not support Epetra_MultiVector inputs that do not"
	      " have constant stride.";
	    throw std::runtime_error (os.str());
	  }
	LDA = A.Stride();
      }

      messenger_ptr pMessenger_;
      node_tsqr_ptr pNodeTsqr_;
      tsqr_ptr pTsqr_;
    };

  } // namespace Trilinos
} // namespace TSQR

#endif // __TSQR_Trilinos_TsqrAdaptor_Epetra_MultiVector_hpp
