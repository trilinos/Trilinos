#ifndef __TSQR_Trilinos_TsqrAdaptor_Epetra_MultiVector_hpp
#define __TSQR_Trilinos_TsqrAdaptor_Epetra_MultiVector_hpp

/// \file TsqrAdaptor_Epetra_MultiVector.hpp
///
/// \warning Users should _not_ include this file directly.  Include
///   "TsqrAdaptor.hpp" instead.  If HAVE_ANASAZI_EPETRA is defined,
///   then this file will be included automatically.  If for some
///   reason you need to include this file directly, be sure to
///   include "TsqrAdaptor.hpp" first.

#include "Epetra_MultiVector.hpp"
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
			       R.stride(), pComm_.get(), contiguousCacheBlocks);
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
	const scalar_type* const pQ_in = Q_in.Values();
	scalar_type* const pQ_out = Q_out.Values();
	pTsqr_->explicit_Q (nrowsLocal, ncols_in, pQ_in, LDQ_in, 
			    ncols_out, pQ_out, LDQ_out, factorOutput, 
			    pComm_.get(), contiguousCacheBlocks);
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
	const scalar_type* const pA_in = A_in.Values();
	scalar_type* const pA_out = A_out.Values();
	pTsqr_->cache_block (nrowsLocal, ncols, pA_out.get(), 
			     pA_in.get(), LDA_in);
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
	const scalar_type* const pA_in = A_in.Values();
	scalar_type* const pA_out = A_out.Values();
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
	    // FIXME (mfh 14 June 2010) Storage of A uses nonconstant
	    // stride internally, but that doesn't necessarily mean we
	    // can't run TSQR.  It depends on what get1dViewNonConst()
	    // returns.  If it's copied and packed into a matrix with
	    // constant stride, then we are free to run TSQR.
	    std::ostringstream os;
	    os << "TSQR does not support Epetra_MultiVector inputs that do not"
	      " have constant stride.";
	    throw std::runtime_error (os.str());
	  }
	LDA = A.Stride();
      }
    };

  } // namespace Trilinos
} // namespace TSQR

#endif // __TSQR_Trilinos_TsqrAdaptor_Epetra_MultiVector_hpp
