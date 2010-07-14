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

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

namespace TSQR {
  namespace Trilinos {

    template< class S, class LO, class GO, class Node >
    class TsqrTpetraAdaptor :
      public TsqrAdaptor< S, LO, GO, Tpetra::MultiVector< S, LO, GO, Node > >
    {
    public:
      typedef Node node_type;
      // mfh 14 Jul 2010: This is because C++ is too stupid to know
      // how to inherit typedefs from a templated base class, when the
      // derived class itself is templated.  Go figure.
      typedef TsqrAdaptor< S, LO, GO, Tpetra::MultiVector< S, LO, GO, Node > > base_type;
      typedef typename base_type::comm_ptr comm_ptr;
      typedef typename base_type::multivector_type multivector_type;
      typedef typename base_type::scalar_type scalar_type;
      typedef typename base_type::local_ordinal_type local_ordinal_type;

      TsqrTpetraAdaptor (const comm_ptr& comm,
			 const Teuchos::ParameterList& plist) : 
	base_type (comm, plist)
      {}

    private:
      virtual void
      fetchDims (const multivector_type& A,
		 local_ordinal_type& nrowsLocal, 
		 local_ordinal_type& ncols, 
		 local_ordinal_type& LDA) const
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

      virtual Teuchos::ArrayRCP< scalar_type > 
      fetchNonConstView (multivector_type& A) const
      {
	return A.get1dViewNonConst();
      }

      virtual Teuchos::ArrayRCP< const scalar_type > 
      fetchConstView (const multivector_type& A) const
      {
	return A.get1dView();
      }
    };

  } // namespace Trilinos
} // namespace TSQR

#endif // __TSQR_Trilinos_TsqrAdaptor_Tpetra_MultiVector_SerialNode_hpp
