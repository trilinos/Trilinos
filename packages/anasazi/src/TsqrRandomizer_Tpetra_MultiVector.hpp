#ifndef __TSQR_Trilinos_Randomizer_Tpetra_MultiVector_hpp
#define __TSQR_Trilinos_Randomizer_Tpetra_MultiVector_hpp

/// \file TsqrRandomizer_Tpetra_MultiVector.hpp
///
/// \warning Users should _not_ include this file directly.  Include
///   "TsqrRandomizer.hpp" instead.  If HAVE_ANASAZI_TPETRA is defined,
///   then this file will be included automatically.  If for some
///   reason you need to include this file directly, be sure to
///   include "TsqrRandomizer.hpp" first.

#include "Tpetra_MultiVector.hpp"
#ifdef HAVE_KOKKOS_TBB
#  include "Kokkos_TBBNode.hpp"
#endif // HAVE_KOKKOS_TBB

#include <stdexcept>
#include <sstream>

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

namespace TSQR {
  namespace Trilinos {

    template< class S, class LO, class GO, class NodeType, class Gen >
    class TpetraRandomizer :
      public Randomizer< S, LO, GO, 
			 Tpetra::MultiVector< S, LO, GO, NodeType >,
			 Gen >
    {
    public:
      // mfh 14 Jul 2010: I'm uncomfortable naming this "MV" because
      // the base class already has MV defined in it.  I don't want to
      // go all circular (MV -> multivector_type -> MV).  TMV here is
      // just a way to make the base_type typedef fit on one line.
      typedef Tpetra::MultiVector< S, LO, GO, NodeType > TMV;

      // The C++ compiler needs help inheriting typedefs, when both
      // the base class and the derived class are templated.
      typedef Randomizer< double, int, int, TMV, Gen > base_type;
      typedef typename base_type::multivector_type multivector_type;
      typedef typename base_type::local_ordinal_type local_ordinal_type;
      typedef typename base_type::scalar_type scalar_type;

      typedef typename base_type::normalgen_ptr normalgen_ptr;
      typedef typename base_type::ordinal_messenger_ptr ordinal_messenger_ptr;
      typedef typename base_type::scalar_messenger_ptr scalar_messenger_ptr;

      /// \brief Constructor
      ///
      TpetraRandomizer (const normalgen_ptr& pGen,
			const ordinal_messenger_ptr& pOrdinalMess,
			const scalar_messenger_ptr& pScalarMess) : 
	base_type (pGen, pOrdinalMess, pScalarMess) 
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
	    os << "The local component of the input matrix has fewer row"
	      "s (" << nrowsLocal << ") than columns (" << ncols << ").  "
	      "TSQR::Trilinos::TpetraRandomizer does not support this case.";
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
	    os << "TSQR::Trilinos::TpetraRandomizer does not support "
	      "Tpetra::MultiVector inputs that do not have constant stride.";
	    throw std::runtime_error (os.str());
	  }
	LDA = A.getStride();
      }

      virtual Teuchos::ArrayRCP< scalar_type > 
      fetchNonConstView (multivector_type& A) const
      {
	// This won't be efficient if the entries of A live in a
	// different memory space (e.g., GPU node), but it will be
	// correct for any Tpetra::MultiVector type.
	return A.get1dViewNonConst();
      }
    };
  } // namespace Trilinos
} // namespace TSQR

#endif // __TSQR_Trilinos_Randomizer_Tpetra_MultiVector_hpp
