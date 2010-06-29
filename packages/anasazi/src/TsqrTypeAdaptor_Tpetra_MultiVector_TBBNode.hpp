#ifndef __TSQR_Trilinos_TsqrTypeAdaptor_Tpetra_MultiVector_TBBNode_hpp
#define __TSQR_Trilinos_TsqrTypeAdaptor_Tpetra_MultiVector_TBBNode_hpp

/// \file TsqrTypeAdaptor_Tpetra_MultiVector_TBBNode.hpp
///
/// \warning Users should _not_ include this file directly.  Include
///   "TsqrTypeAdaptor.hpp" instead.  If HAVE_ANASAZI_TPETRA is
///   defined, then this file will be included automatically.  If for
///   some reason you need to include this file directly, be sure to
///   include "TsqrTypeAdaptor.hpp" first.
///
/// \note All content in this file is disabled unless HAVE_KOKKOS_TBB
///   is defined.

#include "Tpetra_MultiVector.hpp"
#ifdef HAVE_KOKKOS_TBB
#  include "Kokkos_TBBNode.hpp"
#  include "Tsqr_TbbTsqr.hpp"
#endif // HAVE_KOKKOS_TBB

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

namespace TSQR {
  namespace Trilinos {

#ifdef HAVE_KOKKOS_TBB
    template< class S, class LO, class GO >
    class TsqrTypeAdaptor< S, LO, GO, Tpetra::MultiVector< S, LO, GO, Kokkos::TBBNode > > {
    public:
      typedef S scalar_type;
      typedef LO local_ordinal_type;
      typedef GO global_ordinal_type;
      typedef Tpetra::MultiVector< S, LO, GO, Kokkos::TBBNode > multivector_type;

      typedef TSQR::TbbTsqr< LO, S > node_tsqr_type;
      typedef TSQR::Tsqr< LO, S, node_tsqr_type > tsqr_type;
      typedef TsqrFactory< LO, S, node_tsqr_type, tsqr_type > factory_type;
    };
#endif // HAVE_KOKKOS_TBB

  } // namespace Trilinos
} // namespace TSQR


#endif // __TSQR_Trilinos_TsqrTypeAdaptor_Tpetra_MultiVector_TBBNode_hpp
