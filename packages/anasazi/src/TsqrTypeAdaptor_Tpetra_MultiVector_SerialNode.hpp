#ifndef __TSQR_Trilinos_TsqrTypeAdaptor_Tpetra_MultiVector_SerialNode_hpp
#define __TSQR_Trilinos_TsqrTypeAdaptor_Tpetra_MultiVector_SerialNode_hpp

/// \file TsqrTypeAdaptor_Tpetra_MultiVector_SerialNode.hpp
///
/// \warning Users should _not_ include this file directly.  Include
///   "TsqrTypeAdaptor.hpp" instead.  If HAVE_ANASAZI_TPETRA is
///   defined, then this file will be included automatically.  If for
///   some reason you need to include this file directly, be sure to
///   include "TsqrTypeAdaptor.hpp" first.

#include "Tpetra_MultiVector.hpp"
#include "Kokkos_SerialNode.hpp"
#include "Tsqr_SequentialTsqr.hpp"
#include "Tsqr.hpp"

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

namespace TSQR {
  namespace Trilinos {

    template< class S, class LO, class GO >
    class TsqrTypeAdaptor< S, LO, GO, Tpetra::MultiVector< S, LO, GO, Kokkos::SerialNode > > {
    public:
      typedef S scalar_type;
      typedef LO local_ordinal_type;
      typedef GO global_ordinal_type;
      typedef Tpetra::MultiVector< S, LO, GO, Kokkos::SerialNode > multivector_type;

      typedef TSQR::SequentialTsqr< LO, S > node_tsqr_type;
      typedef TSQR::Tsqr< LO, S, node_tsqr_type > tsqr_type;
      typedef SequentialTsqrFactory< local_ordinal_type, scalar_type > factory_type;
    };

  } // namespace Trilinos
} // namespace TSQR


#endif // __TSQR_Trilinos_TsqrTypeAdaptor_Tpetra_MultiVector_SerialNode_hpp
