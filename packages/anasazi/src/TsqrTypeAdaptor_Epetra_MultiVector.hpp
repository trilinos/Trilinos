#ifndef __TSQR_Trilinos_TsqrTypeAdaptor_Epetra_MultiVector_hpp
#define __TSQR_Trilinos_TsqrTypeAdaptor_Epetra_MultiVector_hpp

/// \file TsqrTypeAdaptor_Epetra_MultiVector.hpp
///
/// \warning Users should _not_ include this file directly.  Include
///   "TsqrTypeAdaptor.hpp" instead.  If HAVE_ANASAZI_EPETRA is
///   defined, then this file will be included automatically.  If for
///   some reason you need to include this file directly, be sure to
///   include "TsqrTypeAdaptor.hpp" first.

#include "Epetra_MultiVector.h" // sic (not .hpp)
#include "Tsqr_SequentialTsqr.hpp"
#include "Tsqr.hpp"

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

namespace TSQR {
  namespace Trilinos {

    template<>
    class TsqrTypeAdaptor< double, int, int, Epetra_MultiVector > {
    public:
      typedef double scalar_type;
      typedef int local_ordinal_type;
      typedef int global_ordinal_type;
      typedef Epetra_MultiVector multivector_type;

      typedef TSQR::SequentialTsqr< int, double > node_tsqr_type;
      typedef TSQR::Tsqr< int, double, node_tsqr_type > tsqr_type;
      typedef TsqrFactory< int, double, node_tsqr_type, tsqr_type > factory_type;
    };

  } // namespace Trilinos
} // namespace TSQR


#endif // __TSQR_Trilinos_TsqrTypeAdaptor_Epetra_MultiVector_hpp
