#ifndef __TSQR_Trilinos_TsqrCommFactory_Tpetra_hpp
#define __TSQR_Trilinos_TsqrCommFactory_Tpetra_hpp

#include "Tsqr_TpetraMessenger.hpp"

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

namespace TSQR {
  namespace Trilinos {

    template< class S, class LO, class GO, class Node >
    class TsqrCommFactory< S, LO, GO, Tpetra::MultiVector< S, LO, GO, Node > > : 
      public TsqrCommFactory< S, LO, GO, Tpetra::MultiVector< S, LO, GO, Node > >
    {
    public:
      // MV def looks circular but is not, because of how C++ inherits typedefs.

      typedef Tpetra::MultiVector< S, LO, GO, Node > MV;
      typedef TsqrCommFactory< S, LO, GO, MV > base_type;

      typedef typename base_type::comm_ptr comm_ptr;
      typedef typename base_type::scalar_messenger_ptr scalar_messenger_ptr;
      typedef typename base_type::ordinal_messenger_ptr ordinal_messenger_ptr;

      virtual void
      makeMessengers (const comm_ptr& comm,
		      scalar_messenger_ptr& scalarMessenger,
		      ordinal_messenger_ptr& ordinalMessenger)
      {
	scalarMessenger = Teuchos::rcp (new TSQR::Trilinos::TpetraMessenger< S > (comm));
	ordinalMessenger = Teuchos::rcp (new TSQR::Trilinos::TpetraMessenger< LO > (comm));
      }
    };
  } // namespace Trilinos
} // namespace TSQR

#include "TsqrCommFactory_Tpetra.hpp"
#include "TsqrCommFactory_Epetra.hpp"

#endif // __TSQR_Trilinos_TsqrCommFactory_Tpetra_hpp
