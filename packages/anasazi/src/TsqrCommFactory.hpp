#ifndef __TSQR_Trilinos_TsqrCommFactory_hpp
#define __TSQR_Trilinos_TsqrCommFactory_hpp

#include "Teuchos_RCP.hpp"
#include "TsqrTypeAdaptor.hpp"
#include "Tsqr_MessengerBase.hpp"

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

namespace TSQR {
  namespace Trilinos {

    template< class S, class LO, class GO, class MV >
    class TsqrCommFactory {
    public:
      typedef S scalar_type;
      typedef LO local_ordinal_type;
      typedef GO global_ordinal_type;
      typedef MV multivector_type;

      typedef TsqrTypeAdaptor< S, LO, GO, MV >    type_adaptor;
      typedef typename type_adaptor::comm_type    comm_type;
      typedef Teuchos::RCP< comm_type >           comm_ptr;
      typedef Teuchos::RCP< MessengerBase< S > >  scalar_messenger_ptr;
      typedef Teuchos::RCP< MessengerBase< LO > > ordinal_messenger_ptr;

      virtual void
      makeMessengers (const comm_ptr& comm,
		      scalar_messenger_ptr& scalarMessenger,
		      ordinal_messenger_ptr& ordinalMessenger) = 0;
    };
  } // namespace Trilinos
} // namespace TSQR

#include "TsqrCommFactory_Tpetra.hpp"
#include "TsqrCommFactory_Epetra.hpp"

#endif // __TSQR_Trilinos_TsqrCommFactory_hpp
