#ifndef __TSQR_Trilinos_TsqrCommFactory_Epetra_hpp
#define __TSQR_Trilinos_TsqrCommFactory_Epetra_hpp

#include "Tsqr_EpetraMessenger.hpp"

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

namespace TSQR {
  namespace Trilinos {

    class EpetraCommFactory :
      public CommFactory< double, int, int, Epetra_MultiVector >
    {
    public:
      typedef CommFactory< double, int, int, Epetra_MultiVector > base_type;
      typedef base_type::comm_ptr comm_ptr;
      typedef base_type::scalar_messenger_ptr scalar_messenger_ptr;
      typedef base_type::ordinal_messenger_ptr ordinal_messenger_ptr;

      EpetraCommFactory () {}
      virtual ~EpetraCommFactory () {}

      virtual void
      makeMessengers (const comm_ptr& comm,
		      scalar_messenger_ptr& scalarMessenger,
		      ordinal_messenger_ptr& ordinalMessenger)
      {
	using TSQR::Trilinos::EpetraMessenger;

	scalarMessenger = Teuchos::rcp (new EpetraMessenger< S > (comm));
	ordinalMessenger = Teuchos::rcp (new EpetraMessenger< LO > (comm));
      }
    };
  } // namespace Trilinos
} // namespace TSQR

#endif // __TSQR_Trilinos_TsqrCommFactory_Epetra_hpp
