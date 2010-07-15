#ifndef __TSQR_Trilinos_TsqrCommFactory_Tpetra_hpp
#define __TSQR_Trilinos_TsqrCommFactory_Tpetra_hpp

#include "Tsqr_TpetraMessenger.hpp"

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

namespace TSQR {
  namespace Trilinos {

    template< class S, class LO, class GO, class Node >
    class TpetraCommFactory< S, LO, GO, Tpetra::MultiVector< S, LO, GO, Node > > : 
      public CommFactory< S, LO, GO, Tpetra::MultiVector< S, LO, GO, Node > >
    {
    public:
      // MV def looks circular but is not, because C++ does not pass
      // typedefs from a base class to a derived class when both the
      // base and derived classes are templated.  MV here is just a
      // shorthand so that the base_type typedef fits in one line.
      typedef Tpetra::MultiVector< S, LO, GO, Node > MV;
      typedef CommFactory< S, LO, GO, MV > base_type;

      // C++ doesn't pass typedefs from a templated base class to a
      // templated derived class.  Here, we avoid repeating code from
      // the base class by pulling in its typedefs into the derived
      // class.
      typedef typename base_type::comm_ptr comm_ptr;
      typedef typename base_type::scalar_messenger_ptr scalar_messenger_ptr;
      typedef typename base_type::ordinal_messenger_ptr ordinal_messenger_ptr;

      TpetraCommFactory () {}
      virtual ~TpetraCommFactory () {}

      virtual void
      makeMessengers (const comm_ptr& comm,
		      scalar_messenger_ptr& scalarMessenger,
		      ordinal_messenger_ptr& ordinalMessenger)
      {
	using TSQR::Trilinos::TpetraMessenger;

	scalarMessenger = Teuchos::rcp (new TpetraMessenger< S > (comm));
	ordinalMessenger = Teuchos::rcp (new TpetraMessenger< LO > (comm));
      }
    };
  } // namespace Trilinos
} // namespace TSQR

#endif // __TSQR_Trilinos_TsqrCommFactory_Tpetra_hpp
