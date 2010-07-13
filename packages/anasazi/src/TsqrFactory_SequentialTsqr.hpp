#ifndef __TSQR_Trilinos_TsqrFactory_SequentialTsqr_hpp
#define __TSQR_Trilinos_TsqrFactory_SequentialTsqr_hpp

#include "Tsqr_SequentialTsqr.hpp"
#include "Tsqr.hpp"
#include "Teuchos_ParameterListExceptions.hpp"

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

namespace TSQR {
  namespace Trilinos {

    template< class LO, class S >
    class TsqrFactory< LO, S, SequentialTsqr< LO, S >, Tsqr< LO, S, SequentialTsqr< LO, S > > > {
    public:
      typedef Teuchos::RCP< const Teuchos::Comm<int> > comm_ptr;
      typedef Teuchos::RCP< MessengerBase< S > >       messenger_ptr;
      typedef SequentialTsqr< LO, S >                  node_tsqr_type;
      typedef Tsqr< LO, S, SequentialTsqr< LO, S > >   tsqr_type;
      typedef Teuchos::RCP< node_tsqr_type >           node_tsqr_ptr;
      typedef Teuchos::RCP< tsqr_type >                tsqr_ptr;

      static void
      makeTsqr (const Teuchos::ParameterList& plist,
		const comm_ptr& comm,
		messenger_ptr& messenger,
		node_tsqr_ptr& node_tsqr,
		tsqr_ptr& tsqr)
      {
	messenger = makeMessenger (comm);
	node_tsqr = makeNodeTsqr (plist);
	tsqr = tsqr_ptr (new tsqr_type (*node_tsqr, messenger.get()));
      }

    private:
      TsqrFactory ();
      ~TsqrFactory ();
      TsqrFactory (const TsqrFactory&);
      TsqrFactory& operator= (const TsqrFactory&);

      static node_tsqr_ptr
      makeNodeTsqr (const Teuchos::ParameterList& plist);

      static messenger_ptr 
      makeMessenger (const comm_ptr& comm)
      {
	using TSQR::Trilinos::TrilinosMessenger;
	return messenger_ptr (new TrilinosMessenger< S > (comm));
      }
    };

    template< class LO, class S >
    typename TsqrFactory< LO, S, SequentialTsqr< LO, S >, Tsqr< LO, S, SequentialTsqr< LO, S > > >::node_tsqr_ptr
    TsqrFactory< LO, S, SequentialTsqr< LO, S >, Tsqr< LO, S, SequentialTsqr< LO, S > > >::
    makeNodeTsqr (const Teuchos::ParameterList& plist)
    {
      using Teuchos::Exceptions::InvalidParameter;
      size_t cache_block_size = 0;
      // All this try/catch stuff is because the C++ compiler can't
      // deduce the right two-argument get() function (second argument
      // would be the default).
      try {
	cache_block_size = plist.get< size_t > (std::string("cache-block-size"));
      } catch (InvalidParameter&) {
	cache_block_size = 0;
      }
	
      node_tsqr_ptr node_tsqr (new node_tsqr_type (cache_block_size));
      return node_tsqr;
    }

  }
}

#endif // __TSQR_Trilinos_TsqrFactory_SequentialTsqr_hpp
