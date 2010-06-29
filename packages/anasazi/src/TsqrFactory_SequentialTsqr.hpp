#ifndef __TSQR_Trilinos_TsqrFactory_SequentialTsqr_hpp
#define __TSQR_Trilinos_TsqrFactory_SequentialTsqr_hpp

#include "Tsqr_SequentialTsqr.hpp"
#include "Tsqr.hpp"

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

namespace TSQR {
  namespace Trilinos {

    template< class LO, class S >
    class TsqrFactory< LO, S, SequentialTsqr< LO, S >, Tsqr< LO, S, SequentialTsqr< LO, S > > > {
    public:
      typedef Teuchos::RCP< const Teuchos::Comm<int> > comm_ptr;
      typedef Teuchos::RCP< const MessengerBase< S > > messenger_ptr;
      typedef SequentialTsqr< LO, S >                  node_tsqr_type;
      typedef Tsqr< LO, S, node_tsqr_type >            tsqr_type;
      typedef Teuchos::RCP< node_tsqr_type >           node_tsqr_ptr;
      typedef Teuchos::RCP< tsqr_type >                tsqr_ptr;

      static void
      makeTsqr (const Teuchos::ParameterList& plist,
		messenger_ptr& messenger,
		node_tsqr_ptr& node_tsqr,
		tsqr_ptr& tsqr)
      {
	messenger_ptr messenger (new TrilinosMessenger (comm));
	node_tsqr_ptr node_tsqr = makeNodeTsqr (plist);
	tsqr_ptr tsqr (node_tsqr.get(), messenger.get());
	return tsqr_wrapper_type (messenger, node_tsqr, tsqr);
      }

    private:
      // None of these are implemented: you can't construct, destruct,
      // or copy a TsqrFactory().
      TsqrFactory ();
      ~TsqrFactory ();
      TsqrFactory (const TsqrFactory&);
      TsqrFactory& operator= (const TsqrFactory&);

      static node_tsqr_ptr
      makeNodeTsqr (const Teuchos::ParameterList& plist);
    };

    template< class LO, class S >
    TsqrFactory< LO, S, SequentialTsqr< LO, S >, Tsqr< LO, S, SequentialTsqr< LO, S > > >::node_tsqr_type
    TsqrFactory< LO, S, SequentialTsqr< LO, S >, Tsqr< LO, S, SequentialTsqr< LO, S > > >::
    makeNodeTsqr (const Teuchos::ParameterList& plist)
    {
      const size_t cache_block_size = plist.get("cache_block_size", 0);
      node_tsqr_ptr node_tsqr (new node_tsqr_type (cache_block_size));
      return node_tsqr_ptr;
    }

  }
}

#endif // __TSQR_Trilinos_TsqrFactory_SequentialTsqr_hpp
