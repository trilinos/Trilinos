#ifndef __TSQR_Trilinos_TsqrFactory_TbbTsqr_hpp
#define __TSQR_Trilinos_TsqrFactory_TbbTsqr_hpp

// mfh 28 Jun 2010: This should suffice to pull in the HAVE_KOKKOS_TBB
// define.
#include "Tpetra_MultiVector.hpp"

#ifdef HAVE_KOKKOS_TBB
#  include "TbbTsqr.hpp"
#endif // HAVE_KOKKOS_TBB

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

namespace TSQR {
  namespace Trilinos {

#ifdef HAVE_KOKKOS_TBB
    template< class LO, class S >
    class TsqrFactory< LO, S, TSQR::TBB::TbbTsqr< LO, S >, Tsqr< LO, S, TSQR::TBB::TbbTsqr< LO, S > > > {
    public:
      typedef Teuchos::RCP< const Teuchos::Comm<int> > comm_ptr;
      typedef Teuchos::RCP< const MessengerBase< S > > messenger_ptr;
      typedef TSQR::TBB::TbbTsqr< LO, S >              node_tsqr_type;
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
    TsqrFactory< LO, S, TSQR::TBB::TbbTsqr< LO, S >, Tsqr< LO, S, TbbTsqr< LO, S > > >::node_tsqr_type
    TsqrFactory< LO, S, TSQR::TBB::TbbTsqr< LO, S >, Tsqr< LO, S, TbbTsqr< LO, S > > >::
    makeNodeTsqr (const Teuchos::ParameterList& plist)
    {
      const int num_cores = plist.get("num_cores", 1);
      const size_t cache_block_size = plist.get("cache_block_size", 0);

      node_tsqr_ptr node_tsqr (new node_tsqr_type (num_cores, cache_block_size));
      return node_tsqr;
    }
#endif // HAVE_KOKKOS_TBB

  } // namespace Trilinos
} // namespace TSQR


#endif // __TSQR_Trilinos_TsqrFactory_TbbTsqr_hpp
