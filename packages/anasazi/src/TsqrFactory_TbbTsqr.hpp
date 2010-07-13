#ifndef __TSQR_Trilinos_TsqrFactory_TbbTsqr_hpp
#define __TSQR_Trilinos_TsqrFactory_TbbTsqr_hpp

// mfh 28 Jun 2010: This should suffice to pull in the HAVE_KOKKOS_TBB
// define.
#include "Tpetra_MultiVector.hpp"
#include "Teuchos_ParameterListExceptions.hpp"

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
      typedef Teuchos::RCP< MessengerBase< S > >       messenger_ptr;
      typedef TSQR::TBB::TbbTsqr< LO, S >              node_tsqr_type;
      typedef Tsqr< LO, S, node_tsqr_type >            tsqr_type;
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
	tsqr = tsqr_ptr (new tsqr_type (node_tsqr.get(), messenger.get()));
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
    typename TsqrFactory< LO, S, TSQR::TBB::TbbTsqr< LO, S >, Tsqr< LO, S, TSQR::TBB::TbbTsqr< LO, S > > >::node_tsqr_ptr
    TsqrFactory< LO, S, TSQR::TBB::TbbTsqr< LO, S >, Tsqr< LO, S, TSQR::TBB::TbbTsqr< LO, S > > >::
    makeNodeTsqr (const Teuchos::ParameterList& plist)
    {
      using Teuchos::Exceptions::InvalidParameter;
      int num_cores = 1;
      size_t cache_block_size = 0;
      const std::string cacheBlockSizeParamName ("cacheBlockSize");
      const std::string numCoresParamName ("numCores");

      try {
	num_cores = plist.get<int>(numCoresParamName);
      } catch (InvalidParameter&) {
	num_cores = 1;
      }

      try {
	cache_block_size = plist.get<size_t>(cacheBlockSizeParamName);
      } catch (InvalidParameter&) {
	// Intranode TSQR interprets cache_block_size==0 as "set cache
	// block size to a reasonable positive default value."
	cache_block_size = 0;
      }

      node_tsqr_ptr node_tsqr (new node_tsqr_type (num_cores, cache_block_size));
      return node_tsqr;
    }
#endif // HAVE_KOKKOS_TBB

  } // namespace Trilinos
} // namespace TSQR


#endif // __TSQR_Trilinos_TsqrFactory_TbbTsqr_hpp
