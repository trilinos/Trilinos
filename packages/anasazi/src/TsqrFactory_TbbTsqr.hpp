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
    class TbbTsqrFactory :
      public TsqrFactory< LO, S, TSQR::TBB::TbbTsqr< LO, S >, 
			  Tsqr< LO, S, TSQR::TBB::TbbTsqr< LO, S > > > 
    {
    public:
      // Help C++ pull in the typedefs from the base class.  C++ needs
      // help when both the base and the derived classes are
      // templated.
      typedef TsqrFactory< LO, S, TSQR::TBB::TbbTsqr< LO, S >, Tsqr< LO, S, TSQR::TBB::TbbTsqr< LO, S > > > base_type;

      typedef typename base_type::node_tsqr_type node_tsqr_type;
      typedef typename base_type::node_tsqr_ptr  node_tsqr_ptr;

      TbbTsqrFactory () {}
      virtual ~TbbTsqrFactory () {}

    private:
      /// \param [in] plist Parameter list with the following keys:
      ///   \li "cacheBlockSize": Cache size in bytes.  Default is
      ///       zero, which means that TSQR will guess a reasonable
      ///       cache size.  The size should correspond to that of the
      ///       largest cache that is private to each CPU core, if
      ///       such a private cache exists; alternately, it should
      ///       correspond to the amount of shared cache, divided by
      ///       the number of cores sharing that cache.
      ///   \li "numCores": If node_tsqr_type requires this (TbbTsqr
      ///       does, for Intel Threading Building Blocks intranode
      ///       parallelism), it should be set to the number of CPU
      ///       cores that are to participate in the intranode
      ///       parallel part of TSQR.  If node_tsqr_type does not
      ///       require this, then the parameter is ignored.
      virtual node_tsqr_ptr
      makeNodeTsqr (const Teuchos::ParameterList& plist) const
      {
	using Teuchos::Exceptions::InvalidParameter;
	int numCores = 1;
	size_t cacheBlockSize = 0;

	try {
	  const std::string numCoresParamName ("numCores");
	  numCores = plist.get<int>(numCoresParamName);
	} catch (InvalidParameter&) {
	  numCores = 1;
	}

	try {
	  const std::string cacheBlockSizeParamName ("cacheBlockSize");
	  cacheBlockSize = plist.get<size_t>(cacheBlockSizeParamName);
	} catch (InvalidParameter&) {
	  // Intranode TSQR interprets cacheBlockSize==0 as "set cache
	  // block size to a reasonable positive default value."
	  cacheBlockSize = 0;
	}

	node_tsqr_ptr node_tsqr (new node_tsqr_type (numCores, cacheBlockSize));
	return node_tsqr;
      }
    };
#endif // HAVE_KOKKOS_TBB

  } // namespace Trilinos
} // namespace TSQR


#endif // __TSQR_Trilinos_TsqrFactory_TbbTsqr_hpp
