#ifndef __TSQR_Trilinos_TsqrFactory_SequentialTsqr_hpp
#define __TSQR_Trilinos_TsqrFactory_SequentialTsqr_hpp

#include "Tsqr_SequentialTsqr.hpp"
#include "Tsqr.hpp"
#include "Teuchos_ParameterListExceptions.hpp"

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

namespace TSQR {
  namespace Trilinos {

    /// \class SequentialTsqrFactory
    ///
    /// Subclass of TsqrFactory that knows how to instantiate
    /// SequentialTsqr for intranode TSQR.
    template< class LO, class S >
    class SequentialTsqrFactory :
      public TsqrFactory< LO, S, SequentialTsqr< LO, S >, 
			  Tsqr< LO, S, SequentialTsqr< LO, S > > >
    {
    public:
      typedef TsqrFactory< LO, S, SequentialTsqr< LO, S >, Tsqr< LO, S, SequentialTsqr< LO, S > > > base_type;

      // Pull in the typedefs from the base class.  C++ doesn't do
      // this when both the base and the derived classes are
      // templated.
      typedef typename base_type::node_tsqr_type node_tsqr_type;
      typedef typename base_type::node_tsqr_ptr  node_tsqr_ptr;

      SequentialTsqrFactory () {}
      virtual ~SequentialTsqrFactory () {}

    private:
      /// \param [in] plist Parameter list with the following key:
      ///   \li "cacheBlockSize": Cache size in bytes.  Default is
      ///       zero, which means that TSQR will guess a reasonable
      ///       cache size.  The size should correspond to that of the
      ///       largest cache that is private to each CPU core, if
      ///       such a private cache exists; alternately, it should
      ///       correspond to the amount of shared cache, divided by
      ///       the number of cores sharing that cache.
      virtual node_tsqr_ptr
      makeNodeTsqr (const Teuchos::ParameterList& plist) const
      {
	using Teuchos::Exceptions::InvalidParameter;
	size_t cacheBlockSize = 0;

	// All this try/catch stuff is because the C++ compiler can't
	// deduce the right two-argument get() function (second argument
	// would be the default).
	try {
	  const std::string cacheBlockSizeParamName ("cacheBlockSize");
	  cacheBlockSize = plist.get< size_t > (cacheBlockSizeParamName);
	} catch (InvalidParameter&) {
	  cacheBlockSize = 0;
	}
	
	node_tsqr_ptr node_tsqr (new node_tsqr_type (cacheBlockSize));
	return node_tsqr;
      }
    };

  } // namespace Trilinos
} // namespace TSQR

#endif // __TSQR_Trilinos_TsqrFactory_SequentialTsqr_hpp
