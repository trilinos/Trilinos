#ifndef __TSQR_Trilinos_TsqrFactory_hpp
#define __TSQR_Trilinos_TsqrFactory_hpp

#include "Teuchos_Comm.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RCP.hpp"

#include "Tsqr_MessengerBase.hpp"

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

namespace TSQR {
  namespace Trilinos {

    /// \class TsqrFactory
    ///
    /// Using TSQR requires holding on to three different RCPs:
    /// \li MessengerBase<S>
    /// \li node_tsqr_type
    /// \li tsqr_type
    /// They are returned by makeTsqr().  Pass it a ParameterList
    /// appropriate to the type.  All take a "cache_block_size"
    /// parameter (size of the cache in bytes), and the TBB
    /// node-parallel TSQR variant takes a "num_cores" parameter
    /// (number of cores per MPI process to use in TBB node-parallel
    /// TSQR).
    template< class LO, class S, class NodeTsqrType, class TsqrType >
    class TsqrFactory {
    public:
      typedef Teuchos::RCP< const Teuchos::Comm<int> > comm_ptr;
      typedef Teuchos::RCP< const MessengerBase< S > > messenger_ptr;
      typedef NodeTsqrType                             node_tsqr_type;
      typedef TsqrType                                 tsqr_type;
      typedef Teuchos::RCP< node_tsqr_type >           node_tsqr_ptr;
      typedef Teuchos::RCP< tsqr_type >                tsqr_ptr;

      /// \brief Instantiate and return the objects required by TSQR
      ///
      /// Instantiate and return (through the output arguments) the
      /// three objects required by TSQR.
      ///
      /// \param plist [in] Parameter list with two keys:  
      ///   \li "cache_block_size": Cache size in bytes.  Default is
      ///       zero, which means that TSQR will guess a reasonable
      ///       cache size.  The size should correspond to that of the
      ///       largest cache that is private to each CPU core, if
      ///       such a private cache exists; alternately, it should
      ///       correspond to the amount of shared cache, divided by
      ///       the number of cores sharing that cache.
      ///   \li "num_cores": If node_tsqr_type requires this (TbbTsqr
      ///       does, for Intel Threading Building Blocks intranode
      ///       parallelism), it should be set to the number of CPU
      ///       cores that are to participate in the intranode
      ///       parallel part of TSQR.  If node_tsqr_type does not
      ///       require this, then the parameter is ignored.
      /// \param messenger [out] On output, points to the
      ///   MessengerBase< S > object that TSQR will use for internode
      ///   communication.
      /// \param node_tsqr [out] On output, points to the
      ///   node_tsqr_type object that TSQR will use for the intranode
      ///   part of its computations.
      /// \param tsqr [out] On output, points to the node_tsqr_type
      ///   object that TSQR will use for the internode part of its
      ///   computations.
      static void
      makeTsqr (const Teuchos::ParameterList& plist,
		messenger_ptr& messenger,
		node_tsqr_ptr& node_tsqr,
		tsqr_ptr& tsqr);

    private:
      /// Constructor not implemented, so you can't construct a
      /// TsqrFactory instance; you have to use the static interface.
      TsqrFactory ();
      /// Destructor not implemented, so you can't construct a
      /// TsqrFactory instance; you have to use the static interface.
      ~TsqrFactory ();
      /// Copy constructor not implemented, so you can't construct a
      /// (copy of a) TsqrFactory instance.
      TsqrFactory (const TsqrFactory&);
      /// Assignment operator not implemented, so you can't construct
      /// a (copy of a) TsqrFactory instance.
      TsqrFactory& operator= (const TsqrFactory&);

      /// \brief Instantiate and return TSQR's intranode object
      ///
      /// \param plist [in] Same as the epinonymous input of makeTsqr()
      ///
      /// \return (Smart) pointer to the node_tsqr_type object that
      ///   TSQR will use for the intranode part of its computations
      ///
      /// \note For implementers: this is the interesting method.
      ///   makeTsqr() is "generic," more or less.  makeNodeTsqr() is
      ///   the part that varies a lot for different node_tsqr_type
      ///   types.  This is also the method that reads parameters from
      ///   the plist.  This pattern is the compile-time polymorphism
      ///   equivalent of the "Non-Virtual Interface" (NVI) idiom,
      ///   where the "virtual" methods (here, the methods that vary
      ///   for different template parameters) are private, and the
      ///   "nonvirtual" methods (here, the methods that are the same
      ///   for different template parameters) are part of the public
      ///   interface.
      static node_tsqr_ptr
      makeNodeTsqr (const Teuchos::ParameterList& plist);
    };
  } // namespace Trilinos
} // namespace TSQR

#include "TsqrFactory_SequentialTsqr.hpp"
// mfh 28 Jun 2010: this only does anything if HAVE_KOKKOS_TBB is
// defined.  The file included below includes Tpetra_MultiVector.hpp
// in order to get that.
#include "TsqrFactory_TbbTsqr.hpp"

#endif // __TSQR_Trilinos_TsqrFactory_hpp
