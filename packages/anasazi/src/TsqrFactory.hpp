#ifndef __TSQR_Trilinos_TsqrFactory_hpp
#define __TSQR_Trilinos_TsqrFactory_hpp

#include "Teuchos_Comm.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RCP.hpp"
#include "Tsqr_MessengerBase.hpp"
#include "TsqrTrilinosMessenger.hpp"

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

namespace TSQR {
  namespace Trilinos {

    /// \class TsqrFactory
    /// \brief Base class for TSQR implementation instantiation
    ///
    /// Each child class of TsqrFactory know how to instantiate a
    /// particular TSQR implementation.  TsqrFactory contains all
    /// common functionality for this task.
    ///
    /// \note Unless you need to change the interface between Trilinos
    /// and TSQR, you don't need to do anything with TsqrFactory or
    /// its subclasses.  Just choose the appropriate subclass of
    /// TsqrAdaptor.  TsqrFactory and its subclasses don't have
    /// anything to do with any of the Trilinos multivector classes.
    ///
    /// \note If you have implemented a new intranode TSQR
    /// factorization type, you'll need to create a subclass of
    /// TsqrFactory that knows how to instantiate that intranode TSQR
    /// class.
    ///
    /// \note If you want to change which TSQR implementation is
    /// invoked for a particular multivector (MV) class, change the
    /// factory_type typedef in the specialization of TsqrTypeAdaptor
    /// for MV.
    template< class LO, class S, class NodeTsqrType, class TsqrType >
    class TsqrFactory {
    public:
      typedef Teuchos::RCP< const Teuchos::Comm<int> > comm_ptr;
      typedef Teuchos::RCP< MessengerBase< S > >       messenger_ptr;
      typedef NodeTsqrType                             node_tsqr_type;
      typedef TsqrType                                 tsqr_type;
      typedef Teuchos::RCP< node_tsqr_type >           node_tsqr_ptr;
      typedef Teuchos::RCP< tsqr_type >                tsqr_ptr;

      /// \brief Instantiate and return the objects required by TSQR
      ///
      /// Instantiate and return (through the output arguments) the
      /// three objects required by TSQR.
      ///
      /// \param plist [in] Parameter list (keys depend on the
      ///   subclass; keys are accessed in the subclass' makeNodeTsqr() 
      ///   method)
      /// \param comm_ptr [in] Pointer to the underlying internode
      ///   communication handler.
      /// \param messenger [out] On output, points to the
      ///   MessengerBase< S > object that TSQR will use for internode
      ///   communication.
      /// \param node_tsqr [out] On output, points to the
      ///   node_tsqr_type object that TSQR will use for the intranode
      ///   part of its computations.
      /// \param tsqr [out] On output, points to the node_tsqr_type
      ///   object that TSQR will use for the internode part of its
      ///   computations.
      virtual void
      makeTsqr (const Teuchos::ParameterList& plist,
		const comm_ptr& comm,
		messenger_ptr& messenger,
		node_tsqr_ptr& node_tsqr,
		tsqr_ptr& tsqr) const
      {
	messenger = makeMessenger (comm);
	node_tsqr = makeNodeTsqr (plist);
	tsqr = tsqr_ptr (new tsqr_type (*node_tsqr, messenger.get()));
      }

      TsqrFactory () {}
      virtual ~TsqrFactory () {};

    private:
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
      virtual node_tsqr_ptr
      makeNodeTsqr (const Teuchos::ParameterList& plist) const = 0;

      virtual messenger_ptr 
      makeMessenger (const comm_ptr& comm) const
      {
	using TSQR::Trilinos::TrilinosMessenger;
	return messenger_ptr (new TrilinosMessenger< S > (comm));
      }
    };
  } // namespace Trilinos
} // namespace TSQR

#include "TsqrFactory_SequentialTsqr.hpp"
// mfh 28 Jun 2010: this only does anything if HAVE_KOKKOS_TBB is
// defined.  The file included below includes Tpetra_MultiVector.hpp
// in order to get that.
#include "TsqrFactory_TbbTsqr.hpp"

#endif // __TSQR_Trilinos_TsqrFactory_hpp
