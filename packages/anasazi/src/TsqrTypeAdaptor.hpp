#ifndef __TSQR_Trilinos_TsqrTypeAdaptor_hpp
#define __TSQR_Trilinos_TsqrTypeAdaptor_hpp

#include "Teuchos_RCP.hpp"
#include "AnasaziConfigDefs.hpp"
#include "TsqrFactory.hpp"
#include "Tsqr.hpp"

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

namespace TSQR {
  namespace Trilinos {

    class UndefinedComm {};

    /// \class TsqrTypeAdaptor
    ///
    /// \brief Mapping between multivector class MV and appropriate
    /// {intra,inter}-node TSQR classes.
    ///
    /// TsqrAdaptor has to map a multivector type to two different
    /// classes:
    ///
    /// \li node_tsqr_type, responsible for the intranode part of the
    ///   TSQR computations
    /// \li tsqr_type, responsible for the internode part of the TSQR
    ///   computations
    ///
    /// TsqrTypeAdaptor maps from the multivector type MV, to these
    /// two classes.  It also gives the appropriate TsqrFactory type
    /// to use for constructing a TsqrAdaptor. 
    ///
    /// \note Implementers who want to port TSQR to a new MV class (by
    ///   mapping to an existing TSQR implementation) should first
    ///   specialize a new TsqrTypeAdaptor class for that MV.  They
    ///   should then implement the corresponding TsqrAdaptor class.
    template< class S, class LO, class GO, class MV >
    class TsqrTypeAdaptor {
    public:
      typedef S scalar_type;
      typedef LO local_ordinal_type;
      typedef GO global_ordinal_type;
      typedef MV multivector_type;

      ///
      /// Type representing the intranode part of TSQR
      ///
      typedef TSQR::SequentialTsqr< LO, S > node_tsqr_type;

      ///
      /// Type representing the internode part of TSQR.
      /// Depends on node_tsqr_type.
      ///
      typedef TSQR::Tsqr< LO, S, node_tsqr_type > tsqr_type;

      ///
      /// Type of the TsqrFactory object that knows how to construct
      /// node_tsqr_type and tsqr_type objects.
      ///
      typedef TsqrFactory< LO, S, node_tsqr_type, tsqr_type > factory_type;

      ///
      /// Type of the (raw) communicator object used by the given
      /// multivector type.  Communicator objects are always handled
      /// via Teuchos::RCP.
      ///
      typedef UndefinedComm comm_type;
      typedef Teuchos::RCP< const comm_type > comm_ptr;
    };

  } // namespace Trilinos
} // namespace TSQR

// FIXME (mfh 15 Jul 2010) Need to finish Epetra wrappers
// #ifdef HAVE_ANASAZI_EPETRA
// #  include "TsqrTypeAdaptor_Epetra_MultiVector.hpp"
// #endif // HAVE_ANASAZI_EPETRA

#ifdef HAVE_ANASAZI_TPETRA
#  include "TsqrTypeAdaptor_Tpetra_MultiVector_SerialNode.hpp"
#  include "TsqrTypeAdaptor_Tpetra_MultiVector_TBBNode.hpp"
#endif // HAVE_ANASAZI_TPETRA

#endif // __TSQR_Trilinos_TsqrTypeAdaptor_hpp
