// @HEADER
// *****************************************************************************
//          Kokkos: Node API and Parallel Node Kernels
//
// Copyright 2008 NTESS and the Kokkos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef __TSQR_Trilinos_TsqrFactory_hpp
#define __TSQR_Trilinos_TsqrFactory_hpp

/// \file TsqrFactory.hpp
/// \brief Base class for TSQR implementation instantiation
///
/// \warning TSQR users should _not_ include this file directly.

#include "Tsqr.hpp"
#include "Teuchos_Comm.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RCP.hpp"

namespace TSQR {
  namespace Trilinos {

    /// \class TsqrFactory
    /// \brief Base class for TSQR implementation instantiation
    ///
    /// Each child class of TsqrFactory know how to instantiate a
    /// particular TSQR implementation.  TsqrFactory contains all
    /// common functionality for this task.
    ///
    /// \tparam LO The (local) ordinal type used by TSQR.
    /// \tparam S The Scalar type used by TSQR; the type of the
    ///   entries of the matrices to factor.
    /// \tparam NodeTsqrType The type of the intraprocess part of TSQR.
    /// \tparam DistTsqrType The type of the interprocess part of TSQR.
    ///
    /// \note Unless you need to change the interface between Trilinos
    ///   and TSQR, you don't need to do anything with TsqrFactory or
    ///   its subclasses.  Just choose the appropriate subclass of \c
    ///   \c TsqrAdaptor.  TsqrFactory and its subclasses don't have
    ///   anything to do with any of the Trilinos multivector classes.
    ///
    /// \note If you have implemented a new intraprocess TSQR
    ///   factorization type (NodeTsqrType), you <i>may</i> need to
    ///   create a subclass (not specialization) of TsqrFactory that
    ///   knows how to instantiate that intraprocess TSQR class.
    ///   Alternately, you could write NodeTsqrType so that the
    ///   provided default implementation of makeNodeTsqr works.
    ///
    /// \note If you have implemented a new interprocess TSQR
    ///   factorization type (DistTsqrType), you <i>may</i> need to
    ///   create a subclass (not specialization) of TsqrFactory that
    ///   knows how to instantiate that interprocess TSQR class.
    ///   Alternately, you could write DistTsqrType so that the
    ///   provided default implementation of makeDistTsqr works.
    ///
    /// \note If you want to change which TSQR implementation is
    ///   invoked for a particular multivector (MV) class, you don't
    ///   need to do anything with this class, as long as the
    ///   appropriate subclass of TsqrFactory for the desired
    ///   NodeTsqrType and DistTsqrType exists.  Just change the
    ///   factory_type typedef in the specialization of \c \c
    ///   TsqrTypeAdaptor for the MV class.
    template<class LO, class S, class NodeTsqrType, class DistTsqrType>
    class TsqrFactory {
    public:
      typedef LO local_ordinal_type;
      typedef S scalar_type;
      typedef NodeTsqrType node_tsqr_type;
      typedef DistTsqrType dist_tsqr_type;

      typedef MessengerBase<S> scalar_messenger_type;
      typedef Tsqr<LO, S> tsqr_type;

      /// \brief Instantiate and return the TSQR implementation.
      ///
      /// \param plist [in/out] Parameter list (keys depend on the
      ///   subclass; keys are accessed in the subclass' makeNodeTsqr
      ///   method).  On output: On output: Missing parameters are
      ///   filled in with default values.
      ///
      /// \param nodeTsqr [out] On output, points to the
      ///   node_tsqr_type object that TSQR will use for the
      ///   intraprocess part of its computations.
      ///
      /// \param distTsqr [out] On output, points to the
      ///   dist_tsqr_type object that TSQR will use for the
      ///   interprocess part of its computations.
      ///
      /// \return The node_tsqr_type instance that implements TSQR.
      Teuchos::RCP<tsqr_type>
      makeTsqr (const Teuchos::RCP<Teuchos::ParameterList>& plist,
                Teuchos::RCP<node_tsqr_type>& nodeTsqr,
                Teuchos::RCP<dist_tsqr_type>& distTsqr)
      {
        using Teuchos::RCP;
        using Teuchos::rcp;

        nodeTsqr = makeNodeTsqr (plist);
        distTsqr = makeDistTsqr (plist);
        return rcp (new tsqr_type (nodeTsqr, distTsqr));
      }

      //! Virtual destructor for memory safety of derived classes.
      virtual ~TsqrFactory () = default;

    private:
      /// \brief Instantiate and return TSQR's intraprocess object.
      ///
      /// \param plist [in/out] Same as the epinonymous input of
      ///   makeTsqr.
      ///
      /// \return The node_tsqr_type object that TSQR will use for the
      ///   intraprocess part of its computations.
      ///
      /// \note For implementers: this and makeDistTsqr are the two
      ///   methods to implement.  makeTsqr's implementation is
      ///   "generic"; it does not depend on node_tsqr_type or
      ///   dist_tsqr_type.  The implementation of makeNodeTsqr varies
      ///   for different node_tsqr_type types.  This pattern is the
      ///   compile-time polymorphism equivalent of the "Non-Virtual
      ///   Interface" (NVI) idiom, where the "virtual" methods (here,
      ///   the methods that vary for different template parameters)
      ///   are private, and the "nonvirtual" methods (here, the
      ///   methods that are the same for different template
      ///   parameters) are part of the public interface.
      virtual Teuchos::RCP<node_tsqr_type>
      makeNodeTsqr (const Teuchos::RCP<Teuchos::ParameterList>& plist) const
      {
        return Teuchos::rcp (new node_tsqr_type (plist));
      }

      /// \brief Instantiate and return TSQR's interprocess object.
      ///
      /// \param messenger [in] Object used by TSQR for communicating
      ///   between MPI processes.
      ///
      /// \param plist [in/out] Same as the epinonymous input of
      ///   makeTsqr.
      ///
      /// \return The dist_tsqr_type object that TSQR will use for the
      ///   interprocess part of its computations.
      ///
      /// \note For implementers: this and makeNodeTsqr are the two
      ///   interesting methods.  makeTsqr's implementation is
      ///   "generic"; it does not depend on node_tsqr_type or
      ///   dist_tsqr_type.  The implementation of makeDistTsqr
      ///   varies for different dist_tsqr_type types.
      virtual Teuchos::RCP<dist_tsqr_type>
      makeDistTsqr (const Teuchos::RCP<scalar_messenger_type>& messenger,
                    const Teuchos::RCP<Teuchos::ParameterList>& plist) const
      {
        auto ret = Teuchos::rcp (new dist_tsqr_type (messenger));
        ret->setParameterList (plist);
      }
    };
  } // namespace Trilinos
} // namespace TSQR

#endif // __TSQR_Trilinos_TsqrFactory_hpp
