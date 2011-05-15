// @HEADER
// ***********************************************************************
//
//                 Anasazi: Block Eigensolvers Package
//                 Copyright (2010) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
// @HEADER

#ifndef __TSQR_Trilinos_TsqrFactory_hpp
#define __TSQR_Trilinos_TsqrFactory_hpp

/// \file TsqrFactory.hpp
///
/// \warning TSQR users should _not_ include this file directly.

#include <Tsqr_NodeTsqrFactory.hpp>
#include <Teuchos_Comm.hpp>
#include <Tsqr_MessengerBase.hpp>
#include <Tsqr.hpp>

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
    /// \note If you have implemented a new internode TSQR
    /// factorization type, you'll need to create a subclass of
    /// TsqrFactory that knows how to instantiate that internode TSQR
    /// class.
    ///
    /// \note If you want to change which TSQR implementation is
    /// invoked for a particular multivector (MV) class, you don't
    /// need to do anything with this class, as long as the
    /// appropriate subclass of TsqrFactory for the desired
    /// NodeTsqrType and DistTsqrType exists.  Just change the
    /// factory_type typedef in the specialization of TsqrTypeAdaptor
    /// for the MV class.
    template< class LO, class S, class NodeTsqrType, class DistTsqrType >
    class TsqrFactory {
    public:
      typedef NodeTsqrType                        node_tsqr_type;
      typedef Teuchos::RCP< node_tsqr_type >      node_tsqr_ptr;

      typedef Teuchos::RCP< MessengerBase< S > >  scalar_messenger_ptr;
      typedef DistTsqrType                        dist_tsqr_type;
      typedef Teuchos::RCP< dist_tsqr_type >      dist_tsqr_ptr;

      typedef Tsqr< LO, S, node_tsqr_type > tsqr_type;
      typedef Teuchos::RCP< tsqr_type >     tsqr_ptr;

      /// \brief Instantiate and return TSQR implementation
      ///
      /// Instantiate and return (through the output arguments) the
      /// two TSQR implementation objects.
      ///
      /// \param plist [in] Parameter list (keys depend on the
      ///   subclass; keys are accessed in the subclass' makeNodeTsqr() 
      ///   method)
      /// \param scalar_messenger_ptr [in] Pointer to the underlying
      ///   internode communication handler, as initialized by
      ///   TSQR::Trilinos::CommFactory.
      /// \param node_tsqr [out] On output, points to the
      ///   node_tsqr_type object that TSQR will use for the intranode
      ///   part of its computations.
      /// \param tsqr [out] On output, points to the node_tsqr_type
      ///   object that TSQR will use for the internode part of its
      ///   computations.
      virtual void
      makeTsqr (const Teuchos::ParameterList& plist,
		const scalar_messenger_ptr& messenger,
		tsqr_ptr& tsqr) const
      {
	node_tsqr_ptr nodeTsqr = makeNodeTsqr (plist);
	dist_tsqr_ptr distTsqr = makeDistTsqr (messenger, plist);
	tsqr = Teuchos::rcp (new tsqr_type (nodeTsqr, distTsqr));
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
      /// \note For implementers: this and makeDistTsqr are the two
      ///   interesting methods.  makeTsqr()'s implementation is
      ///   "generic"; it does not depend on node_tsqr_type or
      ///   dist_tsqr_type.  The implementation of makeNodeTsqr()
      ///   varies for different node_tsqr_type types.  This pattern
      ///   is the compile-time polymorphism equivalent of the
      ///   "Non-Virtual Interface" (NVI) idiom, where the "virtual"
      ///   methods (here, the methods that vary for different
      ///   template parameters) are private, and the "nonvirtual"
      ///   methods (here, the methods that are the same for different
      ///   template parameters) are part of the public interface.
      virtual node_tsqr_ptr
      makeNodeTsqr (const Teuchos::ParameterList& plist) const = 0;

      /// \brief Instantiate and return TSQR's internode object
      ///
      /// \param messenger [in] Object used by TSQR for communicating
      ///   between MPI processes
      ///
      /// \param plist [in] Same as the epinonymous input of makeTsqr()
      ///
      /// \return (Smart) pointer to the dist_tsqr_type object that
      ///   TSQR will use for the internode part of its computations
      ///
      /// \note For implementers: this and makeNodeTsqr() are the two
      ///   interesting methods.  makeTsqr()'s implementation is
      ///   "generic"; it does not depend on node_tsqr_type or
      ///   dist_tsqr_type.  The implementation of makeDistTsqr()
      ///   varies for different dist_tsqr_type types.  
      virtual dist_tsqr_ptr
      makeDistTsqr (const scalar_messenger_ptr& messenger,
		    const Teuchos::ParameterList& plist) const = 0;
    };
  } // namespace Trilinos
} // namespace TSQR

#endif // __TSQR_Trilinos_TsqrFactory_hpp
