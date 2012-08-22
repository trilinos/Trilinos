//@HEADER
// ************************************************************************
// 
//          Kokkos: Node API and Parallel Node Kernels
//              Copyright (2008) Sandia Corporation
// 
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
// 
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ************************************************************************
//@HEADER

#ifndef __TSQR_Trilinos_TsqrFactory_hpp
#define __TSQR_Trilinos_TsqrFactory_hpp

/// \file TsqrFactory.hpp
/// \brief Base class for TSQR implementation instantiation
///
/// \warning TSQR users should _not_ include this file directly.

#include <Tsqr_NodeTsqrFactory.hpp>
#include <Teuchos_Comm.hpp>
#include <Tsqr_MessengerBase.hpp>
#include <Tsqr.hpp>


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
    /// \tparam NodeTsqrType The type of the intranode part of TSQR.
    /// \tparam DistTsqrType The type of the internode part of TSQR.
    ///
    /// \note Unless you need to change the interface between Trilinos
    ///   and TSQR, you don't need to do anything with TsqrFactory or
    ///   its subclasses.  Just choose the appropriate subclass of \c
    ///   \c TsqrAdaptor.  TsqrFactory and its subclasses don't have
    ///   anything to do with any of the Trilinos multivector classes.
    ///
    /// \note If you have implemented a new intranode TSQR
    ///   factorization type (NodeTsqrType), you <i>may</i> need to
    ///   create a subclass (not specialization) of TsqrFactory that
    ///   knows how to instantiate that intranode TSQR class.
    ///   Alternately, you could write NodeTsqrType so that the
    ///   provided default implementation of \c makeNodeTsqr() works.
    ///
    /// \note If you have implemented a new internode TSQR
    ///   factorization type (DistTsqrType), you <i>may</i> need to
    ///   create a subclass (not specialization) of TsqrFactory that
    ///   knows how to instantiate that internode TSQR class.
    ///   Alternately, you could write DistTsqrType so that the
    ///   provided default implementation of \c makeDistTsqr() works.
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
      typedef Tsqr<LO, S, node_tsqr_type> tsqr_type;

      /// \brief Instantiate and return the TSQR implementation.
      ///
      /// \param plist [in/out] Parameter list (keys depend on the
      ///   subclass; keys are accessed in the subclass'
      ///   makeNodeTsqr() method).  On output: On output: Missing
      ///   parameters are filled in with default values.
      ///
      /// \param nodeTsqr [out] On output, points to the
      ///   node_tsqr_type object that TSQR will use for the intranode
      ///   part of its computations.
      ///
      /// \param distTsqr [out] On output, points to the
      ///   dist_tsqr_type object that TSQR will use for the internode
      ///   part of its computations.
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

      void
      prepareTsqr

		const Teuchos::RCP<scalar_messenger_type>& messenger,

      //! Virtual destructor for memory safety of derived classes.
      virtual ~TsqrFactory () {};

    private:
      /// \brief Instantiate and return the TSQR's intranode object.
      ///
      /// \param plist [in/out] Same as the epinonymous input of 
      ///   \c makeTsqr().
      ///
      /// \return The node_tsqr_type object that TSQR will use for the
      ///   intranode part of its computations.
      ///
      /// \note For implementers: this and \c makeDistTsqr() are the
      ///   two methods to implement.  makeTsqr()'s implementation is
      ///   "generic"; it does not depend on node_tsqr_type or
      ///   dist_tsqr_type.  The implementation of makeNodeTsqr()
      ///   varies for different node_tsqr_type types.  This pattern
      ///   is the compile-time polymorphism equivalent of the
      ///   "Non-Virtual Interface" (NVI) idiom, where the "virtual"
      ///   methods (here, the methods that vary for different
      ///   template parameters) are private, and the "nonvirtual"
      ///   methods (here, the methods that are the same for different
      ///   template parameters) are part of the public interface.
      virtual Teuchos::RCP<node_tsqr_type>
      makeNodeTsqr (const Teuchos::RCP<Teuchos::ParameterList>& plist) const
      {
	return Teuchos::rcp (new node_tsqr_type (plist));
      }

      /// \brief Instantiate and return TSQR's internode object.
      ///
      /// \param messenger [in] Object used by TSQR for communicating
      ///   between MPI processes.
      ///
      /// \param plist [in/out] Same as the epinonymous input of 
      ///   \c makeTsqr().
      ///
      /// \return The dist_tsqr_type object that TSQR will use for the
      ///   internode part of its computations.
      ///
      /// \note For implementers: this and \c makeNodeTsqr() are the
      ///   two interesting methods.  makeTsqr()'s implementation is
      ///   "generic"; it does not depend on node_tsqr_type or
      ///   dist_tsqr_type.  The implementation of makeDistTsqr()
      ///   varies for different dist_tsqr_type types.
      virtual Teuchos::RCP<dist_tsqr_type>
      makeDistTsqr (const Teuchos::RCP<scalar_messenger_type>& messenger,
		    const Teuchos::RCP<Teuchos::ParameterList>& plist) const
      {
	(void) plist;
	return Teuchos::rcp (new dist_tsqr_type (messenger));
      }
    };
  } // namespace Trilinos
} // namespace TSQR

#endif // __TSQR_Trilinos_TsqrFactory_hpp
