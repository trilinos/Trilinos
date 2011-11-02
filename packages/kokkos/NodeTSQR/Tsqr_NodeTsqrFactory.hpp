//@HEADER
// ************************************************************************
// 
//          Kokkos: Node API and Parallel Node Kernels
//              Copyright (2009) Sandia Corporation
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
// ************************************************************************
//@HEADER

#ifndef __TSQR_NodeTsqrFactory_hpp
#define __TSQR_NodeTsqrFactory_hpp

#include <Kokkos_ConfigDefs.hpp> // HAVE_KOKKOS_TBB

#ifdef HAVE_KOKKOS_TBB
#  include <Kokkos_TBBNode.hpp>
#  include <TbbTsqr.hpp>
#endif // HAVE_KOKKOS_TBB

#ifdef HAVE_KOKKOS_THREADPOOL
#  include <Tsqr_KokkosNodeTsqr.hpp>
#endif // HAVE_KOKKOS_THREADPOOL

#include <Kokkos_SerialNode.hpp>
#include <Tsqr_SequentialTsqr.hpp>

#include <Teuchos_ParameterList.hpp>
#include <Teuchos_ParameterListExceptions.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_ScalarTraits.hpp>
#include <Teuchos_TypeNameTraits.hpp>

#include <stdexcept>


namespace TSQR {

  /// \class NodeTsqrFactory
  /// \brief Factory for creating an instance of the right \c NodeTsqr subclass.
  /// \author Mark Hoemmen
  ///
  /// This class maps from a particular Kokkos Node type, to the
  /// corresponding NodeTsqr subclass.  It lets you construct a
  /// default parameter list for that NodeTsqr subclass, as well as an
  /// instance of the NodeTsqr subclass.  It also provides typedefs
  /// for template metaprogramming.
  template<class Node, class Scalar, class LocalOrdinal>
  class NodeTsqrFactory {
  public:
    /// \typedef node_type
    /// \brief The Kokkos Node type.
    typedef Node node_type;
    typedef Teuchos::RCP<node_type> node_ptr;

    /// \typedef node_tsqr_type
    /// \brief The NodeTsqr subclass corresponding to the Kokkos Node type.
    typedef SequentialTsqr<LocalOrdinal, Scalar> node_tsqr_type;

    /// \brief Default parameter list for intranode TSQR.
    ///
    /// \note The default implementation returns an empty (not null)
    ///   parameter list.  Each specialization for a specific Node
    ///   type redefines this method to return a parameter list
    ///   appropriate for that Node type's TSQR implementation.
    static Teuchos::RCP<const Teuchos::ParameterList>
    getDefaultParameters ()
    {
      return Teuchos::parameterList ("NodeTsqr");
    }

    /// \brief Return a pointer to the intranode TSQR implementation.
    ///
    /// \param node [in/out] Pointer to the Kokkos Node instance.
    ///
    /// \param plist [in/out] Parameter list for configuring the
    ///   NodeTsqr implementation.
    static Teuchos::RCP<node_tsqr_type>
    makeNodeTsqr (const Teuchos::RCP<node_type>& node,
		  const Teuchos::RCP<Teuchos::ParameterList>& plist)
    {
      // This method is implemented with correct behavior for those
      // Kokkos Node types for which we have implemented an intranode
      // TSQR implementation.
      TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error,
				 "TSQR is not supported on your Kokkos Node type "
				 << Teuchos::TypeNameTraits<node_type>::name()
				 << ".");
    }

    /// \brief Prepare the NodeTsqr instance for use by setting its Kokkos Node instance.
    ///
    /// Some NodeTsqr subclasses can't compute anything until they
    /// have a pointer to a Kokkos Node instance.  Call this method
    /// before invoking any computational methods of the NodeTsqr
    /// subclass instance.
    ///
    /// Precondition: ! nodeTsqr.is_null() && ! node.is_null().
    ///
    /// Postcondition: nodeTsqr->ready() == true.
    static void 
    prepareNodeTsqr (const Teuchos::RCP<node_tsqr_type>& nodeTsqr,
		     const Teuchos::RCP<node_type>& node) 
    {
      // This method is implemented with correct behavior for those
      // Kokkos Node types for which we have implemented an intranode
      // TSQR implementation.
      TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error,
				 "TSQR is not supported on your Kokkos Node type "
				 << Teuchos::TypeNameTraits<node_type>::name()
				 << ".");
    }
  };

#ifdef HAVE_KOKKOS_TBB
  //
  // Specialization of NodeTsqrFactory for Kokkos::TBBNode.
  //
  template<class Scalar, class LocalOrdinal>
  class NodeTsqrFactory<Kokkos::TBBNode, Scalar, LocalOrdinal> {
  public:
    typedef Kokkos::TBBNode node_type;
    typedef Teuchos::RCP<node_type> node_ptr;
    typedef TBB::TbbTsqr<LocalOrdinal, Scalar> node_tsqr_type;

    static Teuchos::RCP<const Teuchos::ParameterList>
    getDefaultParameters ()
    {
      using Teuchos::ParameterList;
      using Teuchos::parameterList;
      using Teuchos::RCP;

      RCP<ParameterList> params = parameterList ("NodeTsqr");
      // Create a temporary node_tsqr_type instance in order to get
      // default parameters.  The empty input parameter list will get
      // filled in with default values of missing parameters.
      node_tsqr_type nodeTsqr (params);

      return params;
    }

    static Teuchos::RCP<node_tsqr_type>
    makeNodeTsqr (const Teuchos::RCP<node_type>& node,
		  const Teuchos::RCP<Teuchos::ParameterList>& params)
    {
      (void) node;
      return Teuchos::rcp (new node_tsqr_type (params));
    }

    static void 
    prepareNodeTsqr (const Teuchos::RCP<node_tsqr_type>& nodeTsqr,
		     const Teuchos::RCP<node_type>& node) 
    {
      // TbbTsqr interacts directly with TBB and doesn't use the
      // Kokkos Node.  Thus, the TbbTsqr instance doesn't need to have
      // its Kokkos Node instance set.
      (void) nodeTsqr;
      (void) node;
    }
  };
#endif // HAVE_KOKKOS_TBB

#if defined(HAVE_KOKKOS_THREADPOOL)
  //
  // Specialization of NodeTsqrFactory for Kokkos::TPINode.
  //
  template<class Scalar, class LocalOrdinal>
  class NodeTsqrFactory<Kokkos::TPINode, Scalar, LocalOrdinal> {
  public:
    typedef Kokkos::TPINode node_type;
    typedef Teuchos::RCP<node_type> node_ptr;
    typedef KokkosNodeTsqr<LocalOrdinal, Scalar, node_type> node_tsqr_type;

    static Teuchos::RCP<const Teuchos::ParameterList>
    getDefaultParameters ()
    {
      using Teuchos::ParameterList;
      using Teuchos::parameterList;
      using Teuchos::RCP;

      RCP<ParameterList> params = parameterList ("NodeTsqr");
      // Create a temporary node_tsqr_type instance in order to get
      // default parameters.  The empty input parameter list will get
      // filled in with default values of missing parameters.
      node_tsqr_type nodeTsqr (params);

      return params;
    }

    static Teuchos::RCP<node_tsqr_type>
    makeNodeTsqr (const Teuchos::RCP<node_type>& node,
		  const Teuchos::RCP<Teuchos::ParameterList>& params)
    {
      return Teuchos::rcp (new node_tsqr_type (node, params));
    }

    static void 
    prepareNodeTsqr (const Teuchos::RCP<node_tsqr_type>& nodeTsqr,
		     const Teuchos::RCP<node_type>& node) 
    {
      // KokkosNodeTsqr needs a pointer to the Kokkos Node instance.
      nodeTsqr->setNode (node);
    }
  };
#endif // defined(HAVE_KOKKOS_THREADPOOL)

  //
  // Specialization of NodeTsqrFactory for Kokkos::SerialNode.
  //
  template<class Scalar, class LocalOrdinal>
  class NodeTsqrFactory<Kokkos::SerialNode, Scalar, LocalOrdinal> {
  public:
    typedef Kokkos::SerialNode node_type;
    typedef Teuchos::RCP<node_type> node_ptr;
    typedef SequentialTsqr<LocalOrdinal, Scalar> node_tsqr_type;

    static Teuchos::RCP<const Teuchos::ParameterList>
    getDefaultParameters ()
    {
      using Teuchos::ParameterList;
      using Teuchos::parameterList;
      using Teuchos::RCP;

      RCP<ParameterList> params = parameterList ("NodeTsqr");
      // Create a temporary node_tsqr_type instance in order to get
      // default parameters.  The empty input parameter list will get
      // filled in with default values of missing parameters.
      node_tsqr_type nodeTsqr (params);

      return params;
    }

    static Teuchos::RCP<node_tsqr_type>
    makeNodeTsqr (const Teuchos::RCP<node_type>& node,
		  const Teuchos::RCP<Teuchos::ParameterList>& params)
    {
      (void) node;
      return rcp (new node_tsqr_type (params));
    }

    static void 
    prepareNodeTsqr (const Teuchos::RCP<node_tsqr_type>& nodeTsqr,
		     const Teuchos::RCP<node_type>& node) 
    {
      // SequentialTsqr doesn't need the Kokkos Node instance.
      (void) nodeTsqr;
      (void) node;
    }
  };

} // namespace TSQR

#endif // __TSQR_NodeTsqrFactory_hpp
