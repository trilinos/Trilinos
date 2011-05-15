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

#ifndef __TSQR_TestUtils_hpp
#define __TSQR_TestUtils_hpp
///
/// \file Tsqr_TestUtils.hpp
/// \brief Utilities for testing various TSQR components.
/// \author Mark Hoemmen
///

// The header file included below includes Kokkos_ConfigDefs.hpp, so
// this file doesn't have to.  Avoiding unnecessary includes can save
// build time, since the file doesn't have to be opened and the
// preprocessor doesn't have to look for the ifndef and define
// directives.
#include <Kokkos_DefaultNode.hpp>

namespace Teuchos {
  // Forward declaration of Teuchos::Comm, so that we can use
  // RCP<Comm<int> > as the argument of methods defined in this header
  // file, without needing to include Teuchos_Comm.hpp.
  template<class Ordinal>
  class Comm;
}

namespace TSQR {
  namespace Test {

    /// \fn getValidNodeParameters
    /// \brief Return valid parameter list with default values.
    ///
    /// The returned parameter list corresponds to the given NodeType
    /// (Kokkos Node type), and is meant to be the input argument of
    /// \c getNode().
    template<class NodeType>
    Teuchos::RCP<Teuchos::ParameterList> getValidNodeParameters ();

#ifdef HAVE_KOKKOS_TBB
    //
    // Specialization for TBBNode: 
    // - "Num Threads" (int) option, defaults to -1 for late init.
    //
    template<>
    Teuchos::RCP<Teuchos::ParameterList> 
    getValidNodeParameters<Kokkos::TBBNode> () 
    {
      using Teuchos::ParameterList;
      using Teuchos::parameterList;
      using Teuchos::RCP;

      RCP<ParameterList> plist = parameterList ("TBBNode");
      // -1 tells the task scheduler to control the number of threads.
      plist->set ("Num Threads", -1);
      return plist;
    }
#endif // HAVE_KOKKOS_TBB

#ifdef HAVE_KOKKOS_THREADPOOL
    //
    // Specialization for TPINode: 
    // - "Num Threads" (int) option, defaults to 0.  This number
    //   seems to be taken more seriously than TBBNode's input.
    //
    // - "Verbose" (int) option to print info about number of
    //   threads; defaults to 0.
    //
    template<>
    Teuchos::RCP<Teuchos::ParameterList> 
    getValidNodeParameters<Kokkos::TPINode> () 
    {
      using Teuchos::ParameterList;
      using Teuchos::parameterList;
      using Teuchos::RCP;

      const int numThreads = 8;

      RCP<ParameterList> plist = parameterList ("TPINode");
      plist->set ("Num Threads", numThreads);
      plist->set ("Verbose", 1);
      return plist;
    }
#endif // HAVE_KOKKOS_THREADPOOL

    //
    // Specialization for SerialNode, which takes no parameters.
    //
    template<>
    Teuchos::RCP<Teuchos::ParameterList> 
    getValidNodeParameters<Kokkos::SerialNode> () 
    {
      using Teuchos::ParameterList;
      using Teuchos::parameterList;
      using Teuchos::RCP;

      RCP<ParameterList> plist = parameterList ("SerialNode");
      return plist;
    }

    /// \brief Return a Kokkos Node instance with the given parameters.
    ///
    /// \param plist [in/out] Return value of \c
    ///   getValidNodeParameters() for the given NodeType.  This
    ///   function reserves the right to modify the input parameter
    ///   list (for example, to fill in any missing parameters with
    ///   defaults).  Do not rely on this behavior.
    template<class NodeType>
    Teuchos::RCP<const NodeType>
    getNode (const Teuchos::RCP<Teuchos::ParameterList>& plist)
    {
      return Teuchos::rcp (new NodeType (*plist));
    }

    /// \class Cons
    /// \brief Typedef container enabling iteration over compile-time type list.
    ///
    /// One can use the typedefs in a Cons to "iterate" recursively
    /// over a list of types, that is defined at compile time.
    /// CarType may be any type; these are the "values" in the type
    /// list.  CdrType must be either a Cons or a NullCons.
    ///
    /// The names Cons, Car, and Cdr come from Lisp.  (Don't write
    /// "Lisp" in all caps, unless you are referring to early versions
    /// of the language.)  A cons is a list.  If x is a cons, then
    /// (car x) returns the head of the list, and (cdr x) returns the
    /// rest of the list.
    template<class CarType, class CdrType>
    struct Cons {
      typedef CarType car_type;
      typedef CdrType cdr_type;
    };

    /// \class NullCons
    /// \brief Base case for \c Cons template recursion.
    ///
    /// NullCons doesn't need car_type or cdr_type typedefs.  Classes
    /// that iterate over a Cons type list should define
    /// specializations that make sense for a NullCons, if they want
    /// iteration to work for an empty type list (a NullCons).
    struct NullCons {};

  } // namespace Test
} // namespace TSQR

