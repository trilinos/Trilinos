// @HEADER
// ***********************************************************************
// 
//          Tpetra: Templated Linear Algebra Services Package
//                 Copyright (2008) Sandia Corporation
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
// @HEADER

#ifndef TPETRA_CONFIGDEFS_HPP
#define TPETRA_CONFIGDEFS_HPP

#ifndef __cplusplus
#define __cplusplus
#endif // ifndef __cplusplus

/* this section undefines all the things autotools defines for us that we wish it didn't. */

#ifdef PACKAGE
#undef PACKAGE
#endif // ifdef PACKAGE

#ifdef PACKAGE_NAME
#undef PACKAGE_NAME
#endif // ifdef PACKAGE_NAME

#ifdef PACKAGE_BUGREPORT
#undef PACKAGE_BUGREPORT
#endif // ifdef PACKAGE_BUGREPORT

#ifdef PACKAGE_STRING
#undef PACKAGE_STRING
#endif // ifdef PACKAGE_STRING

#ifdef PACKAGE_TARNAME
#undef PACKAGE_TARNAME
#endif // ifdef PACKAGE_TARNAME

#ifdef PACKAGE_VERSION
#undef PACKAGE_VERSION
#endif // ifdef PACKAGE_VERSION

#ifdef VERSION
#undef VERSION
#endif // ifdef VERSION

// end of undoing autoconf's work section

#include <Tpetra_config.h>
#include <Teuchos_ConfigDefs.hpp>
#include <Kokkos_ConfigDefs.hpp>

//! %Tpetra namespace
namespace Tpetra { 
  // Used in all Tpetra code that explicitly must a type (like a loop index)
  // that is used with the Teuchos::Array[View,RCP] classes.

  //! Size type for Teuchos Array objects.
  typedef Teuchos_Ordinal Array_size_type; 
}

// these make some of the macros in Tpetra_Util.hpp much easier to describe
#ifdef HAVE_TPETRA_THROW_EFFICIENCY_WARNINGS
  #define TPETRA_THROWS_EFFICIENCY_WARNINGS 1
#else
  #define TPETRA_THROWS_EFFICIENCY_WARNINGS 0
#endif

#ifdef HAVE_TPETRA_PRINT_EFFICIENCY_WARNINGS
  #define TPETRA_PRINTS_EFFICIENCY_WARNINGS 1
#else
  #define TPETRA_PRINTS_EFFICIENCY_WARNINGS 0
#endif

#ifdef HAVE_TPETRA_THROW_ABUSE_WARNINGS
  #define TPETRA_THROWS_ABUSE_WARNINGS 1
#else
  #define TPETRA_THROWS_ABUSE_WARNINGS 0
#endif

#ifdef HAVE_TPETRA_PRINT_ABUSE_WARNINGS
  #define TPETRA_PRINTS_ABUSE_WARNINGS 1
#else
  #define TPETRA_PRINTS_ABUSE_WARNINGS 0
#endif


#include <functional>

#ifndef __CUDACC__
// mem management
#include <Teuchos_ArrayView.hpp>
#include <Teuchos_ArrayRCP.hpp>
#include <Teuchos_Array.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_Tuple.hpp>
// traits classes
#include <Teuchos_OrdinalTraits.hpp>
#include <Teuchos_ScalarTraits.hpp>
#include <Teuchos_TypeNameTraits.hpp>
#include <Teuchos_NullIteratorTraits.hpp>
#include <Teuchos_SerializationTraits.hpp>
// comm
#include <Teuchos_Comm.hpp>
#include <Teuchos_CommHelpers.hpp>
// misc
#include <Teuchos_ParameterList.hpp>
#endif

//! Namespace for Tpetra classes and methods
namespace Tpetra {

  /** \brief Global size_t object. 
  
      Set at configure time, this type is intended to support scenarios where the global memory allocation is larger than that of a single node.

      Currently, it is typedefed to size_t.
   */
      
  typedef size_t global_size_t;

  /*! Local versus global allocation of Map elements */
  enum LocalGlobal {
    LocallyReplicated,  /*!< Indicates that map elements are locally replicated across all nodes */
    GloballyDistributed /*!< Indicates that map elements are globally distributed across all nodes */
  };

  /*! Return status of Map lookup */
  enum LookupStatus {
    AllIDsPresent, /*!< Indicates that all queried IDs were present in the Map */
    IDNotPresent   /*!< Indicates that at least one of the specified IDs was not present in the Map */
  };

  /*! Allocation profile for matrix/graph entries */
  enum ProfileType {
    StaticProfile,  /*!< Indicates a single, static allocation (strict and more efficient) */
    DynamicProfile  /*!< Indicates multiple, dynamic allocations (flexibile, but less efficient) */
  };

  /*! Optimize storage option */
  enum OptimizeOption {
    DoOptimizeStorage,   /*!< Indicates that storage should be optimized */
    DoNotOptimizeStorage /*!< Indicates that storage should not be optimized */
  };

  /*!  \brief Tpetra::Combine Mode enumerable type */
  /*! 
    If set to Add, existing values will be summed with new values.
    If set to Insert, new values will be inserted that don't currently exist.
    If set to Replace, existing values will be replaced with new values.

    NOTE: Add and Replace are intended for modifying values that already exist,
    but it will function correctly if those values don't already exist. (i.e.
    zero will be inserted, and then summed with or replaced by the new value.)
    However, performance may suffer. (The same goes for Insert.)
  */
  enum CombineMode {
    ADD,     /*!< Existing values will be summed with new values. */
    INSERT,  /*!< Insert new values that don't currently exist. */
    REPLACE, /*!< Existing values will be replaced with new values. */
    ABSMAX   /*!< Replacment is <tt>max( abs(old_value), abs(new_value) )</tt> */
  };

  enum EPrivateComputeViewConstructor {
    COMPUTE_VIEW_CONSTRUCTOR
  };

  enum EPrivateHostViewConstructor {
    HOST_VIEW_CONSTRUCTOR
  };

  // import Teuchos memory management classes into Tpetra
#ifndef __CUDACC__
  using Teuchos::ArrayRCP;
  using Teuchos::ArrayView;
  using Teuchos::Array;
  using Teuchos::OrdinalTraits;
  using Teuchos::ScalarTraits;
  using Teuchos::RCP;
  using Teuchos::Tuple;
  using Teuchos::Comm;
  using Teuchos::null;

  using Teuchos::outArg;
  using Teuchos::tuple;
  using Teuchos::arcp;
  using Teuchos::rcp;
  using Teuchos::rcpFromRef;
  using Teuchos::rcp_const_cast;
  using Teuchos::av_reinterpret_cast;
  using Teuchos::arcp_reinterpret_cast;

  using Teuchos::typeName;

  using Teuchos::ParameterList;
  using Teuchos::parameterList;
  using Teuchos::sublist;
#endif

  // Tpetra functor objects
  // inspired by SGI-specific project2nd, project1st
  template <class Arg1, class Arg2>
  class firstArg : std::binary_function<Arg1,Arg2,Arg1> {
    public:
    typedef Arg1 first_argument_type;
    typedef Arg2 second_argument_type;
    typedef Arg1 result_type;
    inline Arg1 operator()(const Arg1 &arg1, const Arg2 &arg2) {return arg1;}
  };

  template <class Arg1, class Arg2>
  class secondArg : std::binary_function<Arg1,Arg2,Arg2> {
    public:
    typedef Arg1 first_argument_type;
    typedef Arg2 second_argument_type;
    typedef Arg2 result_type;
    inline Arg2 operator()(const Arg1 &arg1, const Arg2 &arg2) {return arg2;}
  };

} // end of Tpetra namespace


//! Namespace for %Tpetra example classes and methods
namespace TpetraExamples {
}

namespace Tpetra {
  //! Namespace for %Tpetra Reduction/Tranformation Interface
  namespace RTI {
  }
}

namespace Tpetra {
  //! Namespace for external %Tpetra functionality
  namespace Ext {
  }
  //! Collection of matrix-matrix operations.
  namespace MatrixMatrix {
  }
}

namespace Tpetra {
  //! Sweep direction for Gauss-Seidel or Successive Over-Relaxation (SOR).
  enum ESweepDirection {
    Forward = 0,
    Backward,
    Symmetric
  };
}

#endif // TPETRA_CONFIGDEFS_HPP
