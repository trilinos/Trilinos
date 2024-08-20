// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TPETRA_CONFIGDEFS_HPP
#define TPETRA_CONFIGDEFS_HPP

#include "Tpetra_Details_DefaultTypes.hpp"
#include "Teuchos_ConfigDefs.hpp"

namespace Tpetra {
  // Used in all Tpetra code that explicitly must a type (like a loop index)
  // that is used with the Teuchos::Array[View,RCP] classes.

  //! Size type for Teuchos Array objects.
  typedef Teuchos_Ordinal Array_size_type;
}

// these make some of the macros in Tpetra_Util.hpp much easier to describe
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

//#ifndef __CUDACC__
// mem management
#include "Teuchos_Array.hpp" // includes ArrayRCP
#include "Teuchos_RCP.hpp"
#include "Teuchos_Tuple.hpp" // includes ArrayView
// traits classes
#include "Teuchos_OrdinalTraits.hpp"
#include "Teuchos_ScalarTraits.hpp"
#include "Teuchos_TypeNameTraits.hpp"
#include "Teuchos_NullIteratorTraits.hpp"
#include "Teuchos_SerializationTraits.hpp"
// comm
#include "Teuchos_CommHelpers.hpp"
// misc
#include "Teuchos_ParameterList.hpp"
//#endif

namespace Tpetra {

  /// \brief Global size_t object.
  ///
  /// This type is intended to support scenarios where the global
  /// memory allocation is larger than that of a single node.
  ///
  /// Currently, it is typedefed to \c size_t.
  typedef size_t global_size_t;

  /// \brief Enum for local versus global allocation of Map entries.
  ///
  /// \c LocallyReplicated means that the Map's entries are locally
  /// replicated across all processes.
  ///
  /// \c GloballyDistributed means that the Map's entries are globally
  /// distributed across all processes.
  enum LocalGlobal {
    LocallyReplicated,
    GloballyDistributed
  };

  /// \brief Return status of Map remote index lookup (getRemoteIndexList()).
  enum LookupStatus {
    AllIDsPresent, /*!< All queried indices were present in the Map */
    IDNotPresent   /*!< At least one of the specified indices was not present in the Map */
  };


/*! Optimize storage option */
  enum OptimizeOption {
    DoOptimizeStorage,   /*!< Indicates that storage should be optimized */
    DoNotOptimizeStorage /*!< Indicates that storage should not be optimized */
  };

  enum EPrivateComputeViewConstructor {
    COMPUTE_VIEW_CONSTRUCTOR
  };

  enum EPrivateHostViewConstructor {
    HOST_VIEW_CONSTRUCTOR
  };

  /// \class project1st
  /// \brief Binary function that returns its first argument.
  /// \tparam Arg1 Type of the first argument, and type of the
  ///   return value.
  /// \tparam Arg2 Type of the second argument.  It may differ from
  ///   the type of the first argument.
  ///
  /// This function object might be defined in the std namespace; it
  /// is an SGI extension to the STL.  We can't count on it living in
  /// the std namespace, but it's useful, so we define it here.
  ///
  /// If you apply <tt>using namespace std;</tt> to the global
  /// namespace, and if your STL implementation includes project1st,
  /// it will cause collisions with this definition.  We recommend for
  /// this and other reasons not to state <tt>using namespace
  /// std;</tt> in the global namespace.
  template<class Arg1, class Arg2>
  class project1st {
  public:
    typedef Arg1 first_argument_type;
    typedef Arg2 second_argument_type;
    typedef Arg1 result_type;
    Arg1 operator () (const Arg1& x, const Arg2& ) const {
      return x;
    }
  };

  /// \class project2nd
  /// \brief Binary function that returns its second argument.
  /// \tparam Arg1 Type of the first argument.
  /// \tparam Arg2 Type of the second argument, and type of the return
  ///   value.  It may differ from the type of the first argument.
  ///
  /// This function object might be defined in the std namespace; it
  /// is an SGI extension to the STL.  We can't count on it living in
  /// the std namespace, but it's useful, so we define it here.
  ///
  /// If you apply <tt>using namespace std;</tt> to the global
  /// namespace, and if your STL implementation includes project1st,
  /// it will cause collisions with this definition.  We recommend for
  /// this and other reasons not to state <tt>using namespace
  /// std;</tt> in the global namespace.
  template<class Arg1, class Arg2>
  class project2nd {
  public:
    typedef Arg1 first_argument_type;
    typedef Arg2 second_argument_type;
    typedef Arg2 result_type;
    Arg2 operator () (const Arg1& , const Arg2& y) const {
      return y;
    }
  };

} // end of Tpetra namespace


// We include this after the above Tpetra namespace declaration,
// so that we don't interfere with Doxygen's ability to find the
// Tpetra namespace declaration.
#include "Tpetra_CombineMode.hpp"


//! Namespace for %Tpetra example classes and methods
namespace TpetraExamples {
}

namespace Tpetra {
  //! Namespace for external %Tpetra functionality
  namespace Ext {
  }

  /// \brief Distributed sparse matrix-matrix multiply and add.
  ///
  /// This namespace includes functions for computing the sum or product
  /// of two distributed sparse matrices, each of which is represented
  /// as a Tpetra::CrsMatrix.
  namespace MatrixMatrix {
  }

  /// \brief Distributed sparse triple matrix product.
  ///
  /// This namespace includes functions the product of three
  /// distributed sparse matrices, each of which is represented as a
  /// Tpetra::CrsMatrix.
  namespace TripleMatrixMultiply {
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

// For backwards compatibility
namespace KokkosClassic {
  using ::Tpetra::ESweepDirection;
}


#include <Kokkos_Complex.hpp>

// Specializations of Teuchos::SerializationTraits for
// Kokkos::complex<{float,double}>.

namespace Teuchos {
  template<typename Ordinal>
  class SerializationTraits<Ordinal, ::Kokkos::complex<float> >
    : public DirectSerializationTraits<Ordinal, ::Kokkos::complex<float> >
  {};

  template<typename Ordinal>
  class SerializationTraits<Ordinal, ::Kokkos::complex<double> >
    : public DirectSerializationTraits<Ordinal, ::Kokkos::complex<double> >
  {};
} // namespace Teuchos

#endif // TPETRA_CONFIGDEFS_HPP
